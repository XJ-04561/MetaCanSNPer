
import os
import logging
try:
	import CanSNPer2.modules.LogKeeper as LogKeeper

	LOGGER = LogKeeper.createLogger(__name__)
	import CanSNPer2.modules.ErrorFixes as ErrorFixes
	from CanSNPer2.modules.DirectoryLibrary import DirectoryLibrary
	from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions
except:
	import LogKeeper as LogKeeper

	LOGGER = LogKeeper.createLogger(__name__)
	from DatabaseConnection import CanSNPdbFunctions
	from DirectoryLibrary import DirectoryLibrary
	import ErrorFixes as ErrorFixes

from collections.abc import Callable
from threading import Thread, Condition
from subprocess import run, DEVNULL, PIPE, STDOUT, CompletedProcess


class Error(Exception):
	"""docstring for Error"""
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class MauveError(Error):
	"""docstring for MauveE"""
	pass

'''More OOP handling of multiple processes'''
class ThreadGroup:
	threads : list[Thread]
	commands : list[str]
	args : list[list] | list
	kwargs : list[dict] | dict
	threadKwargs : dict

	def __init__(self, target : Callable, args : list[list] | list=None, kwargs : list[dict] | dict=None, n : int=None, **threadKwargs):
		'''
			Has multiple behaviors for arguments and keyword arguments supplied to the target function:
				* args and kwargs are 2D: must be same length on axis 0, n is ignored.
				* args is 1D and kwargs is 2D: must be same length on axis 0, n is ignored.
					Each element of args goes to only one thread.
				* args is 1D and kwargs is 1D (dict): n is ignored, length of args is used.
					Each element of args goes to only one thread.
					Given kwargs is used for every thread.
				* args is 1D or 2D but has length 1 on axis 0: kwargs decides n if 2D else n is used.
					The given 
		'''
		self.target = target

		# Error management for developers
		if args is None and kwargs is None and n is None:
			raise TypeError("""ThreadGroup created without specifying number of threads to be created.
				   The number is hinted from the length of args or kwargs. or specified with n.""")
		elif args is not None and kwargs is not None:
			if len(args) != len(kwargs):
				raise TypeError("""ThreadGroup number of threads is hinted from the length of args or kwargs.
					They are not the same length. (They were {} and {}, respectively)""".format(len(args), len(kwargs)))
			self.args = args
			self.kwargs = kwargs
		elif args is not None:
			self.n = len(args)
			self.args = args
			self.kwargs = [{} for _ in range(self.n)]
		elif kwargs is not None:
			self.n = len(kwargs)
			self.kwargs = kwargs
			self.args = [[] for _ in range(self.n)]
		else:
			self.n = n
			self.args = [[] for _ in range(self.n)]
			self.kwargs = [{} for _ in range(self.n)]
		
		self.threadKwargs = threadKwargs

		self.threads = [Thread(target=self.target, args=args, kwargs=kwargs, **threadKwargs) for args, kwargs in zip(self.args, self.kwargs)]

	def newFinished(self) -> bool:
		'''Returns True if any thread in self.threads is dead. Also returns True if all threads '''
		return any(not t.is_alive() for t in self.threads if t is not None) or all(t is None for t in self.threads)

	def start(self):
		for t in self.threads:
			t.start()

	def waitNext(self, timeout=None) -> list[int]:
		Condition.wait_for(self.newFinished, timeout=timeout)

		# list comprehension to return all dead threads and also replace those same threads with 'None'.
		return [i for i in range(len(self.threads)) if not self.threads[i].is_alive() if self.threads.__setitem__(i, None) is None]

	def finished(self) -> bool:
		return all(t is None for t in self.threads)

# Aligner and Mapper classes to inherit from
class ProcessWrapper:
	softwareName : str
	returncodes : list[list[int]]
	previousErrors : list
	category : str
	
	def start(self, *args):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		pass

	def run(self, command, log : str, pReturncodes : list[int], *args, **kwargs) -> None:
		try:
			n=[i for i in range(len(self.returncodes)) if self.returncodes[i] is pReturncodes]
		except:
			n="?"
		
		LOGGER.debug("Thread {n} - Running command: {command}".format(n=n, command=command))
		p : CompletedProcess = run(command.split() if type(command) is str else command, *args, **kwargs)
		LOGGER.debug("Thread {n} - Returned with exitcode: {exitcode}".format(n=n, exitcode=p.returncode))

		if log is not None:
			LOGGER.debug("Thread {n} - Logging output to: {logpath}".format(n=n, logpath=log))
			with open(log, "w") as logFile:
				logFile.write(p.stdout.read().decode("utf-8"))
				logFile.write("\n")
		pReturncodes.append(p.returncode)
		self.handleRetCode(p.returncode, prefix="Thread {n} - ".format(n=n))

		LOGGER.debug("Thread {n} - Finished!".format(n=n))

	def createCommand(self, *args):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		pass
	
	def handleRetCode(self, returncode : int, prefix : str=""):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		if returncode == 0:
			LOGGER.debug("{prefix}{softwareName} finished with exitcode 0.".format(prefix=prefix, softwareName=self.softwareName))
		else:
			LOGGER.warning("{prefix}WARNING {softwareName} finished with a non zero exitcode: {returncode}".format(prefix=prefix, softwareName=self.softwareName, returncode=returncode))
	
	def preProcess(self):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		pass

	def hickups(self):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		pass

	def fixable(self):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		pass

	def planB(self):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		pass

class IndexingWrapper(ProcessWrapper):
	Lib : DirectoryLibrary
	queryName : str
	commandTemplate : str # Should contain format tags for {target}, {ref}, {output}, can contain more.
	outFormat : str
	# inFormat : str
	logFormat : str
	kwargs : dict
	boolFlags : list[str] = []
	valueFlags : list[str] = []
	format : str

	returncodes : list[list[int]]
	threadGroup : ThreadGroup
	solutions : dict
	
	def __init__(self, lib : DirectoryLibrary, database : CanSNPdbFunctions, outputTemplate : str, kwargs : dict={}):
		self.Lib = lib
		self.database = database
		self.queryName = os.path.splitext(os.path.basename(self.Lib.getQuery()))[0] ## get name of file and remove ending
		self.outputTemplate = outputTemplate
		# self.inFormat = inFormat
		self.kwargs = kwargs
		self.returncodes = [[] for _ in range(len(lib.getReferences()))]
		self.threadGroup = None
		self.solutions : ErrorFixes.SolutionContainer = ErrorFixes.get(self.softwareName)
		
		if not os.path.exists(os.path.join(self.Lib.tmpDir, self.softwareName)):
			os.mkdir(os.path.join(self.Lib.tmpDir, self.softwareName))

		self.formatDict = {
			"tmpDir" : self.Lib.tmpDir,
			"refDir" : self.Lib.refDir,
			"query" : self.Lib.query,
			"queryName" : self.queryName,
			# "inFormat" : self.inFormat,
			"options" : " ".join([flag for flag in self.boolFlags if flag in self.kwargs and self.kwargs[flag] is True]+
			[x for flag in self.valueFlags if flag in self.kwargs for x in [flag, self.kwargs[flag]]]) # Creates list of ["--flag1", "arg1", "--flag2", "arg2", ..., "--flagN", "argN"]
		}
	
	def start(self) -> list[tuple[tuple[str,str],str]]:
		'''Starts alignment processes in new threads. Returns information of the output of the processes, but does not ensure the processes have finished.'''

		commands, logs, outputs = self.createCommand()
		self.threadGroup = ThreadGroup(self.run, args=zip(commands, logs, self.returncodes), stdout=PIPE, stderr=STDOUT)

		self.threadGroup.start()
		
		return outputs
	
	def waitNext(self, timeout=None):
		return self.threadGroup.waitNext(timeout=timeout)

	def wait(self, timeout=5):
		while not self.threadGroup.finished():
			self.threadGroup.waitNext(timeout=timeout)
	
	def createCommand(self) -> tuple[list[str], list[str], list[tuple[tuple[str,str],str]]]:
		
		logs = []
		commands = []
		outputs = []
		for ref, refPath in self.Lib.getReferences():
			refName, refFormat = os.path.splitext(os.path.basename(refPath))

			output = self.outputTemplate.format(target=self.queryName, ref=refName, format=self.outFormat)
			logfile = os.path.join(self.Lib.logDir, output + ".{software}.log".format(software=self.softwareName))
			output = os.path.join(self.Lib.tmpDir, self.category, output)

			os.mkdir(os.path.join(self.Lib.tmpDir, self.softwareName))

			'''Not every command needs all information, but the format function is supplied with a dictionary that has
			everything that could ever be needed.'''
			
			self.formatDict["refName"] = refName
			self.formatDict["refPath"] = refPath
			self.formatDict["indexPath"] = self.Lib.indexed[(self.queryName, refName)] if (self.queryName, refName) in self.Lib.indexed else None
			self.formatDict["SNPs"] = self.Lib.SNPs[refName] if refName in self.Lib.SNPs else None
			self.formatDict["output"] = output
			self.formatDict["logfile"] = logfile

			command = self.commandTemplate.format(self.formatDict)

			'''
			Adds on all the extra flags and arguments defined in the dict self.kwargs.
			For flags that do not require a followup argument, set the value of the keyword as True.
			'''
			command = " ".join([command] + ["--{key} {value}".format(key=key, value=value if value is not True else "") for key, value in self.kwargs])

			commands.append(command)
			logs.append(logfile)
			outputs.append(((self.queryName, refName), output))
		return commands, logs, outputs

	def finished(self):
		if self.threadGroup is None:
			raise ValueError("{classType}.threadGroup is not initialized so it cannot be status checked with {classType}.finished()".format(classType=type(self).__name__))
		return self.threadGroup.finished()
	
	def hickups(self):
		'''Checks whether any process finished with a non-zero exitcode at the latest run.'''
		return any(e[-1]!=0 for e in self.returncodes)

	def fixable(self):
		'''Checks whether there is a known or suspected solution available for any errors that occured. If tried once,
		then will not show up again for that process.'''
		return any(e[-1] in self.solutions for e in self.returncodes)

	def planB(self):
		'''Runs suggested solutions for non-zero exitcodes.'''
		current = [e[-1] if e[-1] not in e[:-1] else 0 for e in self.returncodes ]
		errors = set(current)

		for e in errors:
			if e in self.solutions:
				failedThreads = [i for i, e2 in enumerate(current) if e == e2]
				self.solutions[e](self, failedThreads)


'''
	Aligner, Mapper, and Caller classes to inherit from
'''

class Aligner(IndexingWrapper):
	category = "Indexes"
	pass
	
class Mapper(IndexingWrapper):
	category = "Indexes"
	pass

class SNPCaller(IndexingWrapper):
	category = "SNPs"
	pass
