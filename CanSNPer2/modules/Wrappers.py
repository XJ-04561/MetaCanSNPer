
import os
import logging
import CanSNPer2.modules.LogKeeper as LogKeeper

LOGGER = LogKeeper.createLogger(__name__)

from collections.abc import Callable
from threading import Thread, Condition
from subprocess import run, DEVNULL, PIPE, STDOUT, CompletedProcess

import CanSNPer2.modules.ErrorFixes as ErrorFixes
from CanSNPer2.modules.DirectoryLibrary import DirectoryLibrary

class Error(Exception):
	"""docstring for Error"""
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class MauveError(Error):
	"""docstring for MauveE"""
	pass

'''Declarations used only for typehinting'''
class CanSNPer2:
	logdir : str
	refdir : str
	get_references : Callable
	keep_going : bool
class Aligner:
	tmpdir : str
	query : str


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
		return any(not t.is_alive() for t in self.threads if t is not None)

	def start(self):
		for t in self.threads:
			t.start()

	def waitNext(self, timeout=5) -> list[int]:
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
	
	def start(self, *args):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		pass

	def run(self, command, log : str, pReturncodes : list[int], *args, **kwargs) -> None:
		p : CompletedProcess = run(command.split() if type(command) is str else command, *args, **kwargs)
		if log is not None:
			with open(log, "w") as logFile:
				logFile.write(p.stdout.read().decode("utf-8"))
				logFile.write("\n")
		pReturncodes.append(p.returncode)
		self.handleRetCode(p.returncode)

	def createCommand(self, *args):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		pass
	
	def handleRetCode(self, returncode : int):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		if returncode == 0:
			LOGGER.debug("{softwareName} finished with exitcode 0.".format(softwareName=self.softwareName))
		else:
			LOGGER.warning("WARNING {softwareName} finished with a non zero exitcode: {returncode}".format(softwareName=self.softwareName, returncode=returncode))
	
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
	logFormat : str
	returncodes : list[list[int]]
	threadGroup : ThreadGroup
	kwargs : dict
	solutions : dict
	
	def __init__(self, lib : DirectoryLibrary, outputTemplate : str, kwargs : dict={}):
		self.Lib = lib
		self.query_name = os.path.basename(self.Lib.getQuery()).rsplit(".",1)[0] ## get name of file and remove ending
		self.outputTemplate = outputTemplate
		self.kwargs = kwargs
		self.returncodes = [[] for _ in range(len(lib.getReferences()))]
		self.threadGroup = None
		self.solutions = ErrorFixes.Indexers[self.softwareName]
	
	def start(self):
		''' '''

		commands, logs, outputs = self.createCommand()
		self.threadGroup = ThreadGroup(self.run, args=zip(commands, logs, self.returncodes), stdout=PIPE, stderr=STDOUT)

		self.threadGroup.start()
		
		return outputs
	
	def wait(self, timeout=5):
		while not self.threadGroup.finished():
			self.threadGroup.waitNext(timeout=timeout)
	
	def createCommand(self) -> tuple[list[str], list[str], list[str]]:
		
		logs = []
		commands = []
		outputs = []
		for ref in self.Lib.getReferences():
			output = self.outputTemplate.format(ref=ref)
			logfile = output + "{software}.log".format(software=self.softwareName)

			command = self.commandTemplate.format(target=self.Lib.getQuery(), ref=ref, output=output)

			'''
			Adds on all the extra flags and arguments defined in the dict self.kwargs.
			For flags that do not require a followup argument, set the value of the keyword as True.
			'''
			command = " ".join([command] + ["--{key} {value}".format(key=key, value=value if value is not True else "") for key, value in self.kwargs])

			commands.append(command)
			logs.append(logfile)
			outputs.append(output)
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
			


			




# Aligner, Mapper, and Caller classes to inherit from
class Aligner(IndexingWrapper):
	pass
	
class Mapper(IndexingWrapper):
	pass

class SNPCaller(IndexingWrapper):
	pass

#
#	Implemented Aligners, Mappers, and Callers are defined here.
#

'''
	All that is needed to create a new implementation is to inherit from the correct software type (Listed above) and set
	the two class attributes accordingly. See 'Timeout' and 'Sleep' for a minimalist example.
'''

class Timeout(Aligner):
	softwareName = "timeout"
	commandTemplate = "timeout {ref} {target}"

class Sleep(Aligner):
	softwareName = "sleep"
	commandTemplate = "sleep {ref} {target}"

class ProgressiveMauve(Aligner):
	softwareName = "progressiveMauve"
	commandTemplate = "progressiveMauve {ref} {target}"

	def handleRetCode(self, returncode):
		if returncode == 0:
			LOGGER.debug("progressiveMauve finished with exitcode 0.")
		elif returncode == 11:
			LOGGER.warning("WARNING progressiveMauve finished with a exitcode: {returncode}".format(returncode=returncode))
			LOGGER.debug("This progressiveMauve error is showing up for bad genomes containing short repetitive contigs or sequence contains dashes.".format(returncode=returncode))
		elif returncode == -6:
			LOGGER.warning("Input sequence is not free of gaps, replace gaps with N and retry!!")
		else:
			LOGGER.warning("WARNING progressiveMauve finished with a non zero exitcode: {returncode}\nThe script will terminate when all processes are finished read {log} for more info".format(log=self.logFile,returncode=returncode))

Aligners : dict[Aligner] = {
	"progressiveMauve" : ProgressiveMauve
}
Mappers : dict[Mapper] = {

}
Callers : dict[SNPCaller] = {

}

