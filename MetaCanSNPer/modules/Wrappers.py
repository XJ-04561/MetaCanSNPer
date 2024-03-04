
import os
import logging
try:
	import MetaCanSNPer.modules.LogKeeper as LogKeeper

	LOGGER = LogKeeper.createLogger(__name__)
	import MetaCanSNPer.modules.ErrorFixes as ErrorFixes
	from MetaCanSNPer.modules.DirectoryLibrary import DirectoryLibrary
	from MetaCanSNPer.modules.Databases import DatabaseReader
	from MetaCanSNPer.modules.VCFhandler import CreateVCF
except:
	import LogKeeper as LogKeeper

	LOGGER = LogKeeper.createLogger(__name__)
	import ErrorFixes as ErrorFixes
	from DirectoryLibrary import DirectoryLibrary
	from Databases import DatabaseReader
	from VCFhandler import CreateVCF

from collections.abc import Callable
from threading import Thread, Condition
from subprocess import run, DEVNULL, PIPE, STDOUT, CompletedProcess


'''
	All that is needed to create a new implementation is to inherit from the correct software type ('Aligner' in this case) and set
	the two class attributes accordingly. See 'Timeout' and 'Sleep' for a minimalist example.
'''

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
	Lib : DirectoryLibrary
	database : DatabaseReader
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
	queryName : str
	commandTemplate : str # Should contain format tags for {target}, {ref}, {output}, can contain more.
	outFormat : str
	# inFormat : str
	logFormat : str
	flags : list[str]
	format : str

	returncodes : list[list[int]]
	outputs : list[tuple[tuple[str,str],str]]
	"""[((QUERY_NAME, REFERENCE_NAME), OUTPUT_PATH), ...]"""
	threadGroup : ThreadGroup
	solutions : dict

	
	def __init__(self, lib : DirectoryLibrary, database : DatabaseReader, outputTemplate : str, flags : list[str]=[]):
		self.Lib = lib
		self.database = database
		self.queryName = self.Lib.queryName ## get name of file and remove ending
		self.outputTemplate = outputTemplate
		
		self.flags = flags
		self.returncodes = [[] for _ in range(len(lib.getReferences()))]
		self.threadGroup = None
		self.solutions : ErrorFixes.SolutionContainer = ErrorFixes.get(self.softwareName)
		
		if not os.path.exists(os.path.join(self.Lib.tmpDir, self.softwareName)):
			os.mkdir(os.path.join(self.Lib.tmpDir, self.softwareName))

		self.formatDict = {
			"tmpDir" : self.Lib.tmpDir,
			"refDir" : self.Lib.refDir,
			"query" : " ".join(self.Lib.query),
			"queryName" : self.queryName,
			"options" : " ".join(self.flags)
		}
	
	def start(self) -> list[tuple[tuple[str,str],str]]:
		'''Starts alignment processes in new threads. Returns information of the output of the processes, but does not ensure the processes have finished.'''

		commands, logs, outputs = self.createCommand()
		self.threadGroup = ThreadGroup(self.run, args=zip(commands, logs, self.returncodes), stdout=PIPE, stderr=STDOUT)

		self.threadGroup.start()

		self.outputs = outputs
		return outputs
	
	def waitNext(self, timeout=None):
		return self.threadGroup.waitNext(timeout=timeout)

	def wait(self, timeout=5):
		while not self.threadGroup.finished():
			self.threadGroup.waitNext(timeout=timeout)
	
	def createCommand(self) -> tuple[list[str], list[str], list[tuple[tuple[str,str],str]]]:
		
		outDir = self.Lib.tmpDir.forceFind(self.softwareName)

		logs = []
		commands = []
		outputs = []
		for _, refName, _, _, _ in self.database.references:
			refPath = self.Lib.references[refName]
			'''Not every command needs all information, but the format function is supplied with a dictionary that has
			everything that could ever be needed.'''

			self.formatDict["refName"] = refName
			self.formatDict["refPath"] = refPath
			self.formatDict["mapPath"] = self.Lib.maps[(self.queryName, refName)] if (self.queryName, refName) in self.Lib.maps else None
			self.formatDict["alignmentPath"] = self.Lib.alignments[(self.queryName, refName)] if (self.queryName, refName) in self.Lib.alignments else None
			self.formatDict["SNPs"] = self.Lib.SNPs[refName] if refName in self.Lib.SNPs else None
			
			output = self.outputTemplate.format(self.formatDict)
			logfile = outDir > output+".log"
			output = outDir > output

			self.formatDict["output"] = output

			command = self.commandTemplate.format(self.formatDict)

			commands.append(command)
			logs.append(logfile)
			outputs.append(((self.queryName, refName), output))
		return commands, logs, outputs

	def updateWhileWaiting(self, outputDict : dict):
		while not self.finished():
			finished = self.waitNext()
			for i in finished:
				key, path = self.output[i]
				if self.returncodes[i][-1] == 0:
					outputDict[key] = path
					LOGGER.info("Finished running {software} {n}/{N}.".format(software=self.softwareName, n=len(outputDict), N=len(self.output)))

	def finished(self):
		if self.threadGroup is None:
			raise ValueError("{classType}.threadGroup is not initialized so it cannot be status checked with {classType}.finished()".format(classType=type(self).__name__))
		return self.threadGroup.finished()
	
	def hickups(self):
		'''Checks whether any process finished with a non-zero exitcode at the latest run. Returns True if process has not ran yet.'''
		return len(self.returncodes[0])==0 or any(e[-1]!=0 for e in self.returncodes)

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
	
	def preProcess(self, force : bool=False):
		# Create VCF files that contain the to-be called SNPs

		SNPFiles = {}
		for _, genome, _, _, _ in self.database.references:
			refPath = self.Lib.references[genome]
			accession = open(refPath, "r").readline()[1:].split()[0]
			filename = "{ref}.vcf".format(ref=os.path.join(self.Lib.refDir.writable, os.path.basename(refPath)))

			if force is True or not os.path.exists(filename):
				vcfFile = CreateVCF(filename, referenceFile=refPath)
				SNPs = self.database.SNPsByGenome[genome]
				
				for snpID, pos, ref, alt in SNPs:
					# CHROM has to be the same as the accession id that is in the reference file.
					vcfFile.add(CHROM=accession, POS=pos, ID=snpID, REF="N", ALT="A,T,C,G")
				vcfFile.close()
				
			SNPFiles[genome] = filename
		self.Lib.setSNPfiles(SNPFiles)
