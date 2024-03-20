
import os, shutil
from collections.abc import Callable
from threading import Semaphore, Condition, current_thread, main_thread, ThreadError
from threading import Thread as _Thread
from subprocess import run, DEVNULL, PIPE, STDOUT, CompletedProcess
from VariantCallFixer import openVCF

import MetaCanSNPer.modules.LogKeeper as LogKeeper
import MetaCanSNPer.modules.ErrorFixes as ErrorFixes
from MetaCanSNPer.modules.DirectoryLibrary import DirectoryLibrary
from MetaCanSNPer.modules.Hooks import Hooks
from MetaCanSNPer.modules.Databases import DatabaseReader
from MetaCanSNPer.Globals import *
import MetaCanSNPer.Globals as Globals

LOGGER = LogKeeper.createLogger(__name__)

class ThreadGroup: pass

class Thread(_Thread):

	group : ThreadGroup
	exception : Exception

	# def __new__(cls, group : ThreadGroup=None, **kwargs):
	# 	obj = super().__new__(cls)
	# 	obj.group = group
	# 	obj.exception = None
	# 	return obj
	
	def __init__(self, group : ThreadGroup=None, **kwargs):
		if "daemon" not in kwargs:
			kwargs["daemon"] = True
		
		super().__init__(**kwargs)
		self.group = group
		self.exception = None
	
	def run(self, *args, **kwargs):
		try:
			super().run(*args, **kwargs)
		except Exception as e:
			self.exception = e

"""OOP handling of multiple processes"""
class ThreadGroup:
	threads : list[Thread]
	commands : list[str]
	args : list[list] | list
	kwargs : list[dict] | dict
	threadKwargs : dict
	returncodes : dict[int]

	def __init__(self, target : Callable, args : list[list] | list=None, kwargs : list[dict] | dict=None, n : int=None, names : list=None, daemons : list[bool]|bool=True):
		"""
		"""
		self.target = target
		self.returncodes = {}
		
		# Error management for developers
		if n is not None:
			self.args = [args or []] * n
			self.kwargs = [kwargs or {}] * n
		elif len(args) == len(kwargs):
			self.n = len(args)
			self.args = args
			self.kwargs = kwargs
		else:
			raise TypeError("ThreadGroup created without specifying number of threads to be created.\nThe number is hinted from the length of args or kwargs. or specified with n.")

		try:
			assert type(names) not in [str, bytes, type(None)]
			assert len(names) == self.n
			self.names = names
		except:
			self.names = list(range(1, self.n+1))

		try:	len(daemons);	self.daemons = daemons
		except:	self.daemons = [daemons] * len(self.args)
		
		if any(self.n != x for x in map(len, [self.kwargs, self.names, self.daemons])):
			raise ValueError("ThreadGroup entries for all kwargs except target must be contiguous, scalar, or None. The exception is for `args` and `kwargs` if `n` is given.")

		# def targetWrapper(target, *args, **kwargs):
		# 	try:
		# 		target(*args, **kwargs)
		# 	except Exception as e:
		# 		object.__setattr__(current_thread(), "exception", e)

		self.threads = []
		for args, kwargs, name, daemon in zip(self.args, self.kwargs, self.names, self.daemons):
			self.threads.append(Thread(group=self, target=target, args=args, kwargs=kwargs, name=name, daemon=daemon))


	def __iter__(self):
		return iter(self.threads)
	
	def __getitem__(self, key):
		return self.threads[key]

	# def newFinished(self) -> bool:
	# 	"""Returns True if any thread in self.threads is dead. Also returns True if all threads """
	# 	return any(not t.is_alive() for t in self.threads if t is not None) or all(t is None for t in self.threads)

	def start(self):
		for t in self.threads:
			t.start()

	def results(self) -> list[int]:
		return [t.name for t in self.threads if not t.is_alive() and t.exception is None]

	def finished(self) -> bool:
		return not any(t.is_alive() for t in self.threads)
	

# Aligner and Mapper classes to inherit from
class ProcessWrapper:
	Lib : DirectoryLibrary
	database : DatabaseReader
	hooks : Hooks
	softwareName : str
	history : dict[int, list[int]] # History
	ignoredErrors : set
	category : str
	semaphore : Semaphore
	solutions : ErrorFixes.SolutionContainer
	skip : set
	threadGroup : ThreadGroup

	def __init__(self, lib : DirectoryLibrary, database : DatabaseReader, outputTemplate : str, hooks : Hooks=Hooks()):
		self.Lib = lib
		self.database = database
		self.queryName = self.Lib.queryName ## get name of file and remove ending
		self.outputTemplate = outputTemplate
		self.hooks = hooks
		self.threadGroup = None
		self.semaphore = Semaphore() # Starts at 1.
		self.history = {}
		self.skip = set()
		self.outputs = {}

	def createCommand(self, *args, **kwargs):
		"""Not implemented for the template class, check the wrapper of the
		specific software you are intending to use."""
		raise NotImplementedError(f"Called `.createCommand` on a object of a class which has not implemented it. Object class is `{type(self)}`")

	def run(self, command, log : str, *args, **kwargs) -> None:
		""""""
		
		try:
			this : Thread = current_thread()
			threadN, TG = int(this.name), this.group
			LOGGER.debug(f"{self.category} Thread {threadN} - Running command: {command}")
			self.hooks.trigger(f"{self.category}Progress", {"progress" : 0.0, "threadN" : threadN})

			p : CompletedProcess = run(command.split() if type(command) is str else command, *args, **kwargs)
			LOGGER.debug(f"{self.category} Thread {threadN} - Returned with exitcode: {p.returncode}")

			self.hooks.trigger(f"{self.category}Progress", {"progress" : 1.0, "threadN" : threadN})

			TG.returncodes[threadN] = p.returncode
			self.handleRetCode(p.returncode, prefix=f"Thread {threadN} - ")

			if log is not None:
				if p.returncode == 0:
					LOGGER.debug(f"{self.softwareName} Thread {threadN} - Logging output to: {log}")
				else:
					LOGGER.error(f"{self.softwareName} Thread {threadN} - Child process encountered an issue, logging output to: {log}")

				with open(log, "w") as logFile:
					logFile.write(p.stdout.decode("utf-8"))
					logFile.write("\n")

			LOGGER.debug(f"{self.category} Thread {threadN} - Finished!")
			self.hooks.trigger(f"{self.category}Finished", {"threadN":threadN, "thread":this})
			self.semaphore.release()
		except Exception as e:
			e.add_note(f"{self.category} Thread {threadN}")
			LOGGER.exception(e)

	def start(self) -> list[tuple[tuple[str,str],str]]:
		'''Starts processes in new threads. Returns information of the output of the processes, but does not ensure the processes have finished.'''

		commands, logs, outputs = self.createCommand()
		zipped, names = [], []
		for i, (c, l) in enumerate(zip(commands, logs)):
			i += 1
			if i not in self.skip:
				if i not in self.history:
					self.history[i] = []
				zipped.append((c,l))
				names.append(i)
				self.outputs[i] = outputs[i-1]
		
		self.threadGroup = ThreadGroup(self.run, args=zipped, kwargs=[{"stdout":PIPE, "stderr":PIPE}]*len(commands), names=names, daemons=True)

		self.threadGroup.start()

		return self.outputs
	
	def waitNext(self, timeout=None):
		self.semaphore.acquire(timeout=timeout)
		return map(int, self.threadGroup.results())

	def wait(self, timeout=5):
		while not self.finished() and all(t.exception is None for t in self.threadGroup):
			self.waitNext(timeout=timeout)
		if any(t.exception is not None for t in self.threadGroup):
			raise ThreadError(f"{self.softwareName!r} crashed/failed to start. Check logs for more details.")

	def updateWhileWaiting(self, outputDict : dict):
		finished = set()
		while not self.finished() and all(t.exception is None for t in self.threadGroup):
			for i in set(self.waitNext(timeout=1)).difference(finished):
				finished.add(i)
				key, path = self.outputs[i]
				try:
					if self.threadGroup.returncodes[i] == 0:
						outputDict[key] = path
					LOGGER.info(f"Finished running {self.softwareName} {len(outputDict)}/{len(self.outputs)}.")
				except:
					raise ThreadError(f"{self.softwareName!r} crashed/failed to start. Check logs for more details.")

	def finished(self):
		if self.threadGroup is None:
			raise ThreadError("{classType}.threadGroup is not initialized so it cannot be status checked with {classType}.finished()".format(classType=type(self).__name__))
		return self.threadGroup.finished()
	
	def canRun(self):
		'''Checks whether any process finished with a non-zero exitcode at the latest run. Returns True if process has not ran yet.'''
		# Expects that None!=0 is evaluated as True
		return self.threadGroup is None
	
	def hickups(self):
		'''Checks whether any process finished with a non-zero exitcode at the latest run. Returns True if process has not ran yet.'''
		# Expects that None!=0 is evaluated as True
		return any(e!=0 for e in self.threadGroup.returncodes.values())

	def fixable(self):
		'''Checks whether there is a known or suspected solution available for any errors that occured. If tried once,
		then will not show up again for that process.'''
		return any(e in self.solutions for i,e in self.threadGroup.returncodes.items() if e not in self.history[i][:-1])

	def planB(self):
		'''Runs suggested solutions for non-zero exitcodes.'''
		errors = {e:[] for e in set(self.threadGroup.returncodes.values())}
		previousUnsolved = []
		for i, e in self.threadGroup.returncodes.items():
			if e in self.history[i]:
				previousUnsolved.append((i,e))
			else:
				errors[e].append(i)
				self.history[i].append(e)
		
		if len(previousUnsolved) > 0:
			for i, e in previousUnsolved:
				LOGGER.error(f"Thread {i} of {self.softwareName} ran command: {self.threadGroup.args[i][0]!r} and returned exitcode {self.threadGroup.returncodes[i]} even after applying a fix for this exitcode.")
			raise ChildProcessError(f"{self.softwareName} process(es) returned non-zero value(s). A solution was attempted but the process returned the same exitcode. Check logs for details.")

		failed = [(e,errors[e]) for e in errors if e not in self.solutions and e != 0 and e not in self.ignoredErrors]
		if len(failed) != 0:
			for e, threads in failed:
				for i in threads:
					LOGGER.error(f"Thread {i} of {self.softwareName} ran command: {self.threadGroup.args[i][0]!r} and returned exitcode {self.threadGroup.returncodes[i]}.")
			raise ChildProcessError(f"{self.softwareName} process(es) returned non-zero value(s) which do not have an implemented fix. Check logs for details.")
		
		for e in errors:
			if e == 0 or e in self.ignoredErrors:
				for i in errors[e]:
					self.skip.add(i)
			else:
				self.solutions[e](self, errors[e])
		self.threadGroup = None

	def handleRetCode(self, returncode : int, prefix : str=""):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		if returncode == 0:
			LOGGER.info(f"{prefix} {self.softwareName} finished with exitcode 0.")
		else:
			LOGGER.warning(f"{prefix} {self.softwareName} finished with a non zero exitcode: {returncode}")

	def preProcess(self, *args, **kwargs):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		pass

	def displayOutcomes(self, out=print):
		msg = [
			f"{self.softwareName!r}: Processes finished with exit codes for which there are no implemented solutions.",
			f"{'QUERY':<30}|{'REFERENCE':<58} = {'EXITCODE':<8}"
		]
		if self.threadGroup is None:
			for (key, path), i in zip(self.outputs, sorted(self.history.keys())):
				e = self.history[i][-1]
				msg.append(f"{self.Lib.queryName:<30}|{key:<58} = {e:^8}")
		else:
			for (key, path), i in zip(self.outputs, sorted(self.threadGroup.returncodes.keys())):
				e = self.threadGroup.returncodes[i]
				msg.append(f"{self.Lib.queryName:<30}|{key:<58} = {e:^8}")
		
		out("\n".join(msg))

class IndexingWrapper(ProcessWrapper):
	queryName : str
	commandTemplate : str # Should contain format tags for {target}, {ref}, {output}, can contain more.
	outFormat : str
	inFormat : str
	logFormat : str
	flags : list[str]
	format : str

	outputs : list[tuple[str,str]]
	"""[((refName, queryName), outputPath), ...]"""
	threadGroup : ThreadGroup
	solutions : ErrorFixes.SolutionContainer

	
	def __init__(self, lib : DirectoryLibrary, database : DatabaseReader, outputTemplate : str, hooks : Hooks=Hooks(), flags : list[str]=[]):
		super().__init__(lib=lib, database=database, outputTemplate=outputTemplate, hooks=hooks)
		
		self.flags = flags
		self.threadGroup = None
		self.solutions : ErrorFixes.SolutionContainer = ErrorFixes.get(self.softwareName)(self)
		
		if not os.path.exists(self.Lib.tmpDir > self.softwareName):
			pMakeDirs(self.Lib.tmpDir.writable > self.softwareName)

		self.formatDict = {
			"tmpDir" : self.Lib.tmpDir,
			"refDir" : self.Lib.refDir,
			"query" : self.Lib.query,
			"queryName" : self.queryName,
			"options" : " ".join(self.flags),
			"outFormat" : self.outFormat
		}
	
	def createCommand(self) -> tuple[list[str], list[str], list[tuple[tuple[str,str],str]]]:
		
		outDir = self.Lib.tmpDir.create(self.softwareName, purpose="w")

		logs = []
		commands = []
		outputs = []
		for _, refName, _, _, _ in self.database.references:
			refPath = self.Lib.references[refName]
			'''Not every command needs all information, but the format function is supplied with a dictionary that has
			everything that could ever be needed.'''

			self.formatDict["refName"] = refName
			self.formatDict["refPath"] = refPath
			self.formatDict["mapPath"] = self.Lib.maps[refName]
			self.formatDict["alignmentPath"] = self.Lib.alignments[refName]
			self.formatDict["targetSNPs"] = self.Lib.targetSNPs[refName]
			
			output = outDir > self.outputTemplate.format(**self.formatDict)
			logfile = output+".log"

			self.formatDict["output"] = output
			self.formatDict["logFile"] = logfile

			LOGGER.info(f"Created command (if --dry-run has been specified, this will not be the true command ran):\n{self.commandTemplate.format(**self.formatDict)}")
			if not Globals.DRY_RUN:
				command = self.commandTemplate.format(**self.formatDict)
			else:
				if shutil.which("sleep") is not None:
					command = "sleep {}".format(random.randint(1, 5))
				elif shutil.which("timeout") is not None:
					command = "timeout {}".format(random.randint(1, 5))
				else:
					raise MissingDependancy("To perform a --dry-run, either the command 'sleep SECONDS' or 'timeout SECONDS' is needed.")
				open(output, "w").close()

			commands.append(command)
			logs.append(logfile)
			outputs.append((refName, output))
		return commands, logs, outputs

#
#	Aligner, Mapper, and SNPCaller classes to inherit from
#

class Aligner(IndexingWrapper):
	"""All that is needed to create a new implementation is to inherit from the correct software type ('Aligner' in this case) and set
	the two class attributes accordingly.
	"""
	category = "Aligners"
	
class Mapper(IndexingWrapper):
	"""All that is needed to create a new implementation is to inherit from the correct software type ('Mapper' in this case) and set
	the two class attributes accordingly.
	"""
	category = "Mappers"

class SNPCaller(IndexingWrapper):
	"""All that is needed to create a new implementation is to inherit from the correct software type ('SNPCaller' in this case) and set
	the two class attributes accordingly.
	"""
	category = "SNPCallers"
	
	def preProcess(self, force : bool=False):
		# Create VCF files that contain the to-be called SNPs

		SNPFiles = {}
		for _, genome, _, _, _ in self.database.references:
			refPath = self.Lib.references[genome]
			accession = open(refPath, "r").readline()[1:].split()[0]
			filename = f"{self.Lib.refDir.writable > pName(refPath)}.vcf"

			if force is True or not pExists(filename):
				vcfFile = openVCF(filename, "w", referenceFile=refPath)
				
				for snpID, pos, ref, alt in self.database.SNPsByGenome[genome]:
					# CHROM has to be the same as the accession id that is in the reference file.
					vcfFile.add(CHROM=accession, POS=pos, ID=snpID, REF="N", ALT="A,T,C,G")
				vcfFile.close()
				
			SNPFiles[genome] = filename
		self.Lib.settargetSNPs(SNPFiles)

