
import os, shutil
from collections.abc import Callable
from threading import Thread, Semaphore, Condition, current_thread, main_thread
from subprocess import run, DEVNULL, PIPE, STDOUT, CompletedProcess
from VariantCallFixer import openVCF

import MetaCanSNPer.modules.LogKeeper as LogKeeper
import MetaCanSNPer.modules.ErrorFixes as ErrorFixes
from MetaCanSNPer.modules.DirectoryLibrary import DirectoryLibrary
from MetaCanSNPer.modules.Databases import DatabaseReader
from MetaCanSNPer.Globals import *
import MetaCanSNPer.Globals as Globals

LOGGER = LogKeeper.createLogger(__name__)


'''OOP handling of multiple processes'''
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

		def targetWrapper(target, *args, **kwargs):
			try:
				target(*args, **kwargs)
			except Exception as e:
				object.__setattr__(current_thread(), "exception", e)

		self.threads = [Thread(target=targetWrapper, args=(self.target,)+args, kwargs=kwargs, **threadKwargs) for args, kwargs in zip(self.args, self.kwargs)]
		
		for t in self.threads:
			object.__setattr__(t, "exception", None)

	def __iter__(self):
		return iter(self.threads)

	# def newFinished(self) -> bool:
	# 	"""Returns True if any thread in self.threads is dead. Also returns True if all threads """
	# 	return any(not t.is_alive() for t in self.threads if t is not None) or all(t is None for t in self.threads)

	def start(self):
		for t in self.threads:
			t.start()

	def results(self) -> list[int]:
		return [i for i in range(len(self.threads)) if not self.threads[i].is_alive() and self.threads[i].exception is None]

	def finished(self) -> bool:
		return not any(t.is_alive() for t in self.threads)

# Aligner and Mapper classes to inherit from
class ProcessWrapper:
	Lib : DirectoryLibrary
	database : DatabaseReader
	softwareName : str
	returncodes : list[list[int]]
	previousErrors : list
	category : str
	semaphore : Semaphore

	def __init__(self):
		self.semaphore = Semaphore() # Starts at 1.

	def createCommand(self, *args, **kwargs):
		"""Not implemented for the template class, check the wrapper of the
		specific software you are intending to use."""
		raise NotImplementedError(f"Called `.createCommand` on a object of a class which has not implemented it. Object class is `{type(self)}`")

	def run(self, command, log : str, pReturncodes : list[int], *args, **kwargs) -> None:
		
		try:
			n=[i for i in range(len(self.returncodes)) if self.returncodes[i] is pReturncodes]
		except:
			n="?"
		
		try:
			LOGGER.debug(f"Thread {n} - Running command: {command}")
			p : CompletedProcess = run(command.split() if type(command) is str else command, *args, **kwargs)
			LOGGER.debug(f"Thread {n} - Returned with exitcode: {p.returncode}")
			pReturncodes.append(p.returncode)
			self.handleRetCode(p.returncode, prefix=f"Thread {n} - ")

			if log is not None:
				if p.returncode == 0:
					LOGGER.debug(f"Thread {n} - Logging output to: {log}")
				else:
					LOGGER.error(f"Thread {n} - Child process encountered an issue, logging output to: {log}")

				with open(log, "w") as logFile:
					logFile.write(p.stdout.decode("utf-8"))
					logFile.write("\n")

			LOGGER.debug(f"Thread {n} - Finished!")
			self.semaphore.release()
		except Exception as e:
			e.add_note(f"Thread {n}")
			LOGGER.exception(e)

	def start(self) -> list[tuple[tuple[str,str],str]]:
		'''Starts processes in new threads. Returns information of the output of the processes, but does not ensure the processes have finished.'''

		commands, logs, outputs = self.createCommand()
		self.threadGroup = ThreadGroup(self.run, args=list(zip(commands, logs, self.returncodes)), kwargs=[{"stdout":PIPE, "stderr":STDOUT}]*len(commands), daemon=True)
		
		self.threadGroup.start()

		self.outputs = outputs
		return outputs
	
	def waitNext(self, timeout=None):
		self.semaphore.acquire(timeout=timeout)
		return self.threadGroup.results()

	def wait(self, timeout=5):
		while not self.finished() and all(t.exception is None for t in self.threadGroup):
			self.waitNext(timeout=timeout)
		if any(t.exception is not None for t in self.threadGroup):
			raise ChildProcessError(f"{self.softwareName!r} crashed/failed to start. Check logs for more details.")

	def updateWhileWaiting(self, outputDict : dict):
		finished = set()
		while not self.finished() and all(t.exception is None for t in self.threadGroup):
			for i in set(self.waitNext(timeout=1)).difference(finished):
				finished.add(i)
				key, path = self.outputs[i]
				if self.returncodes[i][-1] == 0:
					outputDict[key] = path
					LOGGER.info(f"Finished running {self.softwareName} {len(outputDict)}/{len(self.outputs)}.")
		if any(t.exception is not None for t in self.threadGroup):
			raise ChildProcessError(f"{self.softwareName!r} crashed/failed to start. Check logs for more details.")

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

	def handleRetCode(self, returncode : int, prefix : str=""):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		if returncode == 0:
			LOGGER.info("{prefix}{softwareName} finished with exitcode 0.".format(prefix=prefix, softwareName=self.softwareName))
		else:
			LOGGER.warning("{prefix}WARNING {softwareName} finished with a non zero exitcode: {returncode}".format(prefix=prefix, softwareName=self.softwareName, returncode=returncode))

	def preProcess(self, *args, **kwargs):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		pass

	def displayOutcomes(self, out=print):
		msg = [
			f"{self.softwareName!r}: Processes finished with exit codes for which there are no implemented solutions.",
			f"{'QUERY':<58}{'REFERENCE':<58}={'EXITCODE':>4}"
		]
		for (key, path), e in zip(self.outputs, self.returncodes):
			msg.append(f"{self.Lib.queryName!r:<30}{key!r:<60}={e[-1]:>4}")
		
		out("\n".join(msg))

class IndexingWrapper(ProcessWrapper):
	queryName : str
	commandTemplate : str # Should contain format tags for {target}, {ref}, {output}, can contain more.
	outFormat : str
	inFormat : str
	logFormat : str
	flags : list[str]
	format : str

	returncodes : list[list[int]]
	outputs : list[tuple[str,str]]
	"""[((refName, queryName), outputPath), ...]"""
	threadGroup : ThreadGroup
	solutions : ErrorFixes.SolutionContainer

	
	def __init__(self, lib : DirectoryLibrary, database : DatabaseReader, outputTemplate : str, flags : list[str]=[]):
		super().__init__()
		self.Lib = lib
		self.database = database
		self.queryName = self.Lib.queryName ## get name of file and remove ending
		self.outputTemplate = outputTemplate
		
		self.flags = flags
		self.returncodes = [[] for _ in range(len(database.references))]
		self.threadGroup = None
		self.solutions : ErrorFixes.SolutionContainer = ErrorFixes.get(self.softwareName)(self)
		
		if not os.path.exists(self.Lib.tmpDir > self.softwareName):
			os.mkdir(self.Lib.tmpDir.writable > self.softwareName)

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
			
			output = outDir > self.outputTemplate.format(self.formatDict)
			logfile = output+".log"

			self.formatDict["output"] = output
			self.formatDict["logFile"] = logfile

			LOGGER.info(f"Created command (if --dry-run has been specified, this will not be the true command ran):\n{self.commandTemplate.format(self.formatDict)}")
			if not Globals.DRY_RUN:
				command = self.commandTemplate.format(self.formatDict)
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
