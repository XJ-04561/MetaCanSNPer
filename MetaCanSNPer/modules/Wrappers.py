
import os, shutil
from typing import Any, Iterable, Callable
from threading import Semaphore, ThreadError

import MetaCanSNPer.modules.LogKeeper as LogKeeper
import MetaCanSNPer.modules.ErrorFixes as ErrorFixes
from MetaCanSNPer.modules.DirectoryLibrary import DirectoryLibrary
from MetaCanSNPer.modules.Hooks import Hooks
from MetaCanSNPer.modules.Databases import DatabaseReader
from MetaCanSNPer.modules.Commands import Command
from MetaCanSNPer.Globals import *
import MetaCanSNPer.Globals as Globals

LOGGER = LogKeeper.createLogger(__name__)

# Aligner and Mapper classes to inherit from
class ProcessWrapper:
	Lib : DirectoryLibrary
	database : DatabaseReader
	hooks : Hooks
	softwareName : str
	history : dict[int, list[int]] # History
	ignoredErrors : set
	category : str
	outputs : dict[str, str]
	semaphore : Semaphore
	solutions : ErrorFixes.SolutionContainer
	skip : set
	command : Command
	_hooksList : dict

	def __init__(self, lib : DirectoryLibrary, database : DatabaseReader, outputTemplate : str, out : dict[str,str]={}, hooks : Hooks=Hooks()):
		self.Lib = lib
		self.database = database
		self.queryName = self.Lib.queryName ## get name of file and remove ending
		self.outputTemplate = outputTemplate
		self.hooks = hooks
		self.command = None
		self.semaphore = Semaphore() # Starts at 1.
		self.history = {}
		self.skip = set()
		self.outputs = out
		self._hooksList = {}

	def formatCommands(self, *args, **kwargs):
		"""Not implemented for the template class, check the wrapper of the
		specific software you are intending to use."""
		raise NotImplementedError(f"Called `.formatCommands` on a object of a class which has not implemented it. Object class is `{type(self)}`")

	def updateOutput(self, eventInfo, outputs : dict[int,str]):
		"""Not Implemented in the base class."""
		pass

	def createCommand(self):
		""""""
		names, commands, logs, outputs = self.formatCommands()
		for name in names:
			if name not in self.history:
				self.history[name] = []
				self.outputs[name] = None
		
		self.hooks.removeHook(f"{self.category}ProcessFinished", self._hooksList.get("ProcessFinished"))
		self._hooksList["ProcessFinished"] = self.hooks.addHook(f"{self.category}ProcessFinished", target=self.updateOutput, args=[dict(zip(names, outputs))])

		self.command = Command(commands, self.category, self.hooks, logFiles=logs, names=names)

	def start(self) -> list[tuple[tuple[str,str],str]]:
		"""Starts processes in new threads. Returns information of the output of the processes, but does not ensure the processes have finished."""
		self.createCommand()
		self.command.start()
	
	def run(self) -> list[tuple[tuple[str,str],str]]:
		"""Runs processes in this thread."""
		self.createCommand()
		self.command.run()
	
	def waitNext(self, timeout=None):
		self.semaphore.acquire(timeout=timeout)

	def wait(self, timeout=5):
		finished = 0
		while not self.finished():
			self.waitNext(timeout=timeout)
			if finished < (finished := len(self.command.returncodes)):
				LOGGER.info(f"Finished running {self.softwareName} {finished}/{len(self.command)}.")
		if any(r is None for r in self.command.returncodes.values()):
			raise ThreadError(f"{self.softwareName!r} crashed/failed to start. Check logs for more details.")

	def finished(self):
		if self.command is None:
			raise ThreadError("{classType}.command is not initialized so it cannot be status checked with {classType}.finished()".format(classType=type(self).__name__))
		return len(self.command.returncodes) == len(self.command)
	
	def canRun(self):
		'''Checks whether any process finished with a non-zero exitcode at the latest run. Returns True if process has not ran yet.'''
		# Expects that None!=0 is evaluated as True
		return self.command is None
	
	def hickups(self):
		'''Checks whether any process finished with a non-zero exitcode at the latest run. Returns True if process has not ran yet.'''
		# Expects that None!=0 is evaluated as True
		return self.command.returncodes == {} or any(e!=0 for e in self.command.returncodes.values())

	def fixable(self):
		'''Checks whether there is a known or suspected solution available for any errors that occured. If tried once,
		then will not show up again for that process.'''
		return any(e in self.solutions for i,e in self.command.returncodes.items() if e not in self.history[i][:-1])

	def planB(self):
		'''Runs suggested solutions for non-zero exitcodes.'''
		errors = {e:[] for e in set(self.command.returncodes.values())}
		previousUnsolved = []
		for i, e in self.command.returncodes.items():
			if e in self.history[i]:
				previousUnsolved.append((i,e))
			else:
				errors[e].append(i)
				self.history[i].append(e)
		
		if len(previousUnsolved) > 0:
			for i, e in previousUnsolved:
				LOGGER.error(f"Thread {i} of {self.softwareName} ran command: {self.command.args[i][0]!r} and returned exitcode {self.command.returncodes[i]} even after applying a fix for this exitcode.")
			raise ChildProcessError(f"{self.softwareName} process(es) returned non-zero value(s). A solution was attempted but the process returned the same exitcode. Check logs for details.")

		failed = [(e,errors[e]) for e in errors if e not in self.solutions and e != 0 and e not in self.ignoredErrors]
		if len(failed) != 0:
			for e, threads in failed:
				for i in threads:
					LOGGER.error(f"Thread {i} of {self.softwareName} ran command: {self.command.args[i][0]!r} and returned exitcode {self.command.returncodes[i]}.")
			raise ChildProcessError(f"{self.softwareName} process(es) returned non-zero value(s) which do not have an implemented fix. Check logs for details.")
		
		for e in errors:
			if e == 0 or e in self.ignoredErrors:
				for i in errors[e]:
					self.skip.add(i)
			else:
				self.solutions[e](self, errors[e])
		self.command = None

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
		if self.command is None:
			for i in sorted(self.history.keys()):
				e = self.history[i][-1] if self.history[i] != [] else ""
				key, _ = self.outputs[i] if i in self.outputs else ("", None)
				msg.append(f"{self.Lib.queryName:<30}|{key:<28} = {e:^8}")
		else:
			for i in sorted(self.command.returncodes.keys()):
				e = self.command.returncodes.get(i, "")
				key, _ = self.outputs[i] if i in self.outputs else ("", None)
				msg.append(f"{self.Lib.queryName:<30}|{key:<28} = {e:^8}")
		
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
	solutions : ErrorFixes.SolutionContainer

	
	def __init__(self, lib : DirectoryLibrary, database : DatabaseReader, outputTemplate : str, out : dict={}, hooks : Hooks=Hooks(), flags : list[str]=[]):
		super().__init__(lib=lib, database=database, outputTemplate=outputTemplate, hooks=hooks)
		
		self.flags = flags
		self.command = None
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
	
	def updateOutput(self, eventInfo, outputs: dict[int, str]):
		
		try:
			name = eventInfo["threadN"]
			assert self.command[name] is eventInfo["object"]
			self.semaphore.release()
			if self.command.returncodes[name] not in self.solutions:
				if self.command.returncodes[name] == 0:
					self.outputs[self.database.references[name][1]] = outputs[name]
					self.hooks.trigger(f"{self.category}Progress", {"threadN" : name, "progress" : 1.0})
					self.hooks.trigger(f"{self.category}Finished", {"threadN" : name})
				else:
					self.hooks.trigger(f"{self.category}Progress", {"threadN" : name, "progress" : None})
					self.hooks.trigger(f"{self.category}Finished", {"threadN" : name})
		except (AssertionError) as e:
			e.add_note(f'<{self.command[eventInfo["threadN"]]!r} is {eventInfo["object"]!r} = {self.command[eventInfo["threadN"]] is eventInfo["object"]}>')
			LOGGER.exception(e, stacklevel=logging.DEBUG)
		except Exception as e:
			LOGGER.exception(e, stacklevel=logging.DEBUG)

	def formatCommands(self) -> tuple[list[Any], list[str], list[str], list[tuple[tuple[str,str],str]]]:
		
		outDir = self.Lib.tmpDir.create(self.softwareName, purpose="w")

		names = []
		logs = []
		commands = []
		outputs = []
		for i, (_, refName, _, _, _) in enumerate(self.database.references):
			if i in self.skip: continue

			refPath = self.Lib.references[refName]
			'''Not every command needs all information, but the format function is supplied with a dictionary that has
			everything that could ever be needed.'''

			self.formatDict["refName"] = refName
			self.formatDict["refPath"] = refPath
			self.formatDict["mapPath"] = self.Lib.maps[refName]
			self.formatDict["alignmentPath"] = self.Lib.alignments[refName]
			self.formatDict["targetSNPs"] = self.Lib.targetSNPs[refName]
			
			output = outDir > self.outputTemplate.format(**self.formatDict)
			logfile = self.Lib.resultDir > f"{self.softwareName}_{refName}.log"

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

			names.append(i)
			commands.append(command)
			logs.append(logfile)
			outputs.append((refName, output))
		return names, commands, logs, outputs

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
