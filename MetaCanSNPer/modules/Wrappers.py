
import os, shutil, re
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

illegalPattern = re.compile(r"[^\w_ \-\.]")

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

	def __init__(self, lib : DirectoryLibrary, database : DatabaseReader, outputTemplate : str, out : dict[str,str]={}, hooks : Hooks=Hooks(), settings : dict[str,Any]={}):
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
		self.settings = settings

	def formatCommands(self, *args, **kwargs):
		"""Not implemented for the template class, check the wrapper of the
		specific software you are intending to use."""
		raise NotImplementedError(f"Called `.formatCommands` on a object of a class which has not implemented it. Object class is `{type(self)}`")

	def updateOutput(self, eventInfo, outputs : dict[int,str]):
		"""Not Implemented in the base class."""
		pass

	def createCommand(self):
		""""""
		names, commands, outputs = self.formatCommands()
		for i in range(len(names))[::-1]:
			name, outFile = outputs[i]
			if name not in self.history:
				self.history[name] = []
				self.outputs[name] = None
			if self.settings.get("saveTemp") is True and pExists(outFile):
				self.history[name].append(0)
				self.outputs[name] = outFile
				self.skip.add(name)

				self.hooks.trigger(f"{self.category}Progress", {"threadN" : name, "progress" : 1.0})
				self.hooks.trigger(f"{self.category}Finished", {"threadN" : name})

				names.pop(i), commands.pop(i), outputs.pop(i)
			
		self.hooks.removeHook(f"{self.category}ProcessFinished", self._hooksList.get("ProcessFinished"))
		self._hooksList["ProcessFinished"] = self.hooks.addHook(f"{self.category}ProcessFinished", target=self.updateOutput, args=[dict(outputs)])

		self.command = Command(commands, self.category, self.hooks, logDir=self.Lib.resultDir.create("SoftwareLogs"), names=names)

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
		return any(e!=0 for e in self.command.returncodes.values())

	def fixable(self):
		'''Checks whether there is a known or suspected solution available for any errors that occured. If tried once,
		then will not show up again for that process.'''
		return any(e in self.solutions for name,e in self.command.returncodes.items() if e not in self.history[name][:-1])

	def planB(self):
		'''Runs suggested solutions for non-zero exitcodes.'''
		errors = {e:[] for e in set(self.command.returncodes.values())}
		previousUnsolved = []
		for name, e in self.command.returncodes.items():
			if e in self.history[name]:
				previousUnsolved.append((name,e))
			else:
				errors[e].append(name)
				self.history[name].append(e)
		
		if len(previousUnsolved) > 0:
			for name, e in previousUnsolved:
				LOGGER.error(f"Thread {name} of {self.softwareName} ran command: {self.command.args[name][0]!r} and returned exitcode {self.command.returncodes[name]} even after applying a fix for this exitcode.")
			raise ChildProcessError(f"{self.softwareName} process(es) returned non-zero value(s). A solution was attempted but the process returned the same exitcode. Check logs for details.")

		failed = [(e,errors[e]) for e in errors if e not in self.solutions and e != 0 and e not in self.ignoredErrors]
		if len(failed) != 0:
			for e, threads in failed:
				for name in threads:
					LOGGER.error(f"Thread {name} of {self.softwareName} ran command: {self.command.args[name][0]!r} and returned exitcode {self.command.returncodes[name]}.")
			raise ChildProcessError(f"{self.softwareName} process(es) returned non-zero value(s) which do not have an implemented fix. Check logs for details.")
		
		for e in errors:
			if e == 0 or e in self.ignoredErrors:
				for name in errors[e]:
					self.skip.add(name)
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
			for name in sorted(self.history.keys()):
				e = self.history[name][-1] if self.history[name] != [] else ""
				key = self.database.references[name][1]
				msg.append(f"{self.Lib.queryName:<30}|{key:<28} = {e:^8}")
		else:
			for name in sorted(self.command.returncodes.keys()):
				e = self.command.returncodes.get(name, "")
				key = self.database.references[name][1]
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

	
	def __init__(self, lib : DirectoryLibrary, database : DatabaseReader, outputTemplate : str, out : dict[str,str], hooks : Hooks=Hooks(), flags : list[str]=[], settings : dict[str,Any]={}):
		super().__init__(lib=lib, database=database, outputTemplate=outputTemplate, out=out, hooks=hooks, settings=settings)
		
		self.flags = flags
		self.command = None
		self.solutions : ErrorFixes.SolutionContainer = ErrorFixes.get(self.softwareName)(self)

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
			genome = eventInfo["threadN"]
			assert self.command is eventInfo["Command"]
			self.semaphore.release()
			if self.command.returncodes[genome] not in self.solutions:
				if self.command.returncodes[genome] == 0:
					outFile = outputs[genome]
					
					if self.settings.get("saveTemp") is True:
						for fileName in os.listdir(self.Lib.tmpDir.find(f".{self.softwareName}", purpose="w") / illegalPattern.sub("-", genome)):
							try:
								os.rename(
									self.Lib.tmpDir.find(f".{self.softwareName}", purpose="w") / illegalPattern.sub("-", genome) / fileName,
									(self.Lib.tmpDir / self.softwareName ).create(illegalPattern.sub("-", genome), purpose="w") / fileName
								)
							except:
								pass

					self.outputs[genome] = outFile
					self.hooks.trigger(f"{self.category}Progress", {"threadN" : genome, "progress" : 1.0})
					self.hooks.trigger(f"{self.category}Finished", {"threadN" : genome})
				else:
					if self.settings.get("saveTemp") is True:
						genome, outFile = outputs[genome]
						for dPath in self.Lib.tmpDir / f".{self.softwareName}" / illegalPattern.sub("-", genome):
							if pExists(dPath):
								try:
									shutil.rmtree(dPath)
								except:
									pass
					self.hooks.trigger(f"{self.category}Progress", {"threadN" : genome, "progress" : None})
					self.hooks.trigger(f"{self.category}Finished", {"threadN" : genome})
		except (AssertionError) as e:
			e.add_note(f'<{self.command!r} is {eventInfo["Command"]!r} = {self.command is eventInfo["Command"]}>')
			LOGGER.exception(e, stacklevel=logging.DEBUG)
		except Exception as e:
			LOGGER.exception(e, stacklevel=logging.DEBUG)

	def formatCommands(self) -> tuple[list[Any], list[str], list[tuple[tuple[str,str],str]]]:
		
		if self.settings.get("saveTemp") is True:
			for dPath in self.Lib.tmpDir / f".{self.softwareName}":
				if pExists(dPath):
					try:
						shutil.rmtree(dPath)
					except:
						pass

			outDirTmp = self.Lib.tmpDir.create(f".{self.softwareName}", purpose="w")
			outDirFinal = self.Lib.tmpDir.create(self.softwareName, purpose="w")
		else:
			outDirTmp = self.Lib.tmpDir.create(self.softwareName, purpose="w")
			outDirFinal = self.Lib.tmpDir.create(self.softwareName, purpose="w")

		names = []
		commands = []
		outputs = []
		for _, refName, _, _, _ in self.database.references:
			if refName in self.skip: continue

			"""Not every command needs all information, but the format function is supplied with a dictionary that has
			everything that could ever be needed."""

			self.formatDict["refName"] = refName
			self.formatDict["refPath"] = self.Lib.references[refName]
			self.formatDict["mapPath"] = self.Lib.maps[refName]
			self.formatDict["alignmentPath"] = self.Lib.alignments[refName]
			self.formatDict["targetSNPs"] = self.Lib.targetSNPs[refName]
			
			output = outDirTmp.create(illegalPattern.sub("-", refName), purpose="w") / self.outputTemplate.format(**self.formatDict)

			self.formatDict["output"] = output

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

			names.append(refName)
			commands.append(command)
			outputs.append((refName, outDirFinal.create(illegalPattern.sub("-", refName), purpose="w") / self.outputTemplate.format(**self.formatDict)))
		return names, commands, outputs

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
