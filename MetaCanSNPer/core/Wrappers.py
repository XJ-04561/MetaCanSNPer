
import os, shutil, re
from typing import Any, Iterable, Callable, TypeVar
from threading import Semaphore, ThreadError

import MetaCanSNPer.core.ErrorFixes as ErrorFixes
from MetaCanSNPer.core.DirectoryLibrary import DirectoryLibrary
from MetaCanSNPer.core.Hooks import Hooks
from MetaCanSNPer.modules.Database import MetaCanSNPerDatabase
from MetaCanSNPer.core.Commands import Command
from MetaCanSNPer.Globals import *
import MetaCanSNPer.Globals as Globals

_T = TypeVar("_T")
illegalPattern = re.compile(r"[^\w_ \-\.]")

class MissingDependancy(Exception): pass

# Aligner and Mapper classes to inherit from
class ProcessWrapper(Logged):
	
	queryName : str
	commandTemplate : str # Should contain format tags for {target}, {ref}, {output}, can contain more.
	outFormat : str
	inFormat : str
	logFormat : str
	flags : list[str]
	format : str

	Lib : DirectoryLibrary
	database : MetaCanSNPerDatabase
	hooks : Hooks
	softwareName : str = None
	history : dict[int, list[int]] # History
	ignoredErrors : set
	category : str
	outputs : dict[str, str]
	semaphore : Semaphore
	solutions : ErrorFixes.SolutionContainer
	skip : set
	command : Command

	@ClassProperty
	def subclasses(self):
		if isinstance(self, type):
			return {subClass.__name__.lower():subClass for subClass in self.__subclasses__()} | {subClass.softwareName.lower():subClass for subClass in self.__subclasses__() if hasattr(subClass, "softwareName")}
		else:
			return {subClass.__name__.lower():subClass for subClass in type(self).__subclasses__()} | {subClass.softwareName.lower():subClass for subClass in type(self).__subclasses__() if hasattr(subClass, "softwareName")}

	_hooksList : dict

	def __init__(self, lib : DirectoryLibrary, database : MetaCanSNPerDatabase, outputTemplate : str, out : dict[str,str], hooks : Hooks=Hooks(), flags : list[str]=[], settings : dict[str,Any]={}):
		
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
		self.flags = flags
		self.settings = settings
		
		self.solutions : ErrorFixes.SolutionContainer = ErrorFixes.get(self.softwareName)(self)

		self.formatDict = {
			"tmpDir" : self.Lib.tmpDir.writable,
			"query" : self.Lib.query,
			"queryName" : self.queryName,
			"options" : " ".join(self.flags),
			"outFormat" : self.outFormat
		}
	
	@classmethod
	def get(cls : _T, name : str):
		selectedClass : _T = cls.subclasses.get(name.lower())
		if not selectedClass:
			nt = "\n\t"
			raise MissingDependency(f"Missing dependency: {name!r}\nIf this program is not available to you, consider using one of the following instead:\n{nt.join(map(*this.softwareName, filter(*this.softwareName, cls.__subclasses__())))}")
		return selectedClass

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

				self.hooks.trigger(f"{self.category}Progress", {"name" : name, "value" : 1.0})
				self.hooks.trigger(f"{self.category}Finished", {"name" : name})

				names.pop(i), commands.pop(i), outputs.pop(i)
			
		self.hooks.removeHook(f"{self.category}Finished", self._hooksList.get("Finished"))
		self._hooksList["ProcessFinished"] = self.hooks.addHook(f"{self.category}Finished", target=self.updateOutput, args=[dict(outputs)])

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
				self.LOG.info(f"Finished running {self.softwareName} {finished}/{len(self.command)}.")
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
				self.LOG.error(f"Thread {name} of {self.softwareName} ran command: {self.command.args[name][0]!r} and returned exitcode {self.command.returncodes[name]} even after applying a fix for this exitcode.")
			raise ChildProcessError(f"{self.softwareName} process(es) returned non-zero value(s). A solution was attempted but the process returned the same exitcode. Check logs for details.")

		failed = [(e,errors[e]) for e in errors if e not in self.solutions and e != 0 and e not in self.ignoredErrors]
		if len(failed) != 0:
			for e, threads in failed:
				for name in threads:
					self.LOG.error(f"Thread {name} of {self.softwareName} ran command: {self.command.args[name][0]!r} and returned exitcode {self.command.returncodes[name]}.")
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
			self.LOG.info(f"{prefix} {self.softwareName} finished with exitcode 0.")
		else:
			self.LOG.warning(f"{prefix} {self.softwareName} finished with a non zero exitcode: {returncode}")

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
				msg.append(f"{self.Lib.queryName:<30}|{name:<28} = {e:^8}")
		else:
			for name in sorted(self.command.returncodes.keys()):
				e = self.command.returncodes.get(name, "")
				msg.append(f"{self.Lib.queryName:<30}|{name:<28} = {e:^8}")
		
		out("\n".join(msg))
	
	def updateOutput(self, eventInfo, outputs: dict[str, str]):
		
		try:
			genome = eventInfo["name"]
			if self.command is not eventInfo.get("Command"):
				self.LOG.debug(f'{self.command is not eventInfo.get("Command") =}')
				return
			self.semaphore.release()
			if self.command.returncodes[genome] not in self.solutions:
				if self.command.returncodes[genome] == 0:
					outFile = outputs[genome]
					
					if self.settings.get("saveTemp") is True:
						for fileName in os.listdir(self.Lib.tmpDir.find(f".{self.softwareName}", purpose="w") / illegalPattern.sub("-", genome)):
							try:
								os.rename(
									(self.Lib.tmpDir / f".{self.softwareName}" / illegalPattern.sub("-", genome)).find(fileName),
									(self.Lib.tmpDir / self.softwareName / illegalPattern.sub("-", genome)).writable / fileName
								)
							except:
								pass

					self.outputs[genome] = outFile
					self.hooks.trigger(f"{self.category}Progress", {"name" : genome, "value" : 1.0})
					self.hooks.trigger(f"{self.category}Finished", {"name" : genome})
				else:
					if self.settings.get("saveTemp") is True:
						outFile = outputs[genome]
						for dPath in self.Lib.tmpDir / f".{self.softwareName}" / illegalPattern.sub("-", genome):
							if pExists(dPath):
								try:
									shutil.rmtree(dPath, ignore_errors=True)
								except:
									pass
					self.hooks.trigger(f"{self.category}Progress", {"name" : genome, "value" : None})
					self.hooks.trigger(f"{self.category}Finished", {"name" : genome})
		except Exception as e:
			self.LOG.exception(e, stacklevel=logging.DEBUG)

	def formatCommands(self) -> tuple[list[Any], list[str], list[tuple[tuple[str,str],str]]]:
		
		if self.settings.get("saveTemp") is True:
			for dPath in self.Lib.tmpDir / f".{self.softwareName}":
				if pExists(dPath):
					try:
						shutil.rmtree(dPath, ignore_errors=True)
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
		for _, refName, *_ in self.database.references:
			if refName in self.skip: continue

			"""Not every command needs all information, but the format function is supplied with a dictionary that has
			everything that could ever be needed."""

			self.formatDict["refName"] = refName
			self.formatDict["refPath"] = self.Lib.references.get(refName)
			self.formatDict["mapPath"] = self.Lib.maps.get(refName)
			self.formatDict["alignmentPath"] = self.Lib.alignments.get(refName)
			self.formatDict["targetSNPs"] = self.Lib.targetSNPs.get(refName)
			
			output = outDirTmp.create(illegalPattern.sub("-", refName), purpose="w") / self.outputTemplate.format(**self.formatDict)

			self.formatDict["output"] = output

			self.LOG.info(f"Created command (if --dry-run has been specified, this will not be the true command ran):\n{self.commandTemplate.format(**self.formatDict)}")
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

class Aligner(ProcessWrapper):
	"""All that is needed to create a new implementation is to inherit from the correct software type ('Aligner' in this case) and set
	the two class attributes accordingly.
	"""
	category = "Aligners"
	
class Mapper(ProcessWrapper):
	"""All that is needed to create a new implementation is to inherit from the correct software type ('Mapper' in this case) and set
	the two class attributes accordingly.
	"""
	category = "Mappers"

class SNPCaller(ProcessWrapper):
	"""All that is needed to create a new implementation is to inherit from the correct software type ('SNPCaller' in this case) and set
	the two class attributes accordingly.
	"""
	category = "SNPCallers"
