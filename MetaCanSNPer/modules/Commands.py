



import re
from subprocess import Popen, PIPE, CompletedProcess
from typing import TextIO
from threading import Thread

from MetaCanSNPer.modules.LogKeeper import createLogger

LOGGER = createLogger(__name__)

from MetaCanSNPer.modules.Hooks import Hooks

parallelPattern = re.compile("\s*[&]\s*")
sequentialPattern = re.compile("\s*[;]\s*")
pipePattern = re.compile("\s*[|]\s*")
dumpPattern = re.compile("\s*[>]\s*")
argsPattern = re.compile("""(['].*?[']|["].*?["]|\S+)""", flags=re.MULTILINE)
whitePattern = re.compile("\s*")

class SequentialCommands: pass
class PipeCommands: pass
class DumpCommands: pass
class Commands: pass

class Command:

	commands : Commands
	logFile : TextIO
	raw : str
	category : str
	hooks : Hooks

	def __init__(self, string : str=None, commands : list[str]=None, category=None, hooks : Hooks=None, logFile : TextIO=None):
		if commands is not None:
			self.raw = " & ".join(commands)
		else:
			self.raw = string
		self.logFile = logFile
		self.category = category
		self.hooks = hooks

		def parallelFinished(eventInfo, self : Command):
			try:
				i = self.commands._list.index(eventInfo["object"])
			except:
				return
			self.hooks.trigger(f"{self.category}ProcessFinished", {"threadN" : i, "Command" : self})

		self._hook = self.hooks.addHook(f"SequentialCommands{self.category}Finished", target=parallelFinished, args=[self])

		self.commands = Commands(self.raw, category, hooks, logFile=logFile)
	
	def __del__(self):
		self.hooks.removeHook(f"SequentialCommands{self.category}Finished", self._hook)

	def __iter__(self):
		return iter(self._list)
	
	def __str__(self):
		return self.raw
	
	def __format__(self, format_spec):
		return self.raw.__format__(format_spec)

	def __getitem__(self, key):
		return self._list[key]
	
	def start(self):
		self.commands.start()

	def run(self):
		self.commands.run()

class Commands(Command):
	pattern : re.Pattern = parallelPattern
	nextType = SequentialCommands
	hooks : Hooks
	category : str
	_list : list[SequentialCommands]

	def __init__(self, string, category : str, hooks : Hooks, logFile : TextIO=None):
		self.raw = string
		self.logFile = logFile
		self.category = category
		self.hooks = hooks
		_list = [[]]
		
		for c in argsPattern.split(string.strip()):
			if self.pattern.fullmatch(c):
				_list.append([])
			else:
				_list[-1].append(c)
		
		self._list = [self.nextType("".join(l), logFile=self.logFile) for l in _list]
	
	def start(self) -> list[CompletedProcess]:
		processes = []
		
		for sc in self._list: # TODO Not actually parallel YET
			sc.start()
			processes.extend(p)
			if processes[-1].returncode != 0:
				return processes
		return processes
	
	def run(self) -> list[CompletedProcess]:
		processes = []
		
		for sc in self._list: # TODO Not actually parallel YET
			p = sc.run()
			processes.extend(p)
			if processes[-1].returncode != 0:
				return processes
		return processes
	

class DumpCommands(Commands):
	pattern = dumpPattern
	nextType = lambda string : list(filter(lambda s : whitePattern.fullmatch(s) is None, argsPattern.split(string)))
	command : list[str]
	outFile : TextIO

	def __init__(self, string, logFile : TextIO=None):
		_list = [[]]
		for c in argsPattern.split(string.strip()):
			if self.pattern.fullmatch(c):
				_list.append([])
			else:
				_list[-1].append(c)
		
		command = _list[0]
		self.command = list(filter(lambda s : whitePattern.fullmatch(s) is None, command))

		if len(_list) == 1:
			self.outFile = None
		elif len(_list) == 2:
			if len(_list[1]) != 1:
				LOGGER.exception(ValueError(f"Output can not be dumped to multiple filenames. Filenames given: {_list[1]}"))
				raise ValueError(f"Output can not be dumped to multiple filenames. Filenames given: {_list[1]}")
			self.outFile = open(_list[1][0], "ab")
		else:
			LOGGER.exception(ValueError(f"Output dumped more or less than once using '>' in one command. Command: {'>'.join(map(''.join, _list))}"))
			raise ValueError(f"Output dumped more or less than once using '>' in one command. Command: {'>'.join(map(''.join, self._list))}")
		
	def run(self, stdin=None, stdout=PIPE, stderr=PIPE, **kwargs) -> Popen:
		p : Popen = Popen(self.command, stdin=stdin, stdout=self.outFile or stdout, stderr=stderr, **kwargs)
		p.start()
		return p

class PipeCommands(Commands):
	pattern = pipePattern
	nextType = DumpCommands
	_list : list[DumpCommands]

	def run(self, **kwargs) -> CompletedProcess:
		first = object()
		first.__setattr__("stdout", None)
		processes : list[Popen] = [first]
		for i, dc in enumerate(self._list[:-1]):
			p : Popen = dc.run(stdin=processes[i].stdout, stdout=PIPE, stderr=PIPE, **kwargs)
			processes.append(p)
		processes.append( self._list[-1].run(stdin=processes[-1].stdout, stdout=PIPE, stderr=PIPE, **kwargs))

		processes[-1].wait()
		processes.pop(0) # Remove dummy object
		for i, p in enumerate(processes):
			print(f"[{self.raw.split()[0]}]", file=self.logFile)
			print(f"[{self.raw.split()[0]}.STDERR]", file=self.logFile)
			print(p.stderr, file=self.logFile)
			if i+1 == len(processes):
				print(f"[{self.raw.split()[0]}.STDOUT]", file=self.logFile)
				print(p.stdout, file=self.logFile)
			if p.returncode != 0:
				return processes[:i+1]
		return processes

class SequentialCommands(Commands):
	pattern = sequentialPattern
	nextType = PipeCommands
	_list : list[PipeCommands]
	processes : list[Popen]
	thread : Thread

	def start(self : SequentialCommands) -> list[CompletedProcess]:
		self.processes = []
		def runInSequence(self : SequentialCommands):
			for pc in self._list:
				p = pc.run()
				self.processes.extend(p)
				if self.processes[-1].returncode != 0:
					break
			self.hooks.trigger(f"SequentialCommands{self.category}Finished", {"object" : self})
		self.thread = Thread(target=runInSequence, args=[self], daemon=True)
		self.thread.start()

	def run(self) -> list[CompletedProcess]:
		self.processes = []
		for pc in self._list:
			p = pc.run()
			self.processes.extend(p)
			if self.processes[-1].returncode != 0:
				return self.processes
		return self.processes
		
class Commands(Command):
	nextType = SequentialCommands