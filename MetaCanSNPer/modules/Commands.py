



import re, shutil, sys
from subprocess import Popen, PIPE, CompletedProcess
from typing import TextIO, BinaryIO
from threading import Thread

from MetaCanSNPer.modules.LogKeeper import createLogger

LOGGER = createLogger(__name__)

from MetaCanSNPer.modules.Hooks import Hooks

def bPrint(*strings, sep=b" ", end=b"\n", file : BinaryIO=None, encoding : str="utf-8"):
	file.write(sep.join(map(lambda s : s.encode("utf-8"), strings)) + end)

parallelPattern = re.compile("\s*[&]\s*")
sequentialPattern = re.compile("\s*[;]\s*")
pipePattern = re.compile("\s*[|]\s*")
dumpPattern = re.compile("\s*[>]\s*")
argsPattern = re.compile("""(['].*?[']|["].*?["]|\S+)""", flags=re.MULTILINE)
whitePattern = re.compile("\s*")

class ParallelCommands: pass
class SequentialCommands: pass
class PipeCommands: pass
class DumpCommands: pass
class Commands: pass


class Command:
	"""Triggers the event `f"{self.category}ProcessFinished"` on each finished parallel command"""


	commands : ParallelCommands

	def __init__(self, string : str=None, category=None, hooks : Hooks=None, logFiles : list[str]|str=None):
		try:
			self.raw = string
			self.logFiles = [logFiles] if type(logFiles) is str else logFiles
			self.category = category
			self.hooks = hooks
			self.returncodes = {}
			self.exceptions = {}

			def parallelFinished(eventInfo, self : Command):
				try:
					i = self.commands._list.index(eventInfo["object"])
					assert id(eventInfo["object"]) == id(self.commands._list[i]) # Unsure if needed, but ensures they are the same exact object.
					self.returncodes[i] = None if len(eventInfo["object"].processes) == 0 else eventInfo["object"].processes[-1].returncode
					self.hooks.trigger(f"{self.category}ProcessFinished", {"threadN" : i, "Command" : self} | eventInfo)
				except AssertionError:
					LOGGER.error(f"`assert id(eventInfo[\"object\"]) == id(self.commands._list[i])` did not pass.\n\t{eventInfo['object']=}\n\t{self.commands._list[i]=}")
					return
				except Exception as e:
					e.add_note(f"'parallelFinished' event exception.")
					LOGGER.exception(e)
					return

			self._hook = self.hooks.addHook(f"SequentialCommands{self.category}Finished", target=parallelFinished, args=[self])

			self.commands = ParallelCommands(string, category, hooks, logFiles=self.logFiles)
		except Exception as e:
			e.add_note(f"{type(self).__name__} failed to initialize.")
			LOGGER.exception(e)
			raise e
	
	def __len__(self):
		return len(self.commands._list)
	
	def __getitem__(self, key):
		return iter(self.commands)
	
	def __getitem__(self, key):
		return self.commands[key]
	
	def __repr__(self):
		return f"<{type(self).__name__} {hex(id(self))} raw={self.raw!r}>"
	
	def __del__(self):
		self.hooks.removeHook(f"SequentialCommands{self.category}Finished", self._hook)
	
	def start(self):
		self.commands.start()

	def run(self):
		self.commands.run()

class Commands:
	"""Only meant to be inherited from"""

	pattern : re.Pattern
	nextType : Commands
	hooks : Hooks
	category : str
	_list : list[Commands]
	logFile : TextIO
	raw : str

	def __init__(self, args, category : str, hooks : Hooks, logFile : TextIO=None):
		try:
			self.raw = "".join(args).strip()
			self.logFile = logFile
			self.category = category
			self.hooks = hooks
			_list = [[]]
			
			for c in args:
				if self.pattern.fullmatch(c):
					_list.append([])
				else:
					_list[-1].append(c)
			
			self._list = [self.nextType(l, category, hooks, logFile=self.logFile) for l in _list]
		except Exception as e:
			e.add_note(f"{type(self).__name__} failed to initialize.")
			LOGGER.exception(e)
			raise e
	
	def __iter__(self):
		return iter(self._list)
	
	def __str__(self):
		return self.raw
	
	def __repr__(self):
		return f"<{type(self).__name__} {hex(id(self))} raw={self.raw!r}>"

	def __getitem__(self, key):
		return self._list[key]

class DumpCommands(Commands):
	"""Takes string partitioning command and file to put stdout into by ">".
	If no ">" is used, then whatever `stdout=` is provided to the `.run()`-method
	is where the process stdout will be sent."""
	pattern = dumpPattern
	command : list[str]
	outFile : TextIO
	_list : list[str]

	@staticmethod
	def nextType(args, *overFlowArgs, **kwargs) -> list[str]:
		return args

	def __init__(self, args, category : str, hooks : Hooks, logFile : TextIO=None):
		try:
			super().__init__(args, category, hooks, logFile)
		
			command = self._list[0]
			self.command = list(filter(lambda s : whitePattern.fullmatch(s) is None, command))

			if len(self._list) == 1:
				self.outFile = None
			elif len(self._list) == 2:
				if len(self._list[1]) != 1:
					LOGGER.exception(ValueError(f"Output can not be dumped to multiple filenames. Filenames given: {self._list[1]}"))
					raise ValueError(f"Output can not be dumped to multiple filenames. Filenames given: {self._list[1]}")
				self.outFile = open(self._list[1][0], "ab")
			else:
				LOGGER.exception(ValueError(f"Output dumped more or less than once using '>' in one command. Command: {'>'.join(map(''.join, self._list))}"))
				raise ValueError(f"Output dumped more or less than once using '>' in one command. Command: {'>'.join(map(''.join, self._list))}")
		except Exception as e:
			e.add_note(f"{type(self).__name__} failed to initialize.")
			LOGGER.exception(e)
			raise e
		
	def run(self, stdin=None, stdout=PIPE, stderr=PIPE, **kwargs) -> Popen:
		LOGGER.debug(f"Starting {self!r}")
		ex = shutil.which(self.command[0])
		if ex is None:
			LOGGER.exception(FileNotFoundError(2, "Could not find command/executable", f"{self.command[0]}"))
			raise FileNotFoundError(2, "Could not find command/executable", f"{self.command[0]}")
		else:
			self.command[0] = ex
		p : Popen = Popen(self.command, stdin=stdin, stdout=self.outFile or stdout, stderr=stderr, **kwargs)
		LOGGER.debug(f"Started {self!r}")
		return p

class PipeCommands(Commands):
	"""Takes string dividing commands by "|" and runs the commands simultaneously,
	but piping to each other. Waits for the last command to finish before
	.run() returns."""
	pattern = pipePattern
	nextType = DumpCommands
	_list : list[DumpCommands]

	def run(self, **kwargs) -> CompletedProcess:
		if len(self._list) == 0: return []

		processes : list[Popen] = [self._list[0].run(stdout=PIPE, stderr=PIPE, **kwargs)]
		for i, dc in enumerate(self._list[1:]):
			processes.append( dc.run(stdin=processes[i].stdout, stdout=PIPE, stderr=PIPE, **kwargs))

		processes[-1].wait()
		
		for i, p in enumerate(processes):
			bPrint(f"{p.args!r}\n{p.args!r}[EXITCODE]\n{p.returncode}\n{p.args!r}[STDERR]", file=self.logFile)
			self.logFile.write(p.stderr.read()) if p.stderr.readable() else bPrint("", file=self.logFile)
			if p == processes[-1] and self._list[i].outFile is None:
				bPrint(f"{p.args!r}[STDOUT]", file=self.logFile)
				self.logFile.write(p.stdout.read()) if p.stdout.readable() else bPrint("", file=self.logFile)
			elif p == processes[-1]:
				bPrint(f"{p.args!r}[STDOUT]\n# Output dumped to: {self._list[i].outFile}", file=self.logFile)
			if p.returncode != 0:
				return processes[:i+1]
		return processes

class SequentialCommands(Commands):
	"""Runs commands partitioned by ";" in sequence.
	* If `.start()` is used, a thread is created to run the commands concurrently.
	  The thread triggers a 'SequentialCommands[CATEGORY]Finished' event with the
	  SequentialCommands object itself provided in the eventInfo dict as
	  the 'object' key.
	* If `.run()` is used, the method only returns once all commands have finished.

	Triggers the event `f"SequentialCommands{self.category}Finished"` to its
	SequentialCommands.hooks that is taken from the constructor for
	SequentialCommands when a command has finished, irrespective of using
	`.start()` or `.run()`.
	"""
	pattern = sequentialPattern
	nextType = PipeCommands
	_list : list[PipeCommands]
	processes : list[Popen]
	thread : Thread

	def start(self : SequentialCommands) -> None:
		try:
			self.processes = []
			def runInSequence(self : SequentialCommands):
				try:
					for pc in self._list:
						p = pc.run()
						self.processes.extend(p)
						if self.processes[-1].returncode != 0:
							break
					self.hooks.trigger(f"SequentialCommands{self.category}Finished", {"object" : self})
				except Exception as e:
					if type(e) is FileNotFoundError:
						self.hooks.trigger(f"ReportError", {"exception" : e})
					try:
						e.add_note(f"<In Thread running: {self.raw!r}>")
						LOGGER.exception(e)
					except:
						LOGGER.exception(e)
					self.hooks.trigger(f"SequentialCommands{self.category}Finished", {"object" : self})
			self.thread = Thread(target=runInSequence, args=[self], daemon=True)
			self.thread.start()
		except Exception as e:
			e.add_note(f"{type(self).__name__} failed to `.start()`.")
			LOGGER.exception(e)
			raise e

	def run(self) -> list[CompletedProcess]:
		self.processes = []
		for pc in self._list:
			p = pc.run()
			self.processes.extend(p)
			if self.processes[-1].returncode != 0:
				self.hooks.trigger(f"SequentialCommands{self.category}Finished", {"object" : self})
				return self.processes
		self.hooks.trigger(f"SequentialCommands{self.category}Finished", {"object" : self})
		return self.processes

class ParallelCommands(Commands):
	pattern : re.Pattern = parallelPattern
	nextType = SequentialCommands
	_list : list[SequentialCommands]

	def __init__(self, string, category : str, hooks : Hooks, logFiles : list[str]=None):
		
		try:
			self.raw = string
			self.category = category
			self.hooks = hooks
			_list = [[]]
			
			for c in argsPattern.split(string.strip()):
				if self.pattern.fullmatch(c):
					_list.append([])
				else:
					_list[-1].append(c)
			
			if logFiles is None:
				self.logFiles = [None]*len(_list)
			else:
				self.logFiles = [open(l, "ab") for l in logFiles]
			
			self._list = [self.nextType(l, category, hooks, logFile=logFile) for l, logFile in zip(_list, self.logFiles)]
		except Exception as e:
			e.add_note(f"{type(self).__name__} failed to initialize.")
			LOGGER.exception(e)
			raise e

	def __del__(self):
		for l in self.logFiles:
			try:
				l.close()
			except:
				pass
	
	def start(self) -> list[CompletedProcess]:
		for sc in self._list:
			sc.start()
	
	def run(self) -> list[CompletedProcess]:
		processes = []
		
		for sc in self._list:
			p = sc.run()
			processes.extend(p)
			if processes[-1].returncode != 0:
				return processes
		return processes