



import re, shutil, sys
from subprocess import Popen, PIPE, CompletedProcess
from typing import TextIO, BinaryIO, Iterable
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
whitePattern = re.compile("\s*")
argsPattern = re.compile(r"(['][^']*?['])|([\"][^\"]*?[\"])|(\S+)", flags=re.MULTILINE+re.DOTALL)
quotePattern = re.compile(r"['\" ]*")

class ParallelCommands: pass
class SequentialCommands: pass
class PipeCommands: pass
class DumpCommands: pass
class Commands: pass


class Command:
	"""Triggers the event `f"{self.category}ProcessFinished"` on each finished parallel command"""


	commands : ParallelCommands

	def __init__(self, args : str|list, category=None, hooks : Hooks=None, logFiles : list[str]|str=None, names : Iterable=None):
		try:
			self.raw = args if type(args) is str else " & ".join(args)
			self.logFiles = [logFiles] if type(logFiles) is str else logFiles
			self.category = category
			self.hooks = hooks
			self.returncodes = {}
			self.exceptions = {}
			
			def parallelFinished(eventInfo, self : Command):
				try:
					for name, com in self.commands:
						if eventInfo["object"] is com:
							self.returncodes[name] = None if len(eventInfo["object"].processes) == 0 else eventInfo["object"].processes[-1].returncode
							self.hooks.trigger(f"{self.category}ProcessFinished", {"threadN" : name, "Command" : self} | eventInfo)
							return
				except Exception as e:
					e.add_note("'parallelFinished' event exception.")
					LOGGER.exception(e)
					return

			self._hook = self.hooks.addHook(f"SequentialCommands{self.category}Finished", target=parallelFinished, args=[self])

			self.commands = ParallelCommands(self.raw, category, hooks, logFiles=self.logFiles, names=names)
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
				LOGGER.debug(f"{self.pattern}.fullmatch({c}) -> {self.pattern.fullmatch(c)}")
				if c is None:
					continue
				elif c.strip() == "":
					continue
				elif self.pattern.fullmatch(c):
					_list.append([])
				else:
					if c.startswith("'"):
						_list[-1].append(c.strip("'"))
					elif c.startswith("\""):
						_list[-1].append(c.strip("\""))
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
			LOGGER.debug(f"{command=}, {self.command=}")

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
			LOGGER.debug(f"{command=}, {self.command=}, {self.outFile=}")
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
			if p.wait(0.5) is None: # .wait returns the returncode. Wait a few millisec so that prior processes don't close pipe before they terminate.
				LOGGER.error(f"Section of pipe closed. Returncodes of commands in pipe: {[p.returncode for p in processes]}\nPipeCommands in question: {self}")
				break

		for i, p in enumerate(processes):
			bPrint(f"{p.args!r}\n{p.args!r}[EXITCODE]\n{p.returncode}\n{p.args!r}[STDERR]", file=self.logFile)
			self.logFile.write(p.stderr.read()) if p.stderr.readable() else bPrint("", file=self.logFile)
			if p == processes[-1] and self._list[i].outFile is None:
				bPrint(f"{p.args!r}[STDOUT]", file=self.logFile)
				self.logFile.write(p.stdout.read()) if p.stdout.readable() else bPrint("", file=self.logFile)
			elif p == processes[-1]:
				bPrint(f"{p.args!r}[STDOUT]\n# Output dumped to: {self._list[i].outFile}", file=self.logFile)
		for i, p in enumerate(processes):
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

	def __init__(self, args : str, category : str, hooks : Hooks, logFiles : list[str]=None, names : Iterable=None):
		
		try:
			self.category = category
			self.hooks = hooks
			self.logFiles = []
			if names is None:
				def _names():
					n = 0
					while True:
						yield n
						n += 1
				names = _names()
			else:
				names = iter(names)
			if logFiles is None:
				def _logFiles():
					while True:
						yield None
				logFiles = _logFiles()
			else:
				logFiles = iter(map(lambda f: open(f, "ab"), logFiles))
			
			logFile = next(logFiles)
			name = next(names)
			_list = {name:[]}
			self._logFiles = {name:logFile}
			self.raw = args
			for c in argsPattern.split(self.raw.strip()):
				if c is None:
					continue
				elif c.strip() == "":
					continue
				elif self.pattern.fullmatch(c):
					name = next(names)
					_list[name] = []
					_logFiles[name] = next(logFiles)
				else:
					if c.startswith("'"):
						_list[name].append(c.strip("'"))
					elif c.startswith("\""):
						_list[name].append(c.strip("\""))
					else:
						_list[name].append(c)
			
			self._list = {name:self.nextType(command, category, hooks, logFile=self._logFiles[name]) for name, command in _list.items()}
		except Exception as e:
			e.add_note(f"{type(self).__name__} failed to initialize.")
			LOGGER.exception(e)
			raise e

	def __del__(self):
		for l in self.logFiles.values():
			try:
				l.close()
			except:
				pass
	
	def start(self) -> list[CompletedProcess]:
		for sc in self._list.values():
			sc.start()
	
	def run(self) -> list[CompletedProcess]:
		processes = []
		
		for sc in self._list.values():
			p = sc.run()
			processes.extend(p)
			if processes[-1].returncode != 0:
				return processes
		return processes