



import re, shutil, sys, os
from subprocess import Popen, PIPE, CompletedProcess
from typing import TextIO, BinaryIO, Iterable, Any
from threading import Thread
from PseudoPathy import Path

from MetaCanSNPer.modules.LogKeeper import createLogger

LOGGER = createLogger(__name__)

from MetaCanSNPer.modules.Hooks import Hooks
from MetaCanSNPer.Globals import SOFTWARE_NAME

def bPrint(*strings, sep=b" ", end=b"\n", file : BinaryIO=None, encoding : str="utf-8"):
	file.write(sep.join(map(lambda s : s.encode("utf-8"), strings)) + end)

parallelPattern = re.compile("\s*[&]\s*")
sequentialPattern = re.compile("\s*[;]\s*")
pipePattern = re.compile("\s*[|]\s*")
dumpPattern = re.compile("\s*[>]\s*")
whitePattern = re.compile("\s*")
argsPattern = re.compile(r"(['][^']*?['])|([\"][^\"]*?[\"])|(\S+)", flags=re.MULTILINE+re.DOTALL)
quotePattern = re.compile(r"['\" ]*")
illegalPattern = re.compile(r"[^\w_ \-\.]")

class ParallelCommands: pass
class SequentialCommands: pass
class PipeCommands: pass
class DumpCommands: pass
class Commands: pass


class Command:
	"""Triggers the event `f"{self.category}ProcessFinished"` on each finished parallel command"""


	commands : ParallelCommands

	def __init__(self, args : str|list, category=None, hooks : Hooks=None, logDir : Path=Path("."), names : Iterable=None):
		try:
			self.raw = args if type(args) is str else " & ".join(args)
			self.category = category
			self.hooks = hooks
			self.names = names
			self.returncodes = {}
			self.exceptions = {}
			if self.raw.strip() != "":
				def parallelFinished(eventInfo, self : Command):
					try:
						for name in self.commands:
							if eventInfo["object"] is self.commands[name]:
								self.returncodes[name] = None if len(eventInfo["object"].returncodes) == 0 else eventInfo["object"].returncodes[-1]
								self.hooks.trigger(f"{self.category}ProcessFinished", {"threadN" : name, "Command" : self})
								return
					except Exception as e:
						e.add_note("'parallelFinished' event exception.")
						LOGGER.exception(e)
						return

				self._hook = self.hooks.addHook(f"SequentialCommands{self.category}Finished", target=parallelFinished, args=[self])
				
				self.commands = ParallelCommands(self.raw, category, hooks, logDir=logDir, names=names)
			else:
				self.commands = None
			
		except Exception as e:
			e.add_note(f"{type(self).__name__} failed to initialize.")
			LOGGER.exception(e)
			raise e
	
	def __len__(self):
		return len(self.commands._list) if self.commands is not None else 0
	
	def __iter__(self):
		return iter(self.commands) if self.commands is not None else iter([])
	
	def __getitem__(self, key):
		if self.commands is not None:
			return self.commands[key]
		else:
			raise KeyError("Key {key!r} not associated with a command.")
	
	def __repr__(self):
		return f"<{type(self).__name__} {hex(id(self))} raw={self.raw!r}>"
	
	def __del__(self):
		try:
			self.hooks.removeHook(f"SequentialCommands{self.category}Finished", self._hook)
		except:
			pass
	
	def start(self):
		if self.commands is not None:
			self.commands.start()

	def run(self):
		if self.commands is not None:
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

	def __init__(self, args : list[str], category : str, hooks : Hooks, logDir : Path=Path(".")):
		try:
			self.raw = "".join(args).strip()
			self.category = category
			self.hooks = hooks
			_list = [[]]
			
			for c in args:
				LOGGER.debug(f"{self.pattern}.fullmatch({c}) -> {self.pattern.fullmatch(c)}")
				if c is None:
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
			
			self._list = [self.nextType(l, category, hooks, logDir=logDir) for l in _list if len(l) > 0]
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
	outFileName : str
	_list : list[str]

	@staticmethod
	def nextType(args, *overFlowArgs, **kwargs) -> list[str]:
		return list(filter(lambda s : whitePattern.fullmatch(s) is None, args))

	def __init__(self, args, category : str, hooks : Hooks, logDir : Path=Path(".")):
		try:
			super().__init__(args, category, hooks, logDir=logDir)
		
			self.command = self._list[0]
			LOGGER.debug(f"{self.command=}")
			if len(self.command) > 1 and self.command[1].isalnum():
				logFile = logDir.writable > f"{self.command[0]}_{self.command[1]}_{SOFTWARE_NAME}.log"
			else:
				logFile = logDir.writable > f"{self.command[0]}_{SOFTWARE_NAME}.log"
			
			try:
				self.logFile = open(logFile, "wb")
			except:
				LOGGER.warning(f"Failed to create {logFile=}")
				self.logFile = open(os.devnull, "wb")
			
			self.outFile = None
			self.outFileName = None

			if len(self._list) == 2:
				if len(self._list[1]) != 1:
					LOGGER.exception(ValueError(f"Output can not be dumped to multiple filenames. Filenames given: {self._list[1]}"))
					raise ValueError(f"Output can not be dumped to multiple filenames. Filenames given: {self._list[1]}")
				
				self.outFileName = self._list[1][0]
				if not os.path.exists(self._list[1][0]):
					self.outFile = open(self.outFileName+".tmp", "wb")
			elif len(self._list) > 2:
				LOGGER.exception(ValueError(f"Output dumped more or less than once using '>' in one command. Command: {'>'.join(map(''.join, self._list))}"))
				raise ValueError(f"Output dumped more or less than once using '>' in one command. Command: {'>'.join(map(''.join, self._list))}")
			LOGGER.debug(f"{self.command=}, {self.outFile=}, {self.logFile=}")
		except Exception as e:
			e.add_note(f"{type(self).__name__} failed to initialize.")
			LOGGER.exception(e)
			raise e
		
	def run(self, stdin=None, stdout=None, stderr=None, **kwargs) -> Popen:
		LOGGER.debug(f"Starting {self!r}")
		ex = shutil.which(self.command[0])
		if ex is None:
			LOGGER.exception(FileNotFoundError(2, "Could not find command/executable", f"{self.command[0]}"))
			raise FileNotFoundError(2, "Could not find command/executable", f"{self.command[0]}")
		else:
			self.command[0] = ex
		p : Popen = Popen(self.command, stdin=stdin, stdout=self.outFile or stdout or open(os.devnull, "wb"), stderr=stderr or self.logFile, **kwargs)
		LOGGER.debug(f"Started {self!r}")
		return p

	def __del__(self):
		try:
			self.outFile.close()
			os.rename(self.outFile.name, self.outFileName)
		except:
			pass
		try:
			self.logFile.close()
		except:
			pass

class PipeCommands(Commands):
	"""Takes string dividing commands by "|" and runs the commands simultaneously,
	but piping to each other. Waits for the last command to finish before
	.run() returns."""
	pattern = pipePattern
	nextType = DumpCommands
	_list : list[DumpCommands]

	def run(self, **kwargs) -> CompletedProcess:
		if len(self._list) == 0: return []

		processes : list[Popen] = []
		lastSTDOUT = None
		for i, dc in enumerate(self._list[:-1]):
			processes.append( dc.run(stdin=lastSTDOUT, stdout=PIPE, **kwargs))
			try:
				lastSTDOUT.close()
			except:
				pass
			lastSTDOUT = processes[i].stdout
		processes.append( self._list[-1].run(stdin=lastSTDOUT, **kwargs))
		try:
			lastSTDOUT.close()
		except:
			pass

		processes[-1].wait()

		for p in processes[::-1]:
			if p.returncode is None:
				LOGGER.error(f"Section of pipe closed. Returncodes of commands in pipe: {[p.returncode for p in processes]}\nPipeCommands in question: {self}")

		for i, p in enumerate(processes):
			if p.returncode not in [0, None]:
				return p.returncode
		try:
			self[-1].outFile.close()
			os.rename(self[-1].outFile.name, self[-1].outFileName)
		except:
			pass
		return processes[-1].returncode

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
	returncodes : list[int]
	thread : Thread

	def start(self : SequentialCommands) -> None:
		try:
			self.returncodes = []
			def runInSequence(self : SequentialCommands):
				try:
					for pc in self._list:
						returncode = pc.run()
						self.returncodes.append(returncode)
						if returncode != 0:
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
			returncode = pc.run()
			self.processes.append(returncode)
			if returncode != 0:
				self.hooks.trigger(f"SequentialCommands{self.category}Finished", {"object" : self})
				return self.processes
		self.hooks.trigger(f"SequentialCommands{self.category}Finished", {"object" : self})
		return self.processes

class ParallelCommands(Commands):
	pattern : re.Pattern = parallelPattern
	nextType = SequentialCommands
	_list : dict[Any,SequentialCommands]

	def __init__(self, args : str, category : str, hooks : Hooks, logDir : Path=Path("."), names : Iterable=None):
		
		try:
			self.category = category
			self.hooks = hooks
			if names is None:
				def _names():
					n = 0
					while True:
						yield n
						n += 1
				names = _names()
			else:
				names = iter(names)
			
			name = next(names)
			_list = {name:[]}
			self.raw = args
			for c in argsPattern.split(self.raw.strip()):
				if c is None:
					continue
				elif self.pattern.fullmatch(c):
					name = next(names)
					_list[name] = []
				else:
					if c.startswith("'"):
						_list[name].append(c.strip("'"))
					elif c.startswith("\""):
						_list[name].append(c.strip("\""))
					else:
						_list[name].append(c)
			
			self._list = {name:self.nextType(command, category, hooks, logDir=logDir > illegalPattern.sub("-", str(name))) for name, command in _list.items()}
		except Exception as e:
			e.add_note(f"{type(self).__name__} failed to initialize.")
			LOGGER.exception(e)
			raise e
	
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