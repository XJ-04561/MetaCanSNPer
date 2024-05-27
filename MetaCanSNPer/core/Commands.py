


from subprocess import Popen, PIPE, CompletedProcess
from threading import Thread, Condition

from MetaCanSNPer.core.Hooks import Hooks, GlobalHooks
from MetaCanSNPer.Globals import *

def bPrint(*strings, sep=b" ", end=b"\n", file : BinaryIO=None, encoding : str="utf-8"):
	file.write(sep.join(map(lambda s : s.encode("utf-8"), strings)) + end)

parallelPattern = re.compile(r"\s*[&]\s*")
sequentialPattern = re.compile(r"\s*([;])\s*|\s*([&][&])\s*|\s*([|][|])\s*")
pipePattern = re.compile(r"\s*[|]\s*")
dumpPattern = re.compile(r"\s*[>]\s*")
whitePattern = re.compile(r"\s*")
argsPattern = re.compile(r"(['][^']*?['])|([\"][^\"]*?[\"])|(\S+)", flags=re.MULTILINE+re.DOTALL)
quotePattern = re.compile(r"['\" ]*")
illegalPattern = re.compile(r"[^\w_ \-\.]")

class Command(Logged):
	"""Triggers the event `f"{self.category}Finished"` on each finished parallel command"""

	runningCondition : Condition
	commands : "ParallelCommands"

	def __init__(self, args : str|list, category=None, hooks : Hooks=GlobalHooks, logDir : Path=Path("."), names : Iterable=None):
		try:
			self.raw = args if type(args) is str else " & ".join(args)
			self.category = category
			self.hooks = hooks
			self.names = names
			self.returncodes = {}
			self.exceptions = {}
			self.runningCondition = Condition()
			if self.raw.strip() != "":
				self._hook = self.hooks.addHook(f"{self.category}Finished", target=self.parallelFinished)
				
				self.commands = ParallelCommands(self.raw, category, hooks=hooks, logDir=logDir, names=names)
			else:
				self.commands = None
			
		except Exception as e:
			e.add_note(f"{type(self).__name__} failed to initialize.")
			self.LOG.exception(e)
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

	def parallelFinished(self, eventInfo):
		for name in self.commands:
			if eventInfo["instance"] is self.commands[name]:
				self.returncodes[name] = eventInfo["value"]
				self.hooks.trigger(f"{self.category}Finished", {"name" : name, "instance" : self, "value" : 3})
				break

	def start(self):
		self.LOG.info(f"Starting {self}")
		self.runningCondition.acquire(False)
		if self.commands is not None:
			self.commands.start()

	def run(self):
		self.LOG.info(f"Running {self}")
		self.runningCondition.acquire(False)
		if self.commands is not None:
			self.commands.run()
	
	def wait(self, timeout=None):
		if self.commands is not None:
			self.commands.wait(timeout=timeout)

class Commands(Logged):
	"""Only meant to be inherited from"""

	pattern : re.Pattern
	nextType : "Commands"
	hooks : Hooks
	category : str
	_list : list["Commands"]
	logFile : TextIO
	raw : str

	def __init__(self, args : list[str], category : str, hooks : Hooks=GlobalHooks, logDir : Path=Path(".")):
		try:
			self.raw = "".join(args).strip()
			self.category = category
			self.hooks = hooks
			_list = [[]]
			self.separators = []
			
			for c in args:
				self.LOG.debug(f"{self.pattern}.fullmatch({c}) -> {self.pattern.fullmatch(c)}")
				if c is None:
					continue
				elif self.pattern.fullmatch(c):
					_list.append([])
					self.separators.append(self.pattern.fullmatch(c).group().strip())
				else:
					if c.startswith("'"):
						_list[-1].append(c.strip("'"))
					elif c.startswith("\""):
						_list[-1].append(c.strip("\""))
					else:
						_list[-1].append(c)
			self.separators.append("&&")
			
			self._list = [self.nextType(l, category, hooks, logDir=logDir) for l in _list if len(l) > 0]
		except Exception as e:
			e.add_note(f"{type(self).__name__} failed to initialize.")
			self.LOG.exception(e)
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

	def __init__(self, args, category : str, hooks : Hooks=GlobalHooks, logDir : Path=Path(".")):
		try:
			super().__init__(args, category, hooks, logDir=logDir)
		
			self.command = self._list[0]
			self.LOG.debug(f"{self.command=}")
			if len(self.command) > 1 and self.command[1].isalnum():
				logFile = logDir.writable / f"{self.command[0]}_{self.command[1]}_{SOFTWARE_NAME}.log"
			else:
				logFile = logDir.writable / f"{self.command[0]}_{SOFTWARE_NAME}.log"
			
			try:
				self.logFile = open(logFile, "wb")
			except:
				self.LOG.warning(f"Failed to create {logFile=}")
				self.logFile = DEV_NULL_BYTES
			
			self.outFile = None
			self.outFileName = None

			if len(self._list) == 2:
				if len(self._list[1]) > 1:
					self.LOG.exception(ValueError(f"Output can not be dumped to multiple filenames. Filenames given: {self._list[1]}"))
					raise ValueError(f"Output can not be dumped to multiple filenames. Filenames given: {self._list[1]}")
				
				self.outFileName = self._list[1][0]
				if not os.path.exists(self._list[1][0]):
					self.outFile = open(self.outFileName+".tmp", "wb")
			elif len(self._list) > 2:
				self.LOG.exception(ValueError(f"Output dumped more or less than once using '>' in one command. Command: {'>'.join(map(''.join, self._list))}"))
				raise ValueError(f"Output dumped more or less than once using '>' in one command. Command: {'>'.join(map(''.join, self._list))}")
			self.LOG.debug(f"{self.command=}, {self.outFile=}, {self.logFile=}")
		except Exception as e:
			e.add_note(f"{type(self).__name__} failed to initialize.")
			self.LOG.exception(e)
			raise e
		
	def run(self, stdin=None, stdout=None, stderr=None, **kwargs) -> Popen:
		self.LOG.info(f"Running {self}")
		ex = shutil.which(self.command[0])
		if ex is None:
			self.LOG.exception(FileNotFoundError(2, "Could not find command/executable", f"{self.command[0]}"))
			raise FileNotFoundError(2, "Could not find command/executable", f"{self.command[0]}")
		else:
			self.command[0] = ex
		p : Popen = Popen(self.command, stdin=stdin, stdout=self.outFile or stdout or DEV_NULL_BYTES, stderr=stderr or self.logFile, **kwargs)
		self.LOG.debug(f"Started {self!r}")
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
		self.LOG.info(f"Running {self}")
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
				self.LOG.error(f"Section of pipe closed. Returncodes of commands in pipe: {[p.returncode for p in processes]}\nPipeCommands in question: {self}")

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

	def start(self : "SequentialCommands") -> None:
		self.LOG.info(f"Starting {self}")
		try:
			self.returncodes = []
			def runInSequence(self : SequentialCommands):
				try:
					for i, pc in enumerate(self._list):
						returncode = pc.run()
						self.returncodes.append(returncode)
						if self.separators[i] == ";":
							pass
						elif returncode != 0 and self.separators[i] == "&&":
							break
						elif returncode == 0 and self.separators[i] == "||":
							break
					self.hooks.trigger(f"SequentialCommands{self.category}Finished", {"name" : None, "value" : self.returncodes[-1] if self.returncodes else None, "instance" : self})
				except Exception as e:
					if type(e) is FileNotFoundError:
						self.hooks.trigger(f"ReportError", {"exception" : e})
					try:
						e.add_note(f"<In Thread running: {self.raw!r}>")
						self.LOG.exception(e)
					except:
						self.LOG.exception(e)
					self.hooks.trigger(f"SequentialCommands{self.category}Finished", {"name" : None, "value" : self.returncodes[-1] if self.returncodes else None, "instance" : self})
			self.thread = Thread(target=runInSequence, args=[self], daemon=True)
			self.thread.start()
		except Exception as e:
			e.add_note(f"{type(self).__name__} failed to `.start()`.")
			self.LOG.exception(e)
			raise e

	def run(self) -> list[CompletedProcess]:
		self.LOG.info(f"Running {self}")
		self.processes = []
		for pc in self._list:
			returncode = pc.run()
			self.processes.append(returncode)
			if returncode != 0:
				self.hooks.trigger(f"SequentialCommands{self.category}Finished", {"instance" : self})
				return self.processes
		self.hooks.trigger(f"SequentialCommands{self.category}Finished", {"instance" : self})
		return self.processes
	
	def wait(self, timeout=None):
		
		self.thread.join(timeout=timeout)

class ParallelCommands(Commands):
	pattern : re.Pattern = parallelPattern
	nextType = SequentialCommands
	_list : dict[Any,SequentialCommands]
	_thread : Thread

	def __init__(self, args : str, category : str, hooks : Hooks=GlobalHooks, logDir : Path=Path("."), names : Iterable=None):
		
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
			
			self._list = {name:self.nextType(command, category, hooks, logDir=logDir / illegalPattern.sub("-", str(name))) for name, command in _list.items()}
		except Exception as e:
			e.add_note(f"{type(self).__name__} failed to initialize.")
			self.LOG.exception(e)
			raise e

	def start(self):
		self.LOG.info(f"Starting {self}")
		self._thread = Thread(target=self.run)
		self._thread.start()
		self.LOG.info(f"Started {self}")
		
	
	def run(self):
		self.LOG.info(f"Running {self}")
		for sc in self._list.values():
			sc.start()
		for sc in self._list.values():
			sc.wait()
	
	def wait(self, timeout=None):

		nt = "\n\t"
		self.LOG.info(f"Waiting for ({id(self):X}):\n\t{nt.join(map(format, self._list))}")
		self._thread.join()
		self.LOG.info(f"Done waiting for ({id(self):X})")