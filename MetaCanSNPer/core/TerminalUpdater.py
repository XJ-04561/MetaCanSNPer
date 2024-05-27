
from MetaCanSNPer.core.Hooks import Hooks
from threading import Thread, Condition
from time import sleep
from timeit import default_timer as timer
from typing import TextIO, Iterable
from functools import cached_property, partial
from MetaCanSNPer.Globals import *
from colors import black, red, green, yellow, blue, magenta, cyan, white
import textwrap

_NOT_FOUND = object()

# print(sys.stdout.encoding)
# if sys.stdout.encoding == "utf-8":
# 	SQUARE = "\u2588"
# 	HALF_SQUARE = "\u258C"
# else:

SQUARE = "="
HALF_SQUARE = ":"

class Printer:

	LOCK : Lock = Lock()
	ANSI_MATCH = re.compile("\u001b.*?m|\x1b.*?m")
	def ANSI_REMOVE(self, string):
		return self.ANSI_MATCH.sub("", string)

	def __init__(self, out=sys.stdout):
		self.out = out
		self.last = 0
	def __call__(self, msg):
		with self.LOCK:
			print("\b"*self.last, end=msg, flush=True, file=self.out)
			self.last = len(self.ANSI_REMOVE(msg))
	def clear(self):
		with self.LOCK:
			print("\b"*self.last+" "*self.last, end="\b"*self.last, flush=True, file=self.out)
			self.last = 0

class HitchableDict(dict):

	onSet : Callable = None

	def __new__(cls, *args, **kwargs):
		obj = super().__new__(cls, *args, **kwargs)
		return obj

	def __init__(self, data, onSet=None, **kwargs):
		super().__init__(data, **kwargs)
		self.onSet = onSet

	def __setitem__(self, key, value):
		super().__setitem__(key, value)
		if self.onSet is not None:
			self.onSet()

# Taken from django @ https://github.com/django/django/blob/main/django/core/management/color.py
@cache
def supportsColor():
	"""
	From django @ https://github.com/django/django/blob/main/django/core/management/color.py
	Return True if the running system's terminal supports color,
	and False otherwise.
	"""
	import sys, os
	try:
		import colorama # type: ignore

		# Avoid initializing colorama in non-Windows platforms.
		colorama.just_fix_windows_console()
	except (
		AttributeError,  # colorama <= 0.4.6.
		ImportError,  # colorama is not installed.
		# If just_fix_windows_console() accesses sys.stdout with
		# WSGIRestrictedStdout.
		OSError,
	):
		HAS_COLORAMA = False
	else:
		HAS_COLORAMA = True

	def vt_codes_enabled_in_windows_registry():
		"""
		Check the Windows Registry to see if VT code handling has been enabled
		by default, see https://superuser.com/a/1300251/447564.
		"""
		try:
			# winreg is only available on Windows.
			import winreg
		except ImportError:
			return False
		else:
			try:
				reg_key = winreg.OpenKey(winreg.HKEY_CURRENT_USER, "Console")
				reg_key_value, _ = winreg.QueryValueEx(reg_key, "VirtualTerminalLevel")
			except FileNotFoundError:
				return False
			else:
				return reg_key_value == 1

	# isatty is not always implemented, #6223.
	is_a_tty = hasattr(sys.stdout, "isatty") and sys.stdout.isatty()

	return is_a_tty and (
		sys.platform != "win32"
		or (HAS_COLORAMA and getattr(colorama, "fixed_windows_console", False))
		or "ANSICON" in os.environ
		or
		# Windows Terminal supports VT codes.
		"WT_SESSION" in os.environ
		or
		# Microsoft Visual Studio Code's built-in terminal supports colors.
		os.environ.get("TERM_PROGRAM") == "vscode"
		or vt_codes_enabled_in_windows_registry()
	)

class Indicator(Logged):
	
	running : bool = False

	_symbols : list[str]
	sep : str
	_borders : tuple[str,str]
	crashSymbol : str
	finishSymbol : str
	_threads : dict[str,float]
	condition : Condition
	finishedThreads : set
	names : list[str]
	rowTemplate : str
	rowLock : Lock

	backspaces : str
	whitespaces : str

	sepLength : int
	borderLength : int
	innerLength : int
	length : int
	
	@property
	def terminalWidth(self):
		if ISATTY:
			return os.get_terminal_size()[0]
		else:
			return 80
	
	@property
	def threads(self) -> HitchableDict[str,float]:
		return self._threads
	
	@threads.setter
	def threads(self, threads : HitchableDict[str,float]):
		with self.rowLock:
			self._threads = threads
			self.names = sorted(threads.keys())
			self.shortKeys = [name if len(name) < self.length else name[:self.length-2]+"..." for name in self.names]
			self.N = len(self.names)
			self.entries = list(map(lambda _:" "*self.innerLength, range(self.N)))
			self.finishedThreads.intersection_update(self._threads)

	def __init__(self, threads : HitchableDict, symbols : tuple[str], length : int=15, message : str="", sep : str=" ", borders : tuple[str,str]=("[", "]"), crashSymbol="X", finishSymbol=SQUARE, out=sys.stdout, preColor : str=None, partition : str=None, crashColor : str=None, skippedColor : str=None, finishColor : str=None, postColor : str=None, progColor : str=None):
		
		self.condition = Condition()
		self.rowLock = Lock()

		self.preColor		= partial(yellow, bg="black")
		self.progColor		= partial(magenta, bg="black")
		self.partition		= partial(white, bg="black")
		self.crashColor		= partial(red, bg="black")
		self.skippedColor	= partial(red, bg="blue")
		self.finishColor	= partial(green, bg="black")
		self.postColor		= partial(cyan, bg="black")

		self.n = 0

		self.message = message
		self.out = out

		self.symbols		= symbols
		self.sep			= sep
		self.borders		= borders
		self.crashSymbol	= crashSymbol
		self.finishSymbol	= finishSymbol

		self.length = length

		self.finishedThreads = set()
		self.threads = threads
	
	@property
	def symbols(self) -> tuple[str]:
		return self._symbols
	@symbols.setter
	def symbols(self, value):
		self._symbols = tuple(value)
		
	@property
	def borders(self) -> tuple[str,str]:
		return self._borders
	@borders.setter
	def borders(self, value):
		self._borders = tuple(value)

	@Default["length", "borderLength"]
	def innerLength(self):
		return self.length - self.borderLength

	@Default["sep"]
	def sepLength(self):
		return len(self.sep)

	@Default["borders"]
	def borderLength(self):
		return len(self.borders[0]) + len(self.borders[1])

	@Default["borderLength", "sepLength"]
	def outerLength(self):
		return self.borderLength + self.sepLength
	
	@cache
	def createRowTemplate(self, width : int, N : int) -> tuple[str, str, str]:
		"""createRowTemplate(self, width : int) -> backspaces, whitespaces, rowTemplate
		"""
		firstRow = f"{self.message}: {{time:<{width-2-len(self.message)}}}"
		spacerRow = " " * width
		
		maxCols = (width+self.sepLength-2) // (self.length+self.sepLength)
		
		namesList = []
		for cols in itertools.batched(map(lambda i: f"{{names[{i}]}}", range(N)), maxCols):
			namesList.append( " "+self.sep.join(cols).ljust(width-1))
		barsList = []
		for cols in itertools.batched(map(lambda i: f"{self.borders[0]}{{bars[{i}]}}{self.borders[1]}", range(N)), maxCols):
			barsList.append( " "+self.sep.join(cols).ljust(width-1))
		
		rowTemplate = []
		for namesRow, barsRow in zip(namesList, barsList):
			rowTemplate.append(namesRow)
			rowTemplate.append(barsRow)
			rowTemplate.append(spacerRow)
		
		return firstRow + "".join(rowTemplate)

	@property
	def rowTemplate(self) -> str:
		return self.createRowTemplate(self.terminalWidth, self.N)

	@property
	def rowGenerator(self) -> Generator[tuple[int,str],None,None]:
		"""Progress special cases:
		None	-	Service crashed
		2		-	Never ran/skipped
		3		-	Completed
		0<->1	-	Running
		<0		-	Not started
		1.0		-	Postprocessing
		"""
		raise NotImplementedError("rowGenerator not implemented in the base class")

	def run(self):
		
		startTime = timer()
		self.running = True
		
		flushPrint = Printer(self.out)
		
		while self.running:
			with self.rowLock:
				flushPrint(self.rowTemplate.format(time=formatTimestamp(timer()-startTime), names=self.shortKeys, bars=tuple(self.rowGenerator)))
			sleep(0.1)
			self.condition.acquire(timeout=0.40)
		
		with self.rowLock:
			if None in self.threads.values():
				flushPrint(self.rowTemplate.format(time=f"{red('Failed!')} {formatTimestamp(timer()-startTime)}", names=self.shortKeys, bars=tuple(self.rowGenerator)))
			elif self.finishedThreads.issuperset(self.threads):
				flushPrint(flushPrint(f"{self.message} {green('Done!')} {formatTimestamp(timer()-startTime)}"))
			else:
				flushPrint(self.rowTemplate.format(time=f"{yellow('Interrupted!')} {formatTimestamp(timer()-startTime)}", names=self.shortKeys, bars=tuple(self.rowGenerator)))
	
	def checkDone(self):
		if all(v is None or v == 2 or v == 3 for v in self.threads.values()):
			self.running = False

	def kill(self):
		self.running = False
		try:
			self.condition.notify_all()
		except RuntimeError:
			pass

class LoadingBar(Indicator):

	def __init__(self, threads, length=15, fill=SQUARE, halfFill=HALF_SQUARE, background=" ", **kwargs):
		super().__init__(threads, [fill, halfFill, background], length=length, **kwargs)
		self.innerLength = length - self.borderLength

	@property
	def rowGenerator(self) -> Generator[str,None,None]:
		"""Progress special cases:
		None	-	Service crashed
		2		-	Never ran/skipped
		3		-	Completed
		0<->1	-	Running
		<0		-	Not started
		1.0		-	Postprocessing
		"""
		crashString = self.crashColor(self.crashSymbol * self.innerLength)
		finishString = self.finishColor(self.symbols[0] * self.innerLength)
		skippedString = self.skippedColor(self.symbols[0] * self.innerLength)
		
		for name, prog in zip(self.names, map(self.threads.get, self.names)):
			if prog is None:
				yield crashString
			elif prog == 2:
				yield skippedString
			elif name in self.finishedThreads:
				yield finishString
			elif 0 <= prog < 1.0:
				fillLength = int(self.innerLength*2*prog)
				fillLength, halfBlock = fillLength//2, fillLength%2
				
				yield self.progColor(f"{self.symbols[0]*fillLength}{self.symbols[1]*halfBlock}")+self.partition(f"{self.symbols[2]*(self.innerLength - fillLength - halfBlock)}")
			elif prog < 0:
				yield self.preColor(self.symbols[2]*self.n + self.symbols[0] + self.symbols[2]*(self.innerLength - self.n - 1))
			elif prog == 1:
				yield self.postColor(self.symbols[0]*self.n + self.symbols[2] + self.symbols[0]*(self.innerLength - self.n - 1))
			else:
				yield crashString
		self.n = (self.n+1) % self.innerLength

class Spinner(Indicator):

	symbols : list[str] = ["|", "/", "-", "\\"]

	@overload
	def __init__(self, threads: HitchableDict, symbols: tuple[str]=["|", "/", "-", "\\"], length: int = 10,
			  message: str = "", sep: str = " ", borders: tuple[str, str] = ("[", "]"), crashSymbol="X",
			  finishSymbol=SQUARE, out=sys.stdout, preColor: str = None, partition: str = None,
			  crashColor: str = None, skippedColor: str = None, finishColor: str = None,
			  postColor: str = None, progColor: str = None): ...
	def __init__(self, *args, **kwargs):
		if len(args) < 2 and "symbols" not in kwargs:
			args = (*args, self.symbols)
		super().__init__(*args, **kwargs)

	@property
	def rowGenerator(self) -> Generator[tuple[int,str],None,None]:
		"""Progress special cases:
		None	-	Service crashed
		2		-	Never ran/skipped
		3		-	Completed
		0<->1	-	Running
		<0		-	Not started
		1.0		-	Postprocessing
		"""
		NSymbols = len(self.symbols)
		crashString = self.crashColor(self.crashSymbol)
		finishString = self.finishColor(self.finishSymbol)
		skippedString = self.skippedColor(self.finishSymbol)
		
		for name, prog in zip(self.names, map(self.threads.get, self.names)):
			if prog is None:
				yield crashString
			elif prog == 2:
				yield skippedString
			elif name in self.finishedThreads:
				yield finishString
			elif 0 <= prog < 1:
				yield self.progColor(self.symbols[self.n])
			elif prog < 0:
				yield self.preColor("." if self.n%2 else ":")
			elif prog == 1:
				yield self.postColor("." if self.n%2 else ":")
			else:
				yield crashString
		self.n = (self.n+1) % NSymbols

class TextProgress(Indicator):
	
	symbols : tuple[str]=[".", ",", ":", "|", "I", "H", "#"]

	@overload
	def __init__(self, threads: HitchableDict, symbols: tuple[str]=[".", ",", ":", "|", "I", "H", "#"], length: int = 3,
			  message: str = "", sep: str = " ", borders: tuple[str, str] = ("[", "]"), crashSymbol="X",
			  finishSymbol=SQUARE, out=sys.stdout, preColor: str = None, partition: str = None,
			  crashColor: str = None, skippedColor: str = None, finishColor: str = None,
			  postColor: str = None, progColor: str = None): ...
	def __init__(self, *args, **kwargs):
		if len(args) < 2 and "symbols" not in kwargs:
			args = (*args, self.symbols)
		super().__init__(*args, **kwargs)

	@property
	def rowGenerator(self) -> Generator[tuple[int,str],None,None]:
		"""Progress special cases:
		None	-	Service crashed
		2		-	Never ran/skipped
		3		-	Completed
		0<->1	-	Running
		<0		-	Not started
		1.0		-	Postprocessing
		"""
		NSymbols = len(self.symbols)
		crashString = self.crashColor(self.crashSymbol)
		finishString = self.finishColor(self.finishSymbol)
		skippedString = self.skippedColor(self.finishSymbol)
		
		for name, prog in zip(self.names, map(self.threads.get, self.names)):
			if prog is None:
				yield crashString
			elif prog == 2:
				yield skippedString
			elif name in self.finishedThreads:
				yield finishString
			elif 0 <= prog < 1:
				yield self.progColor(self.symbols[int(NSymbols*prog)])
			elif prog < 0:
				yield self.preColor("."*self.n + ":"*(self.n<self.innerLength) + "."*(self.innerLength-self.n))
			elif prog == 1:
				yield self.postColor("o"*self.n + "O"*(self.n<self.innerLength) + "o"*(self.innerLength-self.n))
			else:
				yield crashString
		self.n = (self.n+1) % (self.innerLength + 1)

class TerminalUpdater(Logged):
	
	threads : HitchableDict
	thread : Thread
	hooks : Hooks
	out : TextIO
	printer : Indicator = None

	def __init__(self, message, category, hooks : Hooks, threadNames : Iterable, out=sys.stdout, printer : Indicator=Spinner):
		
		self.startTime = timer()
		
		self.message = message
		self.category = category
		self.hooks = hooks
		self.threadNames = threadNames
		self.threads = HitchableDict(map(lambda name : (name,-1.0), threadNames), onSet=self.updatePrinter)
		self.out = out
		
		self.running = True

		"""Progress special cases:
		None	-	Service crashed
		2		-	Never ran/skipped
		3		-	Completed
		0<->1	-	Running
		<0		-	Not started
		1.0		-	Postprocessing
		"""
		
		# self.hooks.addHook(f"{self.category}Initialized", self.initializedCallback) # Not needed, is already initialized by default.
		self.hooks.addHook(f"{self.category}Skipped", self.skippedCallback)
		self.hooks.addHook(f"{self.category}Starting", self.startingCallback)
		self.hooks.addHook(f"{self.category}Progress", self.progressCallback)
		self.hooks.addHook(f"{self.category}PostProcess", self.postProcessCallback)
		self.hooks.addHook(f"{self.category}Finished", self.finishedCallback)
		self.hooks.addHook(f"{self.category}Failed", self.failedCallback)

		if printer is not None:
			self.setPrinter(printer)

	def __enter__(self):
		self.start()
		return self

	def __exit__(self, *args):
		self.stop()

	def initializedCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["name"] in self.threads:
			self.threads[eventInfo["name"]] = -1.0

	def skippedCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["name"] in self.threads:
			self.threads[eventInfo["name"]] = 2
			self.printer.finishedThreads.add(eventInfo["name"])

	def startingCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["name"] in self.threads:
			self.threads[eventInfo["name"]] = 0.0

	def progressCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["name"] in self.threads:
			self.threads[eventInfo["name"]] = eventInfo["value"]

	def postProcessCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["name"] in self.threads:
			self.threads[eventInfo["name"]] = 1.0

	def finishedCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["name"] in self.threads:
			self.threads[eventInfo["name"]] = 3
			self.printer.finishedThreads.add(eventInfo["name"])

	def failedCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["name"] in self.threads:
			self.threads[eventInfo["name"]] = None
	
	def setPrinter(self, printer : Indicator=Spinner, *args, **kwargs):
		self.printer = printer(self.threads, *args, **({"out":self.out, "message":self.message} | kwargs))

	def start(self, *args, **kwargs):
		self.thread = Thread(target=self.mainLoop, args=args, kwargs=kwargs, daemon=True)
		self.thread.start()
	
	def wait(self, timeout=3):
		self.thread.join(timeout=timeout)
	
	def kill(self):
		"""Does not wait for thread to stop."""
		self.printer.running = False
		try:
			self.printer.condition.notify_all()
		except RuntimeError:
			pass

	def stop(self):
		"""Waits for thread to stop."""
		self.printer.running = False
		try:
			self.printer.condition.notify_all()
		except RuntimeError:
			pass
		self.thread.join()

	def mainLoop(self, *args, **kwargs):
		try:
			self.printer.run(*args, **kwargs)
		except Exception as e:
			self.LOG.exception(e)

	def updatePrinter(self):
		try:
			self.printer.threads = self.threads
		except:
			pass