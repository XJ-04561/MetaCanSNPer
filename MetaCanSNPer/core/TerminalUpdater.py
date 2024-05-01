
from MetaCanSNPer.core.Hooks import Hooks
from threading import Thread, Condition
from time import sleep
from timeit import default_timer as timer
from sys import stdout
from typing import TextIO, Iterable
from functools import cached_property
from MetaCanSNPer.Globals import *
from MetaCanSNPer.core.LogKeeper import createLogger
import textwrap

LOGGER = createLogger(__name__)
_NOT_FOUND = object()

class HashCachedProperty:
	__cache : dict
	attrname : str = None
	def __init__(self, *watch):
		self.__cache = {}
		self.watch = watch

	def __set_name__(self, owner, name):
		if self.attrname is None:
			self.attrname = name
		elif name != self.attrname:
			raise TypeError(
				"Cannot assign the same cached_property to two different names "
				f"({self.attrname!r} and {name!r})."
			)

	def __call__(self, func):
		self.func = func
		return self
	
	def __get__(self, instance, owner=None):
		try:
			if (key := tuple(map(instance.__getattribute__, self.watch))) in self._cache:
				return self._cache[key]
			else:
				self.__cache[key] = out = self.func(instance)
				return out
		except TypeError: # Unhashable
			LOGGER.error("Attempted to fetch cached value for property with \
				dependencies that were in part or whole unhashable. Dependencies: \
				" + ", ".join([f"{instance.__name__}.{name}={instance.__getattribute__(name)!r}" for name in self.watch]))
			return self.func


class HitchableDict(dict):

	onSet : Callable = None

	def __new__(cls, *args, **kwargs):
		obj = super().__new__(cls, *args, **kwargs)
		return obj

	def __init__(self, data, onSet=None, **kwargs):
		super().__init__(data, **kwargs)
		self.onSet = onSet

	def __setitem__(self, key, value):
		dict.__setitem__(key, value)
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

class Indicator:
	
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

	backspaces : str
	whitespaces : str

	sepLength : int
	borderLength : int
	innerLength : int
	length : int
	
	@property
	def threads(self) -> HitchableDict[str,float]:
		return self._threads
	
	@threads.setter
	def threads(self, threads : HitchableDict[str,float]):
		self._threads = threads
		self.names = sorted(threads.keys())
		self.shortKeys = [textwrap.shorten(name, self.length, placeholder="..") for name in self.names]
		self.N = len(self.names)
		self.entries = list(map(lambda _:" "*self.innerLength, range(self.N)))
		self.createRowTemplate(os.get_terminal_size()[0])
		self.finishedThreads.intersection_update(self._threads)

	def __init__(self, threads : HitchableDict, symbols : tuple[str], length : int=1, message : str="", sep : str=" ", borders : tuple[str,str]=("[", "]"), crashSymbol=None, finishSymbol="\u2588", out=stdout, preColor : str=None, partition : str=None, crashColor : str=None, skippedColor : str=None, finishColor : str=None, postColor : str=None, progColor : str=None):
		
		self.preColor		= preColor or ("\u001b[33;40m" if supportsColor() else "")
		self.progColor		= progColor or ("\u001b[35;40m" if supportsColor() else "")
		self.partition		= partition or ("\u001b[37;40m" if supportsColor() else "")
		self.crashColor		= crashColor or ("\u001b[31;40m" if supportsColor() else "")
		self.skippedColor	= skippedColor or ("\u001b[31;44m" if supportsColor() else "")
		self.finishColor	= finishColor or ("\u001b[32;40m" if supportsColor() else "")
		self.postColor		= postColor or ("\u001b[36;40m" if supportsColor() else "")

		self.message = message
		self.condition = Condition()
		self.out = out

		self.symbols		= symbols
		self.sep			= sep
		self.borders		= borders
		self.crashSymbol	= crashSymbol
		self.finishSymbol	= finishSymbol

		self.length = length

		self.threads = threads
		self.finishedThreads = set()
	
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

	@HashCachedProperty("length", "borderLength")
	def innerLength(self):
		return self.length - self.borderLength

	@HashCachedProperty("sep")
	def sepLength(self):
		return len(self.sep)

	@HashCachedProperty("borders")
	def borderLength(self):
		return len(self.borders[0]) + len(self.borders[1])

	@HashCachedProperty("borderLength", "sepLength")
	def outerLength(self):
		return self.borderLength + self.sepLength
	
	@cache
	def createRowTemplate(self, width : int) -> tuple[str, str, str]:
		"""createRowTemplate(self, width : int) -> backspaces, whitespaces, rowTemplate
		"""
		firstRow = self.message+" {time:<"+str(width-1-len(self.message))+"}"
		entriesPerRow = max(width-self.sepLength, self.length) // (self.length+self.sepLength)
		spacerRow = " " * width
		whiteSep : str = " " * self.sepLength
		sepReplace = f"{self.borders[0]}{whiteSep}{self.borders[1]}", f"{self.borders[0]}{self.sep}{self.borders[1]}"
		emptyEntry = " "*self.length

		rowTemplate = []
		
		namesList = textwrap.wrap(whiteSep.join(map(lambda i: f"{'{'}names[{i}]{'}'}", range(self.N))), width)
		barsList = textwrap.wrap(whiteSep.join(map(lambda i: f"{self.borders[0]}{'{'}bars[{i}]{'}'}{self.borders[1]}", range(self.N))), width)
		
		for namesRow, barsRow in zip(namesList, barsList):
			rowTemplate.append(namesRow.center(width))
			rowTemplate.append(barsRow.center(width).replace(*sepReplace))
			rowTemplate.append(spacerRow)
		rowTemplate = firstRow + "".join(rowTemplate)
		
		backspaces = "\b" * len(rowTemplate)
		whitespaces = " " * len(self.backspaces)
		
		return backspaces, whitespaces, rowTemplate

	def rowGenerator(self) -> Generator[tuple[int,str],None,None]:
		"""Progress special cases:
		None	-	Service crashed
		1 : int	-	Never ran/skipped
		1.0		-	Completed
		0<->1	-	Running
		<0		-	Not started
		>1		-	Postprocessing
		"""
		raise NotImplementedError("rowGenerator not implemented in the base class")

	def run(self):
		
		startTime = timer()
		self.running = True
		
		generator = self.rowGenerator()
		while self.running:
			self.backspaces, self.whitespaces, self.rowTemplate = self.createRowTemplate(os.get_terminal_size()[0])
			for i in range(len(self.names)):
				try:
					self.entries[i] = next(generator)
				except IndexError:
					pass
			
			s = timer()-startTime
			h, m, s, ms = s // 3600, (s % 3600) // 60, s % 60, s % 1
			if all(p is None or p == 1 for p in self.threads.values()):
				break
			print(self.backspaces, end=self.rowTemplate.format(time=f"{h:0>2d}:{m:0>2d}:{s:0>2d},{ms:0<3d}", names=self.names, bars=self.entries), flush=True, file=self.out)
			self.condition.acquire(timeout=0.2)
		
		if None in self.threads.values():
			print(self.backspaces, end=self.rowTemplate.format(time=f"{h:0>2d}:{m:0>2d}:{s:0>2d},{ms:0<3d} Failed!", names=self.names, bars=self.entries), flush=True, file=self.out)
		elif all(map((1).__eq__, self.threads.values())):
			print(self.backspaces, end=self.rowTemplate.format(time=f"{h:0>2d}:{m:0>2d}:{s:0>2d},{ms:0<3d} Done!", names=self.names, bars=self.entries), flush=True, file=self.out)
		else:
			print(self.backspaces, end=self.rowTemplate.format(time=f"{h:0>2d}:{m:0>2d}:{s:0>2d},{ms:0<3d} Interrupted!", names=self.names, bars=self.entries), flush=True, file=self.out)
		print(flush=True, file=self.out)
	
	def checkDone(self):
		if all(v is not None or v >= 1.0 for v in self.threads.values()):
			self.running = False

	def kill(self):
		self.running*=0
		self.condition.release()

class LoadingBar(Indicator):

	def __init__(self, threads, length=10, fill="\u2588", halfFill="\u258C", background=" ", **kwargs):
		super().__init__(threads, [fill, halfFill, background], **kwargs)
		self.innerLength = length - self.borderLength

	def rowGenerator(self) -> Generator[str,None,None]:
		"""Progress special cases:
		None	-	Service crashed
		1 : int	-	Never ran/skipped
		1.0		-	Completed
		0<->1	-	Running
		<0		-	Not started
		>1		-	Postprocessing
		"""
		crashString = self.crashColor+self.crashSymbol*self.innerLength
		finishString = self.finishColor+self.symbols[0]*self.innerLength
		skippedString = self.skippedColor+self.symbols[0]*self.innerLength
		n = 0
		while self.running:
			for name, prog in zip(self.names, map(self.threads.get, self.names)):
				if prog is None:
					yield crashString
				elif isinstance(prog, int) and prog == 1:
					yield skippedString
				elif name in self.finishedThreads:
					yield finishString
				elif 0 <= prog <= 1.0:
					fillLength = int(self.innerLength*2*prog)
					fillLength, halfBlock = fillLength//2, fillLength%2
					
					yield f"{self.progColor}{self.symbols[0]*fillLength}{self.symbols[1]*halfBlock}{self.partition}{self.symbols[2]*(self.innerLength - fillLength - halfBlock)}"
				elif prog < 0:
					yield self.preColor+self.symbols[2]*n + self.symbols[0] + self.symbols[2]*(self.innerLength - n - 1)
				elif prog > 1:
					yield self.postColor+self.symbols[0]*n + self.symbols[2] + self.symbols[0]*(self.innerLength - n - 1)
				else:
					yield crashString
			n = (n+1) % self.innerLength

class Spinner(Indicator):

	symbols : list[str] = ["|", "/", "-", "\\"]

	def rowGenerator(self) -> Generator[tuple[int,str],None,None]:
		"""Progress special cases:
		None	-	Service crashed
		1 : int	-	Never ran/skipped
		1.0		-	Completed
		0<->1	-	Running
		<0		-	Not started
		>1		-	Postprocessing
		"""
		NSymbols = len(self.symbols)
		crashString = self.crashColor + self.crashSymbol
		finishString = self.finishColor + self.finishSymbol
		skippedString = self.skippedColor + self.finishSymbol
		n = 0
		while self.running:
			for name, prog in zip(self.names, map(self.threads.get, self.names)):
				if prog is None:
					yield crashString
				elif isinstance(prog, int) and prog == 1:
					yield skippedString
				elif name in self.finishedThreads:
					yield finishString
				elif 0 <= prog <= 1:
					yield self.progColor + self.symbols[n]
				elif prog < 0:
					yield self.preColor + "." if n%2 else ":"
				elif 1 < prog:
					yield self.postColor + "." if n%2 else ":"
				else:
					yield crashString
			n = (n+1) % NSymbols

class TextProgress(Indicator):
	def rowGenerator(self) -> Generator[tuple[int,str],None,None]:
		"""Progress special cases:
		None	-	Service crashed
		1 : int	-	Never ran/skipped
		1.0		-	Completed
		0<->1	-	Running
		<0		-	Not started
		>1		-	Postprocessing
		"""
		NSymbols = len(self.symbols)
		crashString = self.crashColor + self.crashSymbol
		finishString = self.finishColor + self.finishSymbol
		skippedString = self.skippedColor + self.finishSymbol
		n = 0
		while self.running:
			for name, prog in zip(self.names, map(self.threads.get, self.names)):
				if prog is None:
					yield crashString
				elif isinstance(prog, int) and prog == 1:
					yield skippedString
				elif name in self.finishedThreads:
					yield finishString
				elif 0 <= prog <= 1:
					yield self.progColor + self.symbols[int(NSymbols*prog)]
				elif prog < 0:
					yield self.preColor + "."*n + ":"*(n<self.innerLength) + "."*(self.innerLength-n)
				elif 1 < prog:
					yield self.postColor + "o"*n + "O"*(n<self.innerLength) + "o"*(self.innerLength-n)
				else:
					yield crashString
			n = (n+1) % (self.innerLength + 1)

class TerminalUpdater:
	
	threads : HitchableDict
	thread : Thread
	hooks : Hooks
	out : TextIO
	printer : Indicator = None

	def __init__(self, message, category, hooks : Hooks, threadNames : Iterable, out=stdout, printer : Indicator=Spinner):
		
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
		1 : int	-	Never ran/skipped
		1.0		-	Completed
		0<->1	-	Running
		<0		-	Not started
		>1		-	Postprocessing
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
			self.threads[eventInfo["name"]] = 1
			self.printer.finishedThreads.add(eventInfo["name"])

	def startingCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["name"] in self.threads:
			self.threads[eventInfo["name"]] = 0.0

	def progressCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["name"] in self.threads and "progress" in eventInfo:
			self.threads[eventInfo["name"]] = eventInfo["progress"]

	def postProcessCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["name"] in self.threads:
			self.threads[eventInfo["name"]] = 1.1

	def finishedCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["name"] in self.threads:
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
		self.printer.condition.release()

	def stop(self):
		"""Waits for thread to stop."""
		self.printer.running = False
		self.printer.condition.release()
		self.thread.join()

	def mainLoop(self, *args, **kwargs):
		try:
			self.printer.run(*args, **kwargs)
		except Exception as e:
			LOGGER.exception(e)
			print("Exception occured in print-loop", flush=True, file=self.out)

	def updatePrinter(self):
		try:
			self.printer.threads = self.threads
			for name, value in self.threads.items():
				if name in self.printer.finishedThreads and value != 1:
					self.printer.finishedThreads.remove(name)
		except:
			pass