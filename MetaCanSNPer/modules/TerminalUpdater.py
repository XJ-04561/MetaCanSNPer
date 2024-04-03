
from MetaCanSNPer.modules.Hooks import Hooks
from threading import Thread, Condition
from time import sleep
from timeit import default_timer as timer
from sys import stdout
from typing import TextIO, Iterable
from functools import cached_property
from MetaCanSNPer.Globals import *
from MetaCanSNPer.modules.LogKeeper import createLogger

LOGGER = createLogger(__name__)

class HitchableDict(dict):

	onSet : Callable = None

	def __new__(cls, *args, **kwargs):
		obj = super().__init__(cls, *args, **kwargs)
		return obj

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

	symbols : list[str]
	sep : str
	borders : tuple[str,str]
	crashSymbol : str
	finishSymbol : str
	_threads : dict[str,float]
	condition : Condition

	backspaces : str
	whitespaces : str

	sepLength : int
	borderLength : int
	innerLength : int
	
	@property
	def threads(self) -> HitchableDict[str,float]:
		return self._threads
	
	@threads.setter
	def threads(self, threads : HitchableDict[str,float]):
		self._threads = threads
		self.keys = sorted(threads.keys())
		self.N = len(self.keys)
		self.createRowTemplate()

	def __init__(self, threads : HitchableDict, symbols : list[str], message="", sep=" ", borders=("[", "]"), crashSymbol=None, finishSymbol="\u2588", out=stdout, time: bool=True, preColor : str=None, partition : str=None, crashColor : str=None, skippedColor : str=None, finishColor : str=None, postColor : str=None, progColor : str=None):
		
		self.preColor		= preColor or ("\u001b[33;40m" if supportsColor() else "")
		self.progColor		= progColor or ("\u001b[35;40m" if supportsColor() else "")
		self.partition		= partition or ("\u001b[37;40m" if supportsColor() else "")
		self.crashColor		= crashColor or ("\u001b[31;40m" if supportsColor() else "")
		self.skippedColor	= skippedColor or ("\u001b[31;44m" if supportsColor() else "")
		self.finishColor	= finishColor or ("\u001b[32;40m" if supportsColor() else "")
		self.postColor		= postColor or ("\u001b[36;40m" if supportsColor() else "")

		self.message = message
		self.condition = Condition()
		self.time = self.time
		self.out = out

		self.symbols		= symbols
		self.innerLength 	= len(self.symbols[0])
		self.sep			= sep
		self.borders		= borders
		self.crashSymbol	= crashSymbol
		self.finishSymbol	= finishSymbol
		
		self.sepLength = len(self.sep)
		self.borderLength = len(self.borders[0]) + len(self.borders[1])

		self.threads = threads
	
	def createRowTemplate(self):
		if self.time is True:
			self.rowTemplate = self.sep.join([f"{self.borders[0]}"+"{targets["+str(i)+"]}"+f"{self.borders[1]}" for i in range(self.N)]) + " t = {seconds:>15.3f}"
			self.backspaces = "\b" * ((self.innerLength+self.borderLength)*self.N + self.sepLength*max(0, self.N-1) + 20)
		else:
			self.rowTemplate = self.sep.join(f"{self.borders[0]}"+"{finishSymbol}"+f"{self.borders[1]}")
			self.backspaces = "\b" * ((self.innerLength+self.borderLength)*self.N + self.sepLength*max(0, self.N-1))
		self.whitespaces = self.backspaces.replace("\b", " ")
		self.entries = list(map(lambda _:" "*self.innerLength, range(self.N)))

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
		N = self.N
		
		print(self.message, end=self.whitespaces, flush=True, file=self.out)
		self.printedLength = len(self.message)
		
		generator = self.rowGenerator()
		while self.running:
			for i in range(len(self.keys)):
				try:
					self.entries[i] = next(generator)
				except IndexError:
					pass
			print(self.backspaces, end=self.rowTemplate.format(seconds=timer()-startTime, targets=self.entries), flush=True, file=self.out)
			self.condition.acquire(timeout=0.2)
	
	def checkDone(self):
		if all(v is not None or v >= 1.0 for v in self.threads.values()):
			self.running = False

	def kill(self):
		self.running*=0
		self.condition.release()

	def cleanup(self):
		print(self.backspaces, end="\b"*self.printedLength, flush=True, file=self.out)

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
			for prog in map(self.threads.get, self.keys):
				if prog is None:
					yield crashString
				elif isinstance(prog, int) and prog == 1:
					yield skippedString
				elif prog == 1.0:
					yield finishString
				elif 0 <= prog < 1.0:
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
			for prog in map(self.threads.get, self.keys):
				if prog is None:
					yield crashString
				elif isinstance(prog, int) and prog == 1:
					yield skippedString
				elif prog == 1.0:
					yield finishString
				elif 0 <= prog < 1:
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
			for prog in map(self.threads.get, self.keys):
				if prog is None:
					yield crashString
				elif isinstance(prog, int) and prog == 1:
					yield skippedString
				elif prog == 1.0:
					yield finishString
				elif 0 <= prog < 1:
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
		self.threads = HitchableDict(map(lambda name : (name,0.0), threadNames))
		self.threads.onSet
		self.out = out
		self.finishedThreads = set()
		
		self.running = True

		"""Progress special cases:
		None	-	Service crashed
		1.0		-	Completed
		0<->1	-	Running
		<0		-	Not started
		>1		-	Postprocessing
		"""
		self.hooks.addHook(f"{self.category}Initialized", self.initializedCallback)
		self.hooks.addHook(f"{self.category}Skipped", self.skippedCallback)
		self.hooks.addHook(f"{self.category}Starting", self.startingCallback)
		self.hooks.addHook(f"{self.category}Progress", self.progressCallback)
		self.hooks.addHook(f"{self.category}PostProcess", self.postProcessCallback)
		self.hooks.addHook(f"{self.category}Finished", self.finishedCallback)
		self.hooks.addHook(f"{self.category}Failed", self.failedCallback)

		if printer is not None:
			self.setPrinter(printer)

	def initializedCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["thread"] in self.threads:
			self.threads[eventInfo["thread"]] = -1.0

	def skippedCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["thread"] in self.threads:
			self.threads[eventInfo["thread"]] = 1

	def startingCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["thread"] in self.threads:
			self.threads[eventInfo["thread"]] = 0.0

	def progressCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["thread"] in self.threads:
			self.threads[eventInfo["thread"]] = eventInfo["progress"]

	def postProcessCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["thread"] in self.threads:
			self.threads[eventInfo["thread"]] = 1.1

	def finishedCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["thread"] in self.threads:
			self.threads[eventInfo["thread"]] = 1.0

	def failedCallback(self, eventInfo : dict[str,Any]):
		if eventInfo["thread"] in self.threads:
			self.threads[eventInfo["thread"]] = None
	
	def setPrinter(self, printer : Indicator=Spinner, *args, **kwargs):
		self.printer = printer(self.threads, *args, **({"out":self.out, "message":self.message} | kwargs))

	def start(self, *args, **kwargs):
		self.thread = Thread(target=self.mainLoop, args=args, kwargs=kwargs, daemon=True)
		self.thread.start()
	
	def wait(self, timeout=3):
		self.thread.join(timeout=timeout)
	
	def kill(self):
		"""Does not wait for thread to stop."""
		self.running *= False
		self.printer.condition.release()

	def stop(self):
		"""Waits for thread to stop."""
		self.running *= False
		self.printer.condition.release()
		self.thread.join()
	
	def cleanup(self):
		self.printer.cleanup()

	def mainLoop(self, *args, **kwargs):
		try:
			self.printer.run(*args, **kwargs)
		except Exception as e:
			LOGGER.exception(e)
			print("Exception occured in print-loop", flush=True, file=self.out)
