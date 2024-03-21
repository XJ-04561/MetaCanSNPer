
from MetaCanSNPer.modules.Hooks import Hooks
from threading import Thread
from time import sleep
from sys import stdout
from typing import TextIO
from functools import cached_property
from MetaCanSNPer.modules.LogKeeper import createLogger

LOGGER = createLogger(__name__)

class TerminalUpdater:
	
	threads : dict
	thread : Thread
	hooks : Hooks
	out : TextIO
	def __init__(self, message, category, hooks : Hooks, nThreads, out=stdout):
		
		if nThreads < 1:
			raise ValueError(f"'nThreads' must have a value over zero.")
		
		self.message = message
		self.category = category
		self.hooks = hooks
		self.threads = {key:0.0 for key in range(nThreads)}
		self.out = out
		self.finishedThreads = set()
		
		self.running = True

		self.hooks.addHook(f"{self.category}Progress", self.updateProgress)
		self.hooks.addHook(f"{self.category}Finished", self.stopLoadingbar)

		self.setPrintFunc(self.showLoadingSymbol)
	
	def setPrintFunc(self, printFunc, args=[], kwargs={}):
		self.thread = Thread(target=printFunc, args=args, kwargs=kwargs, daemon=True)

	def start(self):
		self.thread.start()
	
	def wait(self, timeout=3):
		self.thread.join(timeout=timeout)

	def deadmans(self, timeout=3):
		self.thread.join(timeout=timeout)
		self.running *= False
		self.thread.join()
	
	def kill(self):
		self.running *= False

	def stop(self):
		self.running *= False
		self.thread.join()
	
	def updateProgress(self, eventInfo : dict):
		self.threads[eventInfo["threadN"]] *= 0.0
		self.threads[eventInfo["threadN"]] += eventInfo["progress"]
		
	def stopLoadingbar(self, eventInfo : dict):
		self.finishedThreads.add(eventInfo["threadN"])
		if self.finishedThreads.issuperset(self.threads):
			self.running *= False

	# Taken from django @ https://github.com/django/django/blob/main/django/core/management/color.py
	@cached_property
	def supportsColor(self):
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

	"""Print Functions"""

	def showLoadingSymbol(self, symbols : list[str]=("|", "/", "-", "\\"), sep=" ", borders=("[", "]"), crashSymbol="X"):

		try:
			keys = sorted(self.threads.keys())
			sepLength = len(sep)
			N = len(keys)
			borderLength = len(borders[0]) + len(borders[1])
			backspaces = "\b" * ((len(symbols[0])+borderLength)*N + sepLength*max(0, N-1))
			m = len(symbols)
			n = [0 for _ in range(N)]

			if self.supportsColor:
				borders = ("\u001b[37;40m"+borders[0], "\u001b[37;40m"+borders[1]+"\u001b[0m")
				symbols = list(map(lambda x : "\u001b[36;40m"+x, symbols))
				crashSymbol = f"\u001b[31;40m{crashSymbol}"

			print(f"{self.message} ... ", end=backspaces.replace("\b", " "), flush=True, file=self.out)
			while self.running:
				msg = ""
				for i, key in enumerate(self.threads):
					prog = self.threads[key]
					if self.running:
						msg += f"{borders[0]}{symbols[n[i]] if prog is not None else crashSymbol}{borders[1]}{sep if i < N-1 else ''}"
						n[i]=(n[i]+1)%m
					else:
						print(backspaces+backspaces.replace("\b", " "), end=backspaces, flush=True, file=self.out)
						return print("Done!" if self.finishedThreads.issuperset(self.threads) else "Failed!", flush=True, file=self.out)
				print(backspaces, end=msg, flush=True, file=self.out)
				sleep(0.2)
			print(backspaces+backspaces.replace("\b", " "), end=backspaces, flush=True, file=self.out)
			print("Done!" if self.finishedThreads.issuperset(self.threads) else "Failed!", flush=True, file=self.out)
		except Exception as e:
			LOGGER.exception(e)
			print("", flush=True, file=self.out)

	def showLoadingMiniBars(self, symbols : list[str]= [".", "_", "\u2584", "#", "\u2588"], sep=" ", borders=("[", "]"), crashSymbol : str="/"):

		try:
			keys = sorted(self.threads.keys())
			sepLength = len(sep)
			N = len(keys)
			borderLength = len(borders[0]) + len(borders[1])
			backspaces = "\b" * ((len(symbols[0])+borderLength)*N + sepLength*max(0, len(self.threads)-1))
			
			if self.supportsColor:
				borders = ("\u001b[37;40m"+borders[0], "\u001b[37;40m"+borders[1]+"\u001b[0m")
				symbols = list(map(lambda x : "\u001b[36;40m"+x, symbols))
				crashSymbol = f"\u001b[31;40m{crashSymbol}"

			print(f"{self.message} ... ", end=backspaces.replace("\b", " "), flush=True, file=self.out)
			while self.running:
				msg = ""
				for i, key in enumerate(keys):
					prog = self.threads[key]
					if self.running:
						if prog is not None:
							msg += f"{borders[0]}{symbols[int(N*prog)]}{borders[1]}{sep if i < N-1 else ''}"
						else:
							msg += f"{borders[0]}{crashSymbol}{borders[1]}{sep if i < N-1 else ''}"
					else:
						print(backspaces+backspaces.replace("\b", " "), end=backspaces, flush=True, file=self.out)
						return print("Done!", flush=True, file=self.out)
				print(backspaces, end=msg, flush=True, file=self.out)
				sleep(0.5)
			print(backspaces+backspaces.replace("\b", " "), end=backspaces, flush=True, file=self.out)
			print("Done!" if self.finishedThreads.issuperset(self.threads) else "Failed!", flush=True, file=self.out)
		except Exception as e:
			LOGGER.exception(e)
			print("", flush=True, file=self.out)

	def showLoadingBar(self, length=10, borders=("[", "]"), fill="\u2588", halfFill="\u258C", background=" ", sep=" ", partition="", crashSymbol : str="/"):

		try:
			innerLength = length - len(borders[0]) - len(borders[1]) - len(partition)
			keys = sorted(self.threads.keys())
			sepLength = len(sep)
			N = len(keys)
			backspaces = "\b" * (length * N + sepLength * max(0, N - 1) + len(partition)*N)
			
			crashColor = ""
			if self.supportsColor:
				borders = ("\u001b[37;40m"+borders[0]+"\u001b[32;40m", "\u001b[37;40m"+borders[1]+"\u001b[0m")
				partition = "\u001b[31;40m"
				crashColor = "\u001b[31;40m"

			print(f"{self.message} ... ", end=backspaces.replace("\b", " "), flush=True, file=self.out)
			while self.running:
				msg = ""
				for i, key in enumerate(keys):
					prog = self.threads[key]
					if self.running:
						if prog is not None:
							fillLength = int(innerLength*2*prog)
							fillLength, halfBlock = fillLength//2, fillLength%2
							emptyLength = innerLength - fillLength - halfBlock
							
							msg += f"{borders[0]}{fill*fillLength}{partition}{halfFill*halfBlock}{background*emptyLength}{borders[1]}{sep if i < N-1 else ''}"
						else:
							msg += f"{borders[0]}{crashColor}{innerLength*crashSymbol}{borders[1]}{sep if i < N-1 else ''}"
					else:
						print(backspaces+backspaces.replace("\b", " "), end=backspaces, flush=True, file=self.out)
						return print("Done!", flush=True, file=self.out)
				print(backspaces, end=msg, flush=True, file=self.out)
				sleep(0.6)
			print(backspaces+backspaces.replace("\b", " "), end=backspaces, flush=True, file=self.out)
			print("Done!" if self.finishedThreads.issuperset(self.threads) else "Failed!", flush=True, file=self.out)
		except Exception as e:
			LOGGER.exception(e)
			print("", flush=True, file=self.out)
