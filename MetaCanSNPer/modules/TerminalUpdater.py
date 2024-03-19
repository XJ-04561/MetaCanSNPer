
from MetaCanSNPer.modules.Hooks import Hooks
from threading import Thread
from time import sleep
from sys import stdout
from typing import TextIO

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
		self.threads = {key+1:0.0 for key in range(nThreads)}
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
	
	def wait(self):
		self.thread.join(timeout=3)
	
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

	"""Print Functions"""

	def showLoadingSymbol(self, symbols : list[str]=["|", "/", "-", "\\"], sep=" ", borders=("[", "]")):

		keys = sorted(self.threads.keys())
		sepLength = len(sep)
		N = len(self.threads)
		borderLength = len(borders[0]) + len(borders[1])
		backspaces = "\b" * ((1+borderLength)*N + sepLength*max(0, N-1))
		m = len(symbols)
		n = [0 for _ in range(len(self.threads))]
		print(f"{self.message} ... ", end="", flush=True, file=self.out)
		while self.running:
			for i, key in enumerate(keys):
				prog = self.threads[key]
				if self.running:
					print(f"{borders[0]}{symbols[n[i]]}{borders[1]}", end="" if key == keys[-1] else sep, flush=True, file=self.out)
					n[i]=(n[i]+1)%m
				else:
					backspaces = "\b" * ((1+borderLength+sepLength)*i)
					print(backspaces+backspaces.replace("\b", " ")+backspaces, end="", flush=True, file=self.out)
					return print("Done!", flush=True, file=self.out)
			sleep(0.2)
			print(backspaces, end="", flush=True, file=self.out)
		print(backspaces.replace("\b", " ")+backspaces, end="", flush=True, file=self.out)
		print("Done!" if self.finishedThreads.issuperset(self.threads) else "Failed!", flush=True, file=self.out)

	def showLoadingMiniBars(self, symbols : list[str]= [".", "_", "\u2584", "", "\u2588"], sep=" ", borders=("[", "]")):

		keys = sorted(self.threads.keys())
		sepLength = len(sep)
		N = len(self.threads)
		borderLength = len(borders[0]) + len(borders[1])
		backspaces = "\b" * ((1+borderLength)*N + sepLength*max(0, len(self.threads)-1))
		print(f"{self.message} ... ", end="", flush=True, file=self.out)
		while self.running:
			for i, key in enumerate(keys):
				prog = self.threads[key]
				if self.running:
					print(f"{borders[0]}{symbols[int(N*prog)]}{borders[1]}", end="" if key == keys[-1] else sep, flush=True, file=self.out)
				else:
					backspaces = "\b" * ((1+borderLength+sepLength)*i)
					print(backspaces+backspaces.replace("\b", " ")+backspaces, end="", flush=True, file=self.out)
					return print("Done!", flush=True, file=self.out)
			sleep(0.5)
			print(backspaces, end="", flush=True, file=self.out)
		print(backspaces.replace("\b", " ")+backspaces, end="", flush=True, file=self.out)
		print("Done!" if self.finishedThreads.issuperset(self.threads) else "Failed!", flush=True, file=self.out)

	def showLoadingBar(self, length=10, borders=("[", "]"), fill="\u2588", halfFill="\u258C", background=" ", sep=" "):

		keys = sorted(self.threads.keys())
		innerLength = length - len(borders[0]) - len(borders[1])
		sepLength = len(sep)
		backspaces = "\b" * (length * len(self.threads) + sepLength * max(0, len(self.threads) - 1))
		print(f"{self.message} ... ", end="", flush=True, file=self.out)
		while self.running:
			for i, key in enumerate(keys):
				prog = self.threads[key]
				if self.running:
					fillLength = int(innerLength*2*prog)
					fillLength, halfBlock = fillLength//2, fillLength%2
					emptyLength = innerLength - fillLength - halfBlock
					
					print(f"{borders[0]}{fill*fillLength}{halfFill*halfBlock}{background*emptyLength}{borders[1]}", end="" if key == keys[-1] else sep, flush=True, file=self.out)
				else:
					backspaces = "\b" * (length+sepLength) * i
					print(backspaces+backspaces.replace("\b", " ")+backspaces, end="", flush=True, file=self.out)
					return print("Done!", flush=True, file=self.out)
			sleep(0.6)
			print(backspaces, end="", flush=True, file=self.out)
		print(backspaces.replace("\b", " ")+backspaces, end="", flush=True, file=self.out)
		print("Done!" if self.finishedThreads.issuperset(self.threads) else "Failed!", flush=True, file=self.out)
