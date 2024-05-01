
import sqlite3, logging
from types import FunctionType, MethodType
from collections import defaultdict
from threading import Thread, _DummyThread
from typing import Callable
from time import sleep

from This import this
from PseudoPathy import Path, DirectoryPath
from PseudoPathy.PathShortHands import *

NULL_LOGGER = logging.Logger("NULL_LOGGER", 100)

class ThreadDescriptor:

	threads : dict[int,Thread] = defaultdict(list)
	func : FunctionType | MethodType
	thread : Thread = _DummyThread()
	threads : defaultdict[int,list[Thread]]
	def __init__(self, func):
		self.func = func
		self.threads = defaultdict(list)
	
	def __call__(self, *args, **kwargs):
		self.thread = Thread(target=self.func, args=args, kwargs=kwargs, daemon=True)
		self.threads[id(getattr(self.func, "__self__", self))].append(self.thread)
		self.thread.start()
		return self
	def __repr__(self):
		return f"<{self.__class__.__name__}({self.func}) at {hex(id(self))}>"
	
	def wait(self):
		self.thread.join()
	
class ThreadMethod(ThreadDescriptor):

	func : MethodType

class ThreadFunction(ThreadDescriptor):
	
	func : FunctionType
	
	def __get__(self, instance, owner=None):
		return ThreadMethod(self.func.__get__(instance, owner=owner))
	def __set_name__(self, instance, name):
		type(self.func).__set_name__(self, instance, name)
		self.func.__set_name__(self, instance, name)

def threadDescriptor(func):
	"""The actual decorator to use."""
	return ThreadFunction(func)

class ReportHook:

	totalBlocks : int = None
	def __init__(self, reportHook):
		self.reportHook = reportHook
	
	def __call__(self, block, blockSize, totalSize):
		if self.totalBlocks is None:
			self.totalBlocks = (totalSize // blockSize) + 1
		
		self.reportHook(block / self.totalBlocks)

class Job:

	filename : str
	out : Path
	reportHook : Callable = lambda prog : None
	_queueConnection : sqlite3.Connection
	LOGGER : logging.Logger = NULL_LOGGER
	def __init__(self, filename, reportHook=None, out=Path("."), conn=None, *, logger=LOGGER):
		self.filename = filename
		self.reportHook = reportHook or self.reportHook
		self.out = out
		self._queueConnection = conn
		self.LOGGER = logger
	
	def reserveQueue(self):
		try:
			self._queueConnection.execute("INSERT OR FAIL INTO queueTable (name) VALUES (?);", [self.filename])
		except sqlite3.IntegrityError:
			return False
		else:
			return True

	def isListed(self):
		return self._queueConnection.execute("SELECT CASE WHEN EXISTS(SELECT 1 FROM queueTable WHERE name = ?) THEN TRUE ELSE FALSE;", [self.filename]).fetchone()[0]
	def isQueued(self):
		return self._queueConnection.execute("SELECT CASE WHEN EXISTS(SELECT 1 FROM queueTable WHERE name = ? AND progress < 0.0) THEN TRUE ELSE FALSE;", [self.filename]).fetchone()[0]
	def isDownloading(self):
		return self._queueConnection.execute("SELECT CASE WHEN EXISTS(SELECT 1 FROM queueTable WHERE name = ? AND progress >= 0.0 AND progress < 1.0) THEN TRUE ELSE FALSE;", [self.filename]).fetchone()[0]
	def isDone(self):
		return self._queueConnection.execute("SELECT CASE WHEN EXISTS(SELECT 1 FROM queueTable WHERE name = ? AND progress = 1.0) THEN TRUE ELSE FALSE;", [self.filename]).fetchone()[0]
	def isPostProcess(self):
		return self._queueConnection.execute("SELECT CASE WHEN EXISTS(SELECT 1 FROM queueTable WHERE name = ? AND progress > 1.0) THEN TRUE ELSE FALSE;", [self.filename]).fetchone()[0]
	def isDead(self):
		return self._queueConnection.execute("SELECT CASE WHEN EXISTS(SELECT 1 FROM queueTable WHERE name = ? AND timestamp + 10.0 < julianday(CURRENT_TIMESTAMP)) THEN TRUE ELSE FALSE;", [self.filename]).fetchone()[0]
	
	def updateProgress(self, name, prog : float):
		self._queueConnection.execute("UPDATE queueTable SET progress = ?, timestamp = julianday(CURRENT_TIMESTAMP) WHERE name = ?;", [prog, name])
		self.reportHook(prog)

	def getProgress(self):
		if (ret := self._queueConnection.execute("SELECT progress FROM queueTable WHERE name = ?;", [self.filename]).fetchone()) is None:
			return None
		else:
			return ret[0]
	
	def updateLoop(self, timeStep : float=0.25):
		
		while not self.isDone():
			if self.isDead():
				self.reportHook(None)
				break
			self.reportHook(self.getProgress())
			sleep(timeStep)
		else:
			self.reportHook(1.0)
	
	def run(self, sources):

		from urllib.request import urlretrieve
		reportHook = ReportHook(self.reportHook)
		outFile = self.out / self.filename
		for sourceName, sourceLink in sources:
			try:
				(outFile, msg) = urlretrieve(sourceLink.format(filename=self.filename), filename=outFile, reporthook=reportHook) # Throws error if 404
				return outFile, sourceName
			except Exception as e:
				self.LOGGER.info(f"Couldn't download from source={sourceName}, url: {sourceLink.format(filename=self.filename)!r}")
				self.LOGGER.exception(e, stacklevel=logging.DEBUG)
		self.reportHook(None)
		self.LOGGER.error(f"No database named {self.filename!r} found online. Sources tried: {', '.join(map(*this[0] + ': ' + this[1], sources))}")
		return None, None


class Downloader:

	directory : DirectoryPath = DirectoryPath(".")
	SOURCES : tuple[tuple[str]] = ()
	reportHook : Callable = lambda prog : None
	"""Function that takes arguments: block, blockSize, totalSize"""
	timeStep : float = 0.2
	database : Path
	jobs : list

	_queueConnection : sqlite3.Connection
	_threads : list[Thread]= []
	LOGGER : logging.Logger = NULL_LOGGER

	def __init__(self, directory=directory, *, logger=None, reportHook=None):

		if pAccess(directory, "rw"):
			self.directory = directory
		elif not pExists(directory) and pMakeDirs(directory, "rw"):
			self.directory = directory
		else:
			raise PermissionError(f"Missing read and/or write permissions in directory: {directory}")
		self.database = f"{self.__module__}_QUEUE.db"
		self.jobs = []
		
		self._queueConnection = sqlite3.connect(self.directory / self.database)
		self._queueConnection.execute("CREATE TABLE IF NOT EXISTS queueTable (name TEXT UNIQUE, progress DECIMAL default -1.0, modified INTEGER DEFAULT julianday(CURRENT_TIMESTAMP));")
		if logger is not None:
			self.LOGGER = logger
		if reportHook is not None:
			self.reportHook = reportHook

	def addSources(self, *sources : str):
		self.SOURCES = self.SOURCES + sources

	def wait(self):
		threadList = self.download.threads[id(self)]
		threadSet = set(threadList)
		for t in threadSet:
			try:
				t.join()
				self.download.threads[id(self)].remove(t)
			except:
				pass
	
	@threadDescriptor
	def download(self, filename : str, reportHook=reportHook, *, logger=None) -> None:
		
		job = Job(filename, reportHook=reportHook, out=self.directory, conn=self._queueConnection, logger=logger or self.LOGGER)
		self.jobs.append(job)
		if job.reserveQueue():
			job.run()
		elif job.isDone():
			reportHook(int(1))
		else:
			while not job.reserveQueue():
				job.updateLoop(timeStep=self.timeStep)
				if job.reserveQueue():
					break
			else:
				job.run(self.SOURCES)
