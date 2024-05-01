
from MetaCanSNPer.Globals import *
import MetaCanSNPer.Globals as Globals
import MetaCanSNPer.core.Hooks as Hooks
from MetaCanSNPer.core.LogKeeper import createLogger
from collections import defaultdict
from threading import Thread, _DummyThread
import sqlite3

LOGGER = createLogger(__name__)

NULL_LOGGER = logging.Logger("NULL_LOGGER", 100)

def gunzip(filepath : str):
	'''gunzip's given file. Only necessary for software that requires non-zipped data.'''
	import gzip
	outFile = open(filepath.rstrip(".gz"), "wb")

	outFile.write(gzip.decompress(open(filepath, "rb").read()))
	outFile.close()

class URL:
	string : str
	patterns : tuple[re.Pattern]
	@classmethod
	def format(self, query : str|tuple[str]=None):
		if isinstance(query, tuple):
			query = {name:value for pat, q in zip(self.patterns, query) for name, value in pat.match(q).groupdict().items()}
		else:
			{"query" : query}

		return self.string.format(**query)

class ThreadDescriptor:

	threads : dict[int,Thread] = defaultdict(list)
	func : FunctionType | MethodType
	thread : Thread = _DummyThread()
	threads : defaultdict[int,list[Thread]]
	def __init__(self, func):
		self.func = func
		self.threads = defaultdict(list)
	
	def __call__(self, *args, **kwargs):
		self.thread = Thread(target=self._funcWrapper, args=args, kwargs=kwargs, daemon=True)
		self.threads[id(getattr(self.func, "__self__", self))].append(self.thread)
		self.thread.start()
		return self
	def __repr__(self):
		return f"<{self.__class__.__name__}({self.func}) at {hex(id(self))}>"
	
	def _funcWrapper(self, *args, **kwargs):
		try:
			self.func(*args, **kwargs)
		except Exception as e:
			e.add_note(f"This exception occured in a thread running the following function call: {printCall(self.func, args, kwargs)}")
			LOGGER.exception(e)
			self.func.__self__

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

	query : Any|Iterable
	filename : str
	out : Path
	reportHook : Callable = lambda prog : None
	_queueConnection : sqlite3.Connection
	LOGGER : logging.Logger = NULL_LOGGER
	def __init__(self, query, filename, reportHook=None, out=Path("."), conn=None, *, logger=LOGGER):
		self.query = query
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
	
	def run(self, sources, postProcess : Callable=None):

		try:
			from urllib.request import urlretrieve
			reportHook = ReportHook(self.reportHook)
			for sourceName, sourceLink in sources:
				try:
					(outFile, msg) = urlretrieve(sourceLink.format(query=self.query), filename=self.out / self.filename, reporthook=reportHook) # Throws error if 404
					if postProcess is not None:
						self.reportHook(2.0)
						postProcess(outFile)
					return outFile, sourceName
				except Exception as e:
					self.LOGGER.info(f"Couldn't download from source={sourceName}, url: {sourceLink.format(filename=self.filename)!r}")
					self.LOGGER.exception(e, stacklevel=logging.DEBUG)
			self.reportHook(None)
			self.LOGGER.error(f"No database named {self.filename!r} found online. Sources tried: {', '.join(map(*this[0] + ': ' + this[1], sources))}")
			return None, None
		except Exception as e:
			self.reportHook(None)
			raise e


class Downloader:

	directory : DirectoryPath = DirectoryPath(".")
	SOURCES : tuple[tuple[str]] = ()
	reportHook : Callable = lambda prog : None
	postProcess : Callable=None
	"""Function that takes arguments: block, blockSize, totalSize"""
	timeStep : float = 0.2
	database : Path
	jobs : list
	

	_queueConnection : sqlite3.Connection
	_threads : list[Thread]= []
	LOGGER : logging.Logger = NULL_LOGGER

	def __init__(self, directory=directory, *, reportHook=None, logger=None):

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
	def download(self, query, filename : str, reportHook=reportHook) -> None:
		
		job = Job(query, filename, reportHook=reportHook, out=self.directory, conn=self._queueConnection, logger=self.LOGGER)
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
				job.run(self.SOURCES, postProcess=self.postProcess)


# Implementations

class DownloaderReportHook:

	category : str
	hooks : Hooks
	name : str
	def __init__(self, category : str, hooks : Hooks, name : str):
		self.category = category
		self.hooks = hooks
		self.name = name

	def __call__(self, prog):
		self.hooks.trigger(self.category+"Progress", {"progress" : prog, "name" : self.name})

class NCBI_URL(URL):
	string = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{idHead}/{id1_3}/{id4_6}/{id7_9}/{genome_id}_{assembly}/{genome_id}_{assembly}_genomic.fna.gz"
	patterns = [
		re.compile(r"(?P<idHead>\w{3})_(?P<id1_3>\d{3})(?P<id4_6>\d{3})(?P<id7_9>\d{3})"),
		re.compile(".*")
	]

class DatabaseDownloader(Downloader):
	SOURCES = [
		("MetaCanSNPer-data", "https://github.com/XJ-04561/MetaCanSNPer-data/raw/master/database/{query}"), # MetaCanSNPer
		("CanSNPer2-data", "https://github.com/FOI-Bioinformatics/CanSNPer2-data/raw/master/database/{query}") # Legacy CanSNPer
	]
	LOGGER = LOGGER

class ReferenceDownloader(Downloader):
	SOURCES = [
		("NCBI", NCBI_URL)
	]
	LOGGER = LOGGER
	postProcess = gunzip