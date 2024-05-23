
from MetaCanSNPer.Globals import *
import MetaCanSNPer.Globals as Globals
import MetaCanSNPer.core.Hooks as Hooks
from MetaCanSNPer.core.LogKeeper import createLogger
from collections import defaultdict
from threading import Thread, _DummyThread, Lock, Semaphore, Condition, current_thread
import sqlite3
from queue import Queue, Empty as EmptyQueueException
from urllib.request import urlretrieve, HTTPError

LOGGER = Globals.LOGGER.getChild("Downloader")

def correctDatabase(filename, finalFilename):
	from MetaCanSNPer.modules.Database import MetaCanSNPerDatabase, NoChromosomesInDatabase
	database = MetaCanSNPerDatabase(filename, "w")
	database.fix()
	if not database.valid:
		e = database.exception
		if not isinstance(e, NoChromosomesInDatabase):
			e.add_note(f"Tables hash: {database.tablesHash}\nIndexes Hash: {database.indexesHash}")
			raise e
	database.close()
	os.rename(filename, finalFilename)

def gunzip(filepath : str, outfile : str):
	'''gunzip's given file. Only necessary for software that requires non-zipped data.'''
	import gzip
	open(outfile, "wb").write(
		gzip.decompress(open(filepath, "rb").read())
	)

class URL:
	string : str
	patterns : tuple[re.Pattern]
	@classmethod
	def format(self, query : str|tuple[str]=None) -> str:
		if isinstance(query, tuple):
			query = {name:value for pat, q in zip(self.patterns, query) for name, value in pat.match(q).groupdict().items()}
		else:
			{"query" : query}

		return self.string.format(**query)


class DatabaseThread:

	LOG : logging.Logger = LOGGER.getChild("DatabaseThread")

	queue : Queue[list[str,list,Lock, list]]
	queueLock : Lock
	running : bool
	
	filename : str
	_thread : Thread
	_connection : sqlite3.Connection

	def __init__(self, filename : str):
		self.running = True
		self.queue = Queue()
		self.queueLock = Lock()
		self.queueLock.acquire()
		self.filename = filename
		self._thread = Thread(target=self.mainLoop, daemon=True)
		self._thread.start()

	def mainLoop(self):
		try:
			self._connection = sqlite3.connect(self.filename, autocommit=True)
			while self.running:
				try:
					string, params, lock, results = self.queue.get(timeout=15)
					try:
						results.extend(self._connection.execute(string, params).fetchall())
					except Exception as e:
						results.append(e)
					lock.release()
					self.queue.task_done()
				except EmptyQueueException:
					pass
				except Exception as e:
					self.LOG.exception(e)
		except Exception as e:
			self.LOG.exception(e)
			try:
				results.append(e)
				lock.release()
			except:
				pass
			for _ in range(self.queue.unfinished_tasks):
				string, params, lock, results = self.queue.get(timeout=15)
				lock.release()
		self._connection.close()

	def execute(self, string : str, params : list=[]):
		lock = Lock()
		lock.acquire()
		results = []
		self.queue.put([string, params, lock, results])
		lock.acquire()
		
		if results and isinstance(results[-1], Exception):
			raise results[-1]
		
		return results
	
	def executemany(self, *statements : tuple[str, list]):
		fakeLock = lambda :None
		fakeLock.release = lambda :None

		with self.queueLock:
			results = [[] for _ in len(statements)]
			for i, statement in enumerate(statements[:-1]):
				self.queue.put([*statement, fakeLock, results[i]])
			lock = Lock()
			lock.acquire()
			self.queue.put([*statements[-1], lock, results[-1]])
		lock.acquire()
		
		if any(r and isinstance(r[-1], Exception) for r in results):
			raise next(filter(lambda r:r and isinstance(r[-1], Exception), results))[-1]
		
		return results

	def __del__(self):
		self.running = False

class ThreadDescriptor:

	LOG : logging.Logger = LOGGER.getChild("ThreadDescriptor")

	func : FunctionType | MethodType
	owner : "Downloader"
	def __init__(self, func, owner=None):
		self.func = func
		self.owner = owner
	
	def __call__(self, *args, **kwargs):
		
		self.LOG.info(f"Starting thread for {printCall(self.func, args, kwargs)}")

		thread = Thread(target=self._funcWrapper, args=args, kwargs=kwargs, daemon=True)
		self.owner._threads.append(thread)
		thread.start()

		self.LOG.info(f"Started thread for {printCall(self.func, args, kwargs)}")
		return thread
	def __repr__(self):
		return f"<{self.__class__.__name__}({self.func}) at {hex(id(self))}>"
	
	def __set_name__(self, owner, name):
		if hasattr(self.func, "__set_name__"):
			self.func.__set_name__(self, owner, name)
		self.__name__ = name
		self.owner = owner

	def _funcWrapper(self, *args, **kwargs):
		try:
			self.func(*args, **kwargs)
		except Exception as e:
			current_thread().exception = e
			e.add_note(f"This exception occured in a thread running the following function call: {printCall(self.func, args, kwargs)}")
			self.LOG.exception(e)

	def wait(self):
		self.owner.wait()
	
class ThreadMethod(ThreadDescriptor):

	func : MethodType

class ThreadFunction(ThreadDescriptor):
	
	func : FunctionType
	
	def __get__(self, instance, owner=None):
		if instance is not None:
			return ThreadMethod(self.func.__get__(instance), instance)
		else:
			return ThreadMethod(self.func.__get__(owner), owner)

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

	LOG : logging.Logger = LOGGER.getChild("Job")

	_queueConnection : DatabaseThread

	query : Any|Iterable
	filename : str
	out : Path = Path(".")
	reportHook : Callable = None

	def __init__(self, query, filename, conn, reportHook=None, out=None, *, logger=None):
		
		if logger is not None:
			self.LOG = logger
		self._queueConnection = conn

		self.query = query
		self.filename = filename
		if reportHook:
			self.reportHook = reportHook
		if out:
			self.out = out
	
	def reserveQueue(self):
		try:
			self._queueConnection.execute("INSERT OR FAIL INTO queueTable (name) VALUES (?);", [self.filename])
		except sqlite3.IntegrityError:
			return False
		else:
			return True

	def isListed(self):
		return self._queueConnection.execute("SELECT CASE WHEN EXISTS(SELECT 1 FROM queueTable WHERE name = ?) THEN TRUE ELSE FALSE END;", [self.filename])[0][0]
	def isQueued(self):
		return self._queueConnection.execute("SELECT CASE WHEN EXISTS(SELECT 1 FROM queueTable WHERE name = ? AND progress < 0.0) THEN TRUE ELSE FALSE END;", [self.filename])[0][0]
	def isDownloading(self):
		return self._queueConnection.execute("SELECT CASE WHEN EXISTS(SELECT 1 FROM queueTable WHERE name = ? AND progress >= 0.0 AND progress < 1.0) THEN TRUE ELSE FALSE END;", [self.filename])[0][0]
	def isDone(self):
		return self._queueConnection.execute("SELECT CASE WHEN EXISTS(SELECT 1 FROM queueTable WHERE name = ? AND progress = 1.0) THEN TRUE ELSE FALSE END;", [self.filename])[0][0]
	def isPostProcess(self):
		return self._queueConnection.execute("SELECT CASE WHEN EXISTS(SELECT 1 FROM queueTable WHERE name = ? AND progress > 1.0) THEN TRUE ELSE FALSE END;", [self.filename])[0][0]
	def isDead(self):
		return self._queueConnection.execute("SELECT CASE WHEN EXISTS(SELECT 1 FROM queueTable WHERE name = ? AND modified + 10.0 < UNIXEPOCH()) THEN TRUE ELSE FALSE END;", [self.filename])[0][0]
	
	def updateProgress(self, prog : float):
		self._queueConnection.execute("UPDATE queueTable SET progress = ?, modified = UNIXEPOCH() WHERE name = ?;", [prog, self.filename])
		self.reportHook(prog)

	def getProgress(self):
		ret = self._queueConnection.execute("SELECT progress FROM queueTable WHERE name = ?;", [self.filename])
		if not ret:
			return None
		else:
			return ret[0][0]
	
	def updateLoop(self, timeStep : float=0.25):
		
		while not self.isDone():
			if self.isDead():
				self.reportHook(None)
				break
			self.reportHook(self.getProgress())
			sleep(timeStep)
		else:
			self.reportHook(1.0)
	
	def run(self, sources : list[tuple[str, type[URL]]], postProcess : Callable=None):

		try:
			reportHook = ReportHook(self.updateProgress)
			outFile = "N/A"
			for sourceName, sourceLink in sources:
				try:
					link = sourceLink.format(query=self.query)
					filename = link.rsplit("/")[-1]
					(outFile, msg) = urlretrieve(link, filename=self.out / filename, reporthook=reportHook) # Throws error if 404

					if postProcess is not None:
						self.updateProgress(2.0)
						postProcess(outFile, self.out / self.filename)
					else:
						os.rename(outFile, self.out / self.filename)
					self.updateProgress(1.0)
					return outFile, sourceName

				except HTTPError as e:
					self.LOG.info(f"Couldn't download from source={sourceName}, url: {sourceLink.format(query=self.query)}, due to {e.args}")
					# self.LOG.exception(e, stacklevel=logging.DEBUG)
				except Exception as e:
					e.add_note(f"This occurred while processing {outFile} downloaded from {sourceLink.format(query=self.query)}")
					self.LOG.exception(e)
					raise e
			self.updateProgress(None)
			self.LOG.error(f"No database named {self.filename!r} found online. Sources tried: {', '.join(map(*this[0] + ': ' + this[1], sources))}")
			raise e
		except Exception as e:
			self.updateProgress(None)
			raise e

class Downloader:

	LOG : logging.Logger = LOGGER.getChild("Downloader")
	SOURCES : tuple[tuple[str]] = ()

	directory : DirectoryPath = DirectoryPath(".")
	reportHook : Callable=None
	postProcess : Callable=None
	"""Function that takes arguments: block, blockSize, totalSize"""
	timeStep : float = 0.2
	database : Path = f"{__module__}_QUEUE.db"
	jobs : list
	hooks : Hooks = cached_property(lambda self:Hooks())

	_queueConnection : DatabaseThread
	_threads : list[Thread]= []

	def __init__(self, directory=directory, *, reportHook=None, logger=None, hooks=None, threads=None):

		if pAccess(directory, "rw"):
			self.directory = directory
		elif not pExists(directory) and pMakeDirs(directory, "rw"):
			self.directory = directory
		else:
			raise PermissionError(f"Missing read and/or write permissions in directory: {directory}")
		self.jobs = []
		
		self._queueConnection = DatabaseThread(self.directory / self.database)
		self._queueConnection.execute("CREATE TABLE IF NOT EXISTS queueTable (name TEXT UNIQUE, progress DECIMAL DEFAULT -1.0, modified INTEGER DEFAULT (UNIXEPOCH()));")
		if logger:
			self.LOG = logger
		if reportHook:
			self.reportHook = reportHook
		if hooks:
			self.hooks = hooks
		self._threads = threads if threads is not None else []
		self.LOG.info(f"Created {type(self)} object at 0x{id(self):0>16X}")

	def __del__(self):
		del self._queueConnection

	def addSources(self, *sources : str):
		self.SOURCES = self.SOURCES + sources

	def wait(self):
		threads = self._threads.copy()
		for t in threads:
			try:
				t.join()
				self._threads.remove(t)
			except:
				pass
			if hasattr(t, "exception"):
				raise t.exception
	
	@threadDescriptor
	def download(self, query, filename : str, reportHook : "DownloaderReportHook"=None) -> None:
		
		if reportHook:
			pass
		elif self.reportHook:
			reportHook = self.reportHook
		else:
			reportHook = DownloaderReportHook(getattr(self, "__name__", getattr(self.__class__, "__name__")), self.hooks, filename)
		
		job = Job(query, filename, conn=self._queueConnection, reportHook=reportHook, out=self.directory, logger=self.LOG)
		self.jobs.append(job)

		if job.reserveQueue():
			job.run(self.SOURCES, postProcess=self.postProcess)
		elif job.isDone():
			reportHook(int(1))
		else:
			while not job.reserveQueue():
				job.updateLoop(timeStep=self.timeStep)
				if job.isDone():
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
		re.compile(r"(?P<genome_id>(?P<idHead>\w{3})_(?P<id1_3>\d{3})(?P<id4_6>\d{3})(?P<id7_9>\d{3}).*)"),
		re.compile(r"(?P<assembly>.*)")
	]

class DatabaseDownloader(Downloader):
	SOURCES = [
		("MetaCanSNPer-data", "https://github.com/XJ-04561/MetaCanSNPer-data/raw/master/database/{query}"), # MetaCanSNPer
		("CanSNPer2-data", "https://github.com/FOI-Bioinformatics/CanSNPer2-data/raw/master/database/{query}") # Legacy CanSNPer
	]
	LOG = LOGGER.getChild("DatabaseDownloader")
	postProcess = staticmethod(correctDatabase)

class ReferenceDownloader(Downloader):
	SOURCES = [
		("NCBI", NCBI_URL)
	]
	LOG = LOGGER.getChild("ReferenceDownloader")
	postProcess = staticmethod(gunzip)