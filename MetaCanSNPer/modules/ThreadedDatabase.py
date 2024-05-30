
from threading import Thread, Lock
import sqlite3
from queue import Queue, Empty as EmptyQueueException
from MetaCanSNPer.Globals import Logged

class CursorLike:
	def __init__(self, data : list):
		self.data = data
		self.dataIterator = iter(data)
	def __iter__(self):
		for row in self.data:
			yield row
	def __next__(self):
		return next(self.dataIterator)
	def fetchone(self):
		if self.data:
			return self.data[0]
		else:
			return None
	def fetchall(self):
		return self.data

class DatabaseThread(Logged):

	OPEN_DATABASES = {}
	queue : Queue[list[str,list,Lock, list]]
	queueLock : Lock
	running : bool
	
	filename : str
	_thread : Thread
	@property
	def _connection(self) -> "DatabaseThread":
		return self

	def __new__(cls, filename : str, factory=sqlite3.Connection):
		if (filename, factory) in cls.OPEN_DATABASES and cls.OPEN_DATABASES[filename, factory]._thread.is_alive():
			return cls.OPEN_DATABASES[filename, factory]
		else:
			cls.OPEN_DATABASES[filename, factory] = super().__new__(cls)
			return cls.OPEN_DATABASES[filename, factory]

	def __init__(self, filename : str, factory=sqlite3.Connection):
		self.running = True
		self.queue = Queue()
		self.queueLock = Lock()
		self.queueLock.acquire()
		self.filename = filename
		self._factory = factory
		self._thread = Thread(target=self.mainLoop, daemon=True)
		self._thread.start()

	def mainLoop(self):
		try:
			_connection = sqlite3.connect(self.filename, factory=self._factory)
			while self.running:
				try:
					string, params, lock, results = self.queue.get(timeout=15)
					try:
						results.extend(_connection.execute(string, params).fetchall())
					except Exception as e:
						self.LOG.exception(e)
						results.append(e)
					try:
						lock.release()
					except:
						pass
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
		_connection.close()

	def execute(self, string : str, params : list=[]):
		lock = Lock()
		lock.acquire()
		results = []
		self.queue.put([string, params, lock, results])
		lock.acquire()
		
		if results and isinstance(results[-1], Exception):
			raise results[-1]
		
		return CursorLike(results)
	
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