
import textwrap, gzip
from threading import Thread, Lock, ThreadError, Semaphore, Condition
from urllib.request import urlretrieve
from typing import Callable, Any
from time import sleep

import MetaCanSNPer.core.LogKeeper as LogKeeper
from MetaCanSNPer.core.Hooks import Hooks, DummyHooks, urlretrieveReportHook
from MetaCanSNPer.Globals import *
import MetaCanSNPer.Globals as Globals

LOGGER = LogKeeper.createLogger(__name__)

class DownloadFailed(Exception):
	pass

class Job:
	"""Each query into the DownloadQueue is represented by a "job"-container.
	This container contains everything needed for the queue to work, but also
	for the dispatcher to send information back to the queuer."""

	id : int
	data : dict
	target : Callable
	args : list
	kwargs : dict[str,Any]
	
	def __init__(self, id=None, target=lambda : None, args : list=[], kwargs : dict={}, data : dict={}):
		self.id = id
		self.target = target
		self.args = args
		self.kwargs = kwargs
		self.data = data
	
	def __getattr__(self, key):
		return self.data[key]
	
	def get(self, key):
		return self.__dict__.get(key) or self.data.get(key)

THREAD_COUNT_LIMIT = 10

class DummyThread:
	def is_alive(self):
		return False

class WorkerQueue:
	""""""

	queue : list[int]
	active : int
	jobs : dict[int,Job]
	finished : set[int]
	manager : Thread
	hooks : Hooks
	RUNNING : bool
	workers : list[Thread]

	def __init__(self, hooks : Hooks=None):
		
		self.hooks = hooks or Hooks()
		self.semaphore = Semaphore(THREAD_COUNT_LIMIT)
		self.queue = []
		self.jobs = {}
		self.finished = set()
		self.workers = []
		self.RUNNING = True

	def __hash__(self):
		return self.id

	def __del__(self):
		self.RUNNING = False

	def start(self):
		self.manager = Thread(target=self.mainLoop, daemon=True)
		self.manager.start()
	
	def mainLoop(self):
		try:
			while self.RUNNING:
				if self.semaphore.acquire(timeout=2):
					self.active = self.queue.pop(0)
					job = self.jobs[self.active]
					self.workers = list(filter(Thread.is_alive, self.workers))
					try:
						self.workers.append( thread := Thread(target=job.target, args=job.args, kwargs=job.kwargs|{"data":job.data}))
						thread.start()
					except:
						pass
					finally:
						self.finished.add(job.id)
		except Exception as e:
			LOGGER.exception(e)
			self.hooks.trigger("DownloaderCrashed", {"object" : self})

	def push(self, job : Job, id : int=None) -> int:
		if job.id is not None:
			pass
		elif id is not None:
			job.id = id
		else:
			n = 0
			while (id := int("{}{:0>10}".format(*random.random().as_integer_ratio()))) in self.jobs:
				if n > 10000: raise StopIteration("Couldn't find a unique id for job.")
				n += 1
			job.id = id
		
		self.jobs[job.id] = job
		self.queue.append(job.id)
		self.semaphore.release()

		return job.id
	
	def restart(self):
		self.RUNNING = False
		try:
			self.manager.join(timeout=7)
			if self.manager.is_alive():
				raise ThreadError("Could not restart the DownloadQueue due to previous worker not finishing.")
			del self.manager
		except AttributeError:
			pass
		except RuntimeError:
			del self.manager
		finally:
			self.RUNNING = Globals.RUNNING
			self.manager = Thread(target=self.mainLoop, daemon=True)
			self.manager.start()
		

class DownloadQueue(WorkerQueue):
	""""""

	worker : WorkerQueue = WorkerQueue()
	created : list

	def __init__(self, hooks : Hooks=None):
		
		self.hooks = hooks or Hooks()
		self.semaphore = Semaphore(0)
		self.created = []
		try:
			self.worker.start()
		except RuntimeError:
			pass # Already running

		self.hooks.addHook("DownloadReferencesFinished", target=self.finishedCallback)
		self.hooks.addHook("DownloadReferencesFailed", target=self.crashedCallback)
	
	def __del__(self):
		for id in self.created:
			try:
				del self.finished[id]
			except:
				pass

	def __enter__(self):
		return self

	def __exit__(self, *args):
		for id in self.created:
			try:
				del self.worker.finished[id]
			except:
				pass

	def __iter__(self):

		for i in range(len(self.created)):
			self.semaphore.acquire()
			for jobID in self.created:
				if jobID in self.worker.finished:
					break
			yield jobID
			self.created.remove(jobID)

		return

	@property
	def jobs(self):
		return self.worker.jobs

	def finishedCallback(self, eventInfo):
		if eventInfo["name"] in map(lambda jobID : self.worker.jobs[jobID].data["name"], self.created):
			self.semaphore.release()
	
	def crashedCallback(self, eventInfo):
		if eventInfo["name"] in map(lambda jobID : self.worker.jobs[jobID].data["name"], self.created):
			self.semaphore.release()

	def download(self, genbank_id : str, refseq_id : str, assembly_name : str, dst : DirectoryPath, filename : str=None, source : str="genbank", force=False, hooks : Hooks=DummyHooks()):
		"""Returns ID of job. Returns -1 if the file to be downloaded already exists and force==False."""
		filename = filename or f"{assembly_name}.fna"
		if dst.find(filename, purpose="r") is not None and not force:
			LOGGER.debug(f"File: {filename!r} already exists, not queueing for download.")
			self.hooks.trigger("DownloadReferencesSkipped", {"name" : assembly_name})
			return -1
		try:
			n1,n2,n3 = textwrap.wrap(genbank_id.split("_")[-1].split(".")[0],3)  ## Get the three number parts
		except ValueError:
			LOGGER.warning(f"Could not download {genbank_id}")
			self.hooks.trigger("DownloadReferencesFailed", {"name" : assembly_name})
			return
		
		link = NCBI_FTP_LINK.format(source=SOURCED[source], n1=n1, n2=n2, n3=n3, genome_id=genbank_id, assembly=assembly_name)
		job = Job(target=DownloadQueue._download, kwargs={"filename":dst.writable / filename, "link":link, "force":force, "hooks":hooks}, data={"name":assembly_name})
	
		self.created.append(self.worker.push(job))
		return self.created[-1]

	@staticmethod
	def _download(link : str=None, filename : DirectoryPath=None, force=False, hooks : Hooks=DummyHooks(), data:dict={}) -> str | None:
		'''Download genomes from refseq or genbank on request. Default kwarg of force=False makes the download not take place if file already exists.
			ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/00x/xxx/GCF_00000xxxx.x_ASMxxxv1/GCF_00000xxxx.x_ASMxxxv1_genomic.fna.gz
		'''
		try:
			hooks.trigger("DownloadReferencesStarting", {"name" : data["name"]})
			if force or not os.path.exists(filename) and not os.path.exists(f"{filename}.gz"):
				LOGGER.debug(f"Downloading: {link} -> {filename}")
				reportHook = urlretrieveReportHook("DownloadReferences", hooks, name=data["name"], steps=100)
				urlretrieve(link, f"{filename}.gz", reporthook=lambda *args : reportHook.send(args))
				reportHook.close()
			
			if os.path.exists(f"{filename}.gz") or force:
				hooks.trigger("DownloadReferencesPostProcess", {"name" : data["name"]})
				LOGGER.debug("Unzipping Reference file for {f}.".format(f=os.path.basename(filename).strip("_genomic.fna.gz")))
				DownloadQueue.gunzip(f"{filename}.gz")
				try:
					os.remove(f"{filename}.gz")
				except:
					pass
			hooks.trigger("DownloadReferencesFinished", {"name" : data["name"]})
			LOGGER.debug("Reference file for {f} already exists!".format(f=filename.split("/")[-1].strip("_genomic.fna.gz")))
		except Exception as e:
			LOGGER.exception(e)
			hooks.trigger("DownloadReferencesFailed", {"name" : data["name"]})
		
	@staticmethod
	def gunzip(filename : str, dst : str=None, wait=True):
		'''gunzip's given file. Only necessary for software that requires non-zipped data.'''
		if dst is not None:
			outFile = open(dst, "wb")
		else:
			outFile = open(filename.rstrip(".gz"), "wb")

		if wait is True:
			outFile.write(gzip.decompress(open(filename, "rb").read()))
			outFile.close()
		else:
			t = Thread( target=lambda f : (f.write(gzip.decompress(open(filename, "rb").read())), f.close()), args=[outFile], daemon=True)
			t.start()
			return t
		return None
