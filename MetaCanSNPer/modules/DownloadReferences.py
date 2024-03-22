
import os, textwrap, gzip, sys
from threading import Thread, Lock, ThreadError, Semaphore
from urllib.request import urlretrieve
from typing import Callable, Any
from subprocess import Popen
from time import sleep, time
import random
from PseudoPathy import DirectoryPath

import MetaCanSNPer.modules.LogKeeper as LogKeeper
from MetaCanSNPer.modules.Hooks import Hooks
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

class WorkerQueue:
	""""""

	queue : list[int]
	active : int
	jobs : dict[int,Job]
	finished : set[int]
	worker : Thread
	lock : Lock
	RUNNING : bool
	created : list[int]

	def __init__(self, hooks : Hooks=Hooks()):
		
		self.hooks = hooks
		self.semaphore = Semaphore()
		self.created = []
		self.queue = []
		self.jobs = {}
		self.finished = set()
		self.RUNNING = True
		self.lock = Lock()
		self.worker = Thread(target=self.mainLoop, daemon=True)
		self.worker.start()

		self.hooks.trigger("downloadFinished", target=self.releaseLock)
		self.hooks.trigger("downloadCrashed", target=lambda eventInfo : self.start())

	def __hash__(self):
		return self.id

	def __del__(self):
		self.RUNNING = False

	def releaseLock(self, eventInfo):
		if eventInfo["job"].id in self.created:
			self.semaphore.release()

	def start(self):
		self.worker = Thread(target=self.mainLoop, daemon=True)
		self.worker.start()
	
	def mainLoop(self):
		try:
			while self.RUNNING:
				if len(self.queue) == 0:
					sleep(1)
				else:
					self.active = self.queue.pop(0)
					job = self.jobs[self.active]
					job.target(*job.args, **job.kwargs)
					self.hooks.trigger("downloadFinished", {"job" : job})
		except Exception as e:
			LOGGER.exception(e)
			self.hooks.trigger("downloaderCrashed", {"object" : self})

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
		self.created.append(job.id)
		self.semaphore.acquire(False)
		self.queue.append(job.id)

		return job.id


	def wait(self, timeout=None):
		
		for i, jobID in enumerate(self.created):
			while jobID not in self.finished:
				if self.semaphore.acquire(timeout=timeout):
					pass
				else:
					if not self.worker.is_alive():
						self.created = self.created[i:]
						raise False
		self.created = []
		raise True
	
	def __iter__(self, timeout=None):

		for i, jobID in enumerate(self.created):
			while jobID not in self.finished:
				if self.semaphore.acquire(timeout=timeout):
					yield jobID
				else:
					if not self.worker.is_alive():
						self.created = self.created[i:]
						raise StopIteration
		self.created = []
		raise StopIteration
		

class DownloadQueue(WorkerQueue):
	""""""

	queue : list[int] = []
	jobs : dict[int,Job] = {}
	finished : set[int] = set()
	worker : Thread = None
	lock : Lock
	RUNNING : bool = Globals.RUNNING

	created : list

	def __init__(self, hooks : Hooks=Hooks()):
		
		self.hooks = hooks
		self.semaphore = Semaphore()
		self.created = []
		try:
			self.worker.start()
		except RuntimeError:
			pass # Already running

		self.hooks.trigger("downloadFinished", target=self.releaseLock)
		self.hooks.trigger("downloadCrashed", target=lambda eventInfo : self.restart())
	
	def __del__(self):
		for id in self.created:
			del self.jobs[id]

	def releaseLock(self, eventInfo):
		if eventInfo["job"].id in self.created:
			self.semaphore.release()

	def download(self, genbank_id : str, refseq_id : str, assembly_name : str, dst : DirectoryPath, filename : str=None, source : str="genbank", force=False, stdout=sys.stdout):
		"""Returns ID of job. Returns -1 if the file to be downloaded already exists and force==False."""
		filename = dst > (filename or f"{assembly_name}.fna")
		if pExists(filename) and not force:
			LOGGER.debug(f"File: {filename!r} already exists, not queueing for download.")
			return -1
		try:
			n1,n2,n3 = textwrap.wrap(genbank_id.split("_")[-1].split(".")[0],3)  ## Get the three number parts
		except ValueError:
			LOGGER.warning(f"Could not download {genbank_id}")
			return
		
		link = NCBI_FTP_LINK.format(source=SOURCED[source], n1=n1, n2=n2, n3=n3, genome_id=genbank_id, assembly=assembly_name)
		job = Job(target=DownloadQueue._download, kwargs={"filename":filename, "link":link, "force":force, "stdout":stdout})
		id = self.push(job)
		return id

	@staticmethod
	def _download(link : str=None, filename : DirectoryPath=None, force=False, stdout=sys.stdout) -> str | None:
		'''Download genomes from refseq or genbank on request. Default kwarg of force=False makes the download not take place if file already exists.
			ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/00x/xxx/GCF_00000xxxx.x_ASMxxxv1/GCF_00000xxxx.x_ASMxxxv1_genomic.fna.gz
		'''
		print(f"Downloading {filename!r} ... ", end="", flush=True, file=stdout)

		if os.path.exists(f"{filename}.gz") and not force:
			LOGGER.debug("Unzipping Reference file for {f}.".format(f=os.path.basename(filename).strip("_genomic.fna.gz")))
			DownloadQueue.gunzip(f"{filename}.gz")

		elif os.path.exists(filename) and not force:
			LOGGER.debug("Reference file for {f} already exists!".format(f=filename.split("/")[-1].strip("_genomic.fna.gz")))

		elif not os.path.exists(filename) or force:
			LOGGER.debug(f"Downloading: {link} <-> {filename}")
			urlretrieve(link, f"{filename}.gz")

			LOGGER.debug("Unzipping Reference file: '{f}'".format(f=os.path.basename(filename).strip("_genomic.fna.gz")))
			DownloadQueue.gunzip(f"{filename}.gz")
		
		print("Done!", flush=True, file=stdout)

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
	
	@staticmethod
	def restart():
		DownloadQueue.RUNNING = False
		try:
			DownloadQueue.worker.join(timeout=7)
			if DownloadQueue.worker.is_alive():
				raise ThreadError("Could not restart the DownloadQueue due to previous worker not finishing.")
		except AttributeError:
			DownloadQueue.worker = Thread(target=WorkerQueue.mainLoop, args=[DownloadQueue], daemon=True)
		except RuntimeError:
			del DownloadQueue.worker
		
		DownloadQueue.RUNNING = Globals.RUNNING
		DownloadQueue.worker = Thread(target=WorkerQueue.mainLoop, args=[DownloadQueue], daemon=True)
		DownloadQueue.worker.start()

# Start worker thread
DownloadQueue.worker = Thread(target=WorkerQueue.mainLoop, args=[DownloadQueue], daemon=True)