
import os, textwrap
from threading import Thread, Lock
from urllib.request import urlretrieve
from typing import Callable, Any
from subprocess import Popen
from time import sleep, time
import random
from PseudoPathy import DirectoryPath

import MetaCanSNPer.modules.LogKeeper as LogKeeper
from MetaCanSNPer.Globals import *
import MetaCanSNPer.Globals as Globals

LOGGER = LogKeeper.createLogger(__name__)

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
	jobs : dict[int,Job]
	finished : set[int]
	worker : Thread
	lock : Lock
	RUNNING : bool

	def __init__(self):
		
		self.queue = []
		self.jobs = {}
		self.finished = {}
		self.running = True
		self.lock = Lock()
		self.worker = Thread(target=self.mainLoop, daemon=True)
		self.worker.start()

	def mainLoop(self):
		while self.RUNNING:
			if len(self.queue) == 0:
				sleep(1)
			else:
				taskID = self.queue.pop(0)
				task = self.jobs[taskID]
				task.target(*task.args, **task.kwargs)

				self.lock.release()
				self.lock.acquire(False)

	def __del__(self):
		self.running *= False
	
	def push(self, job : Job, id : int=None):
		if id is None and job.id is None:
			n = 0
			while (id := int("{}{:0>10}".format(*random.random().as_integer_ratio()))) in self.jobs:
				if n > 10000: raise StopIteration("Couldn't find a unique id for job.")
				n += 1
		self.jobs[id] = job
		self.queue.append(id)

	#
		# 	IM HERE
		#

	def waitFor(self, *args, timeout=300):
		querySet = set(args)
		N = len(querySet)
		
		for i in range(N):
			startOfJob = time()
			if i >= len(querySet):
				continue # Prevents more than one job being finished in one release-cycle and then waiting forever.

			LOGGER.debug(f"Waiting for download {1+N-len(querySet)}/{N}")
			finished : set
			while len(finished) == 0 and time() - startOfJob < timeout:
				self.lock.acquire(timeout=5)
				finished = querySet.intersection(self.finished)
			querySet.difference_update(finished)

			if len(finished) == 0:
				yield -1
				raise StopAsyncIteration(f"WorkerQueue.__iter__: Timeout(={timeout}) was reached.")
			for f in list(finished):
				yield f
	
	def __iter__(self, timeout=300):
		return self.iter(*self.queue, timeout=timeout)
		

class DownloadQueue: pass

class DownloadQueue:
	""""""

	queue : list[int] = []
	jobs : dict[Job] = {}
	finished : set[int] = set()
	worker : Thread = Thread(target=WorkerQueue.mainLoop, args=[DownloadQueue], daemon=True)
	lock : Lock = Lock()
	RUNNING : bool = Globals.RUNNING

	def __init__(self) -> None:
		try:
			self.worker.start()
		except RuntimeError:
			pass # Already running
	
	def __del__(self):
		pass

	def download(self, genbank_id : str, refseq_id : str, assembly_name : str, dst : DirectoryPath, filename : str=None, source : str="genbank", force=False):
		filename = dst > (filename or f"{assembly_name}.fna")
		job = Job()
		self.queue.append({"target":DownloadQueue._download, "args":[genbank_id, refseq_id, assembly_name], "kwargs":{"filename":filename, "source":source, "force":force}})
		return job

	@staticmethod
	def _download(genbank_id : str, refseq_id : str, assembly_name : str, filename : DirectoryPath, source : str="genbank", force=False) -> str | None:
		'''Download genomes from refseq or genbank on request. Default kwarg of force=False makes the download not take place if file already exists.
			ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/00x/xxx/GCF_00000xxxx.x_ASMxxxv1/GCF_00000xxxx.x_ASMxxxv1_genomic.fna.gz
		'''
		if os.path.exists(f"{filename}.gz") and not force:
			LOGGER.debug("Unzipping Reference file for {f}.".format(f=os.path.basename(filename).strip("_genomic.fna.gz")))
			DownloadQueue.gunzip(f"{filename}.gz")

		elif os.path.exists(filename) and not force:
			LOGGER.debug("Reference file for {f} already exists!".format(f=filename.split("/")[-1].strip("_genomic.fna.gz")))

		elif not os.path.exists(filename) or force:
			try:
				n1,n2,n3 = textwrap.wrap(genbank_id.split("_")[-1].split(".")[0],3)  ## Get the three number parts
			except ValueError:
				LOGGER.warning(f"Could not download {genbank_id}")
				return
			
			link = NCBI_FTP_LINK.format(source=SOURCED[source], n1=n1, n2=n2, n3=n3, genome_id=genbank_id, assembly=assembly_name)
			# link = f"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GC{SOURCED[source]}/{n1}/{n2}/{n3}/{genbank_id}_{assembly_name}/{genbank_id}_{assembly_name}_genomic.fna.gz"

			LOGGER.debug(f"Downloading: {link} <-> {filename}")
			urlretrieve(link, f"{filename}.gz")

			LOGGER.debug("Unzipping Reference file: '{f}'".format(f=os.path.basename(filename).strip("_genomic.fna.gz")))
			DownloadQueue.gunzip(f"{filename}.gz")

	@staticmethod
	def gunzip(filename, dst=None, wait=True):
		'''gunzip's given file. Only necessary for software that requires non-zipped data.'''
		p = Popen(["gunzip", "-f", filename])
		bnamePath, ext = os.path.splitext(filename)
		if dst is not None:
			try:
				os.rename(bnamePath, dst)
			except FileExistsError:
				LOGGER.debug("'{}' File already exists, skip!".format(dst))
		if wait is True:
			p.wait()
