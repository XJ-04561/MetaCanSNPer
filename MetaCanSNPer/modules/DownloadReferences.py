
import os, textwrap
from threading import Thread, Condition
from urllib.request import urlretrieve
from typing import Callable
from subprocess import Popen

from MetaCanSNPer.modules.LogKeeper import createLogger
from MetaCanSNPer.Globals import *

LOGGER = createLogger(__name__)

class _DownloadQueue:
	_queue : list[dict[str,Callable|list|dict]] = []
	_history : list[dict[str,Callable|list|dict]] = []
	worker : Thread
	running : bool
	def __init__(self):
		self._queue = []
		self._history = []
		self.worker = Thread(target=self.mainLoop)

	def mainLoop(self):
		while self.running:
			Condition.wait_for(predicate=lambda : len(self._queue) > 0, timeout=3)
			if len(self._queue) == 0: continue
			
			task = self._queue.pop(0)
			if task in self._history:
				continue
			else:
				self._history.append(task)
			task["target"](*task["args"], **task["kwargs"])

	def __del__(self):
		self.running = False

	def download(self, genbank_id : str, refseq_id : str, assembly_name : str, dst : str, source : str="genbank", force=False):
		self._queue.append({"target":_DownloadQueue.download, "args":[genbank_id, refseq_id, assembly_name, dst], "kwargs":{"source":source, "force":force}})

	@staticmethod
	def _download(genbank_id : str, refseq_id : str, assembly_name : str, dst : str, filename : str=None, source : str="genbank", force=False) -> str | None:
		'''Download genomes from refseq or genbank on request. Default kwarg of force=False makes the download not take place if file already exists. Does not unzip gzipped downloads.
			ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/00x/xxx/GCF_00000xxxx.x_ASMxxxv1/GCF_00000xxxx.x_ASMxxxv1_genomic.fna.gz
		'''
		retrieved_file = os.path.join(f"{dst}", filename or f"{assembly_name}.fna")
		if os.path.exists(retrieved_file+".gz") and not force:
			LOGGER.debug("Unzipping Reference file for {f}.".format(f=os.path.basename(retrieved_file).strip("_genomic.fna.gz")))
			_DownloadQueue.gunzip(retrieved_file+".gz")
		elif os.path.exists(retrieved_file) and not force:
			try:
				n1,n2,n3 = textwrap.wrap(genbank_id.split("_")[-1].split(".")[0],3)  ## Get the three number parts
			except ValueError:
				LOGGER.warning("Could not download {refid}".format(refid=genbank_id))
				return
			link = NCBI_FTP_LINK.format(source=SOURCED[source], n1=n1, n2=n2, n3=n3, genome_id=genbank_id, assembly=assembly_name)
		
			LOGGER.debug("Downloading: "+ link +" <-> "+ retrieved_file)
			urlretrieve(link, retrieved_file+".gz")
			LOGGER.debug("Unzipping Reference file: '{f}'".format(f=os.path.basename(retrieved_file).strip("_genomic.fna.gz")))
			_DownloadQueue.gunzip(retrieved_file+".gz")
		else:
			LOGGER.debug("Reference file for {f} already exists!".format(f=retrieved_file.split("/")[-1].strip("_genomic.fna.gz")))
		
		return retrieved_file

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

'''
This should run on every import. It makes sure that a script only uses one object to do their downloads. Which
enables sequential queueing of downloads without overwriting each other.

Expected behavior is for downloads to be occurring in parallel, but that downloads of the same materials do not
occur more than once, even if force=True
'''

DownloadQueue = _DownloadQueue()