
import os
import textwrap
from urllib.request import urlretrieve
import logging
try:
	from CanSNPer2.modules.LogKeeper import createLogger

	LOGGER = createLogger(__name__)
except:
	from LogKeeper import createLogger

	LOGGER = createLogger(__name__)
SOURCED = {"refseq":"F", "genbank": "A"}
REFERENCE_FILENAME_FORMAT = "{genome_id}_{assembly}_genomic.fna"
NCBI_FTP_LINK = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GC{source}/{n1}/{n2}/{n3}/{genome_id}_{assembly}/{genome_id}_{assembly}_genomic.fna.gz"

## import standard python libraries for subprocess and multiprocess
from subprocess import Popen,PIPE,STDOUT

def download(genbank_id : str, refseq_id : str, assembly_name : str, dst : str, source : str="genbank", force=False) -> str | None:
	'''Download genomes from refseq or genbank on request. Default kwarg of force=False makes the download not take place if file already exists. Does not unzip gzipped downloads.
		ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/00x/xxx/GCF_00000xxxx.x_ASMxxxv1/GCF_00000xxxx.x_ASMxxxv1_genomic.fna.gz
	'''
	retrieved_file = os.path.join("{directory}", (REFERENCE_FILENAME_FORMAT)).format(genome_id=genbank_id,assembly=assembly_name,directory=dst)
	if os.path.exists(retrieved_file+".gz") and not force:
		LOGGER.debug("Unzipping Reference file for {f}.".format(f=os.path.basename(retrieved_file).strip("_genomic.fna.gz")))
		gunzip(retrieved_file+".gz")
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
		gunzip(retrieved_file+".gz")
	else:
		LOGGER.debug("Reference file for {f} already exists!".format(f=retrieved_file.split("/")[-1].strip("_genomic.fna.gz")))
	
	return retrieved_file

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