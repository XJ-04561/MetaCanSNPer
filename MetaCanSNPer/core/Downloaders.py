
from MetaCanSNPer.Globals import *
from GeekyGadgets.Downloader import Downloader
from GeekyGadgets.URL import URL_TEMPLATE, URL

def correctDatabase(filename, finalFilename):
	from MetaCanSNPer.core.Database import MetaCanSNPerDatabase, NoChromosomesInDatabase
	database = MetaCanSNPerDatabase(filename, "w")
	database.fix()
	if not database.valid:
		e = database.exception
		if not isinstance(e, NoChromosomesInDatabase):
			e.add_note(f"Tables hash: {database.tablesHash}\nIndexes Hash: {database.indexesHash}")
			raise e
	database.close()
	if filename != finalFilename:
		os.rename(filename, finalFilename)

def gunzip(filepath : str, outfile : str):
	'''gunzip's given file. Only necessary for software that requires non-zipped data.'''
	import gzip
	
	with gzip.open(filepath, "rb") as zipped:
		rawfile = open(outfile, "wb")
		for row in zipped:
			rawfile.write(row)
		rawfile.close()



class NCBI_URL(URL_TEMPLATE):
	url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{idHead}/{id1_3}/{id4_6}/{id7_9}/{genome_id}_{assembly}/{genome_id}_{assembly}_genomic.fna.gz"
	formatters = (
		re.compile(r"(?P<genome_id>(?P<idHead>\w{3})_(?P<id1_3>\d{3})(?P<id4_6>\d{3})(?P<id7_9>\d{3}).*)")
	)
	@classmethod
	def complexFormat(cls, query) -> URL:
		genomeID, assemblyName = query
		return cls.url.format(**cls.formatters[0].match(genomeID).groupdict(), assembly=assemblyName)

class FOI_METACANSNPER_URL(URL_TEMPLATE):
	url = "https://github.com/FOI-Bioinformatics/MetaCanSNPer-data/raw/master/database/{query}"

class FOI_CANSNPER2_URL(URL_TEMPLATE):
	url = "https://github.com/FOI-Bioinformatics/CanSNPer2-data/raw/master/database/{query}"



class DatabaseDownloader(Downloader):
	SOURCES = [
		("MetaCanSNPer-data", FOI_METACANSNPER_URL), # MetaCanSNPer
		("CanSNPer2-data", FOI_CANSNPER2_URL) # Legacy CanSNPer
	]
	postProcess = staticmethod(correctDatabase)

class ReferenceDownloader(Downloader):
	SOURCES = [
		("NCBI", NCBI_URL)
	]
	postProcess = staticmethod(gunzip)