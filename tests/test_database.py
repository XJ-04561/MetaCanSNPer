
from MetaCanSNPer import *
from MetaCanSNPer.CommandLineParser import initializeMainObject
from MetaCanSNPer.core.DirectoryLibrary import DirectoryLibrary
from MetaCanSNPer.modules.Downloader import DatabaseDownloader, DownloaderReportHook
from MetaCanSNPer.modules.Database import MetaCanSNPerDatabase

def test_database():
	
	from MetaCanSNPer.modules.Database import ReferencesTable, Genome, GenomeID

	assert ReferencesTable.Genome == Genome
	assert ReferencesTable.GenomeID == GenomeID

# import pytest
# @pytest.mark.skip
def test_database_download():

	organism = "francisella_tularensis"
	Lib = DirectoryLibrary(organism, ["FSC458.fq"])
	databaseName = f"{organism}.db"
	hooks = Lib.hooks

	from urllib.request import urlretrieve, HTTPError
	for sourceName, sourceLink in DatabaseDownloader.SOURCES:
		try:
			(outFile, msg) = urlretrieve(sourceLink.format(query=databaseName), filename=Lib.databaseDir.writable / databaseName) # Throws error if 404
			print(f"Success from {sourceName}!", outFile)
			break
		except Exception as e:
			pass
	else:
		raise FileNotFoundError()

	for directory in Lib.databaseDir:
		if databaseName in directory:
			if MetaCanSNPerDatabase(directory / databaseName, "r", organism=organism).valid is True:
				hooks.trigger("DatabaseDownloaderProgress", {"name" : databaseName, "progress" : int(1)})
				databasePath = directory / databaseName
				break
	else:
		RH = DownloaderReportHook("DatabaseDownloader", hooks, databaseName)
		
		databasePath = Lib.databaseDir.writable / databaseName
		DD = DatabaseDownloader(Lib.databaseDir.writable)

		DD.download(databaseName, databaseName, reportHook=RH)
		
		DD.wait()
	
	database = MetaCanSNPerDatabase(databasePath, "r", organism=organism)
	
	database.clearIndexes()

	if not database.valid:
		database.close()
		database = MetaCanSNPerDatabase(databasePath, "w", organism=organism)
		database.fix()
		database.close()
		database = MetaCanSNPerDatabase(databasePath, "r", organism=organism)

	print(database.exception)

	assert database.valid

# def test_init():
	
# 	class NameSpace:
# 		settingsFile = None
# 		silent = True
# 		query = ["FSC458.fq"]
# 		organism = "francisella_tularensis"
# 		database = None
# 		sessionName = None

# 	mObj : MetaCanSNPer = initializeMainObject(NameSpace)
# 	assert mObj.databasePath is not None

if __name__ == "__main__":
	test_database_download()