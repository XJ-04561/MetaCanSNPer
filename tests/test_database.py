
from MetaCanSNPer.Globals import *
from MetaCanSNPer.CommandLineParser import initializeMainObject
from MetaCanSNPer.core.DirectoryLibrary import DirectoryLibrary
from MetaCanSNPer.modules.Downloader import DatabaseDownloader, DownloaderReportHook
from MetaCanSNPer.modules.Database import MetaCanSNPerDatabase

def test_database():
	
	from MetaCanSNPer.modules.Database import ReferencesTable, Genome, GenomeID

	assert ReferencesTable.Genome == Genome
	assert ReferencesTable.GenomeID == GenomeID

import pytest
@pytest.mark.skip
def test_database_download():

	organism = "francisella_tularensis"
	Lib = DirectoryLibrary(organism, ["FSC458.fq"])
	databaseName = f"{organism}.db"
	hooks = Lib.hooks

	for directory in Lib.databaseDir:
		print(f"Checking {directory / databaseName}")
		if os.path.exists(directory / databaseName):
			print("Exists!")
			if MetaCanSNPerDatabase(directory / databaseName, "r", organism=organism).valid:
				print("Is valid!")
				hooks.trigger("DatabaseDownloaderProgress", {"name" : databaseName, "value" : int(1)})
				databasePath = directory / databaseName
				break
	else:
		print("Me is here!")
		outDir = Lib.databaseDir.writable
		databasePath = outDir / databaseName
		DD = DatabaseDownloader(outDir, hooks=hooks)

		ret = DD.download(databaseName, databaseName)
		print(f"This is the thread! {ret} {ret.is_alive()=}")
		print(f"This is the Database Thread! {DD._queueConnection._thread.is_alive()=}")

		DD.wait()
		print(f"Is the thread still alive? {ret} {ret.is_alive()=}")
		import time
		time.sleep(0.5)
		print(f"Is the thread *STILL* alive? {ret} {ret.is_alive()=}")
		print(f"Is the Database Thread still alive?! {DD._queueConnection._thread.is_alive()=}")
		print(DD._threads)
		from pprint import pprint
		pprint(dict(map(lambda x:(x, x.is_alive()), DD._threads)))
	
	database = MetaCanSNPerDatabase(databasePath, "r", organism=organism)

	if not database.valid:
		database.reopen("w")
		database.fix()
		database.reopen("r")

	print(database.exception)

	assert database.valid

# if __name__ == "__main__":
# 	test_database_download()
