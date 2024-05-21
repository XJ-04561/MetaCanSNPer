
from MetaCanSNPer.core.DirectoryLibrary import DirectoryLibrary
from MetaCanSNPer.modules.Downloader import DatabaseDownloader, DownloaderReportHook
from MetaCanSNPer.modules.Database import MetaCanSNPerDatabase


def test_database_download():

	organism = "francisella_tularensis"
	Lib = DirectoryLibrary(organism, ["FSC458.fq"])
	databaseName = f"{organism}.db"
	hooks = Lib.hooks

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
	
	if not database.valid:
		database.close()
		database = MetaCanSNPerDatabase(databasePath, "w", organism=organism)
		database.fix()
		database.close()
		database = MetaCanSNPerDatabase(databasePath, "r", organism=organism)

	assert database.valid
	
if __name__ == "__main__":
	test_database_download()