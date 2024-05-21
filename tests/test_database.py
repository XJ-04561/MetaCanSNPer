
from MetaCanSNPer import *
from MetaCanSNPer.CommandLineParser import initializeMainObject


def test_init():
	
	class NameSpace:
		settingsFile = None
		silent = True
		query = ["FSC458.fq"]
		organism = "francisella_tularensis"
		database = None
		sessionName = None

	mObj : MetaCanSNPer = initializeMainObject(NameSpace)
	assert mObj.databasePath is not None