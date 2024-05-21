
from MetaCanSNPer.core.DirectoryLibrary import *


def test_init():
	
	DL = DirectoryLibrary("francisella_tularensis", ["FSC458.fq.gz"])

	assert DL.organism == "francisella_tularensis"
	assert DL.query == ("./FSC458.fq.gz",)
	assert DL.queryName == "FSC458"