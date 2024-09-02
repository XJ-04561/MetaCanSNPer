
from MetaCanSNPer.core.DirectoryLibrary import *
import os

def test_init():

	os.makedirs(os.path.splitext(__file__)[0], exist_ok=True)
	os.chdir(os.path.splitext(__file__)[0])
	
	DL = DirectoryLibrary("francisella_tularensis", ["FSC458.fq.gz"])

	assert DL.organism == "francisella_tularensis"
	assert DL.query == (os.path.join(".", "FSC458.fq.gz"),)
	assert DL.queryName == "FSC458"