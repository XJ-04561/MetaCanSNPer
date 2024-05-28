
from MetaCanSNPer.core import MetaCanSNPer
from MetaCanSNPer.CLI import initializeMainObject, NameSpace, parser, runJob, handleOptions, saveResults, separateCommands
import pytest

@pytest.mark.skip
def test_init():

	filenames = "\n".join([
		"RAW_SEQUENCE_DATAFILE_R1.fq",
		"RAW_SEQUENCE_DATAFILE_R2.fq"
	])

	argString = f"""MetaCanSNPer
	--query
	{filenames}
	--organism francisella_tularensis
	--dry-run
	--saveTemp
	--debug
	--mapper
	minimap2
	--snpCaller
	gatk_Mutect2
	--mapperOptions
	-x
	sr"""

	argv = list(filter(None, argString.split()))
	
	argsDict = separateCommands(argv)

	args = parser.parse_args(argsDict["args"], namespace=NameSpace)

	handleOptions(args)

	mObj : MetaCanSNPer = initializeMainObject(args)
	assert mObj.databasePath is not None

	runJob(mObj, args, argsDict)

	saveResults(mObj, args)

	assert mObj.SNPresults or mObj.SNPresults == {}



if __name__ == "__main__":
	test_init()