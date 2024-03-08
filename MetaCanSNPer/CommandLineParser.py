#!/usr/bin/env python3

import logging, sys, argparse

## import MetaCanSNPer specific modules
from MetaCanSNPer.Globals import *
from MetaCanSNPer.Globals import __version__
import MetaCanSNPer.modules.LogKeeper as LogKeeper
from MetaCanSNPer.modules.MetaCanSNPer import MetaCanSNPer
import MetaCanSNPer.Globals as Globals

LOGGER = LogKeeper.createLogger(__name__)

def separateCommands(argv : list[str]) -> dict[str,list[str]]:
	order = [(0, "args")]
	try:
		order.append((argv.index("--mapperOptions"), "--mapperOptions"))
	except:
		pass
	try:
		order.append((argv.index("--alignerOptions"), "--alignerOptions"))
	except:
		pass
	try:
		order.append((argv.index("--snpCallerOptions"), "--snpCallerOptions"))
	except:
		pass
	order.append((len(argv), None))
	order.sort(key=lambda t : t[0])

	out = {}
	for i in range(len(order)-1):
		out[order[i][1]] = argv[order[i][0]+1:order[i+1][0]]
	
	return out

def createParser():
	"""Initiate MetaCanSNPer command line argument parser"""
	parser = argparse.ArgumentParser(description="MetaCanSNPer", usage="""MetaCanSNPer --query RAW_SEQUENCE_DATAFILE.* [RAW_SEQUENCE_DATAFILE_2.*] --database DATABASE_FILE.db --mapper MAPPER_COMMAND --snpCaller SNPCALLER_COMMAND \\
								  \t--mapperOptions [Flags as they would be passed to the mapper] \\
								  \t--snpCallerOptions [Flags as they would be passed to the snpCaller]
								  Examples:
								  MetaCanSNPer --query RAW_SEQUENCE_DATAFILE.fq --database DATABASE_FILE.db --mapper minimap2 --snpCaller gatk_Mutect2 \\
								  \t--mapperOptions -x ava-ont
								  MetaCanSNPer --query RAW_SEQUENCE_DATAFILE_R1.fq RAW_SEQUENCE_DATAFILE_R2.fq --database DATABASE_FILE.db --mapper minimap2 --snpCaller gatk_Mutect2 \\
								  \t--mapperOptions -x sr
								  MetaCanSNPer --query SEQUENCE_ASSEMBLY.fna --database DATABASE_FILE.db --aligner progressiveMauve --snpCaller ParseXMFA2
								  """)

	parser.add_argument("--version", action="store_true", help=argparse.SUPPRESS)
	parser.add_argument("--list", action="store_true", help="To list implemented software and exit.")

	# The 'if True:' structures are used to minimize and expand sections in an IDE.
	requiredArguments = parser.add_argument_group("Required arguments")
	if True:
		requiredArguments.add_argument("--query", nargs="+",	metavar="query",		help="Raw sequence data file supported by the intended Aligner/Mapper.")
		requiredArguments.add_argument("-d", "--database",		metavar="database",		help="Filename of CanSNP database to be used.")
		mapOrAlign = requiredArguments.add_mutually_exclusive_group(required=True)
		if True:
			mapOrAlign.add_argument("--mapper",					metavar="mapper",		help="Name of installed and supported Mapper software.")
			mapOrAlign.add_argument("--aligner",				metavar="aligner",		help="Name of installed and supported Alignment software.")
		requiredArguments.add_argument("--snpCaller",			metavar="snpCaller",	help="Name of installed and supported SNP Calling software.")

	optionalArguments = parser.add_argument_group("Optional arguments")
	if True:
		optionalArguments.add_argument("--saveTemp",			action="store_true",	help="Don't dispose of temporary directories/files.")
		optionalArguments.add_argument("--settingsFile",		metavar="settingsFile",	help="Path to .TOML file containing settings for MetaCanSNPer. Check the 'defaultConfig.toml' to see what can be included in a settings file.")
		
		# Not used by the argparser, but is used for the help-page and for splitting the argv
		mapperOptions = optionalArguments.add_argument("--mapperOptions",		metavar="Mapper options",		help=MAPPER_OPTIONS_EXPLAINER)
		alignerOptions = optionalArguments.add_argument("--alignerOptions",		metavar="Aligner options",		help=ALIGNER_OPTIONS_EXPLAINER)
		snpCallerOptions = optionalArguments.add_argument("--snpCallerOptions",	metavar="SNP Caller options",	help=SNP_CALLER_OPTIONS_EXPLAINER)

	directoryOptions = parser.add_argument_group("Directory Options")
	if True:
		directoryOptions.add_argument("-W", "--workDir",		metavar="DIR", default=None, help="Work directory")
		directoryOptions.add_argument("-U", "--userDir",		metavar="DIR", default=None, help="User directory")
		directoryOptions.add_argument("-I", "--installDir",		metavar="DIR", default=None, help="Installation directory")
		directoryOptions.add_argument("-Q", "--targetDir",		metavar="DIR", default=None, help="Target (Query) directory")
		directoryOptions.add_argument("-T", "--tmpDir",			metavar="DIR", default=None, help="Temporary directory")
		directoryOptions.add_argument("-R", "--refDir",			metavar="DIR", default=None, help="References directory")
		directoryOptions.add_argument("-D", "--databaseDir",	metavar="DIR", default=None, help="Databases directory")
		directoryOptions.add_argument("-O", "--outDir",			metavar="DIR", default=None, help="Output directory")
		directoryOptions.add_argument("-S", "--sessionName",	metavar="DIR", default=None, help="Session Name/Directory")

	debugOptions = parser.add_argument_group("Logging and debug options")
	if True:
		debugOptions.add_argument("--verbose",	action="store_true",	help="Verbose output")
		debugOptions.add_argument("--debug",	action="store_true",	help="Debug output")
		debugOptions.add_argument("--suppress",	action="store_true",	help="Suppress warnings")
		debugOptions.add_argument("--dry-run",	action="store_true",	help="Don't run the processes of the mapper/aligner/snpCaller, just run a randomised (1 - 5 sec) `sleep` call.")
	

	return parser

def main():
	
	argsDict = separateCommands(sys.argv)
	
	parser = createParser()

	if len(sys.argv)<2:
		parser.print_help()
		parser.exit()

	args = parser.parse_args(argsDict["args"])

	if args.version:
		print(f"MetaCanSNPer - version {__version__}")
		exit()

	elif args.list:
		from MetaCanSNPer.modules.Wrappers import Mapper, Aligner, SNPCaller
		print("\nMappers:")
		for mapper in Mapper.__subclasses__():			print(f"\t{mapper.softwareName}")
		print("\nAligners:")
		for aligner in Aligner.__subclasses__():		print(f"\t{aligner.softwareName}")
		print("\nSNPCallers:")
		for snpCaller in SNPCaller.__subclasses__():	print(f"\t{snpCaller.softwareName}")
		exit()

	if args["dry-run"]:
		Globals.DRY_RUN = args["dry-run"] # Don't run the processes of the mapper/aligner/snpCaller, just run a randomised `sleep` call
	if args.saveTemp:
		Globals.PPGlobals.DISPOSE = False # Don't dispose of temporary directories/files.
	
	if args.debug:
		Globals.LOGGER_FILEHANDLER.setLevel(logging.DEBUG)
	elif args.verbose:
		Globals.LOGGER_FILEHANDLER.setLevel(logging.INFO)
	elif args.supress:
		Globals.LOGGER_FILEHANDLER.setLevel(logging.ERROR)
	else:
		pass # The default logging level for the logging package is logging.WARNING

	mObj = MetaCanSNPer(settings=args, settingsFile=args["settingsFile"])

	mObj.setQuery(args.query)
	mObj.setDatabase(args.database)

	if args.sessionName is not None: mObj.setSessionName(args.sessionName)

	if "mapper" in args:
		mObj.createMap(softwareName=args.mapper, kwargs=argsDict["--mapperOptions"] if "--mapperOptions" in argsDict else {})
	if "aligner" in args:
		mObj.createAlignment(softwareName=args.aligner, kwargs=argsDict["--alignerOptions"] if "--alignerOptions" in argsDict else {})

	mObj.callSNPs(softwareName=args.snpCaller, kwargs=argsDict["--snpCallerOptions"] if "--snpCallerOptions" in argsDict else {})

	mObj.saveResults()
	mObj.saveSNPdata()

	print("Done!")

if oname=="__main__":
	main()
	
