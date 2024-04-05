#!/usr/bin/env python3

from timeit import default_timer as timer
startTime = timer()
import logging, sys, argparse, traceback
from threading import Thread
from time import sleep
from typing import Callable
import re

## import MetaCanSNPer specific modules
from MetaCanSNPer.Globals import *
from MetaCanSNPer.Globals import __version__
import MetaCanSNPer.modules.LogKeeper as LogKeeper
from MetaCanSNPer.modules.MetaCanSNPer import MetaCanSNPer
from MetaCanSNPer.modules.TerminalUpdater import TerminalUpdater, Spinner, LoadingBar, TextProgress
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

"MetaCanSNPer --query RAW_SEQUENCE_DATAFILE_R1.fq RAW_SEQUENCE_DATAFILE_R2.fq --database DATABASE_FILE.db --dry-run --debug --mapper minimap2 --snpCaller gatk_Mutect2 --mapperOptions -x sr"

def createParser():
	"""Initiate MetaCanSNPer command line argument parser"""
	parser = argparse.ArgumentParser(description="MetaCanSNPer", usage="""MetaCanSNPer --query RAW_SEQUENCE_DATAFILE.* [RAW_SEQUENCE_DATAFILE_2.*] \\
--database DATABASE_FILE.db --mapper MAPPER_COMMAND --snpCaller SNPCALLER_COMMAND \\
\t--mapperOptions [Flags as they would be passed to the mapper] \\
\t--snpCallerOptions [Flags as they would be passed to the snpCaller]
Examples:
MetaCanSNPer --query RAW_SEQUENCE_DATAFILE.fq --database DATABASE_FILE.db \\
\t--mapper minimap2 --snpCaller gatk_Mutect2 \\
\t--mapperOptions -x ava-ont
MetaCanSNPer --query RAW_SEQUENCE_DATAFILE_R1.fq RAW_SEQUENCE_DATAFILE_R2.fq \\
\t--database DATABASE_FILE.db \\
\t--mapper minimap2 --snpCaller gatk_Mutect2 \\
\t--mapperOptions -x sr
MetaCanSNPer --query SEQUENCE_ASSEMBLY.fna --database DATABASE_FILE.db \\
\t--aligner progressiveMauve --snpCaller ParseXMFA2
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
		optionalArguments.add_argument("--saveTemp",			action="store_true",					help="Don't dispose of temporary directories/files.")
		optionalArguments.add_argument("--settingsFile",		metavar="settingsFile",					help="Path to .TOML file containing settings for MetaCanSNPer. Check the 'defaultConfig.toml' to see what can be included in a settings file.")
		optionalArguments.add_argument("--downloadTimeout",		metavar="downloadTimeout", default=300,	help="Maximum time to take when downloading a reference genome file. If exceeded by the program it will throw an error and stop the current run. A value of -1 means no maximum time.")
		
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
		debugOptions.add_argument("--silent",	action="store_true",	help="Disables printing to terminal except for any error messages which might appear.")
		debugOptions.add_argument("--dry-run",	action="store_true",	help="Don't run the processes of the mapper/aligner/snpCaller, just run a randomised (1 - 5 sec) `sleep` call.")
	

	return parser

def handleOptions(args : argparse.Namespace):
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

	if args.dry_run:
		Globals.DRY_RUN = args.dry_run # Don't run the processes of the mapper/aligner/snpCaller, just run a randomised `sleep` call
	if args.saveTemp:
		Globals.PPGlobals.DISPOSE = False # Don't dispose of temporary directories/files.
	
	if args.debug:
		Globals.LOGGER_FILEHANDLER.setLevel(logging.DEBUG)
	elif args.verbose:
		Globals.LOGGER_FILEHANDLER.setLevel(logging.INFO)
	elif args.suppress:
		Globals.LOGGER_FILEHANDLER.setLevel(logging.ERROR)
	else:
		pass # The default logging level for the logging package is logging.WARNING
	
	Globals.DOWNLOAD_TIMEOUT = args.downloadTimeout

	if args.installDir is not None:
		import PseudoPathy.Globals
		PseudoPathy.Globals.PROGRAM_DIRECTORY = args.installDir

def initializeMainObject(args):
	if args.silent:
		TerminalUpdater = lambda *args, **kwargs: open(os.devnull, "r")

	mObj = MetaCanSNPer(settings=vars(args), settingsFile=args.settingsFile)
	
	if not args.silent: print(f"Checking query {args.query}")
	mObj.setQuery(args.query)

	with TerminalUpdater(f"Checking database {args.database}:", category="DownloadDatabase", hooks=mObj.hooks, threadNames=[os.path.basename(args.database)], printer=LoadingBar):
		mObj.setDatabase(args.database)


	if args.sessionName is not None: mObj.setSessionName(args.sessionName)

	with TerminalUpdater(f"Checking Reference Genomes:", category="DownloadReferences", hooks=mObj.hooks, threadNames=[mObj.database.ReferenceTable.get(DB.Assembly)], printer=LoadingBar):
		mObj.setReferenceFiles()
	
	return mObj

def runJob(mObj : MetaCanSNPer, args : argparse.Namespace, argsDict : dict):
	if args.silent:
		TerminalUpdater : TerminalUpdater = lambda *args, **kwargs: open(os.devnull, "r")
	
	genomes = mObj.database.ReferenceTable.all(DB.Genome)
	
	if args.mapper is not None:
		with TerminalUpdater(f"Creating Mappings:", category="Mappers", hooks=mObj.hooks, threadNames=[genomes], printer=Spinner):
			mObj.createMap(softwareName=args.mapper, flags=argsDict.get("--mapperOptions", {}))

	if args.aligner is not None:
		with TerminalUpdater(f"Creating Alignments:", category="Aligners", hooks=mObj.hooks, threadNames=[genomes], printer=Spinner):
			mObj.createAlignment(softwareName=args.aligner, flags=argsDict.get("--alignerOptions", {}))
	
	with TerminalUpdater(f"Creating Alignments:", category="Aligners", hooks=mObj.hooks, threadNames=[genomes], printer=Spinner):
		mObj.callSNPs(softwareName=args.snpCaller, flags=argsDict.get("--snpCallerOptions", {}))

def saveResults(mObj : MetaCanSNPer, args : argparse.Namespace):
	if not args.silent:
		print("Saving results...", end="", flush=True)

	mObj.saveSNPdata()
	mObj.saveResults()

	if not args.silent:
		print(" Done!", flush=True)

def main():
	
	argsDict = separateCommands(sys.argv)
	
	parser = createParser()

	if len(sys.argv)<2:
		parser.print_help()
		parser.exit()

	args : argparse.Namespace = parser.parse_args(argsDict["args"])


	handleOptions(args)
	
	try:
		mObj = initializeMainObject(args)

		runJob(mObj, args, argsDict)

		saveResults(mObj, args)

		print(f"{SOFTWARE_NAME} finished in {timer() - startTime:.3f} seconds! Results exported to: {mObj.Lib.resultDir}")
	except Exception as e:
		LOGGER.exception(e)
		for exc in mObj.exceptions:
			print(f"{type(exc).__name__}: {exc}", file=sys.stderr)

		if not args.silent: print(f"{SOFTWARE_NAME} ended before completing query. ", end="")
		print("Exception occured: \n", file=sys.stderr)

		if args.debug:
			pattern = re.compile(r"(\[[\w.]+\] \d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3} - ERROR: .*?)(?:\[[\w.]+\] \d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3} - \w+?:|\$)", flags=re.DOTALL+re.MULTILINE)
			for i, err in enumerate(pattern.finditer(open(LOGGER_FILEHANDLER.baseFilename, "r").read())):
				print(f"Exception [{i}]\n{err.group(1)}", file=sys.stderr)
		else:
			print(f"{type(e).__name__}: "+str(e), file=sys.stderr)
		exit(1)
