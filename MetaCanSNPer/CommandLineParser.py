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
from MetaCanSNPer.core.MetaCanSNPer import MetaCanSNPer
from MetaCanSNPer.core.TerminalUpdater import TerminalUpdater, Spinner, LoadingBar, TextProgress
import MetaCanSNPer.Globals as Globals
import MetaCanSNPer as package

LOGGER = LOGGER.getChild(__name__.split(".")[-1])

class NameSpace(argparse.Namespace):

	listSoftware : bool
	version : bool

	query : list[str]
	organism : str
	database : str

	mapper : str
	aligner : str
	snpCaller : str
	
	settingsFile : str

	workDir : str
	userDir : str
	installDir : str
	targetDir : str
	tmpDir : str
	refDir : str
	databaseDir : str
	outDir : str

	sessionName : str
	
	saveTemp : bool
	debug : bool
	verbose : bool
	suppress : bool
	silent : bool
	dry_run : bool

	@overload
	def __init__(self,
				query : list[str],
				organism : str,
				database : str = None,
				mapper : str = None,
				aligner : str = None,
				snpCaller : str = None,
				settingsFile : str = None,
				workDir : str = None,
				userDir : str = None,
				installDir : str = None,
				targetDir : str = None,
				tmpDir : str = None,
				refDir : str = None,
				databaseDir : str = None,
				outDir : str = None,
				sessionName : str = None,
				listSoftware : bool = False,
				version : bool = False,
				saveTemp : bool = False,
				debug : bool = False,
				verbose : bool = False,
				suppress : bool = False,
				silent : bool = False,
				dry_run : bool = False
				): ...

	__init__ = argparse.Namespace.__init__

	def __getitem__(self, name):
		return getattr(self, name)
	def __setitem__(self, name, value):
		return setattr(self, name, value)
	def __iter__(self):
		return iter(self._get_kwargs)

parser = argparse.ArgumentParser(prog=__package__, description=package.__doc__, usage="""MetaCanSNPer --query RAW_SEQUENCE_DATAFILE.* [RAW_SEQUENCE_DATAFILE_2.*] \\
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
parser.add_argument("--list", dest="listSoftware", action="store_true", help="To list implemented software and exit.")

# The 'if True:' structures are used to minimize and expand sections in an IDE.
requiredArguments = parser.add_argument_group("Required arguments")
if True:
	requiredArguments.add_argument("--query", nargs="+",	metavar="query",		required=True, help="Raw sequence data file supported by the intended Aligner/Mapper.")
	requiredArguments.add_argument("--organism",			metavar="organism",		required=True, help="Name of organism queried. (Use \"_\" in place of spaces)")
	
	mapOrAlign = requiredArguments.add_mutually_exclusive_group(required=True)
	if True:
		mapOrAlign.add_argument("--mapper",					metavar="mapper",		help="Name of installed and supported Mapper software.")
		mapOrAlign.add_argument("--aligner",				metavar="aligner",		help="Name of installed and supported Alignment software.")
	requiredArguments.add_argument("--snpCaller",			metavar="snpCaller",	help="Name of installed and supported SNP Calling software.")

optionalArguments = parser.add_argument_group("Optional arguments")
if True:
	optionalArguments.add_argument("-d", "--database",	metavar="database",								help="Filename of CanSNP database to be used.")
	optionalArguments.add_argument("--saveTemp",		action="store_true",							help="Don't dispose of temporary directories/files.")
	optionalArguments.add_argument("--settingsFile",	metavar="settingsFile",							help="Path to .TOML file containing settings for MetaCanSNPer. Check the 'defaultConfig.toml' to see what can be included in a settings file.")

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

def separateCommands(argv : list[str]) -> dict[str,list[str]]:
	
	import itertools
	extraOptionFlags = ["args", "--mapperOptions", "--alignerOptions", "--snpCallerOptions"]
	workList = argv.copy()
	workList[0] = "args"
	proDict = []
	while workList:
		extraOptionFlags.remove(workList[0])
		mode, *args = itertools.takewhile(lambda x:x not in extraOptionFlags, workList)
		proDict.append((mode, args))
		workList = workList[len(args)+1:]
	
	return dict(proDict)

"MetaCanSNPer --query RAW_SEQUENCE_DATAFILE_R1.fq RAW_SEQUENCE_DATAFILE_R2.fq --database DATABASE_FILE.db --dry-run --debug --mapper minimap2 --snpCaller gatk_Mutect2 --mapperOptions -x sr"


def handleOptions(args : NameSpace):
	if args.version:
		print(f"MetaCanSNPer - version {__version__}")
		exit()

	elif args.listSoftware:
		from MetaCanSNPer.core.Wrappers import Mapper, Aligner, SNPCaller
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
		logging.basicConfig(level=logging.DEBUG)
	elif args.verbose:
		logging.basicConfig(level=logging.INFO)
	elif args.suppress:
		logging.basicConfig(level=logging.ERROR)
	else:
		pass # The default logging level for the logging package is logging.WARNING

def initializeMainObject(args : NameSpace) -> MetaCanSNPer:

	from MetaCanSNPer.modules.Database import ReferencesTable
	mObj = MetaCanSNPer(args.organism, args.query, settings=vars(args), settingsFile=args.settingsFile)

	database = args.database or mObj.databaseName
	
	with TerminalUpdater(f"Checking database {database}:", category="DatabaseDownloader", hooks=mObj.hooks, threadNames=[database], printer=LoadingBar, out=sys.stdout if ISATTY else DEV_NULL):
		mObj.setDatabase(args.database)

	if args.sessionName is not None: mObj.setSessionName(args.sessionName)

	with TerminalUpdater(f"Checking Reference Genomes:", category="ReferenceDownloader", hooks=mObj.hooks, threadNames=list(map(lambda a:f"{a}.fna", mObj.database[ReferencesTable.AssemblyName])), printer=LoadingBar, out=sys.stdout if ISATTY else DEV_NULL):
		mObj.setReferenceFiles()
	
	return mObj

def runJob(mObj : MetaCanSNPer, args : NameSpace, argsDict : dict):
	
	from MetaCanSNPer.modules.Database import ReferencesTable
	genomes = list(mObj.database[ReferencesTable.Genome])
	
	if args.mapper is not None:
		with TerminalUpdater(f"Creating Mappings:", category="Mappers", hooks=mObj.hooks, threadNames=genomes, printer=Spinner, out=sys.stdout if ISATTY else DEV_NULL):
			mObj.createMap(softwareName=args.mapper, flags=argsDict.get("--mapperOptions", {}))

	if args.aligner is not None:
		with TerminalUpdater(f"Creating Alignments:", category="Aligners", hooks=mObj.hooks, threadNames=genomes, printer=Spinner, out=sys.stdout if ISATTY else DEV_NULL):
			mObj.createAlignment(softwareName=args.aligner, flags=argsDict.get("--alignerOptions", {}))
	
	with TerminalUpdater(f"Calling SNPs:", category="SNPCallers", hooks=mObj.hooks, threadNames=genomes, printer=Spinner, out=sys.stdout if ISATTY else DEV_NULL):
		mObj.callSNPs(softwareName=args.snpCaller, flags=argsDict.get("--snpCallerOptions", {}))

def saveResults(mObj : MetaCanSNPer, args : argparse.Namespace):

	with TerminalUpdater(f"Saving Results:", category="SavingResults", hooks=mObj.hooks, threadNames=[mObj.queryName], printer=Spinner, out=sys.stdout if ISATTY else DEV_NULL):
		mObj.saveSNPdata()
		outDir = mObj.saveResults()
		mObj.hooks.trigger("SavingResultsFinished", {"name" : mObj.queryName})
	return outDir

def main(argVector : list[str]=sys.argv) -> int:
	
	# mainParser = argparse.ArgumentParser(prog=__package__, description=package.__doc__)
	# mainParser.add_argument("Mode", choices=["mainParser"], type=str.capitalize)

	argsDict = separateCommands(argVector)

	if len(argVector) < 2:
		parser.print_help()
		parser.exit()

	args : NameSpace = parser.parse_args(argsDict["args"], namespace=NameSpace())

	handleOptions(args)

	try:
		mObj = initializeMainObject(args)

		runJob(mObj, args, argsDict)

		outDir = saveResults(mObj, args)
	except Exception as e:
		LOGGER.exception(e)

		if "mObj" in globals():
			for exc in mObj.exceptions:
				print(f"{type(exc).__name__}: {exc}", file=sys.stderr)

		if not args.silent or not ISATTY: print(f"{SOFTWARE_NAME} ended before completing query. ", end="")
		print("Exception occurred: \n", file=sys.stderr)

		if args.debug:
			pattern = re.compile(r"(\[[\w.]+\] \d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3} - ERROR: .*?)(?:\[[\w.]+\] \d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3} - \w+?:|\$)", flags=re.DOTALL+re.MULTILINE)
			for i, err in enumerate(pattern.finditer(open(LOGGING_FILEPATH, "r").read())):
				print(f"Exception [{i}]\n{err.group(1)}", file=sys.stderr)
		else:
			print(f"{type(e).__name__}: "+str(e), file=sys.stderr)
		exit(1)
	else:
		if not args.silent or not ISATTY:
			print(f"{SOFTWARE_NAME} finished in {timer() - startTime:.3f} seconds! Results exported to: {outDir}")
	
	return 0