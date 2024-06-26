#!/usr/bin/env python3

from timeit import default_timer as timer
startTime = timer()
import logging, sys, argparse, traceback
from threading import Thread, Event
from time import sleep
from typing import Callable
import re

## import MetaCanSNPer specific modules
from MetaCanSNPer.Globals import *
from MetaCanSNPer.Globals import __version__
from MetaCanSNPer.core.MetaCanSNPer import MetaCanSNPer
from MetaCanSNPer.core.TerminalUpdater import TerminalUpdater, Spinner, LoadingBar, TextProgress
from MetaCanSNPer.core.Hooks import Hooks, GlobalHooks
import MetaCanSNPer.Globals as Globals
import MetaCanSNPer as package
from MetaCanSNPer.modules.SubSampling import splitFastq, subSampleName, SUB_SAMPLE_NAMES

LOGGER = LOGGER.getChild(__name__.split(".")[-1])
PARSE_INT_PATTERN = re.compile(r"(\d+(?:[.]\d+)?|(?:\d+)?[.]\d+)(?:([*]{2}|[\^])(\d+(?:[.]\d+)?|(?:\d+)?[.]\d+))?")
def parseInt(string):
	"""Parses strings of integers in either of the following three forms: N, N_e**N_p, N_e^N_p."""
	m = PARSE_INT_PATTERN.match(string)
	if m is None:
		return int(string)
	
	exponent, operator, power = m.groups()
	if power is not None:
		return round(float(exponent) ** float(power))
	return int(exponent)

class NameSpace(argparse.Namespace):

	listSoftware : bool = False
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

	sessionName : str = None
	
	subSampled : bool = property(lambda self: any([self.reads, self.coverage, self.dilution, self.bytes]))
	subSampleType : str = property(lambda self: next(filter(self.get, ["reads", "coverage", "dilution", "bytes"]), None))
	reads : list[int,int] = None
	coverage : list[int,int,int] = None
	dilution : list[int,int] = None
	bytes : list[int,int] = None

	saveTemp : bool = False
	debug : bool = False
	verbose : bool = False
	suppress : bool = False
	silent : bool = False
	dryRun : bool = False

	@overload
	def __init__(self,
		query : list[str],
		organism : str,
		database : str = None,
		mapper : str = None,
		aligner : str = None,
		snpCaller : str = None,
		reads : list[int,int] = None,
		coverage : list[int,int,int] = None,
		dilution : list[int,int] = None,
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
		dryRun : bool = False
		): ...

	__init__ = argparse.Namespace.__init__

	def get(self, name):
		return getattr(self, name, None)
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
	requiredArguments.add_argument("--query", nargs="+",	metavar=("FILE", "FILES"), type=FilePath, required=True, help="Raw sequence data file supported by the intended Aligner/Mapper.")
	requiredArguments.add_argument("--organism",			metavar="NAME",		required=True, help="Name of organism queried. (Use \"_\" in place of spaces)")

servicesArguments = parser.add_argument_group("Choosing Software to run", description="If no software is given, a "
											  "default will be used from your personal default flags or from a "
											  "settingsFile specified by the '--settingsFile' flag.")
if True:
	servicesArguments.add_argument("--mapper",				metavar="NAME",		help="Name of installed and supported Mapper software.")
	servicesArguments.add_argument("--aligner",				metavar="NAME",		help="Name of installed and supported Alignment software.")
	servicesArguments.add_argument("--snpCaller",			metavar="NAME",	help="Name of installed and supported SNP Calling software.")

optionalArguments = parser.add_argument_group("Optional arguments")
if True:
	optionalArguments.add_argument("-d", "--database",	metavar="FILE", type=FilePath,			help="Filename of CanSNP database to be used.")
	optionalArguments.add_argument("--saveTemp",		action="store_true",					help="Don't dispose of temporary directories/files.")
	optionalArguments.add_argument("--settingsFile",	metavar="FILE", type=FilePath,							help="Path to .TOML file containing settings for MetaCanSNPer. Check the 'defaultConfig.toml' to see what can be included in a settings file.")
	subSamplingArguments = optionalArguments.add_argument_group("Sub Sampling", description="Run sub samples of the given query, sub sampling by one of three modes. Values can be of either type N, N_e**N_p, N_e^N_p.")
	exclusiveSSArguments = subSamplingArguments.add_mutually_exclusive_group()
	exclusiveSSArguments.add_argument("--reads", nargs=2, metavar=("N", "M"), type=parseInt, default=None, help="Run N sub samples each consisting of M number of reads.")
	exclusiveSSArguments.add_argument("--coverage", nargs=3, metavar=("N", "M", "COVERAGE"), type=parseInt, default=None, help="Run N sub samples each with enough reads to get ~M coverage during mapping. COVERAGE is the expected coverage the query file would achieve on its own.")
	exclusiveSSArguments.add_argument("--dilute", nargs=2, metavar=("N", "M"), type=parseInt, default=None, help="Run N sub samples each diluted by a factor of M (`sampleBytes * (1/M)`).")
	exclusiveSSArguments.add_argument("--bytes", nargs=2, metavar=("N", "M"), type=parseInt, default=None, help="Run N sub samples each containing roughly M bytes.")

	# Not used by the argparser, but is used for the help-page and for splitting the argv
	mapperOptions = optionalArguments.add_argument("--mapperOptions",		metavar="FLAGS",		help=MAPPER_OPTIONS_EXPLAINER)
	alignerOptions = optionalArguments.add_argument("--alignerOptions",		metavar="FLAGS",		help=ALIGNER_OPTIONS_EXPLAINER)
	snpCallerOptions = optionalArguments.add_argument("--snpCallerOptions",	metavar="FLAGS",		help=SNP_CALLER_OPTIONS_EXPLAINER)

directoryOptions = parser.add_argument_group("Directory Options")
if True:
	directoryOptions.add_argument("-W", "--workDir",		metavar="DIRECTORY", type=DirectoryPath, default=None, help="Work directory")
	directoryOptions.add_argument("-U", "--userDir",		metavar="DIRECTORY", type=DirectoryPath, default=None, help="User directory")
	directoryOptions.add_argument("-I", "--installDir",		metavar="DIRECTORY", type=DirectoryPath, default=None, help="Installation directory")
	directoryOptions.add_argument("-Q", "--targetDir",		metavar="DIRECTORY", type=DirectoryPath, default=None, help="Target (Query) directory")
	directoryOptions.add_argument("-T", "--tmpDir",			metavar="DIRECTORY", type=DirectoryPath, default=None, help="Temporary directory")
	directoryOptions.add_argument("-R", "--refDir",			metavar="DIRECTORY", type=DirectoryPath, default=None, help="References directory")
	directoryOptions.add_argument("-D", "--databaseDir",	metavar="DIRECTORY", type=DirectoryPath, default=None, help="Databases directory")
	directoryOptions.add_argument("-O", "--outDir",			metavar="DIRECTORY", type=DirectoryPath, default=None, help="Output directory")
	directoryOptions.add_argument("-S", "--sessionName",	metavar="NAME", type=str, default=None, help="Session Name/Directory")

debugOptions = parser.add_argument_group("Logging and debug options")
if True:
	debugOptions.add_argument("--verbose",	action="store_true",	help="Verbose output")
	debugOptions.add_argument("--debug",	action="store_true",	help="Debug output")
	debugOptions.add_argument("--suppress",	action="store_true",	help="Suppress warnings")
	debugOptions.add_argument("--silent",	action="store_true",	help="Disables printing to terminal except for any error messages which might appear.")
	debugOptions.add_argument("--dryRun",	action="store_true",	help="Don't run the processes of the mapper/aligner/snpCaller, just run a randomised (1 - 5 sec) `sleep` call.")

def checkDependencies(args : NameSpace):

	import shutil
	from MetaCanSNPer.core.Wrappers import Aligner, Mapper, SNPCaller
	requiredDeps = []
	optionalDeps = ["samtools"]

	if args.mapper:
		requiredDeps.extend(Mapper.get(args.mapper).dependencies)
	if args.aligner:
		requiredDeps.extend(Aligner.get(args.aligner).dependencies)
	if args.snpCaller:
		requiredDeps.extend(SNPCaller.get(args.snpCaller).dependencies)

	missed = []
	for dep in requiredDeps:
		if not shutil.which(dep) and not shutil.which(dep+".exe"):
			missed.append(dep)
	if len(missed) == 1:
		raise MissingDependency(f"Missing required dependency: {missed[0]}.")
	elif missed:
		nt = "\n\t"
		raise MissingDependency(f"Missing required dependencies:\n{nt.join(missed)}.")
	
	missed = []
	for dep in optionalDeps:
		if not shutil.which(dep) and not shutil.which(dep+".exe"):
			missed.append(dep)
	if len(missed) == 1:
		print(f"Missing optional dependency: {missed[0]}.")
		while (string := input("You may still run without it, do you want to run anyway [Y/N]? ").strip().lower()) not in ["y", "n"]: pass
		if string == "n":
			exit(2)
		else:
			pass
	elif missed:
		nt = "\n\t"
		print(f"Missing optional dependencies:\n{nt.join(missed)}.")
		while (string := input("You may still run without them, do you want to run anyway [Y/N]? ").strip().lower()) not in ["y", "n"]: pass
		if string == "n":
			exit(2)
		else:
			pass

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

	if args.dryRun:
		Globals.DRY_RUN = args.dryRun # Don't run the processes of the mapper/aligner/snpCaller, just run a randomised `sleep` call
	if args.saveTemp:
		Globals.PPGlobals.DISPOSE = False # Don't dispose of temporary directories/files.
	
	if args.debug:
		logging.basicConfig(level=logging.DEBUG)
		LOGGING_FILEHANDLER.setLevel(logging.DEBUG)
		Globals.DEBUG = True
		Globals.MAX_DEBUG = True
		Globals.SQLOOPGlobals.MAX_DEBUG = True
	elif args.verbose:
		logging.basicConfig(level=logging.INFO)
	elif args.suppress:
		logging.basicConfig(level=logging.ERROR)
	else:
		pass # The default logging level for the logging package is logging.WARNING

@overload
def initializeData(args : NameSpace, /) -> list[tuple[str]]: ...
@overload
def initializeData( *, organism : str, query : list[FilePath], reads : list[int, int]) -> list[tuple[str]]: ...
@overload
def initializeData( *, organism : str, query : list[FilePath], coverage : list[int, int, int]) -> list[tuple[str]]: ...
@overload
def initializeData( *, organism : str, query : list[FilePath], dilution : list[int, int]) -> list[tuple[str]]: ...
@overload
def initializeData( *, organism : str, query : list[FilePath], bytes : list[int, int]) -> list[tuple[str]]: ...
def initializeData(args : NameSpace|None=None, /, **kwargs) -> list[FileList[FilePath]]:
	
	args = NameSpace(**(dict(args._get_kwargs()) | kwargs))
	if not args.subSampled:
		return [FileList(args.query)]

	organism = args.organism
	query = args.query
	
	from MetaCanSNPer.core.DirectoryLibrary import DirectoryLibrary

	DL = DirectoryLibrary(organism, query)
	if args.saveTemp:
		outDir = DL.dataDir.create("SubSampling").create(DL.queryName).create(f"{SUB_SAMPLE_NAMES[args.subSampleType]}-{'-'.join(map(shortNumber, args[args.subSampleType]))}")
	else:
		outDir = PseudoPathyFunctions.createTempDir(f"{SUB_SAMPLE_NAMES[args.subSampleType]}-{'-'.join(map(shortNumber, args[args.subSampleType]))}")
	
	with TerminalUpdater(f"Creating Sub-samples:", category="SplitFastq", names=[DL.queryName], hooks=GlobalHooks, printer=LoadingBar, length=65, out=sys.stdout if ISATTY else DEV_NULL) as TU:
		
		newFiles = splitFastq(args[args.subSampleType][0], DL.query, outDir=outDir, hooks=TU.hooks, **{args.subSampleType:args[args.subSampleType][1:]})

	return newFiles

@overload
def initializeMainObjects(args : NameSpace, /) -> tuple[str,list[MetaCanSNPer]]: ...
@overload
def initializeMainObjects( *, organism : str, query : list[FilePath], queryFiles : list[FileList[str]], database : FilePath|None=None, subSampled : bool, saveTemp : bool, settingsFile : FilePath) -> tuple[str,list[MetaCanSNPer]]: ...
def initializeMainObjects(args : NameSpace=None, /, *, organism : str|None=None, query : list[FilePath]|None=None, queryFiles : list[FileList[str]]|None=None, database : FilePath|None=None, subSampled : bool|None=None, saveTemp : bool|None=None, settingsFile : FilePath|None=None) -> tuple[str,list[MetaCanSNPer]]:

	from MetaCanSNPer.modules.Database import ReferencesTable
	from MetaCanSNPer.core.Hooks import GlobalHooks, Hooks
	from MetaCanSNPer.core.DirectoryLibrary import DirectoryLibrary
	
	organism = organism or args.organism
	query : FileList[FilePath] = FileList(query or args.query)
	queryName = query.name
	subSampled = queryFiles is not None
	queryFiles = queryFiles if subSampled else []
	settings = vars(args)
	
	if args is not None:
		database = database or args.database
		saveTemp = saveTemp if saveTemp is not None else args.saveTemp
		settingsFile = settingsFile or args.settingsFile
	
	if not args.subSampled:
		instances = [MetaCanSNPer(organism, query, settings=settings)]
		groupSessionName = instances[0].sessionName
	else:
		sampleType = args.subSampleType
		sampleINTs = args[sampleType]
		tmpDir = None

		timeString = time.strftime('-%Y-%m-%d_%H-%M-%S', time.localtime())
		groupSessionName = subSampleName(args.sessionName or queryName, sampleType, *sampleINTs) + timeString

		DL = DirectoryLibrary(organism, query)

		if args.saveTemp:
			timeString = ""
			if settings.get("tmpDir") is None:
				tmpDir = DL.userCacheDir / groupSessionName
		
		instances : list[MetaCanSNPer] = [
			MetaCanSNPer(
				organism, query,
				settings=settings, settingsFile=args.settingsFile,
				hooks=Hooks() if i > 0 else GlobalHooks,
				sessionName=(subSessionName := subSampleName(args.sessionName or queryName, sampleType, *sampleINTs, index=i+1) + timeString),
				tmpDir=tmpDir / subSessionName if tmpDir is not None else None)
			for i, query in enumerate(queryFiles)
		]
	mObj = instances[0]

	database = database or mObj.databaseName
	
	with TerminalUpdater(f"Checking database {database!r}:", category="DatabaseDownloader", names=[database], printer=LoadingBar, length=30, out=sys.stdout if ISATTY else DEV_NULL) as TU:
		mObj.setDatabase(database, sequential=True)
		GlobalHooks.trigger("DatabaseDownloaderPostProcess", {"name" : database, "value" : 1.0})
		for obj in instances[1:]:
			obj.databaseName = mObj.databaseName
			obj.Lib.database = obj.database = mObj.database
		GlobalHooks.trigger("DatabaseDownloaderFinished", {"name" : database, "value" : 3})

	refFiles = [f"{assemblyName}.fna" for *_, assemblyName in mObj.database.references]
	with TerminalUpdater(f"Checking Reference Genomes:", category="ReferenceDownloader", names=refFiles, printer=LoadingBar, length=30, out=sys.stdout if ISATTY else DEV_NULL) as TU:
		mObj.setReferenceFiles(sequential=True)
		for refFile in refFiles:
			GlobalHooks.trigger("ReferenceDownloaderProgress", {"name" : refFile, "value" : 1.0})
		for instance in instances[1:]:
			instance.setReferenceFiles()
		for refFile in refFiles:
			GlobalHooks.trigger("ReferenceDownloaderProgress", {"name" : refFile, "value" : 3})
	
	return groupSessionName, instances

def runAndUpdate(mObj, failedEvent, func, lock, softwareName, flags, TU, categoryName, jobs):
	try:
		func(mObj, softwareName=softwareName, flags=flags)
		with lock:
			for name, value in TU.threads.items():
				TU.hooks.trigger(f"{categoryName}Progress", {"name" : name, "value" : value + 1 / jobs}, block=True)
		return 0
	except Exception as e:
		LOGGER.exception(e)
		failedEvent.set()
		for name, value in TU.threads.items():
			TU.hooks.trigger(f"{categoryName}Failed", {"name" : name, "value" : None})
		return 1

def runJobs(instances, func, args, argsDict, category, categoryName, names, message):
	jobs = len(instances)
	commonLock = Lock()
	LocalHooks = Hooks()
	
	with TerminalUpdater(message, category=categoryName, hooks=LocalHooks, names=names, printer=LoadingBar, length=30, out=sys.stdout if ISATTY else DEV_NULL) as TU:
		for name in names:
			LocalHooks.trigger(f"{categoryName}Starting", {"name" : name, "value" : 0})
		failedEvent = Event()
		for mObj in instances:
			if 0 != runAndUpdate(mObj, failedEvent, func, commonLock, args[category] or mObj.settings[category], argsDict.get(f"--{category}Options", {}), TU, categoryName, jobs):
				raise ChildProcessError(f"{categoryName} process failed.")
		for name in names:
			LocalHooks.trigger(f"{categoryName}Finished", {"name" : name, "value" : 3})

def runPrograms(instances : list[MetaCanSNPer], args : NameSpace, argsDict : dict):
	
	from MetaCanSNPer.modules.Database import ReferencesTable
	
	genomes = list(instances[0].database[ReferencesTable.Genome])
	queryFormat = instances[0].query[0].lower()
	if args.subSampled:
		if args.mapper or any(queryFormat.endswith(ext) for ext in ["fastq", "fq", "fastq.gz", "fq.gz"]):
			runJobs(instances, MetaCanSNPer.createMap, args, argsDict, "mapper", "Mappers", genomes, "Creating Mappings:")
		
		if args.aligner or any(queryFormat.endswith(ext) for ext in ["fasta", "fna", "fasta.gz", "fna.gz"]):
			runJobs(instances, MetaCanSNPer.createAlignment, args, argsDict, "aligner", "Aligners", genomes, "Creating Alignments:")

		runJobs(instances, MetaCanSNPer.callSNPs, args, argsDict, "snpCaller", "SNPCallers", genomes, "Calling SNPs:")
	else:
		mObj = instances[0]
		
		if args.mapper or any(queryFormat.endswith(ext) for ext in ["fastq", "fq", "fastq.gz", "fq.gz"]):
			with TerminalUpdater(f"Creating Mappings:", category="Mappers", hooks=mObj.hooks, names=genomes, printer=Spinner, out=sys.stdout if ISATTY else DEV_NULL):
				
				mObj.createMap(softwareName=args.mapper or mObj.settings["mapper"], flags=argsDict.get("--mapperOptions", {}))
		
		if args.aligner or any(queryFormat.endswith(ext) for ext in ["fasta", "fna", "fasta.gz", "fna.gz"]):
			with TerminalUpdater(f"Creating Alignments:", category="Aligners", hooks=mObj.hooks, names=genomes, printer=Spinner, out=sys.stdout if ISATTY else DEV_NULL):
				
				mObj.createAlignment(softwareName=args.aligner or mObj.settings["aligner"], flags=argsDict.get("--alignerOptions", {}))
		
		with TerminalUpdater(f"Calling SNPs:", category="SNPCallers", hooks=mObj.hooks, names=genomes, printer=Spinner, out=sys.stdout if ISATTY else DEV_NULL):
			
			mObj.callSNPs(softwareName=args.snpCaller or mObj.settings["snpCaller"], flags=argsDict.get("--snpCallerOptions", {}))

def saveResults(instances : list[MetaCanSNPer], args : NameSpace, sessionName : str) -> Path:
	
	if not args.subSampled:
		mObj = instances[0]

		with TerminalUpdater(f"Saving Results:", category="SavingResults", hooks=mObj.hooks, names=[mObj.queryName], printer=Spinner, out=sys.stdout if ISATTY else DEV_NULL):
			mObj.hooks.trigger("SavingResultsStarting", {"name" : mObj.queryName, "value" : 0})
			mObj.saveSNPdata()
			outDir = mObj.saveResults()
			mObj.hooks.trigger("SavingResultsFinished", {"name" : mObj.queryName, "value" : 3})

		return outDir
	else:
		from MetaCanSNPer.core.DirectoryLibrary import DirectoryLibrary

		DL = DirectoryLibrary(args.organism, args.query, sessionName=sessionName)
		realOutDir = DL.resultDir
		jobs = len(instances)
		LocalHooks = Hooks()
		name = DL.queryName

		with TerminalUpdater(f"Saving Results:", category="SavingResults", hooks=LocalHooks, names=[name], printer=LoadingBar, length=65, out=sys.stdout if ISATTY else DEV_NULL):
			LocalHooks.trigger("SavingResultsStarting", {"name" : name, "value" : 0})
			for i, mObj in enumerate(instances):
				mObj.saveSNPdata()
				outDir = mObj.saveResults()
				for filename in map(FilePath, os.listdir(outDir)):
					os.rename(
						outDir / filename,
						realOutDir / subSampleName(filename, args.subSampleType, *args[args.subSampleType], index=i+1)
					)
				shutil.rmtree(outDir, ignore_errors=True)
				LocalHooks.trigger("SavingResultsProgress", {"name" : name, "value" : (i+1) / jobs})
			LocalHooks.trigger("SavingResultsFinished", {"name" : name, "value" : 3})
		
		return DirectoryPath(realOutDir)

def main(argVector : list[str]=sys.argv) -> int:
	
	# mainParser = argparse.ArgumentParser(prog=__package__, description=package.__doc__)
	# mainParser.add_argument("Mode", choices=["mainParser"], type=str.capitalize)
	
	argsDict = separateCommands(argVector)

	if len(argVector) < 2:
		parser.print_help()
		parser.exit()

	args : NameSpace = parser.parse_args(argsDict["args"], namespace=NameSpace())

	print(f"\nRunning {SOFTWARE_NAME}...\n", file=sys.stderr)

	if not args.dryRun:
		checkDependencies(args)
	
	handleOptions(args)

	filenames = initializeData(args)

	sessionName, instances = initializeMainObjects(args, queryFiles=filenames)

	runPrograms(instances, args, argsDict)

	outDir = saveResults(instances, args, sessionName)

	print(f"Results exported to:\n\t{outDir}", file=sys.stderr)

	return 0

from functools import wraps
@wraps(main)
def _main_wrapper(mainFunc, *args, **kwargs):
	errno = 1
	try:
		errno = mainFunc(*args, **kwargs)
	except Exception as e:
		LOGGER.exception(e)

		message = str(e)
		indentation = "    "
		rowLength = 80 - len(indentation)
		rows = bool(len(message) % rowLength) + len(message) // rowLength
		message = f"\n{indentation}".join(message[i*rowLength:(i+1)*rowLength] for i in range(rows))
		if Globals.DEBUG:
			Globals.LOGGING_ERRORHANDLER.flush()
			print("".join(Globals.LOGGING_ERRORMESSAGES), file=sys.stderr)
			print(f"{SOFTWARE_NAME} ended before completing query. Exceptions that occurred are listed above.", file=sys.stderr)
		else:
			print(f"{SOFTWARE_NAME} ended before completing query. ", file=sys.stderr)
			print(f"Due to `{e.__class__.__name__}`:\n{message}", file=sys.stderr)
		print(flush=True, file=sys.stderr)
		errno = getattr(e, "errno", 1)
	else:
		print(f"{SOFTWARE_NAME} finished in {timer() - startTime:.3f} seconds!", flush=True, file=sys.stderr)
	finally:
		exit(errno)

main = _main_wrapper.__get__(main)