#!/usr/bin/env python3

import logging, sys, argparse, traceback
from threading import Thread
from time import sleep
from typing import Callable

## import MetaCanSNPer specific modules
from MetaCanSNPer.Globals import *
from MetaCanSNPer.Globals import __version__
import MetaCanSNPer.modules.LogKeeper as LogKeeper
from MetaCanSNPer.modules.MetaCanSNPer import MetaCanSNPer, Hooks
import MetaCanSNPer.Globals as Globals

LOGGER = LogKeeper.createLogger(__name__)


def showLoadingSymbol(running : bool, progress : list[float], symbols : list[str]=["\b|", "\b/", "\b-", "\b\\"], sep=" ", borders=("[", "]")):
	if type(progress) is float:
		progress = [progress]

	sepLength = len(sep)
	N = len(symbols)
	backspaces = "\b" * (len(progress)+borders[0]*len(progress)+borders[1]*len(progress)+sepLength*min(0, len(progress)-1))
	n = [0 for _ in range(len(progress))]
	while running:
		print(backspaces, end="", flush=True)
		for i, prog in enumerate(progress):
			if running:
				n[i]=(n[i]+1)%len(symbols)
				print(borders[0]+symbols[n[i]]+borders[1]+sep, end="", flush=True)
			else:
				backspaces = "\b" * (i+borders[0]*i+borders[1]*i+sepLength*min(0, i-1))
				print(backspaces, end="", flush=True)
				print(backspaces.replace("\b", " "), end="", flush=True)
				print(backspaces, end="", flush=True)
				return
		print("\b \b", end="", flush=True)
		sleep(0.2)

def showLoadingMiniBars(running : bool, progress : list[float], symbols : list[str]= [".", "_", "\u2584", "", "\u2588"], sep=" ", borders=("[", "]")):
	if type(progress) is float:
		progress = [progress]

	sepLength = len(sep)
	N = len(progress)
	backspaces = "\b" * (len(progress)+borders[0]*len(progress)+borders[1]*len(progress)+sepLength*min(0, len(progress)-1))
	while running:
		print(backspaces, end="", flush=True)
		for i, prog in enumerate(progress):
			if running:
				print(borders[0]+symbols[int(N*prog)]+borders[1]+sep, end="", flush=True)
			else:
				backspaces = "\b" * (i+borders[0]*i+borders[1]*i+sepLength*min(0, i-1))
				print(backspaces, end="", flush=True)
				print(backspaces.replace("\b", " "), end="", flush=True)
				print(backspaces, end="", flush=True)
				return
		print("\b \b", end="", flush=True)
		sleep(0.5)

def showLoadingBar(running, progress : float|list[float], length=10, border=("[", "]"), fill="\u2588", halfFill="\u258C", background=" ", sep=" "):
	if type(progress) is float:
		progress = [progress]
	innerLength = length - len(border[0]) - len(border[1])
	sepLength = len(sep)
	backspaces = "\b" * (length * len(progress) + sepLength * min(0, len(progress) - 1))
	while running:
		print(f"{backspaces}", end="", flush=True)
		for i, prog in enumerate(progress):
			if running:
				fillLength = int(innerLength*2*prog)
				fillLength, halfBlock = fillLength//2, fillLength%2
				emptyLength = innerLength - fillLength - halfBlock
				print(f"{border[0]}{fill*fillLength}{halfFill*halfBlock}{background*emptyLength}{border[1]}{sep}", end="", flush=True)
			else:
				backspaces = "\b" * (length * i + sepLength * min(0, i - 1))
				print(backspaces, end="", flush=True)
				print(backspaces.replace("\b", " "), end="", flush=True)
				print(backspaces, end="", flush=True)
				return
				
			print("\b \b"*sepLength, end="", flush=True)
		sleep(0.6)

def updateTerminal(msg, category, hooks : Hooks, nThreads):


	def updateProgress(eventInfo : dict, progress : float, threads : dict[int,float]):
		threads[eventInfo["threadN"]] = eventInfo["progress"]
		progress -= progress
		progress += sum(threads.values()) / len(threads)

	print(f"{msg} ... ", end="", flush=True)
	notDone = True
	progress = 0.0
	threads = {key+1:0.0 for key in range(nThreads)}

	hooks.addHook(f"{category}Progress", updateProgress, args=[progress, threads])
	t = Thread(target=showLoadingBar, args=[notDone, progress], daemon=True)
	t.start()

	def stopLoadingbar(eventInfo : dict, t : Thread, notDone : bool, finishedThreads : dict):
		finishedThreads[eventInfo["threadN"]] = 1
		if sum(finishedThreads.values()) == len(finishedThreads):
			notDone *= False
			print("\bDone!", flush=True)
	
	hooks.addHook(f"{category}Finished", stopLoadingbar, args=[t, notDone, {}])

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
		debugOptions.add_argument("--dry-run",	action="store_true",	help="Don't run the processes of the mapper/aligner/snpCaller, just run a randomised (1 - 5 sec) `sleep` call.")
	

	return parser

def main():
	
	argsDict = separateCommands(sys.argv)
	
	parser = createParser()

	if len(sys.argv)<2:
		parser.print_help()
		parser.exit()

	args : argparse.Namespace = parser.parse_args(argsDict["args"])

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

	flags = dict(args._get_kwargs())
	try:
		mObj = MetaCanSNPer(settings=flags, settingsFile=args.settingsFile)
		
		mObj.setQuery(flags["query"])
		mObj.setDatabase(flags["database"])
		
		if flags["sessionName"] is not None: mObj.setSessionName(flags["sessionName"])
		
		if flags.get("mapper") is not None:
			updateTerminal("Creating maps ... ", "Mappers", mObj.hooks, len(mObj.database.references))

			mObj.createMap(softwareName=flags["mapper"], flags=argsDict.get("--mapperOptions", {}))

		if flags.get("aligner") is not None:
			updateTerminal("Creating Alignments ... ", "Aligners", mObj.hooks, len(mObj.database.references))

			mObj.createAlignment(softwareName=flags["aligner"], flags=argsDict.get("--alignerOptions", {}))
		
		updateTerminal("Calling SNPs", "SNPCallers", mObj.hooks, len(mObj.database.references))

		mObj.callSNPs(softwareName=flags["snpCaller"], flags=argsDict.get("--snpCallerOptions", {}))

		print(f"{SOFTWARE_NAME} finished! Results exported to: {mObj.Lib.resultDir}")
	except Exception as e:
		LOGGER.exception(e)
		print(f"{SOFTWARE_NAME} ended before completing query. Exception that caused it:")
		print()
		if args.debug:
			raise e
		else:
			string : list = traceback.format_exc().split("\n")
			output = []
			for i in range(len(string[1:])):
				row : str = string[i+1]
				if row.strip("\r")[:1] in [" ", "\t"]:
					continue
				else:
					for j in range(i, len(string[1:])):
						output.append(string[j+1])
					break
			print("\n".join(output))


if oname=="__main__":
	main()
	
