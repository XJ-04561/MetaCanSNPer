#!/usr/bin/env python3
"""
MetaCanSNPer
"""

oname = __name__
__name__ 		= "MetaCanSNPer"
__version__ 	= "3.0.0"
__author__ 		= "Fredrik Sörensen"
__credits__ 	= ["Fredrik Sörensen", "David Sundell"]
__license__ 	= "GPLv3"
__maintainer__ 	= "FOI bioinformatics group"
__email__ 		= ["bioinformatics@foi.se", "fredrik.sörensen@foi.se"]
__date__ 		= "2024-02-27"
__status__ 		= "Prototype"


import logging
try:
	## import MetaCanSNPer specific modules
	import MetaCanSNPer.modules.LogKeeper as LogKeeper
	from MetaCanSNPer.modules.DirectoryLibrary import DirectoryLibrary
	from MetaCanSNPer.modules.MetaCanSNPer import MetaCanSNPer
except:
	## import MetaCanSNPer specific modules
	import modules.LogKeeper as LogKeeper
	from modules.DirectoryLibrary import DirectoryLibrary
	from modules.MetaCanSNPer import MetaCanSNPer

LOGGER = LogKeeper.createLogger(__name__)

import os,sys
## Basic config file

"""MetaCanSNPer settings"""
import argparse

DIRECTORY_OPTIONS = ["workDir", "userDir", "installDir", "targetDir", "tmpDir", "refDir", "databaseDir", "outDir", "sessionName"]

INDEXING_OPTIONS_EXPLAINER = """
To provide flags/arguments for the chosen Mapper or Aligner, provide them
directly after the '--indexerOptions' flag, only interrupted by the end of the
command call or the corresponding flag for SNP Caller options.
"""
SNP_CALLER_OPTIONS_EXPLAINER = """
To provide flags/arguments for the chosen SNP Caller, provide them directly
after the '--snpCallerOptions' flag, only interrupted by the end of the
command call or the corresponding flag for Mapper or Aligner options.
"""

def main():
	"""Initiate CanSNPer2 object"""
	parser = argparse.ArgumentParser(description="MetaCanSNPer")

	# The 'if True:' structures are used to minimize and expand sections in an IDE.

	requiredArguments = parser.add_argument_group("Required arguments")
	if True:
		requiredArguments.add_argument("-q", "--query",		metavar="query",		help="Raw sequence data file supported by the intended Aligner/Mapper.")
		requiredArguments.add_argument("-d", "--database",	metavar="database",		help="Filename of CanSNP database to be used.")
		requiredArguments.add_argument("--snpCaller",		metavar="snpCaller",	help="Name of installed and supported SNP Calling software.")
	
	mapOrAlign = parser.add_mutually_exclusive_group(required=True)
	if True:
		mapOrAlign.add_argument("--mapper",					metavar="mapper",		help="Name of installed and supported mapper software.")
		mapOrAlign.add_argument("--aligner",				metavar="aligner",		help="Name of installed and supported alignment software.")

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

	indexingOptions = parser.add_argument_group(title="--indexerOptions : Indexing options", description=INDEXING_OPTIONS_EXPLAINER)

	snpCallerOptions = parser.add_argument_group(title="--snpCallerOptions : SNP Caller options", description=SNP_CALLER_OPTIONS_EXPLAINER)

	debugOptions = parser.add_argument_group("Logging and debug options")
	if True:
		debugOptions.add_argument("--verbose",			action="store_const", const=logging.INFO,					help="Verbose output")
		debugOptions.add_argument("--debug",			action="store_const", const=logging.DEBUG,					help="Debug output")
		debugOptions.add_argument("--supress",	action="store_const", const=logging.ERROR,	default=logging.WARNING,	help="Supress warnings")
	
	parser.add_argument("--version", action="store_true", help=argparse.SUPPRESS)

	args, remaining = parser.parse_known_args()
	if len(sys.argv)==1:
		parser.print_help()
		parser.exit()
	"""
		Setup logging and debug options
	"""
	if args.version:
		print("MetaCanSNPer - version {version}".format(version=__version__))
		exit()

	givenDirectories = {d:args[d] for d in DIRECTORY_OPTIONS if d in args}
	DL = DirectoryLibrary(**givenDirectories)

	mObj = MetaCanSNPer(lib=DL)

	mObj.setQuery(args.query)
	mObj.setDatabase(args.database)
	mObj.connectDatabase()

	if args.sessionName is not None: mObj.setSessionName(args.sessionName)

	indexerArgs, snpCallerArgs = [], []
	if len(remaining) > 0:
		if remaining[0] not in ["--indexerOptions", "--snpCallingOptions"]:
			LOGGER.critical("Unknown flags provided not as belonging to a Mapper, Aligner, or SNP Caller.")
			exit()
		current = indexerArgs
		for f in remaining:
			## I AM HERE


	if "mapper" in args:
		mObj.createIndex(softwareName=args.mapper)
	elif "aligner" in args:
		mObj.createIndex(softwareName=args.aligner)

	mObj.callSNPs(softwareName=args.snpCaller)

if oname=="__main__":
	main()
