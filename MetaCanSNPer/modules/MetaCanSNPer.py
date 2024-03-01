'''
MetaCanSNPer module: A toolkit for SNP-typing using NGS data.
Copyright (C) 2019 David Sundell @ FOI bioinformatics group
'''

import os
import time
import logging
try:
	## import MetaCanSNPer specific modules
	import MetaCanSNPer.modules.LogKeeper as LogKeeper
	from MetaCanSNPer.modules.Databases import DatabaseReader
	from MetaCanSNPer.CanSNPerTree import __version__
	from MetaCanSNPer.modules.DirectoryLibrary import DirectoryLibrary
	from MetaCanSNPer.modules.Wrappers import Aligner, Mapper, SNPCaller, IndexingWrapper
	import MetaCanSNPer.modules.Aligners as Aligners
	import MetaCanSNPer.modules.Mappers as Mappers
	import MetaCanSNPer.modules.SNPCallers as SNPCallers
except:
	## import MetaCanSNPer specific modules
	import LogKeeper as LogKeeper
	from Databases import DatabaseReader
	from DirectoryLibrary import DirectoryLibrary
	from Wrappers import Aligner, Mapper, SNPCaller, IndexingWrapper
	import Aligners
	import Mappers
	import SNPCallers

LOGGER = LogKeeper.createLogger(__name__)

try:
	import tomllib as toml
except ModuleNotFoundError:
	try:
		import tomli as toml
	except:
		raise ModuleNotFoundError("TOML-reading module not found. For Python 3.5> it should be included as 'tomllib', if older than 2.5, then install 'tomli'.")

import random
random.seed()

SOFTWARE_NAME = "MetaCanSNPer"

# OS Alibis
pSep = os.path.sep
pJoin = os.path.join
pIsAbs = lambda *args, **kwargs: os.path.isabs(os.path.expanduser(*args, **kwargs))
pExpUser = os.path.expanduser
pAbs = os.path.abspath
pNorm = os.path.normpath
pDirName = os.path.dirname

def loadFlattenedTOML(filename):
	tmp : dict[str,dict] = toml.load(open(filename, "rb"))
	# Settings hierarchy looks like this: ["Category"]["Flag"] -> Value
	# Flatten hierarchy so flags are easily accessible
	settings = {}
	for flags in tmp.values():
		for flag, value in flags.items():
			settings[flag] = value
	return settings

class MetaCanSNPer:
	outputTemplate = "{0[refName]}_{0[queryName]}.{0[outFormat]}"
	databaseName : str
	database : DatabaseReader
	Lib : DirectoryLibrary
	settings : dict
	sessionName : str

	"""docstring for MetaCanSNPer"""
	def __init__(self, lib : DirectoryLibrary=None, database=None, settings : dict={}, settingsFile : str=None, sessionName : str=None):
		self.startTime = time.localtime()
		
		if lib is None:
			self.Lib = DirectoryLibrary(settings=self.settings)
		else:
			self.Lib = lib

		self.settings : dict = loadFlattenedTOML(os.path.join(self.Lib.installDir, "defaultFlags.toml"))
		if settingsFile is not None:
			if pIsAbs(settingsFile):
				self.settings.update( loadFlattenedTOML(settingsFile))
			else:
				self.settings.update( loadFlattenedTOML(self.Lib.targetDir[settingsFile]))
		self.settings.update(settings)

		if database is not None:
			self.setDatabase(database=database)
	
		self.sessionName = sessionName

		# Need for this is unknown, might remove later.
		self.summarySet = set()
		self.calledGenome = {}

	'''Database functions'''

	def setDatabase(self, database : str):
		if database in self.Lib.databaseDir:
			path = self.Lib.databaseDir[database]
			self.databaseName = path
			self.connectDatabase()
		else:
			LOGGER.warning("Database not found: '{}'".format( database))

	def connectDatabase(self):
		if self.databaseName is None:
			LOGGER.error("Database not specified. Can't connect unless a valid database is set through MetaCanSNPer.setDatabase.")
			return 1
		try:
			self.database = DatabaseReader(self.databaseName)
			return 0
		except Exception as e:
			LOGGER.error(e)
		return 1

	'''MetaCanSNPer set directories'''
	
	def setTargetDir(self, path : str):						self.Lib.setTargetDir(path)
	def setRefDir(self, genomeName : str, path : str):		self.Lib.setRefDir(genomeName, path)
	def setDatabaseDir(self, path : str):					self.Lib.setDatabaseDir(path)
	def setTmpDir(self, path : str):						self.Lib.setTmpDir(path)
	def setOutDir(self, path : str):						self.Lib.setOutDir(path)

	# Setting the same docstring for each 'set*Dir' function
	setTargetDir.__doc__ = setRefDir.__doc__ = setDatabaseDir.__doc__ = setTmpDir.__doc__ = setOutDir.__doc__ = """
	'path' can be relative or absolute. Check the docstring for the DirectoryLibrary to see where the relative
	path will start from.
	"""

	'''MetaCanSNPer set values'''

	def setQuery(self, query : str):
		self.Lib.setQuery(query)
		self.queryName, self.queryFormat = os.path.splitext(os.path.basename(self.Lib.query)) ## get name of file and remove ending
		if self.sessionName is None:
			self.Lib.setSessionName("Sample-{queryName}-{dateYYYYMMDD}".format(queryName=self.queryName, dateYYYYMMDD="{:0>4}.{:0>2}.{:0>2}.{:>2}{:>2},{:>2}".format(*(self.startTime[:6]))))

	def setSessionName(self, name):
		self.sessionName = name
		self.Lib.setSessionName(name)

	'''MetaCanSNPer get functions'''

	def getReferences(self):
		'''Fetch names of reference genomes in connected database. Download any reference genomes not present locally.'''
		return self.database.references

	def getQuery(self):
		return self.Lib.query

	'''Indexing methods'''

	def getSoftwareClass(self, softwareName):
		for types in [Aligners, Mappers, SNPCallers]:
			if types.get(softwareName) is not None:
				return types.get(softwareName)
		return None

	def runSoftware(self, softwareClass : IndexingWrapper, outputDict : dict={}, kwargs : dict={}):
		'''Align sequences using subprocesses.'''

		LOGGER.info("Creating software wrapper for '{}' of type '{}'".format(softwareClass.softwareName, softwareClass.__name__))
		software : IndexingWrapper = softwareClass(self.Lib, self.database, self.outputTemplate, kwargs=kwargs)

		# Check that error did not occur.
		while software.hickups():	
			LOGGER.info("Checking for pre-processing for {software}.".format(software=software.softwareName))
			software.preProcess()

			LOGGER.info("Starting {software}.".format(software=software.softwareName))
			output = software.start()
			
			# Run the alignment just like before, but hopefully fixed.
			software.updateWhileWaiting(outputDict)
			LOGGER.info("Finished running all instances of {software}.".format(software=software.softwareName))

			if software.fixable():
				LOGGER.info("{software} failed. Fix for the specific error exists. Implementing and running again.".format(software=software.softwareName))
				software.planB()
			else:
				LOGGER.error("Software error, no fix implemented. Returning empty list of outputs")
				for key, path in output:
					del outputDict[key]
				return outputDict
		
		return outputDict
	
	def createMap(self, softwareName : str, kwargs : dict={}):
		''''''
		MapperType : Mapper = Mappers.get(softwareName)

		LOGGER.info("Loading References from database.")
		self.database.references
		LOGGER.info("Loaded a total of {n} References.".format(n=len(self.database.references)))
		
		self.runSoftware(MapperType, outputDict=self.Lib.SNPs, kwargs=kwargs)

	def createAlignment(self, softwareName : str, kwargs : dict={}):
		''''''
		AlignerType : Aligner = Aligners.get(softwareName)

		LOGGER.info("Loading References from database.")
		self.database.references
		LOGGER.info("Loaded a total of {n} References.".format(n=len(self.database.references)))
		
		self.runSoftware(AlignerType, outputDict=self.Lib.SNPs, kwargs=kwargs)

	def callSNPs(self, softwareName : str, kwargs : dict={}):
		''''''
		SNPCallerType : SNPCaller = SNPCallers.get(softwareName)

		LOGGER.info("Loading References from database.")
		self.database.references
		LOGGER.info("Loading SNPs from database.")
		self.database.SNPsByGenome
		LOGGER.info("Loaded a total of {n} SNPs.".format(n=sum(len(SNPs) for SNPs in self.database.SNPsByGenome.values())))
		
		self.runSoftware(SNPCallerType, outputDict=self.Lib.SNPs, kwargs=kwargs)

	def traverseTree(self):
		'''Depth-first tree search.'''
		nodeID = 2
		snpID = self.database.nodes[nodeID]
		path = [snpID]

		while nodeID in self.database.tree:
			for childID in self.database.tree[nodeID]:
				childSNPID = self.database.nodes[childID]
				if self.database.SNPsByID[childSNPID][2]  == self.SNPresults[childSNPID]:
					nodeID = childID
					break
			if nodeID != childID:
				# Node is not a leaf, but no children are matches for the target.
				break
			path.append(childSNPID)
		
		return path

	'''Functions'''

	def saveResults(self, where : str=None):
		if where is None:
			pathFile = open(os.path.join(self.Lib.resultDir, self.Lib.queryName+"_path.tsv"), "w")
			finalFile = open(os.path.join(self.Lib.resultDir, self.Lib.queryName+"_final.tsv"), "w")
		else:
			pathFile = open(os.path.join(where, self.Lib.queryName+"_path.tsv"), "w")
			finalFile = open(os.path.join(where, self.Lib.queryName+"_final.tsv"), "w")

		path = self.traverseTree()

		pathFile.write("\t".join(path))
		finalFile.write(path[-1])
		
		pathFile.close()
		finalFile.close()


	def saveSNPdata(self, where : str=None):
		""""""
		if where is None:
			called = open(self.Lib.resultDir > self.Lib.queryName+"_snps.tsv", "w")
			notCalled = open(self.Lib.resultDir > self.Lib.queryName+"_not_called.tsv", "w")
			noCoverage = open(self.Lib.resultDir > self.Lib.queryName+"_no_coverage.tsv", "w")
			unique = open(self.Lib.resultDir > self.Lib.queryName+"_unique.tsv", "w")
		else:
			called = open(where > self.Lib.queryName+"_snps.tsv", "w")
			notCalled = open(where > self.Lib.queryName+"_not_called.tsv", "w")
			noCoverage = open(where > self.Lib.queryName+"_no_coverage.tsv", "w")
			unique = open(where > self.Lib.queryName+"_unique.tsv", "w")

		header = "\t".join(["Name","Reference","Pos","Ancestral base","Derived base", "Target base\n"])
		entryTemplate = "{name}\t{reference}\t{pos}\t{ancestral}\t{derived}\t{target}\n"

		called.write(header)
		notCalled.write(header)
		noCoverage.write(header)
		unique.write(header)

		for genomeID, genome, genbankID, refseqID, assemblyName in self.database.references:
			'''Print SNPs to tab separated file'''
			for snpID, position, ancestral, derived in self.database.SNPsByGenome[genome]:
				N : str = self.SNPresults[snpID]
				entry = entryTemplate.format(name=snpID, reference=genome, pos=position,
				 							 ancestral=ancestral, derived=derived, target=N)
				if derived == N:
					called.write(entry)
				elif ancestral == N:
					notCalled.write(entry)
				elif N.isalpha():
					unique.write(entry)
				else:
					noCoverage.write(entry)
		called.close()
		notCalled.close()
		noCoverage.close()
		unique.close()

	def readQueriesFrom(self, queryFile : str):
		'''If query input is a text file, parse file'''
		with open(queryFile, "r") as f:
			query = [query.strip() for query in f if len(query.strip())]
		return query

	def cleanup(self):
		'''Remove files in temporary folder'''
		LOGGER.info("Clean up temporary files... ")
		del self.Lib
		LOGGER.info("Done!")
		return
	