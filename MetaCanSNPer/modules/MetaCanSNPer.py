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

class Error(Exception):
	"""docstring for Error"""
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)


class MetaCanSNPerError(Error):
	"""docstring for MauveE"""
	pass

class MetaCanSNPer:
	outputTemplate = "{ref}_{target}.{format}"
	databaseName : str
	database : DatabaseReader
	Lib : DirectoryLibrary[str]
	settings : dict
	sessionName : str

	"""docstring for MetaCanSNPer"""
	def __init__(self, lib : DirectoryLibrary=None, tmpDir=None, refDir="References", database="CanSNPer.fdb", settingsFile : str=None, sessionName : str=None, **kwargs):
		self.startTime = time.localtime()
		
		self.setDatabase(database)
		if lib is None:
			self.Lib = DirectoryLibrary(tmpDir=tmpDir, refDir=os.path.join(refDir, os.path.splitext(database)[0]))
		else:
			self.Lib = lib
		settings : dict = toml.load(open(self.Lib.get("defaultFlags.toml", "installDir") if settingsFile is None else self.Lib.get(settingsFile, "workDir")), "rb")
		# Settings hierarchy looks like this: ["Category"]["Flag"] -> Value

		# Flatten hierarchy so flags are easily accessible
		self.settings = {}
		for flags in settings.values():
			for flag, value in flags.items():
				if flag in kwargs:
					self.settings[flag] = kwargs[flag]
				else:
					self.settings[flag] = value

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

	def connectDatabase(self):
		self.database = DatabaseReader(self.databaseName)

	'''MetaCanSNPer set directories'''
	
	def setTargetDir(self, path : str):		self.Lib.setTargetDir(path)
	def setRefDir(self, path : str):		self.Lib.setRefDir(path)
	def setDatabaseDir(self, path : str):	self.Lib.setDatabaseDir(path)
	def setTmpDir(self, path : str):		self.Lib.setTmpDir(path)
	def setOutDir(self, path : str):		self.Lib.setOutDir(path)

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

	'''MetaCanSNPer get functions'''

	def getReferences(self):
		'''Fetch names of reference genomes in connected database. Download any reference genomes not present locally.'''
		return self.database.references

	def getQuery(self):
		return self.Lib.query

	'''Indexing methods'''

	def createIndex(self, softwareName : str, kwargs : dict={}):
		'''Align sequences using subprocesses.'''

		LOGGER.debug("Checking for Aligner named '{}'".format( softwareName))
		indexerType : type = Aligners.get(softwareName) # Fetch the object class of the specified aligner software.
		if indexerType is None:
			LOGGER.debug("Checking for Mapper named '{}'".format( softwareName))
			indexerType : type = Mappers.get(softwareName) # Fetch the object class of the specified Mapper software.
		else:
			# No Aligner or Mapper found that is implemented yet.
			LOGGER.error("No software defined for name '{}'".format( softwareName))
			raise NotImplementedError("No software defined for name '{}'".format( softwareName))
		
		LOGGER.info("Checking if all references in the database have been downloaded. Downloads them if they haven't.")

		LOGGER.info("Creating Indexer '{}' of type '{}'".format( softwareName, indexerType.__name__))
		indexer : IndexingWrapper= indexerType(self.Lib, self.database, self.outputTemplate, kwargs=kwargs)

		LOGGER.info("Checking for pre-processing for indexer.")
		indexer.preProcess()

		LOGGER.info("Starting indexing.")
		output = indexer.start()

		# Start aligning
		while not indexer.finished():
			finished = indexer.waitNext()
			for i in finished:
				key, path = output[i]
				if key not in self.Lib.indexed and indexer.returncodes[i][-1] == 0:
					self.Lib.indexed[key] = path

		# Check that error did not occur.
		while indexer.hickups():
			if not indexer.fixable():
				# Error does not have a known solution that has worked during this run.
				return []
			else:
				# Error has a possible fix, implement the fix.
				indexer.planB()
			
			# Run the alignment just like before, but hopefully fixed.
			while not indexer.finished():
				finished = indexer.waitNext()
				for i in finished:
					key, path = output[i][0]
					if indexer.returncodes[i][-1] == 0:
						self.Lib.indexed[key] = path
	

	def callSNPs(self, softwareName : str, kwargs : dict={}):
		''''''
		SNPCallerType : SNPCaller = SNPCallers.get(softwareName)

		LOGGER.info("Loading SNPs from database.")
		references = {}
		for genome, reference in self.Lib.getReferences().items():
			SNPs, snpList = self.database.get_snps(reference=genome)
			references[genome] = SNPs
			LOGGER.info("Loaded {n} SNPs for genome {genome}.".format(n=len(snpList), genome=genome))
		LOGGER.info("Loaded a total of {n} SNPs.".format(n=sum(len(SNPs) for SNPs in references.values())))
		
		LOGGER.info("Creating SNPCaller '{}' of type '{}'".format( softwareName, SNPCallerType.__name__))
		snpCaller : SNPCaller = SNPCallerType(self.Lib, self.database, self.outputTemplate, kwargs=kwargs)
		
		LOGGER.info("Checking for pre-processing for SNPCaller.")
		snpCaller.preProcess()

		# Start aligning
		output = snpCaller.start()


		while not snpCaller.finished():
			finished = snpCaller.waitNext()
			for i in finished:
				key, path = output[i]
				if key not in self.Lib.SNPs and snpCaller.returncodes[i][-1] == 0:
					self.Lib.SNPs[key] = path

		# Check that error did not occur.
		while snpCaller.hickups():
			if not snpCaller.fixable():
				# Error does not have a known solution that has worked during this run.
				return []
			else:
				# Error has a possible fix, implement the fix.
				snpCaller.planB()
			
			# Run the alignment just like before, but hopefully fixed.
			while not snpCaller.finished():
				finished = snpCaller.waitNext()
				for i in finished:
					key, path = output[i][0]
					if snpCaller.returncodes[i][-1] == 0:
						self.Lib.SNPs[key] = path
		
		resultData = self.Lib.getSNPdata()

		# Interpret the SNP calls

		self.SNPresults = {}
		for snpID, (pos, anc, der) in self.database.SNPsByID.items():
			self.SNPresults[snpID] = resultData[pos]

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

	def saveResults(self, session : str=None):
		path = self.traverseTree()

		pathFile = open(os.path.join(self.Lib.resultDir, self.Lib.queryName+"_path.tsv"), "w")
		finalFile = open(os.path.join(self.Lib.resultDir, self.Lib.queryName+"_final.tsv"), "w")

		pathFile.write("\t".join(path))
		finalFile.write(path[-1])
		
		pathFile.close()
		finalFile.close()


	def saveSNPdata(self, session : str=None):
		""""""
		called = open(os.path.join(self.Lib.resultDir, self.Lib.queryName+"_snps.tsv"), "w")
		notCalled = open(os.path.join(self.Lib.resultDir, self.Lib.queryName+"_not_called.tsv"), "w")
		noCoverage = open(os.path.join(self.Lib.resultDir, self.Lib.queryName+"_no_coverage.tsv"), "w")
		unique = open(os.path.join(self.Lib.resultDir, self.Lib.queryName+"_unique.tsv"), "w")

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
	