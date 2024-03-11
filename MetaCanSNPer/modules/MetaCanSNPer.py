'''
MetaCanSNPer module: A toolkit for SNP-typing using NGS data.
Copyright (C) 2024 Fredrik Sörensen @ Umeå University
'''

import os, time
import tomllib as toml
from VariantCallFixer.Functions import getSNPdata

## import MetaCanSNPer specific modules
from MetaCanSNPer.Globals import *
import MetaCanSNPer.Globals as Globals
import MetaCanSNPer.modules.LogKeeper as LogKeeper
from MetaCanSNPer.modules.Databases import DatabaseReader
from MetaCanSNPer.modules.DirectoryLibrary import DirectoryLibrary
from MetaCanSNPer.modules.Wrappers import Aligner, Mapper, SNPCaller, IndexingWrapper
import MetaCanSNPer.modules.Aligners as Aligners
import MetaCanSNPer.modules.Mappers as Mappers
import MetaCanSNPer.modules.SNPCallers as SNPCallers

LOGGER = LogKeeper.createLogger(__name__)

def loadFlattenedTOML(filename):
	with open(filename, "rb") as f:
		tmp : dict[str,dict] = toml.load(f)
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
	SNPresults : dict

	"""docstring for MetaCanSNPer"""
	def __init__(self, lib : DirectoryLibrary=None, database=None, settings : dict={}, settingsFile : str=None, sessionName : str=None):

		LOGGER.info("Initializing MetaCanSNPer object.")
		self.startTime = time.localtime()
		self.sessionName = sessionName

		if lib is None:
			self.Lib = DirectoryLibrary(settings={k:(v if type(v) is not list else tuple(v)) for k,v in settings.items()})
		else:
			self.Lib = lib

		self.settings = self.Lib.settings

		if settingsFile is not None:
			if pIsAbs(settingsFile):
				settingsFlag = loadFlattenedTOML(settingsFile)
			else:
				settingsFlag = loadFlattenedTOML(self.Lib.targetDir[settingsFile])
			self.Lib.updateSettings({flag:(settingsFlag[flag] if type(settingsFlag[flag]) is not list else tuple(settingsFlag[flag])) for flag in set(settingsFlag).difference(self.Lib.settings)})
		if pExists(self.Lib.installDir > "defaultFlags.toml"):
			defaultFlags = loadFlattenedTOML(self.Lib.installDir > "defaultFlags.toml")
			self.Lib.updateSettings({flag:(defaultFlags[flag] if type(defaultFlags[flag]) is not list else tuple(defaultFlags[flag])) for flag in set(defaultFlags).difference(self.Lib.settings)})
		else:
			with open(self.Lib.installDir > "defaultFlags.toml", "w") as f:
				f.write(DEFAULT_TOML_TEMPLATE)

		if database is not None:
			self.setDatabase(database=database)
		
		self.SNPresults = {}

	'''Database functions'''

	def setDatabase(self, database : str):
		LOGGER.debug(f"Setting database to:{database}")
		if (path := self.Lib.databaseDir[database]) is not None:
			self.databaseName = path
			self.Lib.updateSettings({"organism":pName(database)})
			self.connectDatabase()
			self.setReferences(self.database.references)
		else:
			LOGGER.warning(f"Database not found: '{database}'")
			raise FileNotFoundError(f"Database not found: '{database}'")

	def connectDatabase(self):
		LOGGER.debug(f"Connecting to database:{self.databaseName}")
		if self.databaseName is None:
			LOGGER.error("Database not specified. Can't connect unless a valid database is set through MetaCanSNPer.setDatabase.")
			exit(1)
		try:
			self.database = DatabaseReader(self.databaseName)
			return 0
		except Exception as e:
			LOGGER.error(e)
			print(e)
			exit(1)

	'''MetaCanSNPer set directories'''
	
	def setTargetDir(self, path : str):						self.Lib.setTargetDir(path)
	def setRefDir(self, organism : str, path : str):		self.Lib.setRefDir(organism, path)
	def setDatabaseDir(self, path : str):					self.Lib.setDatabaseDir(path)
	def setTmpDir(self, path : str):						self.Lib.setTmpDir(path)
	def setOutDir(self, path : str):						self.Lib.setOutDir(path)

	# Setting the same docstring for each 'set*Dir' function
	setTargetDir.__doc__ = setRefDir.__doc__ = setDatabaseDir.__doc__ = setTmpDir.__doc__ = setOutDir.__doc__ = """
	`path` can be relative or absolute. Check the docstring for the DirectoryLibrary to see where the relative
	path will start from.
	"""

	'''MetaCanSNPer set values'''

	def setQuery(self, query : list[str]):
		self.Lib.setQuery(query)
		self.queryName = self.Lib.queryName ## get name of file and remove ending
		if self.sessionName is None:
			self.Lib.setSessionName("Sample-{queryName}-{dateYYYYMMDD}".format(queryName=self.queryName, dateYYYYMMDD="{:0>4}.{:0>2}.{:0>2}.{:>2}{:>2},{:>2}".format(*(self.startTime[:6]))))

	def setSessionName(self, name):
		self.sessionName = name
		self.Lib.setSessionName(name)

	def setReferences(self, references : list[str,str,str,str,str]):
		self.Lib.setReferences(references=references)

	'''MetaCanSNPer get functions'''

	def getReferences(self):
		'''Fetch names of reference genomes in connected database. Download any reference genomes not present locally.'''
		LOGGER.debug(f"Fetching references from database.")
		return self.database.references

	def getQuery(self):
		return self.Lib.query

	'''Indexing methods'''

	def getSoftwareClass(self, softwareName):
		for types in [Aligners, Mappers, SNPCallers]:
			if types.get(softwareName) is not None:
				return types.get(softwareName)
		return None

	def runSoftware(self, softwareClass : IndexingWrapper, outputDict : dict={}, flags : list=[]):
		'''Align sequences using subprocesses.'''

		LOGGER.info("Creating software wrapper for '{}' of type '{}'".format(softwareClass.softwareName, softwareClass.__name__))
		software : IndexingWrapper = softwareClass(self.Lib, self.database, self.outputTemplate, flags=flags)

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
				msg = [
					f"{softwareClass.softwareName!r}: Processes finished with exit codes for which there are no implemented solutions.",
					f"{'QUERY':<58}{'REFERENCE':<58}={'EXITCODE':>4}"
				]
				for (key, path), e in zip(output, software.returncodes):
					del outputDict[key]
					msg.append(f"{self.queryName!r:<30}{key!r:<60}={e[-1]:>4}")
				LOGGER.error("\n".join(msg))
				raise ChildProcessError("\n".join(msg))
		
		return outputDict
	
	def createMap(self, softwareName : str, flags : list=[]):
		''''''
		LOGGER.info(f"Creating map using:{softwareName}")
		MapperType : Mapper = Mappers.get(softwareName)

		LOGGER.info("Loading References from database.")
		self.database.references
		LOGGER.info("Loaded a total of {n} References.".format(n=len(self.database.references)))
		
		self.runSoftware(MapperType, outputDict=self.Lib.maps, flags=flags)

	def createAlignment(self, softwareName : str, flags : list=[]):
		''''''
		LOGGER.info(f"Creating alignment using:{softwareName}")
		AlignerType : Aligner = Aligners.get(softwareName)

		LOGGER.info("Loading References from database.")
		self.database.references
		LOGGER.info("Loaded a total of {n} References.".format(n=len(self.database.references)))
		
		self.runSoftware(AlignerType, outputDict=self.Lib.alignments, flags=flags)

	def callSNPs(self, softwareName : str, flags : list=[]):
		''''''
		LOGGER.info(f"Calling SNPs using:{softwareName}")
		SNPCallerType : SNPCaller = SNPCallers.get(softwareName)

		LOGGER.info("Loading References from database.")
		self.database.references
		LOGGER.info("Loading SNPs from database.")
		self.database.SNPsByGenome
		LOGGER.info("Loaded a total of {n} SNPs.".format(n=sum(len(SNPs) for SNPs in self.database.SNPsByGenome.values())))
		
		self.runSoftware(SNPCallerType, outputDict=self.Lib.resultSNPs, flags=flags)

		for genome, filePath in self.Lib.resultSNPs:
			getSNPdata(filePath, out=self.SNPresults)
	
	def traverseTree(self):
		'''Depth-first tree search.'''
		LOGGER.info(f"Traversing tree to get genotype called.")
		node = self.database.tree
		path = [node.nodeID]
		
		for child in node.children:
			if child.nodeID == 2: continue
			childSNPID = self.database.nodes[child.nodeID]
			pos, anc, der = self.database.SNPsByID[childSNPID]
			if der == self.SNPresults[pos]:
				node = child
				break

		while node.nodeID not in path:
			path.append(node.nodeID)
			for child in node.children:
				childSNPID = self.database.nodes[child.nodeID]
				pos, anc, der = self.database.SNPsByID[childSNPID]
				if der == self.SNPresults[pos]:
					node = child
					break
		
		LOGGER.info(f"Traversed:{path}")
		return path

	'''Functions'''

	def saveResults(self, dst : str=None):
		LOGGER.debug(f"open({dst or self.Lib.resultDir!r} > {self.Lib.queryName!r}+'_path.tsv', 'w')")
		LOGGER.debug(f"open({dst or self.Lib.resultDir!r} > {self.Lib.queryName!r}+'_final.tsv', 'w')")
		pathFile = open((dst or self.Lib.resultDir) > self.Lib.queryName+"_path.tsv", "w")
		finalFile = open((dst or self.Lib.resultDir) > self.Lib.queryName+"_final.tsv", "w")

		path = self.traverseTree()

		pathFile.write("\t".join(path))
		finalFile.write(path[-1])
		
		pathFile.close()
		finalFile.close()


	def saveSNPdata(self, where : str=None):
		""""""
		LOGGER.debug(f"open({where or self.Lib.resultDir!r} > {self.Lib.queryName!r}+'_snps.tsv', 'w')")
		LOGGER.debug(f"open({where or self.Lib.resultDir!r} > {self.Lib.queryName!r}+'_not_called.tsv', 'w')")
		LOGGER.debug(f"open({where or self.Lib.resultDir!r} > {self.Lib.queryName!r}+'_no_coverage.tsv', 'w')")
		LOGGER.debug(f"open({where or self.Lib.resultDir!r} > {self.Lib.queryName!r}+'_unique.tsv', 'w')")
		called = open((where or self.Lib.resultDir) > self.Lib.queryName+"_snps.tsv", "w")
		notCalled = open((where or self.Lib.resultDir) > self.Lib.queryName+"_not_called.tsv", "w")
		noCoverage = open((where or self.Lib.resultDir) > self.Lib.queryName+"_no_coverage.tsv", "w")
		unique = open((where or self.Lib.resultDir) > self.Lib.queryName+"_unique.tsv", "w")

		header = "\t".join(["Name","Reference","Pos","Ancestral base","Derived base", "Target base\n"])
		entryTemplate = "{name}\t{reference}\t{pos}\t{ancestral}\t{derived}\t{target}\n"

		called.write(header)
		notCalled.write(header)
		noCoverage.write(header)
		unique.write(header)

		for genomeID, genome, genbankID, refseqID, assemblyName in self.database.references:
			'''Print SNPs to tab separated file'''
			for snpID, position, ancestral, derived in self.database.SNPsByGenome[genome]:
				N, *args = self.SNPresults[snpID]
				entry = "\t".join([snpID, genome, position, ancestral, derived, N]) + "\n"
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
		LOGGER.debug(f"Reading queries from: {queryFile!r}")
		with open(queryFile, "r") as f:
			query = [query.strip() for query in f if len(query.strip())]
		return query
	
	def __del__(self):
		Globals.RUNNING = False

	def cleanup(self):
		'''Remove files in temporary folder'''
		LOGGER.info("Clean up temporary files... ")
		del self.Lib
		LOGGER.info("Done!")
		return
	