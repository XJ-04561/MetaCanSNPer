'''
MetaCanSNPer module: A toolkit for SNP-typing using NGS data.
Copyright (C) 2024 Fredrik Sörensen @ Umeå University
'''

import os, time
import tomllib as toml
from VariantCallFixer.Functions import getSNPdata
from PseudoPathy import Path

## import MetaCanSNPer specific modules
from MetaCanSNPer.Globals import *
import MetaCanSNPer.Globals as Globals
import MetaCanSNPer.modules.LogKeeper as LogKeeper
from MetaCanSNPer.modules.DirectoryLibrary import DirectoryLibrary
from MetaCanSNPer.modules.Hooks import Hooks
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
	outputTemplate = "{refName}_{queryName}.{outFormat}"
	databasePath : Path
	databaseName : str
	database : DatabaseReader
	Lib : DirectoryLibrary
	settings : dict
	sessionName : str
	SNPresults : dict
	exceptions : list[Exception]

	"""docstring for MetaCanSNPer"""
	def __init__(self, lib : DirectoryLibrary=None, database : str=None, settings : dict={}, settingsFile : str=None, sessionName : str=None):

		LOGGER.info("Initializing MetaCanSNPer object.")
		self.startTime = time.localtime()
		self.sessionName = sessionName
		self.exceptions = []
		self.hooks = Hooks()
		self.hooks.addHook("ReportError", target=lambda eventInfo : self.exceptions.append(eventInfo["exception"]))
		
		if lib is None:
			self.Lib = DirectoryLibrary(settings={k:(v if type(v) is not list else tuple(v)) for k,v in settings.items()})
		else:
			self.Lib = lib
		
		self.settings = self.Lib.settings

		if settingsFile is not None:
			if pIsAbs(settingsFile):
				settingsFlag = loadFlattenedTOML(settingsFile)
			elif self.Lib.targetDir[settingsFile] is not None:
				settingsFlag = loadFlattenedTOML(self.Lib.targetDir[settingsFile])
			else:
				raise FileNotFoundError(f"Could not found the settingsFile {settingsFile!r} in any of {self.Lib.targetDir}")
			self.Lib.updateSettings({flag:(settingsFlag[flag] if type(settingsFlag[flag]) is not list else tuple(settingsFlag[flag])) for flag in set(settingsFlag).difference(self.Lib.settings)})
		if "defaultFlags.toml" in self.Lib.commonGroups.shared:
			defaultFlags = loadFlattenedTOML(self.Lib.commonGroups.shared.find("defaultFlags.toml", "r"))
			self.Lib.updateSettings({flag:(defaultFlags[flag] if type(defaultFlags[flag]) is not list else tuple(defaultFlags[flag])) for flag in set(defaultFlags).difference(self.Lib.settings)})
		else:
			if self.Lib.commonGroups.shared.writable is None:
				LOGGER.error(f"Could not write to any of the directories in {self.Lib.commonGroups.shared}")
				raise PermissionError("Could not find a directory to write program-essential files to.")
			with open(self.Lib.commonGroups.shared.writable / "defaultFlags.toml", "w") as f:
				f.write(DEFAULT_TOML_TEMPLATE)

		self.databasePath = None
		self.databaseName = None

		if database is not None:
			self.setDatabase(database=database)
		
		self.SNPresults = {}

	'''Database functions'''

	def setDatabase(self, database : str, silent : bool=False):
		LOGGER.debug(f"Setting database to:{database}")
		self.databaseName = os.path.basename(database)
		if (path := self.Lib.databaseDir.find(database, purpose="r")) is not None:
			LOGGER.info(f"Found database {database!r} in path {path!r}")
		else:
			if not silent: print(f"Downloading database {self.databaseName} ... ", end="", flush=True)
			if (path := downloadDatabase(self.databaseName, dst=self.Lib.databaseDir.forceFind("", "w") / database)) is None:
				if not silent: print("Failed!", flush=True)
				LOGGER.error(f"Database not found locally or online: {database!r}\nLocal directories checked: {self.Lib.databaseDir}")
				raise FileNotFoundError(f"Database not found: {database!r}")
			else:
				if not silent: print("Done!", flush=True)
			LOGGER.info(f"Found database {database!r} online and downloaded to path {path!r}")
		self.database : DatabaseReader = DatabaseReader(path)
		try:
			self.database.validateDatabase(self.database.checkDatabase())
		except:
			LOGGER.info(f"Database {path!r} is not up to date/does not have correct schema.")
			self.database.close()
			if (path := self.Lib.databaseDir.find(self.databaseName, purpose="w")) is None:
				if (path := self.Lib.databaseDir.find(purpose="w")) is not None:
					CanSNPDB.Commands.download(databaseName=self.databaseName, outDir=path)
				else:
					raise PermissionError(f"Can't find a valid database: {database!r} or find a writeable directory in which to create one.")
			self.database : DatabaseWriter = DatabaseWriter(path / self.databaseName)
			code = self.database.checkDatabase()
			try:
				self.database.validateDatabase(code)
			except IsLegacyCanSNPer2:
				self.Lib.setReferences(self.database.ReferenceTable.get(DB.ALL))
			except:
				pass
			self.database.rectifyDatabase(code)
			self.database.validateDatabase(code)
			self.database.commit()
			self.database.close()
			self.database : DatabaseReader = DatabaseReader(path)


		self.databasePath = path / self.databaseName
		self.Lib.updateSettings({"organism":pName(self.databaseName)}) # TODO : Get organism name from safer source than filename
		self.Lib.references = None
		self.connectDatabase()

	def connectDatabase(self):
		LOGGER.debug(f"Connecting to database:{self.databasePath}")
		if self.databasePath is None:
			raise NameError("Database not specified. Can't connect unless a valid database is set through MetaCanSNPer.setDatabase.")
		
		self.database = DatabaseReader(self.databasePath)

		LOGGER.debug(f"Connected to database:{self.databasePath}")
		return 0

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
		LOGGER.debug(f"Query in DirectoryLibrary: {self.Lib.query}")
		self.queryName = self.Lib.queryName ## get name of file and remove ending
		if self.sessionName is None:
			self.setSessionName("Sample-{queryName}-{dateYYYYMMDD}".format(queryName=self.queryName, dateYYYYMMDD="{:0>4}-{:0>2}-{:0>2}-{:0>2}.{:0>2}.{:0<3}".format(*(self.startTime[:6]))))

	def setSessionName(self, name):
		self.sessionName = name
		self.Lib.setSessionName(name)

	def setReferenceFiles(self, references : list[str,str,str,str,str]=None, silent : bool=False):
		self.Lib.setReferences(references=references or self.database.references, silent=silent)

	'''MetaCanSNPer get functions'''

	def getReferenceFiles(self) -> dict[str,Path]:
		return self.Lib.references

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

		LOGGER.info(f"Creating software wrapper for {softwareClass.softwareName!r} of type {softwareClass.__name__!r}")
		software : IndexingWrapper = softwareClass(self.Lib, self.database, self.outputTemplate, out=outputDict, hooks=self.hooks, flags=flags, settings=self.settings)

		# Check that error did not occur.
		while software.canRun():
			LOGGER.info(f"Checking for pre-processing for {software.softwareName}.")
			software.preProcess()

			LOGGER.info(f"Starting {software.softwareName}.")
			software.start()
			
			LOGGER.info(f"Waiting for {software.softwareName} to finish.")
			software.wait()
			LOGGER.info(f"Finished running all instances of {software.softwareName}.")

			if software.hickups():
				if software.fixable():
					LOGGER.info(f"{software.softwareName} failed. Fix for the specific exitcodes exist. Implementing and running again.")
					software.planB()
				else:
					software.displayOutcomes(out=LOGGER.error)
					raise ChildProcessError(f"{software.softwareName} returned with a non-zero exitcode for which there is no implemented solution.")
		
		return outputDict
	
	def createMap(self, softwareName : str, flags : list=[]):
		''''''
		LOGGER.info(f"Creating map using:{softwareName}")
		MapperType : Mapper = Mappers.get(softwareName)

		if self.Lib.references is None:
			LOGGER.error("References not set.")
			raise FileNotFoundError("References not set. Can be set with MetaCanSNPer.setReferences")
		
		LOGGER.info(f"Loaded a total of {sum(1 for r in self.database.references)} References.")
		
		self.runSoftware(MapperType, outputDict=self.Lib.maps, flags=flags)

		LOGGER.info(f"Result of Mapping in: {self.Lib.maps}")

	def createAlignment(self, softwareName : str, flags : list=[]):
		''''''
		LOGGER.info(f"Creating alignment using:{softwareName}")
		AlignerType : Aligner = Aligners.get(softwareName)

		if self.Lib.references is None:
			LOGGER.error("References not set.")
			raise FileNotFoundError("References not set. Can be set with MetaCanSNPer.setReferences")
		
		LOGGER.info(f"Loaded a total of {sum(1 for r in self.database.references)} References.")
		
		self.runSoftware(AlignerType, outputDict=self.Lib.alignments, flags=flags)

		LOGGER.info(f"Result of Alignment in: {self.Lib.alignments}")

	def callSNPs(self, softwareName : str, flags : list=[]):
		''''''
		LOGGER.info(f"Calling SNPs using:{softwareName}")
		SNPCallerType : SNPCaller = SNPCallers.get(softwareName)

		if self.Lib.references is None:
			LOGGER.error("References not set.")
			raise FileNotFoundError("References not set. Can be set with MetaCanSNPer.setReferences")
		
		LOGGER.info("Loading SNPs from database.")
		self.Lib.setTargetSNPs()
		LOGGER.info(f"Loaded a total of {len(self.database.SNPTable)} SNPs.")
		
		self.runSoftware(SNPCallerType, outputDict=self.Lib.resultSNPs, flags=flags)

		LOGGER.info(f"Result of SNPCalling in: {self.Lib.resultSNPs}")

		for genome, filePath in self.Lib.resultSNPs:
			for pos, (chromosome, ref) in getSNPdata(filePath, values=["CHROM", "REF"]):
				(nodeID,) = self.database.SNPTable.first(DB.NodeID, Position=pos, Chromosome=chromosome)
				if (pos, genome) not in self.SNPresults:
					self.SNPresults[nodeID] = {}
				self.SNPresults[nodeID][pos] = ref
	
	def traverseTree(self):
		'''Depth-first tree search.'''
		LOGGER.info(f"Traversing tree to get genotype called.")
		node = self.database.tree
		paths : list[list[Branch]] = []
		nodeScores = {node.nodeID:0}

		award = (1, -1, 0)

		paths.append([node])
		while paths[-1] != []:
			paths.append([])
			for node in paths[-2]:
				for child in node.children:
					if child.nodeID in nodeScores: continue
					for childSNPID, (pos, anc, der) in self.database.SNPsByNode(child.nodeID):
						nodeScores[child.nodeID] = nodeScores[node.nodeID]
						(called, *_) = self.SNPresults[childSNPID]
						if der == called:
							nodeScores[child.nodeID] += award[0]
						elif anc == called:
							nodeScores[child.nodeID] += award[1]
						else:
							nodeScores[child.nodeID] += award[2]
					paths[-1].append(child)
			if paths[-1] == []:
				paths = paths[:-1]
				LOGGER.info(f"Finished traversed tree.")
				return max(nodeScores.items(), key=lambda nodeTupe: nodeTupe[1]), nodeScores

	'''Functions'''

	def saveResults(self, dst : str=None):
		LOGGER.debug(f"open({dst or self.Lib.resultDir.writable!r} / {self.Lib.queryName!r}+'_final.tsv', 'w')")
		with open((dst or self.Lib.resultDir.writable) / self.Lib.queryName+"_final.tsv", "w") as finalFile:
			(finalNodeID, score), scores = self.traverseTree()
			
			finalFile.write(f"{self.database.nodeName(finalNodeID):<20}{score}\n")
			if self.settings.get("debug"):
				finalFile.write("\n")
				for nodeID in scores:
					finalFile.write(f"{self.database.nodeName(nodeID):<20}{scores[nodeID]}\n")

	def saveSNPdata(self, dst : str=None):
		""""""

		header = "Name\tReference\tPos\tAncestral base\tDerived base\tTarget base\n"

		if self.settings.get("debug"):
			LOGGER.debug(f"open({dst or self.Lib.resultDir.writable!r} / {self.Lib.queryName!r}+'_snps.tsv', 'w')")
			LOGGER.debug(f"open({dst or self.Lib.resultDir.writable!r} / {self.Lib.queryName!r}+'_not_called.tsv', 'w')")
			LOGGER.debug(f"open({dst or self.Lib.resultDir.writable!r} / {self.Lib.queryName!r}+'_no_coverage.tsv', 'w')")
			LOGGER.debug(f"open({dst or self.Lib.resultDir.writable!r} / {self.Lib.queryName!r}+'_unique.tsv', 'w')")
			called = open((dst or self.Lib.resultDir.writable) / self.Lib.queryName+"_snps.tsv", "w")
			notCalled = open((dst or self.Lib.resultDir.writable) / self.Lib.queryName+"_not_called.tsv", "w")
			noCoverage = open((dst or self.Lib.resultDir.writable) / self.Lib.queryName+"_no_coverage.tsv", "w")
			unique = open((dst or self.Lib.resultDir.writable) / self.Lib.queryName+"_unique.tsv", "w")
			
			called.write(header)
			notCalled.write(header)
			noCoverage.write(header)
			unique.write(header)
		else:
			LOGGER.debug(f"open({dst or self.Lib.resultDir.writable!r} / {self.Lib.queryName!r}+'_snps.tsv', 'w')")
			LOGGER.debug(f"open({dst or self.Lib.resultDir.writable!r} / {self.Lib.queryName!r}+'_not_called.tsv', 'w')")
			called = open((dst or self.Lib.resultDir.writable) / self.Lib.queryName+"_snps.tsv", "w")
			noCoverage = unique = notCalled = open((dst or self.Lib.resultDir.writable) / self.Lib.queryName+"_not_called.tsv", "w")
			
			called.write(header)
			noCoverage.write(header)

		for genomeID, genome, genbankID, refseqID, assemblyName in self.database.references:
			'''Print SNPs to tab separated file'''
			for snpID, position, ancestral, derived in self.database.SNPsByGenome[genome]:
				N, *args = self.SNPresults[snpID]
				entry = f"{snpID}\t{genome}\t{position}\t{ancestral}\t{derived}\t{N}\n"
				if derived == N:
					called.write(entry)
				elif ancestral == N:
					notCalled.write(entry)
				elif N.isalpha():
					unique.write(entry)
				else:
					noCoverage.write(entry)
		called.close()
		noCoverage.close()
		try:
			notCalled.close()
			unique.close()
		except:
			pass

	def readQueriesFrom(self, queryFile : str):
		'''If query input is a text file, parse file'''
		LOGGER.debug(f"Reading queries from: {queryFile!r}")
		with open(queryFile, "r") as f:
			query = [query.strip() for query in f if len(query.strip())]
		return query
	
	# def __del__(self):
	# 	Globals.RUNNING = False

	def cleanup(self):
		'''Remove files in temporary folder'''
		LOGGER.info("Clean up temporary files... ")
		del self.Lib
		LOGGER.info("Done!")
		return
	