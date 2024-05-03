'''
MetaCanSNPer module: A toolkit for SNP-typing using NGS data.
Copyright (C) 2024 Fredrik Sörensen @ Umeå University
'''

import tomllib as toml
from VariantCallFixer.Functions import getSNPdata

## import MetaCanSNPer specific modules
from MetaCanSNPer.Globals import *
import MetaCanSNPer.Globals as Globals
import MetaCanSNPer.core.LogKeeper as LogKeeper
from MetaCanSNPer.core.DirectoryLibrary import DirectoryLibrary
from MetaCanSNPer.core.Hooks import Hooks, DummyHooks
from MetaCanSNPer.core.Wrappers import Aligner, Mapper, SNPCaller, IndexingWrapper
import MetaCanSNPer.core.Aligners as Aligners
import MetaCanSNPer.core.Mappers as Mappers
import MetaCanSNPer.core.SNPCallers as SNPCallers
from MetaCanSNPer.core.TerminalUpdater import TerminalUpdater

from MetaCanSNPer.modules.Database import (MetaCanSNPerDatabase, DatabaseError, verifyDatabase, correctDatabase,
										   Parent, NodeID, GenoType, Position, Ancestral, Derived, SNPReference,
										   Date, ChromID, Chromosome, GenomeID, Genome, Strain, GenbankID,
										   RefseqID, Assembly,
										   Branch)
from MetaCanSNPer.modules.Downloader import DatabaseDownloader, DownloaderReportHook, ReferenceDownloader



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

class MetaCanSNPer(DescribeAnnotations):
	outputTemplate = "{refName}_{queryName}.{outFormat}"
	organism : str = NotSet
	databasePath : Path = NotSet
	databaseName : str = NotSet
	database : MetaCanSNPerDatabase = NotSet
	Lib : DirectoryLibrary
	settings : dict
	sessionName : str = NotSet
	SNPresults : dict
	exceptions : list[Exception]

	queryName : str

	"""docstring for MetaCanSNPer"""
	def __init__(self, lib : DirectoryLibrary=None, database : str=None, settings : dict={}, settingsFile : str=None, sessionName : str=None):

		LOGGER.info("Initializing MetaCanSNPer object.")
		self.startTime = time.localtime()
		self.sessionName = sessionName
		self.exceptions = []
		self.hooks = Hooks()
		self.hooks.addHook("ReportError", target=lambda eventInfo : self.exceptions.append(eventInfo["exception"]))
		
		if lib is None:
			self.Lib = DirectoryLibrary(settings={k:(v if type(v) is not list else tuple(v)) for k,v in settings.items()}, hooks=self.hooks)
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
			self.connectDatabase()
		
		self.SNPresults = {}

	'''Database functions'''

	def setOrganism(self, organism : str):
		self.Lib.updateSettings({"organism" : organism})
		del self.databaseName, self.database, self.databasePath
		
	def setDatabase(self, databaseName : str=organism + ".db"):

		self.Lib.references = None
		if databaseName is Globals._NOT_SET:
			databaseName = f"{self.organism}.db"
		LOGGER.debug(f"Setting database to:{databaseName}")
		self.databaseName = databaseName
		
		for directory in self.Lib.databaseDir:
			if databaseName in directory:
				if MetaCanSNPerDatabase(directory / databaseName, "r").valid is True:
					self.hooks.trigger("DatabaseDownloaderProgress", {"name" : databaseName, "progress" : int(1)})
					self.databasePath = directory / databaseName
					break
		else:
			LOGGER.info(f"No valid database={databaseName} found. Looking for writeable versions or directories that can be updated or downloaded to, respectively.")
			RH = DownloaderReportHook("DatabaseDownloader", self.hooks, databaseName)
			
			self.databasePath = self.Lib.databaseDir.writable / databaseName
			DD = DatabaseDownloader(self.Lib.databaseDir.writable)

			DD.download(databaseName, databaseName, reportHook=RH)
			
			DD.wait()
		
		self.database = MetaCanSNPerDatabase(self.databasePath, "r")
		LOGGER.info(f"Database {self.databaseName} loaded from: {self.databasePath}!")
	
	def setReferenceFiles(self, references : Iterable[tuple[int,str,str,str,str,str]]=database.references):

		if references is NotSet:
			assert self.database is not NotSet, "Database not yet set"
			references = self.database.references

		directory = self.Lib.refDir.writable
		DD = ReferenceDownloader(directory)

		self.Lib.references.clear()

		for genomeID, genome, strain, genbankID, refseqID, assemblyName in references:
			DD.download((genbankID, assemblyName), f"{assemblyName}.fna", reportHook=DownloaderReportHook("ReferenceDownloader", self.hooks, genome))
			
			self.Lib.references[genome] = directory / f"{assemblyName}.fna"
			self.Lib.references[assemblyName] = directory / f"{assemblyName}.fna"
		DD.wait()

	'''MetaCanSNPer set directories'''
	
	def setTargetDir(self, path : str):						self.Lib.setTargetDir(path)
	def setRefDir(self, organism : str, path : str):		self.Lib.setRefDir(organism, path)
	def setDatabaseDir(self, path : str):					self.Lib.setDatabaseDir(path)
	def setTmpDir(self, path : str):						self.Lib.setTmpDir(path)
	def setOutDir(self, path : str):						self.Lib.setOutDir(path)
	"""
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

	def setSessionName(self, name="Sample-" + queryName + "{dateYYYYMMDD}"):
		self.sessionName = name
		self.Lib.setSessionName(name)

	@property
	def queryName(self):
		return self.Lib.queryName

	'''MetaCanSNPer get functions'''

	@property
	def referenceFiles(self) -> dict[str,Path]:
		'''Fetch files of reference genomes.'''
		return self.Lib.references

	@property
	def references(self):
		'''Fetch names of reference genomes in connected database.'''
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
		LOGGER.info(f"Loaded a total of {len(self.database.SNPs)} SNPs.")
		
		self.runSoftware(SNPCallerType, outputDict=self.Lib.resultSNPs, flags=flags)

		LOGGER.info(f"Result of SNPCalling in: {self.Lib.resultSNPs}")

		for genome, filePath in self.Lib.resultSNPs:
			for pos, (chromosome, ref) in getSNPdata(filePath, values=["CHROM", "REF"]):
				(nodeID,) = self.database[NodeID][Position==pos, Chromosome==chromosome]
				if (pos, genome) not in self.SNPresults:
					self.SNPresults[nodeID] = {}
				self.SNPresults[nodeID][pos] = ref
	
	def traverseTree(self):
		'''Depth-first tree search.'''
		LOGGER.info(f"Traversing tree to get genotype called.")
		rootNode = self.database.tree
		paths : list[list[Branch]] = []
		nodeScores = {rootNode.node : 0}

		award = (1, -1, 0)

		paths.append([rootNode])
		while paths[-1] != []:
			paths.append([])
			for parent in paths[-2]:
				for child in parent.children:
					nodeScores[child.node] = nodeScores[parent.node]
					for childSNPID, pos, anc, der, *_ in self.database.SNPsByNode[child.node]:
						(called, *_) = self.SNPresults[childSNPID]
						if der == called:
							nodeScores[child.node] += award[0]
						elif anc == called:
							nodeScores[child.node] += award[1]
						else:
							nodeScores[child.node] += award[2]
					paths[-1].append(child)
		
		paths = paths[:-1]
		LOGGER.info(f"Finished traversing tree.")
		return max(nodeScores.items(), key=lambda x: x[1]), nodeScores

	'''Functions'''

	def saveResults(self, dst : str=None):
		LOGGER.debug(f"open({dst or self.Lib.resultDir.writable!r} / {self.Lib.queryName!r}+'_final.tsv', 'w')")
		with open((dst or self.Lib.resultDir.writable) / self.Lib.queryName+"_final.tsv", "w") as finalFile:
			(finalNodeID, score), scores = self.traverseTree()
			
			finalFile.write("{:<20}{score}\n".format(*self.database[GenoType][NodeID == finalNodeID], score=score))
			if self.settings.get("debug"):
				finalFile.write("\n")
				for nodeID in scores:
					finalFile.write("{:<20}{score}\n".format(*self.database[GenoType][NodeID == nodeID], score=scores[nodeID]))

	def saveSNPdata(self, dst : str=None):
		""""""

		header = "Name\tReference\tChromosome\tPosition\tAncestral base\tDerived base\tTarget base\n"

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
			for nodeID, position, ancestral, derived, chromosome in self.database[NodeID, Position, Ancestral, Derived, Chromosome][GenomeID == genomeID]:
				N, *args = self.SNPresults[position, chromosome]
				entry = f"{nodeID}\t{genome}\t{chromosome}\t{position}\t{ancestral}\t{derived}\t{N}\n"
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
	