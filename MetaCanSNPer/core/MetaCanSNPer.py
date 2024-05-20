'''
MetaCanSNPer module: A toolkit for SNP-typing using NGS data.
Copyright (C) 2024 Fredrik Sörensen @ Umeå University
'''

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

from MetaCanSNPer.modules.Database import *
from MetaCanSNPer.modules.Downloader import DatabaseDownloader, DownloaderReportHook, ReferenceDownloader


class MetaCanSNPer:
	
	LOG : logging.Logger = LOGGER
	hooks = Hooks()

	outputTemplate = "{refName}_{queryName}.{outFormat}"
	organism : str

	settings : dict = DEFAULT_SETTINGS.copy()
	Lib : DirectoryLibrary
	database : MetaCanSNPerDatabase = NotSet

	databasePath : Path
	databaseName : str
	
	referenceFiles = property(lambda self:self.Lib.references)
	references = property(lambda self:self.database.references, doc="""Fetch names of reference genomes in connected database.""")

	query = property(lambda self:self.Lib.query)
	queryName : str = Default["Lib.query"](lambda self:self.Lib.queryName)
	sessionName : str = Default["Lib.query"](lambda self:self.Lib.sessionName)
	
	SNPresults : dict = Default["Lib.query"](lambda self:dict())
	exceptions : list[Exception]


	@overload
	def __init__(self, /, organism : str, query : list[str]|str, *, lib : DirectoryLibrary=None, database : str=None, settings : dict={}, settingsFile : str=None, sessionName : str=None): ...

	def __init__(self, /, organism : str, query : list[str]|str, **kwargs):

		self.LOG.info(f"Initializing MetaCanSNPer object at 0x{id(self):0>16X}.")
		self.exceptions = []
		self.hooks.addHook("ReportError", target=lambda eventInfo : self.exceptions.append(eventInfo["exception"]))

		if "settingsFile" in kwargs and os.path.exists(kwargs["settingsFile"]):
			self.settings |= loadFlattenedTOML(kwargs["settingsFile"])
		elif "settingsFile" in kwargs:
			raise FileNotFoundError(f"Could not found the settingsFile {kwargs['settingsFile']!r}")
			
		self.settings |= kwargs.get("settings", {})

		self.Lib = DirectoryLibrary(organism, query, settings=kwargs.get("settings"), hooks=self.hooks)

		self.LOG = self.LOG.getChild(f"[{self.Lib.queryName}]")

		self.organism = organism

	def setOrganism(self, organism : str):
		self.Lib.organism = organism
		
	def setDatabase(self, databaseName : str=None):
		"""databaseName Defaults to `{organism}.db`"""

		if databaseName is None:
			self.databaseName = databaseName = f"{self.organism}.db"
		
		self.LOG.debug(f"Setting database to:{databaseName}")
		
		for directory in self.Lib.databaseDir:
			if databaseName in directory:
				if MetaCanSNPerDatabase(directory / databaseName, "r", organism=self.organism).valid is True:
					self.hooks.trigger("DatabaseDownloaderProgress", {"name" : databaseName, "progress" : int(1)})
					self.databasePath = directory / databaseName
					break
		else:
			self.LOG.info(f"No valid database={databaseName} found. Looking for writeable versions or directories that can be updated or downloaded to, respectively.")
			RH = DownloaderReportHook("DatabaseDownloader", self.hooks, databaseName)
			
			self.databasePath = self.Lib.databaseDir.writable / databaseName
			DD = DatabaseDownloader(self.Lib.databaseDir.writable)

			DD.download(databaseName, databaseName, reportHook=RH)
			
			DD.wait()
		
		self.database = MetaCanSNPerDatabase(self.databasePath, "r", organism=self.organism)
		self.LOG.info(f"Database {self.databaseName} loaded from: {self.databasePath}!")
	
	def setReferenceFiles(self, references : Iterable[tuple[int,str,str,str,str,str]]=database.references):

		if references is NotSet:
			assert hasattr(self, "database"), "Database not yet set. Setting references requires either a table of references as an argument or that a database is connected."
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
	
	def setTargetDir(self, path : str):		setattr(self.Lib, "setTargetDir", path)
	def setRefDir(self, path : str):		setattr(self.Lib, "setRefDir", path)
	def setDatabaseDir(self, path : str):	setattr(self.Lib, "setDatabaseDir", path)
	def setTmpDir(self, path : str):		setattr(self.Lib, "setTmpDir", path)
	def setOutDir(self, path : str):		setattr(self.Lib, "setOutDir", path)
	"""
	`path` can be relative or absolute. Check the docstring for the DirectoryLibrary to see where the relative
	path will start from.
	"""

	'''MetaCanSNPer set values'''

	def setSessionName(self, name):
		self.sessionName = name
		self.Lib.sessionName = name


	'''Indexing methods'''

	def getSoftwareClass(self, softwareName):
		for types in [Aligners, Mappers, SNPCallers]:
			if types.get(softwareName) is not None:
				return types.get(softwareName)
		return None

	def runSoftware(self, softwareClass : IndexingWrapper, outputDict : dict={}, flags : list=[]):
		'''Align sequences using subprocesses.'''

		self.LOG.info(f"Creating software wrapper for {softwareClass.softwareName!r} of type {softwareClass.__name__!r}")
		software : IndexingWrapper = softwareClass(self.Lib, self.database, self.outputTemplate, out=outputDict, hooks=self.hooks, flags=flags, settings=self.settings)

		# Check that error did not occur.
		while software.canRun():
			self.LOG.info(f"Checking for pre-processing for {software.softwareName}.")
			software.preProcess()

			self.LOG.info(f"Starting {software.softwareName}.")
			software.start()
			
			self.LOG.info(f"Waiting for {software.softwareName} to finish.")
			software.wait()
			self.LOG.info(f"Finished running all instances of {software.softwareName}.")

			if software.hickups():
				if software.fixable():
					self.LOG.info(f"{software.softwareName} failed. Fix for the specific exitcodes exist. Implementing and running again.")
					software.planB()
				else:
					software.displayOutcomes(out=self.LOG.error)
					raise ChildProcessError(f"{software.softwareName} returned with a non-zero exitcode for which there is no implemented solution.")
		
		return outputDict
	
	def createMap(self, softwareName : str, flags : list=[]):
		''''''
		self.LOG.info(f"Creating map using:{softwareName}")
		MapperType : Mapper = Mappers.get(softwareName)

		if self.Lib.references is None:
			self.LOG.error("References not set.")
			raise FileNotFoundError("References not set. Can be set with MetaCanSNPer.setReferences")
		
		self.LOG.info(f"Loaded a total of {sum(1 for r in self.database.references)} References.")
		
		self.runSoftware(MapperType, outputDict=self.Lib.maps, flags=flags)

		self.LOG.info(f"Result of Mapping in: {self.Lib.maps}")

	def createAlignment(self, softwareName : str, flags : list=[]):
		''''''
		self.LOG.info(f"Creating alignment using:{softwareName}")
		AlignerType : Aligner = Aligners.get(softwareName)

		if self.Lib.references is None:
			self.LOG.error("References not set.")
			raise FileNotFoundError("References not set. Can be set with MetaCanSNPer.setReferences")
		
		self.LOG.info(f"Loaded a total of {sum(1 for r in self.database.references)} References.")
		
		self.runSoftware(AlignerType, outputDict=self.Lib.alignments, flags=flags)

		self.LOG.info(f"Result of Alignment in: {self.Lib.alignments}")

	def callSNPs(self, softwareName : str, flags : list=[]):
		''''''
		self.LOG.info(f"Calling SNPs using:{softwareName}")
		SNPCallerType : SNPCaller = SNPCallers.get(softwareName)

		if self.Lib.references is None:
			self.LOG.error("References not set.")
			raise FileNotFoundError("References not set. Can be set with MetaCanSNPer.setReferences")
		
		self.LOG.info("Loading SNPs from database.")
		self.Lib.setTargetSNPs()
		self.LOG.info(f"Loaded a total of {len(self.database.SNPs)} SNPs.")
		
		self.runSoftware(SNPCallerType, outputDict=self.Lib.resultSNPs, flags=flags)

		self.LOG.info(f"Result of SNPCalling in: {self.Lib.resultSNPs}")

		for genome, filePath in self.Lib.resultSNPs:
			for pos, (chromosome, ref) in getSNPdata(filePath, values=["CHROM", "REF"]):
				(nodeID,) = self.database[NodeID, Position==pos, Chromosome==chromosome]
				if (pos, genome) not in self.SNPresults:
					self.SNPresults[nodeID] = {}
				self.SNPresults[nodeID][pos] = ref
	
	def traverseTree(self):
		'''Depth-first tree search.'''
		self.LOG.info(f"Traversing tree to get genotype called.")
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
		self.LOG.info(f"Finished traversing tree.")
		return max(nodeScores.items(), key=lambda x: x[1]), nodeScores

	'''Functions'''

	def saveResults(self, dst : str=None):
		self.LOG.debug(f"open({dst or self.Lib.resultDir.writable!r} / {self.Lib.queryName!r}+'_final.tsv', 'w')")
		with open((dst or self.Lib.resultDir.writable) / self.Lib.queryName+"_final.tsv", "w") as finalFile:
			(finalNodeID, score), scores = self.traverseTree()
			
			finalFile.write("{:<20}{score}\n".format(*self.database[GenoType, NodeID == finalNodeID], score=score))
			if self.settings.get("debug"):
				finalFile.write("\n")
				for nodeID in scores:
					finalFile.write("{:<20}{score}\n".format(*self.database[GenoType, NodeID == nodeID], score=scores[nodeID]))

	def saveSNPdata(self, dst : str=None):
		""""""

		header = "Name\tReference\tChromosome\tPosition\tAncestral base\tDerived base\tTarget base\n"

		if self.settings.get("debug"):
			self.LOG.debug(f"open({dst or self.Lib.resultDir.writable!r} / {self.Lib.queryName!r}+'_snps.tsv', 'w')")
			self.LOG.debug(f"open({dst or self.Lib.resultDir.writable!r} / {self.Lib.queryName!r}+'_not_called.tsv', 'w')")
			self.LOG.debug(f"open({dst or self.Lib.resultDir.writable!r} / {self.Lib.queryName!r}+'_no_coverage.tsv', 'w')")
			self.LOG.debug(f"open({dst or self.Lib.resultDir.writable!r} / {self.Lib.queryName!r}+'_unique.tsv', 'w')")
			called = open((dst or self.Lib.resultDir.writable) / self.Lib.queryName+"_snps.tsv", "w")
			notCalled = open((dst or self.Lib.resultDir.writable) / self.Lib.queryName+"_not_called.tsv", "w")
			noCoverage = open((dst or self.Lib.resultDir.writable) / self.Lib.queryName+"_no_coverage.tsv", "w")
			unique = open((dst or self.Lib.resultDir.writable) / self.Lib.queryName+"_unique.tsv", "w")
			
			called.write(header)
			notCalled.write(header)
			noCoverage.write(header)
			unique.write(header)
		else:
			self.LOG.debug(f"open({dst or self.Lib.resultDir.writable!r} / {self.Lib.queryName!r}+'_snps.tsv', 'w')")
			self.LOG.debug(f"open({dst or self.Lib.resultDir.writable!r} / {self.Lib.queryName!r}+'_not_called.tsv', 'w')")
			called = open((dst or self.Lib.resultDir.writable) / self.Lib.queryName+"_snps.tsv", "w")
			noCoverage = unique = notCalled = open((dst or self.Lib.resultDir.writable) / self.Lib.queryName+"_not_called.tsv", "w")
			
			called.write(header)
			noCoverage.write(header)

		for genomeID, genome, genbankID, refseqID, assemblyName in self.database.references:
			'''Print SNPs to tab separated file'''
			for nodeID, position, ancestral, derived, chromosome in self.database[NodeID, Position, AncestralBase, DerivedBase, Chromosome, GenomeID == genomeID]:
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
	