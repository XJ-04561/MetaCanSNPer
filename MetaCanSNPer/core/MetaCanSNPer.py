'''
MetaCanSNPer module: A toolkit for SNP-typing using NGS data.
Copyright (C) 2024 Fredrik Sörensen @ Umeå University
'''

from VariantCallFixer.Functions import getSNPdata

## import MetaCanSNPer specific modules
from MetaCanSNPer.Globals import *
import MetaCanSNPer.Globals as Globals
from MetaCanSNPer.core.DirectoryLibrary import DirectoryLibrary
from MetaCanSNPer.core.Hooks import Hooks, DummyHooks, GlobalHooks
from MetaCanSNPer.core.Wrappers import Aligner, Mapper, SNPCaller, ProcessWrapper
import MetaCanSNPer.software.Aligners as Aligners
import MetaCanSNPer.software.Mappers as Mappers
import MetaCanSNPer.software.SNPCallers as SNPCallers
from MetaCanSNPer.core.TerminalUpdater import TerminalUpdater

from MetaCanSNPer.core.Database import MetaCanSNPerDatabase, TreeTable, NodeID, Position, GenomeID, AncestralBase, DerivedBase, Chromosome, Genotype, ChromosomeID, Genome
from MetaCanSNPer.core.Downloaders import DatabaseDownloader, DownloaderReportHook, ReferenceDownloader


class MetaCanSNPer(Logged):
	
	hooks = GlobalHooks

	outputTemplate = "{refName}_{queryName}.{outFormat}"
	organism : str

	settings : dict
	Lib : DirectoryLibrary
	database : MetaCanSNPerDatabase

	databasePath : Path
	databaseName : str = Default["organism"](lambda self:f"{self.organism}.db")
	
	referenceFiles = property(lambda self:self.Lib.references)
	references = property(lambda self:self.database.references, doc="""Fetch names of reference genomes in connected database.""")

	query = property(lambda self:self.Lib.query, lambda self, value:setattr(self.Lib, "query", value), lambda self:delattr(self.Lib, "query"), DirectoryLibrary.query.__doc__)
	queryName : str = Default["Lib.query"](lambda self:self.Lib.queryName)
	sessionName : str = Default["Lib.query"](lambda self:self.Lib.sessionName)
	
	SNPresults : dict = Default["Lib.query"](lambda self:dict())
	exceptions : list[Exception]

	@overload
	def __init__(self, /, organism : str, query : list[str]|str, *, lib : DirectoryLibrary=None, database : str=None,
			  	hooks : Hooks=None, settings : dict={}, settingsFile : str=None, sessionName : str=None): ...

	def __init__(self, /, organism : str, query : list[str]|str, **kwargs):

		self.LOG.info(f"Initializing MetaCanSNPer object at 0x{id(self):0>16X}.")
		self.exceptions = []
		self.hooks = kwargs.get("hooks", self.hooks)
		self.hooks.addHook("ReportError", target=lambda eventInfo : self.exceptions.append(eventInfo["exception"]))
		self.settings = DEFAULT_SETTINGS.copy()
		
		if kwargs.get("settingsFile"):
			if not os.path.exists(kwargs["settingsFile"]):
				raise FileNotFoundError(f"Could not find the settingsFile {kwargs['settingsFile']!r}")
			self.settings |= loadTOML(kwargs["settingsFile"])
			
		for flag, value in kwargs.get("settings", {}).items():
			if isinstance(value, bool) or value:
				self.settings[flag] = value

		self.Lib = DirectoryLibrary(organism, query, settings=self.settings, hooks=self.hooks)

		self.LOG = self.LOG.getChild(f"[{self.Lib.queryName}]")

		self.organism = organism

		requiredDeps = []

		if not self.settings.get("mapper") and not self.settings.get("aligner"):
			raise ValueError("No aligner or mapper has been provided through flags or settings-files.")
		elif not self.settings.get("snpCaller"):
			raise ValueError("No SNP caller has been provided through flags or settings-files.")
		
		if not self.settings.get("dry_run"):
			if self.settings.get("mapper"):
				requiredDeps.extend(Mapper.get(self.settings["mapper"]).dependencies)
			if self.settings.get("aligner"):
				requiredDeps.extend(Aligner.get(self.settings["aligner"]).dependencies)
			if self.settings.get("snpCaller"):
				requiredDeps.extend(SNPCaller.get(self.settings["snpCaller"]).dependencies)

			missed = []
			for dep in requiredDeps:
				if not shutil.which(dep) and not shutil.which(dep+".exe"):
					missed.append(dep)
			if len(missed) == 1:
				raise MissingDependency(f"Missing required dependency: {missed[0]}.")
			elif missed:
				nt = "\n\t"
				raise MissingDependency(f"Missing required dependencies:\n{nt.join(missed)}.")

	def setOrganism(self, organism : str):
		self.Lib.organism = self.organism = organism
		
	def setDatabase(self, databaseName : str=None):
		"""databaseName Defaults to `{organism}.db`"""

		if databaseName is not None:
			self.databaseName = databaseName
		else:
			databaseName = self.databaseName
		
		self.LOG.debug(f"Setting database to:{databaseName}")
		
		for directory in self.Lib.databaseDir:
			if os.path.exists(directory / databaseName):
				self.LOG.debug(f"{databaseName} in {directory}")
				if MetaCanSNPerDatabase(directory / databaseName, "r", organism=self.organism).valid:
					self.LOG.debug(f"{directory / databaseName} is valid")
					self.hooks.trigger("DatabaseDownloaderSkipped", {"name" : databaseName, "value" : 2})
					self.databasePath = directory / databaseName
					break
			else:
				self.LOG.debug(f"{databaseName} not in {directory}")
		else:
			self.LOG.info(f"No valid database={databaseName} found. Looking for writeable versions or directories that can be updated or downloaded to, respectively.")
			
			outDir = self.Lib.databaseDir.writable
			self.databasePath = outDir / databaseName
			DD = DatabaseDownloader(outDir, hooks=self.hooks)

			DD.download(databaseName, databaseName)
			
			DD.wait()
		
		self.Lib.database = self.database = MetaCanSNPerDatabase(self.databasePath, "r", organism=self.organism)
		self.LOG.info(f"Database {self.databaseName} loaded from: {self.databasePath}!")
		if Globals.MAX_DEBUG:
			nt = "\n\t"
			self.LOG.debug(f"Database {self.databaseName} has attributes:\n{nt.join([name+' = '+str(getattr(self.database, name)) for name in dir(self.database) if not name.startswith('_')])}")
	
	def setReferenceFiles(self, references : Iterable[tuple[int,str,str,str,str,str]]=None):

		if not references:
			assert hasattr(self, "database"), "Database not yet set. Setting references requires either a table of references as an argument or that a database is connected."
			references = self.database.references

		directory = self.Lib.refDir.writable
		DD = ReferenceDownloader(directory, hooks=self.hooks)

		del self.Lib.references
		
		for genomeID, genome, strain, genbankID, refseqID, assemblyName in references:
			filename = f"{assemblyName}.fna"
			if self.Lib.references[genome]:
				self.LOG.debug(f"Reference genome already exists! {self.Lib.references[genome]}")
				self.hooks.trigger("ReferenceDownloaderSkipped", {"name" : filename, "value" : 2})
				continue
			
			DD.download((genbankID, assemblyName), filename)
			
			self.Lib.references[genome] = directory / filename
		
		DD.wait()

		if not self.database.valid:
			self.database.reopen("w")
			self.database.fix()
			if not self.database.valid:
				raise self.database.exception
			self.database.reopen("r")

	'''MetaCanSNPer set directories'''
	
	def setTargetDir(self, path : str):		setattr(self.Lib, "targetDir", path)
	def setRefDir(self, path : str):		setattr(self.Lib, "refDir", path)
	def setDatabaseDir(self, path : str):	setattr(self.Lib, "databaseDir", path)
	def setTmpDir(self, path : str):		setattr(self.Lib, "tmpDir", path)
	def setOutDir(self, path : str):		setattr(self.Lib, "outDir", path)
	"""
	`path` can be relative or absolute. Check the docstring for the DirectoryLibrary to see where the relative
	path will start from.
	"""

	'''MetaCanSNPer set values'''

	def setSessionName(self, name):
		self.Lib.sessionName = self.sessionName = name

	'''Indexing methods'''

	def runSoftware(self, softwareClass : ProcessWrapper, outputDict : dict={}, flags : list=[]):
		'''Align sequences using subprocesses.'''

		self.LOG.info(f"Creating software wrapper for {softwareClass.softwareName!r} of type {softwareClass.__name__!r}")
		software : ProcessWrapper = softwareClass(self.Lib, self.database, self.outputTemplate, out=outputDict, hooks=self.hooks, flags=flags, settings=self.settings)

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
		MapperType : Mapper = Mapper.get(softwareName)

		if self.Lib.references is None:
			self.LOG.error("References not set.")
			raise FileNotFoundError("References not set. Can be set with MetaCanSNPer.setReferences")
		
		self.LOG.info(f"Loaded a total of {sum(1 for r in self.database.references)} References.")
		
		self.runSoftware(MapperType, outputDict=self.Lib.maps, flags=flags)

		self.LOG.info(f"Result of Mapping in: {self.Lib.maps}")

	def createAlignment(self, softwareName : str, flags : list=[]):
		''''''
		self.LOG.info(f"Creating alignment using:{softwareName}")
		AlignerType : Aligner = Aligner.get(softwareName)

		if self.Lib.references is None:
			self.LOG.error("References not set.")
			raise FileNotFoundError("References not set. Can be set with MetaCanSNPer.setReferences")
		
		self.LOG.info(f"Loaded a total of {sum(1 for r in self.database.references)} References.")
		
		self.runSoftware(AlignerType, outputDict=self.Lib.alignments, flags=flags)

		self.LOG.info(f"Result of Alignment in: {self.Lib.alignments}")

	def callSNPs(self, softwareName : str, flags : list=[]):
		''''''
		LOG = self.LOG.getChild("callSNPs")
		LOG.info(f"Calling SNPs using:{softwareName}")
		SNPCallerType : SNPCaller = SNPCaller.get(softwareName)

		if self.Lib.references is None:
			LOG.error("References not set.")
			raise FileNotFoundError("References not set. Can be set with MetaCanSNPer.setReferences")
		
		LOG.info("Loading SNPs from database.")
		from MetaCanSNPer.modules.CanSNP2VCF import CanSNP2VCF
		CanSNP2VCF(self.Lib)
		LOG.info(f"Loaded a total of {len(self.database.SNPs)} SNPs.")
		
		self.runSoftware(SNPCallerType, outputDict=self.Lib.resultSNPs, flags=flags)

		LOG.info(f"Result of SNPCalling in: {self.Lib.resultSNPs}")

		for genome, filePath in self.Lib.resultSNPs.items():
			for nodeID in self.database[NodeID, Genome==genome]:
				if nodeID not in self.SNPresults:
					self.SNPresults[nodeID] = {}
			if Globals.DRY_RUN:
				continue
			for (chrom, pos), ref in getSNPdata(filePath, key=["CHROM", "POS"], values="REF"):
				LOG.debug(f"({chrom}, {pos}), {ref}")
				chromID = list(self.database[ChromosomeID, Chromosome==chrom, Genome==genome])[0]
				nodeID = self.database[NodeID, Position==pos, ChromosomeID==chromID]
				self.SNPresults[nodeID][pos] = ref
		LOG.info("Got nodes: " + ", ".join(map(str, self.SNPresults)))
	
	def traverseTree(self):
		'''Depth-first tree search.'''
		self.LOG.info(f"Traversing tree to get Genotype called.")
		rootNode = self.database.tree
		paths = []
		nodeScores = {rootNode.node : 0}

		award = (1, -1, 0)

		miscCalls = []
		paths.append([rootNode])
		while paths[-1] != []:
			paths.append([])
			for parent in paths[-2]:
				for child in parent.children:
					if Globals.DRY_RUN:
						continue
					
					nodeScores[child.node] = nodeScores[parent.node]
					for nodeID, pos, anc, der, *_ in self.database.SNPsByNode[child.node]:
						base = self.SNPresults[nodeID].get(pos, "-")
						
						if der == base:
							nodeScores[child.node] += award[0]
						elif anc == base:
							nodeScores[child.node] += award[1]
						elif base == "-":
							pass
						else:
							miscCalls.append((nodeID, pos, anc, der, base))
					paths[-1].append(child)
		if miscCalls:
			self.LOG.warning("These nodes had non-recognized bases:\n\t NODE_ID, POS, ANCESTRAL, DERIVED, TARGET\n\t" + "\n\t".join(map(str, miscCalls)))
		paths = paths[:-1]
		self.LOG.info(f"Finished traversing tree.")
		return max(nodeScores.items(), key=lambda x: x[1]), nodeScores

	'''Functions'''

	def saveResults(self, dst : str=None):
		
		outDir = dst or self.Lib.resultDir.writable
		self.LOG.debug(f"open({outDir!r} / {self.Lib.queryName!r}+'_final.tsv', 'w')")
		with open((outDir) / self.Lib.queryName+"_final.tsv", "w") as finalFile:
			(finalNodeID, score), scores = self.traverseTree()
			res = self.database[Genotype, TreeTable, NodeID == finalNodeID]
			
			finalFile.write("{:<20}{score}\n".format(res, score=score))
			if self.settings.get("debug"):
				finalFile.write("\n")
				for nodeID in scores:
					finalFile.write("{:<20}{score}\n".format(self.database[Genotype, TreeTable, NodeID == nodeID], score=scores[nodeID]))
		return outDir

	def saveSNPdata(self, dst : str=None):
		""""""

		outDir = dst or self.Lib.resultDir.writable
		header = "Name\tReference\tChromosome\tPosition\tAncestral base\tDerived base\tTarget base\n"

		self.LOG.debug(f"open({outDir!r} / {self.Lib.queryName!r}+'_snps.tsv', 'w')")
		self.LOG.debug(f"open({outDir!r} / {self.Lib.queryName!r}+'_not_called.tsv', 'w')")
		self.LOG.debug(f"open({outDir!r} / {self.Lib.queryName!r}+'_no_coverage.tsv', 'w')")
		called = open((outDir) / self.Lib.queryName+"_snps.tsv", "w")
		notCalled = open((outDir) / self.Lib.queryName+"_not_called.tsv", "w")
		noCoverage = open((outDir) / self.Lib.queryName+"_no_coverage.tsv", "w")
		
		called.write(header)
		notCalled.write(header)
		noCoverage.write(header)

		for genomeID, genome, strain, genbankID, refseqID, assemblyName in self.database.references:
			'''Print SNPs to tab separated file'''
			for nodeID, position, ancestral, derived, chromosome in self.database[NodeID, Position, AncestralBase, DerivedBase, Chromosome, GenomeID == genomeID]:
				if Globals.DRY_RUN:
					continue
				N = self.SNPresults[nodeID].get(position, "-")
				genotype = self.database[Genotype, NodeID==nodeID]
				if Globals.MAX_DEBUG: self.LOG.debug(f"Got {genotype=} and base={N!r} for {nodeID=}")

				entry = f"{genotype}\t{genome}\t{chromosome}\t{position}\t{ancestral}\t{derived}\t{N}\n"
				if N == derived or N == ancestral:
					called.write(entry)
				elif N.isalpha():
					notCalled.write(entry)
				else:
					noCoverage.write(entry)
		
		called.close()
		noCoverage.close()
		notCalled.close()
			
		return outDir