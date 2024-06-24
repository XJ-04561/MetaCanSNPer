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
import MetaCanSNPer.core.Aligners as Aligners
import MetaCanSNPer.core.Mappers as Mappers
import MetaCanSNPer.core.SNPCallers as SNPCallers
from MetaCanSNPer.core.TerminalUpdater import TerminalUpdater

from MetaCanSNPer.modules.Database import MetaCanSNPerDatabase, TreeTable, NodeID, Position, GenomeID, AncestralBase, DerivedBase, Chromosome, Genotype, ChromosomeID, Genome, Parent
from MetaCanSNPer.modules.Downloader import DatabaseDownloader, DownloaderReportHook, ReferenceDownloader


class MetaCanSNPer(Logged):
	
	hooks : Hooks = GlobalHooks

	outputTemplate = "{refName}_{queryName}.{outFormat}"
	organism : str = property(lambda self:self.Lib.organism,
							  lambda self, value: setattr(self.Lib, "organism", value),
							  lambda self: delattr(self.Lib, "organism"))

	settings : dict
	Lib : DirectoryLibrary
	database : MetaCanSNPerDatabase

	databasePath : Path
	databaseName : str = Default["organism"](lambda self:f"{self.organism}.db")
	
	referenceFiles = property(lambda self:self.Lib.references)
	references = property(lambda self:self.database.references,
						  doc="""Fetch names of reference genomes in connected database.""")

	query : FileList = property(lambda self:self.Lib.query,
								lambda self, value:setattr(self.Lib, "query", value),
								lambda self:delattr(self.Lib, "query"), DirectoryLibrary.query.__doc__)
	queryName : str = Default["Lib.query"](lambda self:self.Lib.queryName)
	sessionName : str = property(lambda self:self.Lib.sessionName,
								 lambda self, value: setattr(self.Lib, "sessionName", value),
								 lambda self: delattr(self.Lib, "sessionName"))
	
	SNPresults : dict = Default["Lib.query"](lambda self:defaultdict(dict))
	exceptions : list[Exception]

	@overload
	def __init__(self, /, organism : str, query : Iterable[str]|str, *, lib : DirectoryLibrary=None, database : str=None,
			  	hooks : Hooks=None, settings : dict={}, settingsFile : str=None, sessionName : str=None): ...

	def __init__(self, /, organism : str, query : Iterable[str]|str, **kwargs):

		self.LOG.info(f"Initializing MetaCanSNPer object at 0x{id(self):0>16X}.")
		self.exceptions = []
		self.hooks = kwargs.get("hooks", self.hooks)
		self.hooks.addHook("ReportError", target=lambda eventInfo : self.exceptions.append(eventInfo["exception"]))
		self.settings = DEFAULT_SETTINGS.copy()
		
		if kwargs.get("settingsFile"):
			if not os.path.exists(kwargs["settingsFile"]):
				raise FileNotFoundError(f"Could not find the settingsFile {kwargs['settingsFile']!r}")
			self.settings |= loadFlattenedTOML(kwargs["settingsFile"])
			
		for flag, value in kwargs.pop("settings", {}).items():
			if isinstance(value, bool) or value:
				self.settings[flag] = value

		self.Lib = DirectoryLibrary(organism, query, settings=self.settings, hooks=self.hooks, **kwargs)

		self.LOG = self.LOG.getChild(f"[{self.Lib.queryName}]")

		self.organism = organism

		requiredDeps = []

		# if not self.settings.get("mapper") and not self.settings.get("aligner"):
		# 	raise ValueError("No aligner or mapper has been provided through flags or settings-files.")
		# elif not self.settings.get("snpCaller"):
		# 	raise ValueError("No SNP caller has been provided through flags or settings-files.")
		
		if not self.settings.get("dryRun"):
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
		
	def setDatabase(self, databaseName : str=None, sequential : bool=False):
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
			DD = DatabaseDownloader(outDir, hooks=self.hooks, sequential=sequential)

			DD.download(databaseName, databaseName)
			
			DD.wait()
		
		self.Lib.database = self.database = MetaCanSNPerDatabase(self.databasePath, "r", organism=self.organism)
		self.LOG.info(f"Database {self.databaseName} loaded from: {self.databasePath}!")
		if Globals.MAX_DEBUG:
			nt = "\n\t"
			self.LOG.debug(f"Database {self.databaseName} has attributes:\n{nt.join([name+' = '+str(getattr(self.database, name)) for name in dir(self.database) if not name.startswith('_')])}")
	
	def setReferenceFiles(self, references : Iterable[tuple[int,str,str,str,str,str]]=None, sequential : bool=False):

		if not references:
			assert hasattr(self, "database"), "Database not yet set. Setting references requires either a table of references as an argument or that a database is connected."
			references = self.database.references

		directory = self.Lib.refDir.writable
		DD = ReferenceDownloader(directory, hooks=self.hooks, sequential=sequential)

		del self.Lib.references
		
		for genomeID, genome, strain, genbankID, refseqID, assemblyName in references:
			filename = f"{assemblyName}.fna"
			if self.Lib.references[genome]:
				self.LOG.debug(f"Reference genome already exists! {self.Lib.references[genome]}")
				if not sequential:
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
			if Globals.DRY_RUN:
				continue
			vcfFile = openVCF(filePath, "r")
			for chrom, pos, counts in vcfFile["CHROM", "POS", "AD"]:
				depth = sum(counts)
				LOG.debug(f"(CHROM, POS, Allele-readDepth (A, T, C, G) ) {chrom}, {pos}, {depth} = {' + '.join(map(str, counts))}")
				
				chromID = next(self.database[ChromosomeID, Chromosome==chrom, Genome==genome], None)
				nodeID, derivedBase, ancestralBase = self.database[NodeID, DerivedBase, AncestralBase, Position==pos, ChromosomeID==chromID]
				
				match derivedBase:
					case "A":
						derivedReads = counts[0]
					case "T":
						derivedReads = counts[1]
					case "C":
						derivedReads = counts[2]
					case "G":
						derivedReads = counts[3]
					case _:
						derivedReads = 0
				match ancestralBase:
					case "A":
						ancestralReads = counts[0]
					case "T":
						ancestralReads = counts[1]
					case "C":
						ancestralReads = counts[2]
					case "G":
						ancestralReads = counts[3]
					case _:
						ancestralReads = 0
				self.SNPresults[nodeID][pos] = (derivedReads, ancestralReads, counts, depth)
		
		LOG.info("Got nodes: " + ", ".join(map(str, self.SNPresults)))
	
	def traverseTree(self) -> tuple[int,dict[int,list[int,int,list[int,int,int,int],int]]]:
		'''Depth-first tree search.'''
		self.LOG.info(f"Traversing tree to get Genotype called.")
		
		rootNode = self.database.tree
		assert next(rootNode.children, False)

		nodeScores = {rootNode.node : [0,0,[0,0,0,0],0]}
		
		paths = [[rootNode]]
		while paths[-1] != []:
			paths.append([])
			for parent in paths[-2]:
				for child in parent.children:
					if Globals.DRY_RUN:
						continue
					
					nodeScores[child.node] = [0,0,[0,0,0,0],0]
					for nodeID, pos, *_ in self.database.SNPsByNode[child.node]:
						derived, ancestral, counts, depth = self.SNPresults[nodeID][pos]
						nodeScores[child.node][0] += derived
						nodeScores[child.node][1] += ancestral
						nodeScores[child.node][2] = tuple(n+o for n, o in zip(counts, nodeScores[child.node][2]))
						nodeScores[child.node][3] += depth
					paths[-1].append(child)
		
		paths = paths[:-1]
		
		averageDepth = sum(x[3] for x in nodeScores.values()) / (len(nodeScores)-1)
		for nodeID, (derived, ancestral, counts, depth) in reversed(list(nodeScores.items())):
			if derived == max(counts) and derived > averageDepth * 0.25:
				self.LOG.info(f"Finished traversing tree.")
				return nodeID, nodeScores
		else:
			self.LOG.info(f"Finished traversing tree.")
			return rootNode.node, nodeScores

	def saveResults(self, dst : str=None):
		
		outDir = dst or self.Lib.resultDir.writable
		self.LOG.debug(f"open({outDir!r} / {self.Lib.queryName!r}+'_final.tsv', 'w')")
		with open((outDir) / self.Lib.queryName+"_final.tsv", "w") as finalFile:
			finalNodeID, scores = self.traverseTree()
			
			genotype = self.database[Genotype, TreeTable, NodeID == finalNodeID]
			derived, ancestral, counts, depth = scores[finalNodeID]
			
			if derived > 1:
				logRatio = format(math.log10(derived)/math.log10(depth), ".3f")
			else:
				logRatio = "NaN"
			finalFile.write(f"{'Genotype':<20} {'log_10(Called)/log_10(Depth)':>30} {'Called':>20} {'Ancestral':>20} {'Non-Canonical Bases':>20} {'Depth':>20}\n")
			finalFile.write(f"{genotype:<20} {logRatio:>30} {derived:>20} {ancestral:>20} {depth-(derived+ancestral):>20} {depth:>20}\n\n")
			for nodeID, (derived, ancestral, counts, depth) in scores.items():
				genotype = self.database[Genotype, TreeTable, NodeID == nodeID]
				
				if derived > 1:
					logRatio = format(math.log10(derived)/math.log10(depth), ".3f")
				else:
					logRatio = "NaN"
				finalFile.write(f"{genotype:<20} {logRatio:>30} {derived:>20} {ancestral:>20} {depth-(derived+ancestral):>20} {depth:>20}\n")
		
		# Export graph to GraphML format.

		HEADER = ("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
		"<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"  \n"
		"    xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
		"    xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n"
		"  <key id=\"genotype\" for=\"node\" attr.name=\"Variant\" attr.type=\"string\"/>\n"
		"  <key id=\"called\" for=\"edge\" attr.name=\"Called SNP\" attr.type=\"int\">\n"
		"    <default>NaN</default>\n"
		"  </key>\n"
		"  <key id=\"ancestral\" for=\"edge\" attr.name=\"Ancestral Bases\" attr.type=\"int\">\n"
		"    <default>NaN</default>\n"
		"  </key>\n"
		"  <key id=\"nonCanon\" for=\"edge\" attr.name=\"Non-Canonical Bases\" attr.type=\"int\">\n"
		"    <default>NaN</default>\n"
		"  </key>\n"
		"  <key id=\"depth\" for=\"edge\" attr.name=\"Read Depth\" attr.type=\"int\">\n"
		"    <default>NaN</default>\n"
		"  </key>\n"
		"  <key id=\"ratio\" for=\"edge\" attr.name=\"Ratio Called/Depth\" attr.type=\"float\">\n"
		"    <default>NaN</default>\n"
		"  </key>\n"
		"  <key id=\"logRatio\" for=\"edge\" attr.name=\"Ratio log_10(Called)/log_10(Depth)\" attr.type=\"float\">\n"
		"    <default>NaN</default>\n"
		"  </key>\n"
		f"  <graph id=\"{self.queryName}\" edgedefault=\"directed\">\n")
		nodeEntries = []
		edges = []
		for nodeID, (derived, ancestral, counts, depth) in scores.items():
			parent, genotype = self.database[Parent, Genotype, NodeID == nodeID]
			assert isinstance(genotype, str), f"{nodeID=} should have a genotype string, but instead had {genotype=}"
			nodeEntries.append(
				f"    <node id=\"{nodeID}\">\n"
				f"      <data id=\"genotype\">{genotype}</data>\n"
				"    </node>\n"
			)
			
			edges.append(
				f"    <edge id=\"{genotype}\" source=\"{parent}\" target=\"{nodeID}\">\n"
				f"      <data id=\"called\">{derived}</data>\n"
				f"      <data id=\"ancestral\">{ancestral}</data>\n"
				f"      <data id=\"nonCanon\">{depth-(derived+ancestral)}</data>\n"
				f"      <data id=\"depth\">{depth}</data>\n"
				f"      <data id=\"ratio\">{derived/depth if depth > 0 else 'NaN'}</data>\n"
				f"      <data id=\"logRatio\">{math.log10(derived)/math.log10(depth) if derived > 1 else 'NaN'}</data>\n"
				f"    </edge>\n"
			)
		FOOTER = ("  </graph>\n"
		"</graphml>\n")
		with open((outDir) / self.Lib.queryName+".graphml.xml", "w") as graphFile:
			graphFile.write(HEADER)
			for nodeEntry in nodeEntries:
				graphFile.write(nodeEntry)
			for edge in edges:
				graphFile.write(edge)
			graphFile.write(FOOTER)
		
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
				*_, counts, depth = self.SNPresults[nodeID][position]
				genotype = self.database[Genotype, NodeID==nodeID]
				if Globals.MAX_DEBUG: self.LOG.debug(f"Got {genotype=} and baseCounts={counts} for {nodeID=}")
				
				N = 'ATCG'[max(range(4), key=counts.__getitem__)]
				entry = f"{genotype}\t{genome}\t{chromosome}\t{position}\t{ancestral}\t{derived}\t{N}\t{', '.join(c+'='+str(n) for c,n in zip('ATCG', counts))}\n"
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