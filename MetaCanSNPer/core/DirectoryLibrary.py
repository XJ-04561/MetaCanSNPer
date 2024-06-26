

import logging.handlers
from MetaCanSNPer.Globals import *
import MetaCanSNPer.Globals as Globals

from MetaCanSNPer.core.Hooks import Hooks
from MetaCanSNPer.modules.Database import MetaCanSNPerDatabase

'''Container and handler of directories and files'''
class DirectoryLibrary(SoftwareLibrary, Logged):
	'''
		This structure allows for defining file managing as well and making sure that directory and file path
		information is simply passed around by handing around the same 'Library'.
	'''

	SOFTWARE_NAME = SOFTWARE_NAME
	
	organism : str

	database : MetaCanSNPerDatabase
	settings : dict = cached_property(lambda self : dict())
	hooks : Hooks = cached_property(lambda self : Hooks())

	query : PathList[FilePath]
	"""The data file being queried. Can be multiple files in the case of Illumina and other datasets split into parts."""
	queryName : str = cached_property(lambda self : self.query.name)
	"""Defaults to a sequence alignment of the query files involved."""
	sessionName : str = cached_property(lambda self : f"Sample-{self.queryName}-{time.strftime('%Y-%m-%d_%H-%M-%S', time.localtime())}")
	"""Defaults to 'Sample-[QUERY_NAME]-[CURRENT_DATE]'"""
	
	targetDir : PathGroup
	refDir : PathGroup
	SNPDir : PathGroup
	databaseDir : PathGroup
	tmpDir : PathGroup
	outDir : PathGroup
	resultDir : DirectoryPath
	logDir : DirectoryPath

	targetSNPs : dict[str,Path]
	references : dict[str,Path]
	
	@property
	def query(self):
		try:
			return self.__dict__["query"]
		except:
			raise AttributeError("Attribute not yet set 'query'")
	@query.setter
	def query(self, value : str|list):
		if isinstance(value, str):
			value = [value]

		queryList = [FilePath(self.targetDir.find(q)) for q in value]
		
		if None in queryList:
			raise FileNotFoundError(f"Query files not found: {', '.join(filter(lambda x:x not in self.targetDir, value))}")
		
		self.__dict__["query"] = FileList(queryList)

	@Default["workDir", "userDir"]
	def targetDir(self) -> DirectoryGroup:
		return self.workDir | self.userDir
	
	@Default["dataDir", "organism"]
	def refDir(self) -> DirectoryGroup:
		return self.dataDir / "References" / self.organism
	
	@Default["dataDir", "organism"]
	def SNPDir(self) -> DirectoryGroup:
		return self.dataDir / "SNPs" / self.organism
	
	@Default["dataDir"]
	def databaseDir(self) -> DirectoryGroup:
		return self.dataDir / "Databases"
	
	@Default["userCacheDir"]
	def tmpDir(self) -> DirectoryPath:
		if self.settings.get("saveTemp"):
			return self.userCacheDir.writable / f"{self.organism}_{self.queryName}"
		else:
			return PseudoPathyFunctions.createTempDir(f"{self.organism}_{self.queryName}", dir=self.userCacheDir.writable)
	
	@Default["targetDir", "userDir", "SOFTWARE_NAME"]
	def outDir(self) -> DirectoryGroup:
		return DirectoryPath(self.workDir) | (self.userDir / self.SOFTWARE_NAME)
	
	@Default["outDir", "sessionName"]
	def resultDir(self) -> DirectoryPath:
		"""Should not be overriden, instead look to instance.outputDir and instance.sessionName separately.
		This will automatically change to reflect those two values. ([OUTPUT_DIR]/[SESSION_NAME]/)"""
		return self.outDir / self.sessionName

	@Default["outDir", "sessionName"]
	def logDir(self) -> DirectoryPath:
		return UniqueDirectoryPath(self.userCacheDir.writable / self.sessionName)
	
	@logDir.setter
	def logDir(self, value):
		if isinstance(value, (PathGroup, Path)):
			self.__dict__["logDir"] = value
		else:
			self.__dict__["logDir"] = Path(value)
	
	@Default["refDir", "database"]
	def targetSNPs(self) -> dict[str,Path]:
		"""{GENOME_NAME : TARGET_SNPS_FILE_PATH}"""
		return {genome:self.SNPDir.find(f"{genome}.vcf") for _, genome, *_ in self.database.references}
	
	@Default["refDir", "database"]
	def references(self) -> dict[str,Path]:
		"""{GENOME_NAME : REFERENCE_GENOME_FILE_PATH}"""
		return {genome:self.refDir.find(f"{assemblyName}.fna") or self.refDir.find(f"{assemblyName}.fasta") for _, genome, strain, genbankID, refseqID, assemblyName in self.database.references}
	
	maps : dict[str,Path]			= cached_property(lambda self : dict())
	"""{GENOME_NAME : MAPPED_QUERY_FILE_PATH}"""
	alignments : dict[str,Path]		= cached_property(lambda self : dict())
	"""{GENOME_NAME : ALIGNED_QUERY_FILE_PATH}"""
	resultSNPs : dict[str,Path]		= cached_property(lambda self : dict())
	"""{GENOME_NAME : CALLED_SNPS_FILE_PATH}"""
	
	def __init__(self, organism : str, query : list[str]|str, settings : dict={}, **kwargs):
		"""Directories that can be passed as kwargs:
		workDir, targetDir, tmpDir, refDir, databaseDir, outDir

		Needs only a settings dictionary implementing all keys specified in the 'defaultFlags.toml' that should be in
		the installation directory, but keyword arguments can be used to force use directories.
		"""

		self.LOG.info("Creating DirectoryLibrary object.")
		
		self.query = query
		self.organism = organism
		self.settings |= settings

		super().__init__(**{name:value for name, value in kwargs.items() if value is not None})
		self.LOG = self.LOG.getChild(f"[{self.sessionName}]")
		
		self.targetDir.create(purpose="w")
		self.refDir.create(purpose="w")
		self.SNPDir.create(purpose="w")
		self.databaseDir.create(purpose="w")
		self.tmpDir.create(purpose="w")
		self.outDir.create(purpose="w")
		self.logDir.create(purpose="w")

	if Globals.MAX_DEBUG:
		def __getattribute__(self, name):
			super().__getattribute__("LOG").debug(f"[0x{id(self):0>16X}] Getting {name!r}")
			ret = super().__getattribute__(name)
			try:
				if hasattr(ret, "__str__"):
					super().__getattribute__("LOG").debug(f"[0x{id(self):0>16X}] Got {name!r} as {ret} ({ret!r})")
				else:
					super().__getattribute__("LOG").debug(f"[0x{id(self):0>16X}] Got {name!r} as {ret!r}")
			except:
				super().__getattribute__("LOG").debug(f"[0x{id(self):0>16X}] Got {name!r} [Could not format]")
			return ret

		def __setattr__(self, name, value):
			try:
				if hasattr(value, "__str__"):
					super().__getattribute__("LOG").debug(f"[0x{id(self):0>16X}] Setting {name!r} to: {value}({value!r})")
				else:
					super().__getattribute__("LOG").debug(f"[0x{id(self):0>16X}] Setting {name!r} to: {value!r}")
			except:
				super().__getattribute__("LOG").debug(f"[0x{id(self):0>16X}] Setting {name!r} to: [Could not format]")
			super().__setattr__(name, value)