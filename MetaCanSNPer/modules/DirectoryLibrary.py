

import os, random
from PseudoPathy import MinimalPathLibrary, PathLibrary, PathGroup, Path, DirectoryPath, FilePath, DisposablePath
from PseudoPathy.Functions import createTemp
import PseudoPathy.Globals

import MetaCanSNPer.modules.LogKeeper as LogKeeper
from MetaCanSNPer.modules.DownloadReferences import DownloadQueue
from MetaCanSNPer.modules.FileNameAlignment import align as fileNameAlign
from MetaCanSNPer.Globals import *

LOGGER = LogKeeper.createLogger(__name__)
PseudoPathy.Globals.LOGGER = LOGGER

class FileList(list):
	def __new__(cls, *data):
		obj = super(FileList, cls).__new__(cls, *data)
		return obj
	
	def __str__(self):
		return " ".join(self)

'''Container and handler of directories and files'''
class DirectoryLibrary(PathLibrary):
	'''
		This structure allows for defining file managing as well and making sure that directory and file path
		information is simply passed around by handing around the same 'Library'.
	'''

	targetDir : PathGroup
	refDir : PathGroup
	databaseDir : PathGroup
	tmpDir : PathGroup
	outDir : PathGroup
	resultDir : DirectoryPath
	logDir : DirectoryPath

	maps : MinimalPathLibrary
	alignments : MinimalPathLibrary
	targetSNPS : MinimalPathLibrary
	resultSNPs : MinimalPathLibrary
	references : MinimalPathLibrary
	"""Current `MinimalPathLibrary` of reference files:
		{GENOME_NAME : REFERENCE_PATH}"""

	# Non-Pathy attributes. Must be initialized using object.__setattr__

	query : FileList[FilePath]
	sessionName : str
	queryName : str
	chromosomes : list[str]
	
	def __init__(self, settings : dict, reference : str=None, sessionName="Unnamed-Session", **kwargs):
		'''Directories that can be passed as kwargs:
		workDir, targetDir, tmpDir, refDir, databaseDir, outDir

		Needs only a settings dictionary implementing all keys specified in the 'defaultFlags.toml' that should be in
		the installation directory, but keyword arguments can be used to force use directories.
		'''

		LOGGER.debug("Creating: {}".format(self))

		if (workDir := kwargs.pop("workDir", None)) is not None: os.chdir(workDir)

		super(DirectoryLibrary, self).__init__(self, **{name:kwargs[name] for name in {"workDir", "userDir", "installDir"}.intersection(kwargs)})

		self.settings = settings
		self.sessionName = sessionName

		# Use kwargs as a first, and settings as a second.

		LOGGER.debug(f"Work dir:    {self.workDir!r}")
		LOGGER.debug(f"Install dir: {self.installDir!r}")
		LOGGER.debug(f"User dir:    {self.userDir!r}")

		super(DirectoryLibrary, self).__init__()
		
		PseudoPathy.Globals.DISPOSE = self.settings.get("saveTemp", True)
		
		self.setTargetDir(kwargs.pop("targetDir", None) or self.settings.get("targetDir"))
		if reference is not None: self.setRefDir(reference, refDir=kwargs.pop("refDir", None) or self.settings.get("refDir"))
		self.setDatabaseDir(kwargs.pop("databaseDir", None) or self.settings.get("databaseDir"))
		self.setTmpDir(kwargs.pop("tmpDir", None) or self.settings.get("tmpDir"))
		self.setOutDir(kwargs.pop("outDir", None) or self.settings.get("outDir"))
		# If argument is None, then a default set of directories are used. Non-abs paths are appended to a set of Paths (workDir and userDir mostly)
		
		LOGGER.debug(str(self))
		
		object.__setattr__(self, "query", [])
		object.__setattr__(self, "queryName", [])
		object.__setattr__(self, "sessionName", "")
		object.__setattr__(self, "chromosomes", [])

		self.maps = MinimalPathLibrary()
		self.alignments = MinimalPathLibrary()
		self.targetSNPS = MinimalPathLibrary()
		self.resultSNPs = MinimalPathLibrary()

	'''Set-Pathing'''

	def setTargetDir(self, targetDir : str=None):
		if targetDir is None:
			self.targetDir = self.commonGroups.withBackup
		elif pIsAbs(targetDir):
			self.targetDir = Path(targetDir)
		else:
			self.targetDir = self.commonGroups.W > targetDir
	
	def setRefDir(self, reference : str, refDir : str=None):
		"""Reference directory is set differently than the rest of the directory groups, as it needs to both be
		specific to the .fna file containing directory, but can also look for that directory in more directories than
		those hardcoded by adding it through the --refDir flag or as 'refDir = ' in a .TOML."""

		if pIsAbs(reference):
			self.refDir = DirectoryPath(reference, purpose="r")
		else:
			if refDir is None:
				self.refDir = self.commonGroups.shared > SOFTWARE_NAME+"-Data" > "References"
			elif pIsAbs(refDir):
				self.refDir = DirectoryPath(refDir, reference, purpose="r")
			else:
				self.refDir = self.workDir > refDir
	
	def setDatabaseDir(self, databaseDir : str=None):
		if databaseDir is None:
			self.databaseDir = self.commonGroups.shared > SOFTWARE_NAME+"-Data" > "Databases"
		elif pIsAbs(databaseDir):
			self.databaseDir = PathGroup(databaseDir, purpose="r")
		else:
			self.databaseDir = DirectoryPath(self.workDir, purpose="r") > databaseDir

	def setTmpDir(self, tmpDir : str=None):
		"""Will by default create randomized unique temporary directories in the userDirectory as a first option, and in
		the workDirectory as a second. If given an absolute path temporaryDirectory, all other directories created in
		that temporary directory will be removed after the reference to that new directory is trash-collected. If given
		a relative path, the path will first be looked for in the userDirectory and secondly in the workDirectory,
		removes child paths once they are trash collected the same way as for absolute paths."""
		if tmpDir is None:
			self.tmpDir = createTemp(self.commonGroups.personal.writeable, prefix="tmp-[MetaCanSNPer]")
		else:
			tmpDir = DirectoryPath(tmpDir, purpose="w")
			if pIsAbs(tmpDir):
				self.tmpDir = createTemp(tmpDir, prefix="tmp-[MetaCanSNPer]")
			else:
				self.tmpDir = createTemp(self.commonGroups.personal > tmpDir, prefix="tmp-[MetaCanSNPer]")
	
	def setOutDir(self, outDir : str=None):
		if outDir is None:
			self.outDir = self.commonGroups.locals > SOFTWARE_NAME+"-Results"
			self.outDir.defaultPurpose = "w"
		elif pIsAbs(outDir):
			self.outDir = Path(outDir, purpose="w")
		else:
			self.outDir = self.commonGroups.locals > outDir
			self.outDir.defaultPurpose = "w"

		self.resultDir, self.logDir = self.outDir.create(self.sessionName)
		if self.resultDir is None:
			raise PermissionError(f"No directory with writing permissions found! Looked through: {self.outDir}")

	'''Set-Values'''
	
	def setSessionName(self, name):
		self.sessionName = name

	def setQuery(self, query : list[str]):
		self.query = FileList
		for q in query:
			if pIsAbs(q):
				self.query.append(FilePath(q))
				self.access(q, mode="r")
			else:
				self.query.append(self.targetDir.find(q))
		self.queryName = fileNameAlign(*map(os.path.basename, self.query))
	
	def setReferences(self, references : list[str,str,str,str,str]):
		for genome, strain, genbank_id, refseq_id, assembly_name in references:
			filename = DownloadQueue.download(genbank_id, refseq_id, assembly_name, dst=self.refDir.writable, filename=f"{assembly_name}.fna")
			if not os.path.exists(filename):
				msg = f"Could not download reference genome: [{genbank_id=}, {refseq_id=}, {assembly_name=}]"
				LOGGER.error(msg)
				raise FileNotFoundError(msg)
			self.references[genome] = filename
			self.chromosomes = [open(filename, "r").readline()[1:].split()[0]]

	def setMaps(self, maps : dict[str,str]):
		if type(maps) is MinimalPathLibrary:
			self.alignments = maps
		else:
			for r,path in maps.items():
				self.access(path, mode="r")
				self.maps[r] = path

	def setAlignments(self, alignments : dict[str,str]):
		if type(alignments) is MinimalPathLibrary:
			self.alignments = alignments
		else:
			for r,path in alignments.items():
				self.access(path, mode="r")
				self.alignments[r] = path
	
	def setTargetSNPs(self, targetSNPs : dict[str,str]):
		if type(targetSNPs) is MinimalPathLibrary:
			self.targetSNPs = targetSNPs
		else:
			for r,path in targetSNPs.items():
				self.access(path, mode="r")
				self.targetSNPs[r] = path

	def setResultSNPs(self, resultSNPs : dict[str,str]):
		if type(resultSNPs) is MinimalPathLibrary:
			self.resultSNPs = resultSNPs
		else:
			for r,path in resultSNPs.items():
				self.access(path, mode="r")
				self.resultSNPs[r] = path


if __name__ == "__main__":
	p1 = DirectoryPath(os.path.abspath("."))
	p2 = DirectoryPath(pExpUser("~"))
	print("p1 = ", p1)
	print("p2 = ", p2)

	dg1 = p1 | p2
	print("p1 | p2 = ", dg1)

	p3 = p2 > "Documents"
	print("p2 > \"Documents\" = ", p3)

	p4 = p1 + "-Data"
	print("p1 + \"-Data\" = ", p4)

	dg2 = p1 | p2 | p3 | p4
	print("p1 | p2 | p3 | p4 = ", dg2)

	print("dg1 | dg2 = ", dg1 | dg2)

	ex1 = p2 > "myReads.fq"+".log"
	print(ex1)