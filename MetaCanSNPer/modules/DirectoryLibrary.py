

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

	settings : dict
	query : FileList[FilePath]
	sessionName : str
	queryName : str
	chromosomes : list[str]
	
	def __init__(self, settings : dict, sessionName="Unnamed-Session", **kwargs):
		"""Directories that can be passed as kwargs:
		workDir, targetDir, tmpDir, refDir, databaseDir, outDir

		Needs only a settings dictionary implementing all keys specified in the 'defaultFlags.toml' that should be in
		the installation directory, but keyword arguments can be used to force use directories.
		"""

		LOGGER.info("Creating DirectoryLibrary object.")

		super(DirectoryLibrary, self).__init__(self)

		object.__setattr__(self, "settings", {})
		object.__setattr__(self, "sessionName", sessionName)

		self.updateSettings(settings | kwargs)

		LOGGER.debug(f"Work dir:    {self.workDir!r}")
		LOGGER.debug(f"Install dir: {self.installDir!r}")
		LOGGER.debug(f"User dir:    {self.userDir!r}")
		
		PseudoPathy.Globals.DISPOSE = self.settings.get("saveTemp", True)
		
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
		LOGGER.debug(f"Set targetDir to:\n{self.targetDir}")
	
	def setRefDir(self, organism : str, refDir : str=None):
		"""Reference directory is set differently than the rest of the directory groups, as it needs to both be
		specific to the .fna file containing directory, but can also look for that directory in more directories than
		those hardcoded by adding it through the --refDir flag or as 'refDir = ' in a .TOML."""

		if pIsAbs(organism):
			self.refDir = DirectoryPath(organism, purpose="r")
		else:
			if refDir is None:
				self.refDir = (self.commonGroups.shared > (SOFTWARE_NAME+"-Data")) > "References"
			elif pIsAbs(refDir):
				self.refDir = DirectoryPath(refDir, organism, purpose="r")
			else:
				self.refDir = self.workDir > refDir
		LOGGER.debug(f"Set refDir to:\n{self.refDir}")
	
	def setDatabaseDir(self, databaseDir : str=None):
		print(databaseDir)
		if databaseDir is None:
			self.databaseDir = (self.commonGroups.shared > (SOFTWARE_NAME+"-Data")) > "Databases"
		elif pIsAbs(databaseDir):
			self.databaseDir = DirectoryPath(databaseDir, purpose="r")
		else:
			self.databaseDir = self.workDir > databaseDir
		LOGGER.debug(f"Set databaseDir to:\n{self.databaseDir}")

	def setTmpDir(self, tmpDir : str=None):
		"""Will by default create randomized unique temporary directories in the userDirectory as a first option, and in
		the workDirectory as a second. If given an absolute path temporaryDirectory, all other directories created in
		that temporary directory will be removed after the reference to that new directory is trash-collected. If given
		a relative path, the path will first be looked for in the userDirectory and secondly in the workDirectory,
		removes child paths once they are trash collected the same way as for absolute paths."""
		if tmpDir is None:
			self.tmpDir = createTemp(self.commonGroups.personal.writable > SOFTWARE_NAME, prefix="tmp-[MetaCanSNPer]")
		else:
			tmpDir = DirectoryPath(tmpDir, purpose="w")
			if pIsAbs(tmpDir):
				self.tmpDir = createTemp(tmpDir, prefix="tmp-[MetaCanSNPer]")
			else:
				self.tmpDir = createTemp(self.commonGroups.personal > tmpDir, prefix="tmp-[MetaCanSNPer]")
		LOGGER.debug(f"Set tmpDir to:\n{self.tmpDir}")
	
	def setOutDir(self, outDir : str=None):
		if outDir is None:
			self.outDir = self.commonGroups.locals > SOFTWARE_NAME+"-Results"
			self.outDir.defaultPurpose = "w"
		elif pIsAbs(outDir):
			self.outDir = Path(outDir, purpose="w")
		else:
			self.outDir = self.commonGroups.locals > outDir
			self.outDir.defaultPurpose = "w"
		LOGGER.debug(f"Set outDir to:\n{self.outDir}")

		self.resultDir = self.logDir = self.outDir.create(self.sessionName)
		LOGGER.debug(f"Set resultDir and logDir to:\n{self.resultDir}")
		if self.resultDir is None:
			raise PermissionError(f"No directory with writing permissions found! Looked through: {self.outDir}")

	'''Set-Values'''
	
	def updateSettings(self, settings : dict):
		new = dict(set(settings.items()).difference(self.settings.items()))
		
		LOGGER.debug(f"Updating settings with:\n{new}")
		self.settings.update(settings)
		
		if new.get("workDir") is not None:		os.chdir(self.settings.get("workDir"))
		if new.get("userDir") is not None:		self.userDir = self.settings.get("userDir")
		if new.get("installDir") is not None:	self.installDir = self.settings.get("installDir")

		if "targetDir" in new:				self.setTargetDir(self.settings.get("targetDir"))
		if "organism" in new:				self.setRefDir(self.settings.get("organism"), refDir=self.settings.get("refDir"))
		elif "refDir" in new and \
			"reference" in self.settings:	self.setRefDir(self.settings.get("reference"), refDir=self.settings.get("refDir"))
		if "databaseDir" in new:			self.setDatabaseDir(self.settings.get("databaseDir"))
		if "outDir" in new:					self.setOutDir(self.settings.get("outDir"))
		if "tmpDir" in new:					self.setTmpDir(self.settings.get("tmpDir"))

	def setSessionName(self, name):
		self.sessionName = name
		LOGGER.debug(f"Setting sessionName to:{self.sessionName}")

	def setQuery(self, query : list[str]):
		self.query = FileList()
		for q in query:
			if pIsAbs(q):
				self.query.append(FilePath(q))
				self.access(q, mode="r")
			else:
				self.query.append(self.targetDir.find(q))
				if self.query[-1] is None:
					LOGGER.error(f"Query file {q!r} could not be found.")
					raise FileNotFoundError(f"Query file {q!r} could not be found.")
		LOGGER.debug(f"Setting query to:{self.query}")
		self.queryName = fileNameAlign(*[pName(q) for q in self.query])
		LOGGER.debug(f"Setting queryName to:{self.queryName}")
	
	def setReferences(self, references : list[str,str,str,str,str]):
		self.references = MinimalPathLibrary()
		LOGGER.info(f"Downloading references:{references}")
		for genome, strain, genbank_id, refseq_id, assembly_name in references:
			filename = DownloadQueue.download(genbank_id, refseq_id, assembly_name, dst=self.refDir.writable, filename=f"{assembly_name}.fna")
			self.references[genome] = filename
		n=1
		for genome, strain, genbank_id, refseq_id, assembly_name in references:
			LOGGER.info(f"Waiting for downloads to complete {n}/{len(references)}")
			filename = DownloadQueue.waitNext()
			if not os.path.exists(filename):
				msg = f"Could not download reference genome: [{genbank_id=}, {refseq_id=}, {assembly_name=}] -> {filename!r}"
				LOGGER.error(msg)
				raise FileNotFoundError(msg)
			self.chromosomes = [open(filename, "r").readline()[1:].split()[0]]
			n+=1
		LOGGER.debug(f"Finished downloading {len(references)} reference genomes!")
		
		LOGGER.debug(f"Setting references to:{self.references}")

	def setMaps(self, maps : dict[str,str]):
		if type(maps) is MinimalPathLibrary:
			self.maps = maps
		else:
			self.maps = MinimalPathLibrary()
			for r,path in maps.items():
				self.access(path, mode="r")
				self.maps[r] = path
		LOGGER.debug(f"Setting maps to:{self.maps}")

	def setAlignments(self, alignments : dict[str,str]):
		if type(alignments) is MinimalPathLibrary:
			self.alignments = alignments
		else:
			self.alignments = MinimalPathLibrary()
			for r,path in alignments.items():
				self.access(path, mode="r")
				self.alignments[r] = path
		LOGGER.debug(f"Setting alignments to:{self.alignments}")
	
	def setTargetSNPs(self, targetSNPs : dict[str,str]):
		if type(targetSNPs) is MinimalPathLibrary:
			self.targetSNPs = targetSNPs
		else:
			self.targetSNPs = MinimalPathLibrary()
			for r,path in targetSNPs.items():
				self.access(path, mode="r")
				self.targetSNPs[r] = path
		LOGGER.debug(f"Setting targetSNPs to:{self.targetSNPs}")

	def setResultSNPs(self, resultSNPs : dict[str,str]):
		if type(resultSNPs) is MinimalPathLibrary:
			self.resultSNPs = resultSNPs
		else:
			self.resultSNPs = MinimalPathLibrary()
			for r,path in resultSNPs.items():
				self.access(path, mode="r")
				self.resultSNPs[r] = path
		LOGGER.debug(f"Setting resultSNPs to:{self.resultSNPs}")


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