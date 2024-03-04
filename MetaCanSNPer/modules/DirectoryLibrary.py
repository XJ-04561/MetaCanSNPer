

import os
from functools import cached_property, cache
import shutil
import random
random.seed()
import logging
try:
	import MetaCanSNPer.modules.LogKeeper as LogKeeper
	from MetaCanSNPer.modules.DownloadReferences import DownloadQueue
	from MetaCanSNPer.modules.VCFhandler import CreateVCF, ReadVCF
	from MetaCanSNPer.modules.FileNameAlignment import align as fileNameAlign
except:
	import LogKeeper as LogKeeper
	from DownloadReferences import DownloadQueue
	from VCFhandler import openVCF, RowDict


LOGGER = LogKeeper.createLogger(__name__)
PERMS_LOOKUP = {"r":"read", "w":"write", "x":"execute"}
PERMS_LOOKUP_OS = {"r":os.R_OK, "w":os.W_OK, "x":os.X_OK}
MAX_RESULTS = 10**9
SOFTWARE_NAME = "MetaCanSNPer"

# OS alibis
pSep = os.path.sep
pJoin = os.path.join
pIsAbs = lambda *args, **kwargs: os.path.isabs(os.path.expanduser(*args, **kwargs))
pExpUser = os.path.expanduser
pAbs = os.path.abspath
pNorm = os.path.normpath
pDirName = os.path.dirname
pMakeDirs = lambda x : [os.makedirs(x), LOGGER.debug("os.makedirs('{}')".format(x))][0]
pIsFile = os.path.isfile

''' Functions '''

def findValidPath(paths, mode="r") -> str:
	for d in paths:
		if all(os.access(d, PERMS_LOOKUP_OS[c]) for c in mode):
			return d
	return None

def newRandomName(path : str, ext : str=None):
	newName = "tmp_{}{}{}".format(*random.random().as_integer_ratio(), "."+ext if ext is not None else "")
	n=0
	while os.path.exists(pJoin(path, newName)):
		newName = "tmp_{}{}{}".format(*random.random().as_integer_ratio(), "."+ext if ext is not None else "")
		if n>1000000-1:
			raise FileExistsError("Could not generate random file/directory name for path='{path}' after {n} attempts".format(path=path, n=n))
			break
	return Path(path, newName)

class Path(str):
	def __init__(self) -> None:
		LOGGER.debug("<Path '{}'>".format(self))

	def __new__(cls, *paths):
		obj = super(Path, cls).__new__(cls, pJoin(*paths))
		return obj

	def __add__(self, right):
		return Path(str.__add__(self.rstrip(pSep), right)) if not pIsFile(right) else FilePath(str.__add__(self.rstrip(pSep), right))

	def __gt__(self, right):
		return Path(self, right) if not pIsFile(right) else FilePath(self, right)
	
	def __lt__(self, left):
		return Path(left, self)

	def __or__(self, right):
		if type(right) is DirectoryGroup:
			return DirectoryGroup(self, *right._roots, purpose=right.defaultPurpose)
		else:
			return DirectoryGroup(self, right)
	
	def __ror__(self, left):
		if type(left) is DirectoryGroup:
			return DirectoryGroup(*left._roots, self, purpose=left.defaultPurpose)
		else:
			return DirectoryGroup(left, self)

class FilePath(Path):
	def __init__(self) -> None:
		LOGGER.debug("<FilePath '{}'>".format(self))
	def __lt__(self, right):
		raise TypeError("Attempted to append path to a filepath, this behovior is not currently supported (Allowed). Attempted to combine '{}' with '{}'".format(self, right))

class DisposablePath(Path):
	def __del__(self):
		try:
			LOGGER.debug("shutil.rmtree('{}', ignore_errors=True)".format(self))
			shutil.rmtree(self, ignore_errors=True)
		except:
			LOGGER.warning("Failed to remove content of directory: '{}'".format(self))

class DirectoryGroup:
	_roots : list[Path]

	def __init__(self, *paths : tuple[str], purpose="r", selfCleaning=False):
		self.defaultPurpose = purpose
		for p in paths:
			if not os.path.exists(p):
				try:
					os.makedirs(p)
				except:
					pass
		self._roots = [p if type(p) is Path else Path(p) for p in paths]
		self.selfCleaning = selfCleaning

	def __or__(self, right):
		if type(right) is DirectoryGroup:
			return DirectoryGroup(*self._roots, *right._roots, purpose=right.defaultPurpose)
		else:
			return DirectoryGroup(*self._roots, right, purpose=self.defaultPurpose)
	
	def __ror__(self, left):
		if type(left) is DirectoryGroup:
			return DirectoryGroup(*left._roots, *self._roots, purpose=left.defaultPurpose+self.defaultPurpose)
		else:
			return DirectoryGroup(left, *self._roots, purpose=self.defaultPurpose)
		
	def __ior__(self, right):
		if type(right) is DirectoryGroup:
			self._roots += right._roots
		else:
			self._roots.append(right)
	
	def __gt__(self, right):
		return DirectoryGroup(*[r > right for r in self._roots], purpose=self.defaultPurpose)
	
	def __ge__(self, right):
		self._roots = [r > right for r in self._roots]
	
	def __contains__(self, path : str):
		for r in self._roots:
			if os.path.exists(pJoin(r, path)):
				return True
		return False
	
	def __iter__(self):
		return self._roots

	def __str__(self) -> str:
		return "<DirectoryGroup>\n" + "\n".join(["  {}:d{}{}{}".format(r, *[c if os.access(r, PERMS_LOOKUP_OS[c]) else "-" for c in "rwx"]) for r in self._roots]) + "\n</DirectoryGroup>"

	def __getitem__(self, path : str, purpose:str=None) -> str:
		if purpose is None:
			purpose = self.defaultPurpose
		for r in self._roots:
			if os.path.exists(pJoin(r, path)):
				if all(os.access(pJoin(r, path), PERMS_LOOKUP_OS[p]) for p in purpose):
					return Path(pJoin(r, path))
		return None
	
	def find(self, path : str, purpose:str=None):
		'''Looks for path in the group of directories and returns first found path.'''
		return self.__getitem__(self, path, purpose=purpose if purpose is not None else self.defaultPurpose)

	def forceFind(self, path : str, purpose:str=None):
		'''Looks for path in the group of directories and returns first found path. Will try to create and return path
		in the group if it does not currently exist.'''
		return self.__getitem__(self, path, purpose=purpose) or self.create(path=path)

	
	@property
	def writable(self): return self.__getitem__(self, "", purpose="w")
	@property
	def readable(self): return self.__getitem__(self, "", purpose="r")
	
	def create(self, path : str):
		'''Should not be used to create files, only directories!'''
		for r in self._roots:
			if os.access(r, os.W_OK):
				os.makedirs(pJoin(r, pDirName(path)))
				return DisposablePath(r, pDirName(path)) if self.selfCleaning else Path(r, pDirName(path))
		return None
		

'''Container and handler of directories and files'''
class DirectoryLibrary:
	'''
		This structure allows for defining file managing as well and making sure that directory and file path
		information is simply passed around by handing around the same 'Library'.
	'''

	__lib__ : dict

	workDir : Path
	installDir : Path
	userDir : Path

	targetDir : DirectoryGroup
	refDir : DirectoryGroup
	databaseDir : DirectoryGroup
	tmpDir : DirectoryGroup
	outDir : DirectoryGroup
	resultDir : Path
	logDir : Path
	query : list[FilePath]
	
	sessionName : str
	queryName : str

	maps : dict[tuple[str, str], FilePath]
	'''Dictionary keys are tuples of (reference, query). Dictionary values are the absolute paths to the file.'''
	alignments : dict[tuple[str, str], FilePath]
	'''Dictionary keys are tuples of (reference, query). Dictionary values are the absolute paths to the file.'''
	targetSNPS : dict[str,FilePath]
	'''Dictionary keys are references' genome names. Dictionary values are the absolute paths to the file.'''
	SNPs : dict[tuple[str, str], FilePath]
	'''Dictionary keys are tuples of (reference, query). Dictionary values are the absolute paths to the file.'''
	

	baseDirs = [
		"workDir", "installDir", "userDir",
		"targetDir", "refDir", "tmpDir",
		"outDir", "resultDir", "logDir",
		"databaseDir"
	]

	@cached_property
	def workDir(self):		return Path(os.curdir)
	@cached_property
	def installDir(self):	return Path(pNorm(pJoin(pDirName(__file__), "..")))
	@cached_property
	def userDir(self):		return Path(pExpUser("~"))
	@property
	def resultDir(self):	return self.outDir.forceFind(self.sessionName)
	@property
	def logDir(self):		return self.outDir.forceFind(self.sessionName)
	@property
	def queryName(self):	return fileNameAlign(*self.query)
	
	def __init__(self, settings : dict, reference : str=None, sessionName="Unnamed-Session", **kwargs):
		'''Directories that can be passed as kwargs:
		workDir, targetDir, tmpDir, refDir, databaseDir, outDir

		Needs only a settings dictionary implementing all keys specified in the 'defaultFlags.toml' that should be in
		the installation directory, but keyword arguments can be used to force use directories.
		'''

		LOGGER.debug("Creating: {}".format(self))

		self.settings = settings
		self.sessionName = sessionName

		LOGGER.debug("Work dir:    '{}'".format(self.workDir))
		LOGGER.debug("Install dir: '{}'".format(self.installDir))
		LOGGER.debug("User dir:    '{}'".format(self.userDir))
		
		# Use kwargs as a first, and settings as a second.
		workDir = kwargs.get("workDir")
		targetDir = kwargs.get("targetDir") or self.settings.get("targetDir")
		refDir = kwargs.get("refDir") or self.settings.get("refDir")
		databaseDir = kwargs.get("databaseDir") or self.settings.get("databaseDir")
		tmpDir = kwargs.get("tmpDir") or self.settings.get("tmpDir")
		outDir = kwargs.get("outDir") or self.settings.get("outDir")

		if workDir is not None: os.chdir(workDir)
		self.setTargetDir(targetDir)
		if reference is not None: self.setRefDir(reference, refDir=refDir)
		self.setDatabaseDir(databaseDir)
		self.setTmpDir(tmpDir)
		self.setOutDir(outDir)
		# If argument is None, then a default set of directories are used. Non-abs paths are appended to a set of Paths (workDir and userDir mostly)
		
		LOGGER.debug(str(self))
		
		self.query = []

		self.indexed = {}
		self.maps = {}
		self.targetSNPS = {}
		self.SNPs = {}

	def __getitem__(self, key):
		return self.__getattribute__(key)

	def __repr__(self):
		''''''
		return "<{}.{} at 0x{:0>16}>".format(__name__, type(self).__name__, hex(id(self))[2:])
	
	def __str__(self):
		''''''
		# Works nicely, don't question it.
		return "Directories in Library at 0x{:0>16}:\n".format(hex(id(self))[2:])+"\n".join(["  {:<20}='{}{}{}' '{}'".format(d,*["?" if self[d] not in self.perms else c if self.perms[self[d]][c] else " " for c in "rwx" ], self[d]) for d in (a for a in dir(self) if a in self.baseDirs)])

	'''Set-Pathing'''

	def setTargetDir(self, targetDir : str=None):
		if targetDir is None:
			self.targetDir = DirectoryGroup(self.workDir, self.installDir, self.userDir)
		elif pIsAbs(targetDir):
			self.targetDir = DirectoryGroup(targetDir)
		else:
			self.targetDir = DirectoryGroup(self.workDir) > targetDir
	
	def setRefDir(self, reference : str, refDir : str=None):
		"""Reference directory is set differently than the rest of the directory groups, as it needs to both be
		specific to the .fna file containning directory, but can also look for that directory in more directories than
		those hardcoded by adding it through the --refDir flag or as 'refDir = ' in a .TOML."""

		self.refDir = DirectoryGroup(purpose="r")
		if refDir is not None and pIsAbs(refDir):
			self.refDir |= refDir
		elif refDir is not None:
			self.refDir |= self.workDir > refDir
		else:
			self.refDir |= self.installDir > SOFTWARE_NAME+"-Data" > "References"
			self.refDir |= self.userDir > SOFTWARE_NAME+"-Data" > "References"
			self.refDir |= self.workDir > SOFTWARE_NAME+"-Data" > "References"

		if pIsAbs(reference):
			self.refDir = DirectoryGroup(reference, purpose="r")
		else:
			self.refDir >= reference
	
	def setDatabaseDir(self, databaseDir : str=None):
		if databaseDir is None:
			self.databaseDir = DirectoryGroup(self.installDir, self.userDir, self.workDir, purpose="r") > SOFTWARE_NAME+"-Data" > "Databases"
		elif pIsAbs(databaseDir):
			self.databaseDir = DirectoryGroup(databaseDir, purpose="r")
		else:
			self.databaseDir = DirectoryGroup(self.workDir, purpose="r") > databaseDir

	def setTmpDir(self, tmpDir : str=None):
		"""Will by default create randomized unique temporary directories in the userDirectory as a first option, and in
		the workDirectory as a second. If given an absolute path temporaryDirectory, all other directories created in
		that temporary directory will be removed after the reference to that new directory is trash-collected. If given
		a relative path, the path will first be looked for in the userDirectory and secondly in the workDirectory,
		removes child paths once they are trash collected the same way as for absolute paths."""
		if tmpDir is None:
			self.tmpDir = DirectoryGroup(DisposablePath(newRandomName(self.userDir)), DisposablePath(newRandomName(self.workDir)), purpose="w")
		elif pIsAbs(tmpDir):
			self.tmpDir = DirectoryGroup(tmpDir, purpose="w", selfCleaning=True if not self.settings["saveTemp"] else False)
		else:
			self.tmpDir = DirectoryGroup(self.userDir, self.workDir, purpose="w", selfCleaning=True if not self.settings["saveTemp"] else False) > tmpDir
	
	def setOutDir(self, outDir : str=None):
		if outDir is None:
			self.outDir = DirectoryGroup(self.workDir, self.userDir, purpose="w")
		elif pIsAbs(outDir):
			self.outDir = DirectoryGroup(outDir, purpose="w")
		else:
			self.outDir = DirectoryGroup(self.workDir, self.userDir, purpose="w") > outDir


	'''Set-Values'''
	
	def setSessionName(self, name):
		self.sessionName = name

	def setQuery(self, query : list[str]):
		self.query = []
		for q in query:
			if pIsAbs(q):
				self.query.append(FilePath(q))
				self.access(q, mode="r")
			else:
				self.query.append(self.targetDir[q])
	
	def setReferences(self, references : list[str,str,str,str,str]):
		self.references = {}
		for genome, strain, genbank_id, refseq_id, assembly_name in references:
			filename = DownloadQueue.download(genbank_id, refseq_id, assembly_name, dst=self.refDir.writable)
			if not os.path.exists(filename):
				msg = "Could not download reference genome: {genbank_id='{genbank_id}', refseq_id='{refseq_id}' assembly_name='{assembly_name}'}".format(genbank_id=genbank_id, refseq_id=refseq_id, assembly_name=assembly_name)
				LOGGER.error(msg)
				raise FileNotFoundError(msg)
			self.references[genome] = filename

	def setIndexes(self, indexes : dict[tuple[str,str],str]):
		for (r,q),iP in indexes.items():
			self.access(iP)
		self.indexes = indexes
	
	def setTargetSNPs(self, targetSNPs):
		self.targetSNPs = targetSNPs

	def setSNPs(self, SNPs : dict[str,str]):
		for (r,q), path in SNPs.items():
			self.access(path, mode="r")
		self.SNPs = SNPs
	
	'''Get values'''

	def getReferences(self) -> dict[str,(str, str, str, str, str)]:
		'''Returns current list of references as a dictionary:
			{GENOME_NAME : (REFERENCE_PATH, STRAIN, GENBANK_ID, REFSEQ_ID, ASSEMBLY_NAME)}'''
		
		return self.references
	
	def getSNPdata(self) -> dict[int,tuple[str,None]]:
		'''getSNPdata() -> {POS : CALLED}
		'''
		SNPs = {}
		for (refName, queryName), path in self.SNPs.items():
			if path.lower().endswith(".vcf"):
				reader = openVCF(path, "r")
				for entry in reader:
					SNPs[entry.POS] = (entry.REF, None) # Called Base is in the "REF" field
					# More fields might be added instead of None.
		
		return SNPs

	def access(self, path, mode : str="rwx", create : bool=False):
		'''Throws appropriate errors if access is not possible.'''
		
		if not os.path.exists(path):
			if create:
				LOGGER.debug("Path does not exist, creating it instead: '{}'".format(path))
				os.makedirs(path)
			else:
				LOGGER.error("Path does not exist: '{}'".format(path))
				raise FileNotFoundError("Path does not exist: '{}'".format(path))
		
		if not all(os.access(path, PERMS_LOOKUP_OS[c]) for c in mode):
			LOGGER.error("Missing {perms} permissions for path: '{path}'".format(path=path, perms="+".join([PERMS_LOOKUP[c] for c in mode if not self.perms[path][c]])))
			raise PermissionError("Missing {perms} permissions for path: '{path}'".format(path=path, perms="+".join([PERMS_LOOKUP[c] for c in mode if not self.perms[path][c]])))
		
		return True
	
	def accessible(self, path, mode : str="rwx", create : bool=False):
		'''Returns True/False for the accessibility question.'''
		
		if not os.path.exists(path):
			if create:
				LOGGER.debug("Path does not exist, creating it instead: '{}'".format(path))
				os.makedirs(path)
			else:
				LOGGER.debug("Path does not exist: '{}'".format(path))
				return False
		
		if not all(os.access(path, PERMS_LOOKUP_OS[c]) for c in mode):
			LOGGER.debug("Missing {perms} permissions for path: '{path}'".format(path=path, perms="+".join([PERMS_LOOKUP[c] for c in mode if not self.perms[path][c]])))
			return False
		
		return True

if __name__ == "__main__":
	p1 = Path(os.path.abspath("/"))
	p2 = Path(pExpUser("~"))
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