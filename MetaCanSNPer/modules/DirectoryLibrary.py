

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
except:
	import LogKeeper as LogKeeper
	from DownloadReferences import DownloadQueue
	from VCFhandler import CreateVCF, ReadVCF


LOGGER = LogKeeper.createLogger(__name__)
PERMS_LOOKUP = {"r":"read", "w":"write", "x":"execute"}
PERMS_LOOKUP_OS = {"r":os.R_OK, "w":os.W_OK, "x":os.X_OK}
MAX_RESULTS = 10**9
SOFTWARE_NAME = "MetaCanSNPer"

# OS alibis
pJoin = os.path.join
pIsAbs = os.path.isabs
pExpUser = os.path.expanduser
pAbs = os.path.abspath
pNorm = os.path.normpath
pDirName = os.path.dirname

''' Functions '''

def findValidPath(self, paths, mode="r") -> str:
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
	return pJoin(path, newName)

class DirectoryGroup:
	_roots : list[str]

	def __init__(self, *paths : tuple[str], purpose="r"):
		self.defaultPurpose = purpose
		for p in paths:
			if not os.path.exists(p):
				try:
					os.makedirs(p)
				except:
					pass
		self._roots = [*paths]
	
	def __contains__(self, path : str):
		for r in self._roots:
			if os.path.exists(pJoin(r, path)):
				return True
		return False

	def __str__(self) -> str:
		return "{}\t: {}".format(self._roots, [c if os.access(r, PERMS_LOOKUP_OS[c]) else " " for r in self._roots for c in "rwx"])

	def __getitem__(self, path : str, purpose:str=None) -> str:
		if purpose is None:
			purpose = self.defaultPurpose
		for r in self._roots:
			if os.path.exists(pJoin(r, path)):
				if os.access(pJoin(r, path), PERMS_LOOKUP_OS[purpose]):
					return pJoin(r, path)
		return None
	
	def find(self, path : str, purpose:str=None):
		'''Looks for path in the group of directories and returns first found path.'''
		return self.__getitem__(self, path, purpose=purpose if purpose is not None else self.defaultPurpose)

	def forceFind(self, path : str, purpose:str=None):
		'''Looks for path in the group of directories and returns first found path. Will try to create and return path
		in the group if it does not currently exist.'''
		return self.__getitem__(self, path, purpose=purpose if purpose is not None else self.defaultPurpose) or self.create(path=path)

	
	@property
	def writable(self): return self.__getitem__(self, "", purpose="w")
	@property
	def readable(self): return self.__getitem__(self, "", purpose="r")
	
	def create(self, path : str):
		'''Should not be used to create files, only directories!'''
		for r in self._roots:
			if os.access(r, os.W_OK):
				os.makedirs(pJoin(r, pDirName(path)))
				return pJoin(r, pDirName(path))
		return None
		

'''Container and handler of directories and files'''
class DirectoryLibrary:
	'''
		This structure allows for defining file managing as well and making sure that directory and file path
		information is simply passed around by handing around the same 'Library'.
	'''

	__lib__ : dict

	workDir : str
	installDir : str
	userDir : str

	targetDir : DirectoryGroup
	refDir : DirectoryGroup
	databaseDir : DirectoryGroup
	tmpDir : DirectoryGroup
	outDir : DirectoryGroup
	resultDir : str
	logDir : str
	query : list[str]
	
	sessionName : str
	queryName : str

	indexed : dict[tuple[str, str], str]
	'''Dictionary keys are tuples of (reference, query). Dictionary values are the absolute paths to the file.'''
	maps : dict[tuple[str, str], str]
	'''Dictionary keys are tuples of (reference, query). Dictionary values are the absolute paths to the file.'''
	SNPs : dict[tuple[str, str], str]
	'''Dictionary keys are tuples of (reference, query). Dictionary values are the absolute paths to the file.'''


	baseDirs = [
		"workDir", "installDir", "userDir",
		"targetDir", "refDir", "tmpDir",
		"outDir", "resultDir", "logDir",
		"databaseDir"
	]

	@cached_property
	def workDir(self):		return os.curdir
	@cached_property
	def installDir(self):	return pNorm(pJoin(pDirName(__file__), ".."))
	@cached_property
	def userDir(self):		return pExpUser("~")
	@property
	def resultDir(self):	return self.outDir.forceFind(self.sessionName)
	@property
	def logDir(self):		return self.outDir.forceFind(self.sessionName)
	@cached_property
	def queryName(self):	return os.path.splitext(os.path.basename(self.query))[0]
	
	def __init__(self, settings : dict, sessionName="Unnamed-Session", **kwargs):
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
		refDir = kwargs.get("refDir")
		databaseDir = kwargs.get("databaseDir") or self.settings.get("databaseDir")
		tmpDir = kwargs.get("tmpDir") or self.settings.get("tmpDir")
		outDir = kwargs.get("outDir") or self.settings.get("outDir")

		# Use pre-defined directories if kwargs and settings did not have specifications.
		if workDir is not None: os.chdir(workDir)
		self.targetDir = DirectoryGroup(*[self.workDir, self.installDir, self.userDir] if targetDir is None else targetDir)
		self.refDir = DirectoryGroup(*[pJoin(d, SOFTWARE_NAME+"-Data", "References", refDir) if "refDir" not in self.settings else pJoin(d, self.settings["refDir"]) for d in [self.installDir, self.userDir, self.workDir]] if "refDir" not in self.settings or not pIsAbs(pExpUser(refDir)) else refDir)
		self.databaseDir = DirectoryGroup(*[pJoin(d, SOFTWARE_NAME+"-Data", "Databases") for d in [self.installDir, self.userDir, self.workDir]] if databaseDir is None else databaseDir)
		self.tmpDir = DirectoryGroup(*[newRandomName(d) for d in [self.userDir, self.workDir]] if tmpDir is None else tmpDir, purpose="w")
		self.outDir = DirectoryGroup(*[pJoin(d, SOFTWARE_NAME+"-Results", sessionName) for d in [self.workDir, self.userDir]] if outDir is None else pJoin(outDir, sessionName), purpose="w")

		LOGGER.debug(str(self))
		
		self.query = None

		self.indexed = {}
		self.maps = {}
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

	def __del__(self):
		if os.path.exists(self.tmpDir) and (not self.settings["KeepTemporaryFiles"] if "KeepTemporaryFiles" in self.settings else True):
			shutil.rmtree(self.tmpDir)

	'''Set-Pathing'''

	def setTargetDir(self, targetDir : str):
		if pIsAbs(targetDir):
			self.targetDir = DirectoryGroup(targetDir, purpose="r")
		else:
			self.targetDir = DirectoryGroup(*[pJoin(d, targetDir) for d in [self.workDir, self.installDir, self.userDir]], purpose="r")
	
	def setRefDir(self, refDir : str):
		if pIsAbs(refDir):
			self.refDir = DirectoryGroup(refDir, purpose="r")
		else:
			self.refDir = DirectoryGroup(*[pJoin(d, refDir) for d in [self.installDir, self.userDir, self.workDir]], purpose="r")
	
	def setDatabaseDir(self, databaseDir : str):
		if pIsAbs(databaseDir):
			self.databaseDir = DirectoryGroup(databaseDir, purpose="r")
		else:
			self.databaseDir = DirectoryGroup(*[pJoin(d, databaseDir) for d in [self.installDir, self.userDir, self.workDir]], purpose="r")

	def setTmpDir(self, tmpDir : str):
		if pIsAbs(tmpDir):
			self.tmpDir = DirectoryGroup(tmpDir, purpose="w")
		else:
			self.tmpDir = DirectoryGroup(*[pJoin(d, tmpDir) for d in [self.userDir, self.workDir, self.installDir]], purpose="w")
	
	def setOutDir(self, outDir : str=SOFTWARE_NAME+"-Results"):
		if pIsAbs(outDir):
			self.outDir = DirectoryGroup(outDir, purpose="w")
		else:
			self.outDir = DirectoryGroup(*[pJoin(d, outDir) for d in [self.workDir, self.userDir]], purpose="w")

	def setSessionName(self, name):
		self.sessionName = name
	
	'''Set-Values'''

	def setQuery(self, query : str):
		if pIsAbs(query):
			self.query = query
			self.access(self.query, mode="r")
		else:
			self.query = self.targetDir[query]
	
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
		for (q,r),iP in indexes.items():
			self.access(iP)
		self.indexes = indexes
	
	def setSNPs(self, SNPs : dict[str,str]):
		for r,SNP in SNPs.items():
			self.access(SNP, mode="r")
		self.SNPs = SNPs
	
	def getReferences(self) -> dict[str,(str, str, str, str, str)]:
		'''Returns current list of references as a dictionary:
			{GENOME_NAME : (REFERENCE_PATH, STRAIN, GENBANK_ID, REFSEQ_ID, ASSEMBLY_NAME)}'''
		
		return self.references
	
	def getSNPdata(self) -> dict[tuple[str,str],tuple[int,str,str,list[str],int,dict[str,int|str|list[float]|bool],dict[str,dict[str,int|str|list[int,int]]]]]:
		'''getSNPdata() -> {POS : CALLED}
		'''
		SNPs = {}
		for (reference, query), path in self.SNPs.items():
			reader = ReadVCF(filename=path)
			for entry in reader:
				SNPs[entry.POS] = entry.REF # Called Base is in the "REF" field
		
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
	DL = DirectoryLibrary({})

	print(DL)