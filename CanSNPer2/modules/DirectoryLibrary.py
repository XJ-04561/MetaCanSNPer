

import os
from functools import cached_property, cache
import shutil
import random
random.seed()
import logging
try:
	import CanSNPer2.modules.LogKeeper as LogKeeper
	from CanSNPer2.modules.DownloadReferences import DownloadQueue
except:
	import LogKeeper as LogKeeper
	from DownloadReferences import DownloadQueue
import vcf as vcf

LOGGER = LogKeeper.createLogger(__name__)
PERMS_LOOKUP = {"r":"read", "w":"write", "x":"execute"}
PERMS_LOOKUP_OS = {"r":os.R_OK, "w":os.W_OK, "x":os.X_OK}
MAX_RESULTS = 10**9

''' Functions '''

def findValidPath(self, paths, mode="r") -> str:
	for d in paths:
		if all(os.access(d, PERMS_LOOKUP_OS[c]) for c in mode):
			return d
	return None


class DirectoryGroup:
	_roots : list[str]

	def __init__(self, *paths : tuple[str], purpose="r"):
		self.defaultPurpose = purpose
		for p in paths:
			if not os.path.exists(p):
				os.makedirs(p)
		self._roots = [*paths]
	
	def __contains__(self, path : str):
		for r in self._roots:
			if os.path.exists(os.path.join(r, path)):
				return True
		return False

	def __str__(self) -> str:
		return "{}\t: {}".format(self._roots, [c if os.access(r, PERMS_LOOKUP_OS[c]) else " " for r in self._roots for c in "rwx"])

	def __getitem__(self, path : str, purpose=None) -> str:
		if purpose is None:
			purpose = self.defaultPurpose
		for r in self._roots:
			if not os.path.exists(os.path.join(r, path)):
				self.create(path)
			
			if os.access(os.path.join(r, path), PERMS_LOOKUP_OS[purpose]):
				return os.path.join(r, path)
		return None
	
	def create(self, path : str):
		'''Should not be used to create files, only directories!'''
		for r in self._roots:
			if os.access(r, os.W_OK):
				os.makedirs(os.path.join(r, os.path.dirname(path)))
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

	@cached_property
	def workDir(self):
		return os.curdir
	@cached_property
	def installDir(self):
		return os.path.dirname(__file__)
	@cached_property
	def userDir(self):
		return os.path.expanduser("~")

	targetDir : DirectoryGroup
	refDir : DirectoryGroup
	databaseDir : DirectoryGroup
	tmpDir : DirectoryGroup
	outDir : DirectoryGroup
	resultDir : str
	logDir : str

	@property
	def resultDir(self):
		return self.outDir["Results"]
	@property
	def logDir(self):
		return self.outDir["Logs"]

	dirs = [
		"workDir", "installDir", "userDir",
		"targetDir", "refDir", "tmpDir",
		"outDir", "resultDir", "logDir",
		"databaseDir"
	]

	query : list[str]
	indexed : dict[tuple[str, str], str]
	maps : dict[tuple[str, str], str]
	SNPs : dict[tuple[str, str], str]

	'''Dictionary keys are tuples of (reference, query).
	Dictionary values are the absolute paths to the file.'''

	__cache : dict[str]
	
	def __init__(self, settings : dict, workDir : str=None, installDir : str=None, userDir : str=None,
			     targetDir : str=None, tmpDir : str=None, refDir : str=os.path.join("References", "Francisella_Tularensis"),
				 databaseDir : str="Databases", outDir : str=None):
		'''Needs only a settings dictionary implementing all keys specified in the 'defaultFlags.toml' that should be in
		the installation directory, but keyword arguments can be used to force use directories.
		'''

		LOGGER.debug("Creating: {}".format(self))

		self.settings = settings

		if workDir is not None:
			os.chdir(self.workDir)
		
		LOGGER.debug("Work dir:    '{}'".format(self.workDir))
		LOGGER.debug("Install dir: '{}'".format(self.installDir))
		LOGGER.debug("User dir:    '{}'".format(self.userDir))
		
		# Use a list of directories in order of priority, unless being forced via flags
		self.targetDir = DirectoryGroup(*[self.workDir, self.installDir, self.userDir] if targetDir is None else targetDir)
		self.refDir = DirectoryGroup(*[os.path.join(d, "References") for d in [self.installDir, self.userDir, self.workDir]] if refDir is None else refDir)
		self.databaseDir = DirectoryGroup(*[os.path.join(d, "Databases") for d in [self.installDir, self.userDir, self.workDir]] if databaseDir is None else databaseDir)
		self.tmpDir = DirectoryGroup(*[self.newRandomName(d) for d in [self.userDir, self.workDir, self.installDir]] if tmpDir is None else tmpDir)
		self.outDir = DirectoryGroup(*[self.workDir, self.userDir] if outDir is None else outDir)

		LOGGER.debug(str(self))
		
		self.query = None

		self.indexed = {}
		self.maps = {}
		self.SNPs = {}

		self.__cache = {}
	
	def __getitem__(self, key):
		return self.__getattribute__(key)

	def __repr__(self):
		''''''
		return "<{}.{} at 0x{:0>16}>".format(__name__, type(self).__name__, hex(id(self))[2:])
	
	def __str__(self):
		''''''
		# Works nicely, don't question it.
		return "Directories in Library at 0x{:0>16}:\n".format(hex(id(self))[2:])+"\n".join(["  {:<20}='{}{}{}' '{}'".format(d,*["?" if self[d] not in self.perms else c if self.perms[self[d]][c] else " " for c in "rwx" ], self[d]) for d in (a for a in dir(self) if a in self.dirs)])

	def __del__(self):
		if os.path.exists(self.tmpDir) and (not self.settings["KeepTemporaryFiles"] if "KeepTemporaryFiles" in self.settings else True):
			shutil.rmtree(self.tmpDir)
	
	def _permCheck(self, dr):
		self.perms[dr] = {}
		self.perms[dr]["r"] = os.access(dr, os.R_OK)
		self.perms[dr]["w"] = os.access(dr, os.W_OK)
		self.perms[dr]["x"] = os.access(dr, os.X_OK)
	
	@staticmethod
	def newRandomName(path : str, ext : str=None):
		newName = "tmp_{}{}{}".format(*random.random().as_integer_ratio(), "."+ext if ext is not None else "")
		n=0
		while os.path.exists(os.path.join(path, newName)):
			newName = "tmp_{}{}{}".format(*random.random().as_integer_ratio(), "."+ext if ext is not None else "")
			if n>1000000-1:
				raise FileExistsError("Could not generate random file/directory name for path='{path}' after {n} attempts".format(path=path, n=n))
				break
		return os.path.join(path, newName)

	'''Set-functions'''
	def setTargetDir(self, targetDir : str):
		if os.path.isabs(targetDir):
			self.targetDir = DirectoryGroup(targetDir, purpose="w")
		else:
			self.targetDir = DirectoryGroup(*[os.path.join(d, targetDir) for d in [self.workDir, self.installDir, self.userDir]], purpose="r")
	
	def setRefDir(self, refDir : str):
		if os.path.isabs(refDir):
			self.refDir = DirectoryGroup(refDir, purpose="w")
		else:
			self.refDir = DirectoryGroup(*[os.path.join(d, refDir) for d in [self.installDir, self.userDir, self.workDir]], purpose="r")
	
	def setDatabaseDir(self, databaseDir : str):
		if os.path.isabs(databaseDir):
			self.databaseDir = DirectoryGroup(databaseDir, purpose="w")
		else:
			self.databaseDir = DirectoryGroup(*[os.path.join(d, databaseDir) for d in [self.installDir, self.userDir, self.workDir]], purpose="r")

	def setTmpDir(self, tmpDir : str):
		if os.path.isabs(tmpDir):
			self.tmpDir = DirectoryGroup(tmpDir, purpose="w")
		else:
			self.tmpDir = DirectoryGroup(*[os.path.join(d, tmpDir) for d in [self.userDir, self.workDir, self.installDir]], purpose="w")
	
	def setOutDir(self, outDir : str):
		if os.path.isabs(outDir):
			self.outDir = DirectoryGroup(outDir, purpose="w")
		else:
			self.outDir = DirectoryGroup(*[os.path.join(d, outDir) for d in [self.workDir, self.userDir]], purpose="w")
	

	def setQuery(self, query : str, abs : bool=True):
		if os.path.isabs(query):
			self.query = query
			self.access(self.query, mode="r")
		else:
			self.query = self.targetDir[query]
	
	def setReferences(self, references : list[str,str,str,str,str]):
		self.references = {}
		for genome, strain, genbank_id, refseq_id, assembly_name in references:
			filename = DownloadReferences.download(genbank_id, refseq_id, assembly_name, dst=self.refDir)
			if not os.path.exists(filename):
				msg = "Could not download reference genome: {genbank_id='{genbank_id}', refseq_id='{refseq_id}' assembly_name='{assembly_name}'}".format(genbank_id=genbank_id, refseq_id=refseq_id, assembly_name=assembly_name)
				LOGGER.error(msg)
				raise FileNotFoundError(msg)
			self.references[genome] = filename, strain, genbank_id, refseq_id, assembly_name

	def setIndexes(self, indexes : dict[tuple[str,str],str]):
		for (q,r),iP in indexes.items():
			self.access(iP)
		self.indexes = indexes
	
	def setSNPs(self, SNPs : dict[str,str]):
		for r,SNP in SNPs.items():
			self.access(SNP, mode="r")
		self.SNPs = SNPs
	
	'''Get-functions'''

	def get(self, filename : str, hint : str=None):
		'''Gets the absolute path to a file/directory found in one of the directories used by the Library. Hint can be
		provided to limit search to only one of the directories in the Library.
		If a path can not be found, it will be returned. This should be the case for a path which already is absolute.'''

		if filename in self.__cache:
			LOGGER.debug("Returned cached path '{path}' for filename '{filename}'".format(path=self.__cache[filename], filename=filename))
			return self.__cache[filename]

		if os.path.isabs(filename):
			return filename
		elif hint is not None:
			# directory hint has been provided
			for path, dirs, files in os.walk(self[hint]):
				for f in files:
					if f == filename:
						self.__cache[filename] = os.path.join(path, f)
						return os.path.join(path, f)
				for d in dirs:
					if d == filename:
						self.__cache[filename] = os.path.join(path, d)
						return os.path.join(path, d)
					
			LOGGER.info("Looked for path: '{}' in {} but could not find it.\n{}".format(filename, hint, str(self)))
		else:
			# Without a hint, each directory is searched.
			for dr in {self[drVar] for drVar in self.dirs}:
				for path, dirs, files in os.walk(self[dr]):
					for f in files:
						if f == filename:
							return os.path.join(path, f)
					for d in dirs:
						if d == filename:
							return os.path.join(path, d)
		
			LOGGER.info("Looked for path: '{}' in all directories but could not find it.\n{}".format(filename, str(self)))
	
	def getReferences(self) -> dict[str,(str, str, str, str, str)]:
		'''Returns current list of references as a dictionary:
			{GENOME_NAME : (REFERENCE_PATH, STRAIN, GENBANK_ID, REFSEQ_ID, ASSEMBLY_NAME)}'''
		
		return self.references
	
	def gatherSNPdata(self) -> dict[tuple[str,str],tuple[int,str,str,list[str],int,dict[str,int|str|list[float]|bool],dict[str,dict[str,int|str|list[int,int]]]]]:
		'''Returns: {
			(reference, query) : (POS, ID, REF, ALT, QUAL, INFO, samples),
			...
		}
		Where:
			POS		= number	, ID		= name	, REF		= seq,
			ALT		= [seq, ]	, QUAL		= number, INFO		= {"XY" : value, ...},
			samples	= {sampleName : {"XY" : value, }, }
		'''
		SNPs = {}
		for (reference, query), path in self.SNPs.items():
			reader = vcf.Reader(filename=path)
			SNPs[(reference, query)] = [(entry.POS+1, entry.ID, entry.REF, entry.ALT, entry.QUAL, entry.INFO, {sample.sample:{key:sample.data.__getattribute__(key) for key in dir(sample.data) if key.isupper()} for sample in entry.samples}) for entry in reader]
		
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
		self._permCheck(path)
		if not all(self.perms[path][c] for c in mode):
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
		self._permCheck(path)
		if not all(self.perms[path][c] for c in mode):
			LOGGER.debug("Missing {perms} permissions for path: '{path}'".format(path=path, perms="+".join([PERMS_LOOKUP[c] for c in mode if not self.perms[path][c]])))
			return False
		
		return True

if __name__ == "__main__":
	DL = DirectoryLibrary({})

	print(DL)