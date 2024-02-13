

import os
import shutil
import random
random.seed()
import logging
try:
	import CanSNPer2.modules.LogKeeper as LogKeeper
	import CanSNPer2.modules.DownloadReferences as DownloadReferences
except:
	import LogKeeper as LogKeeper
	import DownloadReferences as DownloadReferences

LOGGER = LogKeeper.createLogger(__name__)
PERMS_LOOKUP = {"r":"read", "w":"write", "x":"execute"}

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

	targetDir : str
	refDir : str
	tmpDir : str
	outDir : str
	resultDir : str
	logDir : str
	databaseDir : str

	dirs = [
		"workDir", "installDir", "userDir",
		"targetDir", "refDir", "tmpDir",
		"outDir", "resultDir", "logDir",
		"databaseDir"
	]

	query : list[str]
	
	# Dictionary keys are tuples of (reference, query)
	# Dictionary values are the absolute paths to the file.
	indexed : dict[tuple[str, str], str]
	maps : dict[tuple[str, str], str]
	SNPs : dict[tuple[str, str], str]

	__cache : dict[str]
	
	def __init__(self, settings : dict, workDir : str=None, installDir : str=None, userDir : str=None,
			     tmpDir : str=None, refDir : str=os.path.join("References", "Francisella_Tularensis"),
				 databaseDir : str="Databases", outDir : str=None):
		'''Needs only a settings dictionary implementing all keys specified in the 'defaultFlags.toml' that should be in
		the installation directory, but keyword arguments can be used to force use directories.
		
		workDir		-	Will be used to find queries/targets, and is the second directory to be searched for references
		
		installDir	-	Will be used to store temporarily generated outputs, data, and converted data, and is the first
						place to look for references.

		userDir		-	Will be used as a backup place for temporary storage if permissions don't allow the installDir,
						and is the third place to look for references.
		'''

		LOGGER.debug("Creating: {}".format(self))

		self.settings = settings

		if workDir is None:
			self.workDir = os.getcwd()
		else:
			self.accessible(workDir)
			os.chdir(self.workDir)
		
		# Not really meant to be overriden, but implemented just in case.
		if installDir is None:
			self.installDir, _ = __file__.rsplit(os.path.sep, 1)
		else:
			self.accessible(installDir)
			self.installDir = installDir
		
		# Not really meant to be overriden, but implemented just in case.
		if userDir is None:
			self.userDir = os.path.expanduser("~")
		else:
			self.accessible(userDir)
			self.userDir = userDir

		LOGGER.debug("Work dir:    '{}'".format(self.workDir))
		LOGGER.debug("Install dir: '{}'".format(self.installDir))
		LOGGER.debug("User dir:    '{}'".format(self.userDir))

		self.perms = {}
		for d in [self.workDir, self.installDir, self.userDir]:
			self._permCheck(d)
		
		LOGGER.debug(str(self))

		self.setTargetDir(self.workDir)
		self.setRefDir(refDir)
		self.setDatabaseDir()
		self.setTmpDir(tmpDir)
		self.setOutDir(outDir)
		
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

	def setOutDir(self, outDir : str, abs=False):
		if abs is True:
			self.access(outDir, create=True)
		elif outDir is None:
			if not all(self.perms[self.workDir].values()):
				LOGGER.error("No valid working directory or output directory given.\n" + str(self))
				raise PermissionError("No valid working directory or output directory given.\n" + str(self))
			else:
				self.outDir = os.path.join(self.workDir, "MetaCanSNPer")
				if not os.path.exists(self.outDir):
					os.mkdir(self.outDir)
		else:
			self.access(os.path.join(self.workDir, outDir))
			rOutDir = self.get(outDir, "workDir")
			self._permCheck(rOutDir)
			if all(self.perms[rOutDir].values()):
				self.outDir = rOutDir
			else:
				raise PermissionError("Given outDir='{}' does not have read, write, and execute permissions. It has '{}{}{}'".format(rOutDir, *[c if self.perms[rOutDir][c] else " " for c in "rwx"]))
		self.setResultDir()
		self.setLogDir()
		
	def setResultDir(self, name : str="Results", abs=False):
		'''This method assumes that self.access(self.outDir) does not throw any errors.'''
		if abs is True:
			self.resultDir = name
			self.access(self.resultDir, create=True)
		elif "MergeNewResults" in self.settings and self.settings["MergeNewResults"]:
			self.resultDir = os.path.join(self.outDir, name)
			self.access(self.resultDir, create=True)
		else:
			n=1
			while os.path.exists(os.path.join(self.outDir, "{}-{}".format(name, str(n)))):
				n+=1
			self.resultDir = os.path.join(self.outDir, "{}-{}".format(name, str(n)))
			os.mkdir(self.resultDir)
	
	def setLogDir(self, name : str="Logs", abs=False):
		if abs is True:
			self.logDir = name
		else:
			self.logDir = os.path.join(self.outDir, name)

		self.access(self.logDir, create=True)

	def setTargetDir(self, targetDir : str, abs=True):
		'''targetDir must be an absolute reference. abs parameter exists only for consistency among sister methods.'''
		if self.access(targetDir, mode="r") is True:
			self.targetDir = targetDir

	def setRefDir(self, refDir : str, abs=False):
		if abs is True:
			if self.access(refDir, mode="rx", create=True) is True:
				self.refDir = refDir
		refOrder = []
		for d in [self.installDir, self.workDir, self.userDir]:
			p = os.path.join(d, refDir)
			if self.accessible(p, mode="rx", create=True) is True:
				refOrder.append(p)

		if refOrder == [] and os.path.exists(refDir):
			# refDir is determined to be an absolute path
			self.refDir = refDir
		elif refOrder == []:
			LOGGER.warning("Nowhere to read references from unless references are absolute paths to readable directories. No '{}' directory with reading and execution permission in installDir, workDir, or userDir:\n".format(refDir) + str(self))
			raise PermissionError("Nowhere to read references from unless references are absolute paths to readable directories. No '{}' directory with reading and execution permission in installDir, workDir, or userDir:\n".format(refDir) + str(self))
		else:
			self.refDir = refOrder[0]

	def setDatabaseDir(self, databaseDir : str, abs : bool=False):
		if abs is True:
			if self.access(databaseDir, mode="rx") is True:
				self.databaseDir = databaseDir
		databaseOrder = []
		for d in [self.installDir, self.workDir, self.userDir]:
			p = os.path.join(d, databaseDir)
			if self.accessible(p, mode="rx") is True:
				databaseOrder.append(p)

		if databaseOrder == [] and os.path.exists(databaseDir):
			# databaseDir is determined to be an absolute path
			self.databaseDir = databaseDir
		elif databaseOrder == []:
			LOGGER.warning("Nowhere to read databases from unless databases are absolute paths to readable directories. No '{}' directory with reading and execution permission in installDir, workDir, or userDir:\n".format(databaseDir) + str(self))
			raise PermissionError("Nowhere to read databases from unless databases are absolute paths to readable directories. No '{}' directory with reading and execution permission in installDir, workDir, or userDir:\n".format(databaseDir) + str(self))
		else:
			self.databaseDir = databaseOrder[0]
		

	def setTmpDir(self, tmpDir : str=None):
		'''If no place for a temporary directory is given, will find a nice place to put one, and then create a
		random-name temporary directory there. Temporary directory behavior is determined by the settings/flags.
		 * Note that tmpDir is not the name of the temporary directory, but the directory in which the temporary
		 directory will be placed.'''
		if tmpDir is None:
			tmpOrder = [self.installDir, self.userDir, self.workDir]
			tmpOrder = [dr for dr in tmpOrder if all(self.perms[dr].values())]
			
			if tmpOrder == []:
				LOGGER.error("Nowhere to put temporary files. Missing read+write+execute permission in any of installDir, workDir, or userDir:\n" + str(self))
				raise PermissionError("Missing write/read permissions in userDir, installDir, and workDir. Specify your workDir if you are intending to use a different workDir than the shell you are currently in.")
			
			#	Raise-protected
			self.tmpDir = self.newRandomName(tmpOrder[0])
		else:
			self.tmpDir = self.newRandomName(tmpDir)
		
		self.access(self.tmpDir, create=True)
	
	def setQuery(self, query : str, abs : bool=True):
		if abs is not True:
			self.query = self.get(query, "targetDir")
			self.access(self.query, mode="r")
		else:
			self.query = query
			self.access(self.query, mode="r")
	
	'''Get-functions'''

	def create(self, *paths): # Unsure if needed.
		os.makedirs( os.path.join(*paths))

	def get(self, filename : str, hint : str=None):
		'''Gets the absolute path to a file/directory found in one of the directories used by the Library. Hint can be
		provided to limit search to only one of the directories in the Library.
		If a path can not be found, it will be returned. This should be the case for a path which already is absolute.'''

		if filename in self.__cache:
			LOGGER.debug("Returned cached path '{path}' for filename '{filename}'".format(path=self.__cache[filename], filename=filename))
			return self.__cache[filename]

		if hint is not None:
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

		return filename
	
	def getReferences(self, references : list[str,str,str,str,str]):
		'''Returns list of references.'''
		references = {}
		for genome, strain, genbank_id, refseq_id, assembly_name in references:
			filename = DownloadReferences.download(genbank_id, refseq_id, assembly_name, dst=self.Lib.refDir)
			if not os.path.exists(filename):
				msg = "Could not download reference genome: {genbank_id='{genbank_id}', refseq_id='{refseq_id}' assembly_name='{assembly_name}'}".format(genbank_id=genbank_id, refseq_id=refseq_id, assembly_name=assembly_name)
				LOGGER.error(msg)
				raise FileNotFoundError(msg)
			references[genome] = filename
		return references
	
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