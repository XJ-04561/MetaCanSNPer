

import os
import shutil
import random
random.seed()
from logging import Logger

'''Container and handler of directories and files'''
class DirectoryLibrary:
	workDir : str
	installDir : str
	userDir : str

	targetDir : str
	refDir : str
	tmpDir : str
	logDir : str
	outDir : str
	dirs = ["targetDir", "refDir", "tmpDir"]

	tmps : list[str]
	references : dict[list[str]]
	targets : list[str]

	__cache : dict[str]
	
	def __init__(self, settings : dict, workDir : str=None, installDir : str=None, userDir : str=None, logger : Logger=Logger,
			     tmpDir : str=None, refDir : str="References", outDir : str=None, **kwargs):
		'''Needs no arguments to initialize, but keyword arguments can be used to force use directories. Keyword
		arguments not named in the method head are caught and ignored by **kwargs.
		
		workDir		-	Will be used to find queries/targets, and is the second directory to be searched for references
		
		installDir	-	Will be used to store temporarily generated outputs, data, and converted data, and is the first
						place to look for references.

		userDir		-	Will be used as a backup place for temporary storage if permissions don't allow the installDir,
						and is the third place to look for references.
		'''

		self.settings = settings
		self.logger = logger

		if workDir is None:
			self.workDir = os.getcwd()
		elif not os.path.exists(workDir):
			self.workDir = workDir
			logger.info("Creating working directory {wd}".format(wd=self.workDir))
			os.makedirs(self.workDir)
			os.chdir(self.workDir)
		else:
			os.chdir(self.workDir)
		if installDir is None:
			self.installDir, _ = __file__.rsplit(os.path.sep)
		if userDir is None:
			self.userDir = os.path.expanduser("~")

		self.logger.debug("Work dir:    '{dr}'".format(self.workDir))
		self.logger.debug("Install dir: '{dr}'".format(self.installDir))
		self.logger.debug("User dir:    '{dr}'".format(self.userDir))

		self.perms = {}

		for d in [self.workDir, self.installDir, self.userDir]:
			self.__permCheck__(self, d)
		
		self.logger.debug("Perms for directories:\n" + "\n".join(["{}={}{}{}".format(path,*[c if perms[c] else " " for c in "rwx"]) for path, perms in self.perms]))

		self.setTargetDir(self.workDir)
		self.setRefDir(refDir)
		self.createNewTemporary(tmpDir)
		self.setOutDir(outDir)
		
		self.__cache = {}
	
	def __getitem__(self, key):
		return self.__getattribute__(key)

	def __str__(self):
		return "\n".join(["{:<20}='{}{}{}' '{}'".format(d,*[c if self.perms[self[d]][c] else " " for c in "rwx"], self[d]) for d in self.dirs])

	def __del__(self):
		if os.path.exists(self.tmpDir):
			shutil.rmtree(self.tmpDir)
	
	def __permCheck__(self, dr):
		self.perms[dr]["r"] = os.access(dr, os.R_OK)
		self.perms[dr]["w"] = os.access(dr, os.W_OK)
		self.perms[dr]["x"] = os.access(dr, os.X_OK)
	
	@staticmethod
	def newRandomName(path : str, ext : str=None):
		newName = "tmp_{}{}{}".format(*random.random().as_integer_ratio(), "."+ext if ext is not None else "")
		n=0
		while os.path.exists(os.path.sep.join([path, newName])):
			newName = "tmp_{}{}{}".format(*random.random().as_integer_ratio(), "."+ext if ext is not None else "")
			if n>1000000-1:
				raise FileExistsError("Could not generate random file/directory name for path='{path}' after {n} attempts".format(path=path, n=n))
				break
		return os.path.sep.join([path, newName])

	def setOutDir(self, outDir : str):
		if outDir is None:
			if not all(self.perms[self.workDir].values()):
				self.logger.error("No valid working directory or output directory given.\n" + str(self))
				raise PermissionError("No valid working directory or output directory given.\n" + str(self))
			else:
				self.outDir = os.path.sep.join([self.workDir, "CanSNPer"])
				os.mkdir(self.outDir)
		else:
			self.__permCheck__(outDir)
			if all(self.perms[outDir].values()):
				self.outDir = outDir
			else:
				raise PermissionError("Given outDir='{}' does not have read, write, and execute permissions. It has '{}{}{}'".format(outDir, *[c if self.perms[outDir][c] else " " for c in "rwx"]))
		self.setResultDir()
		self.setLogDir()
		
	def setResultDir(self, name : str="Results"):
		if self.settings["FileManagement"]["MergeNewResults"]:
			self.resultDir = os.path.sep.join([self.outDir, name])
			if not os.path.exists(self.resultDir):
				os.mkdir(self.resultDir)
		else:
			n=1
			while os.path.exists("{}-{}".format(name, str(n))):
				n+1
			self.resultDir = os.path.sep.join([self.outDir, "{}-{}".format(name, str(n))])
			os.mkdir(self.resultDir)
	
	def setLogDir(self, name : str="Logs"):
		
		self.logDir = os.path.sep.join([self.outDir, name])
		if os.path.exists(name):
			os.mkdir(self.logDir)

	def setTargetDir(self, targetDir : str):
		if not all(self.perms[targetDir].values()):
			targetDir = None
			self.logger.warning("No read+write permission in working directory, targets/queries must be absolute paths.: '{wd}'".format(wd=targetDir))
		else:
			self.targetDir = targetDir

	def setRefDir(self, refDir : str):
		for d in [self.workDir, self.installDir, self.userDir]:
			if os.path.exists(os.path.sep.join([d, refDir])):
				self.__permCheck__(self, os.path.sep.join([d, refDir]))
		refOrder = [os.path.sep.join([d, refDir]) for d in [self.installDir, self.workDir, self.userDir]]
		refOrder = [dr for dr in self.refOrder if self.perms[dr]["r"] and self.perms[dr]["x"]]
		if refOrder == []:
			self.logger.warning("Nowhere to read references from unless references are absolute paths to readable directories. No '{}' directory with reading and execution permission in installDir, workDir, or userDir:\n".format(refDir) + str(self))
			raise PermissionError("Nowhere to read references from unless references are absolute paths to readable directories. No '{}' directory with reading and execution permission in installDir, workDir, or userDir:\n".format(refDir) + str(self))
		else:
			self.refDir = refOrder[0]

	def createNewTemporary(self, tmpDir : str=None):
		if tmpDir is None:
			tmpOrder = [self.installDir, self.userDir, self.workDir]
			tmpOrder = [dr for dr in tmpOrder if all(self.perms[dr].values())]
			
			if tmpOrder == []:
				self.logger.error("Nowhere to put temporary files. Missing read+write permission in any of installDir, workDir, or userDir:\n" + str(self))
				raise PermissionError("Missing write/read permissions in userDir, installDir, and workDir. Specify your workDir if you are intending to use a different workDir than the shell you are currently in.")
			
			#	Raise-protected
			self.tmpDir = self.newRandomName(tmpOrder[0])
		else:
			self.tmpDir = self.newRandomName(tmpDir)
		
		self.__permCheck__(self.tmpDir)
		if not all(self.perms[self.tmpDir].values()):
			msg = "Missing crucial permissions in created temporary folder '{dr}' missing {perms} permission.".format(
				dr=self.tmpDir,
				perms=" and ".join([key for key, value in self.perms[self.tmpDir].items() if value]))
			self.logger.error(msg)
			raise PermissionError(msg)

		self.logger.info("Creating new temporary folder '{newDir}'".format(newDir=self.tmpDir))
		try:
			os.mkdir(self.tmpDir)
		except Exception as e:
			self.logger.error("Failed to create temporar directory: '{}'".format(self.tmpDir))
			raise e
		
		return self.tmpDir
	
	def get(self, filename : str, hint : str=None):
		'''Gets the absolute path to a file/directory found in one of the directories used by the Library. Hint can be
		provided to limit search to only one of the directories in the Library.
		If a path can not be found, it will be returned. This should be the case for a path which already is absolute.'''

		if filename in self.__cache:
			self.logger.debug("Returned cached path '{path}' for filename '{filename}'".format(path=self.__cache[filename], filename=filename))
			return self.__cache[filename]

		if hint is not None:
			# directory hint has been provided
			for path, dirs, files in os.walk(self[hint]):
				for f in files:
					if f == filename:
						self.__cache[filename] = os.path.sep.join([path, f])
						return os.path.sep.join([path, f])
					
			self.logger.info("Looked for file: '{}' in {} but could not find it.\n{}".format(filename, hint, str(self)))
		else:
			# Without a hint, each directory is searched.
			for dr in {self[drVar] for drVar in self.dirs}:
				for path, dirs, files in os.walk(self[dr]):
					for f in files:
						if f == filename:
							return os.path.sep.join([path, f])
		
			self.logger.info("Looked for file: '{}' in all directories but could not find it.\n{}".format(filename, str(self)))

		return filename