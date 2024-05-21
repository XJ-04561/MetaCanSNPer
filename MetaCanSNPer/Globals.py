

"""
MetaCanSNPer
"""

oname = __name__
__name__ 		= "MetaCanSNPer"
__version__ 	= "0.1.0"
__author__ 		= "Fredrik Sörensen"
__credits__ 	= ["Fredrik Sörensen", "David Sundell"]
__license__ 	= "GPLv3"
__maintainer__ 	= "FOI bioinformatics group"
__email__ 		= ["bioinformatics@foi.se", "fredrik.sorensen@foi.se"]
__date__ 		= "2024-02-27"
__status__ 		= "Prototype"

## Common Globals

SOFTWARE_NAME = "MetaCanSNPer"
DRY_RUN = False
RUNNING = True
from This import this
from PseudoPathy.PathShortHands import *
import PseudoPathy.Globals as PPGlobals
import VariantCallFixer.Globals as VCFGlobals
from VariantCallFixer import openVCF
from PseudoPathy import *
from collections import namedtuple

from types import FunctionType, MethodType
import random, logging, re, time, os, sys, shutil, itertools
from functools import cache, cached_property
from typing import Iterable, Callable, Any, Generator, Literal, AnyStr, TextIO, BinaryIO, Self, overload
from time import sleep
from collections import defaultdict, OrderedDict
import tomllib as toml
from appdirs import user_log_dir, user_config_dir, site_config_dir
from tempfile import NamedTemporaryFile

LOG_DIR = user_log_dir(SOFTWARE_NAME)
pMakeDirs(LOG_DIR)
with NamedTemporaryFile(prefix=time.strftime("MetaCanSNPer-(%Y-%m-%d)-(%H-%M-%S)-[", time.localtime()), suffix="].log", dir=LOG_DIR, delete=False) as f:
	LOGGING_FILEPATH = f.name
LOGGER_FILEHANDLER = logging.FileHandler(LOGGING_FILEPATH)
LOGGER_FILEHANDLER.setFormatter(logging.Formatter("[%(name)s] %(asctime)s - %(levelname)s: %(message)s"))
LOGGER = logging.Logger("MetaCanSNPer")
LOGGER.addHandler(LOGGER_FILEHANDLER)

DEV_NULL = open(os.devnull, "w")
ISATTY = sys.stdout.isatty()

## Default .toml

DEFAULT_TOML_TEMPLATE = """[FileManagement]
saveTemp = false
referenceFormats = [".fna", ".fasta"]

[Software]
# mapper = minimap2
# aligner = progressiveMauve
# snpCaller = gatk_Mutect2

[Directories]
# workDir = 
# userDir = 
# installDir = 
# targetDir = 
# refDir = 
# databaseDir = 
# tmpDir = 
# outDir = 
# sessionName = 
"""

def loadFlattenedTOML(filename):
	with open(filename, "rb") as f:
		tmp : dict[str,dict] = toml.load(f)
	# Settings hierarchy looks like this: ["Category"]["Flag"] -> Value
	# Flatten hierarchy so flags are easily accessible
	settings = {}
	for flags in tmp.values():
		for flag, value in flags.items():
			settings[flag] = value
	return settings
_configDir = Path(user_config_dir(SOFTWARE_NAME))
_configDir.create(purpose="rw")
if "defaults.toml" not in _configDir:
	with open(_configDir / "defaults.toml", "w") as f:
		f.write(DEFAULT_TOML_TEMPLATE)
DEFAULT_SETTINGS = loadFlattenedTOML(_configDir / "defaults.toml")


PPGlobals.LOGGER = LOGGER.getChild("PseudoPathy")
VCFGlobals.LOGGER = LOGGER.getChild("VariantCallFixer")

class Number: pass
Number = int|float

class NoneStr:
	def __mult__(self, other): return self
	def __add__(self, other): return self
	def __getattribute__(self, name): return self
	def __str__(self): return self

def printCall(func, args, kwargs):
	return f"{getattr(func, '__qualname__', getattr(func, '__name__', func))}({', '.join(itertools.chain(map(str, args), map(lambda keyval : str(keyval[0])+"="+str(keyval[1]), kwargs.items())))})"

random.seed()
from tempfile import NamedTemporaryFile

_NULL_KEY = object()

class UninitializedError(AttributeError):
	def __init__(self, obj=None, name=None, **kwargs):
		if obj is not None:
			objName = repr(type(obj).__name__)
		else:
			objName = "Object"
		if name is not None:
			name = repr(name) + " "
		else:
			name = ""
		super().__init__(f"Attribute {name}of {objName} was accessed, but has yet to be set.", **kwargs)

class InitCheckDescriptor:

	def __init__(cls, className, bases, namespace):
		cls.names = {}

	def __set_name__(self, owner, name):
		object.__getattribute__(self, "names")[id(owner)] = name

	def __set__(self, instance, value):
		instance.__dict__[object.__getattribute__(self, "names")[id(type(instance))]] = value

	def __get__(self, instance, owner=None):
		return instance.__dict__.get(object.__getattribute__(self, "names").get(id(owner), id(_NULL_KEY)), self)

	def __delete__(self, instance):
		pass

	def __mult__(self, other): return self
	def __add__(self, other): return self
	def __sub__(self, other): return self
	def __getattribute__(self, name): return self

	def __bool__(self): return False

class NotSet(metaclass=InitCheckDescriptor): pass


## ArgParser Globals

DIRECTORY_OPTIONS = ["workDir", "userDir", "installDir", "targetDir", "tmpDir", "refDir", "databaseDir", "outDir", "sessionName"]

MAPPER_OPTIONS_EXPLAINER = """
To provide flags/arguments for the chosen Mapper, provide them
directly after the '--indexerOptions' flag, only interrupted by the end of the
command call or the corresponding flag for an Aligner or SNP Caller options.
"""
ALIGNER_OPTIONS_EXPLAINER = """
To provide flags/arguments for the chosen Aligner, provide them
directly after the '--indexerOptions' flag, only interrupted by the end of the
command call or the corresponding flag for a Mapper or SNP Caller options.
"""
SNP_CALLER_OPTIONS_EXPLAINER = """
To provide flags/arguments for the chosen SNP Caller, provide them directly
after the '--snpCallerOptions' flag, only interrupted by the end of the
command call or the corresponding flag for Mapper or Aligner options.
"""

class MissingDependancy(Exception): pass

## ErrorFixes Globals
class Aligner: pass
class Mapper: pass
class SNPCaller: pass


class Default:
	
	SIZE_LIMIT = 10000

	deps : tuple = None
	data : dict
	name : str

	fget : Callable
	fset : Callable
	fdel : Callable

	def __init__(self, fget=None, fset=None, fdel=None, doc=None, deps=None):
		self.data = {}
		if hasattr(fget, "__annotations__"):
			self.__annotations__ = fget.__annotations__
		self.fget = fget
		self.fset = fset
		self.fdel = fdel
		self.__doc__ = doc or fget.__doc__ or getattr(self, "__doc__", None)
		if deps is not None:
			self.deps = deps
		
	def __call__(self, fget=None, fset=None, fdel=None, doc=None):
		self.__init__(fget, fset, fdel, doc=doc)
		return self
	
	def __class_getitem__(cls, deps):
		"""Calls Default(None) and adds the keys provided as the names of attributes upon which this value depends
		before returning. This is useful for creating attributes which have default values which are meant to be
		dependent on other attributes of the same object. When getting the same attribute repeatedly, new attribute
		value instances will not be created, the first one is returned until one of the dependency attributes are
		changed."""
		return cls(deps=deps)
	
	def __set_name__(self, owner, name):
		self.name = name

	def __get__(self, instance, owner=None):
		if instance is None:
			return self
		if self.name in instance.__dict__:
			return instance.__dict__[self.name]
		
		idSpec = None if self.deps is None else tuple(id(getattr(instance, attrName, None)) for attrName in self.deps)
		if id(instance) not in self.data:
			self.data[id(instance)] = (idSpec, self.fget(instance))
			return self.data[id(instance)][1]
		elif self.data[id(instance)][0] != idSpec:
			self.data[id(instance)] = (idSpec, self.fget(instance))
			return self.data[id(instance)][1]
		else:
			return self.data[id(instance)][1]
	
	def __set__(self, instance, value):
		if self.fset is None:
			instance.__dict__[self.name] = value
		else:
			self.fset(instance, value)
	
	def __delete__(self, instance, owner=None):
			if self.fdel is not None:
				self.fdel(instance)
			else:
				instance.__dict__.pop(self.name, None)
				if id(instance) in self.data:
					del self.data[id(instance)]
	
	def setter(self, fset):
		self.fset = fset
		return self
	
	def deleter(self, fdel):
		self.fdel = fdel
		return self