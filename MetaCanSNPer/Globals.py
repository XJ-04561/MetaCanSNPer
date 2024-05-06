

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
from PseudoPathy import MinimalPathLibrary, PathLibrary, PathGroup, Path, DirectoryPath, FilePath, PathList
from PseudoPathy.Library import CommonGroups
from collections import namedtuple

from types import FunctionType, MethodType
import random, logging, re, time, os, sys, shutil, itertools
from functools import cache, cached_property
from typing import Iterable, Callable, Any, Generator, Literal, AnyStr, TextIO, BinaryIO, Self
from time import sleep

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

class _IterTimerInit:
    n : int|None
    def __init__(self, n=None):
        self.n = n
    def __iter__(self):
        return iter(time.localtime()[:self.n])

def replTimeMatch(timer, m):
    return format(timer[m.group(0)[0]], f">0{len(m.group(0))}")

timeLetterPattern = re.compile(r"(Y+|M+|D+|h+|m+|s+)")

class Time(namedtuple("Time", ["Y", "M", "D", "h", "m", "s"], defaults=_IterTimerInit(6))):
    def __format__(self, fs : str):
        timeLetterPattern.sub(fs, replTimeMatch.__get__(self, type(self)))

LOG_DIR = None

## LogKeeper Globals
with NamedTemporaryFile(prefix=time.strftime("MetaCanSNPer-(%Y-%m-%d)-(%H-%M-%S)-", time.localtime()), suffix=".log", dir=LOG_DIR, delete=False) as f:
    LOGGING_FILEPATH = f.name
LOGGER_FILEHANDLER = logging.FileHandler(LOGGING_FILEPATH)
LOGGER_FILEHANDLER.setFormatter(logging.Formatter("[%(name)s] %(asctime)s - %(levelname)s: %(message)s"))

PPGlobals.LOGGER.addHandler(LOGGER_FILEHANDLER)
VCFGlobals.LOGGER.addHandler(LOGGER_FILEHANDLER)

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
