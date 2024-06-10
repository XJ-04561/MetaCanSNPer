

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

SOFTWARE_NAME = "MetaCanSNPer"
DRY_RUN = False
RUNNING = True
DIRECTORY_OPTIONS = ["workDir", "userDir", "installDir", "targetDir", "tmpDir", "refDir", "databaseDir", "outDir", "sessionName"]
DEFAULT_SETTINGS : "Config"

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

DEFAULT_TOML_TEMPLATE = """[FileManagement]
saveTemp = false
referenceFormats = [".fna", ".fasta"]

[Software]
mapper = "minimap2"
aligner = "progressiveMauve"
snpCaller = "gatk_Mutect2"

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

from PseudoPathy.PathShortHands import *
import PseudoPathy.Globals as PPGlobals
import VariantCallFixer.Globals as VCFGlobals
import SQLOOP.Globals as SQLOOPGlobals
from VariantCallFixer import openVCF
from PseudoPathy import *
from PseudoPathy.Paths import FileList
from collections import namedtuple

from types import FunctionType, MethodType
import random, logging, re, time, os, sys, shutil, itertools, logging.handlers as logHandlers
from functools import cache, cached_property
from typing import Iterable, Callable, Any, Generator, Literal, AnyStr, TextIO, BinaryIO, overload
from time import sleep
from collections import defaultdict, OrderedDict

from GeekyGadgets import (
	DEV_NULL, DEV_NULL_BYTES, ISATTY, Logged, Default, ClassProperty, CachedClassProperty, Config, loadTOML,
	forceHash, getAttrChain, Hooks, Hook, DummyHooks, callFormat, timeFormat, LimitedDict, MissingDependency,
	this
)
import GeekyGadgets.MonkeyPatch

from appdirs import user_log_dir, user_config_dir, site_config_dir
from threading import Lock
random.seed()

SQLOOPGlobals.MAX_DEBUG = MAX_DEBUG = False

LOGGING_FILEPATH = UniqueFilePath(DirectoryPath(user_log_dir(SOFTWARE_NAME)).writable, time.strftime("MetaCanSNPer-%Y-%m-%d--%H-%M-%S.log", time.localtime()))
LOGGING_FILEHANDLER = logging.FileHandler(LOGGING_FILEPATH)
logging.basicConfig(handlers=[LOGGING_FILEHANDLER], format="[%(name)s] %(asctime)s - %(levelname)s: %(message)s", level=logging.DEBUG)

Logged.LOG = LOGGER = logging.Logger("MetaCanSNPer", level=logging.DEBUG)

## Default .toml



_configDir = Path(user_config_dir(SOFTWARE_NAME))
_configDir.create(purpose="rw")
if "defaults.toml" not in _configDir:
	with open(_configDir / "defaults.toml", "w") as f:
		f.write(DEFAULT_TOML_TEMPLATE)
DEFAULT_SETTINGS = loadTOML(_configDir / "defaults.toml")