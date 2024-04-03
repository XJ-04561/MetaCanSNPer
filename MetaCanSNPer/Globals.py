

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
from PseudoPathy.PathShortHands import *
import PseudoPathy.Globals as PPGlobals
import VariantCallFixer.Globals as VCFGlobals
from VariantCallFixer import openVCF
from PseudoPathy import MinimalPathLibrary, PathLibrary, PathGroup, Path, DirectoryPath, FilePath, PathList
from PseudoPathy.Library import CommonGroups
from MetaCanSNPerDatabases import DatabaseReader, downloadDatabase, Branch, DatabaseWriter, openDatabase, IsLegacyCanSNPer2, updateFromLegacy
from MetaCanSNPerDatabases import Commands
from MetaCanSNPerDatabases import Columns as DB
import MetaCanSNPerDatabases as CanSNPDB
from MetaCanSNPerDatabases import Globals as DBGlobals
import random, logging, re, time, os, sys, shutil
from functools import cache, cached_property
from typing import Iterable, Callable, Any, Generator, Literal, AnyStr, TextIO, BinaryIO, Self
class Number: pass
Number = int|float
# class Comparisons:
#     def __init__(self, left, right=None, orAnd=None):
#         self.left, self.right, self.orAnd = left, right, orAnd
#     def __eq__(self, other):
#         if self.orAnd:
#             return other == self.left or other == self.right
#         else:
#             return other == self.left and other == self.right
#     def __or__(self, other):
#         return Comparisons(self, other, orAnd=True)
#     def __and__(self, other):
#         return Comparisons(self, other, orAnd=False)

# class Below(Comparisons):
#     def __eq__(self, other):
#         return other < self.left

# class Above(Comparisons):
#     def __eq__(self, other):
#         return other > self.left
    

random.seed()
from tempfile import NamedTemporaryFile


LOG_DIR = None
for root in CommonGroups().locals:
    path = os.path.join(root, f"{SOFTWARE_NAME}-Results", "Logs")
    if pAccess(path, "rwx"):
        LOG_DIR = path
        break
if LOG_DIR is None:
    for root in CommonGroups().locals:
        path = os.path.join(root, f"{SOFTWARE_NAME}-Results", "Logs")
        if pBackAccess(path, "w"):
            if pMakeDirs(path):
                LOG_DIR = path
                break
del root, path
## LogKeeper Globals
with NamedTemporaryFile(prefix=time.strftime("MetaCanSNPer-(%Y-%m-%d)-(%H-%M-%S)-", time.localtime()), suffix=".log", dir=LOG_DIR, delete=False) as f:
    LOGGING_FILEPATH = f.name
LOGGER_FILEHANDLER = logging.FileHandler(LOGGING_FILEPATH)
LOGGER_FILEHANDLER.setFormatter(logging.Formatter("[%(name)s] %(asctime)s - %(levelname)s: %(message)s"))

PPGlobals.LOGGER.addHandler(LOGGER_FILEHANDLER)
VCFGlobals.LOGGER.addHandler(LOGGER_FILEHANDLER)
DBGlobals.LOGGER.addHandler(LOGGER_FILEHANDLER)

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

## DownloadReferences Globals

DOWNLOAD_TIMEOUT = 300
SOURCED = {"refseq":"F", "genbank": "A"}
"""```{"refseq":"F", "genbank": "A"}```"""
NCBI_FTP_LINK = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GC{source}/{n1}/{n2}/{n3}/{genome_id}_{assembly}/{genome_id}_{assembly}_genomic.fna.gz"
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
