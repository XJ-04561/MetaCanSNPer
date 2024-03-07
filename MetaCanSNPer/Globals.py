

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
from PseudoPathy.PathShortHands import *
import random
random.seed()

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

## Database Globals

TABLE_NAME_SNP_ANNOTATION = "snp_annotation"
TABLE_NAME_REFERENCES = "snp_references"
TABLE_NAME_NODES = "nodes"
TABLE_NAME_TREE = "tree"

SNP_COLUMN_NODE_ID      = "node_id"
SNP_COLUMN_SNP_ID       = "snp_id"
SNP_COLUMN_POSITION     = "position"
SNP_COLUMN_ANCESTRAL    = "ancestral_base"
SNP_COLUMN_DERIVED      = "derived_base"
SNP_COLUMN_REFERENCE    = "reference"
SNP_COLUMN_DATE         = "date"
SNP_COLUMN_GENOME_ID    = "genome_i"

REFERENCE_COLUMN_GENOME_ID    = "id"
REFERENCE_COLUMN_GENOME       = "genome"
REFERENCE_COLUMN_STRAIN       = "strain"
REFERENCE_COLUMN_GENBANK      = "genbank_id"
REFERENCE_COLUMN_REFSEQ       = "refseq_id"
REFERENCE_COLUMN_ASSEMBLY     = "assembly_name"

NODE_COLUMN_ID                  = "id"
NODE_COLUMN_NAME                = "name"

TREE_COLUMN_PARENT              = "parent"
TREE_COLUMN_CHILD               = "child"
TREE_COLUMN_RANK                = "rank_i"

## DownloadReferences Globals

SOURCED = {"refseq":"F", "genbank": "A"}
NCBI_FTP_LINK = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GC{source}/{n1}/{n2}/{n3}/{genome_id}_{assembly}/{genome_id}_{assembly}_genomic.fna.gz"
