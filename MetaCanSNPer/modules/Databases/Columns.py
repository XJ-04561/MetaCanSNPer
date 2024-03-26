
from MetaCanSNPer.modules.Databases._Constants import *


class ColumnFlag(int): pass

ALL				= ColumnFlag(0)
nodeID			= ColumnFlag(1)
snpID			= ColumnFlag(2)
genomeID		= ColumnFlag(3)
position		= ColumnFlag(4)
ancestral		= ColumnFlag(5)
derived			= ColumnFlag(6)
snpReference	= ColumnFlag(7)
date			= ColumnFlag(8)
genome			= ColumnFlag(9)
strain			= ColumnFlag(10)
genbankID		= ColumnFlag(11)
refseqID		= ColumnFlag(12)
assembly		= ColumnFlag(13)
treeParent		= ColumnFlag(14)
treeChild		= ColumnFlag(15)
treeRank		= ColumnFlag(16)
# Experimental
genoType		= ColumnFlag(17)
chromosome		= ColumnFlag(18)

NAMES = [
    "*",
	nodeID,
	snpID,
	genomeID,
	position,
	ancestral,
	derived,
	snpReference,
	date,
	genome,
	strain,
	genbankID,
	refseqID,
	assembly,
	treeParent,
	treeChild,
	treeRank,
	genoType,
	chromosome
]

LOOKUP = {
    TABLE_NAME_NODES : {
        ALL			: "*",
        nodeID		: NODE_COLUMN_ID,
        genoType	: NODE_COLUMN_NAME
	},
    TABLE_NAME_REFERENCES : {
        ALL			: "*",
        genomeID	: REFERENCE_COLUMN_GENOME_ID,
		genome		: REFERENCE_COLUMN_GENOME,
		strain		: REFERENCE_COLUMN_STRAIN,
		genbankID	: REFERENCE_COLUMN_GENBANK,
		refseqID	: REFERENCE_COLUMN_REFSEQ,
		assembly	: REFERENCE_COLUMN_ASSEMBLY
	},
    TABLE_NAME_SNP_ANNOTATION : {
        ALL				: "*",
        nodeID			: SNP_COLUMN_NODE_ID,
		snpID			: SNP_COLUMN_SNP_ID,
		position		: SNP_COLUMN_POSITION,
		ancestral		: SNP_COLUMN_ANCESTRAL,
		derived			: SNP_COLUMN_DERIVED,
		snpReference	: SNP_COLUMN_REFERENCE,
		date			: SNP_COLUMN_DATE,
		genomeID		: SNP_COLUMN_GENOME_ID
	},
    TABLE_NAME_TREE : {
        ALL : "*",
        treeParent	: TREE_COLUMN_PARENT,
		treeChild	: TREE_COLUMN_CHILD,
		treeRank	: TREE_COLUMN_RANK
	}
}

RELATIONSHIPS = {
	TABLE_NAME_NODES : {
		snpID			: (TABLE_NAME_SNP_ANNOTATION, nodeID),
		position		: (TABLE_NAME_SNP_ANNOTATION, nodeID),
		ancestral		: (TABLE_NAME_SNP_ANNOTATION, nodeID),
		derived			: (TABLE_NAME_SNP_ANNOTATION, nodeID),
		snpReference	: (TABLE_NAME_SNP_ANNOTATION, nodeID),
		date			: (TABLE_NAME_SNP_ANNOTATION, nodeID),
		genomeID		: (TABLE_NAME_SNP_ANNOTATION, nodeID),
        
		treeParent	: (TABLE_NAME_TREE, nodeID),
		treeChild	: (TABLE_NAME_TREE, nodeID),
		treeRank	: (TABLE_NAME_TREE, nodeID)
	},
    TABLE_NAME_REFERENCES : {
        nodeID			: (TABLE_NAME_SNP_ANNOTATION, genomeID),
		snpID			: (TABLE_NAME_SNP_ANNOTATION, genomeID),
		position		: (TABLE_NAME_SNP_ANNOTATION, genomeID),
		ancestral		: (TABLE_NAME_SNP_ANNOTATION, genomeID),
		derived			: (TABLE_NAME_SNP_ANNOTATION, genomeID),
		snpReference	: (TABLE_NAME_SNP_ANNOTATION, genomeID),
		date			: (TABLE_NAME_SNP_ANNOTATION, genomeID)
	},
    TABLE_NAME_SNP_ANNOTATION : {
        genome		: (TABLE_NAME_REFERENCES, genomeID),
		strain		: (TABLE_NAME_REFERENCES, genomeID),
		genbankID	: (TABLE_NAME_REFERENCES, genomeID),
		refseqID	: (TABLE_NAME_REFERENCES, genomeID),
		assembly	: (TABLE_NAME_REFERENCES, genomeID),
        
		genoType	: (TABLE_NAME_NODES, nodeID),
        
		treeParent	: (TABLE_NAME_TREE, nodeID),
		treeChild	: (TABLE_NAME_TREE, nodeID)
	}
}