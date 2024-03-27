
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
	ALL,
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

LOOKUP : dict[str,dict[ColumnFlag, str]] = {
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
        ALL			: "*",
        treeParent	: TREE_COLUMN_PARENT,
		treeChild	: TREE_COLUMN_CHILD,
		treeRank	: TREE_COLUMN_RANK
	}
}

UNIQUES = set(LOOKUP[TABLE_NAME_NODES]).difference(LOOKUP[TABLE_NAME_REFERENCES], LOOKUP[TABLE_NAME_SNP_ANNOTATION], LOOKUP[TABLE_NAME_TREE])
UNIQUELOOKUP = {
    col : [table for table in LOOKUP if col in LOOKUP[table]][0] for col in UNIQUES
}
COMMONS = set(LOOKUP[TABLE_NAME_NODES]).intersection(LOOKUP[TABLE_NAME_REFERENCES], LOOKUP[TABLE_NAME_SNP_ANNOTATION], LOOKUP[TABLE_NAME_TREE])
COMMONLOOKUP = {
    col : {table for table in LOOKUP if col in LOOKUP[table]} for col in COMMONS
}

RELATIONS : dict[tuple[str,str],ColumnFlag] = {
    (TABLE_NAME_NODES, TABLE_NAME_SNP_ANNOTATION) : nodeID,
    (TABLE_NAME_SNP_ANNOTATION, TABLE_NAME_NODES) : nodeID,
    
	(TABLE_NAME_NODES, TABLE_NAME_TREE) : nodeID,
    (TABLE_NAME_TREE, TABLE_NAME_NODES) : nodeID,
    
	(TABLE_NAME_REFERENCES, TABLE_NAME_SNP_ANNOTATION) : genomeID,
    (TABLE_NAME_SNP_ANNOTATION, TABLE_NAME_REFERENCES) : genomeID,
    
	(TABLE_NAME_SNP_ANNOTATION, TABLE_NAME_REFERENCES) : genomeID,
    (TABLE_NAME_REFERENCES, TABLE_NAME_SNP_ANNOTATION) : genomeID,
    
	(TABLE_NAME_SNP_ANNOTATION, TABLE_NAME_NODES) : nodeID,
    (TABLE_NAME_NODES, TABLE_NAME_SNP_ANNOTATION) : nodeID,
    
	(TABLE_NAME_SNP_ANNOTATION, TABLE_NAME_TREE) : nodeID,
    (TABLE_NAME_TREE, TABLE_NAME_SNP_ANNOTATION) : nodeID
}

RELATIONSHIPS : dict[str,dict[ColumnFlag, str]]= {
	TABLE_NAME_NODES : {
		snpID			: TABLE_NAME_SNP_ANNOTATION,
		position		: TABLE_NAME_SNP_ANNOTATION,
		ancestral		: TABLE_NAME_SNP_ANNOTATION,
		derived			: TABLE_NAME_SNP_ANNOTATION,
		snpReference	: TABLE_NAME_SNP_ANNOTATION,
		date			: TABLE_NAME_SNP_ANNOTATION,
		genomeID		: TABLE_NAME_SNP_ANNOTATION,
        
		treeParent	: TABLE_NAME_TREE,
		treeChild	: TABLE_NAME_TREE,
		treeRank	: TABLE_NAME_TREE
	},
    TABLE_NAME_REFERENCES : {
        nodeID			: TABLE_NAME_SNP_ANNOTATION,
		snpID			: TABLE_NAME_SNP_ANNOTATION,
		position		: TABLE_NAME_SNP_ANNOTATION,
		ancestral		: TABLE_NAME_SNP_ANNOTATION,
		derived			: TABLE_NAME_SNP_ANNOTATION,
		snpReference	: TABLE_NAME_SNP_ANNOTATION,
		date			: TABLE_NAME_SNP_ANNOTATION
	},
    TABLE_NAME_SNP_ANNOTATION : {
        genome		: TABLE_NAME_REFERENCES,
		strain		: TABLE_NAME_REFERENCES,
		genbankID	: TABLE_NAME_REFERENCES,
		refseqID	: TABLE_NAME_REFERENCES,
		assembly	: TABLE_NAME_REFERENCES,
        
		genoType	: TABLE_NAME_NODES,
        
		treeParent	: TABLE_NAME_TREE,
		treeChild	: TABLE_NAME_TREE
	}
}