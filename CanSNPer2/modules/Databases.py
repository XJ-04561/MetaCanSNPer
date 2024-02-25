
from functools import cached_property
import sqlite3

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

class DatabaseReader:
    _connection : sqlite3.Connection

    references : list[tuple[int,str,str,str,str]]
    SNPs : list[tuple[str,int,str,str]]

    def __init__(self, database : str):
        self._connection = sqlite3.connect("file:{}?mode=ro".format(database), uri=True)

    def __del__(self):
        try:
            self._connection.close()
        except:
            pass
    
    def commit(self):
        self._connection.commit()

    def genomeID(self, genome) -> int:
        return self._connection.execute(
            "SELECT ? FROM ? WHERE ? = ?;",
            [
                REFERENCE_COLUMN_GENOME_ID,
                TABLE_NAME_REFERENCES,
                REFERENCE_COLUMN_GENOME,
                genome
            ]).fetchone()

    @cached_property
    def references(self) -> list[tuple[int,str,str,str,str]]:
        ''''''
        return self._connection.execute(
            "SELECT ?, ?, ?, ?, ? FROM ? ORDER BY ? ASC;",
            [
                REFERENCE_COLUMN_GENOME_ID,
                REFERENCE_COLUMN_GENOME,
                REFERENCE_COLUMN_GENBANK,
                REFERENCE_COLUMN_REFSEQ,
                REFERENCE_COLUMN_ASSEMBLY,
                TABLE_NAME_REFERENCES,
                REFERENCE_COLUMN_GENOME_ID
            ]).fetchall()

    def _getSNPsByGenomeId(self, genomeID : int) -> sqlite3.Cursor:
        return self._connection.execute(
            "SELECT ?, ?, ?, ? FROM ? WHERE ? = ? ORDER BY ? ASC;",
            [
                SNP_COLUMN_SNP_ID,
                SNP_COLUMN_POSITION,
                SNP_COLUMN_ANCESTRAL,
                SNP_COLUMN_DERIVED,
                TABLE_NAME_SNP_ANNOTATION,
                SNP_COLUMN_GENOME_ID,
                genomeID,
                SNP_COLUMN_POSITION
            ])
    
    @cached_property
    def SNPs(self) -> sqlite3.Cursor:
        return self._connection.execute(
            "SELECT ?, ?, ?, ? FROM ? ORDER BY ? ASC;",
            [
                SNP_COLUMN_SNP_ID,
                SNP_COLUMN_POSITION,
                SNP_COLUMN_ANCESTRAL,
                SNP_COLUMN_DERIVED,
                TABLE_NAME_SNP_ANNOTATION,
                SNP_COLUMN_POSITION
            ])

    @cached_property
    def SNPsByGenome(self) -> dict[str,list[tuple[str,int,str,str]]]:
        return {genome:self._getSNPsByGenomeId(genome_id) for genome_id, genome, _, _, _ in self.references}

    @cached_property
    def SNPsByID(self) -> dict[str,list[tuple[str,int,str,str]]]:
        return {snpID:(pos, anc, der) for snpID, pos, anc, der in self.SNPs}
    
    @cached_property
    def nodes(self) -> dict[str,list[tuple[str,int,str,str]]]:
        nodes = self._connection.execute(
            "SELECT ?, ? FROM ?;",
            [
                NODE_COLUMN_ID,
                NODE_COLUMN_NAME,
                TABLE_NAME_NODES
            ])
        return dict(nodes)

    @cached_property
    def tree(self) -> dict[int,list[int]]:
        """{nodeID:[child1, child2, ...]}"""
        parameters = [
            TREE_COLUMN_CHILD,
            TABLE_NAME_TREE,
            TREE_COLUMN_PARENT
        ]
        ret = {node:self._connection.execute("SELECT ? FROM ? WHERE ? = ?;", parameters+[node]).fetchall() for node in self.nodes}

        # Exception due to data structure
        del ret[1]
        ret[2] = [i for i in ret[2] if i != 2]
        
        return ret
        


class DatabaseWriter(DatabaseReader):
    _connection : sqlite3.Connection

    def __init__(self, database : str):
        self._connection = sqlite3.connect(database)

    def addSNP(self, nodeID, position, ancestral, derived, reference, date, genomeID):
        self._connection.execute("INSERT (?,?,?,?,?,?,?) INTO ?;", [nodeID, position, ancestral, derived, reference, date, genomeID, TABLE_NAME_SNP_ANNOTATION])
    
    def addReference(self, genomeID, genome, strain, genbank, refseq, assemblyName):
        self._connection.execute("INSERT (?,?,?,?,?,?) INTO ?;", [genomeID, genome, strain, genbank, refseq, assemblyName, TABLE_NAME_REFERENCES])