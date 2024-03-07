
from functools import cached_property
import sqlite3

from MetaCanSNPer.Globals import *

class Branch:


    _connection : sqlite3.Connection
    nodeID : int
    parameters = [
        TREE_COLUMN_CHILD,
        TABLE_NAME_TREE,
        TREE_COLUMN_PARENT
    ]

    def __init__(self, connection : sqlite3.Connection, nodeID : int):
        self._connection = connection
        self.nodeID = nodeID
    
    def children(self):
        return [Branch(self._connection, childID) for childID in self._connection.execute("SELECT ? FROM ? WHERE ? = ?;", self.parameters+[self.nodeID]).fetchall()]

class DatabaseReader:
    _connection : sqlite3.Connection

    references : list[tuple[int,str,str,str,str]]
    SNPs : list[tuple[str,int,str,str]]

    def __init__(self, database : str):
        self._connection = sqlite3.connect(f"file:{database}?mode=ro", uri=True)

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
        '''{SNP_ID : (POSITION, ANCESTRAL, DERIVED)}'''
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
        return Branch(self._connection, self._connection.execute("SELECT ? FROM ? WHERE ? = ?;", parameters+[2]).fetch())
        


class DatabaseWriter(DatabaseReader):
    _connection : sqlite3.Connection

    def __init__(self, database : str):
        self._connection = sqlite3.connect(database)

    def addSNP(self, nodeID, position, ancestral, derived, reference, date, genomeID):
        self._connection.execute("INSERT (?,?,?,?,?,?,?) INTO ?;", [nodeID, position, ancestral, derived, reference, date, genomeID, TABLE_NAME_SNP_ANNOTATION])
    
    def addReference(self, genomeID, genome, strain, genbank, refseq, assemblyName):
        self._connection.execute("INSERT (?,?,?,?,?,?) INTO ?;", [genomeID, genome, strain, genbank, refseq, assemblyName, TABLE_NAME_REFERENCES])