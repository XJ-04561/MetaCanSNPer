
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
	
	@property
	def children(self):
		return [Branch(self._connection, childID) for (childID,) in self._connection.execute(f"SELECT {TREE_COLUMN_CHILD} FROM {TABLE_NAME_TREE} WHERE {TREE_COLUMN_PARENT} = ?", [self.nodeID]).fetchall()]

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
		return self._connection.execute(f"SELECT {REFERENCE_COLUMN_GENOME_ID} FROM {TABLE_NAME_REFERENCES} WHERE {REFERENCE_COLUMN_GENOME} = ?", [genome]).fetchone()

	@cached_property
	def references(self) -> list[tuple[int,str,str,str,str]]:
		''''''
		return self._connection.execute( f"SELECT {REFERENCE_COLUMN_GENOME_ID}, {REFERENCE_COLUMN_GENOME}, {REFERENCE_COLUMN_GENBANK}, {REFERENCE_COLUMN_REFSEQ}, {REFERENCE_COLUMN_ASSEMBLY} FROM {TABLE_NAME_REFERENCES} ORDER BY {REFERENCE_COLUMN_GENOME_ID} ASC").fetchall()

	def _getSNPsByGenomeId(self, genomeID : int) -> sqlite3.Cursor:
		return self._connection.execute(f"SELECT {SNP_COLUMN_SNP_ID}, {SNP_COLUMN_POSITION}, {SNP_COLUMN_ANCESTRAL}, {SNP_COLUMN_DERIVED} FROM {TABLE_NAME_SNP_ANNOTATION} WHERE {SNP_COLUMN_GENOME_ID} = ? ORDER BY {SNP_COLUMN_POSITION} ASC", [genomeID])
	
	@cached_property
	def SNPs(self) -> sqlite3.Cursor:
		return self._connection.execute(f"SELECT {SNP_COLUMN_SNP_ID}, {SNP_COLUMN_POSITION}, {SNP_COLUMN_ANCESTRAL}, {SNP_COLUMN_DERIVED} FROM {TABLE_NAME_SNP_ANNOTATION} ORDER BY {SNP_COLUMN_POSITION} ASC")

	@cached_property
	def SNPsByGenome(self) -> dict[str,list[tuple[str,int,str,str]]]:
		return {genome:self._getSNPsByGenomeId(genome_id) for genome_id, genome, _, _, _ in self.references}

	@cached_property
	def SNPsByID(self) -> dict[str,list[tuple[str,int,str,str]]]:
		'''{SNP_ID : (POSITION, ANCESTRAL, DERIVED)}'''
		return {snpID:(pos, anc, der) for snpID, pos, anc, der in self.SNPs}
	
	@cached_property
	def nodes(self) -> dict[str,list[tuple[str,int,str,str]]]:
		return dict(self._connection.execute(f"SELECT {NODE_COLUMN_ID}, {NODE_COLUMN_NAME} FROM {TABLE_NAME_NODES}"))
	
	def node(self, nodeID):
		for (nodeSNPName, ) in self._connection.execute(f"SELECT {NODE_COLUMN_NAME} FROM {TABLE_NAME_NODES} WHERE {NODE_COLUMN_ID} = ?", [nodeID]):
			return nodeSNPName

	@cached_property
	def tree(self) -> Branch:
		"""{nodeID:[child1, child2, ...]}"""
		for (nodeID,) in self._connection.execute(f"SELECT {TREE_COLUMN_CHILD} FROM {TABLE_NAME_TREE} WHERE {TREE_COLUMN_PARENT} = ?", [2]):
			if nodeID != 2:
				return Branch(self._connection, nodeID)
		


class DatabaseWriter(DatabaseReader):
	_connection : sqlite3.Connection

	def __init__(self, database : str):
		self._connection = sqlite3.connect(database)

	def addSNP(self, nodeID, position, ancestral, derived, reference, date, genomeID):
		self._connection.execute("INSERT (?,?,?,?,?,?,?) INTO ?", [nodeID, position, ancestral, derived, reference, date, genomeID, TABLE_NAME_SNP_ANNOTATION])
	
	def addReference(self, genomeID, genome, strain, genbank, refseq, assemblyName):
		self._connection.execute("INSERT (?,?,?,?,?,?) INTO ?", [genomeID, genome, strain, genbank, refseq, assemblyName, TABLE_NAME_REFERENCES])