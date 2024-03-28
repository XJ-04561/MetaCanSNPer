
from MetaCanSNPer.modules.Databases.Globals import *
import MetaCanSNPer.modules.Databases.Globals as Globals
import MetaCanSNPer.modules.Databases.Columns as Columns
from MetaCanSNPer.modules.Databases.Columns import ColumnFlag
from MetaCanSNPer.modules.Databases._Constants import *

from MetaCanSNPer.modules.Databases.Functions import generateTableQuery

class TableDefinitionMissmatch(Exception): pass

class Table:

	_conn : sqlite3.Connection
	_tableName : str
	_columns : list[str]
	_types : list[tuple[str]]
	_appendRows : list[str]

	def __init__(self, conn : sqlite3.Connection, mode : Literal["r","w","a"]):
		self._conn = conn
		self._mode = mode

	def __len__(self):
		self._conn.execute(f"SELECT COUNT(*) FROM {self._tableName};").fetchone()[0]
	
	def __repr__(self):
		return object.__repr__(self)[:-1] + f" rows={len(self)} columns={self._columns}>"

	def create(self) -> bool:
		queryString = [f"{name} {' '.join(colType)}" for name, colType in zip(self._columns, self._types)]
		queryString += self._appendRows

		self._conn.execute(f"CREATE TABLE IF NOT EXISTS {self._tableName} (\n\t\t{',\n\t\t'.join(queryString)}\n);")
	
	def recreate(self) -> bool:
		try:
			self._conn.execute(f"ALTER TABLE {self._tableName} RENAME TO {self._tableName}2;")
		except:
			# Table doesn't exist
			pass
		
		queryString = [f"{name} {' '.join(colType)}," for name, colType in zip(self._columns, self._types)]
		queryString += self._appendRows

		self._conn.execute(f"CREATE TABLE {self._tableName} (\n\t\t{',\n\t\t'.join(queryString)}\n);")
		
		try:
			self._conn.execute(f"INSERT INTO {self._tableName} SELECT * FROM {self._tableName}2;")
		except:
			pass

	@overload
	def get(self, *columnsToGet : ColumnFlag, orderBy : ColumnFlag|tuple[ColumnFlag,Literal["DESC","ASC"]]|list[tuple[ColumnFlag,Literal["DESC","ASC"]]]=[], nodeID : int=None, snpID : str=None, genomeID : int=None, position : int=None, ancestral : Literal["A","T","C","G"]=None, derived : Literal["A","T","C","G"]=None, snpReference : str=None, date : str=None, genome : str=None, strain : str=None, genbankID : str=None, refseqID : str=None, assembly : str=None, chromosome : str=None) -> Generator[tuple[Any],None,None]:
		pass
	
	@final
	def get(self, *select : ColumnFlag, orderBy : ColumnFlag|tuple[ColumnFlag,Literal["DESC","ASC"]]|list[tuple[ColumnFlag,Literal["DESC","ASC"]]]=[], **where : Any) -> Generator[tuple[Any],None,None]:
		for row in self._conn.execute(*generateTableQuery(*select, orderBy=orderBy, **where)):
			yield row

	@overload
	def first(self, *columnsToGet : ColumnFlag, nodeID : int=None, snpID : str=None, genomeID : int=None, position : int=None, ancestral : Literal["A","T","C","G"]=None, derived : Literal["A","T","C","G"]=None, snpReference : str=None, date : str=None, genome : str=None, strain : str=None, genbankID : str=None, refseqID : str=None, assembly : str=None, chromosome : str=None) -> tuple[Any]:
		pass
	
	@final
	def first(self, *select : ColumnFlag, orderBy : ColumnFlag|tuple[ColumnFlag,Literal["DESC","ASC"]]|list[tuple[ColumnFlag,Literal["DESC","ASC"]]]=[], **where : Any) -> tuple[Any]:
		for row in self.get(*select, orderBy=orderBy, **where):
			return row
	
	@overload
	def all(self, *columnsToGet : ColumnFlag, nodeID : int=None, snpID : str=None, genomeID : int=None, position : int=None, ancestral : Literal["A","T","C","G"]=None, derived : Literal["A","T","C","G"]=None, snpReference : str=None, date : str=None, genome : str=None, strain : str=None, genbankID : str=None, refseqID : str=None, assembly : str=None, chromosome : str=None) -> list[tuple[Any]]:
		pass
	
	@final
	def all(self, *select : ColumnFlag, orderBy : ColumnFlag|tuple[ColumnFlag,Literal["DESC","ASC"]]|list[tuple[ColumnFlag,Literal["DESC","ASC"]]]=[], **where : Any) -> list[tuple[Any]]:
		return list(self.get(*select, orderBy=orderBy, **where))

Table.get.__doc__ = Table.first.__doc__ = Table.all.__doc__ = generateTableQuery.__doc__

class SNPTable(Table):

	_tableName = TABLE_NAME_SNP_ANNOTATION
	_columns = [
		SNP_COLUMN_NODE_ID,
		SNP_COLUMN_SNP_ID,
		SNP_COLUMN_POSITION,
		SNP_COLUMN_ANCESTRAL,
		SNP_COLUMN_DERIVED,
		SNP_COLUMN_REFERENCE,
		SNP_COLUMN_DATE,
		SNP_COLUMN_GENOME_ID
	]
	_types = [
		SNP_COLUMN_NODE_ID_TYPE,
		SNP_COLUMN_SNP_ID_TYPE,
		SNP_COLUMN_POSITION_TYPE,
		SNP_COLUMN_ANCESTRAL_TYPE,
		SNP_COLUMN_DERIVED_TYPE,
		SNP_COLUMN_REFERENCE_TYPE,
		SNP_COLUMN_DATE_TYPE,
		SNP_COLUMN_GENOME_ID_TYPE
	]
	_appendRows = [
		f"FOREIGN KEY ({SNP_COLUMN_GENOME_ID}) REFERENCES {TABLE_NAME_REFERENCES} ({REFERENCE_COLUMN_GENOME_ID})"
	]

class ReferenceTable(Table):

	_tableName = TABLE_NAME_REFERENCES
	_columns = [
		REFERENCE_COLUMN_GENOME_ID,
		REFERENCE_COLUMN_GENOME,
		REFERENCE_COLUMN_STRAIN,
		REFERENCE_COLUMN_GENBANK,
		REFERENCE_COLUMN_REFSEQ,
		REFERENCE_COLUMN_ASSEMBLY
	]
	_types = [
		REFERENCE_COLUMN_GENOME_ID_TYPE,
		REFERENCE_COLUMN_GENOME_TYPE,
		REFERENCE_COLUMN_STRAIN_TYPE,
		REFERENCE_COLUMN_GENBANK_TYPE,
		REFERENCE_COLUMN_REFSEQ_TYPE,
		REFERENCE_COLUMN_ASSEMBLY_TYPE
	]
	_appendRows = []

class NodeTable(Table):

	_tableName = TABLE_NAME_NODES
	_columns = [
		NODE_COLUMN_ID,
		NODE_COLUMN_NAME
	]
	_types = [
		NODE_COLUMN_ID_TYPE,
		NODE_COLUMN_NAME_TYPE
	]
	_appendRows = []

class TreeTable(Table):

	_tableName = TABLE_NAME_TREE
	_columns = [
		TREE_COLUMN_PARENT,
		TREE_COLUMN_CHILD,
		TREE_COLUMN_RANK,
	]
	_types = [
		TREE_COLUMN_PARENT_TYPE,
		TREE_COLUMN_CHILD_TYPE,
		TREE_COLUMN_RANK_TYPE
	]
	_appendRows = [
		f"FOREIGN KEY ({TREE_COLUMN_PARENT}) REFERENCES {TABLE_NAME_NODES} ({NODE_COLUMN_ID})"
		f"FOREIGN KEY ({TREE_COLUMN_CHILD}) REFERENCES {TABLE_NAME_NODES} ({NODE_COLUMN_ID})"
		f"FOREIGN KEY ({TREE_COLUMN_RANK}) REFERENCES {TABLE_NAME_RANKS} ({RANKS_COLUMN_ID})"
		f"unique ({TREE_COLUMN_PARENT}, {TREE_COLUMN_CHILD})"
	]

class RankTable(Table):
	
	_tableName = TABLE_NAME_RANKS
	_columns = [
		RANKS_COLUMN_ID,
		RANKS_COLUMN_RANK
	]
	_types = [
		RANKS_COLUMN_ID_TYPE,
		RANKS_COLUMN_RANK_TYPE
	]
	_appendRows = []

class GenomesTable(Table):
	
	_tableName = TABLE_NAME_RANKS
	_columns = [		
		GENOMES_COLUMN_ID,
		GENOMES_COLUMN_NAME
	]
	_types = [
		GENOMES_COLUMN_ID_TYPE,
		GENOMES_COLUMN_NAME_TYPE
	]
	_appendRows = [
		f"FOREIGN KEY ({GENOMES_COLUMN_ID}) REFERENCES {TABLE_NAME_NODES} ({NODE_COLUMN_ID})"
	]