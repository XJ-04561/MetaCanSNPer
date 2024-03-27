
from MetaCanSNPer.modules.Databases.Globals import *
from MetaCanSNPer.modules.Databases import *
from MetaCanSNPer.modules.Databases.Functions import generateTableQuery

class TableDefinitionMissmatch(Exception): pass

class Table:

	_conn : sqlite3.Connection
	_tableName : str
	_columns : list[str]
	_types : list[tuple[str]]
	_appendRows : list[str]

	def __init__(self, conn : sqlite3.Connection, mode : str):
		self._conn = conn
		if mode.lower() in ["r", "w", "a"]:
			self._mode = mode.lower()
		else:
			ValueError(f"Table argument 'mode' must be either 'r', 'a', or 'w', not {mode!r}.")
		
		try:
			pragma = self._conn.execute(f"PRAGMA table_info({self._tableName});").fetchall()
			if len(pragma) != len(self._columns):
				raise TableDefinitionMissmatch(f"Column count does not match between database and {SOFTWARE_NAME}. {SOFTWARE_NAME} had columns: {self._columns} and Database had PRAGMA: {pragma}")
			for (_, name, dataType, _, _, isPrimaryKey), (colName, colType) in zip(pragma, zip(self._columns, self._types)):
				if name != colName: raise TableDefinitionMissmatch(f"Column name in database ({name}) not consistent with name in {SOFTWARE_NAME} ({colName})")
				if dataType != colType: raise TableDefinitionMissmatch(f"Column data type in database ({dataType}) not consistent with data type in {SOFTWARE_NAME} ({colType})")
				if isPrimaryKey and colType not in colType: raise TableDefinitionMissmatch(f"Column is PRIMARY KEY in database but should not be for compliance with {SOFTWARE_NAME}.")
				if not isPrimaryKey and colType in colType: raise TableDefinitionMissmatch(f"Column is not PRIMARY KEY in database but should be for compliance with {SOFTWARE_NAME}.")
		except Exception as e:
			LOGGER.exception(e)
			if self._mode == "r" or type(e) is not TableDefinitionMissmatch: 
				raise e
			
			if self._mode == "a": self._conn.execute(f"ALTER TABLE {self._tableName} RENAME TO {self._tableName}2;")
			
			queryString = [f"CREATE TABLE {self._tableName} ("]
			
			for name, colType in zip(self._columns, self._types):
				queryString.append( f"{name} {' '.join(colType)},")
			queryString += self._appendRows
			queryString.append(");")

			self._conn.execute("\n".join(queryString))

			if self._mode == "a":
				self._conn.execute(f"INSERT INTO {self._tableName} SELECT * FROM {self._tableName}2;")

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