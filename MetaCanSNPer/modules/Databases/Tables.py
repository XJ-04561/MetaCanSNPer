

from MetaCanSNPer.modules.Databases.Globals import *
from MetaCanSNPer.modules.Databases import *

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
	def get(self, *columnsToGet : ColumnFlag, nodeID : int=None, snpID : str=None, genomeID : int=None, position : int=None, ancestral : Literal["A","T","C","G"]=None, derived : Literal["A","T","C","G"]=None, snpReference : str=None, date : str=None, genome : str=None, strain : str=None, genbankID : str=None, refseqID : str=None, assembly : str=None, chromosome : str=None):
		pass
	
	@final
	def get(self, *select : ColumnFlag, orderBy : ColumnFlag|tuple[ColumnFlag,Literal["DESC","ASC"]]|list[tuple[ColumnFlag,Literal["DESC","ASC"]]]=[], **where : Any):
		
		T4 = [TABLE_NAME_NODES, TABLE_NAME_REFERENCES, TABLE_NAME_SNP_ANNOTATION, TABLE_NAME_TREE]
		T4.remove(self._tableName)
		
		if all(col in Columns.LOOKUP[self._tableName] for col in select):
			# All values needed are in this table.
			if Columns.ALL in select:
				selection = "*"
			else:
				selection = ", ".join([Columns.LOOKUP[self._tableName][col] for col in select])
			
			source = self._tableName

			vars = list(where.keys())
			
			conditions = " AND ".join([f"{Columns.LOOKUP[self._tableName][col]} = ?" for col in vars])

			
			if type(orderBy) is ColumnFlag:
				orderBy = [(orderBy, False)]
			elif type(orderBy) is tuple:
				orderBy = [orderBy]
				
			order = [f"{Columns.LOOKUP[self._tableName][col]} {direction}" for col, direction in orderBy]
			
			keyColumn = f" ORDER BY {', '.join(order)}" if order != [] else ""

			params = [where[col] for col in vars]

		elif all(col in Columns.RELATIONSHIPS[self._tableName] or col in Columns.LOOKUP[self._tableName] for col in select):
			# Selection requires joining with other table.
			thisTable = Columns.LOOKUP[self._tableName]
			relations = Columns.RELATIONSHIPS[self._tableName]

			if Columns.ALL in select:
				selection = "*"
				tables = []
				for table, col in map(lambda col : thisTable.get(col) or relations.get(col), select):
					if table is None: continue
					if table not in tables:
						tables.append(table)
			else:
				selection = []
				tables = []
				for table, col in map(lambda col : thisTable.get(col) or relations.get(col), select):
					selection.append(f"{table}.{Columns.LOOKUP[table][col]}")
					if table not in tables:
						tables.append(table)
				selection = ", ".join(selection)

			source = ", ".join(tables)

			vars = list(where.keys())
			
			conditions = []
			params = []
			for col in vars:
				if col in thisTable:
					conditions.append(f"{self._tableName}.{Columns.LOOKUP[self._tableName][col]} = ?")
					params.append(where[col])
				else:
					otherTable, commonCol = relations[col]
					conditions.append(f"{self._tableName}.{Columns.LOOKUP[self._tableName][commonCol]} = {otherTable}.{Columns.LOOKUP[otherTable][commonCol]}")
					conditions.append(f"{otherTable}.{Columns.LOOKUP[otherTable][col]} = ?")
					params.append(where[col])
			conditions = " AND ".join(conditions)

			
			if type(orderBy) is ColumnFlag:
				orderBy = [(orderBy, "DESC")]
			elif type(orderBy) is tuple:
				orderBy = [orderBy]
			order = []	
			for col, direction in orderBy:
				if col in thisTable:
					order.append(f"{Columns.LOOKUP[self._tableName][col]} {direction}")
				else:
					order.append(f"{Columns.LOOKUP[relations[col][0]][col]} {direction}")
			
			keyColumn = f" ORDER BY {', '.join(order)}" if order != [] else ""

			params = [where[col] for col in vars]
		else:
			raise NotImplementedError(f"Can't find values={[Columns.NAMES[col] for col in select]} associated with column values=({[key+'='+str(val) for key, val in where]}) in table={self._tableName}")
		
		self._conn.execute(f"SELECT {selection} FROM {source} WHERE {conditions} ORDER BY {keyColumn};", params)

	@overload
	def first(self, *columnsToGet : ColumnFlag, nodeID : int=None, snpID : str=None, genomeID : int=None, position : int=None, ancestral : Literal["A","T","C","G"]=None, derived : Literal["A","T","C","G"]=None, snpReference : str=None, date : str=None, genome : str=None, strain : str=None, genbankID : str=None, refseqID : str=None, assembly : str=None, chromosome : str=None):
		pass
	
	@final
	def first(self, *select : ColumnFlag, orderBy : ColumnFlag|tuple[ColumnFlag,Literal["DESC","ASC"]]|list[tuple[ColumnFlag,Literal["DESC","ASC"]]]=[], **where : Any):
		for row in self.get(*select, orderBy=orderBy, **where):
			return row
	
	@overload
	def all(self, *columnsToGet : ColumnFlag, nodeID : int=None, snpID : str=None, genomeID : int=None, position : int=None, ancestral : Literal["A","T","C","G"]=None, derived : Literal["A","T","C","G"]=None, snpReference : str=None, date : str=None, genome : str=None, strain : str=None, genbankID : str=None, refseqID : str=None, assembly : str=None, chromosome : str=None):
		pass
	
	@final
	def all(self, *select : ColumnFlag, orderBy : ColumnFlag|tuple[ColumnFlag,Literal["DESC","ASC"]]|list[tuple[ColumnFlag,Literal["DESC","ASC"]]]=[], **where : Any):
		for row in self.get(*select, orderBy=orderBy, **where):
			yield row

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
