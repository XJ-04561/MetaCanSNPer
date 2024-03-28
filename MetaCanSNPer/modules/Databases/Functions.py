
from MetaCanSNPer.modules.Databases.Globals import *
import MetaCanSNPer.modules.Databases.Globals as Globals
import MetaCanSNPer.modules.Databases.Columns as Columns
from MetaCanSNPer.modules.Databases.Columns import ColumnFlag
from MetaCanSNPer.modules.Databases._Constants import *

import inspect

def private(func):
	print("\n".join(map(repr, inspect.stack()[0])))
	return func

whitespacePattern = re.compile(r"\s+")

def downloadDatabase(databaseName : str, dst : str) -> str:
	from urllib.request import urlretrieve
	
	for source in SOURCES:
		try:
			(filename, msg) = urlretrieve(source.format(databaseName=databaseName), filename=dst) # Throws error if 404
			return filename
		except Exception as e:
			LOGGER.info(f"Database {databaseName!r} not found/accessible on {source!r}.")
			LOGGER.exception(e, stacklevel=logging.INFO)
	LOGGER.error(f"No database named {databaseName!r} found online. Sources tried: {SOURCES}")
	return None

@cache
def generateTableQuery(self, *select : ColumnFlag, orderBy : ColumnFlag|tuple[ColumnFlag,Literal["DESC","ASC"]]|list[tuple[ColumnFlag,Literal["DESC","ASC"]]]=[], **where : Any) -> Generator[tuple[Any],None,None]:
    """All positional arguments should be `ColumnFlag` objects and they are used to
    determine what information to be gathered from the database.
    
    All keyword arguments (except `orderBy`) are the conditions by which each row
    is selected. For example, if you inted to get the row for a specific genbankID
    then you would use the keyword argument as such: `genbankID="GCA_123123123.1"`.
    
    `orderBy` is used to sort the selected data according to `ColumnFlag` with or
    without a direction string ("DESC" or "ASC"). Input can be a list to use many
    columns for sorting, and any item of the list can be with or without a
    direction specified."""
    query = f"SELECT {', '.join(map(Columns.LOOKUP[self._tableName].__getitem__, select))} FROM {self._tableName}"
    params = []
    if where != {}:
        _tmp = []
        for col, val in where.items():
            _tmp.append(f'{col} = ?')
            params.append(val)
        query += f" WHERE {' AND '.join(_tmp)}"
        del _tmp
    
    if type(orderBy) is ColumnFlag:
        orderBy = [(orderBy, "DESC")]
    elif type(orderBy) is tuple:
        orderBy = [orderBy]
    
    if orderBy != []:
        # Create an ordered list of all "ORDER BY X [DIRECTION]"-statements
        orderBy = [tupe if type(tupe) is tuple else (tupe, "DESC") for tupe in orderBy]
        query += f" ORDER BY {', '.join(map(' '.join, orderBy))}"
	
    return query, params

@cache
def generateQuery(*select : ColumnFlag, orderBy : ColumnFlag|tuple[ColumnFlag,Literal["DESC","ASC"]]|list[tuple[ColumnFlag,Literal["DESC","ASC"]]]=[], **where : Any) -> tuple[str,list[Any]]:
	"""All positional arguments should be `ColumnFlag` objects and they are used to
	determine what information to be gathered from the database.
	
	All keyword arguments (except `orderBy`) are the conditions by which each row
	is selected. For example, if you inted to get the row for a specific genbankID
	then you would use the keyword argument as such: `genbankID="GCA_123123123.1"`.
	
	`orderBy` is used to sort the selected data according to `ColumnFlag` with or
	without a direction string ("DESC" or "ASC"). Input can be a list to use many
	columns for sorting, and any item of the list can be with or without a
	direction specified."""

	# Find out which tables are involved in this query. Also makes some necessary assertions for the rest of the function.
	tables = {Columns.UNIQUELOOKUP[col] for col in Columns.UNIQUES.intersection(select, where, map(lambda x : x if type(x) is ColumnFlag else x[0], orderBy))}
	for col in Columns.COMMONS.intersection(select, where, map(lambda x : x[0], orderBy)):
		assert col in Columns.NAMES, f"No such column `{col}` in `Columns.NAMES`"
		if tables.isdisjoint(Columns.COMMONLOOKUP[col]):
			for table in Columns.COMMONLOOKUP[col]:
				if not tables.isdisjoint(Columns.RELATIONSHIPS[table].values()):
					tables.add(table)
					break
			assert table in tables, f"Not implemented to join two disjoint tables, {table=} could not be joined with any of {tables=}"

	def getTable(col):
		for table in tables:
			if col in Columns.LOOKUP[table]:
				return table

	# If selecting all columns then change the selection string into "*", otherwise create a list of "WHERE" statements
	if Columns.ALL in select:
		selection = "*"
	else:
		selection = [f"{table}.{Columns.LOOKUP[table][col]}" for table, col in zip(map(getTable, select), select)]

	# Simply join all the tables into a comma(+space) separated list
	source = ", ".join(tables)
	
	# Create "WHERE"-statements that are meant to show how the tables are connected, like: table1.colA = table2.colB
	connections = []
	for i, table in enumerate(tables):
		for otherTable in tables.intersection(Columns.RELATIONSHIPS[table].values()):
			col = Columns.RELATIONS[table,otherTable]
			connections.append(f"{table}.{Columns.LOOKUP[table][col]} = {table}.{Columns.LOOKUP[otherTable][col]}")
			break

	params = list(where.values())
	conditions = " AND ".join([f"{table}.{Columns.LOOKUP[table][col]} = ?" for table, col in zip(map(getTable, params), params)])

	# Create an ordered list of all "ORDER BY X [DIRECTION]"-statements
	if type(orderBy) is ColumnFlag:
		orderBy = [(orderBy, "DESC")]
	elif type(orderBy) is tuple:
		orderBy = [orderBy]
	if len(orderBy) > 0:
		orderBy = [tupe if type(tupe) is tuple else (tupe, "DESC") for tupe in orderBy]
		
		andTable = lambda col, direction : (getTable(col), col, direction)

		keyColumn = ", ".join([f"{table}.{Columns.LOOKUP[table][col]} {direction}" for table, col, direction in map(andTable, orderBy)])
	
		return f"SELECT {selection} FROM {source} WHERE {conditions} ORDER BY {keyColumn};", params
	else:
		return f"SELECT {selection} FROM {source} WHERE {conditions};", params