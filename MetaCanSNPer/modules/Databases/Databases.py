
from MetaCanSNPer.modules.Databases.Globals import *
import MetaCanSNPer.modules.Databases.Globals as Globals
import MetaCanSNPer.modules.Databases.Columns as Columns
from MetaCanSNPer.modules.Databases.Columns import ColumnFlag
from MetaCanSNPer.modules.Databases._Constants import *

from MetaCanSNPer.modules.Databases.Tables import SNPTable, ReferenceTable, NodeTable, TreeTable, RankTable, GenomesTable
from MetaCanSNPer.modules.Databases.Tree import Branch
from MetaCanSNPer.modules.Databases.Functions import generateQuery

class Database:

	_connection : sqlite3.Connection
	_mode : str
	filename : str

	def __init__(self, database : sqlite3.Connection):
		self.filename = database.execute("PRAGMA database_list;").fetchone()[2]
		self._connection = database

		self.SNPTable = SNPTable(self._connection, self._mode)
		self.ReferenceTable = ReferenceTable(self._connection, self._mode)
		self.NodeTable = NodeTable(self._connection, self._mode)
		self.TreeTable = TreeTable(self._connection, self._mode)
		self.RankTable = RankTable(self._connection, self._mode)
		self.GenomesTable = GenomesTable(self._connection, self._mode)

		if set(TABLES).isdisjoint([table for (table,) in self._connection.execute("SELECT tableName FROM sqlite_master WHERE type='table';")]):
			# Table is new
			self.SNPTable.create()
			self.ReferenceTable.create()
			self.NodeTable.create()
			self.TreeTable.create()
			self.RankTable.create()
			self.GenomesTable.create()
		elif self.schemaHash != Globals.DATABASE_VERSION_HASH:
			LOGGER.warning(f"Database version does not match the currently set version. (user_version={database.execute('PRAGMA user_version;').fetchone()[0]})")
			if Globals.STRICT:
				raise sqlite3.DatabaseError(f"Database version does not match the currently set version. (user_version={database.execute('PRAGMA user_version;').fetchone()[0]})")
			elif self._mode == "w":
				# Transfer data from old tables into new tables
				self.SNPTable.recreate()
				self.ReferenceTable.recreate()
				self.NodeTable.recreate()
				self.TreeTable.recreate()
				self.RankTable.recreate()
				self.GenomesTable.recreate()
	
	def __enter__(self):
		return self
	
	def __exit__(self):
		try:
			self.commit()
		except:
			pass
		del self
	
	@property
	def __version__(self):
		return DATABASE_VERSIONS.get(self.schemaHash, "Unknown")

	def __del__(self):
		try:
			self._connection.close()
		except:
			pass
	
	def __repr__(self):
		return object.__repr__(self)[:-1] + f" version={self.__version__} schemaHash={self.schemaHash} tables={list(zip(TABLES,map(len, TABLES)))}>"
	
	@overload
	def get(self, *columnsToGet : ColumnFlag, orderBy : ColumnFlag|tuple[ColumnFlag,Literal["DESC","ASC"]]|list[tuple[ColumnFlag,Literal["DESC","ASC"]]]=[], nodeID : int=None, snpID : str=None, genomeID : int=None, position : int=None, ancestral : Literal["A","T","C","G"]=None, derived : Literal["A","T","C","G"]=None, snpReference : str=None, date : str=None, genome : str=None, strain : str=None, genbankID : str=None, refseqID : str=None, assembly : str=None, chromosome : str=None) -> Generator[tuple[Any],None,None]:
		pass

	@final
	def get(self, *select : ColumnFlag, orderBy : ColumnFlag|tuple[ColumnFlag,Literal["DESC","ASC"]]|list[tuple[ColumnFlag,Literal["DESC","ASC"]]]=[], **where : Any) -> Generator[tuple[Any],None,None]:
		for row in self._connection.execute(*generateQuery(*select, orderBy=orderBy, **where)):
			yield row
	
	@property
	def SNPs(self) -> Generator[tuple[str,int,str,str],None,None]:
		return self.SNPTable.get(Columns.ALL)

	@property
	def references(self) -> Generator[tuple[int,str,str,str,str],None,None]:
		return self.ReferenceTable.get(Columns.ALL)

	@property
	def nodes(self) -> Generator[tuple[int,str],None,None]:
		return self.NodeTable.get(Columns.ALL)

	@cached_property
	def tree(self) -> Branch:
		"""{nodeID:[child1, child2, ...]}"""
		for (nodeID,) in self._connection.execute(f"SELECT {TREE_COLUMN_CHILD} FROM {TABLE_NAME_TREE} WHERE {TREE_COLUMN_PARENT} = ?;", [2]):
			if nodeID != 2:
				return Branch(self._connection, nodeID)

	@property
	def schemaHash(self):
		return hashlib.md5(
			whitespacePattern.sub(
				" ",
				"; ".join([
					x[0]
					for x in self._connection.execute(f"SELECT sql FROM sqlite_schema;")
					if type(x) is tuple and x[0] is not None
				])
			).encode("utf-8")
		).hexdigest()

class DatabaseReader(Database):
	
	_mode = "r"


class DatabaseWriter(Database):

	_mode = "w"

	def addSNP(self, nodeID, position, ancestral, derived, reference, date, genomeID):
		self._connection.execute(f"INSERT (?,?,?,?,?,?,?) INTO {TABLE_NAME_SNP_ANNOTATION};", [nodeID, position, ancestral, derived, reference, date, genomeID])
	
	def addReference(self, genomeID, genome, strain, genbank, refseq, assemblyName):
		self._connection.execute(f"INSERT (?,?,?,?,?,?) INTO {TABLE_NAME_REFERENCES};", [genomeID, genome, strain, genbank, refseq, assemblyName])

	def addNode(self, nodeID, genoType):
		self._connection.execute(f"INSERT (?,?) INTO {TABLE_NAME_NODES};", [nodeID, genoType])

	def addBranch(self, parentID, childID, rank):
		self._connection.execute(f"INSERT (?,?,?) INTO {TABLE_NAME_TREE};", [parentID, childID, rank])

	def addRank(self, rankID, rankName):
		self._connection.execute(f"INSERT (?,?) INTO {TABLE_NAME_RANKS};", [rankID, rankName])

	def addGenome(self, genomeID, genomeName):
		self._connection.execute(f"INSERT (?,?) INTO {TABLE_NAME_GENOMES};", [genomeID, genomeName])
	
	def commit(self):
		self._connection.commit()

@overload
def openDatabase(database : str, mode : Literal["r"]) -> DatabaseReader:
	pass

@final
def openDatabase(database : str, mode : Literal["w"]) -> DatabaseWriter:
	match mode:
		case "r":
			if not pExists(database):
				raise FileNotFoundError(f"Database file {database} not found on the system.")
			# Convert to URI acceptable filename
			cDatabase = "/".join(filter(lambda s : s != "", database.replace('?', '%3f').replace('#', '%23').split(os.path.sep)))
			if not cDatabase.startswith("/"): # Path has to be absolute already, and windows paths need a prepended '/'
				cDatabase = "/"+cDatabase
			try:
				return DatabaseReader(sqlite3.connect(f"file:{cDatabase}?immutable=1", uri=True))
			except Exception as e:
				LOGGER.error("Failed to connect to database using URI: "+f"file:{cDatabase}?immutable=1")
				raise e
			
		case "w":
			return DatabaseWriter(sqlite3.connect(database))
