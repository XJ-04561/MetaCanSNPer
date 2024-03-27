
from MetaCanSNPer.modules.Databases.Globals import *
from MetaCanSNPer.modules.Databases.Tables import SNPTable, ReferenceTable, NodeTable, TreeTable
import MetaCanSNPer.modules.Databases.Columns as Columns
from MetaCanSNPer.modules.Databases.Tree import Branch
			

class DatabaseReader:
	_connection : sqlite3.Connection

	filename : str

	def __init__(self, database : str):
		if not pExists(database):
			raise FileNotFoundError(f"Database file {database} not found on the system.")
		self.filename = database
		# Convert to URI acceptable filename
		cDatabase = "/".join(filter(lambda s : s != "", database.replace('?', '%3f').replace('#', '%23').split(os.path.sep)))
		if not cDatabase.startswith("/"): # Path has to be absolute already, and windows paths need a prepended '/'
			cDatabase = "/"+cDatabase
		try:
			self._connection = sqlite3.connect(f"file:{cDatabase}?immutable=1", uri=True)
		except Exception as e:
			LOGGER.error("Failed to connect to database using URI: "+f"file:{cDatabase}?immutable=1")
			raise e
		
		self.SNPTable = SNPTable(self._connection, "r")
		self.ReferenceTable = ReferenceTable(self._connection, "r")
		self.NodeTable = NodeTable(self._connection, "r")
		self.TreeTable = TreeTable(self._connection, "r")

	def __del__(self):
		try:
			self._connection.close()
		except:
			pass
	
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


class DatabaseWriter(DatabaseReader):
	_connection : sqlite3.Connection

	def __init__(self, database : str, mode : str):
		self._connection = sqlite3.connect(database)

		self.SNPTable = SNPTable(self._connection, mode=mode)
		self.ReferenceTable = ReferenceTable(self._connection, mode=mode)
		self.NodeTable = NodeTable(self._connection, mode=mode)
		self.TreeTable = TreeTable(self._connection, mode=mode)

	def addSNP(self, nodeID, position, ancestral, derived, reference, date, genomeID):
		self._connection.execute(f"INSERT (?,?,?,?,?,?,?) INTO {TABLE_NAME_SNP_ANNOTATION};", [nodeID, position, ancestral, derived, reference, date, genomeID])
	
	def addReference(self, genomeID, genome, strain, genbank, refseq, assemblyName):
		self._connection.execute(f"INSERT (?,?,?,?,?,?) INTO {TABLE_NAME_REFERENCES};", [genomeID, genome, strain, genbank, refseq, assemblyName])

	def addNode(self, nodeID, genoType):
		self._connection.execute(f"INSERT (?,?) INTO {TABLE_NAME_NODES};", [nodeID, genoType])

	def addBranch(self, parentID, childID, rank):
		self._connection.execute(f"INSERT (?,?,?) INTO {TABLE_NAME_TREE};", [parentID, childID, rank])