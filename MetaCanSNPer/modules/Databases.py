
from functools import cached_property
import sqlite3
from typing import Generator, Callable, Iterable, Self

from MetaCanSNPer.Globals import *
from MetaCanSNPer.modules.LogKeeper import createLogger

LOGGER = createLogger(__name__)


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
	def children(self) -> list[Self]:
		return [Branch(self._connection, childID) for (childID,) in self._connection.execute(f"SELECT {TREE_COLUMN_CHILD} FROM {TABLE_NAME_TREE} WHERE {TREE_COLUMN_PARENT} = ?", [self.nodeID]).fetchall()]

class DatabaseReader:
	_connection : sqlite3.Connection

	filename : str

	references : list[tuple[int,str,str,str,str]]
	SNPs : list[tuple[str,int,str,str]]

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
		LOGGER.debug(f'"SELECT {SNP_COLUMN_SNP_ID}, {SNP_COLUMN_POSITION}, {SNP_COLUMN_ANCESTRAL}, {SNP_COLUMN_DERIVED} FROM {TABLE_NAME_SNP_ANNOTATION} WHERE {SNP_COLUMN_GENOME_ID} = ? ORDER BY {SNP_COLUMN_POSITION} ASC", [genomeID]')
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
	
	def SNPsByNode(self, nodeID : int) -> Generator[tuple[str,tuple[int,str,str]], None, None]:
		for (nodeSNPName, ) in self._connection.execute(f"SELECT {NODE_COLUMN_NAME} FROM {TABLE_NAME_NODES} WHERE {NODE_COLUMN_ID} = ?", [nodeID]):
			yield (nodeSNPName, self.SNPsByID[nodeSNPName])

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