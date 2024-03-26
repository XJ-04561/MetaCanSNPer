
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
	def children(self) -> Generator[Self,None,None]:
		for (childID,) in self._connection.execute(f"SELECT {TREE_COLUMN_CHILD} FROM {TABLE_NAME_TREE} WHERE {TREE_COLUMN_PARENT} = ?", [self.nodeID]):
			yield Branch(self._connection, childID)

# class Table:
# 	_connection : sqlite3.Connection
# 	tableName : str

# 	def __init__(self, conn : sqlite3.Connection, tableName : str):
# 		self.conn = conn
# 		self.tableName = tableName

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

	def genomeID(self, genome : str) -> int:
		for (genomeID,) in self._connection.execute(f"SELECT {REFERENCE_COLUMN_GENOME_ID} FROM {TABLE_NAME_REFERENCES} WHERE {REFERENCE_COLUMN_GENOME} = ?", [genome]):
			return genomeID
		raise ValueError(f"No genome id found for {genome=}")

	@property
	def references(self) -> Generator[tuple[int,str,str,str,str],None,None]:
		''''''
		for genomeID, genome, genbankID, refseqID, assemblyName in self._connection.execute( f"SELECT {REFERENCE_COLUMN_GENOME_ID}, {REFERENCE_COLUMN_GENOME}, {REFERENCE_COLUMN_GENBANK}, {REFERENCE_COLUMN_REFSEQ}, {REFERENCE_COLUMN_ASSEMBLY} FROM {TABLE_NAME_REFERENCES} ORDER BY {REFERENCE_COLUMN_GENOME_ID} ASC"):
			yield genomeID, genome, genbankID, refseqID, assemblyName

	def _getSNPsByGenomeId(self, genomeID : int) -> Generator[tuple[str,int,str,str],None,None]:
		for snpID, pos, anc, der in self._connection.execute(f"SELECT {SNP_COLUMN_SNP_ID}, {SNP_COLUMN_POSITION}, {SNP_COLUMN_ANCESTRAL}, {SNP_COLUMN_DERIVED} FROM {TABLE_NAME_SNP_ANNOTATION} WHERE {SNP_COLUMN_GENOME_ID} = ? ORDER BY {SNP_COLUMN_POSITION} ASC", [genomeID]):
			yield snpID, pos, anc, der
	
	@property
	def SNPs(self) -> Generator[tuple[str,int,str,str],None,None]:
		for snpID, pos, anc, der in self._connection.execute(f"SELECT {SNP_COLUMN_SNP_ID}, {SNP_COLUMN_POSITION}, {SNP_COLUMN_ANCESTRAL}, {SNP_COLUMN_DERIVED} FROM {TABLE_NAME_SNP_ANNOTATION} ORDER BY {SNP_COLUMN_POSITION} ASC"):
			yield snpID, pos, anc, der

	@property
	def SNPsByGenome(self) -> dict[str,Generator[tuple[str,int,str,str],None,None]]:
		out = {}
		for genome_id, genome, _, _, _ in self.references:
			out[genome] = self._getSNPsByGenomeId(genome_id)
		return out

	def SNPByPos(self, pos : int, genome : str=None) -> list[tuple[str]]|str:
		if genome is None:
			return self._connection.execute(f"SELECT {SNP_COLUMN_SNP_ID} FROM {TABLE_NAME_SNP_ANNOTATION} WHERE {SNP_COLUMN_POSITION} = ?", [pos]).fetchall()
		else:
			if (out := self._connection.execute(f"SELECT {SNP_COLUMN_SNP_ID} FROM {TABLE_NAME_SNP_ANNOTATION} WHERE {SNP_COLUMN_POSITION} = ? AND {SNP_COLUMN_GENOME_ID} = ?", [pos, self.genomeID(genome)]).fetchone()) is not None:
				return out[0]
			else:
				raise ValueError(f"No SNP found at {pos=} in {genome=} with genomeID={self.genomeID(genome)}")

	@cached_property
	def SNPsByID(self) -> dict[str,Iterable[tuple[str,int,str,str]]]:
		'''{SNP_ID : [(POSITION, ANCESTRAL, DERIVED), ...]}'''
		return {snpID:self._connection.execute(f"SELECT {SNP_COLUMN_POSITION}, {SNP_COLUMN_ANCESTRAL}, {SNP_COLUMN_DERIVED} FROM {TABLE_NAME_SNP_ANNOTATION} WHERE {SNP_COLUMN_SNP_ID} = ?", [snpID]).fetchone() for (snpID,) in self._connection.execute(f"SELECT {SNP_COLUMN_SNP_ID} FROM {TABLE_NAME_SNP_ANNOTATION}")}
	
	def nodeName(self, nodeID : int) -> dict[str,list[tuple[str,int,str,str]]]:
		for (nodeName,) in self._connection.execute(f"SELECT {NODE_COLUMN_NAME} FROM {TABLE_NAME_NODES} WHERE {NODE_COLUMN_ID} = ?", [nodeID]):
			return nodeName
		raise ValueError(f"No node found for {nodeID=}")
	
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