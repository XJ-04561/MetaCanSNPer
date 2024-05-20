

import SQLOOP.Globals as Globals
from SQLOOP import *
from SQLOOP.core import *
import argparse, sys

from MetaCanSNPer.Globals import *
from MetaCanSNPer.modules.Downloader import DatabaseDownloader, Downloader, URL, gunzip, DownloaderReportHook
from MetaCanSNPer.core.LogKeeper import createLogger

LOGGER = LOGGER.getChild("Database")

SOFTWARE_NAME = "MetaCanSNPer"

import logging, os

LEGACY_HASH = "7630f33662e27489b7bb7b3b121ca4ff"
LEGACY_VERSION = 0


class Parent(Column):						type=INTEGER
class GenoType(Column):						type=TEXT
class NodeID(Column):						type=INTEGER
class Position(Column):						type=INTEGER
class AncestralBase(Column):				type=VARCHAR(1)
class DerivedBase(Column):					type=VARCHAR(1)
class Citation(Column, name="reference"):	type=VARCHAR(20)
class Date(Column):							type=DATETIME
class ChromosomeID(Column):					type=INTEGER
class Chromosome(Column):					type=VARCHAR(30)
class GenomeID(Column):						type=INTEGER
class Genome(Column):						type=VARCHAR(30)
class Strain(Column):						type=VARCHAR(30)
class GenbankID(Column):					type=VARCHAR(30)
class RefseqID(Column):						type=VARCHAR(30)
class AssemblyName(Column):					type=VARCHAR(30)

class TreeTable(Table, name="tree"): pass
class ChromosomesTable(Table, name="chromosomes"): pass
class ReferencesTable(Table, name="snp_references"): pass
class SNPsTable(Table, name="snp_annotation"): pass

class TreeTable(Table, name="tree"):
	Parent = Parent
	Child = NodeID
	GenoType = GenoType
	options = (
		PRIMARY - KEY (Child),
		UNIQUE (Parent, Child),
		UNIQUE (Child),
		NOT - NULL (Parent),
		NOT - NULL (GenoType)
	)
class ReferencesTable(Table, name="snp_references"):
	GenomeID = GenomeID
	Genome = Genome
	Strain = Strain
	GenbankID = GenbankID
	RefseqID = RefseqID
	AssemblyName = AssemblyName
	options = (
		PRIMARY - KEY (ChromosomeID, Position, NodeID),
		FOREIGN - KEY (ChromosomeID) - REFERENCES (ChromosomesTable, ChromosomeID),
		FOREIGN - KEY (NodeID) - REFERENCES (TreeTable, NodeID),
		UNIQUE (Position)
	)
class ChromosomesTable(Table, name="chromosomes"):
	ChromosomeID = ChromosomeID
	Chromosome = Chromosome
	GenomeID = GenomeID
	options = (
		PRIMARY - KEY (ChromosomeID, GenomeID),
		FOREIGN - KEY (GenomeID) - REFERENCES (ReferencesTable, GenomeID),
		UNIQUE (ChromosomeID),
		NOT - NULL (Chromosome),
		NOT - NULL (GenomeID)
	)
class SNPsTable(Table, name="snp_annotation"):
	NodeID = NodeID
	Position = Position
	AncestralBase = AncestralBase
	DerivedBase = DerivedBase
	Citation = Citation
	Date = Date
	ChromosomeID = ChromosomeID
	options = (
		PRIMARY - KEY (GenomeID),
		UNIQUE (GenbankID),
		UNIQUE (RefseqID),
		UNIQUE (AssemblyName)
	)

"""Indexes"""

class TreeTableByParent(Index):
	table = TreeTable
	Parent = Parent

class SNPsTableByPosition(Index):
	table = SNPsTable
	Position = Position

class SNPsTableByNodeID(Index):
	table = SNPsTable
	NodeID = NodeID

class SNPsTableByChromID(Index):
	table = SNPsTable
	ChromosomeID = ChromosomeID

class ReferencesTableByGenome(Index):
	table = ReferencesTable
	Genome = Genome

class ReferencesTableByAssembly(Index):
	table = ReferencesTable
	AssemblyName = AssemblyName

class CanSNPDatabaseError(DatabaseError): pass
class IsLegacyCanSNPer2(CanSNPDatabaseError): pass
class NoTreeConnectedToRoot(CanSNPDatabaseError): pass

import  MetaCanSNPer.modules.LegacyDatabase as Legacy

class NotLegacyCanSNPer2(Assertion):
	legacyObjects = {}
	@classmethod
	def exception(self, database : "MetaCanSNPerDatabase"=None) -> Exception:
		return IsLegacyCanSNPer2("Database is Legacy CanSNPer2 schema.")
	@classmethod
	def condition(self, database : "MetaCanSNPerDatabase") -> bool:
		return database.tablesHash != LEGACY_HASH
	@classmethod
	def rectify(self, database : "MetaCanSNPerDatabase") -> None:
		import appdirs
		import LegacyDatabase as LO
		from json import loads
		from subprocess import check_output as getOutput

		class NodeMinus1(Column, name="node_id-1"): pass
		class ParentMinus1(Column, name="parent-1"): pass

		class TempTable(Table):
			
			class NodeMinus1(Column, name="node_id-1"): pass
			class ParentMinus1(Column, name="parent-1"): pass
			Child = LO.Child
		
		refDir = PathGroup(appdirs.user_data_dir(SOFTWARE_NAME), appdirs.site_data_dir(SOFTWARE_NAME)) / "References" / database.organism
		
		# References
		LOGGER.info("Updating 'References'-table")
		database(BEGIN - TRANSACTION)
		database(ALTER - TABLE (LO.ReferencesTable) - RENAME - TO - TempTable)
		database(COMMIT)
		database(BEGIN - TRANSACTION)
		database(CREATE - TABLE - sql(ReferencesTable))
		database(INSERT - INTO - ReferencesTable - SELECT * FROM (TempTable))
		database(DROP - TABLE (TempTable))
		database(COMMIT)

		# Chromosomes
		LOGGER.info("Updating 'Chromosomes'-table")
		database(BEGIN - TRANSACTION)
		database(CREATE - TABLE - sql(ChromosomesTable) )
		for i, GenbankID, assembly in database(SELECT (LO.GenomeID, LO.GenbankID, LO.AssemblyName) - FROM (LO.ReferencesTable)):
			try:
				chromosome = loads(getOutput("".split()))["assemblies"][0]["assembly"]["biosample"]["sample_ids"][0]["value"]
				assert len(chromosome) != 0, "No genbank entry found"
				database(INSERT - INTO - ChromosomesTable - VALUES (i, chromosome, i))
			except (AssertionError, KeyError) as e:
				LOGGER.exception(e)
				try:
					chromosome = open(refDir.find(f"{assembly}.fna"), "r").readline()[1:].split()[0]
				except FileNotFoundError:
					LOGGER.warning(f"Couldn't find genome with {GenbankID=} either online or in {refDir!r}.")
					raise UnableToDefineChromosomes(f"Could not find naming for chromosomes in entry with {i=}, {GenbankID=}, and {assembly=}.")
			finally:
				database(INSERT - INTO - ChromosomesTable - VALUES (i, chromosome, i))
					
		database(COMMIT)
		

		# SNPs
		LOGGER.info("Updating 'SNP'-table")
		database(BEGIN - TRANSACTION)
		database(ALTER - TABLE - SNPsTable - RENAME - TO - TempTable)
		database(COMMIT)
		database(BEGIN - TRANSACTION)
		database(CREATE - TABLE (sql(SNPsTable)))
		database(INSERT - INTO - SNPsTable - (NodeID, Position, AncestralBase, DerivedBase, Citation, Date, ChromosomeID) - SELECT (NodeMinus1, LO.Position, LO.AncestralBase, LO.DerivedBase, LO.Citation, LO.Date, LO.GenomeID) - FROM (TempTable))
		database(DROP - TABLE - TempTable)
		database(COMMIT)

		
		# Tree
		LOGGER.info("Updating 'Tree'-table")
		database(BEGIN - TRANSACTION)
		database(ALTER - TABLE - TreeTable - RENAME - TO - TempTable)
		database(COMMIT)
		database(BEGIN - TRANSACTION)
		database(CREATE - TABLE - sql(TreeTable))
		database(COMMIT)
		database(BEGIN - TRANSACTION)
		database(INSERT - INTO - TreeTable - (Parent, NodeID, GenoType) - SELECT (TempTable.ParentMinus1, NULL, LO.NodesTable.Name) - FROM (TempTable, LO.NodesTable) - WHERE (TempTable.Child == LO.NodesTable.ID, TempTable.Child > 1) - ORDER - BY - (TempTable.Child - ASC))
		database(COMMIT)
		database(BEGIN - TRANSACTION)
		database(UPDATE - TreeTable - SET (Parent == 0) - WHERE (NodeID == 1))
		database(DROP - TABLE - LO.NodesTable)
		database(DROP - TABLE - TempTable)

		LOGGER.info("Dropping 'genomes'- and 'rank'-tables")

		class Genomes(Table): pass
		class Rank(Table): pass

		database(DROP - TABLE - Genomes)
		database(DROP - TABLE - Rank)
		
		for index in database.indexes:
			database.createIndex(index)

		database(PRAGMA (user_version = database.CURRENT_VERSION))
		database(COMMIT)

class CanSNPNode(Branch):
	table : Table = TreeTable
	parentCol : Column = Parent
	childCol : Column = NodeID

class MetaCanSNPerDatabase(Database):
	assertions = (NotLegacyCanSNPer2, *Globals.ASSERTIONS)

	LOGGER = createLogger(__name__)
	DATABASE_VERSIONS = {
		LEGACY_HASH	: LEGACY_VERSION, # Legacy CanSNPer
	}

	CURRENT_VERSION = 2

	organism : str

	def __init__(self, filename: str, mode: Globals.Mode, organism : str=None):
		self.organism = organism or pName(filename)
		super().__init__(filename, mode)

	@property
	def tree(self):
		try:
			return CanSNPNode(self, next(self()))
		except StopIteration:
			raise NoTreeConnectedToRoot(f"In database {self.filename!r}")
	
	@cached_property
	def references(self):
		return tuple(self[ALL, ReferencesTable])
	
	@property
	def SNPs(self):
		return self[ALL, SNPsTable]
	
	@cached_property
	def chromosomes(self):
		return tuple(self[ALL, ChromosomesTable])
	
	@cached_property
	def SNPsByReference(self):
		return {genome:tuple(self[ALL, SNPsTable, GenomeID == genomeID]) for genomeID, genome in self[GenomeID, Genome]}
	
	@cached_property
	def SNPsByNode(self):
		return {nodeID:tuple(self[ALL, SNPsTable, NodeID == nodeID]) for nodeID in self[NodeID, TreeTable]}

def correctDatabase(filename):
	database = MetaCanSNPerDatabase(filename, "w")
	database.fix()
	if not database.valid:
		e = database.exception
		e.add_note(f"Tables hash: {database.tablesHash}\nIndexes Hash: {database.indexesHash}")
		raise e
	database.close()

DatabaseDownloader.postProcess = correctDatabase

def loadFromReferenceFile(database : Database, file : TextIO, refDir : str="."):
	file.seek(0)
	if "genome	strain	genbank_id	refseq_id	assembly_name" == file.readline():
		for row in file:
			genome, strain, genbank_id, refseq_id, assembly_name = row.strip().split("\t")
			database(INSERT - INTO - ReferencesTable - VALUES (genome, strain, genbank_id, refseq_id, assembly_name))
			GenomeID = next(database(SELECT (GenomeID) - FROM (ReferencesTable) - WHERE (Genome == genome)))
			try:
				chrom = open(os.path.join(refDir, f"{assembly_name}.fna"), "r").readline()[1:].split()[0]
				database(INSERT - INTO - ChromosomesTable - VALUES (None, chrom, GenomeID))
			except FileNotFoundError as e:
				raise MissingReferenceFile(f"Could not find reference file {os.path.join(refDir, f'{assembly_name}.fna')!r}. The file {file.__name__!r} does not specify chromosomes, and so the reference fasta file is required. To set the directory in which to look for .fna references, use the flag '--refDir'")

	elif "chromosomes	genome	strain	genbank_id	refseq_id	assembly_name" == file.readline():
		for row in file:
			chromosomes, genome, strain, genbank_id, refseq_id, assembly_name = row.strip().split("\t")
			database(INSERT - INTO - ReferencesTable - VALUES (genome, strain, genbank_id, refseq_id, assembly_name))
			GenomeID = next(database(SELECT (GenomeID) - FROM (ReferencesTable) - WHERE (Genome == genome)))
			for chrom in chromosomes.split(";"):
				database(INSERT - INTO - ChromosomesTable - VALUES (None, chrom, GenomeID))
	else:
		ValueError("File is not of accepted format.")

def loadFromTreeFile(database : Database, file : TextIO):
	file.seek(0)
	database(INSERT - INTO - TreeTable - VALUES (0, None, file.readline().strip()))
	for row in file:
		*_, parent, child = row.rstrip(file.newlines).rsplit("\t", 2)
		database(INSERT - INTO - TreeTable - VALUES (SELECT (NodeID) - FROM (TreeTable) - WHERE (GenoType == parent), None, child))

def loadFromSNPFile(database : Database, file : TextIO):
	file.seek(0)
	if "snp_id	strain	reference	genome	position	derived_base	ancestral_base" == file.readline().strip():
		for row in file:
			nodeName, strain, reference, genome, position, ancestral, derived = row.rstrip(file.newlines).split("\t")
			database(INSERT - INTO - SNPsTable - VALUES (SELECT (NodeID) - FROM (TreeTable) - WHERE (GenoType == nodeName), position, ancestral, derived, reference, None, SELECT (ChromosomeID) - FROM (ChromosomesTable) - WHERE (Genome == genome)))
	else:
		ValueError("File is not of accepted format.")

def main():
	def write(databasePath : str=None, SNPFile : str=None, treeFile : str=None, referenceFile : str=None, **kwargs):

		LOGGER.debug(f"{databasePath=}")
		database : MetaCanSNPerDatabase = MetaCanSNPerDatabase(databasePath, "w")

		database.checkDatabase(mode = "r")
		
		print(database)
		
		if SNPFile is not None:
			loadFromReferenceFile(database, referenceFile)
			loadFromTreeFile(database, treeFile)
			loadFromSNPFile(database, SNPFile)
		else:
			LOGGER.warning("No files given to build database from. Created an empty database with the current MetaCanSNPer structure.")

		database.close()

	def update(databasePaths : list[str]=None, refDir : str=".", noCopy=False, **kwargs):

		LOGGER.debug(f"{databasePaths=}")
		for databasePath in databasePaths:
			database : MetaCanSNPerDatabase = MetaCanSNPerDatabase(databasePath, "w")

			database.fix()
			
			if database.valid:
				print(f"Updated {databasePath} succesfully!")
				print(database)
			else:
				raise database.exception

			database.close()

	def download(databaseNames : list[str]=[], outDir : str=Path("."), **kwargs):
		
		from MetaCanSNPer.modules.Downloader import DatabaseDownloader
		from MetaCanSNPer.core.TerminalUpdater import TerminalUpdater
		from MetaCanSNPer.core.Hooks import DownloaderReportHook
		LOGGER.debug(f"{databaseNames=}")
		out = []
		TU = TerminalUpdater("Downloading databases")
		DD = DatabaseDownloader(outDir)
		for databaseName in databaseNames:
			RH = DownloaderReportHook("DatabaseDownloader", TU.hooks, databaseName)
			DD.download(databaseName, outDir / databaseName, reportHook=RH)
				
			#raise DownloadFailed(f"Failed to download {databaseName} to {os.path.join(outDir, databaseName)}.")
			out.append(outDir / databaseName)
		DD.wait()
		TU.stop()
		return out

	def test(database : list[Path]= [], refDir : Path=".", outDir : Path=".", noCopy : bool=False, **kwargs):

		LOGGER.debug(f"{kwargs['database']=}")
		print("Testing Download:")
		databasePaths : list[str] = download(databaseNames=kwargs['database'], outDir=outDir)

		print("Testing Update:")
		update(databasePaths=databasePaths, refDir=refDir, noCopy=noCopy)

		print("Testing Read:")
		for databasePath in databasePaths:
			print(f"  {databasePath.replace(os.path.realpath('.'), '.').replace(os.path.expanduser('~'), '~')}")
			LOGGER.debug(f"{databasePath.replace(os.path.realpath('.'), '.').replace(os.path.expanduser('~'), '~')}")

			database = MetaCanSNPerDatabase(databasePath, "r")
			LOGGER.debug(repr(database))

			print(f"    Arbitrary `.get` from one table only.")
			LOGGER.debug(f"Arbitrary `.get` from one table only.")
			string = []
			for row in database[Position, AncestralBase, DerivedBase, ChromosomeID == 1]:
				string.append(",\t".join(map(str, row)))
			LOGGER.debug("\n".join(string))

			print(f"    Get one entry from a table.")
			LOGGER.debug(f"Get one entry from a table.")
			
			chromID, chromosome, GenomeID = next(iter(database[ChromosomesTable]))
			LOGGER.debug(f"{chromID=}, {chromosome=}, {GenomeID=}")

			GenomeID, genome, strain, genbank, refseq, assembly = next(iter(database[ReferencesTable]))
			LOGGER.debug(f"{GenomeID=}, {genome=}, {strain=}, {genbank=}, {refseq=}, {assembly=}")

			print(f"    Arbitrary `.get` from one table referencing adjacent table.")
			LOGGER.debug(f"Arbitrary `.get` from one table referencing adjacent table.")
			string = []
			for row in database[Position, AncestralBase, DerivedBase, Chromosome == chromosome]:
				string.append(",\t".join(map(str, row)))
			LOGGER.debug("\n".join(string))

			print(f"    Arbitrary `.get` from one table referencing across more than one chained table.")
			LOGGER.debug(f"Arbitrary `.get` from one table referencing across more than one chained table.")
			string = []
			for row in database[Position, AncestralBase, DerivedBase, Genome == genome]:
				string.append(",\t".join(map(str, row)))
			LOGGER.debug("\n".join(string))
	
	###############################################################################################################

	parser = argparse.ArgumentParser(prog="MetaCanSNPer")

	modeGroup : argparse._SubParsersAction = parser.add_subparsers(title="Mode", dest="mode", description="Mode with which to open the database.", metavar="MODES")

	writeParser : argparse.ArgumentParser = modeGroup.add_parser("write",	help="Create a database with or without data. Data for database is given through the appropriate File flags.")
	writeParser.add_argument("databasePath", type=os.path.realpath)
	filesGroup = writeParser.add_argument_group(title="Input Files")
	
	filesGroup.add_argument("--SNPFile",		type=os.path.realpath,	help="If used, make sure that the related references and tree nodes are present in the database or in the other flagged files.")
	filesGroup.add_argument("--referenceFile",	type=os.path.realpath)
	filesGroup.add_argument("--treeFile",		type=os.path.realpath)

	writeParser.add_argument("--refDir", help="Directory where the reference genomes are located. This is only required if your --referenceFile doesn't have a `chromosomes` column.")

	optionalGroup = writeParser.add_argument_group(title="Optional Flags")
	optionalGroup.add_argument("--rectify",	action="store_true", help="If used, will edit the database structure if it doesn't comply with the current set schema. If not used, will continue operations without rectifying, but the program might crash due to the difference in schema.")

	writeParser.set_defaults(func=write)

	updateParser : argparse.ArgumentParser = modeGroup.add_parser("update", help="Update an existing database to follow the current standard schema.")
	updateParser.add_argument("databasePaths", nargs="+", type=os.path.realpath)
	updateParser.add_argument("--refDir")
	updateParser.add_argument("--noCopy", action="store_true")
	updateParser.set_defaults(func=update)
	
	downloadParser : argparse.ArgumentParser = modeGroup.add_parser("download", help="Download a database from one of the internally defined sources.")
	downloadParser.add_argument("databaseNames", nargs="+", type=os.path.basename)
	downloadParser.add_argument("--outDir", default=os.path.realpath("."))
	downloadParser.set_defaults(func=download)

	testParser : argparse.ArgumentParser = modeGroup.add_parser("test", help="Test out the features of MetaCanSNPerDatabase to see if your environment is suitable for using it.")
	testParser.add_argument("database", nargs="+", type=os.path.realpath, default=["francisella_tularensis.db"])
	testParser.add_argument("--refDir")
	testParser.add_argument("--noCopy", action="store_true")
	testParser.add_argument("--outDir", default=os.path.realpath(".")) # CommonGroups.shared / "MetaCanSNPer-Data" / "Databases")
	testParser.set_defaults(func=test)

	parser.add_argument("--version", action="store_true")
	parser.add_argument("--debug", action="store_true")
	parser.add_argument("--info", action="store_true")
	parser.add_argument("--noLog", action="store_true")

	if len(sys.argv) <= 1:
		parser.print_help()
		exit(0)
	elif "--version" in sys.argv:
		print(f"MetaCanSNPerDatabase v. {MetaCanSNPerDatabase.CURRENT_VERSION}")
		exit(0)
	elif all(mode not in sys.argv for mode in ["write", "update", "download", "test"]):
		print("No mode chosen check usage to see which mode is appropriate for your intended use.", file=sys.stderr)
		parser.print_help()
		exit(1)

	args : argparse.Namespace = parser.parse_args(sys.argv[2:])
	
	if args.debug:
		LOGGER.setLevel(logging.DEBUG)
	elif args.info:
		LOGGER.setLevel(logging.INFO)
	if args.noLog:
		LOGGER.disabled = True

	try:
		args.func(**dict(args._get_kwargs()))
	except Exception as e:
		LOGGER.exception(e)
		print(f"{type(e).__name__}:", e, file=sys.stderr)
		exit(1)

	print("Done!")

	exit(0)

