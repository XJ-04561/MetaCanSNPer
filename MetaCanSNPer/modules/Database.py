

import SQLOOP.Globals as Globals
from SQLOOP import *
from SQLOOP.core import *
import argparse, sys

from MetaCanSNPer.Globals import *
from MetaCanSNPer.modules.Downloader import DatabaseDownloader, DatabaseThread

LOGGER = LOGGER.getChild(__name__.split(".")[-1])


class Parent(Column):						type=INTEGER
class Genotype(Column):						type=TEXT
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

class TreeTable(Table, name="tree"):
	Parent = Parent
	Child = NodeID
	Genotype = Genotype
	constraints = (
		PRIMARY - KEY (Child),
	)
NewTreeTable = type(Table)("NewTreeTable", (Table,), dict(vars(TreeTable)))
class ReferencesTable(Table, name="snp_references"):
	GenomeID = GenomeID
	Genome = Genome
	Strain = Strain
	GenbankID = GenbankID
	RefseqID = RefseqID
	AssemblyName = AssemblyName
	constraints = (
		PRIMARY - KEY (GenomeID),
		UNIQUE (GenbankID),
		UNIQUE (RefseqID),
		UNIQUE (AssemblyName)
	)
NewReferencesTable = type(Table)("NewReferencesTable", (Table,), dict(vars(ReferencesTable)))
class ChromosomesTable(Table, name="chromosomes"):
	ChromosomeID = ChromosomeID
	Chromosome = Chromosome
	GenomeID = GenomeID
	constraints = (
		PRIMARY - KEY (ChromosomeID),
		# FOREIGN - KEY (GenomeID) - REFERENCES (ReferencesTable, GenomeID)
	)
NewChromosomesTable = type(Table)("NewChromosomesTable", (Table,), dict(vars(ChromosomesTable)))
class SNPsTable(Table, name="snp_annotation"):
	NodeID = NodeID
	Position = Position
	AncestralBase = AncestralBase
	DerivedBase = DerivedBase
	Citation = Citation
	Date = Date
	ChromosomeID = ChromosomeID
	constraints = (
		PRIMARY - KEY (Position),
		# FOREIGN - KEY (ChromosomeID) - REFERENCES (ChromosomesTable, ChromosomeID),
		# FOREIGN - KEY (NodeID) - REFERENCES (TreeTable, NodeID)
	)
NewSNPsTable = type(Table)("NewSNPsTable", (Table,), dict(vars(SNPsTable)))

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
class NoTreeConnectedToRoot(CanSNPDatabaseError): pass
class IsLegacyCanSNPer2(CanSNPDatabaseError): pass
class NoChromosomesInDatabase(CanSNPDatabaseError): pass

class CanSNPNode(Branch):
	table : Table = TreeTable
	parentCol : Column = Parent
	childCol : Column = NodeID

from MetaCanSNPer.modules.LegacyDatabase import NotLegacyCanSNPer2, LEGACY_HASH, LEGACY_VERSION

class HasChromosomes(Assertion, Logged):

	@classmethod
	def exception(self, database : "MetaCanSNPerDatabase"=None) -> Exception:
		return NoChromosomesInDatabase("Database doesn't have chromosomes defined.")
	@classmethod
	def condition(self, database : "MetaCanSNPerDatabase") -> bool:
		try:
			return database(SELECT (COUNT(ALL)) - FROM (ChromosomesTable) - WHERE (Chromosome - IS (NULL))) == 0
		except:
			return False
	@classmethod
	def rectify(self, database : "MetaCanSNPerDatabase") -> None:
		import appdirs
		from json import loads
		from subprocess import check_output as getOutput
		
		refDir = PathGroup([appdirs.user_data_dir(SOFTWARE_NAME), *appdirs.site_data_dir(SOFTWARE_NAME, multipath=True).split(os.pathsep)]) / "References" / database.organism
		
		database(BEGIN - TRANSACTION)

		# Chromosomes
		self.LOG.info("Updating 'Chromosomes'-table")
		database(DELETE - FROM (ChromosomesTable) )

		highestGenomeID = max(database[GenomeID, ReferencesTable])
		j = highestGenomeID + 1
		ref2chromLookup = {}
		for i, genbankID, assembly in database[GenomeID, GenbankID, AssemblyName, ReferencesTable]:
			ref2chromLookup[i] = []
			assemblyFile = refDir.find(f"{assembly}.fna")
			
			### DATASETS SUMMARY DOESN'T HAVE THE .FNA CHROMOSOME NAMES
			# if shutil.which(commandName):
			# 	chromosomes = tuple(map(*this["value"].strip("\"'"), loads(getOutput(f"{commandName} summary genome accession {genbankID} --as-json-lines".split()))["assembly_info"]["biosample"]["sample_ids"]))
			
			if assemblyFile:
				# No genbank entry found
				with open(assemblyFile, "r") as refFile:
					chromosomes = [["N/A", 0]]
					_iter = iter(refFile)
					for row in _iter:
						if row.startswith(">"):
							chromosomes.append([row[1:].split()[0], 0])
						else:
							chromosomes[-1][1] += len(row.strip())
			else:
				self.LOG.error(f"Couldn't find genome with {genbankID=} either online or in {refDir}.")
				raise UnableToDefineChromosomes(f"Can't find fasta file for reference genome with assembly name: {assembly}")

			for [chromosome, length], prevLength in zip(chromosomes[1:], itertools.accumulate(chromosomes[1:-1], lambda x, y:x+y[1], initial=0)):
				database(INSERT - OR - REPLACE - INTO (ChromosomesTable) - (ChromosomeID, Chromosome, GenomeID) - VALUES (j, chromosome, i))
				database(UPDATE (SNPsTable) - SET (chromosome_id = j, position = (Position - prevLength) ) - WHERE (ChromosomeID == i, Position > prevLength, Position - prevLength <= length))
				j += 1
		database(UPDATE (SNPsTable) - SET (chromosome_id = ChromosomeID - highestGenomeID))
		database(UPDATE (ChromosomesTable) - SET (chromosome_id = ChromosomeID - highestGenomeID))
				
		database(COMMIT)

class MetaCanSNPerDatabase(Database, Logged):

	CURRENT_VERSION = 2
	DATABASE_VERSIONS = {
		LEGACY_HASH	: LEGACY_VERSION, # Legacy CanSNPer
	}
	assertions = (NotLegacyCanSNPer2, *Globals.ASSERTIONS, HasChromosomes)
	organism : str

	TreeTable = TreeTable
	ReferencesTable = ReferencesTable
	ChromosomesTable = ChromosomesTable
	SNPsTable = SNPsTable

	TreeTableByParent = TreeTableByParent
	SNPsTableByPosition = SNPsTableByPosition
	SNPsTableByNodeID = SNPsTableByNodeID
	SNPsTableByChromID = SNPsTableByChromID
	ReferencesTableByGenome = ReferencesTableByGenome
	ReferencesTableByAssembly = ReferencesTableByAssembly

	def __init__(self, filename: str, mode: Globals.Mode, organism : str=None):
		self.organism = organism or pName(filename)
		super().__init__(filename, mode, factoryFunc=DatabaseThread)

	@property
	def tree(self):
		try:
			return CanSNPNode(self, first(self[NodeID, TreeTable, Parent == 0]))
		except StopIteration:
			raise NoTreeConnectedToRoot(f"In database {self.filename!r}")
	
	@cached_property
	def references(self):
		return tuple(entry for entry in self[ALL, ReferencesTable])
	
	@property
	def SNPs(self):
		return self[ALL, SNPsTable]
	
	@cached_property
	def chromosomes(self):
		return tuple(entry for entry in self[ALL, ChromosomesTable])
	
	@cached_property
	def SNPsByReference(self):
		return {genome:tuple(entry for entry in self[ALL, SNPsTable, GenomeID == genomeID]) for genomeID, genome in self[GenomeID, Genome]}
	
	@cached_property
	def SNPsByNode(self):
		return {nodeID:tuple(entry for entry in self[ALL, SNPsTable, NodeID == nodeID]) for nodeID in self[NodeID, TreeTable]}

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
		database(INSERT - INTO - TreeTable - VALUES (SELECT (NodeID) - FROM (TreeTable) - WHERE (Genotype == parent), None, child))

def loadFromSNPFile(database : Database, file : TextIO):
	file.seek(0)
	if "snp_id	strain	reference	genome	position	derived_base	ancestral_base" == file.readline().strip():
		for row in file:
			nodeName, strain, reference, genome, position, ancestral, derived = row.rstrip(file.newlines).split("\t")
			database(INSERT - INTO - SNPsTable - VALUES (SELECT (NodeID) - FROM (TreeTable) - WHERE (Genotype == nodeName), position, ancestral, derived, reference, None, SELECT (ChromosomeID) - FROM (ChromosomesTable) - WHERE (Genome == genome)))
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

