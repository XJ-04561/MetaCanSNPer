
from SQLOOP.core import *
from MetaCanSNPer.Globals import *

LEGACY_HASH = 114303753400556555282992836027662767595
LEGACY_VERSION = 0
LOGGER = LOGGER.getChild(__name__.split(".")[-1])

class LO: # NameSpace
	class Parent(Column):						type=INTEGER
	class Child(Column):						type=INTEGER
	class RankI(Column):						type=INTEGER
	class Name(Column):							type=TEXT
	class ID(Column):							type=INTEGER
	class NodeID(Column):						type=INTEGER
	class SNPID(Column, name="snp_id"):			type=TEXT
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

	class NodeMinus1(Column, name="node_id-1"): pass
	class ParentMinus1(Column, name="parent-1"): pass

	class NodesTable(Table, name="nodes"): pass
	class TreeTable(Table, name="tree"): pass
	class ReferencesTable(Table, name="snp_references"): pass
	class SNPsTable(Table, name="snp_annotation"): pass

	class NodesTable(Table, name="nodes"):
		class ID(Column):							type=INTEGER
		class Name(Column):							type=TEXT
	class TreeTable(Table, name="tree"):
		class Parent(Column):						type=INTEGER
		class Child(Column):						type=INTEGER
		class RankI(Column):						type=INTEGER
		class ParentMinus1(Column, name="parent-1"): pass
		class NodeMinus1(Column, name="node_id-1"): pass
	class ReferencesTable(Table, name="snp_references"):
		class ID(Column):							type=INTEGER
		class Genome(Column):						type=VARCHAR(30)
		class Strain(Column):						type=VARCHAR(30)
		class GenbankID(Column):					type=VARCHAR(30)
		class RefseqID(Column):						type=VARCHAR(30)
		class AssemblyName(Column):					type=VARCHAR(30)
	class SNPsTable(Table, name="snp_annotation"):
		class NodeID(Column):						type=INTEGER
		class SNPID(Column, name="snp_id"):			type=TEXT
		class Position(Column):						type=INTEGER
		class AncestralBase(Column):				type=VARCHAR(1)
		class DerivedBase(Column):					type=VARCHAR(1)
		class Citation(Column, name="reference"):	type=VARCHAR(20)
		class Date(Column):							type=DATETIME
		class GenomeID(Column, name="genome_i"):	type=INTEGER


class NotLegacyCanSNPer2(Assertion):

	LOG = LOGGER

	@classmethod
	def exception(self, database : "MetaCanSNPerDatabase"=None) -> Exception:
		from MetaCanSNPer.modules.Database import IsLegacyCanSNPer2
		return IsLegacyCanSNPer2("Database is Legacy CanSNPer2 schema.")
	@classmethod
	def condition(self, database : "MetaCanSNPerDatabase") -> bool:
		return database.tablesHash != LEGACY_HASH
	@classmethod
	def rectify(self, database : "MetaCanSNPerDatabase") -> None:
		import appdirs
		from json import loads
		from subprocess import check_output as getOutput
		from MetaCanSNPer.modules.Database import (
			TreeTable, NewTreeTable, ReferencesTable,
			NewReferencesTable, ChromosomesTable,
			NewChromosomesTable, SNPsTable, NewSNPsTable,
			Parent, Genotype, NodeID, Position, AncestralBase,
			DerivedBase, Citation, Date, ChromosomeID,
			Chromosome, GenomeID, Genome, Strain, GenbankID,
			RefseqID, AssemblyName
		)
		
		refDir = PathGroup(appdirs.user_data_dir(SOFTWARE_NAME), *appdirs.site_data_dir(SOFTWARE_NAME, multipath=True).split(os.pathsep)) / "References" / database.organism
		
		database(BEGIN - TRANSACTION)

		# References
		self.LOG.info("Updating 'References'-table")
		database(CREATE - TABLE - sql(NewReferencesTable))
		database(INSERT - INTO (NewReferencesTable) - SELECT(ALL) - FROM (LO.ReferencesTable))

		# Chromosomes
		self.LOG.info("Updating 'Chromosomes'-table")
		database(CREATE - TABLE - sql(ChromosomesTable) )
		j = 0
		ref2chromLookup = {}
		for i, genbankID, assembly in database(SELECT (LO.ReferencesTable.ID, LO.GenbankID, LO.AssemblyName) - FROM (LO.ReferencesTable)):
			ref2chromLookup[i] = []
			chromosomes = ()
			assemblyFile = refDir.find(f"{assembly}.fna") or Path("?")
			commandName = f"datasets{'.exe' if os.name == 'nt' else ''}"
			
			if shutil.which(commandName):
				chromosomes = tuple(map(*this["value"].strip("\"'"), loads(getOutput(f"{commandName} summary genome accession {genbankID} --as-json-lines".split()))["assembly_info"]["biosample"]["sample_ids"]))
			
			if len(chromosomes) > 0:
				pass # genbank entry found
			elif assemblyFile.exists:
				# No genbank entry found
				chromosomes = map(*this[1:].split()[0], filter(*this.startswith(">"), open(assemblyFile, "r").readline()))
			else:
				self.LOG.exception(UnableToDefineChromosomes(f"Could not find naming for chromosomes in entry with {i=}, {genbankID=}, and {assembly=}."))
				self.LOG.warning(f"Couldn't find genome with {genbankID=} and {assembly=} either online or in {refDir}.")
				chromosomes = (NULL,)

			for chromosome in chromosomes:
				database(INSERT - OR - REPLACE - INTO (ChromosomesTable) - (ChromosomeID, Chromosome, GenomeID) - VALUES (j, chromosome, i))
				ref2chromLookup[i].append(j)
				j += 1
				break # FOR NOW
		
		# SNPs
		# TODO : Get the real chromosome ID based on the position of the SNP and the lengths of the chromosomes
		self.LOG.info("Updating 'SNP'-table")
		database(CREATE - TABLE - sql(NewSNPsTable))
		database(INSERT - INTO (NewSNPsTable) - (NodeID, Position, AncestralBase, DerivedBase, Citation, Date, ChromosomeID) - SELECT (LO.NodeMinus1, LO.Position, LO.AncestralBase, LO.DerivedBase, LO.Citation, LO.Date, LO.SNPsTable.GenomeID) - FROM (LO.SNPsTable))

		
		# Tree
		self.LOG.info("Updating 'Tree'-table")
		database(CREATE - TABLE - sql(NewTreeTable))
		database(INSERT - INTO (NewTreeTable) - (Parent, NodeID, Genotype) - SELECT (LO.TreeTable.ParentMinus1, NULL, LO.NodesTable.Name) - FROM (LO.TreeTable, LO.NodesTable) - WHERE (LO.TreeTable.Child == LO.NodesTable.ID, LO.TreeTable.Child > 1) - ORDER - BY (LO.TreeTable.Child - ASC))
		database(UPDATE (NewTreeTable) - SET (parent = 0) - WHERE (NodeID == 1))


		class Genomes(Table): pass
		class Rank(Table): pass
		database(DROP - TABLE (LO.ReferencesTable))
		database(DROP - TABLE (LO.TreeTable))
		database(DROP - TABLE (LO.NodesTable))
		database(DROP - TABLE (LO.SNPsTable))
		database(DROP - TABLE (Genomes))
		database(DROP - TABLE (Rank))
		
		database(ALTER - TABLE (NewTreeTable) - RENAME - TO (TreeTable))
		database(ALTER - TABLE (NewSNPsTable) - RENAME - TO (SNPsTable))
		database(ALTER - TABLE (NewReferencesTable) - RENAME - TO (ReferencesTable))

		database(DROP - TABLE - IF - EXISTS (NewTreeTable))
		database(DROP - TABLE - IF - EXISTS (NewSNPsTable))
		database(DROP - TABLE - IF - EXISTS (NewReferencesTable))

		for index in database.indexes:
			database.createIndex(index)

		database(PRAGMA (user_version = database.CURRENT_VERSION))
		database(COMMIT)

try:
	from MetaCanSNPer.modules.Database import MetaCanSNPerDatabase, CanSNPDatabaseError, IsLegacyCanSNPer2
except:
	pass