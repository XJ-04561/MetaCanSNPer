
from SQLOOP.core import *
from MetaCanSNPer.core.LogKeeper import createLogger

SOFTWARE_NAME = "MetaCanSNPer"

import logging, os

LEGACY_HASH = "7630f33662e27489b7bb7b3b121ca4ff"
LEGACY_VERSION = 0

# LOGGER = createLogger(__name__)


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

class NodesTable(Table, name="nodes"): pass
class TreeTable(Table, name="tree"): pass
class ReferencesTable(Table, name="snp_references"): pass
class SNPsTable(Table, name="snp_annotation"): pass

class NodesTable(Table, name="nodes"):
	ID = ID
	Name = Name
class TreeTable(Table, name="tree"):
	Parent = Parent
	Child = Child
	RankI = RankI
class ReferencesTable(Table, name="snp_references"):
	ID = ID
	Genome = Genome
	Strain = Strain
	GenbankID = GenbankID
	RefseqID = RefseqID
	AssemblyName = AssemblyName
class SNPsTable(Table, name="snp_annotation"):
	NodeID = NodeID
	SNPID = SNPID
	Position = Position
	AncestralBase = AncestralBase
	DerivedBase = DerivedBase
	Citation = Citation
	Date = Date
	class GenomeID(Column, name="genome_i"):
		type=INTEGER
