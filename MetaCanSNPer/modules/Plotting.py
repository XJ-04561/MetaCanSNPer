

from GeekyGadgets.Illustrative.Graphs import Tree, Branch, Leaf, TreeRoot
from GeekyGadgets.Illustrative.Plots import LinePlot, Histogram, BoxPlot
from GeekyGadgets.Illustrative.Tables import Table
from GeekyGadgets.Illustrative import IllustrativeCollection, IllustrateFigure, IllustrateDocument

from PseudoPathy import FilePath
from GeekyGadgets.Globals import NAN
from GeekyGadgets.Classy import CachedDefault, Default
from GeekyGadgets.Semantics.Markups.HTML import I as Italic, H1 as Head1, P as Paragraph
from GeekyGadgets.Colors import HSV, ColorScale, Color
from GeekyGadgets.TypeHinting import Any, overload, TextIO, Number
from GeekyGadgets.Formatting.SISize import percentEmphasize
from GeekyGadgets.Formatting.Case import TitleCase, SentenceCase
from GeekyGadgets.Formatting.Alignment import alignSignature
from GeekyGadgets.Iterators.Walkers import LeafClimber, BranchWalker
from GeekyGadgets.Math.Stats.Clustering import combineClose, combineDensity
import GeekyGadgets.Iterators as _Iterators
from GeekyGadgets.This import this
from GeekyGadgets.Functions import first
from GeekyGadgets.Logging import ROOT_LOGGER
from math import log10, sqrt, isnan
from functools import cached_property, cache

LOGGER = ROOT_LOGGER.getChild(__name__)

metaCanSNPColorScale = ColorScale(
	HSV(96, 100, 100),
	HSV(215, 100, 100),
	HSV(340, 100, 100),
	keyframes=(0.9, 0.5, 0.1)
)
	
__all__ = (
	"CanSNPTree", "treeGraph", "coverageHistogram", "coverageBoxplot", "variantsTable",
	"nonCanonTable", "badCoverageTable", "interactiveFigure", "figuresDocument"
)

class CanSNPLeaf(Leaf):
	
	name = property(lambda self: getattr(self, "genotype", None))
	properties = [("genotype", "string")]

	@property
	def attributes(self):
		return {
			"coverage" : self.coverage,
			"called" : self.called,
			"nonCanon" : self.nonCanon
		}

	@cached_property
	def hidden(self) -> bool:
		
		if not self.incoming:
			return False
		parent = self.incoming[0].incoming
		if (
			parent.called
			and isinstance(parent.incoming[0].logRatio, Number)
			and 0.5 <= parent.incoming[0].logRatio
		):
			return False
		
		if any(not child.hidden for child in self.children):
			return False
		if isinstance(self.incoming[0].logRatio, Number) and 0.5 <= self.incoming[0].logRatio:
			return False

		return True

	@property
	def coverage(self) -> int|None:
		if not self.incoming:
			return None
		elif not isinstance(self.incoming[0].depth, Number):
			return None
		return self.incoming[0].depth

	@property
	def called(self) -> bool:
		if not self.incoming:
			return False
		elif isinstance(self.incoming[0].logRatio, Number) and self.incoming[0].logRatio >= 0.5:
			return True
		elif not (parent := self.incoming[0].incoming).incoming:
			return False
		
		if not parent.called:
			return False
		
		called = self.incoming[0].called
		prevCalled = parent.incoming[0].called
		if not isinstance(called, Number) or called == 0:
			return False
		elif not isinstance(prevCalled, Number) or prevCalled == 0:
			return False
		elif log10(prevCalled) == 0:
			return False
		elif log10(called) / log10(prevCalled) >= 0.5:
			return True
		else:
			return False
	
	@property
	def callBasis(self):
		
		retLogs = []
		for leaf in LeafClimber(self):
			if not leaf.incoming:
				pass
			elif isinstance(leaf.incoming[0].logRatio, Number) and isnan(leaf.incoming[0].logRatio):
				retLogs.append(0)
			elif isinstance(leaf.incoming[0].logRatio, Number):
				retLogs.append(leaf.incoming[0].logRatio)
			elif leaf.coverage > 1:
				retLogs.append(0)
			else:
				pass
		
		return sum(retLogs) / len(retLogs)

	@property
	def nonCanon(self):
		return self.incoming[0].nonCanon

	def hasNonCanon(self, minRatio : float=0.5) -> bool:
		
		if self.nonCanon <= 0:
			return False
		elif log10(self.incoming[0].depth) == 0:
			return False
		else:
			return minRatio <= log10(self.incoming[0].nonCanon) / log10(self.incoming[0].depth)

	@property
	def color(self) -> str:
		
		if not self.coverage:
			return HSV(0, 50, 50).hex()
		elif self.hasNonCanon():
			return HSV(300, 65, 75).hex()
		elif self.called:
			return HSV(120, 65, 75).hex()
		else:
			return HSV(230, 65, 75).hex()

class CanSNPBranch(Branch):

	properties = [
		("called", "int"), ("ancestral", "int"), ("nonCanon", "int"),
		("depth", "int"), ("ratio", "float"), ("logRatio", "float"),
		# ("color", "string")
	]

	weight = property(lambda self: self.ratio)

	@property
	def color(self):
		if self.outgoing.hasNonCanon():
			return HSV(300, 65, 75).hex()
		else:
			return metaCanSNPColorScale[float(self.weight)].hex()

	@property
	def hidden(self) -> bool:
		return self.outgoing and getattr(self.outgoing, "hidden", False)

class CanSNPTree(Tree):

	name : str = Default(lambda self: self.sessionName)
	nodeClass = CanSNPLeaf
	edgeClass = CanSNPBranch

	root : "CanSNPRoot" = cached_property(lambda self: CanSNPRoot(anchor="OUTER TOP", graph=self))

	properties = [
		("organism", "string"),
		("sessionName", "string"),
		("calledVariants", "list"),
		("reads", "int"),
		("bases", "int"),
		("bytes", "int")
	]

	@CachedDefault["allNodes"]
	def calledNodes(self) -> tuple[CanSNPLeaf]:
		return tuple(filter(lambda node:node.called, self.allNodes))
	@CachedDefault["allNodes"]
	def calledVariants(self) -> tuple[tuple[float,CanSNPLeaf]]:
		
		# paths : list[list[CanSNPBranch]] = [[branch] for branch in self.root.outgoing]
		TOLERANCE = 0.02
		# _ = combineClose(filter(lambda x:isinstance(x, float) and x > 0, map(*this.ratio, BranchWalker(self))))
		clusters = combineDensity(filter(lambda x:isinstance(x, float) and x > 0, map(*this.ratio, BranchWalker(self))))
		
		calledNodes = {}
		for tip in self.endNodes:
			for leaf in LeafClimber(tip):
				
				if leaf.called:
					for cluster in clusters:
						if leaf.incoming[0].ratio in cluster.data:
							fraction = cluster.mean
							break
					else:
						ValueError(f"Ungroupable SNP Call ratio {leaf.incoming[0].ratio}. Identified clusters: {clusters}")
					
					i = 0
					while fraction+i in calledNodes:
						if leaf.isAncestor(calledNodes[fraction+i]):
							if leaf.callBasis > 0.5:
								calledNodes[fraction+i] = leaf
							break
						i += 1
					else:
						if leaf.callBasis > 0.5:
							calledNodes[fraction+i] = leaf
					break
		
		# Filtering away the fractions that are really just sums of other fractions
		nodes = {n:f for f,n in calledNodes.items()}
		for fraction in sorted(calledNodes, key=lambda x:x%1):
			if fraction not in calledNodes:
				continue
			
			for parent in LeafClimber(calledNodes[fraction]):
				if parent in nodes:
					grouping = list(filter(lambda f: parent.isDescendant(calledNodes[f]), calledNodes))
					total = sum(map(lambda x:x%1, grouping))
					if nodes[parent]%1 - total < TOLERANCE*log10(len(grouping)):
						calledNodes.pop(nodes[parent], None)

		return tuple((f%1, v) for f,v in calledNodes.items())
	@CachedDefault["allNodes"]
	def nonCanonNodes(self) -> tuple[CanSNPLeaf]:
		return tuple(filter(lambda node:node.hasNonCanon(), self.allNodes))
	@CachedDefault["allNodes"]
	def noCoverageNodes(self) -> tuple[CanSNPLeaf]:
		return [node for node in self.allNodes if not isinstance(node.coverage, Number) or node.coverage <= 1 or node.coverage < 0.01*self.coverage]
	@CachedDefault["allNodes"]
	def coverage(self) -> float:
		return sum(n.coverage for n in self.allNodes if isinstance(n.coverage, Number)) / len(self.allNodes)
	
	@overload
	@classmethod
	def fromGraphML(cls : "type[CanSNPTree]", file : TextIO, /, *, root : str|None=None): ...
	@overload
	@classmethod
	def fromGraphML(cls : "type[CanSNPTree]", filename : FilePath, /, *, root : str|None=None): ...
	@classmethod
	@cache
	def fromGraphML(cls : "type[CanSNPTree]", file : TextIO|FilePath, /, *, root : str|None=None):
		return super().fromGraphML(file, root=root)

class CanSNPRoot(TreeRoot, CanSNPLeaf):
	id = "0"

def treeGraph(tree : CanSNPTree|list[CanSNPTree], filename : str|None=None, dark : bool=False, background : Color|str=None, **kwargs) -> CanSNPTree:

	if isinstance(tree, CanSNPTree) and not tree.allNodes:
		raise ValueError(f"Tree is empty, can't create graph.")
	elif isinstance(tree, list) and any(not t.allNodes for t in tree):
		raise ValueError(f"Tree is empty, can't create graph.")
	
	figure = tree

	if filename is not None:
		with open(f"{filename}.html", "w") as f:
			IllustrateFigure(figure).illustrateHTML(file=f, background=background)
		with open(f"{filename}.png", "wb") as f:
			IllustrateFigure(figure).illustratePNG(file=f, dark=dark, background=background, **kwargs)

	return figure

def coverageHistogram(*trees : CanSNPTree, filename : str|None=None, dark : bool=False, background : Color|str=None, **kwargs) -> Histogram:

	if any(not t.allNodes for t in trees):
		raise ValueError(f"Tree is empty, can't create graph.")
	
	figure = Histogram(
		*([node.coverage or 0 for node in t.allNodes] for t in trees),
		title="Variant Call Coverage", 
		labels=[t.name for t in trees],
		yTitle="MetaCanSNPs",
		xTitle="Coverage"
	)

	if filename is not None:
		with open(f"{filename}.html", "w") as f:
			IllustrateFigure(figure).illustrateHTML(file=f, background=background)
		with open(f"{filename}.png", "wb") as f:
			IllustrateFigure(figure).illustratePNG(file=f, dark=dark, background=background, **kwargs)

	return figure

def coverageBoxplot(*trees : CanSNPTree|list[CanSNPTree], filename : str|None=None, dark : bool=False, background : Color|str=None, **kwargs) -> BoxPlot:

	if any(not t.allNodes for t in trees):
		raise ValueError(f"Tree is empty, can't create graph.")
	figure = BoxPlot(
		*(
			[node.coverage or 0 for node in t.allNodes]
			for t in trees
		),
		title="Variant Call Coverage", 
		description="Combined coverage of all SNP's pertaining to each variant.", 
		yTitle="Coverage",
		xLabels=tuple(t.name for t in trees)
	)

	if filename is not None:
		with open(f"{filename}.html", "w") as f:
			IllustrateFigure(figure).illustrateHTML(file=f, background=background)
		with open(f"{filename}.png", "wb") as f:
			IllustrateFigure(figure).illustratePNG(file=f, dark=dark, background=background, **kwargs)

	return figure

def variantsTable(tree : CanSNPTree|list[CanSNPTree], filename : str|None=None, dark : bool=False, background : Color|str=None, **kwargs) -> Table:

	if isinstance(tree, CanSNPTree) and not tree.allNodes:
		raise ValueError(f"Tree is empty, can't create graph.")
	elif isinstance(tree, list) and any(not t.allNodes for t in tree):
		raise ValueError(f"Tree is empty, can't create graph.")
	
	figure = Table(
		[
			[node.genotype, percentEmphasize(fraction), node.incoming[0].depth]
			for fraction, node in tree.calledVariants
		],
		columns=["Variant", "Fraction", "Coverage"],
		title="Variants"
	)

	if filename is not None:
		with open(f"{filename}.html", "w") as f:
			IllustrateFigure(figure).illustrateHTML(file=f, background=background)
		with open(f"{filename}.md", "w") as f:
			figure.illustrateMARKDOWN(file=f)
		with open(f"{filename}.png", "wb") as f:
			IllustrateFigure(figure).illustratePNG(file=f, dark=dark, background=background, **kwargs)

	return figure

def nonCanonTable(tree : CanSNPTree|list[CanSNPTree], filename : str|None=None, dark : bool=False, background : Color|str=None, **kwargs) -> Table:

	if isinstance(tree, CanSNPTree) and not tree.allNodes:
		raise ValueError(f"Tree is empty, can't create graph.")
	elif isinstance(tree, list) and any(not t.allNodes for t in tree):
		raise ValueError(f"Tree is empty, can't create graph.")
	
	figure = Table(
		[
			[
				node.genotype,
				percentEmphasize(node.incoming[0].nonCanon/node.incoming[0].depth)
				if node.incoming[0].nonCanon > 0 and node.incoming[0].depth > 0
				else "NaN",
				node.incoming[0].depth
			]
			for node in tree.nonCanonNodes
		],
		columns=["Variant", "Fraction", "Coverage"],
		title="Non-Canon SNPs"
	)

	if filename is not None:
		with open(f"{filename}.html", "w") as f:
			IllustrateFigure(figure).illustrateHTML(file=f, background=background)
		with open(f"{filename}.md", "w") as f:
			figure.illustrateMARKDOWN(file=f)
		with open(f"{filename}.png", "wb") as f:
			IllustrateFigure(figure).illustratePNG(file=f, dark=dark, background=background, **kwargs)

	return figure

def badCoverageTable(tree : CanSNPTree|list[CanSNPTree], filename : str|None=None, dark : bool=False, background : Color|str=None, **kwargs) -> Table:

	if isinstance(tree, CanSNPTree) and not tree.allNodes:
		raise ValueError(f"Tree is empty, can't create graph.")
	elif isinstance(tree, list) and any(not t.allNodes for t in tree):
		raise ValueError(f"Tree is empty, can't create graph.")
	
	figure = Table(
		[
			[
				node.genotype,
				percentEmphasize(node.coverage/tree.coverage)
				if isinstance(node.coverage, Number) and node.coverage > 0 else
				"NaN",
				node.incoming[0].depth
			]
			for node in tree.noCoverageNodes
		],
		columns=["Variant", "% Average", "Coverage"],
		title="Bad Coverage"
	)

	if filename is not None:
		with open(f"{filename}.html", "w") as f:
			IllustrateFigure(figure).illustrateHTML(file=f, background=background)
		with open(f"{filename}.md", "w") as f:
			figure.illustrateMARKDOWN(file=f)
		with open(f"{filename}.png", "wb") as f:
			IllustrateFigure(figure).illustratePNG(file=f, dark=dark, background=background, **kwargs)

	return figure

#
#
#

def variantFractionConsistency(*treeGroups : list[tuple[Number,CanSNPTree]], filename : str|None=None, xUnit : str="Data", legend : bool|None=None, dark : bool=False, background : Color|str=None, **kwargs) -> Table:

	initialVariants = []
	LOGGER.debug(f"Getting initial variant call data from {len(treeGroups)} sets of trees. trees in sets: {', '.join(map(str, map(len, treeGroups)))}")
	for i, trees in enumerate(treeGroups):
		nHighest = max(map(lambda x:x[0],trees))
		fractions = {}
		LOGGER.debug(f"Tree group {i+1}/{len(treeGroups)}")
		biggestSizeTrees = list(filter(lambda x:x[0]==nHighest, trees))
		for n,tree in biggestSizeTrees:
			if n != nHighest:
				continue
			elif not tree.root.children:
				continue
			for f,node in tree.calledVariants:
				if node.genotype not in fractions:
					fractions[node.genotype] = [f]
				else:
					fractions[node.genotype].append(f)
				LOGGER.debug(f"Variant {node.genotype} found in tree {i+1}/{len(biggestSizeTrees)}")
			
		initialVariants.append(
			max((variant for variant in fractions), key=lambda v:sum(fractions[v])/len(fractions[v]))
		)
		LOGGER.debug(f"Variant {initialVariants[0]} determined as initial called variant for tree group {i+1}/{len(treeGroups)}")
	
	figure = LinePlot(
		*(
			[
				(
					n,
					first(
						(
							100
							for f,n in tree.calledVariants
							if n.genotype == initialVariant
						),
						0
					) if tree.root.children
					else 0
				)
				for j,(n,tree) in enumerate(trees)
				if not LOGGER.debug(f"Determining variant call for tree {j+1}/{len(treeGroups)}")
			]
			for i, initialVariant, trees in zip(_Iterators.Count() ,initialVariants, treeGroups)
			if not LOGGER.debug(f"Creating data for tree group {i+1}/{len(treeGroups)}")
		),
		**{
			"title":"Variant Calling Consistency",
			"description":", ".join(map(
				lambda n:str(Italic(SentenceCase(n)[:-1])),
				set(tree.organism for trees in treeGroups for n,tree in trees))),
			"yTitle":r"Variant Fraction %",
			"xTitle":f"{xUnit} / Sampling",
			"yLim":(0, 100),
		}|kwargs
	)
	if legend is True:
		figure.showLegend()
	elif legend is False:
		figure.hideLegend()

	if filename is not None:
		with open(f"{filename}.html", "w") as f:
			IllustrateFigure(figure).illustrateHTML(file=f, background=background)
		with open(f"{filename}.png", "wb") as f:
			IllustrateFigure(figure).illustratePNG(file=f, dark=dark, background=background, **kwargs)

	return figure

#
#
#

def interactiveFigure(tree : CanSNPTree|list[CanSNPTree], filename : str|None=None, dark : bool=False, background : Color|str=None, **kwargs) -> IllustrateFigure:

	if isinstance(tree, CanSNPTree) and not tree.allNodes:
		raise ValueError(f"Tree is empty, can't create graph.")
	elif isinstance(tree, list) and any(not t.allNodes for t in tree):
		raise ValueError(f"Tree is empty, can't create graph.")
	
	figureDocument = IllustrateFigure(IllustrativeCollection(
		treeGraph(tree),
		coverageHistogram(tree),
		coverageBoxplot(tree),
		variantsTable(tree),
		nonCanonTable(tree),
		badCoverageTable(tree),
		title=TitleCase(tree.sessionName),
		description=Italic(tree.organism.replace("_", " ").capitalize())
	))

	if filename is not None:
		with open(f"{filename}.html", "w") as f:
			figureDocument.illustrateHTML(file=f, background=background)
		with open(f"{filename}.png", "wb") as f:
			figureDocument.illustratePNG(file=f, dark=dark, background=background, **kwargs)
	
	return figureDocument

def figuresDocument(tree : CanSNPTree|list[CanSNPTree], filename : str|None=None, dark : bool=False, background : Color|str=None, **kwargs) -> IllustrateFigure:
	
	if isinstance(tree, CanSNPTree) and not tree.allNodes:
		raise ValueError(f"Tree is empty, can't create graph.")
	elif isinstance(tree, list) and any(not t.allNodes for t in tree):
		raise ValueError(f"Tree is empty, can't create graph.")
	
	figureDocument = IllustrateFigure(
		treeGraph(tree),
		coverageHistogram(tree),
		coverageBoxplot(tree),
		title=TitleCase(tree.sessionName)
	)
	figureDocument.header = [
		Head1(TitleCase(tree.sessionName)),
		Paragraph(Italic(tree.organism.replace("_", " ").capitalize()))
	]

	if filename is not None:
		with open(f"{filename}.html", "w") as f:
			figureDocument.illustrateHTML(file=f, background=background)
		with open(f"{filename}.png", "wb") as f:
			figureDocument.illustratePNG(file=f, dark=dark, background=background, **kwargs)
	
	return figureDocument

def tablesDocument(tree : CanSNPTree|list[CanSNPTree], filename : str|None=None, dark : bool=False, background : Color|str=None, **kwargs) -> IllustrateFigure:
	
	if isinstance(tree, CanSNPTree) and not tree.allNodes:
		raise ValueError(f"Tree is empty, can't create graph.")
	elif isinstance(tree, list) and any(not t.allNodes for t in tree):
		raise ValueError(f"Tree is empty, can't create graph.")
	
	tablesDocument = IllustrateFigure(
		variantsTable(tree),
		nonCanonTable(tree),
		badCoverageTable(tree),
		title=TitleCase(tree.sessionName)
	)
	tablesDocument.header = [
		Head1(TitleCase(tree.sessionName)),
		Paragraph(Italic(tree.organism.replace("_", " ").capitalize()))
	]

	if filename is not None:
		with open(f"{filename}.html", "w") as f:
			tablesDocument.illustrateHTML(file=f, background=background)
		with open(f"{filename}.png", "wb") as f:
			tablesDocument.illustratePNG(file=f, dark=dark, background=background, **kwargs)
	
	return tablesDocument