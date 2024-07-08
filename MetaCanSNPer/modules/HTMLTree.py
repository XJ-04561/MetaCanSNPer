
from PseudoPathy import FilePath
from GeekyGadgets.Illustrative import Tree, Branch, Leaf
from GeekyGadgets.TypeHinting import *

__all__ = (
	"CanSNPLeaf", "CanSNPBranch", "CanSNPTree", "drawGraphML"
)

class CanSNPLeaf(Leaf):
	"""Properties:
	genotype : str"""
	
	@property
	def hidden(self) -> bool:
		if not self.incoming:
			return False
		if not self.incoming[0].properties["depth"]:
			return False
		if not isinstance(self.incoming[0].properties["depth"], Number):
			return False
		if any(not child.hidden for child in self.children):
			return False
		if 0.01 < self.incoming[0].properties["ratio"]:
			return False
		return True

	@property
	def color(self) -> str:
		
		if not self.incoming or not self.incoming[0].properties["depth"] or not isinstance(self.incoming[0].properties["depth"], Number):
			return "#303030"
		elif 0.05 < self.incoming[0].properties["nonCanon"] / self.incoming[0].properties["depth"]:
			return "#ff30ff"
		calledSNPs = []
		ratios = []
		node = self
		while node.incoming:
			ratios.append(node.incoming[0].properties["ratio"])
			calledSNPs.append(node.incoming[0].properties["called"])
			node = node.incoming[0].pair[0]
		
		if any(prevCalled > 1 and isinstance(thisNode, Number) and isinstance(prevNode, Number) and 1.1 < thisNode/prevNode for thisNode, prevNode, prevCalled in zip(ratios, ratios[1:], calledSNPs[1:])):
			return "#ff30ff"
		
		return "#20ff20"

class CanSNPBranch(Branch):
	"""Properties:
	called : int
	ancestral : int
	nonCanon : int
	depth : int
	ratio : float
	logRatio : float"""

	weight = property(lambda self: self.properties["ratio"])

class CanSNPTree(Tree):
	nodeClass = CanSNPLeaf
	edgeClass = CanSNPBranch

	weightProp : str = "ratio"

def drawGraphML(filename : FilePath):
	tree = CanSNPTree.fromGraphML(open(filename, "r"))

	tree.illustrate("HTML", filename=filename+".html", nameProp="genotype")
