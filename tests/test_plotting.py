
from MetaCanSNPer.modules.Plotting import (
	interactiveFigure, figuresDocument, tablesDocument, CanSNPTree,
	treeGraph, coverageHistogram, coverageBoxplot, variantsTable, nonCanonTable, badCoverageTable
)
import os, pytest

@pytest.mark.skip
def test_collections():
	
	os.makedirs(os.path.splitext(__file__)[0], exist_ok=True)
	os.chdir(os.path.splitext(__file__)[0])

	tree, = CanSNPTree.fromGraphML("tree.graphml")

	namesAndFuncs = {
		"interactive_figure" : interactiveFigure,
		"figures_document" : figuresDocument,
		"tables_document" : tablesDocument,
	}

	for name, func in namesAndFuncs.items():
		start1 = os.stat(f"{name}.html").st_mtime if os.path.exists(f"{name}.html") else 0
		if start1 > 0:
			continue
		
		func(tree, filename=f"{name}", size=(1280, 720))
		
		assert os.stat(f"{name}.html").st_mtime > start1
		# assert os.stat(f"{name}.png").st_mtime > start2

@pytest.mark.skip
def test_singles():
	
	os.makedirs(os.path.splitext(__file__)[0], exist_ok=True)
	os.chdir(os.path.splitext(__file__)[0])

	trees = [tree for i in ["", "1", "2", "3"] for tree in CanSNPTree.fromGraphML(f"tree{i}.graphml")]
	tree = trees[0]
	
	namesAndFuncs = {
		"tree_graph" : treeGraph,
		"coverage_histogram" : coverageHistogram,
		"coverage_boxplot" : coverageBoxplot,
		"variants_table" : variantsTable,
		"non_canon_table" : nonCanonTable,
		"bad_coverage_table" : badCoverageTable,
	}

	for name, func in namesAndFuncs.items():
		start1 = os.stat(f"{name}.html").st_mtime if os.path.exists(f"{name}.html") else 0
		if start1 > 0:
			continue
		
		if name.startswith("coverage"):
			func(*trees, filename=f"{name}")
		else:
			func(tree, filename=f"{name}")
		
		assert os.stat(f"{name}.html").st_mtime > start1
		# assert os.stat(f"{name}.png").st_mtime > start2

"""
for d in os.listdir("."):
	if not os.path.isdir(d): continue
	front, back = (d.partition(']')[0]+']').split('-', 1)
	for f in os.listdir(d):
		if not f.startswith("tree_") or not f.endswith(".graphml"):
			os.remove(f"./{d}/{f}")
			continue
		n = f.split(".")[0].split("_")[-1]
		os.rename(f"./{d}/{f}", f"./{'-'.join([front, n, back])}.graphml")
	os.rmdir(f"./{d}")
"""

@pytest.mark.skip
def test_format_csv():
	from MetaCanSNPer.modules.Plotting import CanSNPTree
	from GeekyGadgets.Functions import parseNum
	from GeekyGadgets.DataFormats import CSV
	from GeekyGadgets.Formatting.Case import CamelCase
	from timeit import default_timer as timer
	import re
	
	os.makedirs(os.path.splitext(__file__)[0], exist_ok=True)
	os.chdir(os.path.splitext(__file__)[0])

	# FSC148[Bytes-01-25-20000000].graphml

	filePattern = re.compile(r"([^[]+)[[](\w+)-(\d+)-\d+-(\d+)[]][.]graphml")
	
	start = timer()
	data : dict[str,dict[str,list]] = {}
	for filename in os.listdir(os.path.join(".", "trees")):
		m = filePattern.match(filename)
		if m is None:
			print(f"{filename!r} failed to match pattern {filePattern}")
			continue
		genome, unitType, i, size = m.groups()
		
		# if parseNum(i) > 3:
		# 	continue
		size = parseNum(size)
		if unitType not in data:
			data[unitType] = {}
		if genome not in data[unitType]:
			data[unitType][genome] = []
		for tree in CanSNPTree.fromGraphML(open(os.path.join(".", "trees", filename), "r")):
			tree.id = f"{tree.id} | {i}"
			data[unitType][genome].append(
				(
					size,
					tree
				)
			)
	print(f"From GraphML to Python: {timer()-start} s.")

	assert len(set(tree.size for k1 in data for k2 in data[k1] for size, *trees in data[k1][k2] for tree in trees)) == 1

	os.makedirs(os.path.join(".", "treesCSV"), exist_ok=True)

	start = timer()
	for unitType in data:
		for genome in data[unitType]:
			
			csvFile = CSV()
			for size, tree in data[unitType][genome]:
				setattr(tree, CamelCase(unitType), size)
				tree.illustrateCSV(csv=csvFile)
			csvFile.save(os.path.join(".", "treesCSV", f"{genome}_{unitType}.csv"))
	print(f"From Python to CSV: {timer()-start} s.")

	filePattern = re.compile(r"(\w+)_(\w+)[.]csv")
	
	start = timer()
	data = {}
	for filename in os.listdir(os.path.join(".", "treesCSV")):
		m = filePattern.match(filename)
		if m is None:
			print(f"{filename!r} failed to match pattern {filePattern}")
			continue
		genome, unitType = m.groups()

		if unitType not in data:
			data[unitType] = {}
		if genome not in data[unitType]:
			data[unitType][genome] = []
		for tree in CanSNPTree.fromCSV(open(os.path.join(".", "treesCSV", filename), "r")):
			data[unitType][genome].append(
				(
					getattr(tree, CamelCase(unitType)),
					tree
				)
			)
	print(f"From CSV to Python: {timer()-start} s.")
	
	assert len(set(tree.size for k1 in data for k2 in data[k1] for size, *trees in data[k1][k2] for tree in trees)) == 1

	assert False


def test_consistency():
	from MetaCanSNPer.modules.Plotting import variantFractionConsistency, CanSNPTree
	from GeekyGadgets.Functions import parseNum
	from GeekyGadgets.Formatting.Case import CamelCase
	from GeekyGadgets.Logging import setLevel, logTo
	import re

	os.makedirs(os.path.splitext(__file__)[0], exist_ok=True)
	os.chdir(os.path.splitext(__file__)[0])
	
	setLevel(0)
	logTo(filename="test_consistency.log")

	filePattern = re.compile(r"(\w+)_(\w+)[.]csv")
	
	data = {}
	for filename in os.listdir(os.path.join(".", "treesCSV")):
		m = filePattern.match(filename)
		if m is None:
			print(f"{filename!r} failed to match pattern {filePattern}")
			continue
		genome, unitType = m.groups()

		if unitType not in data:
			data[unitType] = {}
		if genome not in data[unitType]:
			data[unitType][genome] = []
		for tree in CanSNPTree.fromCSV(open(os.path.join(".", "treesCSV", filename), "r")):
			data[unitType][genome].append(
				(
					getattr(tree, CamelCase(unitType)),
					tree
				)
			)
	
	assert len(set(tree.size for k1 in data for k2 in data[k1] for size, *trees in data[k1][k2] for tree in trees)) == 1
	
	for unitType in data:
		variantFractionConsistency(*(data[unitType][genome] for genome in data[unitType]), filename=f"variant_fraction_consistency[{unitType}]", size=(1280, 720), labels=list(data[unitType]), xUnit=unitType)
	
# test_consistency()