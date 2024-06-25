
from MetaCanSNPer.Globals import *
from MetaCanSNPer.core.Hooks import *
import hashlib
import gunzip

class DummyIO:
	def write(self, data, **kwargs):
		return len(data)

def equalSampling(N, randomiser : random.Random):
	counts = [1 for _ in range(N)]
	weights = [0 for _ in range(N)]
	choices = list(range(N))
	while True:
		for i in range(N):
			weights[i] = 1 - (counts[i] / sum(counts))
		for i in randomiser.choices(choices, weights, k=100):
			yield i

@overload
def splitFastq(files : int, source : Path|list[Path], *, reads : int, **kwargs) -> list[tuple[str]]: ...
@overload
def splitFastq(files : int, source : Path|list[Path], *, subFactor : int|float, **kwargs) -> list[tuple[str]]: ...
@overload
def splitFastq(files : int, source : Path|list[Path], *, coverage : int, **kwargs) -> list[tuple[str]]: ...
@overload
def splitFastq(files : int, source : Path|list[Path], *,
			   reads : int=None, subFactor : int|float=None, coverage : int=None,
			   outDir : DirectoryPath=None, hooks=GlobalHooks, randomiser=None, steps : int=100) -> list[tuple[str]]: ...
def splitFastq(files : int, source : Path|list[Path], *,
			   reads : int=None, subFactor : int|float=None, coverage : int=None,
			   outDir : DirectoryPath=None, hooks=GlobalHooks, randomiser=None, steps : int=100) -> list[tuple[str]]:

	sourceName = getattr(source, "name", os.path.basename(source) if not isinstance(source, list) else os.path.basename(source[0]))

	if all(filename.endswith(".gz") for filename in source):
		dataOpen = gunzip.gzip.open
	elif not any(filename.endswith(".gz") for filename in source):
		dataOpen = open
	else:
		raise ValueError(f"Files are not consistent in their compression file extensions: {source}")
	
	if randomiser is None:
		randomiser = random.Random()
		randomiser.seed(hashlib.md5("".join(source).encode("utf-8")).digest())
	
	splitNames = []
	for filepath in source:
		name, ext = os.path.basename(filepath)[1:].split(".", 1)
		splitNames.append((os.path.join(outDir or os.path.dirname(filepath), os.path.basename(filepath)[0]+name), ext))

	NLength = len(str(files))
	outNames = [tuple(f"{name}-{str(i+1).zfill(NLength)}-{files}.{ext}" for name, ext in splitNames) for i in range(files)]

	
	dataFiles : list[BinaryIO] = [dataOpen(name, "rb") for name in source]

	outFiles = [[dataOpen(filename, "wb") if not os.path.exists(filename) or os.stat(filename).st_size < 1000 else DummyIO() for filename in source] for source in outNames]

	if all(isinstance(file, DummyIO) for files in outFiles for file in files):
		hooks.trigger("SplitFastqSkipped", {"name" : sourceName, "value" : 2})
		return outNames
	
	nBytes = sum(file.seek(0, 2) for file in dataFiles)
	if reads is not None:
		splitByReads(reads)
	elif coverage is not None:
		splitByCoverage(coverage)
	else:
		splitByDilution(subFactor or files)
		
def splitByDilution(dilutionFactor):
	bytesPerOutfile = nBytes / M
	for file in dataFiles:
		file.seek(0, 0)

	hooks.trigger("SplitFastqStarting", {"name" : sourceName, "value" : 0.0})
	steps : list[float] = [2**6] + [(steps-i)/steps for i in range(steps+1)]
	byteCounts = [0 for _ in range(files)]
	
	for choice in equalSampling(files, randomiser):
		if bytesPerOutfile < sum(byteCounts) / files:
			break

		traversed = sum(of.write(df.readline()) for of, df in zip(outFiles[choice], dataFiles) for _ in range(4))
		
		byteCounts[choice] += traversed

		if traversed == 0: # End of File reached, time to restart
			for df in dataFiles:
				df.seek(0)
		
		if (sum(byteCounts) / files) / bytesPerOutfile > steps[-1]:
			hooks.trigger("SplitFastqProgress", {"name" : sourceName, "value" : min(1.0, (sum(byteCounts) / files) / bytesPerOutfile)})
			steps.pop()

	for files in outFiles:
		for file in files:
			file.close()
	hooks.trigger("SplitFastqFinished", {"name" : sourceName, "value" : 3})
	return outNames