
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
def splitFastq(files : int, source : FilePath|FileList[FilePath], *, reads : int, **kwargs) -> list[tuple[str]]: ...
@overload
def splitFastq(files : int, source : FilePath|FileList[FilePath], *, subFactor : int|float, **kwargs) -> list[tuple[str]]: ...
@overload
def splitFastq(files : int, source : FilePath|FileList[FilePath], *, coverage : int, **kwargs) -> list[tuple[str]]: ...
@overload
def splitFastq(files : int, source : FilePath|FileList[FilePath], *,
			   reads : int=None, subFactor : int|float=None, coverage : int=None,
			   outDir : DirectoryPath=None, hooks=GlobalHooks, randomiser=None, steps : int=100) -> list[tuple[str]]: ...
def splitFastq(files : int, source : FilePath|FileList[FilePath], *,
			   reads : int=None, subFactor : int|float=None, coverage : int=None,
			   outDir : DirectoryPath=None, hooks=GlobalHooks, randomiser=None, steps : int=100) -> list[tuple[str]]:
	if not isinstance(source, list):
		source = FileList([source])
	
	if all(filepath.endswith(".gz") for filepath in source):
		dataOpen = gunzip.gzip.open
	elif not any(filepath.endswith(".gz") for filepath in source):
		dataOpen = open
	else:
		raise ValueError(f"Files are not consistent in their compression file extensions: {source}")
	
	if randomiser is None:
		randomiser = random.Random()
		randomiser.seed(hashlib.md5("".join(source).encode("utf-8")).digest())
	
	splitNames = [((outDir or filepath.directory) / filepath.name, filepath.ext) for filepath in source]
	
	NLength = len(str(files))
	outNames = [tuple(f"{name}[{{splitType}}-{str(i+1).zfill(NLength)}-{files}].{ext}" for name, ext in splitNames) for i in range(files)]
	
	if reads is not None:
		splitByReads(reads, source, outNames, dataOpen, randomiser, hooks=hooks)
	elif coverage is not None:
		splitByCoverage(coverage, source, outNames, dataOpen, randomiser, hooks=hooks)
	else:
		splitByDilution(subFactor or files, source, outNames, dataOpen, randomiser, hooks=hooks)

@overload
def splitByReads(reads : int, source : FilePath, outNames : list[FilePath], dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], hooks : Hooks=GlobalHooks, steps : int=100): ...
@overload
def splitByReads(reads : int, source : FileList[FilePath], outNames : list[list[FilePath]], dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], hooks : Hooks=GlobalHooks, steps : int=100): ...
def splitByReads(reads : int, source : FilePath|FileList[FilePath], outNames : list[FilePath]|list[list[FilePath]], dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], hooks : Hooks=GlobalHooks, steps : int=100):

	hooks.trigger("SplitFastqStarting", {"name" : source.name, "value" : 0.0})
	outNames = [[FilePath(filename.format(splitType="ReadSplit")) for filename in filenames] for filenames in outNames]

	if all(os.path.exists(filename) for filenames in outNames for filename in filenames):
		hooks.trigger("SplitFastqSkipped", {"name" : source.name, "value" : 2})
		return outNames
	
	dataFiles : list[BinaryIO] = [dataOpen(name, "rb") for name in source]
	outFiles = [[dataOpen(filename, "wb") for filename in filenames] for filenames in outNames]
	
	readsWritten = [0 for _ in range(len(outFiles))]
	
	threshold = 0
	for choice in equalSampling(files, randomiser):
		if all(reads <= readsInFile for readsInFile in readsWritten):
			break
		while abs(choice) < len(outFiles): # choice is already satisfied, select next.
			if readsWritten[choice] < reads:
				break # choice now selects a file which needs more data.
			choice -= 1
		else:
			break # All files have the correct size.

		traversedBytes = sum(of.write(df.readline()) for of, df in zip(outFiles[choice], dataFiles) for _ in range(4))

		if traversedBytes == 0:
			for df in dataFiles:
				df.seek(0)
			if 0 == sum(of.write(df.readline()) for of, df in zip(outFiles[choice], dataFiles) for _ in range(4)):
				raise EOFError(f"No reads can be read from files {[f.name for f in dataFiles]}")
		
		readsWritten[choice] += 1
			
		if (sum(readsWritten) / files) / reads > threshold:
			hooks.trigger("SplitFastqProgress", {"name" : source.name, "value" : min(1.0, (sum(readsWritten) / files) / reads)})
			threshold = (int(steps * ((sum(readsWritten) / files) / reads)) + 1) / steps

	for files in outFiles:
		for file in files:
			file.close()
	hooks.trigger("SplitFastqFinished", {"name" : source.name, "value" : 3})
	return outNames

@overload
def splitByCoverage(coverage : int, source : FilePath, outNames : list[FilePath], dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], hooks : Hooks=GlobalHooks, steps : int=100): ...
@overload
def splitByCoverage(coverage : int, source : FileList[FilePath], outNames : list[list[FilePath]], dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], hooks : Hooks=GlobalHooks, steps : int=100): ...
def splitByCoverage(coverage : int, source : FilePath|FileList[FilePath], outNames : list[FilePath]|list[list[FilePath]], dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], hooks : Hooks=GlobalHooks, steps : int=100):

	raise NotImplementedError(f"Not yet implemented!")

	hooks.trigger("SplitFastqStarting", {"name" : source.name, "value" : 0.0})
	outNames = [[FilePath(filename.format(splitType="CoverageSplit")) for filename in filenames] for filenames in outNames]

	if all(os.path.exists(filename) for filenames in outNames for filename in filenames):
		hooks.trigger("SplitFastqSkipped", {"name" : source.name, "value" : 2})
		return outNames
	
	dataFiles : list[BinaryIO] = [dataOpen(name, "rb") for name in source]
	outFiles = [[dataOpen(filename, "wb") for filename in filenames] for filenames in outNames]
	
	nBytes = 0
	for file in dataFiles:
		nBytes += file.seek(0, 2)
		file.seek(0, 0)
	


	bytesPerOutfile = nBytes / coverage
	bytesWritten = [0 for _ in range(len(outFiles))]
	
	threshold = 0
	for choice in equalSampling(files, randomiser):
		if bytesPerOutfile < sum(bytesWritten) / files:
			break
		while abs(choice) < len(outFiles)+1: # choice is already satisfied, select next.
			choice -= 1
			if bytesWritten[choice] < bytesPerOutfile:
				break # choice now selects a file which needs more data.
		else:
			break # All files have the correct size.

		traversed = sum(of.write(df.readline()) for of, df in zip(outFiles[choice], dataFiles) for _ in range(4))
		
		bytesWritten[choice] += traversed

		if traversed == 0: # End of File reached, time to restart
			for df in dataFiles:
				df.seek(0)
		
		if (sum(bytesWritten) / files) / bytesPerOutfile > threshold:
			hooks.trigger("SplitFastqProgress", {"name" : source.name, "value" : min(1.0, (sum(bytesWritten) / files) / bytesPerOutfile)})
			threshold = (int(steps * ((sum(bytesWritten) / files) / bytesPerOutfile)) + 1) / steps

	for files in outFiles:
		for file in files:
			file.close()
	hooks.trigger("SplitFastqFinished", {"name" : source.name, "value" : 3})
	return outNames

@overload
def splitByDilution(dilutionFactor : int, source : FilePath, outNames : list[FilePath], dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], hooks : Hooks=GlobalHooks, steps : int=100): ...
@overload
def splitByDilution(dilutionFactor : int, source : FileList[FilePath], outNames : list[list[FilePath]], dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], hooks : Hooks=GlobalHooks, steps : int=100): ...
def splitByDilution(dilutionFactor : int, source : FilePath|FileList[FilePath], outNames : list[FilePath]|list[list[FilePath]], dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], hooks : Hooks=GlobalHooks, steps : int=100):
	
	hooks.trigger("SplitFastqStarting", {"name" : source.name, "value" : 0.0})
	outNames = [[FilePath(filename.format(splitType="DilutionSplit")) for filename in filenames] for filenames in outNames]

	if all(os.path.exists(filename) for filenames in outNames for filename in filenames):
		hooks.trigger("SplitFastqSkipped", {"name" : source.name, "value" : 2})
		return outNames
	
	dataFiles : list[BinaryIO] = [dataOpen(name, "rb") for name in source]
	outFiles = [[dataOpen(filename, "wb") for filename in filenames] for filenames in outNames]
	
	nBytes = 0
	for file in dataFiles:
		nBytes += file.seek(0, 2)
		file.seek(0, 0)
	
	bytesPerOutfile = nBytes / dilutionFactor
	bytesWritten = [0 for _ in range(len(outFiles))]
	
	threshold = 0
	for choice in equalSampling(files, randomiser):
		if bytesPerOutfile < sum(bytesWritten) / files:
			break
		while abs(choice) < len(outFiles)+1: # choice is already satisfied, select next.
			choice -= 1
			if bytesWritten[choice] < bytesPerOutfile:
				break # choice now selects a file which needs more data.
		else:
			break # All files have the correct size.

		traversed = sum(of.write(df.readline()) for of, df in zip(outFiles[choice], dataFiles) for _ in range(4))
		
		bytesWritten[choice] += traversed

		if traversed == 0: # End of File reached, time to restart
			for df in dataFiles:
				df.seek(0)
		
		if (sum(bytesWritten) / files) / bytesPerOutfile > threshold:
			hooks.trigger("SplitFastqProgress", {"name" : source.name, "value" : min(1.0, (sum(bytesWritten) / files) / bytesPerOutfile)})
			threshold = (int(steps * ((sum(bytesWritten) / files) / bytesPerOutfile)) + 1) / steps

	for files in outFiles:
		for file in files:
			file.close()
	hooks.trigger("SplitFastqFinished", {"name" : source.name, "value" : 3})
	return outNames