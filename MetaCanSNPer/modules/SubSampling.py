
from MetaCanSNPer.Globals import *
from MetaCanSNPer.core.Hooks import *
import hashlib, gunzip, subprocess
		

SUB_SAMPLE_NAMES = {
	"reads" : "Reads",
	"coverage" : "Coverage",
	"dilution" : "Dilution",
	"bases" : "Bases",
	"bytes" : "Bytes"
}

def randomReadGenerator(indices):
	while True:
		choices = list(range(len(indices[0])))
		random.shuffle(choices)
		for choice in choices:
			yield [reads[choice] for reads in indices]

def subSampleName(name : FilePath|DirectoryPath|str, type : Literal["reads","coverage","dilution","bytes"], *N, index : int=None) -> str:
	
	if index is None:
		bracketedID = f"[{SUB_SAMPLE_NAMES[type]}-{'-'.join(map(shortNumber, N))}]"
	else:
		bracketedID = f"[{SUB_SAMPLE_NAMES[type]}-{str(index).zfill(len(str(N[0])))}-{'-'.join(map(shortNumber, N))}]"
	if isinstance(name, FilePath):
		return name.name + bracketedID + "." + name.ext
	elif isinstance(name, DirectoryPath):
		return os.path.basename(name) + bracketedID
	else:
		newName, *ext = os.path.basename(name)[1:].split(".")
		name, ext = os.path.basename(name)[0]+newName, ".".join(ext)
		if ext:
			ext = "." + ext
		
		return name + bracketedID + ext

def readsProgressCallback(outData : list[list[list[list[int,int,int,int,int]|list[int,int,int,int]]]], /, *, reads, **kwargs) -> list[bool]:
	"""Returns `True` if condition has NOT been met."""
	return [sum(map(len, sampleData)) / reads for sampleData in outData]

def coverageProgressCallback(outData : list[list[list[list[int,int,int,int,int]|list[int,int,int,int]]]], /, *, targetCoverage, expectedCoverage, totalReads, **kwargs) -> list[bool]:
	"""Returns `True` if condition has NOT been met."""
	reads = totalReads * targetCoverage/expectedCoverage
	return [sum(map(len, sampleData)) / reads for sampleData in outData]

def dilutionProgressCallback(outData : list[list[list[list[int,int,int,int,int]|list[int,int,int,int]]]], /, *, dilution, totalReads, **kwargs) -> list[bool]:
	"""Returns `True` if condition has NOT been met."""
	reads = totalReads / dilution
	return [sum(map(len, sampleData)) / reads for sampleData in outData]

def basesProgressCallback(outData : list[list[list[list[int,int,int,int,int]|list[int,int,int,int]]]], /, *, bases, **kwargs) -> list[bool]:
	"""Returns `True` if condition has NOT been met."""
	return [sum(itertools.chain(map(lambda x:map(lambda y:y[0], x), sampleData))) / bases for sampleData in outData]

def bytesProgressCallback(outData : list[list[list[list[int,int,int,int,int]|list[int,int,int,int]]]], /, *, bytes, **kwargs) -> list[bool]:
	"""Returns `True` if condition has NOT been met."""
	return [sum(itertools.chain(map(lambda x:map(lambda y:y[2]-y[1], x), sampleData))) / bytes for sampleData in outData]

def getProgressCallback(subSamplingType : str, *args, **kwargs):
	match subSamplingType:
		case "reads":
			return partial(readsProgressCallback, reads=args[0], **kwargs)
		case "coverage":
			return partial(coverageProgressCallback, targetCoverage=args[0], expectedCoverage=args[1], totalReads=args[2], **kwargs)
		case "dilution":
			return partial(dilutionProgressCallback, dilution=args[0], totalReads=args[1], **kwargs)
		case "bases":
			return partial(basesProgressCallback, bases=args[0], **kwargs)
		case "bytes":
			return partial(bytesProgressCallback, bytes=args[0], **kwargs)
	

@overload
def splitFastq(files : int, source : FilePath|FileList[FilePath], *, reads : list[int], **kwargs) -> list[tuple[str]]: ...
@overload
def splitFastq(files : int, source : FilePath|FileList[FilePath], *, coverage : list[int,int], **kwargs) -> list[tuple[str]]: ...
@overload
def splitFastq(files : int, source : FilePath|FileList[FilePath], *, dilution : list[int], **kwargs) -> list[tuple[str]]: ...
@overload
def splitFastq(files : int, source : FilePath|FileList[FilePath], *, bytes : list[int], **kwargs) -> list[tuple[str]]: ...
@overload
def splitFastq(files : int, source : FilePath|FileList[FilePath], *,
			   reads : list[int]=None, dilution : list[int]=None, coverage : list[int,int]=None, bytes : list[int]=None, bases : list[int]=None,
			   outDir : DirectoryPath=None, hooks=GlobalHooks, steps : int=100) -> list[tuple[str]]: ...
def splitFastq(files : int, source : FilePath|FileList[FilePath], *,
			   outDir : DirectoryPath=None, hooks=GlobalHooks, steps : int=100,
			   **kwargs) -> list[tuple[str]]:
	
	LOGGER.info(f"Sub Sampling: From {source}.")

	if isinstance(source, str):
		source = FileList([FilePath(source)])
	elif isinstance(source, Iterable):
		source = FileList(FilePath(name) for name in source)
	
	if all(filepath.endswith(".gz") for filepath in source):
		LOGGER.info(f"Sub Sampling: From and to gzip-data.")
		dataOpen = gunzip.gzip.open
	elif not any(filepath.endswith(".gz") for filepath in source):
		LOGGER.info(f"Sub Sampling: From and to raw-data.")
		dataOpen = open
	else:
		raise ValueError(f"Files are not consistent in their compression file extensions: {source}")
	
	LOGGER.info(f"Sub Sampling: Creating Read Index.")
	readsIndex = []
	for filepath in source:
		with dataOpen(filepath, "rb") as file:
			file : BinaryIO
			
			pos = -1
			readList = []
			while pos != (pos := file.tell()):
				if not file.readline().startswith(b"@"):
					continue
				read = []
				for line in file:
					if not line.strip().isalpha():
						lineSep = line
						break
					read.append(line)
				readLength = sum(map(len, map(str.strip, read)))
				if lineSep.strip() == b"+":
					file.seek(1, readLength)
				else:
					file.seek(1, -len(lineSep))
				readList.append([readLength, pos, file.tell()])

			readsIndex.append(readList)
	LOGGER.info(f"Sub Sampling: Found {', '.join(map(len, readsIndex))} reads for the source file(s).")

	for name in ["reads", "dilution", "coverage", "bytes", "bases"]:
		if name in kwargs:
			LOGGER.info(f"Sub Sampling: Using {values[0]} {SUB_SAMPLE_NAMES[name]}.")
			varName, values = name, kwargs[name]
			progressCallback = getProgressCallback(name, *values, totalReads=len(readsIndex[0]))
			break
	else:
		raise ValueError("No sub sampling information given, check keyword arguments of `splitFastq`.")

	hooks.trigger("SplitFastqStarting", {"name" : source.name, "value" : 0.0})
	outNames = [tuple((outDir or filepath.directory) / subSampleName(filepath, varName, files, *values, index=i+1) for filepath in source) for i in range(files)]

	if all(os.path.exists(filename) for filenames in outNames for filename in filenames):
		hooks.trigger("SplitFastqSkipped", {"name" : source.name, "value" : 2})
		return outNames
	
	dataFiles : list[BinaryIO] = [dataOpen(filepath, "rb") for filepath in source]
	outData : list[list[BinaryIO]] = [[[] for filename in filenames] for filenames in outNames]
	outFiles : list[list[BinaryIO]] = [[dataOpen(filename, "wb") for filename in filenames] for filenames in outNames]
	
	threshold = 0
	
	progressVector = [-1 for _ in outData]
	# First 1/4 of progress
	LOGGER.info(f"Sub Sampling: Read Selection.")
	for readSet in itertools.batched(randomReadGenerator(readsIndex), len(outData)):
		if progressVector == (progressVector := progressCallback(outData)):
			break
		if (1/4) * sum(progressVector) / len(progressVector) >= threshold:
			hooks.trigger("SplitFastqProgress", {"name" : source.name, "value" : min(1.0, (1/4) * sum(progressVector) / len(progressVector))})
			threshold = (int((1/4) * steps * sum(progressVector) / len(progressVector))+1) / steps
		for notDone, sampleData, outReads in zip(progressVector, outData, readSet):
			if notDone:
				for fileData, read in zip(sampleData, outReads):
					fileData.append(read)
	
	# 2/4 of progress
	LOGGER.info(f"Sub Sampling: Read Aggregation.")
	readsAggregates = [[] for _ in outFiles[0]]
	for data, files in zip(outData, outFiles):
		if (1/4) + (1/4) * i / len(outFiles) >= threshold:
			hooks.trigger("SplitFastqProgress", {"name" : source.name, "value" : min(1.0, (1/4) + (1/4) * i / len(outFiles))})
			threshold = (int(steps * ((1/4) + (1/4) * i) / len(outFiles))+1) / steps
		for fileN, file, reads in zip(range(len(files), files, data)):
			readsAggregates[fileN].extend((file, read) for read in reads)

	# 3/4 of progress
	LOGGER.info(f"Sub Sampling: Read Sorting.")
	for i, aggregate in enumerate(readsAggregates):
		if (2/4) + (1/4) * i / len(readsAggregates) >= threshold:
			hooks.trigger("SplitFastqProgress", {"name" : source.name, "value" : min(1.0, (2/4) + (1/4) * i / len(readsAggregates))})
			threshold = (int(steps * ((2/4) + (1/4) * i) / len(readsAggregates))+1) / steps
		aggregate.sort(key=lambda x:x[1][1])

	# 4/4 of progress
	LOGGER.info(f"Sub Sampling: Writing Files.")
	for aggregate, dataFile in zip(readsAggregates, dataFiles):
		if (3/4) + (1/4) * i / len(readsAggregates) >= threshold:
			hooks.trigger("SplitFastqProgress", {"name" : source.name, "value" : min(1.0, (3/4) + (1/4) * i / len(readsAggregates))})
			threshold = (int(steps * ((3/4) + (1/4) * i) / len(readsAggregates))+1) / steps
		for file, read in aggregate:
			dataFile.seek(read[1])
			file.write(dataFile.read(read[2] - read[1]))
	hooks.trigger("SplitFastqProgress", {"name" : source.name, "value" : 1.0})

	for files in outFiles:
		for file in files:
			file.close()
	hooks.trigger("SplitFastqFinished", {"name" : source.name, "value" : 3})

	LOGGER.info(f"Sub Sampling: Finished!")
	
	return outNames
