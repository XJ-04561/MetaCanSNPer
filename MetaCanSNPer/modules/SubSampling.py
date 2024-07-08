
from MetaCanSNPer.Globals import *
import MetaCanSNPer.Globals as Globals
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

def subSampleName(name : FilePath|DirectoryPath|str, type : Literal["reads","coverage","dilution","bases","bytes"], *N, index : int=None) -> str:
	
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

class ProgressTracker(list):
	def __init__(self, n : int, goal : int|float, factor : int|float=1, **kwargs):
		super().__init__(0 for _ in range(n))
		self.n = n
		self.goal = goal
		self.factor = factor
		for name, value in kwargs.items():
			setattr(self, name, value)
	
	def __iter__(self):
		for p in super().__iter__():
			yield self.factor * p / self.goal

	def update(self, other):
		for i, x in enumerate(other):
			self[i] += 1
	
	def mean(self):
		return self.goal * self.norm()
	
	def norm(self):
		return sum(self) / self.n

class ReadsProgressTracker(ProgressTracker): ...

class CoverageProgressTracker(ProgressTracker): ...

class DilutionProgressTracker(ProgressTracker): ...

class BasesProgressTracker(ProgressTracker):
	def update(self, other):
		for i, x in enumerate(other):
			self[i] += sum(map(lambda y:y[0][0], x))

class BytesProgressTracker(ProgressTracker):
	def update(self, other):
		for i, x in enumerate(other):
			self[i] += sum(map(lambda y:y[0][2] - y[0][1], x))

def getProgressCallback(subSamplingType : str, n : int, *args, **kwargs):
	match subSamplingType:
		case "reads":
			return ReadsProgressTracker(n, goal=args[0])
		case "coverage":
			return CoverageProgressTracker(n, goal=args[0], factor=args[1] / args[2])
		case "dilution":
			return DilutionProgressTracker(n, goal=1/args[0], factor=1/args[1])
		case "bases":
			return BasesProgressTracker(n, goal=args[0])
		case "bytes":
			return BytesProgressTracker(n, goal=args[0])

def createReadsIndex(filepath : FilePath, ioFunc : Callable[[str,str],BinaryIO]):
	
	if os.path.exists(f"{filepath.name}_readsIndex.csv"):
		return [[int(x) for x in row.strip().split(",")] for row in open(f"{filepath.name}_readsIndex.csv", "r")]
	
	pos = -1
	readList = []
	with ioFunc(filepath, "rb") as file:
		while pos != (pos := file.tell()):
			if not file.readline().startswith(b"@"):
				continue
			read = []
			for line in file:
				if not line.strip().isalpha():
					qualHeader = line
					break
				read.append(line)
			
			if qualHeader.startswith(b"+"):
				file.seek(sum(map(len, read)), 1)
			elif qualHeader:
				file.seek(-len(qualHeader), 1)
			readList.append([sum(map(len, map(bytes.strip, read))), pos, file.tell()])
	try:
		with open(f"{filepath.name}_readsIndex.csv", "w") as file:
			for row in readList:
				file.write(f"{row[0]},{row[1]},{row[2]}\n")
	except:
		pass
			
	return readList

@overload
def splitFastq(samples : int, source : FilePath|FileList[FilePath], *, reads : list[int], **kwargs) -> list[tuple[str]]: ...
@overload
def splitFastq(samples : int, source : FilePath|FileList[FilePath], *, coverage : list[int,int], **kwargs) -> list[tuple[str]]: ...
@overload
def splitFastq(samples : int, source : FilePath|FileList[FilePath], *, dilution : list[int], **kwargs) -> list[tuple[str]]: ...
@overload
def splitFastq(samples : int, source : FilePath|FileList[FilePath], *, bytes : list[int], **kwargs) -> list[tuple[str]]: ...
@overload
def splitFastq(samples : int, source : FilePath|FileList[FilePath], *,
			   reads : list[int]=None, dilution : list[int]=None, coverage : list[int,int]=None, bytes : list[int]=None, bases : list[int]=None,
			   outDir : DirectoryPath=None, hooks=GlobalHooks, steps : int=100) -> list[tuple[str]]: ...
def splitFastq(samples : int, source : FilePath|FileList[FilePath], *,
			   outDir : DirectoryPath=None, hooks=GlobalHooks, steps : int=100,
			   **kwargs) -> list[tuple[str]]:
	
	Globals.LOGGER.info(f"Sub Sampling: From {source}.")

	if isinstance(source, str):
		source = FileList([FilePath(source)])
	elif isinstance(source, Iterable):
		source = FileList(FilePath(name) for name in source)
	
	if all(filepath.endswith(".gz") for filepath in source):
		Globals.LOGGER.info("Sub Sampling: From and to gzip-data.")
		dataOpen = gunzip.gzip.open
	elif not any(filepath.endswith(".gz") for filepath in source):
		Globals.LOGGER.info("Sub Sampling: From and to raw-data.")
		dataOpen = open
	else:
		raise ValueError(f"Files are not consistent in their compression file extensions: {source}")
	
	Globals.LOGGER.info("Sub Sampling: Creating Read Index.")
	readsIndex = []
	for filepath in source:
		readsIndex.append(createReadsIndex(filepath, dataOpen))
	Globals.LOGGER.info(f"Sub Sampling: Found {', '.join(map(str, map(len, readsIndex)))} reads for the source file(s).")

	for name in ["reads", "dilution", "coverage", "bases", "bytes"]:
		if name in kwargs:
			varName, values = name, kwargs[name]
			Globals.LOGGER.info(f"Sub Sampling: Using {values[0]} {SUB_SAMPLE_NAMES[name]}.")
			progressTracker = getProgressCallback(name, samples, *values, totalReads=len(readsIndex[0]))
			break
	else:
		raise ValueError("No sub sampling information given, check keyword arguments of `splitFastq`.")

	hooks.trigger("SplitFastqStarting", {"name" : source.name, "value" : 0.0})
	outNames = [tuple((outDir or filepath.directory) / subSampleName(filepath, varName, samples, *values, index=i+1) for filepath in source) for i in range(files)]

	if all(os.path.exists(filename) for filenames in outNames for filename in filenames):
		hooks.trigger("SplitFastqSkipped", {"name" : source.name, "value" : 2})
		return outNames
	
	dataFiles : list[BinaryIO] = [dataOpen(filepath, "rb") for filepath in source]
	outData : list[list[BinaryIO]] = [[[] for filename in filenames] for filenames in outNames]
	outFiles : list[list[BinaryIO]] = [[dataOpen(filename, "wb") for filename in filenames] for filenames in outNames]
	
	threshold = 0
	# First 1/4 of progress
	Globals.LOGGER.info(f"Sub Sampling: Read Selection.")
	for readSet in itertools.batched(randomReadGenerator(readsIndex), len(outData)):
		for prog, sampleData, outReads in zip(progressTracker, outData, readSet):
			if prog < 1.0:
				for fileData, read in zip(sampleData, outReads):
					fileData.append(read)
		progressTracker.update(outData)
		if progressTracker.norm() >= 1.0:
			break
		elif steps * (1/4) * progressTracker.norm() >= threshold:
			hooks.trigger("SplitFastqProgress", {"name" : source.name, "value" : min(1.0, (1/4) * progressTracker.norm())})
			threshold += 1
	
	# 2/4 of progress
	Globals.LOGGER.info(f"Sub Sampling: Read Aggregation.")
	readsAggregates = [[] for _ in outFiles[0]]
	for i, data, files in zip(itertools.count(), outData, outFiles):
		if steps * ((1/4) + (1/4) * i / len(outFiles)) >= threshold:
			hooks.trigger("SplitFastqProgress", {"name" : source.name, "value" : min(1.0, (1/4) + (1/4) * i / len(outFiles))})
			threshold += 1
		for fileN, file, reads in zip(itertools.count(), files, data):
			readsAggregates[fileN].extend((file, read) for read in reads)

	# 3/4 of progress
	Globals.LOGGER.info(f"Sub Sampling: Read Sorting.")
	for i, aggregate in enumerate(readsAggregates):
		if steps * ((2/4) + (1/4) * i / len(readsAggregates)) >= threshold:
			hooks.trigger("SplitFastqProgress", {"name" : source.name, "value" : min(1.0, (2/4) + (1/4) * i / len(readsAggregates))})
			threshold += 1
		aggregate.sort(key=lambda x:x[1][1])

	# 4/4 of progress
	Globals.LOGGER.info(f"Sub Sampling: Writing Files.")
	for i, aggregate, dataFile in zip(itertools.count(), readsAggregates, dataFiles):
		if steps * ((3/4) + (1/4) * i / len(readsAggregates)) >= threshold:
			hooks.trigger("SplitFastqProgress", {"name" : source.name, "value" : min(1.0, (3/4) + (1/4) * i / len(readsAggregates))})
			threshold += 1
		for file, read in aggregate:
			dataFile.seek(read[1])
			file.write(dataFile.read(read[2] - read[1]))
	hooks.trigger("SplitFastqProgress", {"name" : source.name, "value" : 1.0})

	for files in outFiles:
		for file in files:
			file.close()
	hooks.trigger("SplitFastqFinished", {"name" : source.name, "value" : 3})

	Globals.LOGGER.info(f"Sub Sampling: Finished!")
	
	return outNames
