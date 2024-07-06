
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
	
	Globals.LOGGER.info(f"Sub Sampling: From {source}.")
	LOGGING_FILEHANDLER.emit(logging.LogRecord(__name__, logging.INFO, __file__, 0, f"Sub Sampling: From {source}.", (), None, splitFastq))

	if isinstance(source, str):
		source = FileList([FilePath(source)])
	elif isinstance(source, Iterable):
		source = FileList(FilePath(name) for name in source)
	
	if all(filepath.endswith(".gz") for filepath in source):
		Globals.LOGGER.info("Sub Sampling: From and to gzip-data.")
		LOGGING_FILEHANDLER.emit(logging.LogRecord(__name__, logging.INFO, __file__, 0, "Sub Sampling: From and to gzip-data.", (), None, splitFastq))
		dataOpen = gunzip.gzip.open
	elif not any(filepath.endswith(".gz") for filepath in source):
		Globals.LOGGER.info("Sub Sampling: From and to raw-data.")
		LOGGING_FILEHANDLER.emit(logging.LogRecord(__name__, logging.INFO, __file__, 0, "Sub Sampling: From and to raw-data.", (), None, splitFastq))
		dataOpen = open
	else:
		raise ValueError(f"Files are not consistent in their compression file extensions: {source}")
	
	Globals.LOGGER.info("Sub Sampling: Creating Read Index.")
	LOGGING_FILEHANDLER.emit(logging.LogRecord(__name__, logging.INFO, __file__, 0, "Sub Sampling: Creating Read Index.", (), None, splitFastq))
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
				
				if lineSep.strip() == b"+":
					file.seek(sum(map(len, read)), 1)
				elif lineSep:
					file.seek(-len(lineSep), 1)
				readList.append([sum(map(len, map(bytes.strip, read))), pos, file.tell()])

			readsIndex.append(readList)
	Globals.LOGGER.info(f"Sub Sampling: Found {', '.join(map(str, map(len, readsIndex)))} reads for the source file(s).")
	LOGGING_FILEHANDLER.emit(logging.LogRecord(__name__, logging.INFO, __file__, 0, f"Sub Sampling: Found {', '.join(map(str, map(len, readsIndex)))} reads for the source file(s).", (), None, splitFastq))

	for name in ["reads", "dilution", "coverage", "bases", "bytes"]:
		if name in kwargs:
			varName, values = name, kwargs[name]
			Globals.LOGGER.info(f"Sub Sampling: Using {values[0]} {SUB_SAMPLE_NAMES[name]}.")
			LOGGING_FILEHANDLER.emit(logging.LogRecord(__name__, logging.INFO, __file__, 0, f"Sub Sampling: Using {values[0]} {SUB_SAMPLE_NAMES[name]}.", (), None, splitFastq))
			progressTracker = getProgressCallback(name, files, *values, totalReads=len(readsIndex[0]))
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
	# First 1/4 of progress
	Globals.LOGGER.info(f"Sub Sampling: Read Selection.")
	LOGGING_FILEHANDLER.emit(logging.LogRecord(__name__, logging.INFO, __file__, 0, f"Sub Sampling: Read Selection.", (), None, splitFastq))
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
	LOGGING_FILEHANDLER.emit(logging.LogRecord(__name__, logging.INFO, __file__, 0, f"Sub Sampling: Read Aggregation.", (), None, splitFastq))
	readsAggregates = [[] for _ in outFiles[0]]
	for i, data, files in zip(itertools.count(), outData, outFiles):
		if steps * ((1/4) + (1/4) * i / len(outFiles)) >= threshold:
			hooks.trigger("SplitFastqProgress", {"name" : source.name, "value" : min(1.0, (1/4) + (1/4) * i / len(outFiles))})
			threshold += 1
		for fileN, file, reads in zip(itertools.count(), files, data):
			readsAggregates[fileN].extend((file, read) for read in reads)

	# 3/4 of progress
	Globals.LOGGER.info(f"Sub Sampling: Read Sorting.")
	LOGGING_FILEHANDLER.emit(logging.LogRecord(__name__, logging.INFO, __file__, 0, f"Sub Sampling: Read Sorting.", (), None, splitFastq))
	for i, aggregate in enumerate(readsAggregates):
		if steps * ((2/4) + (1/4) * i / len(readsAggregates)) >= threshold:
			hooks.trigger("SplitFastqProgress", {"name" : source.name, "value" : min(1.0, (2/4) + (1/4) * i / len(readsAggregates))})
			threshold += 1
		aggregate.sort(key=lambda x:x[1][1])

	# 4/4 of progress
	Globals.LOGGER.info(f"Sub Sampling: Writing Files.")
	LOGGING_FILEHANDLER.emit(logging.LogRecord(__name__, logging.INFO, __file__, 0, f"Sub Sampling: Writing Files.", (), None, splitFastq))
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
	LOGGING_FILEHANDLER.emit(logging.LogRecord(__name__, logging.INFO, __file__, 0, f"Sub Sampling: Finished!", (), None, splitFastq))
	
	return outNames
