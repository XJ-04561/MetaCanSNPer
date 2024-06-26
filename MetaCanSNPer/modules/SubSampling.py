
from MetaCanSNPer.Globals import *
from MetaCanSNPer.core.Hooks import *
import hashlib
import gunzip
		

SUB_SAMPLE_NAMES = {
	"reads" : "Reads",
	"coverage" : "Coverage",
	"dilution" : "Dilution",
	"bytes" : "Bytes"
}

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

def subSampleName(name : FilePath|DirectoryPath|str, type : Literal["reads","coverage","dilution","bytes"], *N, index : int=None) -> str:
	
	if index is None:
		bracketedID = f"[{SUB_SAMPLE_NAMES[type]}-{'-'.join(map(shortNumber, N))}]"
	else:
		bracketedID = f"[{SUB_SAMPLE_NAMES[type]}-{str(index).zfill(len(str(N[0])))}-{'-'.join(map(shortNumber, N))}]"
	if isinstance(name, FilePath):
		return name.name + bracketedID + name.ext
	elif isinstance(name, DirectoryPath):
		return os.path.basename(name) + bracketedID
	else:
		name, *ext = os.path.basename(name)[1:].split(".")
		name, ext = os.path.basename(name)[0]+name, ".".join(ext)
		if ext:
			ext = "." + ext
		
		return name + bracketedID + ext
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
			   reads : list[int]=None, dilution : list[int]=None, coverage : list[int,int]=None, bytes : list[int]=None,
			   outDir : DirectoryPath=None, hooks=GlobalHooks, randomiser=None, steps : int=100) -> list[tuple[str]]: ...
def splitFastq(files : int, source : FilePath|FileList[FilePath], *,
			   reads : list[int]|None=None, dilution : list[int]|None=None, coverage : list[int,int]|None=None, bytes : list[int]|None=None,
			   outDir : DirectoryPath=None, hooks=GlobalHooks, randomiser=None, steps : int=100) -> list[tuple[str]]:
	if not isinstance(source, list):
		source = FileList([source])
	elif isinstance(source, FileList):
		source = FileList(source)
	
	if all(filepath.endswith(".gz") for filepath in source):
		dataOpen = gunzip.gzip.open
	elif not any(filepath.endswith(".gz") for filepath in source):
		dataOpen = open
	else:
		raise ValueError(f"Files are not consistent in their compression file extensions: {source}")
	
	if randomiser is None:
		randomiser = random.Random()
		randomiser.seed(hashlib.md5("".join(source).encode("utf-8")).digest())
	
	if reads is not None:
		splitByReads(files, reads, source, dataOpen, randomiser, outDir=outDir, hooks=hooks)
	elif coverage is not None:
		splitByCoverage(files, coverage, source, dataOpen, randomiser, outDir=outDir, hooks=hooks)
	elif dilution is not None:
		splitByDilution(files, dilution, source, dataOpen, randomiser, outDir=outDir, hooks=hooks)
	else:
		splitByBytes(files, bytes, source, dataOpen, randomiser, outDir=outDir, hooks=hooks)

@overload
def splitByReads(files : int, reads : int, source : FilePath, dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], outDir : DirectoryGroup|None=None, hooks : Hooks=GlobalHooks, steps : int=100): ...
@overload
def splitByReads(files : int, reads : int, source : FileList[FilePath], dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], outDir : DirectoryGroup|None=None, hooks : Hooks=GlobalHooks, steps : int=100): ...
def splitByReads(files : int, reads : int, source : FilePath|FileList[FilePath], dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], outDir : DirectoryGroup|None=None, hooks : Hooks=GlobalHooks, steps : int=100):
	
	
	hooks.trigger("SplitFastqStarting", {"name" : source.name, "value" : 0.0})
	outNames = [tuple((outDir or filepath.directory) / subSampleName(filepath, "reads", files, reads, index=i+1) for filepath in source) for i in range(files)]

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
def splitByCoverage(files : int, coverage : int, source : FilePath, dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], outDir : DirectoryGroup|None=None, hooks : Hooks=GlobalHooks, steps : int=100): ...
@overload
def splitByCoverage(files : int, coverage : int, source : FileList[FilePath], dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], outDir : DirectoryGroup|None=None, hooks : Hooks=GlobalHooks, steps : int=100): ...
def splitByCoverage(files : int, coverage : int, source : FilePath|FileList[FilePath], dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], outDir : DirectoryGroup|None=None, hooks : Hooks=GlobalHooks, steps : int=100):

	raise NotImplementedError(f"Not yet implemented!")

	hooks.trigger("SplitFastqStarting", {"name" : source.name, "value" : 0.0})
	outNames = [tuple((outDir or filepath.directory) / subSampleName(filepath, "coverage", files, coverage, index=i+1) for filepath in source) for i in range(files)]

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
def splitByDilution(files : int, dilutionFactor : int, source : FilePath, dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], outDir : DirectoryGroup|None=None, hooks : Hooks=GlobalHooks, steps : int=100): ...
@overload
def splitByDilution(files : int, dilutionFactor : int, source : FileList[FilePath], dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], outDir : DirectoryGroup|None=None, hooks : Hooks=GlobalHooks, steps : int=100): ...
def splitByDilution(files : int, dilutionFactor : int, source : FilePath|FileList[FilePath], dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], outDir : DirectoryGroup|None=None, hooks : Hooks=GlobalHooks, steps : int=100):
	
	hooks.trigger("SplitFastqStarting", {"name" : source.name, "value" : 0.0})
	outNames = [tuple((outDir or filepath.directory) / subSampleName(filepath, "dilution", files, dilutionFactor, index=i+1) for filepath in source) for i in range(files)]

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

@overload
def splitByBytes(files : int, bytesPerFile : int, source : FilePath, dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], outDir : DirectoryGroup|None=None, hooks : Hooks=GlobalHooks, steps : int=100): ...
@overload
def splitByBytes(files : int, bytesPerFile : int, source : FileList[FilePath], dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], outDir : DirectoryGroup|None=None, hooks : Hooks=GlobalHooks, steps : int=100): ...
def splitByBytes(files : int, bytesPerFile : int, source : FilePath|FileList[FilePath], dataOpen : Callable[[FilePath, str],TextIO], randomiser : Iterator[int], outDir : DirectoryGroup|None=None, hooks : Hooks=GlobalHooks, steps : int=100):
	
	hooks.trigger("SplitFastqStarting", {"name" : source.name, "value" : 0.0})
	outNames = [tuple((outDir or filepath.directory) / subSampleName(filepath, "bytes", files, bytesPerFile, index=i+1) for filepath in source) for i in range(files)]

	if all(os.path.exists(filename) for filenames in outNames for filename in filenames):
		hooks.trigger("SplitFastqSkipped", {"name" : source.name, "value" : 2})
		return outNames
	
	dataFiles : list[BinaryIO] = [dataOpen(name, "rb") for name in source]
	outFiles = [[dataOpen(filename, "wb") for filename in filenames] for filenames in outNames]
	
	bytesWritten = [0 for _ in range(len(outFiles))]
	
	threshold = 0
	for choice in equalSampling(files, randomiser):
		if bytesPerFile < sum(bytesWritten) / files:
			break
		while abs(choice) < len(outFiles)+1: # choice is already satisfied, select next.
			choice -= 1
			if bytesWritten[choice] < bytesPerFile:
				break # choice now selects a file which needs more data.
		else:
			break # All files have the correct size.

		traversed = sum(of.write(df.readline()) for of, df in zip(outFiles[choice], dataFiles) for _ in range(4))
		
		bytesWritten[choice] += traversed

		if traversed == 0: # End of File reached, time to restart
			for df in dataFiles:
				df.seek(0)
		
		if (sum(bytesWritten) / files) / bytesPerFile > threshold:
			hooks.trigger("SplitFastqProgress", {"name" : source.name, "value" : min(1.0, (sum(bytesWritten) / files) / bytesPerFile)})
			threshold = (int(steps * ((sum(bytesWritten) / files) / bytesPerFile)) + 1) / steps

	for files in outFiles:
		for file in files:
			file.close()
	hooks.trigger("SplitFastqFinished", {"name" : source.name, "value" : 3})
	return outNames