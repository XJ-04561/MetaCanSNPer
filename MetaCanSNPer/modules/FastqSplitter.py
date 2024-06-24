
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

def splitFastq(N : int, filenames : FileList[FilePath], hooks=GlobalHooks, randomiser=None, steps : int=100) -> list[tuple[str]]:
	if all(filename.endswith(".gz") for filename in filenames):
		dataOpen = gunzip.gzip.open
	elif not any(filename.endswith(".gz") for filename in filenames):
		dataOpen = open
	else:
		raise ValueError(f"Files are not consistent in their compression file extensions: {filenames}")
	
	if randomiser is None:
		randomiser = random.Random()
		randomiser.seed(hashlib.md5("".join(filenames).encode("utf-8")).digest())
	
	splitNames = []
	for filepath in filenames:
		name, ext = os.path.basename(filepath)[1:].split(".", 1)
		splitNames.append((os.path.join(os.path.dirname(filepath), os.path.basename(filepath)[0]+name), ext))

	outNames = [tuple(f"{name}-{i+1}-{N}.{ext}" for name, ext in splitNames) for i in range(N)]

	
	dataFiles : list[BinaryIO] = [dataOpen(name, "rb") for name in filenames]

	outFiles = [[dataOpen(filename, "wb") if not os.path.exists(filename) or os.stat(filename).st_size < 1000 else DummyIO() for filename in filenames] for filenames in outNames]

	if all(isinstance(file, DummyIO) for files in outFiles for file in files):
		hooks.trigger("SplitFastqSkipped", {"name" : filenames.name, "value" : 2})
		return outNames
	
	nBytes = sum([file.seek(0, 2) for file in dataFiles])
	for file in dataFiles:
		file.seek(0, 0)

	chooser = equalSampling(N, randomiser)

	hooks.trigger("SplitFastqStarting", {"name" : filenames.name, "value" : 0.0})
	steps : list[float] = [2] + [(steps-i)/steps for i in range(steps+1)]
	traversed = 1
	position = 0
	while traversed > 0:
		
		traversed = 0
		choice = next(chooser)
		for of, df in zip(outFiles[choice], dataFiles):
			traversed += of.write(df.readline())
			traversed += of.write(df.readline())
			traversed += of.write(df.readline())
			traversed += of.write(df.readline())
			
		position += traversed
		if position / nBytes > steps[-1]:
			hooks.trigger("SplitFastqProgress", {"name" : filenames.name, "value" : position / nBytes})
			steps.pop()

	chooser.close()
	for files in outFiles:
		for file in files:
			file.close()
	hooks.trigger("SplitFastqFinished", {"name" : filenames.name, "value" : 3})
	return outNames
