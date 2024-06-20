
from MetaCanSNPer.Globals import *
from MetaCanSNPer.core.Hooks import *
import hashlib
import gunzip

def equalSampling(N, randomiser : random.Random):
	counts = [1 for _ in range(N)]
	weights = [0 for _ in range(N)]
	choices = list(range(N))
	while True:
		for i in range(N):
			weights[i] = 1 - (counts[i] / sum(counts))
		for i in randomiser.choices(choices, weights, k=100):
			yield i

def splitFastq(N : int, filenames : FileList, hooks=GlobalHooks, randomiser=None) -> list[tuple[str]]:
	if randomiser is None:
		randomiser = random.Random()
		randomiser.seed(hashlib.md5("".join(filenames)).digest())
	
	nBytes = sum([os.stat(name).st_size for name in filenames])
	dataFiles : list[TextIO] = [gunzip.gzip.open(name, "r") for name in filenames]

	splitNames = [filename[0]+filename[1:].split(".", 1) for filename in filenames]

	outNames = [tuple("".join([name, f".{i+1}.{N}", ext]) for name, ext in splitNames) for i in range(N)]
	outFiles = [[gunzip.gzip.open(filename, "w") for filename in filenames] for filenames in outNames]

	chooser = equalSampling(N, randomiser)

	hooks.trigger("SplitFastqStarting", {"name" : filenames.name, "value" : 0.0})
	traversed = 1
	position = 0
	progress = 0.5
	finalProg = nBytes // 100
	while traversed > 0:
		traversed = 0
		choice = next(chooser)
		for of, df in zip(outFiles[choice], dataFiles):
			entry = df.readline()+df.readline()+df.readline()+df.readline()
			of.write(entry)
			traversed += len(entry)
		position += traversed
		if 0.5 + position // finalProg > progress:
			progress = position // finalProg
			hooks.trigger("SplitFastqProgress", {"name" : filenames.name, "value" : progress / 100})

	chooser.close()
	for files in outFiles:
		for file in files:
			file.close()
	hooks.trigger("SplitFastqFinished", {"name" : filenames.name, "value" : 3})
	return outNames
