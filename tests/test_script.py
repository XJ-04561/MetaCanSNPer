
import os

names = ["LVS", "FSC148", "SCHUS4"]

command = (
	"MetaCanSNPer --query SAMPLES/{name}/raw_data/illumina/{name}_1.fastq.gz SAMPLES/{name}/raw_data/illumina/{name}_2.fastq.gz "
	"--organism francisella_tularensis --bytes 10 {size} "
	"--mapperOptions -x sr")

N = 25
BYTES = 2 * (10 **9)

for f in [1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05, 0.025, 0.01]:
	for name in names:
		print(f"Running {name} at {round(100*f)}% ({BYTES*f} bytes)")
		os.system(command.format(name=name, size=round(BYTES*f)))