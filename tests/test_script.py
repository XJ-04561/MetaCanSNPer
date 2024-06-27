
import os

names = ["LVS", "FSC148", "SCHUS4"]

command = (
	"MetaCanSNPer --query ~/SAMPLES/{name}/raw_data/illumina/{name}_1.fastq ~/SAMPLES/{name}/raw_data/illumina/{name}_2.fastq "
	"--organism francisella_tularensis --bytes {n} {size} "
	"--mapperOptions -x sr")

N = 25
BYTES = 2 * (10 **9)

BANNED = [
	("FSC148", round(BYTES*1.0)),
	("SCHUS4", round(BYTES*1.0)),
	("LVS", round(BYTES*0.8)),
	("FSC148", round(BYTES*0.8)),
	("SCHUS4", round(BYTES*0.8)),
	("FSC148", round(BYTES*0.6)),
	("SCHUS4", round(BYTES*0.6)),
	("FSC148", round(BYTES*0.4)),
]

for f in [1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05, 0.025, 0.01]:
	for name in names:
		if (name, round(BYTES*f)) not in BANNED:
			print(f"Running {name} {N} times, at {round(100*f, 1)}% of {BYTES:.0e} ({round(BYTES*f)} bytes)")
			os.system(command.format(n=N, name=name, size=round(BYTES*f)))
		else:
			print(f"Skipping {name} {N} times, at {round(100*f, 1)}% of {BYTES:.0e} ({round(BYTES*f)} bytes)")