# MetaCanSNPer
MetaCanSNPer: A toolkit for SNP-typing NGS data.
ParseXMFA2: a script/module for looking up positions in an aligned genome from an .XMFA file, without loading the file into memory.

Future planned implementations
* Allow the next process to start when the previous has finished, without waiting for sibling processes to finish.
* Multiple SNPs per node.
* GUI
* Interactive Graph outputs
* Allow for more complex SNP Arbitration (Possibly weighting SNPs)
* Pre-Processing command to determina organisms present in read-data.

Databases supplied can be found at https://github.com/FOI-Bioinformatics/CanSNPer2-data

* Francisella tularensis
* Bacillus anthracis
* Yersinia pestis
 
#### Upcoming databases
* Brucella
* Coxiella


## Installation
As of right now it is only possible to install directly from the repository
```
python setup.py install
```

## Requirements (for manual install conda will install all dependencies)
None.\\
But at least an aligner/mapper and a SNPCaller is needed. (XMFA can be parsed without external software.)

### Supported Aligners
* progressiveMauve

### Supported Mappers
* minimap2

### Supported SNPCallers
* gatk Mutect2

## User guide CanSNPer2 (for custom databases see below)
1. Download pre-built databases from https://github.com/FOI-Bioinformatics/CanSNPer2-data

2. Run genomes
```sh
MetaCanSNPer --query FILE1 [FILE2] --database DOWNLOADED_DATABASE --mapper MAPPER_COMMAND --SNPCaller SNPCALLER_COMMAND
```

For more options MetaCanSNPer --help

## Help Page
```
usage: CommandLineParser.py [-h] [--list] [--query query [query ...]] [-d database] (--mapper mapper | --aligner aligner) [--snpCaller snpCaller] [-s saveTemp]
                            [--settingsFile settingsFile] [--mapperOptions Mapper options] [--alignerOptions Aligner options] [--snpCallerOptions SNP Caller options] [-W DIR]        
                            [-U DIR] [-I DIR] [-Q DIR] [-T DIR] [-R DIR] [-D DIR] [-O DIR] [-S DIR] [--verbose] [--debug] [--supress]

MetaCanSNPer

options:
  -h, --help            show this help message and exit
  --list                To list implemented software and exit.

Required arguments:
  --query query [query ...]
                        Raw sequence data file supported by the intended Aligner/Mapper.
  -d database, --database database
                        Filename of CanSNP database to be used.
  --mapper mapper       Name of installed and supported Mapper software.
  --aligner aligner     Name of installed and supported Alignment software.
  --snpCaller snpCaller
                        Name of installed and supported SNP Calling software.

Optional arguments:
  -s saveTemp, --saveTemp saveTemp
                        Path to .TOML file containing settings for MetaCanSNPer. Check the 'defaultConfig.toml' to see what can be included in a settings file.
  --settingsFile settingsFile
                        Path to .TOML file containing settings for MetaCanSNPer. Check the 'defaultConfig.toml' to see what can be included in a settings file.
  --mapperOptions Mapper options
                        To provide flags/arguments for the chosen Mapper, provide them directly after the '--indexerOptions' flag, only interrupted by the end of the command call    
                        or the corresponding flag for an Aligner or SNP Caller options.
  --alignerOptions Aligner options
                        To provide flags/arguments for the chosen Aligner, provide them directly after the '--indexerOptions' flag, only interrupted by the end of the command call   
                        or the corresponding flag for a Mapper or SNP Caller options.
  --snpCallerOptions SNP Caller options
                        To provide flags/arguments for the chosen SNP Caller, provide them directly after the '--snpCallerOptions' flag, only interrupted by the end of the command   
                        call or the corresponding flag for Mapper or Aligner options.

Directory Options:
  -W DIR, --workDir DIR
                        Work directory
  -U DIR, --userDir DIR
                        User directory
  -I DIR, --installDir DIR
                        Installation directory
  -Q DIR, --targetDir DIR
                        Target (Query) directory
  -T DIR, --tmpDir DIR  Temporary directory
  -R DIR, --refDir DIR  References directory
  -D DIR, --databaseDir DIR
                        Databases directory
  -O DIR, --outDir DIR  Output directory
  -S DIR, --sessionName DIR
                        Session Name/Directory

Logging and debug options:
  --verbose             Verbose output
  --debug               Debug output
  --supress             Supress warnings
```

## About

About this software
===================
Copyright (C) 2024 Fredrik Sörensen @ Umeå University

MetaCanSNPer is implemented as a Python 3 package. It is open source software made available
under the [GPL-3.0 license](LICENSE).

If you experience any difficulties with this software, or you have suggestions, or want
to contribute directly, you have the following options:

- submit a bug report or feature request to the
  [issue tracker](https://github.com/XJ-04561/MetaCanSNPer/issues)
- contribute directly to the source code through the
  [github](https://github.com/XJ-04561/MetaCanSNPer) repository. 'Pull requests' are
  especially welcome.
