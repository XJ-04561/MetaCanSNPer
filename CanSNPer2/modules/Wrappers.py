

import os
from logging import Logger

## import standard python libraries for subprocess and multiprocess
from subprocess import Popen,PIPE,STDOUT

class Error(Exception):
	"""docstring for Error"""
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class MauveError(Error):
	"""docstring for MauveE"""
	pass

'''Declarations used only for typehinting'''
class CanSNPer2:
	logdir : str
	refdir : str
	get_references : function
	keep_going : bool
class Aligner:
	tmpdir : str
	query : str

# Aligner and Mapper classes to inherit from
class Aligner:
	query : str
	queryName : str
	references : str
	logger : Logger
	softwareName : str
	
	def __init__(self, softwareName, **kwargs):
		
		self.softwareName

		for key, item in kwargs.items():
			self.__setattr__(key, item)
	
	def __call__(self, *args):
		''' '''

		commands, logs = self.createCommand(*args)
		self.__run__(commands, logs)

	def __run__(self, commands : list[str], logs : list[str]):
		retvalue=0
		processes = []  #process container
		log_f = {}	  #log container
		error = 0   #Variable for errors
		self.logger.info("Starting {name} on {n} references".format(name=self.softwareName, n=len(commands)))
		for command, log in zip(commands, logs): #Loop through commands,
			self.logger.debug(command)			## In verbose mode print the actual mauve command
			p = Popen(command.split(" "),  stdout=PIPE, stderr=STDOUT)  ##Split command to avoid shell=True pipe stdout and stderr to stdout
			log_f[p.stdout.fileno()] = log	  ## Store the reference to the correct log file
			processes.append(p)
		while processes:
			for p in processes: ## Loop through processes
				exitcode = p.poll()  #get exitcode for each process
				if exitcode is not None: ## if there is no exitcode the program is still running

					## When a process is finished open the log file and write stdout/stderr to log file
					with open(log_f[p.stdout.fileno()], "w") as log:
						print(p.stdout.read().decode('utf-8'),file=log)
					### IF the exitcode is not 0 print a warning and ask user to read potential error messages
					if exitcode == 11:
						self.logger.warning("WARNING progressiveMauve finished with a exitcode: {exitcode}".format(exitcode=exitcode))
						self.logger.debug("This progressiveMauve error is showing up for bad genomes containing short repetative contigs or sequence contains dashes".format(exitcode=exitcode))
						retvalue=11
					elif exitcode == -6:
						self.logger.warning("Input sequence is not free of gaps, replace gaps with N and retry!!")
						retvalue=6
					elif exitcode != 0:
						if not self.keep_going:
							self.logger.error("Error: exitcode-{exitcode}".format(exitcode=exitcode))
						error += 1
						self.logger.warning("WARNING progressiveMauve finished with a non zero exitcode: {exitcode}\nThe script will terminate when all processes are finished read {log} for more info".format(log=log_f[p.stdout.fileno()],exitcode=exitcode))
					## Remove process from container
					processes.remove(p)
		if retvalue == 6 or retvalue == 11: ## This is done if sequence is not free of gaps
			return retvalue
		elif not self.keep_going:
			if error:  ## Error handling regarding mauve subprocesses, stop script if any of them fails
				self.logger.error("Error: {errors} progressiveMauve processes did not run correctly check log files for more information".format(errors=error))
				raise MauveError("Error: {errors} progressiveMauve processes did not run correctly check log files for more information".format(errors=error))
		else:
			retvalue=1
		return retvalue

	def createCommand(self, *args):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		pass

	
class Mapper:
	def __init__():

class Mauve:
	query : str
	queryName : str
	references : list[str]
	logger : Logger
	
	def __init__(self, query, logger : Logger, references=[]):
		self.query = query
		self.query_name = os.path.basename(query).rsplit(".",1)[0] ## get name of file and remove ending
		self.references = references
		self.logger = logger
	
    def __call__(self, query=None, references=None):
        

	def create_command(self, aligner : Aligner | CanSNPer2, query : str=None, references : list[str]=None):
		if query is None:
			query = self.query
			query_name = self.query_name
		else:
			query_name = os.path.basename(query).rsplit(".", 1)[0] ## get name of file and remove ending
			
		if references is None:
			references = self.references
		
		'''Mauve commands'''
		commands =[]	# store execute command
		logs = []	   # store log filepath
		if len(references) == 0:	## If specific references are not given fetch references from the reference folder
			references : list[str] = aligner.get_references()
		for ref in references:	  ## For each reference in the reference folder align to query
			ref_name = ref.rsplit(".",1)[0] ## remove file ending
			
			xmfa_output = "{tmpdir}/{ref}_{target}.xmfa".format(tmpdir=aligner.tmpdir.rstrip("/"),ref=ref_name,target=query_name)
			ref_file = "{refdir}/{ref}".format(refdir=aligner.refdir, ref=ref)
			log_file = "{logdir}/{ref}_{target}.mauve.log".format(logdir=aligner.logdir,ref=ref_name,target=query_name)

			'''Create run command for mauve'''
			command = "{mauve_path}progressiveMauve --output {xmfa} {ref_fasta} {target_fasta}".format(
							mauve_path	  = aligner.mauve_path,
							xmfa			= xmfa_output,
							ref_fasta	   = ref_file,
							target_fasta	= query
			)
			commands.append(command)			## Mauve command
			logs.append(log_file)			   ## Store log files for each alignment
			aligner.xmfa_files.append(xmfa_output) ## Store the path to xmfa files as they will be used later
		return commands,logs
	
	def run_mauve(self, commands : list[str], logs : list[str]) -> int:
		'''Run mauve
			ProgressiveMauve is an alignment program for small genomes
			this function will execute mauve commands as paralell subprocesses

			To run subprocesses securely this functino uses the Popen command includint standard pipes
			stderr will be passed to the stdout pipe to reduce the number of log files required

			This function currently does not support any limitation of number of processes spawned,
			modern OS will however mostly distribute subprocesses efficiently, it will also be limited by the
			number of references supplied. A warning message will be given if the number of input references
			exceeds the number of references in the database
		'''
		retvalue=0
		processes = []  #process container
		log_f = {}	  #log container
		error = 0   #Variable for errors
		self.logger.info("Starting progressiveMauve on {n} references".format(n=len(commands)))
		for i in range(len(commands)): #Loop through commands,
			command = commands[i]
			self.logger.debug(command)			## In verbose mode print the actual mauve command
			p = Popen(command.split(" "),  stdout=PIPE, stderr=STDOUT)  ##Split command to avoid shell=True pipe stdout and stderr to stdout
			log_f[p.stdout.fileno()] = logs[i]	  ## Store the reference to the correct log file
			processes.append(p)
		while processes:
			for p in processes: ## Loop through processes
				exitcode = p.poll()  #get exitcode for each process
				if exitcode is not None: ## if there is no exitcode the program is still running

					## When a process is finished open the log file and write stdout/stderr to log file
					with open(log_f[p.stdout.fileno()], "w") as log:
						print(p.stdout.read().decode('utf-8'),file=log)
					### IF the exitcode is not 0 print a warning and ask user to read potential error messages
					if exitcode == 11:
						self.logger.warning("WARNING progressiveMauve finished with a exitcode: {exitcode}".format(exitcode=exitcode))
						self.logger.debug("This progressiveMauve error is showing up for bad genomes containing short repetative contigs or sequence contains dashes".format(exitcode=exitcode))
						retvalue=11
					elif exitcode == -6:
						self.logger.warning("Input sequence is not free of gaps, replace gaps with N and retry!!")
						retvalue=6
					elif exitcode != 0:
						if not self.keep_going:
							self.logger.error("Error: exitcode-{exitcode}".format(exitcode=exitcode))
						error += 1
						self.logger.warning("WARNING progressiveMauve finished with a non zero exitcode: {exitcode}\nThe script will terminate when all processes are finished read {log} for more info".format(log=log_f[p.stdout.fileno()],exitcode=exitcode))
					## Remove process from container
					processes.remove(p)
		if retvalue == 6 or retvalue == 11: ## This is done if sequence is not free of gaps
			return retvalue
		elif not self.keep_going:
			if error:  ## Error handling regarding mauve subprocesses, stop script if any of them fails
				self.logger.error("Error: {errors} progressiveMauve processes did not run correctly check log files for more information".format(errors=error))
				raise MauveError("Error: {errors} progressiveMauve processes did not run correctly check log files for more information".format(errors=error))
		else:
			retvalue=1
		return retvalue