

import os
from logging import Logger
from threading import Thread, Condition
from subprocess import run, DEVNULL, PIPE

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

'''More OOP handling of multiple processes'''
class ThreadGroup:
	# barrier : Barrier
	threads : list[Thread]
	commands : list[str]

	def __init__(self, target : function, args : list[list]=None, kwargs : list[dict]=None, n : int=None):
		# self.barrier = Barrier(2)
		self.target = target

		# Error management for developers
		if args is None and kwargs is None and n is None:
			raise TypeError("""ThreadGroup created without specifying number of threads to be created.
				   The number is hinted from the length of args or kwargs. or specified with n.""")
		elif args is not None and kwargs is not None:
			if len(args) != len(kwargs):
				raise TypeError("""ThreadGroup number of threads is hinted from the length of args or kwargs.
					They are not the same length. (They were {} and {}, respectively)""".format(len(args), len(kwargs)))
			self.args = args
			self.kwargs = kwargs
		elif args is not None:
			self.n = len(args)
			self.args = args
			self.kwargs = [{} for _ in range(self.n)]
		elif kwargs is not None:
			self.n = len(kwargs)
			self.kwargs = kwargs
			self.args = [[] for _ in range(self.n)]
		else:
			self.n = n
			self.args = [[] for _ in range(self.n)]
			self.kwargs = [{} for _ in range(self.n)]
		self.communicator = [[] for _ in range(self.n)]
		self.threads = [Thread(target=self.target, args=args, kwargs=kwargs) for args, kwargs in zip(self.communicator, self.args, self.kwargs)]

	def newFinished(self):
		return any(not t.is_alive() for t in self.threads if t is not None)

	def start(self):
		for t in self.threads:
			t.start()

	def waitNext(self):
		# self.barrier.wait()
		Condition.wait_for(self.newFinished, timeout=5)
		# list comprehension to return all dead threads and also replace those same threads with 'None'.
		return [i for i in range(len(self.threads)) if not self.threads[i].is_alive() if self.threads.__setitem__(i, None) is None]

	def checkDone(self):
		return all(t is None for t in self.threads)

# Aligner and Mapper classes to inherit from
class ProcessWrapper:
	logger : Logger
	softwareName : str
	returncode : int
	
	def __call__(self, *args):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		pass

	def __run__(self, command, log, *args, **kwargs):
		p = run(command.split() if type(command) is str else command, *args, **kwargs)
		with open(log, "W") as logFile:
			logFile.write(p.stdout.read().decode("utf-8"))
			logFile.write("\n")
		self.returncode = p.returncode
		self.handleRetValue(p.returncode)

	def createCommand(self, *args):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		pass
	
	def runCommand(self, *args):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		pass

	def handleRetValue(self, retvalue):
		'''Not implemented for the template class, check the wrapper of the
		specific software you are intending to use.'''
		pass

class IndexingWrapper(ProcessWrapper):
	query : str
	queryName : str
	references : str
	commandTemplate : str # Should contain format tags for {target}, {ref}, {output}, can contain more.
	outFormat : str
	logFormat : str
	returncode : int
	
	def __init__(self, query : str, references : list[str], logger : Logger, outDir : str=".", kwargs : str=""):
		self.query = query
		self.query_name = os.path.basename(query).rsplit(".",1)[0] ## get name of file and remove ending
		self.references = references
		self.outDir = outDir
		self.logger = logger
		self.kwargs = kwargs
		self.returncode = None
	
	def __call__(self, *args):
		''' '''

		commands, logs, outputs = self.createCommand(*args)
		self.runCommand(commands, logs)
		
		return outputs
	
	def createCommand(self) -> tuple[list[str], list[str], list[str]]:
		outputTemplate = "{tmpdir}/{ref}_{target}.{format}"
		logs = []
		commands = []
		outputs = []
		for ref in self.references:
			output = outputTemplate.format(tmpdir=self.outDir, ref=ref, target=self.query_name, format=self.outFormat)
			logfile = outputTemplate.format(tmpdir=self.outDir, ref=ref, target=self.query_name, format=self.logFormat)

			command = self.commandTemplate.format(target=self.query, ref=ref)
			" ".join([command] + ["--{key} {value}".format(key=key, value=value)])
			
			##
			##	
			##

			commands.append(command)
			logs.append(logfile)
			outputs.append(output)
		return commands, logs, outputs



# Aligner, Mapper, and Caller classes to inherit from
class Aligner(IndexingWrapper):
	pass
	
class Mapper(IndexingWrapper):
	pass

class SNPCaller(IndexingWrapper):
	pass

class Mauve(Aligner):
	softwareName = "progressiveMauve"
	commandTemplate = "progressiveMauve {ref} {target}"

	def handleRetCode(self, returncode):
		if returncode == 0:
			self.logger.debug("progressiveMauve finished with exitcode 0.")
		elif returncode == 11:
			self.logger.warning("WARNING progressiveMauve finished with a exitcode: {returncode}".format(returncode=returncode))
			self.logger.debug("This progressiveMauve error is showing up for bad genomes containing short repetitive contigs or sequence contains dashes.".format(returncode=returncode))
		elif returncode == -6:
			self.logger.warning("Input sequence is not free of gaps, replace gaps with N and retry!!")
		else:
			self.logger.warning("WARNING progressiveMauve finished with a non zero exitcode: {returncode}\nThe script will terminate when all processes are finished read {log} for more info".format(log=self.logFile,returncode=returncode))
		
"""
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
						self.logger.debug("This progressiveMauve error is showing up for bad genomes containing short repetitive contigs or sequence contains dashes.".format(exitcode=exitcode))
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
	"""