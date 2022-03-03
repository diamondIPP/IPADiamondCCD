#!/usr/bin/env python
import time
import subprocess as subp
import multiprocessing as mp
import os, shutil, sys
import ipdb
from optparse import OptionParser
from copy import deepcopy
import numpy as np
from collections import OrderedDict

class ParallelManager:
	def __init__(self, working_dir='.', runlist=[], exec_command='ConverterCaen.py', options=[], num_cores=2, force=False, verb=True):
		self.workind_dir = working_dir
		self.exec_command = exec_command
		self.options = options
		self.verb = verb
		self.runlist = runlist
		self.num_cores = num_cores if num_cores <= int(mp.cpu_count()/2.0) or force else int(mp.cpu_count()/2.0)
		self.num_runs = 0
		# self.settings_list = []
		# self.job_chunks = []
		# self.analysis_processes = {}
		# self.workind_dir = os.getcwd()
		self.queue = {}
		self.queue_running = {}
		self.queue_showing = {}
		self.queue_runs = {}
		self.runs_dic_completed = {}
		self.runs_dic_running = {}

	def SetVariables(self, working_dir, runlist, exec_command, options, num_cores=2, force=False, verb=True):
		self.workind_dir = working_dir
		self.exec_command = exec_command
		self.options = options
		self.verb = verb
		self.runlist = runlist
		self.num_cores = num_cores if num_cores <= int(mp.cpu_count() / 2.0) or force else int(mp.cpu_count() / 2.0) if mp.cpu_count() != 1 else 1

	# def RunParallelAnalysis2(self):
	# 	with open(os.devnull, 'w') as FNULL:
	# 		for jobs in self.job_chunks:
	# 			self.analysis_processes = []
	# 			for it, run in enumerate(jobs):
	# 				if it == len(jobs) - 1:
	# 					print 'Showing output for run', run
	# 					self.analysis_processes.append(subp.Popen(['{wd}/rd42Analysis.py'.format(wd=self.workind_dir), '-w', os.getcwd(), '-s', os.path.abspath(run), '--normal'], bufsize=-1, stdin=subp.PIPE, close_fds=True))
	# 				else:
	# 					self.analysis_processes.append(subp.Popen(['{wd}/rd42Analysis.py'.format(wd=self.workind_dir), '-w', os.getcwd(), '-s', os.path.abspath(run), '--normal', '-q'], bufsize=-1, stdin=subp.PIPE, stdout=FNULL, close_fds=True))
	# 			for job_i in xrange(len(self.analysis_processes)):
	# 				while self.analysis_processes[job_i].poll() is None:
	# 					time.sleep(5)
	# 				self.CloseSubprocess(self.analysis_processes[job_i], stdin=True, stdout=False)
	# 			print 'Done with', jobs

	def RunParallelConversion(self):
		os.chdir(self.workind_dir)
		print('\nRunning', len(self.runlist), 'jobs in parallel for', self.num_cores, 'cores\n')
		self.runs_dic_completed = OrderedDict(list(zip(self.runlist, [False for r in self.runlist])))
		self.runs_dic_running = OrderedDict(list(zip(self.runlist, [False for r in self.runlist])))
		self.queue = {c: None for c in range(self.num_cores)}
		self.queue_running = {c: False for c in range(self.num_cores)}
		self.queue_showing = {c: False for c in range(self.num_cores)}
		self.queue_runs = {c: None for c in range(self.num_cores)}
		first_time = True
		# ipdb.set_trace()
		with open(os.devnull, 'w') as FNULL:
			while not np.array(list(self.runs_dic_completed.values()), '?').all():
				pending = np.bitwise_not(np.bitwise_or(list(self.runs_dic_running.values()), list(self.runs_dic_completed.values())))
				pos_run, num_runs_left = pending.argmax(), pending.sum()
				jobi = self.runlist[pos_run]
				option = self.options[pos_run]
				do_add_queue = not np.array(list(self.queue_running.values()), '?').all() and num_runs_left > 0
				if do_add_queue:
					pos_q, nfree = np.array(list(self.queue_running.values()), '?').argmin(), np.bitwise_not(list(self.queue_running.values())).sum()
					print('\nRunning job', jobi, '...')
					command = [self.workind_dir + '/' + self.exec_command] + option
					if nfree == 1 and not np.array(list(self.queue_showing.values()), '?').any() and self.verb:
						print('\nShowing output for job', jobi)
						self.queue[pos_q] = subp.Popen(command, bufsize=-1, stdin=subp.PIPE, close_fds=True)
						# self.queue[pos_q] = 'blaa'
						self.queue_showing[pos_q] = True
					else:
						self.queue[pos_q] = subp.Popen(command, bufsize=-1, stdin=subp.PIPE, stdout=FNULL, stderr=subp.STDOUT, close_fds=True)
						# self.queue[pos_q] = 'baaf'
					self.queue_running[pos_q] = True
					self.runs_dic_running[jobi] = True
					self.queue_runs[pos_q] = jobi
				if not first_time:
					temp = deepcopy(self.queue_running)
					for p, queue_p in temp.items():
						if queue_p:
							if self.queue[p]:
								temp2 = self.queue[p]
								if temp2.poll() is not None:
									jobj = self.queue_runs[p]
									print('\nJob', jobj, 'completed :). Closing ...\r', end=' ')
									self.CloseSubprocess(self.queue[p], stdin=True, stdout=False)
									self.queue[p] = None
									self.queue_running[p] = False
									if self.queue_showing[p]:
										self.queue_showing[p] = False
									self.runs_dic_running[jobj] = False
									self.runs_dic_completed[jobj] = True
									print('Job', jobj, 'completed :). Closing ... Done\n')
					time.sleep(3)
				else:
					first_time = not np.array(list(self.queue_running.values()), '?').all()

	def GetAvailablePos(self):
		pos = -1
		num_free = 0
		for p, elem in self.queue.items():
			if not elem:
				if pos == -1:
					pos = p
				num_free += 1
		return pos, num_free

	def CloseSubprocess(self, p, stdin=False, stdout=False):
		pid = p.pid
		if stdin:
			p.stdin.close()
		if stdout:
			p.stdout.close()
		time.sleep(1)
		if p.wait() is None:
			print('Could not terminate subprocess... forcing termination')
			p.kill()
			time.sleep(1)
			if p.wait() is None:
				print('Could not kill subprocess... quitting')
				exit(os.EX_SOFTWARE)
		try:
			os.kill(pid, 0)
		except OSError:
			pass
		else:
			print('The subprocess is still running. Killing it with os.kill')
			os.kill(pid, 15)
			try:
				os.kill(pid, 0)
			except OSError:
				pass
			else:
				print('The process does not die... quitting')
				exit(os.EX_SOFTWARE)
		del pid
		p = None

def main():
	parser = OptionParser()
	parser.add_option('-d', '--workingdir', dest='workingdir', type='string', help='full path to the working directory e.g. /home/sandiego/rd42AnalysisTools')
	parser.add_option('-r', '--runlist', dest='runlist', type='string', help='string containing the list of jobs to run. e.g. [job1, job2, job3]')
	parser.add_option('-e', '--exec_command', dest='exec_command', type='string', help='command to execute. e.g. rd42Analysis.py')
	parser.add_option('-o', '--options', dest='options', type='string', help='string containing the list with the options to run. e.g. [-w, /home/sandiego/rd42AnalysisTools, --normal]')
	parser.add_option('-n', '--numcores', dest='numcores', type='int', default=2, help='number of runs to execute in parallel')
	parser.add_option('-f', '--force', dest='force', default=False, action='store_true', help='force to use the specified number of cores')
	parser.add_option('-q', '--quiet', dest='quiet', default=False, action='store_true', help='enables quiet mode: no verbose')

	(options, args) = parser.parse_args()
	workingdir = str(options.workingdir)
	runlisttemp = str(options.runlist)
	runlist = runlisttemp.split('[')[1].split(']')[0].split(',')
	execcommand = str(options.exec_command)
	optionstemp = str(options.options)
	optionslist = optionstemp.split('[')[1].split(']')[0].split(',')
	num = int(options.numcores)
	force = bool(options.force)
	verb = not bool(options.quiet)

	pm = ParallelManager(working_dir=workingdir, runlist=runlist, exec_command=execcommand, options=optionslist, num_cores=num, force=force, verb=verb)


if __name__ == '__main__':
	main()
