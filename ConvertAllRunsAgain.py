#!/usr/bin/env python
import numpy as np
from struct import unpack
import time, os, sys
from optparse import OptionParser
import ipdb
import cPickle as pickle
from Utils import *
from Settings_Caen import Settings_Caen
from Converter_Caen import Converter_Caen
import subprocess as subp
import multiprocessing as mp
from ParallelManager import ParallelManager
# from DataAcquisition import DataAcquisition


class ConvertAllRunsAgain:
	def __init__(self, runsdir='None'):
		self.time0 = time.time()
		self.runsdir = Correct_Path(runsdir)
		runstemp = glob.glob('{d}/*'.format(d=self.runsdir))
		if len(runstemp) < 1:
			ExitMessage('The directory does not have any runs', os.EX_USAGE)
		self.num_cores = 1
		self.runs = []
		self.runs_settings = []
		for run in runstemp:
			settingsfs = glob.glob(run + '/*.settings')
			if len(settingsfs) > 0:
				self.runs.append(run)
				self.runs_settings.append(settingsfs[0])
		self.num_runs = len(self.runs)
		if self.num_runs < 1: ExitMessage('There is not even the required data to convert one run', os.EX_DATAERR)
		self.job_chunks = []
		self.analysis_processes = {}
		self.workind_dir = os.getcwd()
		self.queue = {}
		self.queue_running = {}
		self.queue_showing = {}
		self.queue_runs = {}
		self.runs_dic_completed = {}
		self.runs_dic_running = {}
		self.parallelManager = None

	def DoAll(self, num_cores=2):
		self.time0 = time.time()
		working_dir = os.getcwd()
		self.num_cores = num_cores if num_cores <= int(mp.cpu_count()/2.0) else int(mp.cpu_count()/2.0) if mp.cpu_count() != 1 else 1
		self.job_chunks = [self.runs[i:i + self.num_cores] for i in xrange(0, self.num_runs, self.num_cores)]
		self.num_cores = min(self.num_cores, self.num_runs)
		self.parallelManager = ParallelManager()
		options = [[self.runs_settings[pos], self.runs[pos], '0'] for pos in xrange(self.num_runs)]
		self.parallelManager.SetVariables(working_dir=working_dir, runlist=self.runs, exec_command='Converter_Caen.py', options=options, num_cores=self.num_cores)
		self.parallelManager.RunParallelConversion()
		self.time0 = time.time() - self.time0
		print 'Runs converted in', self.time0, 'seconds,', self.time0/float(self.num_runs + 0.0000001), 'seconds per run. Exiting :D'


def main():
	parser = OptionParser()
	parser.add_option('-d', '--runsdir', dest='runsdir', help='path to folder containing all the run folders to modify')
	(options, args) = parser.parse_args()
	runsdir = str(options.runsdir)
	cara = ConvertAllRunsAgain(runsdir)
	return cara

if __name__ == '__main__':
	cara = main()
