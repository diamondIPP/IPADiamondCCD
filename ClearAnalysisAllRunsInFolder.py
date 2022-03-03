#!/usr/bin/env python
import numpy as np
from struct import unpack
import time, os, sys
from optparse import OptionParser
import ipdb
import pickle as pickle
from Utils import *
from Settings_Caen import Settings_Caen
from ConverterCaen import ConverterCaen
import subprocess as subp
from AnalysisCaenCCD import AnalysisCaenCCD
import ROOT as ro

fit_method = ('Minuit2', 'Migrad', )

class ClearAnalysisAllRunsInFolder:
	def __init__(self, runsdir=''):
		self.time0 = time.time()
		self.runsdir = Correct_Path(runsdir)
		runstemp = glob.glob('{d}/*'.format(d=self.runsdir))
		if len(runstemp) < 1:
			ExitMessage('The directory does not have any runs', os.EX_USAGE)
		# self.num_cores = 1
		self.runs = [runi for runi in runstemp if os.path.isdir(runi)]
		self.runs.sort(key=lambda x: float(x.split('_Pos_')[-1].split('V')[0]) if '_Pos_' in x else -1*float(x.split('_Neg_')[-1].split('V')[0]))
		self.num_runs = len(self.runs)
		if self.num_runs < 1: ExitMessage('There is not even the required data to convert one run', os.EX_DATAERR)

	def DoAll(self):
		self.time0 = time.time()
		self.LoopRuns()
		print('Finished in', time.time() - self.time0, 'seconds')

	def LoopRuns(self):
		# ro.gROOT.SetBatch(True)
		for run in self.runs:
			print('Clearing Analysis in run:', run)
			if os.path.isdir(run):
				root_files = glob.glob('{d}/*.root'.format(d=run))
				if len(root_files) > 0:
					analysis_root_files = glob.glob('{d}/*.analysis*root'.format(d=run))

					if len(analysis_root_files) > 0:
						for anaf in analysis_root_files:
							print('Removing file: ', anaf)
							os.remove(anaf)


def main():
	parser = OptionParser()
	parser.add_option('-d', '--runsdir', dest='runsdir', type='str', default='', help='path to folder containing all the run folders to modify')
	(options, args) = parser.parse_args()
	runsdir = str(options.runsdir)
	carif = ClearAnalysisAllRunsInFolder(runsdir)
	return carif

if __name__ == '__main__':
	carif = main()
