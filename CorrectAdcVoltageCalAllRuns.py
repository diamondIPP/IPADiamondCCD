#!/usr/bin/env python
from optparse import OptionParser
import glob
import shutil
from Utils import *
from Modify_Settings_Caen import Modify_Pickles_Caen


class CorrectAdcVoltageCalAllRuns:
	def __init__(self, runsdir='None'):
		self.runsdir = Correct_Path(runsdir)
		self.runs = glob.glob('{d}/*'.format(d=self.runsdir))

	def MakeAllCalibrationModifications(self, vcaldir=''):
		for run in self.runs:
			if os.path.isdir(run):
				print('Modifying run', run)
				print('Backing up old pickles...', end=' ') ; sys.stdout.flush()
				settingsfs = glob.glob(run + '/*.settings')
				signalsfs =  glob.glob(run + '/*.signal_ch')
				triggersfs =  glob.glob(run + '/*.trigger_ch')
				vetosfs =  glob.glob(run + '/*.veto')
				if len(settingsfs) > 0:
					shutil.copy2(settingsfs[0], settingsfs[0] + '.bkp')
				if len(signalsfs) > 0:
					shutil.copy2(signalsfs[0], signalsfs[0] + '.bkp')
				if len(triggersfs) > 0:
					shutil.copy2(triggersfs[0], triggersfs[0] + '.bkp')
				if len(vetosfs) > 0:
					shutil.copy2(vetosfs[0], vetosfs[0] + '.bkp')
				print('Done')
				modPickleRun = Modify_Pickles_Caen(run)
				modPickleRun.CorrectAdcVoltageCal(vcaldir)
				modPickleRun.SavePickles()
				print('Finish with run', run, '\n')

def main():
	parser = OptionParser()
	parser.add_option('-d', '--runsdir', dest='runsdir', help='path to folder containing all the run folders to modify')
	(options, args) = parser.parse_args()
	runsdir = str(options.runsdir)
	mp = CorrectAdcVoltageCalAllRuns(runsdir)
	return mp

if __name__ == '__main__':
	mp = main()
