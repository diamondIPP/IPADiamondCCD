#!/usr/bin/env python
import os
import shutil
import struct
import subprocess as subp
import sys
import time
from configparser import ConfigParser
from optparse import OptionParser

import ROOT as ro
import numpy as np
import pickle as pickle

from Channel_Caen import Channel_Caen
from Settings_Caen import Settings_Caen
from HV_Control import HV_Control
from Utils import *
from CCD_Caen import CCD_Caen
# from memory_profiler import profile

trig_rand_time = 0.2
wait_time_hv = 7

class VoltageScan:
	def __init__(self, infile='None', vini=5, vend=10, vstep=5, lista=[], listn=[], timebla=1, verbose=False):
		print('Starting CCD program ...')
		self.infile = infile
		self.verb = verbose
		self.vini = vini
		self.vend = vend
		self.vstep = vstep
		self.voltages = np.linspace(self.vini, self.vend, int(round(float(self.vend - self.vini)/float(self.vstep)) + 1), dtype='int32') if lista == [] else np.array(lista, 'int32')
		self.time_sleep = timebla
		self.settings = Settings_Caen(self.infile, self.verb)
		self.settings.ReadInputFile()
		self.settings.Get_Calibration_Constants()
		self.settings.SetOutputFiles()
		self.wd = os.getcwd()

		self.p = {volt: None for volt in self.voltages}
		self.time0 = time.time()
		self.num_events = self.settings.num_events * np.ones(len(self.voltages), 'int32') if listn == [] else np.array(listn, 'int32')
		print(self.num_events)

	def DoVoltageScan(self):

		for it, volt in enumerate(self.voltages):
			self.ResetSettings()
			self.settings.bias = float(volt)
			self.settings.num_events = int(self.num_events[it])
			self.settings.Get_Calibration_Constants()
			self.settings.SetOutputFiles()

			self.p[volt] = CCD_Caen(settingsObj=self.settings)
			self.p[volt].StartHVControl()
			self.p[volt].AdjustBaseLines()
			self.p[volt].SavePickles()
			print('Waiting {s} seconds before starting run...'.format(s=self.time_sleep), end=' ') ; sys.stdout.flush()
			time.sleep(self.time_sleep)
			print('Done')
			written_events = self.p[volt].GetData()
			if self.p[volt].stop_run: print('Run stopped because current is too high')
			self.settings.num_events = written_events
			self.p[volt].SavePickles()
			self.p[volt].settings.MoveBinaryFiles()
			self.p[volt].settings.RenameDigitiserSettings()
			self.p[volt].CloseHVClient()
			if not self.p[volt].settings.simultaneous_conversion:
				self.p[volt].CreateRootFile(files_moved=True)
				while self.p[volt].pconv.poll() is None:
					time.sleep(3)
				self.p[volt].CloseSubprocess('converter', stdin=False, stdout=False)

			print('\nFinished voltage scan for {v} Volts :)\n'.format(v=volt))
			sys.stdout.write('\a\a\a')
			sys.stdout.flush()

			self.p[volt] = None

		print('Finished all voltage scan :D')
		sys.stdout.write('\a\a\a')
		sys.stdout.flush()

	def ResetSettings(self):
		self.settings = Settings_Caen(self.infile, self.verb)
		self.settings.ReadInputFile()
		self.settings.Get_Calibration_Constants()
		self.settings.SetOutputFiles()


def main():
	parser = OptionParser()
	parser.add_option('-i', '--infile', dest='infile', default='', type='string',
	                  help='Input configuration file. e.g. CAENRunConfig.cfg')
	parser.add_option('--start', dest='start', type='int', help='Starting scan voltage', default=0)
	parser.add_option('--stop', dest='stop', type='int', help='Stopping scan voltage', default=0)
	parser.add_option('--step', dest='step', type='int', help='Voltage step between scans', default=0)
	# parser.add_option('--time', dest='time', type='int', help='maximum time in minutes before closing and starting next voltage')
	parser.add_option('-v', '--verbose', dest='verb', default=False, help='Toggles verbose', action='store_true')
	parser.add_option('-a', '--automatic', dest='auto', default=False, help='Toggles automatic conversion and analysis afterwards', action='store_true')
	parser.add_option('-l', '--list', dest='listvals', type='string', help='List of bias voltages to analyse. It should be in brackets. The values should be floats. e.g. [500,550,600]. It overrides "--step, --start, --stop option".', default='[]')
	parser.add_option('-n', '--numevts', dest='numevts', type='string', help='List of number of events to analyse. It should be in brackets. The number of elements has to match the number of bias voltages. e.g. [10000,50000,10000]', default='[]')
	parser.add_option('-t', '--time', dest='timebetween', type='int', help='time in between measurements in seconds.', default=120)

	(options, args) = parser.parse_args()
	infile = str(options.infile)
	auto = bool(options.auto)
	verb = bool(options.verb)
	vini = int(options.start)
	vend = int(options.stop)
	vstep = int(options.step)
	lista = str(options.listvals)
	listn = str(options.numevts)
	print(lista)
	timebla = int(options.timebetween)
	lista1 = []
	if lista != '[]':
		if lista[0] != '[' or lista[-1] != ']':
			ExitMessage('the option -l or --list should start with "[" and end with "]"')
		else:
			lista0 = lista[1:-1].split(',')
			for elem in lista0:
				if IsFloat(elem):
					lista1.append(float(elem))
				else:
					ExitMessage('The entered values in -l or --list are not all integers')
	# tau = int(options.time)
	print(listn)
	listn1 = []
	if listn != '[]':
		if listn[0] != '[' or listn[-1] != ']':
			ExitMessage('the option -n or --numevts should start with "[" and end with "]"')
		else:
			listn0 = listn[1:-1].split(',')
			for elem in listn0:
				if is_int(elem):
					listn1.append(int(elem))
				else:
					ExitMessage('The entered values in -n or --numevts are not all integers')

	if listn1 != '[]' and lista1 != '[]':
		if len(listn1) != len(lista1):
			ExitMessage('The lists have different dimensions. Check')
	vs = VoltageScan(infile, vini, vend, vstep, lista1, listn1, timebla, verb)
	if auto:
		vs.DoVoltageScan()

	return vs

if __name__ == '__main__':
	vs = main()
