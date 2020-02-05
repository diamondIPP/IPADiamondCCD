#!/usr/bin/env python
import numpy as np
from struct import unpack
import time, os, sys
from optparse import OptionParser
import ipdb
import cPickle as pickle
from Utils import *
from Settings_Caen import Settings_Caen
# from DataAcquisition import DataAcquisition


class Modify_Pickles_Caen:
	def __init__(self, outdir='None'):
		self.outdir = Correct_Path(outdir)
		self.settings_path = glob.glob('{d}/*.settings'.format(d=outdir))
		self.signalCh_path = glob.glob('{d}/*.signal_ch'.format(d=outdir))
		self.triggerCh_path = glob.glob('{d}/*.trigger_ch'.format(d=outdir))
		self.vetoCh_path = glob.glob('{d}/*.veto'.format(d=outdir))

		self.settings = pickle.load(open(self.settings_path[0], 'rb')) if len(self.settings_path) == 1 else None
		self.signal_ch = pickle.load(open(self.signalCh_path[0], 'rb')) if len(self.signalCh_path) == 1 else None
		self.trigger_ch = pickle.load(open(self.triggerCh_path[0], 'rb')) if len(self.triggerCh_path) == 1 else None
		self.veto_ch = pickle.load(open(self.vetoCh_path[0], 'rb')) if len(self.vetoCh_path) == 1 else None

	def SavePickles(self):
		if len(self.settings_path) == 1:
			with open(self.settings_path[0], 'wb') as fs:
				pickle.dump(self.settings, fs, pickle.HIGHEST_PROTOCOL)
		if len(self.signalCh_path) == 1:
			with open(self.signalCh_path[0], 'wb') as fs:
				pickle.dump(self.signal_ch, fs, pickle.HIGHEST_PROTOCOL)
		if len(self.triggerCh_path) == 1:
			with open(self.triggerCh_path[0], 'wb') as fs:
				pickle.dump(self.trigger_ch, fs, pickle.HIGHEST_PROTOCOL)
		if len(self.vetoCh_path) == 1:
			with open(self.vetoCh_path[0], 'wb') as fs:
				pickle.dump(self.veto_ch, fs, pickle.HIGHEST_PROTOCOL)

	def ModifyVetoThreshold(self, val):
		self.settings.ac_thr_counts = val
		self.veto_ch.thr_counts = val


def main():
	parser = OptionParser()
	parser.add_option('-d', '--outdir', dest='outdir', help='path to folder containing the pickles')
	(options, args) = parser.parse_args()
	outdir = str(options.outdir)
	mp = Modify_Pickles_Caen(outdir)
	return mp

if __name__ == '__main__':
	mp = main()
