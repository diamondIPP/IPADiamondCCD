#!/usr/bin/env python
import numpy as np
from struct import unpack
import time, os, sys
from optparse import OptionParser
import ipdb
import pickle as pickle
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

		self.voltage_cal_dir = ''

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

	def CorrectAdcVoltageCal(self, vcaldir=''):
		def CheckForDir():
			while self.voltage_cal_dir == '' or not os.path.isdir(self.voltage_cal_dir):
				prompt_text = 'The path for the adc-voltage-cal directory does not exist. Enter a correct path to the directory containing the adc_cal_<ch>_<dcOffsetP>.cal files: ' if self.voltage_cal_dir != '' else 'The adc voltage calibration directory is not present in this run. Please enter it: '
				self.voltage_cal_dir = input(prompt_text)
				self.voltage_cal_dir = Correct_Path(self.voltage_cal_dir)
			print('Set voltage calibration directory:', self.voltage_cal_dir)

		def Get_voltage_calibration():
			if self.voltage_cal_dir != '':
				if os.path.isdir(self.voltage_cal_dir):
					v_adc_calfile = ''
					if 'sig_dc_offset_percent' in list(self.settings.__dict__.keys()):
						v_adc_calfile = 'adc_cal_{c}_{dop}.cal'.format(c=self.settings.sigCh, dop=self.settings.sig_dc_offset_percent)
					else:
						v_adc_calfile = 'adc_cal_{c}_{dop}.cal'.format(c=self.settings.sigCh, dop=45 if self.settings.bias < 0 else -45)
					if os.path.isfile('{d}/{f}'.format(d=self.voltage_cal_dir, f=v_adc_calfile)):
						self.settings.adc_volts_pickle = pickle.load(open('{d}/{f}'.format(d=self.voltage_cal_dir, f=v_adc_calfile), 'rb'))
					else:
						print('Using system default (with: DC Offset Percent=0)')
						self.settings.adc_volts_pickle = {'fit_p0': self.settings.sigRes, 'fit_p1': -1.0}
					print('Using p0:', self.settings.adc_volts_pickle['fit_p0'], 'and p1:', self.settings.adc_volts_pickle['fit_p1'])
				else:
					ExitMessage('Given directory does not exist: {d}'.format(d=self.voltage_cal_dir))
			else:
				ExitMessage('No directory has been se to get the voltage calibrations')

		def Set_Channel(channelObj):
			if 'offseted_adc_to_volts_cal' in list(channelObj.__dict__.keys()):
				channelObj.adc_to_volts_cal = {}
			if 'adc_to_volts_cal' in list(channelObj.__dict__.keys()):
				if channelObj.type in ['signal', 'signal_ch']:
					channelObj.adc_to_volts_cal['p0'] = self.settings.adc_volts_pickle['fit_p0']
					channelObj.adc_to_volts_cal['p1'] = self.settings.adc_volts_pickle['fit_p1']
					print('Set signal channel p0 to:', channelObj.adc_to_volts_cal['p0'], 'and p1 to:', channelObj.adc_to_volts_cal['p1'])
				else:
					if not 'dc_offset_percent' in list(channelObj.__dict__.keys()):
						print('No dc_offset_percent found in the channel object. As it is not signal, assuming it is -45')
						channelObj.dc_offset_percent = -45
					channelObj.adc_to_volts_cal['p1'] = self.settings.sigRes
					channelObj.adc_to_volts_cal['p0'] = np.subtract(np.divide(channelObj.dc_offset_percent, 50, dtype='f8'), 1, dtype='f8')
					channelObj.base_line_u_adcs = np.divide(0, channelObj.adc_to_volts_cal['p1'])
					channelObj.base_line_adcs = RoundInt(np.divide(np.subtract(0, channelObj.adc_to_volts_cal['p0'], dtype='f8'), channelObj.adc_to_volts_cal['p1'], dtype='f8'), 'uint16')
					print('Set', channelObj.type, 'p0 to:', channelObj.adc_to_volts_cal['p0'], ', p1 to:', channelObj.adc_to_volts_cal['p1'], ', base_line_adcs to:', channelObj.base_line_adcs, 'and base_line_u_adcs to:', channelObj.base_line_u_adcs)

		if self.settings:
			if not 'voltage_calib_dir' in list(self.settings.__dict__.keys()):
				print('The adc voltage calibration directory was not defined for this run. ', end=' ') ; sys.stdout.flush()
				if vcaldir != '':
					self.voltage_cal_dir = vcaldir
					print('Using parsed directory:', vcaldir)
				else:
					print('Please enter it:')
			else:
				self.voltage_cal_dir = self.settings.voltage_calib_dir if vcaldir == '' else vcaldir
			CheckForDir()
			Get_voltage_calibration()
			# ipdb.set_trace()
			if self.signal_ch:
				Set_Channel(self.signal_ch)
			if self.trigger_ch:
				Set_Channel(self.trigger_ch)
			if self.veto_ch:
				Set_Channel(self.veto_ch)
		else:
			print('The settings file pickle is currupt for this run')
			return


def main():
	parser = OptionParser()
	parser.add_option('-d', '--outdir', dest='outdir', help='path to folder containing the pickles')
	(options, args) = parser.parse_args()
	outdir = str(options.outdir)
	mp = Modify_Pickles_Caen(outdir)
	return mp

if __name__ == '__main__':
	mp = main()
