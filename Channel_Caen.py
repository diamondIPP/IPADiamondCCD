#!/usr/bin/env python
import visa
import csv
import numpy as np
from struct import unpack
import time, os, sys
from optparse import OptionParser
import progressbar
import ipdb
from pykeyboard import PyKeyboard
from configparser import ConfigParser
import subprocess as subp
import struct
import ROOT as ro
import shutil
from Utils import *
# from DataAcquisition import DataAcquisition


class Channel_Caen:
	def __init__(self, ch=3, ch_type='signal', verb=False):
		self.ch = ch
		self.type = ch_type
		self.verb = verb

		self.base_line_u_adcs = 0
		self.base_line_adcs = 0
		self.sigma_adcs = 10
		self.dc_offset_percent = 0
		self.thr_counts = 0
		self.edge = -1
		self.name = str(ch_type)

		self.adc_to_volts_cal = {'p0': -1, 'p1': np.divide(2, np.subtract(np.power(2, 14, dtype='float64'), 1.0, dtype='float64'), dtype='float64')}

	def Set_Channel(self, settings):
		self.adc_to_volts_cal['p1'] = settings.sigRes
		if self.type == 'signal_ch':
			self.dc_offset_percent = settings.sig_dc_offset_percent if 'sig_dc_offset_percent' in list(settings.__dict__.keys()) else 0
			self.adc_to_volts_cal['p0'] = np.subtract(np.divide(self.dc_offset_percent, 50, dtype='f8'), 1, dtype='f8')
			self.edge = -int(settings.bias/abs(settings.bias)) if settings.bias != 0 else -1
			# Channels 3, 6 and 7 were calibrated with dc voltage and multimeter. The calibration files are on 20180419ch{X}/Runs
			if 'adc_volts_pickle' in list(settings.__dict__.keys()):
				if settings.adc_volts_pickle:
					self.adc_to_volts_cal['p0'] = settings.adc_volts_pickle['fit_p0']
					self.adc_to_volts_cal['p1'] = settings.adc_volts_pickle['fit_p1']
			else:
				if self.ch == 3:
					# With negative bias, signals are positive. With positive bias, signals are negative
					self.adc_to_volts_cal['p0'] = -0.07482169290039371 if settings.bias < 0 else -2.056009761213107
					# self.adc_to_volts_cal['p0'] = 0.035089408942074955 if settings.bias < 0 else -0.02328757136517118
					self.adc_to_volts_cal['p1'] = 0.00013097264906803782 if settings.bias < 0 else 0.00013061281362144022
					# self.adc_to_volts_cal['p1'] = 0.00013089621340339722 if settings.bias < 0 else 0.00013076091653412987
				elif self.ch == 6:
					self.adc_to_volts_cal['p0'] = 0.04183464530415922 if settings.bias < 0 else -0.04426166525468815
					self.adc_to_volts_cal['p1'] = 0.0001294935440001848 if settings.bias < 0 else 0.00012934155437522722
				elif self.ch == 7:
					self.adc_to_volts_cal['p0'] = 0.025512742406801153 if settings.bias < 0 else -0.024895654896994378
					self.adc_to_volts_cal['p1'] = 0.00012875546036321804 if settings.bias < 0 else 0.00013155381396351944
		elif self.type == 'trigger_ch':
			self.dc_offset_percent = -25 if not settings.is_cal_run else 0
			self.adc_to_volts_cal['p0'] = np.subtract(np.divide(self.dc_offset_percent, 50, dtype='f8'), 1, dtype='f8')
			self.base_line_u_adcs = self.Calculate_Universal_ADCs(settings.trig_base_line, self.adc_to_volts_cal['p1']) if not settings.is_cal_run else self.Calculate_Universal_ADCs(settings.trig_cal_base_line, self.adc_to_volts_cal['p1'])
			self.base_line_adcs = self.Volts_to_ADC(settings.trig_base_line) if not settings.is_cal_run else self.Volts_to_ADC(settings.trig_cal_base_line)
			self.thr_counts = settings.trig_thr_counts if not settings.is_cal_run else RoundInt(settings.trig_cal_thr_mv  / 1000. / self.adc_to_volts_cal['p1'])
			self.edge = -1 if not settings.is_cal_run else settings.trig_cal_polarity
		elif self.type == 'veto':
			self.dc_offset_percent = -25
			self.adc_to_volts_cal['p0'] = np.subtract(np.divide(self.dc_offset_percent, 50, dtype='f8'), 1, dtype='f8')
			self.base_line_u_adcs = self.Calculate_Universal_ADCs(settings.ac_base_line, self.adc_to_volts_cal['p1'])
			self.base_line_adcs = self.Volts_to_ADC(settings.ac_base_line)
			self.thr_counts = settings.ac_thr_counts
			self.edge = -1

	def Calculate_Universal_ADCs(self, value_volts, sig_res):
		return np.divide(value_volts, sig_res, dtype='f8')

	def Correct_Base_Line(self, mean_adc, sigma_adc, settings):
		self.sigma_adcs = sigma_adc
		self.base_line_u_adcs = self.Calculate_Universal_ADCs(self.ADC_to_Volts(mean_adc), self.adc_to_volts_cal['p1'])
		self.base_line_adcs = RoundInt(mean_adc, 'uint16')

	def Correct_Threshold(self):
		if self.type == 'trigger_ch':
			self.thr_counts = int(round(max(10 * self.sigma_adcs, self.thr_counts)))
		elif self.type == 'veto':
			self.thr_counts = int(round(max(4 * self.sigma_adcs, self.thr_counts)))

	def ADC_to_Volts(self, adcs):
		return np.add(self.adc_to_volts_cal['p0'], np.multiply(adcs, self.adc_to_volts_cal['p1'], dtype='f8'), dtype='f8')
		# return np.subtract(np.add(np.multiply(adcs, self.adc_to_volts_cal['p1'], dtype='f8'), np.divide(self.dc_offset_percent, 50.0, dtype='f8')), 1, dtype='f8')

	def Volts_to_ADC(self, volts):
		return RoundInt(np.divide(np.subtract(volts, self.adc_to_volts_cal['p0'], dtype='f8'), self.adc_to_volts_cal['p1'], dtype='f8'), 'uint16')
		# return RoundInt(np.multiply(np.divide(1, self.adc_to_volts_cal['p1'], dtype='f8'), np.add(np.subtract(volts, np.divide(self.dc_offset_percent, 50.0, dtype='f8'), dtype='f8'), 1, dtype='f8'), dtype='f8'), 'uint16')

if __name__ == '__main__':
	print('bla')