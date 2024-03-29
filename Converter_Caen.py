#!/usr/bin/env python
import visa
import csv
import numpy as np
from struct import unpack
import time, os, sys
from optparse import OptionParser
import progressbar
import ipdb
# from pykeyboard import PyKeyboard
from ConfigParser import ConfigParser
import subprocess as subp
import struct
import ROOT as ro
import pickle, re
import shutil
from copy import deepcopy
from Utils import *


# from DataAcquisition import DataAcquisition

class Converter_Caen:
	def __init__(self, settings_object='', data_path='', simultaneous_data_conv=True):
		self.settings_object = settings_object

		self.settings_full_path = os.path.abspath(settings_object)
		self.output_dir = '/'.join(self.settings_full_path.split('/')[:-1])
		self.raw_dir = self.output_dir if data_path == '' else data_path
		self.filename = self.settings_full_path.split('/')[-1].split('.settings')[0]
		self.settings = pickle.load(open('{d}/{f}.settings'.format(d=self.output_dir, f=self.filename), 'rb'))
		self.signal_ch = pickle.load(open('{d}/{f}.signal_ch'.format(d=self.output_dir, f=self.filename), 'rb'))
		self.trigger_ch = pickle.load(open('{d}/{f}.trigger_ch'.format(d=self.output_dir, f=self.filename), 'rb'))
		self.is_cal_run = self.settings.is_cal_run if 'is_cal_run' in self.settings.__dict__.keys() else False
		self.doVeto = True if not self.is_cal_run else False
		self.veto_ch = pickle.load(open('{d}/{f}.veto'.format(d=self.output_dir, f=self.filename), 'rb')) if self.doVeto else None

		self.settings.simultaneous_conversion = simultaneous_data_conv #  overrides the flag used while taking data, if it is converted offline
		self.control_hv = self.settings.do_hv_control
		self.r_passive = self.settings.r_passive if 'r_passive' in self.settings.__dict__.keys() else 230e6

		self.signal_path = data_path + '/raw_wave{chs}.dat'.format(chs=self.settings.sigCh) if self.settings.simultaneous_conversion else data_path + '/' + self.filename + '_signal.dat'
		self.trigger_path = data_path + '/raw_wave{cht}.dat'.format(cht=self.settings.trigCh) if self.settings.simultaneous_conversion else data_path + '/' + self.filename + '_trigger.dat'
		self.veto_path = data_path + '/raw_wave{cha}.dat'.format(cha=self.settings.acCh) if self.settings.simultaneous_conversion else data_path + '/' + self.filename + '_veto.dat'
		self.time_path = data_path + '/raw_time.dat' if self.settings.simultaneous_conversion else data_path + '/' + self.filename + '_time.dat'
		self.hv_log_files_path = None
		self.current_hv_log_path = None
		if self.control_hv:
			self.hv_log_files_path = data_path + '/{f}/{d}_CH{ch}'.format(f=self.filename, d=self.settings.hv_supply, ch=self.settings.hv_ch) if self.settings.simultaneous_conversion else '{d}/Runs/{f}/HV_{f}/{s}_CH{ch}'.format(d=self.settings.outdir, f=self.filename, s=self.settings.hv_supply, ch=self.settings.hv_ch)

			print 'HV log files are in:', self.hv_log_files_path

		self.points = self.settings.points
		self.num_events = self.settings.num_events
		self.struct_len = self.settings.struct_len
		self.struct_fmt = self.settings.struct_fmt
		self.adc_res = self.settings.sigRes
		self.adc_offset = 0
		self.sig_offset = self.signal_ch.dc_offset_percent
		self.trig_offset = self.trigger_ch.dc_offset_percent
		self.anti_co_offset = self.veto_ch.dc_offset_percent if self.doVeto else 0
		self.time_res = self.settings.time_res
		self.post_trig_percent = self.settings.post_trig_percent
		self.trig_value = self.settings.ADC_to_Volts(self.settings.GetTriggerValueADCs(self.trigger_ch), self.trigger_ch)
		self.veto_value = self.veto_ch.thr_counts if self.doVeto else 0
		self.dig_bits = self.settings.dig_bits
		self.simultaneous_conversion = self.settings.simultaneous_conversion
		self.time_recal = self.settings.time_calib
		if not self.is_cal_run:
			self.polarity = 1 if self.settings.bias >= 0 else -1
		else:
			self.polarity = 1 if self.settings.pulser_amplitude >=0 else -1

		self.hv_file_name = 'hvfile_{f}.dat'.format(f=self.filename)
		self.hv_dict = None
		self.hv_pos = 0
		self.hv_struct_fmt = self.settings.hv_struct_fmt
		self.hv_struct_len = self.settings.hv_struct_len

		self.trigger_search_window = 0.2e-6
		self.veto_window_around_trigg = 50e-9
		self.peak_pos_estimate = 2e-6
		self.peak_pos_window_low = 1.5e-6
		self.peak_pos_window_high = 3e-6

		self.array_points = np.arange(self.points, dtype=np.dtype('int32'))

		self.raw_file = None
		self.raw_tree = None
		self.eventBra = self.voltBra = self.trigBra = self.vetoBra = self.timeBra = self.vetoedBra = self.badShapeBra = self.badPedBra = None
		self.satEventBra = None
		self.hvVoltageBra = self.hvCurrentBra = None
		self.voltageDiaBra = None
		# self.hourBra = self.minuteBra = self.secondBra = None
		self.hourMinSecBra = None
		# self.timeStampBra = None
		self.time_struct_fmt = '@II'
		self.time_struct_len = struct.calcsize(self.time_struct_fmt)
		try:
			self.time_struct_fmt = self.settings.time_struct_fmt
			self.time_struct_len = self.settings.time_struct_len
		except AttributeError:
			self.time_struct_fmt = '@II'
			self.time_struct_len = struct.calcsize(self.time_struct_fmt)

		self.t0 = time.time()

		self.signal_written_events = self.trigger_written_events = self.anti_co_written_events = None
		self.timestamp_written_events = None
		self.fs = self.ft = self.fa = None
		self.ftime = None
		self.wait_for_data = None

		self.datas = self.datat = self.dataa = None
		self.sigADC = self.trigADC = self.vetoADC = None
		self.sigVolts = self.trigVolts = self.vetoVolts = None
		self.trigPos = None
		self.timeVect = None
		self.vetoed_event = None
		self.bad_shape_event = None
		self.bad_pedstal_event = None
		self.sat_event = None
		self.condition_base_line = None
		self.condition_peak_pos = None
		self.hv_voltage_event = None
		self.hv_current_event = None
		self.hour_event = self.minute_event = self.second_event = None
		self.time_break = int(np.ceil(self.time_recal + 30))
		self.file_hv = None
		self.hv_raw_data = None
		self.hv_struct = None
		self.hv_data = {'event': 0, 'seconds': 0, 'nanoseconds': 0, 'voltage': 0, 'current': 0}

		self.datatime = None
		self.tempTime = ro.TTimeStamp()
		self.tempYear, self.tempMonth, self.tempDay = np.zeros(1, 'int32'), np.zeros(1, 'int32'), np.zeros(1, 'int32')
		self.tempHour, self.tempMinute, self.tempSecond = np.zeros(1, 'int32'), np.zeros(1, 'int32'), np.zeros(1, 'int32')
		# self.hour_min_sec_event = None
		# self.currentTime = None

		self.struct_s = self.struct_t = self.struct_ac = None
		self.struct_time = None

		self.bar = None

	def CheckTimeStampRaw(self):
		if os.path.isfile('{wd}/{r}.root'.format(wd=self.output_dir, r=self.filename)):
			tempfile = ro.TFile('{wd}/{r}.root'.format(wd=self.output_dir, r=self.filename), 'READ')
			temptree = tempfile.Get(self.filename)
			if temptree:
				if temptree.FindLeaf('timeHV'):
					if temptree.GetLeaf('timeHV').GetTypeName() != 'TDatime':
						if not os.path.isfile(self.time_path):
							print 'Extracting timestamp from existing root tree.'
							temptimefile = open(self.time_path, 'wb')
							temptimefile.close()
							leng = temptree.Draw('timeHV.AsDouble()', '', 'goff')
							while leng > temptree.GetEstimate():
								temptree.SetEstimate(leng)
								leng = temptree.Draw('timeHV.AsDouble()', '', 'goff')
							timehv = temptree.GetVal(0)
							timehv = np.array([timehv[i] for i in xrange(leng)], dtype='f8')
							print 'Finished extracting timestamp from existing root tree.'
							timehvseconds = timehv.astype('int32')
							timehvnanoseconds = np.multiply(1e9, np.subtract(timehv, timehvseconds, dtype='f8'), dtype='f8').astype('int32')
							print 'Extracted time in seconds and nanoseconds. Creating binary raw file'
							with open(self.time_path, 'ab') as temptimefile:
								for tsec, tnsec in zip(timehvseconds, timehvnanoseconds):
									datatime = struct.pack(self.time_struct_fmt, int(tsec), int(tnsec))
									temptimefile.write(datatime)
									temptimefile.flush()
							print 'Finished creating raw timestamp file'
				tempfile.Close()
			return

	def SetupRootFile(self):
		if self.simultaneous_conversion:
			print 'Start creating root file simultaneously with data taking'
		else:
			print 'Checking if there is enough data'
			self.CheckSettingsAndBinaries()
			print 'Start creating root file'
		self.raw_file = ro.TFile('{wd}/{r}.root'.format(wd=self.output_dir, r=self.filename), 'RECREATE')
		self.raw_tree = ro.TTree(self.filename, self.filename)
		self.raw_tree.SetAutoFlush(100)
		self.raw_tree.SetAutoSave(-10485760)
		if self.doVeto:
			self.vetoBra = np.zeros(self.points, 'f8')
			self.vetoedBra = np.zeros(1, '?')
		self.eventBra = np.zeros(1, 'I')
		self.voltBra = np.zeros(self.points, 'f8')
		self.trigBra = np.zeros(self.points, 'f8')
		self.timeBra = np.zeros(self.points, 'f8')
		self.badShapeBra = np.zeros(1, dtype=np.dtype('int8'))  # signed char
		self.badPedBra = np.zeros(1, '?')
		self.satEventBra = np.zeros(1, '?')
		# self.timeStampBra = ro.TTimeStamp()
		self.hourMinSecBra = ro.TTimeStamp()
		if self.control_hv:
			self.hvVoltageBra = np.zeros(1, 'f4')
			self.hvCurrentBra = np.zeros(1, 'f4')
			self.voltageDiaBra = np.zeros(1, 'f4')
		self.raw_tree.Branch('event', self.eventBra, 'event/i')
		self.raw_tree.Branch('time', self.timeBra, 'time[{s}]/D'.format(s=self.points))
		self.raw_tree.Branch('voltageSignal', self.voltBra, 'voltageSignal[{s}]/D'.format(s=self.points))
		self.raw_tree.Branch('voltageTrigger', self.trigBra, 'voltageTrigger[{s}]/D'.format(s=self.points))
		if self.doVeto:
			self.raw_tree.Branch('voltageVeto', self.vetoBra, 'voltageVeto[{s}]/D'.format(s=self.points))
			self.raw_tree.Branch('vetoedEvent', self.vetoedBra, 'vetoedEvent/O')
		self.raw_tree.Branch('badShape', self.badShapeBra, 'badShape/B')  # signed char
		self.raw_tree.Branch('badPedestal', self.badPedBra, 'badPedestal/O')
		self.raw_tree.Branch('satEvent', self.satEventBra, 'satEvent/O')
		self.raw_tree.Branch('timeHV', self.hourMinSecBra)
		# self.raw_tree.Branch('timeStamp', self.timeStampBra)
		if self.control_hv:
			self.raw_tree.Branch('voltageHV', self.hvVoltageBra, 'voltageHV/F')
			self.raw_tree.Branch('currentHV', self.hvCurrentBra, 'currentHV/F')
			self.raw_tree.Branch('voltageDia', self.voltageDiaBra, 'voltageDia/F')
		# self.raw_tree.Branch('timeHV', self.hourMinSecBra)

	def GetBinariesNumberWrittenEvents(self):
		self.signal_written_events = int(round(os.path.getsize(self.signal_path) / self.struct_len)) if os.path.isfile(self.signal_path) else 0
		self.trigger_written_events = int(round(os.path.getsize(self.trigger_path) / self.struct_len)) if os.path.isfile(self.trigger_path) else 0
		self.anti_co_written_events = int(round(os.path.getsize(self.veto_path) / self.struct_len)) if self.doVeto and os.path.isfile(self.veto_path) else 0
		self.timestamp_written_events = int(round(os.path.getsize(self.time_path) / self.time_struct_len)) if os.path.isfile(self.time_path) else 0

	def OpenRawBinaries(self):
		self.fs = open(self.signal_path, 'rb')
		self.ft = open(self.trigger_path, 'rb')
		self.fa = open(self.veto_path, 'rb') if self.doVeto and os.path.isfile(self.veto_path) else None
		if os.path.isfile(self.time_path):
			self.ftime = open(self.time_path, 'rb')

	def CreateProgressBar(self, maxVal=1):
		widgets = [
			'Processed: ', progressbar.Counter(),
			' out of {mv} '.format(mv=maxVal), progressbar.Percentage(),
			' ', progressbar.Bar(marker='>'),
			' ', progressbar.Timer(),
			' ', progressbar.ETA()
			# ' ', progressbar.AdaptativeETA(),
			#  ' ', progressbar.AdaptativeTransferSpeed()
		]
		self.bar = progressbar.ProgressBar(widgets=widgets, maxval=maxVal)

	def ConvertEvents(self):
		self.bar.start()
		for ev in xrange(self.num_events):
			self.CheckFilesSizes(ev)
			self.WaitForData(ev)
			self.ReadData(ev)
			self.CheckData()
			self.struct_s = struct.Struct(self.struct_fmt).unpack_from(self.datas)
			self.sigADC = np.array(self.struct_s, 'H')
			self.sigVolts = self.ADC_to_Volts('signal')
			self.struct_t = struct.Struct(self.struct_fmt).unpack_from(self.datat)
			self.trigADC = np.array(self.struct_t, 'H')
			self.trigVolts = self.ADC_to_Volts('trigger')
			self.LookForTime0()
			self.timeVect = np.linspace(-self.trigPos * self.time_res, self.time_res * (self.points - 1 - self.trigPos), self.points, dtype='f8')
			if self.doVeto:
				self.struct_ac = struct.Struct(self.struct_fmt).unpack_from(self.dataa)
				self.vetoADC = np.array(self.struct_ac, 'H')
				self.vetoVolts = self.ADC_to_Volts('veto')
				self.vetoed_event = self.IsEventVetoed()
			self.DefineSignalBaseLineAndPeakPosition()
			self.bad_shape_event = self.IsEventBadShape()
			self.bad_pedstal_event = self.IsPedestalBad()
			self.sat_event = self.IsEventSaturated()
			if self.datatime:
				self.struct_time = struct.Struct(self.time_struct_fmt).unpack_from(self.datatime)
				self.hv_data['seconds'] = int(self.struct_time[0])
				self.hv_data['nanoseconds'] = int(self.struct_time[1])
			self.FillBranches(ev)
			if ev == 10:
				self.raw_tree.OptimizeBaskets()
			self.raw_tree.Fill()
			self.bar.update(ev + 1)

	def CheckFilesSizes(self, ev):
		self.wait_for_data = (self.signal_written_events <= ev) or (self.trigger_written_events <= ev)
		if self.doVeto:
			self.wait_for_data = self.wait_for_data or (self.anti_co_written_events <= ev)
		self.wait_for_data = self.wait_for_data or (self.timestamp_written_events <= ev)

	def WaitForData(self, ev):
		t1 = time.time()
		while self.wait_for_data:
			if self.simultaneous_conversion:
				if time.time() - t1 > self.time_break:
					print 'No data has been saved in file for event {ev} in the past {t} seconds... exiting!'.format(ev=ev, t=self.time_break)
					exit(os.EX_NOINPUT)
				if not self.fs.closed:
					self.fs.close()
				if not self.ft.closed:
					self.ft.close()
				if self.doVeto:
					if not self.fa.closed:
						self.fa.close()
				if not self.ftime.closed:
					self.ftime.close()
				self.GetBinariesNumberWrittenEvents()
				self.CheckFilesSizes(ev)
				if not self.wait_for_data:
					self.OpenRawBinaries()
			else:
				print 'The data is corrupted... exiting'
				exit()

	def ReadData(self, ev):
		self.fs.seek(ev * self.struct_len, 0)
		self.datas = self.fs.read(self.struct_len)
		self.ft.seek(ev * self.struct_len, 0)
		self.datat = self.ft.read(self.struct_len)
		if self.doVeto:
			self.fa.seek(ev * self.struct_len, 0)
			self.dataa = self.fa.read(self.struct_len)
		if self.control_hv:
			if self.ftime:
				self.ftime.seek(ev * self.time_struct_len, 0)
				self.datatime = self.ftime.read(self.time_struct_len)
				self.Read_HV_File()

	def Read_HV_File(self):
		self.struct_time = struct.Struct(self.time_struct_fmt).unpack_from(self.datatime)
		self.hv_data['seconds'] = int(self.struct_time[0])
		self.hv_data['nanoseconds'] = int(self.struct_time[1])
		self.FindLogFilePath(self.hv_data['seconds'], self.hv_data['nanoseconds'])
		temp_hv_dic = self.FindLineInLog(self.hv_data['seconds'], self.hv_data['nanoseconds'])
		self.hv_data['voltage'] = temp_hv_dic['voltage']
		self.hv_data['current'] = temp_hv_dic['current']

		# line = [0, 0]
		# temp_line = ''
		# if os.path.isfile(self.hv_file_name):
		# 	hv_elems = int(round(os.path.getsize(self.hv_file_name) / self.hv_struct_len))
		# 	with open(self.hv_file_name, 'rb') as self.file_hv:
		# 		for pi in xrange(self.hv_pos, hv_elems):
		# 			self.file_hv.seek(pi * self.hv_struct_len, 0)
		# 			self.hv_raw_data = self.file_hv.read(self.hv_struct_len)
		# 			self.hv_struct = struct.Struct(self.hv_struct_fmt).unpack_from(self.hv_raw_data)
		# 			if self.hv_struct[0] <= ev:
		# 				self.hv_pos = pi
		# 			else:
		# 				self.hv_pos = pi - 1
		# 			if self.hv_struct[0] >= ev:
		# 				break
		# 		self.file_hv.seek(self.hv_pos * self.hv_struct_len, 0)
		# 		self.hv_raw_data = self.file_hv.read(self.hv_struct_len)
		# 		self.hv_struct = struct.Struct(self.hv_struct_fmt).unpack_from(self.hv_raw_data)
		# 		self.hv_data['event'], self.hv_data['seconds'], self.hv_data['nanoseconds'], self.hv_data['voltage'], self.hv_data['current'] = ev, self.hv_struct[1], self.hv_struct[2], self.hv_struct[3], self.hv_struct[4]

		# if self.simultaneous_conversion:
		# 	if os.path.isfile(self.hv_file_name):
		# 		with open(self.hv_file_name, 'rb') as self.file_hv:
		# 			self.file_hv.seek(self.hv_pos * self.hv_struct_len, 0)
		# 			self.hv_raw_data = self.file_hv.read(self.hv_struct_len)
		# 			self.hv_struct = struct.Struct(self.hv_struct_fmt).unpack_from(self.hv_raw_data)
		# 			if self.hv_struct[0] > ev:
		# 				self.hv_pos
		# 			temp_line = self.file_hv.readline()
		# 			line = temp_line.split() if temp_line != '' else line
		# self.hv_voltage_event, self.hv_current_event = float(line[0]), float(line[1])
		# del line, temp_line, self.file_hv
		# self.file_hv = None

	def FindLogFilePath(self, timesec, timens):
		list_logs = glob.glob('{d}/*.log'.format(d=self.hv_log_files_path))
		if not list_logs:
			return
		list_logs.sort(key=lambda x: os.path.getmtime(x))
		position = -1
		for it, filet in enumerate(list_logs):
			if os.path.getmtime(filet) >= timesec + timens * 1e-9:
				position = it
				break
		self.current_hv_log_path = list_logs[position]

	def FindLineInLog(self, timesec, timens):
		lines = []
		temptime = time.localtime(timesec + 1e-9 * timens)
		# do_before = True if temptime[5] == 0 else False
		# do_after = True if temptime[5] == 59 else False
		# temptime_before = time.localtime(timesec + 1e-9 * timens - 1)
		# temptime_after = time.localtime(timesec + 1e-9 * timens + 1)
		if self.current_hv_log_path:
			with open('{f}'.format(f=self.current_hv_log_path), 'r') as current_log:
			# while True:
			# 	line = current_log.readline().split()
			# 	if not line:
			# 		break
			# 	line = line.split()
				# if len(line) > 2 and IsFloat(line[1]) and IsFloat(line[2]):
				# 	lines.append(line)
			# lines.append(line.split() for if len(line.split()) > 2 and IsFloat(line.split()[1]) and IsFloat(line.split()[2]))]
				lines = current_log.readlines()
			lines = [line.split() for line in lines if re.match('{h:02d}:{m:02d}:{s:02d}'.format(h=temptime[3], m=temptime[4], s=temptime[5]), line) and len(line.split()) >= 3 and IsFloat(line.split()[1]) and IsFloat(line.split()[2])]
			# lines = [line.split() for line in lines if len(line.split()) >= 3 and IsFloat(line.split()[1]) and IsFloat(line.split()[2])]
			# current_log.close()
		# self.tempYear = temptime[0]
		# self.tempMonth = temptime[1]
		# self.tempDay = temptime[2]
		# self.tempTime.Set(1970, 1, 1, 0, 0, timesec, timens, True, 0)
		# self.tempTime.GetDate(False, 0, self.tempYear, self.tempMonth, self.tempDay)
		# self.tempTime.GetTime(False, 0, self.tempHour, self.tempMinute, self.tempSecond)
		tempdic = {'voltage': self.hv_data['voltage'], 'current': self.hv_data['current']}
		if len(lines) > 0:
			if len(lines) == 1:
				pos = 0
			else:
				lines2 = [abs(timesec + 1e-9 * timens - time.mktime([temptime[0], temptime[1], temptime[2], int(line[0].split(':')[0]), int(line[0].split(':')[1]), int(line[0].split(':')[2]), 0, 0, -1])) for line in lines]
				# lines2 = np.array([ro.TTimeStamp(int(self.tempYear), int(self.tempMonth), int(self.tempDay), int(line[0].split(':')[0]), int(line[0].split(':')[1]), int(line[0].split(':')[2]), 0, False, 0).AsDouble() for line in lines], 'f8')
				# pos = abs(lines2 - self.tempTime.AsDouble()).argmin()
				# pos = np.argmin(lines2)
				pos = lines2.index(min(lines2))
			tempdic['voltage'] = float(lines[pos][1]) if float(lines[pos][1]) != 0 else tempdic['voltage']
			tempdic['current'] = float(lines[pos][2]) if abs(float(lines[pos][2])) < 100e-6 else tempdic['current']
		return tempdic

	def CheckData(self):
		if not self.datas or not self.datat or (self.doVeto and (not self.dataa)):
			print 'No event in signal or trigger files... exiting'
			exit(os.EX_DATAERR)

	def LookForTime0(self):
		guess_pos = int(round(self.points * (100.0 - self.post_trig_percent)/100.0))
		condition_trigg = np.array(np.abs(self.array_points - guess_pos) <= int(round(self.trigger_search_window/self.time_res)), dtype='?')
		condition_no_trigg = np.bitwise_not(condition_trigg, dtype='?')
		# condition_no_trigg = np.array(1 - condition_trigg, dtype='?')
		# mean = np.extract(condition_no_trigg, self.trigVolts).mean()
		# sigma = np.extract(condition_no_trigg, self.trigVolts).std()
		temp_trig_volts = np.copy(self.trigVolts)
		np.putmask(temp_trig_volts, condition_no_trigg, -100 * self.trigger_ch.edge)
		volt_peak_pos = temp_trig_volts.argmin() if self.trigger_ch.edge < 0 else temp_trig_volts.argmax()
		condition_trigg = np.bitwise_and(condition_trigg, np.array(self.array_points <= volt_peak_pos))
		np.putmask(temp_trig_volts, np.bitwise_not(condition_trigg), -100 * self.trigger_ch.edge)
		self.trigPos = np.abs(temp_trig_volts - self.trig_value).argmin()
		del guess_pos, condition_trigg, condition_no_trigg, temp_trig_volts, volt_peak_pos

	def IsEventVetoed(self):
		condition_veto_base_line = np.array(np.abs(self.array_points - self.trigPos) > int(round(self.veto_window_around_trigg / float(self.time_res))), dtype='?')
		condition_search = np.bitwise_not(condition_veto_base_line, dtype='?')
		# condition_search = np.array(1 - condition_veto_base_line, dtype='?')
		meanbl = np.extract(condition_veto_base_line, self.vetoADC).mean()
		# sigma = np.extract(condition_veto_base_line, self.vetoADC).std()
		# vetoValNew = 4 * sigma if self.veto_value < 4 * sigma else self.veto_value
		# veto_event = bool((np.extract(condition_search, self.vetoADC) - mean + vetoValNew).min() <= 0)
		veto_event = bool((np.extract(condition_search, self.vetoADC) - meanbl + self.veto_value).astype('float32').min() <= 0)
		del condition_search, condition_veto_base_line, meanbl
		return veto_event

	def DefineSignalBaseLineAndPeakPosition(self):
		self.condition_base_line = np.array(self.array_points <= self.trigPos, dtype='?')
		# values shaper of the signal indicates it should peak at ~2us. The window is set between 1.5us before and 3us after the peak position designed by the shaper
		self.condition_peak_pos = np.array(np.abs(self.array_points - ((self.peak_pos_estimate + (self.peak_pos_window_high - self.peak_pos_window_low) / 2.)/float(self.time_res) + self.trigPos)) <= ((self.peak_pos_window_high + self.peak_pos_window_low) / 2.)/float(self.time_res), dtype='?')

	# def IsEventBadShape(self):
	# 	# mean = np.extract(self.condition_base_line, self.sigADC).mean()
	# 	sigma = np.extract(self.condition_base_line, self.sigADC).std()
	# 	lim_inf = self.condition_peak_pos.argmax()
	# 	lim_sup = self.points - self.condition_peak_pos[::-1].argmax() - 1
	# 	peak_pos = self.sigADC.argmin() if self.polarity == 1 else self.sigADC.argmax()
	# 	if lim_inf < peak_pos < lim_sup:
	# 		# The event has a good shape
	# 		return 0
	# 	else:
	# 		modified_adc = self.sigADC - sigma if self.polarity == 1 else self.sigADC + sigma
	# 		modified_adc[lim_inf] += 2*sigma if self.polarity == 1 else -2*sigma
	# 		modified_adc[lim_sup] += 2*sigma if self.polarity == 1 else -2*sigma
	# 		peak_pos = modified_adc.argmin() if self.polarity == 1 else modified_adc.argmax()
	# 		if lim_inf < peak_pos < lim_sup:
	# 			# Can't tell if the event has a bad shape
	# 			return -1
	# 		else:
	# 			# Event has bad shape
	# 			return 1

	def IsEventBadShape(self):
		# mean = np.extract(self.condition_base_line, self.sigADC).mean()
		sigma = np.extract(self.condition_base_line, self.sigADC).std()
		lim_inf = self.condition_peak_pos.argmax()
		lim_sup = self.points - self.condition_peak_pos[::-1].argmax() - 1
		peak_pos = RoundInt(self.peak_pos_estimate / self.time_res + self.trigPos)
		sigInf, sigPeak, sigSup = self.sigADC[lim_inf] * (-self.polarity), self.sigADC[peak_pos] * (-self.polarity), self.sigADC[lim_sup] * (-self.polarity)
		# if lim_inf < peak_pos < lim_sup:
		if sigPeak > sigInf and sigPeak > sigSup:
			# The event has a good shape
			return 0
		else:
			sigInf -= 2 * sigma
			sigSup -= 2 * sigma
			sigPeak += 2 * sigma
			if sigPeak > sigInf and sigPeak > sigSup:
				# Can't tell if the event has a bad shape
				return -1
			else:
				# Event has bad shape
				return 1

	def IsPedestalBad(self):
		sigma = np.extract(self.condition_base_line, self.sigADC).std()
		self.adc_res = self.signal_ch.adc_to_volts_cal['p1']
		sigma_volts = sigma * self.adc_res
		diff_volts = abs(np.mean(self.sigADC[0:250]) - np.mean(self.sigADC[(self.trigPos - 250):self.trigPos])) * self.adc_res
		# diff_volts = abs(int(self.sigADC[0]) - int(self.sigADC[self.trigPos])) * self.adc_res
		if sigma_volts >= 10e-3 or diff_volts >= 15e-3:  # if a signal is not flat enough due to previous unstable states
			return True
		condition_base_line_ini = np.bitwise_and(np.less_equal(self.array_points, self.trigPos), np.greater(self.array_points, self.trigPos - 500))
		condition_base_line_end = np.bitwise_and(np.greater_equal(self.array_points, self.trigPos + 4000), np.less(self.array_points, self.trigPos + 4500))
		meanbl_ini = np.extract(condition_base_line_ini, self.sigADC).mean()
		meanbl_end = np.extract(condition_base_line_end, self.sigADC).mean()
		if abs(meanbl_end - meanbl_ini) * self.adc_res > 50e-3:  # if a signal does not return to base line due to multiple effects
			return True
		# when the signal pedestal is well behaved and the signal returns to baseline
		return False

	def IsEventSaturated(self):
		if RoundInt(2 ** 14 -1) in np.extract(np.bitwise_or(self.condition_base_line, self.condition_peak_pos), self.sigADC):
			return True
		if 0 in np.extract(np.bitwise_or(self.condition_base_line, self.condition_peak_pos), self.sigADC):
			return True
		return False

	def FillBranches(self, ev):
		self.eventBra.fill(ev)
		np.putmask(self.timeBra, np.ones(self.points, '?'), self.timeVect)
		np.putmask(self.voltBra, np.ones(self.points, '?'), self.sigVolts)
		np.putmask(self.trigBra, np.ones(self.points, '?'), self.trigVolts)
		if self.doVeto:
			np.putmask(self.vetoBra, np.ones(self.points, '?'), self.vetoVolts)
			self.vetoedBra.fill(self.vetoed_event)
		self.badShapeBra.fill(self.bad_shape_event)
		self.badPedBra.fill(self.bad_pedstal_event)
		self.satEventBra.fill(self.sat_event)
		# tempTime = time.time()
		# tempSec = int(tempTime)
		# tempNanoSec = int(1e9 * abs(tempTime - tempSec))
		# todo
		# self.timeStampBra.Set(1970, 1, 1, 0, 0, tempSec, tempNanoSec, True, 0)
		self.hourMinSecBra.Set(1970, 1, 1, 0, 0, self.hv_data['seconds'], self.hv_data['nanoseconds'], True, 0)
		if self.control_hv:
			self.hvVoltageBra.fill(self.hv_data['voltage'])
			self.hvCurrentBra.fill(self.hv_data['current'])
			self.voltageDiaBra.fill(self.CalculateDiamondVoltage())

	def CalculateDiamondVoltage(self):
		return self.hv_data['voltage'] - self.hv_data['current'] * self.r_passive

	def CloseAll(self):
		self.bar.finish()
		self.raw_file.Write()
		self.raw_file.Close()
		self.fs.close()
		del self.fs
		self.ft.close()
		del self.ft
		if self.doVeto:
			self.fa.close()
			del self.fa
		self.t0 = time.time() - self.t0
		print 'Time creating root tree:', self.t0, 'seconds'
		exit()

	def IsPedestalNotFlat(self, signalADC, points, trigPos, time_res):
		array_points = np.arange(points, dtype=np.dtype('int32'))
		condition_base_line = np.array(array_points - trigPos <= 0, dtype='?')

	def ADC_to_Volts(self, sig_type):
		def ChannelAdcToVolts(adcs, channel):
			if 'adc_to_volts_cal' in channel.__dict__.keys():
				return np.add(channel.adc_to_volts_cal['p0'], np.multiply(adcs, channel.adc_to_volts_cal['p1'], dtype='f8'), dtype='f8')
			else:
				ExitMessage('The channel object does not have "adc_to_volts_cal". Run Modify_Settings_Caen.py first.', os.EX_USAGE)

		adcs, offset = 0, 0
		channel = None
		if sig_type == 'signal':
			adcs = self.sigADC
			channel = self.signal_ch
			# offset = self.sig_offset
			# self.adc_offset = self.signal_ch.adc_to_volts_cal['p0']
			# self.adc_res = self.signal_ch.adc_to_volts_cal['p1']
		elif sig_type == 'trigger':
			adcs = self.trigADC
			channel = self.trigger_ch
			# offset = self.trig_offset
			# self.adc_offset = self.trigger_ch.adc_to_volts_cal['p0']
			# self.adc_res = self.trigger_ch.adc_to_volts_cal['p1']
		elif sig_type == 'veto':
			adcs = self.vetoADC
			channel = self.veto_ch
			# offset = self.anti_co_offset
			# self.adc_offset = self.veto_ch.adc_to_volts_cal['p0']
			# self.adc_res = self.veto_ch.adc_to_volts_cal['p1']
		else:
			print 'Wrong type. Exiting'
			exit()
		# result = np.add(self.adc_offset, np.multiply(self.adc_res, np.add(adcs, np.multiply(2 ** self.dig_bits - 1.0, offset / 100.0 - 0.5, dtype='f8'), dtype='f8'), dtype='f8'), dtype='f8')
		# result = channel.ADC_to_Volts(adcs)
		result = ChannelAdcToVolts(adcs, channel)
		return result

	def CheckSettingsAndBinaries(self):
		self.GetBinariesNumberWrittenEvents()
		if not self.simultaneous_conversion and (self.num_events > 10 and (self.signal_written_events < 10 or self.trigger_written_events < 10) or self.signal_written_events + self.trigger_written_events == 0):
			if os.path.isfile('{wd}/{r}.root'.format(wd=self.output_dir, r=self.filename)):
				os.remove('{wd}/{r}.root'.format(wd=self.output_dir, r=self.filename))
			ExitMessage('It was a flawed run. There are not enough events for analyisis. Exiting', os.EX_NOINPUT)
			exit()

if __name__ == '__main__':
	# first argument is the path to the settings pickle file
	# second argument is the path of the directory that contains the raw data.
	# By default, it assumes simultaneous data conversion. If the conversion is done offline (aka. not simultaneous), then the 3rd parameter has to be given and should be '0'
	if len(sys.argv) < 2:
		print 'Usage is: Converter_Caen.py <settings_pickle_path> <dir_with_raw_data> 0 for offline conversion)'
		exit()
	settings_object = str(sys.argv[1])  # settings pickle path
	if settings_object in ['-h', '--help']:
		print 'Usage is: Converter_Caen.py <settings_pickle_path> <dir_with_raw_data> 0 for offline conversion)'
		exit()
	print 'settings object', settings_object
	if len(sys.argv) > 2:
		data_path = str(sys.argv[2])  # path where the binary data in adcs is. It is a directory path containing the raw files.
		print 'data_path', data_path
	else:
		data_path = ''
		print 'data_path empty ""'
	is_simultaneous_data_conv = True
	if len(sys.argv) > 3:
		if IsInt(str(sys.argv[3])):
			is_simultaneous_data_conv = bool(int(str(sys.argv[3])))
			print 'simultaneous is now', is_simultaneous_data_conv
	converter = Converter_Caen(settings_object=settings_object, data_path=data_path, simultaneous_data_conv=is_simultaneous_data_conv)

	converter.CheckTimeStampRaw()
	converter.SetupRootFile()
	converter.GetBinariesNumberWrittenEvents()
	converter.OpenRawBinaries()
	converter.CreateProgressBar(converter.num_events)
	converter.ConvertEvents()
	converter.CloseAll()




