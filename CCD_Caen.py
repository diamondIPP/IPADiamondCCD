#!/usr/bin/env python
import os
import shutil
import struct
import subprocess as subp
import sys
import time
from ConfigParser import ConfigParser
from optparse import OptionParser

import ROOT as ro
import numpy as np
import cPickle as pickle
import ipdb

from Channel_Caen import Channel_Caen
from Settings_Caen import Settings_Caen
from HV_Control import HV_Control
from Utils import *
# from memory_profiler import profile

trig_rand_time = 10  # for voltage calibration
# trig_rand_time = 0.2  # for system test
wait_time_hv = 7

class CCD_Caen:
	def __init__(self, infile='None', verbose=False, settingsObj=None, iscal=False):
		print 'Starting CCD program ...'
		self.infile = infile
		self.verb = verbose
		self.is_cal_run = iscal
		if self.infile != 'None':
			self.settings = Settings_Caen(self.infile, self.verb, self.is_cal_run)
			self.settings.ReadInputFile()
		elif settingsObj:
			self.settings = settingsObj
		else:
			ExitMessage('No setting file was given, or settings object. Quitting!')
		self.settings.Get_Calibration_Constants()
		self.settings.Get_voltage_calibration()
		self.settings.SetOutputFiles()

		# Create channel objects for signal, trigger and veto
		self.signal_ch = Channel_Caen(self.settings.sigCh, 'signal_ch', self.verb)
		self.signal_ch.Set_Channel(self.settings)
		self.trigger_ch = Channel_Caen(self.settings.trigCh, 'trigger_ch', self.verb)
		self.trigger_ch.Set_Channel(self.settings)
		self.veto_ch = Channel_Caen(self.settings.acCh, 'veto', self.verb) if not self.is_cal_run else None
		if not self.is_cal_run:
			self.veto_ch.Set_Channel(self.settings)

		# declare extra variables that will be used
		self.fs0, self.ft0, self.fv0 = None, None, None
		self.file_time_stamp = None
		self.hv_control = None
		self.stop_run = False
		self.utils = Utils()
		self.RemoveFiles()
		self.t0, self.t1, self.t2 = None, None, None
		self.p, self.pconv = None, None
		self.total_events = 0
		self.written_events_sig, self.written_events_trig, self.written_events_veto = 0, 0, 0
		self.written_events_time = 0
		self.total_events_sig, self.total_events_trig, self.total_events_veto = 0, 0, 0
		self.session_measured_data_sig, self.session_measured_data_trig, self.session_measured_data_veto = 0, 0, 0
		self.total_merged_data_sig, self.total_merged_data_trig, self.total_merged_data_veto = 0, 0, 0
		self.total_merged_data_time = 0
		self.doMerge = False
		self.min_measured_data = 0
		self.min_data_to_write = 0
		self.events_to_write = 0
		self.read_size = 0
		self.sig_written, self.trg_written, self.veto_written = 0, 0, 0
		self.timestamp_written = 0
		self.session_written_events_sig, self.session_written_events_trg, self.session_written_events_veto = 0, 0, 0
		self.fins, self.fint, self.finv = None, None, None
		self.datas, self.datat, self.datav = None, None, None
		self.time_restart = int(np.ceil(self.settings.time_calib + 35))
		self.is_hv_test = False

	def RemoveFiles(self):
		# used, for example, to remove old files that may have stayed due to crashes
		if self.fs0:
			if not self.fs0.closed:
				self.fs0.close()
				del self.fs0
		if self.ft0:
			if not self.ft0.closed:
				self.ft0.close()
				del self.ft0
		if self.fv0:
			if not self.fv0.closed:
				self.fv0.close()
				del self.fv0
		channels = [self.signal_ch.ch, self.trigger_ch.ch, self.veto_ch.ch] if not self.is_cal_run else [self.signal_ch.ch, self.trigger_ch.ch]
		for ch in channels:
			if os.path.isfile('raw_waves{c}.dat'.format(c=ch)):
				os.remove('raw_waves{c}.dat'.format(c=ch))
			if os.path.isfile('waves{c}.dat'.format(c=ch)):
				os.remove('waves{c}.dat'.format(c=ch))
		del channels

	def StartHVControl(self):
		if self.settings.do_hv_control:
			self.hv_control = HV_Control(self.settings)
			print 'Waiting {t} seconds for the HVClient to start... '.format(t=wait_time_hv), ; sys.stdout.flush()
			time.sleep(wait_time_hv)
			print 'Done'
			self.hv_control.CheckVoltage()

	def AdjustBaseLines(self, ntries=5):
		# Read ADCs from files for trigger and veto scintillators
		ntriggers = 10
		t0 = time.time()
		self.CreateEmptyFiles()
		self.CloseFiles()
		self.settings.SetupDigitiser(doBaseLines=True, signal=self.signal_ch, trigger=self.trigger_ch, ac=self.veto_ch)
		self.p = subp.Popen(['{p}/wavedump'.format(p=self.settings.wavedump_path), '{d}/WaveDumpConfig_CCD_BL.txt'.format(d=self.settings.outdir)], bufsize=-1, stdin=subp.PIPE, close_fds=True)
		time.sleep(2)
		if self.p.poll() is None:
			self.GetWaveforms(int(-1 * ntriggers), True, False)
		else:
			print 'Wavedump did not start'
			self.settings.RemoveBinaries()
			self.RemoveFiles()
			if ntries > 0:
				self.AdjustBaseLines(ntries - 1)
			else:
				ExitMessage('There is a problem with Wavedump... exiting')

		if self.total_events != ntriggers:
			print 'Saved', self.total_events, 'which is different to the', ntriggers, 'sent'
			ntriggers = self.total_events

		with open('raw_wave{t}.dat'.format(t=self.trigger_ch.ch), 'rb') as self.ft0:
			triggADCs = np.empty(0, dtype='H')
			for ev in xrange(ntriggers):
				self.ft0.seek(ev * self.settings.struct_len, 0)
				self.datat = self.ft0.read(self.settings.struct_len)
				t = struct.Struct(self.settings.struct_fmt).unpack_from(self.datat)
				triggADCs = np.append(triggADCs, np.array(t, 'H'))
			mean_t = triggADCs.mean()
			std_t = triggADCs.std()

		with open('raw_wave{ac}.dat'.format(ac=self.veto_ch.ch), 'rb') as self.fv0:
			acADCs = np.empty(0, dtype='H')
			for ev in xrange(ntriggers):
				self.fv0.seek(ev * self.settings.struct_len, 0)
				self.datav = self.fv0.read(self.settings.struct_len)
				ac = struct.Struct(self.settings.struct_fmt).unpack_from(self.datav)
				acADCs = np.append(acADCs, np.array(ac, 'H'))
			mean_ac = acADCs.mean()
			std_ac = acADCs.std()

		# clear possible ADCs with non-baseline signals
		for i in xrange(10):
			condition_t = (np.abs(triggADCs - mean_t) < 3 * std_t)
			mean_t = np.extract(condition_t, triggADCs).mean()
			std_t = np.extract(condition_t, triggADCs).std()
			condition_ac = (np.abs(acADCs - mean_ac) < 3 * std_ac)
			mean_ac = np.extract(condition_ac, acADCs).mean()
			std_ac = np.extract(condition_ac, acADCs).std()

		# set channels such that the baselines are near the maximum ADC's leaving space for the scintillator signals. Adjust threshold values
		self.trigger_ch.Correct_Base_Line(mean_adc=mean_t, sigma_adc=std_t, settings=self.settings)
		self.trigger_ch.Correct_Threshold()
		# self.settings.trig_base_line = np.multiply(self.trigger_ch.base_line_u_adcs, self.settings.sigRes, dtype='f8')
		self.settings.trig_base_line = self.trigger_ch.ADC_to_Volts(self.trigger_ch.base_line_adcs)
		self.settings.trig_thr_counts = self.trigger_ch.thr_counts
		self.veto_ch.Correct_Base_Line(mean_adc=mean_ac, sigma_adc=std_ac, settings=self.settings)
		self.veto_ch.Correct_Threshold()
		# self.settings.ac_base_line = np.multiply(self.veto_ch.base_line_u_adcs, self.settings.sigRes, dtype='f8')
		self.settings.ac_base_line = self.veto_ch.ADC_to_Volts(self.veto_ch.base_line_adcs)
		self.settings.ac_thr_counts = self.veto_ch.thr_counts

		del self.ft0, self.datat, t, triggADCs, mean_t, std_t
		del self.fv0, self.datav, ac, acADCs, mean_ac, std_ac
		self.ft0, self.datat= None, None
		self.fv0, self.datav = None, None

	def CreateEmptyFiles(self):
		self.ft0 = open('raw_wave{t}.dat'.format(t=self.trigger_ch.ch), 'wb')
		self.fs0 = open('raw_wave{s}.dat'.format(s=self.signal_ch.ch), 'wb')
		if not self.is_cal_run:
			self.fv0 = open('raw_wave{a}.dat'.format(a=self.veto_ch.ch), 'wb')
		# for timestamp
		self.file_time_stamp = open('raw_time.dat', 'wb')

	def OpenFiles(self, mode='rb'):
		if not self.fs0:
			self.fs0 = open('raw_wave{s}.dat'.format(s=self.signal_ch.ch), mode)
		if not self.ft0:
			self.ft0 = open('raw_wave{t}.dat'.format(t=self.trigger_ch.ch), mode)
		if not self.fv0:
			self.fv0 = open('raw_wave{a}.dat'.format(a=self.veto_ch.ch), mode)

	def CloseFiles(self):
		print 'Closing files'
		if self.ft0:
			self.ft0.close()
			if self.ft0.closed:
				del self.ft0
				self.ft0 = None
		if self.fs0:
			self.fs0.close()
			if self.fs0.closed:
				del self.fs0
				self.fs0 = None
		if self.fv0:
			self.fv0.close()
			if self.fv0.closed:
				del self.fv0
				self.fv0 = None
		if self.file_time_stamp:
			self.file_time_stamp.close()
			if self.file_time_stamp.closed:
				del self.file_time_stamp
				self.file_time_stamp = None

	def GetWaveforms(self, events=1, stdin=False, stdout=False):
		self.t1 = time.time()
		# Negative events is for correcting the baseline of the scintillators and the sigma to set the thresholds
		if events < 0:
			# while self.p.poll() is None:
			time.sleep(1)
			self.p.stdin.write('c')
			self.p.stdin.flush()
			time.sleep(1)
			self.p.stdin.write('s')
			self.p.stdin.flush()
			if self.settings.plot_waveforms:
				# time.sleep(1)
				self.p.stdin.write('P')
				self.p.stdin.flush()
			# time.sleep(1)
			self.p.stdin.write('W')
			self.p.stdin.flush()
			# time.sleep(1)
			for it in xrange(abs(events)):
				time.sleep(0.5)
				self.p.stdin.write('t')
				self.p.stdin.flush()
			self.p.stdin.write('s')
			self.p.stdin.flush()
			time.sleep(1)
			self.p.stdin.write('q')
			self.p.stdin.flush()
			while self.p.poll() is None:
				continue
			if self.settings.do_hv_control: self.stop_run = self.stop_run or self.hv_control.UpdateHVFile()
			self.ConcatenateBinaries();
			self.CloseSubprocess('wave_dump', stdin=stdin, stdout=stdout)
			self.settings.RemoveBinaries()
		else:
			time.sleep(1)
			self.p.stdin.write('c')
			self.p.stdin.flush()
			time.sleep(1)
			self.p.stdin.write('W')
			self.p.stdin.flush()
			if self.settings.plot_waveforms:
				# time.sleep(1)
				self.p.stdin.write('P')
				self.p.stdin.flush()
			# time.sleep(1)
			self.p.stdin.write('s')
			self.p.stdin.flush()
			self.written_events_sig, self.written_events_trig, self.written_events_veto = 0, 0, 0
			# time.sleep(1)
			self.t2 = time.time()
			while self.p.poll() is None:
				if time.time() - self.t1 >= self.settings.time_calib:
					self.p.stdin.write('s')
					self.p.stdin.flush()
					self.settings.RemoveBinaries()
					self.p.stdin.write('c')
					self.p.stdin.flush()
					self.p.stdin.write('q')
					self.p.stdin.flush()
					time.sleep(1)
				elif (self.written_events_sig + self.sig_written >= events) or self.stop_run:
					if self.stop_run: print 'run was stopped'
					self.p.stdin.write('s')
					self.p.stdin.flush()
					self.settings.RemoveBinaries()
					self.p.stdin.write('q')
					self.p.stdin.flush()
					time.sleep(1)
				else:
					if self.settings.random_test and (time.time() - self.t2 > trig_rand_time):
						self.p.stdin.write('t')
						self.p.stdin.flush()
						self.t2 = time.time()
					self.ConcatenateBinaries()
					if self.settings.do_hv_control:
						self.stop_run = self.stop_run or self.hv_control.UpdateHVFile(int(min(self.written_events_sig + self.sig_written, self.settings.num_events)))
					if not self.settings.simultaneous_conversion:
						self.settings.bar.update(int(min(self.written_events_sig + self.sig_written, self.settings.num_events)))
			del self.t1
			self.t1 = None
			self.CloseSubprocess('wave_dump', stdin=stdin, stdout=stdout)
			time.sleep(1)
			del self.t2
			self.t2 = None
			if self.settings.do_hv_control: self.stop_run = self.stop_run or self.hv_control.UpdateHVFile(int(min(self.written_events_sig + self.sig_written, self.settings.num_events)))
		self.total_events_sig = self.CalculateEventsWritten(self.signal_ch.ch)
		self.total_events_trig = self.CalculateEventsWritten(self.trigger_ch.ch)
		self.total_events_veto = self.CalculateEventsWritten(self.veto_ch.ch) if not self.is_cal_run else 0
		if self.total_events_sig == self.total_events_trig:
			self.total_events = self.total_events_sig
		else:
			print 'Written events are of different sizes (signal: {s}, trigger: {t}, veto: {v}). Missmatch!'.format(s=self.total_events_sig, t=self.total_events_trig, v=self.total_events_veto)
			exit()
		del self.total_events_sig, self.total_events_trig, self.total_events_veto
		self.total_events_sig, self.total_events_trig, self.total_events_veto = None, None, None

	def CloseSubprocess(self, pname='wave_dump', stdin=False, stdout=False):
		p = self.p if pname == 'wave_dump' else self.pconv if pname == 'converter' else None
		if not p:
			print 'Something failed! Exiting!'
			exit()
		pid = p.pid
		if stdin:
			p.stdin.close()
		if stdout:
			p.stdout.close()
		if p.wait() is None:
			print 'Could not terminate subprocess... forcing termination'
			p.kill()
			if p.wait() is None:
				print 'Could not kill subprocess... quitting'
				exit()
		try:
			os.kill(pid, 0)
		except OSError:
			pass
		else:
			print 'The subprocess is still running. Killing it with os.kill'
			os.kill(pid, 15)
			try:
				os.kill(pid, 0)
			except OSError:
				pass
			else:
				print 'The process does not die... quitting program'
				exit()
		del p, pid

		if pname == 'wave_dump':
			del self.p
			self.p = None
		elif pname == 'converter':
			del self.pconv
			self.pconv = None

	def ConcatenateBinaries(self):
		self.session_measured_data_sig, self.session_measured_data_trig, self.session_measured_data_veto = 0, 0, 0
		if os.path.isfile('wave{s}.dat'.format(s=self.signal_ch.ch)) and os.path.isfile('wave{t}.dat'.format(t=self.trigger_ch.ch)) and (self.is_cal_run or os.path.isfile('wave{a}.dat'.format(a=self.veto_ch.ch))):
			self.session_measured_data_sig = int(os.path.getsize('wave{s}.dat'.format(s=self.signal_ch.ch)))
			self.session_measured_data_trig = int(os.path.getsize('wave{t}.dat'.format(t=self.trigger_ch.ch)))
			self.session_measured_data_veto = int(os.path.getsize('wave{a}.dat'.format(a=self.veto_ch.ch))) if not self.is_cal_run else 0

		self.total_merged_data_sig = int(os.path.getsize('raw_wave{s}.dat'.format(s=self.signal_ch.ch)))
		self.total_merged_data_trig = int(os.path.getsize('raw_wave{t}.dat'.format(t=self.trigger_ch.ch)))
		self.total_merged_data_veto = int(os.path.getsize('raw_wave{a}.dat'.format(a=self.veto_ch.ch))) if not self.is_cal_run else 0
		self.total_merged_data_time = int(os.path.getsize('raw_time.dat'))
		self.doMerge = (self.session_measured_data_sig + self.sig_written * self.settings.struct_len > self.total_merged_data_sig) and (self.session_measured_data_trig + self.trg_written * self.settings.struct_len > self.total_merged_data_trig)
		if not self.is_cal_run:
			self.doMerge = self.doMerge and (self.session_measured_data_veto + self.veto_written * self.settings.struct_len > self.total_merged_data_veto)
		if self.doMerge:
			# self.OpenFiles(mode='ab')
			self.min_measured_data = min(self.session_measured_data_sig, self.session_measured_data_trig, self.session_measured_data_veto) if not self.is_cal_run else min(self.session_measured_data_sig, self.session_measured_data_trig)
			data_to_write_sig = self.min_measured_data - self.total_merged_data_sig + self.sig_written * self.settings.struct_len
			data_to_write_trg = self.min_measured_data - self.total_merged_data_trig + self.trg_written * self.settings.struct_len
			data_to_write_aco = self.min_measured_data - self.total_merged_data_veto + self.veto_written * self.settings.struct_len if not self.is_cal_run else 0
			self.min_data_to_write = min(data_to_write_sig, data_to_write_trg, data_to_write_aco) if not self.is_cal_run else min(data_to_write_sig, data_to_write_trg)
			del data_to_write_sig, data_to_write_trg, data_to_write_aco
			self.events_to_write = int(np.floor(self.min_data_to_write / float(self.settings.struct_len)))
			self.read_size = self.events_to_write * self.settings.struct_len

			with open('wave{s}.dat'.format(s=self.signal_ch.ch), 'rb') as self.fins:
				self.fins.seek(self.written_events_sig * self.settings.struct_len, 0)
				self.datas = self.fins.read(self.read_size)
			del self.fins
			self.fins = None
			with open('raw_wave{s}.dat'.format(s=self.signal_ch.ch), 'ab') as self.fs0:
				self.fs0.write(self.datas)
				self.fs0.flush()
			del self.fs0, self.datas
			self.fs0, self.datas = None, None

			with open('wave{t}.dat'.format(t=self.trigger_ch.ch), 'rb') as self.fint:
				self.fint.seek(self.written_events_trig * self.settings.struct_len, 0)
				self.datat = self.fint.read(self.read_size)
			del self.fint
			self.fint = None
			with open('raw_wave{t}.dat'.format(t=self.trigger_ch.ch), 'ab') as self.ft0:
				self.ft0.write(self.datat)
				self.ft0.flush()
			del self.ft0, self.datat
			self.ft0, self.datat = None, None

			if not self.is_cal_run:
				with open('wave{a}.dat'.format(a=self.veto_ch.ch), 'rb') as self.finv:
					self.finv.seek(self.written_events_veto * self.settings.struct_len, 0)
					self.datav = self.finv.read(self.read_size)
				del self.finv
				self.finv = None
				with open('raw_wave{a}.dat'.format(a=self.veto_ch.ch), 'ab') as self.fv0:
					self.fv0.write(self.datav)
					self.fv0.flush()
				del self.fv0, self.datav
				self.fv0, self.datav = None, None

			temptime = time.time()
			temptimes = int(temptime)
			temptimens = int(1e9 * (temptime - temptimes))
			datatime = struct.pack(self.settings.time_struct_fmt, temptimes, temptimens)
			with open('raw_time.dat', 'ab') as self.file_time_stamp:
				for ev in xrange(self.events_to_write):
					self.file_time_stamp.write(datatime)
					self.file_time_stamp.flush()
			del self.file_time_stamp
			self.file_time_stamp = None

			self.written_events_sig += int(self.events_to_write)
			self.written_events_trig += int(self.events_to_write)
			self.written_events_veto += int(self.events_to_write) if not self.is_cal_run else 0
			self.written_events_time += int(self.events_to_write)

			del self.events_to_write, self.read_size
			self.events_to_write, self.read_size = None, None

		self.doMerge = False

	def CalculateEventsWritten(self, ch):
		if ch != -1:
			# it is not timestamp file. it is a digitized channel
			return int(round(float(os.path.getsize('raw_wave{c}.dat'.format(c=ch))) / float(self.settings.struct_len)))
		else:
			# it is for timestamp
			return int(round(float(os.path.getsize('raw_time.dat')) / float(self.settings.time_struct_len)))

	# @profile(precision=12)
	def GetData(self, ntries=5):
		self.t0 = time.time()
		self.CreateEmptyFiles()
		self.CloseFiles()
		self.total_events = 0
		print 'Getting {n} events...'.format(n=self.settings.num_events)
		if self.settings.simultaneous_conversion:
			self.CreateRootFile(files_moved=False)
		else:
			self.settings.CreateProgressBar(self.settings.num_events)
			self.settings.bar.start()
		self.settings.SetupDigitiser(doBaseLines=False, signal=self.signal_ch, trigger=self.trigger_ch, ac=self.veto_ch, events_written=self.total_events)
		doBreak = False
		while (self.total_events < self.settings.num_events) and not doBreak:
			print 'Calculating written events'
			self.sig_written = self.CalculateEventsWritten(self.signal_ch.ch)
			self.trg_written = self.CalculateEventsWritten(self.trigger_ch.ch)
			self.veto_written = self.CalculateEventsWritten(self.veto_ch.ch) if not self.is_cal_run else 0
			self.timestamp_written = self.CalculateEventsWritten(-1)
			if not self.stop_run:
				print 'Starting wavedump'
				self.p = subp.Popen(['{p}/wavedump'.format(p=self.settings.wavedump_path), '{d}/WaveDumpConfig_CCD.txt'.format(d=self.settings.outdir)], bufsize=-1, stdin=subp.PIPE, stdout=subp.PIPE, close_fds=True)
				time.sleep(2)
				if self.p.poll() is None:
					self.GetWaveforms(self.settings.num_events, stdin=True, stdout=True)
				else:
					print 'Wavedump did not start'
					if ntries <= 0:
						ExitMessage('There is a problem with Wavedump... exiting')
					else:
						ntries = ntries - 1
			else:
				doBreak = True
				print 'run was stopped. Wavedump will not start (again)'

		self.CloseFiles()
		if not self.settings.simultaneous_conversion:
			print 'Time getting {n} events: {t} seconds'.format(n=self.total_events, t=time.time() - self.t0)
			self.settings.bar.finish()
		else:
			wait_time_if_stopped = 10
			temp_t = time.time()
			while self.pconv.poll() is None:
				time.sleep(2)
				if self.stop_run and time.time() - temp_t > wait_time_if_stopped:
					print 'killing converter because run was stopped'
					self.pconv.kill()
			self.CloseSubprocess('converter', stdin=False, stdout=False)
		return self.total_events

	def CreateRootFile(self, files_moved=False):
		settings_bin_path = os.path.abspath(self.settings.outdir + '/Runs/{f}/{f}.settings'.format(f=self.settings.filename))
		data_bin_path = os.path.abspath(self.settings.outdir + '/Runs/{f}'.format(f=self.settings.filename)) if files_moved else os.getcwd()
		conv_command = ['python', 'Converter_Caen.py', settings_bin_path, data_bin_path]
		if files_moved:
			conv_command.append('0')
		self.pconv = subp.Popen(conv_command, close_fds=True)
		del settings_bin_path

	def CloseHVClient(self):
		if self.settings.do_hv_control:
			self.hv_control.CloseClient()

	def SavePickles(self):
		# save objects of settings, signal_ch, trigger_ch and veto_ch
		with open('{d}/Runs/{f}/{f}.settings'.format(d=self.settings.outdir, f=self.settings.filename), 'wb') as fs:
			pickle.dump(self.settings, fs, pickle.HIGHEST_PROTOCOL)
		with open('{d}/Runs/{f}/{f}.signal_ch'.format(d=self.settings.outdir, f=self.settings.filename), 'wb') as fsig:
			pickle.dump(self.signal_ch, fsig, pickle.HIGHEST_PROTOCOL)
		with open('{d}/Runs/{f}/{f}.trigger_ch'.format(d=self.settings.outdir, f=self.settings.filename), 'wb') as ft:
			pickle.dump(self.trigger_ch, ft, pickle.HIGHEST_PROTOCOL)
		with open('{d}/Runs/{f}/{f}.veto'.format(d=self.settings.outdir, f=self.settings.filename), 'wb') as fv:
			pickle.dump(self.veto_ch, fv, pickle.HIGHEST_PROTOCOL)

	def PrintPlotLimits(self, ti=-5.12e-7, tf=4.606e-6, vmin=-0.7, vmax=0.05):
		print np.double([(tf-ti)/float(self.settings.time_res) +1, ti-self.settings.time_res/2.0,
		                      tf+self.settings.time_res/2.0, (vmax-vmin)/self.settings.sigRes, vmin, vmax])

def main():
	parser = OptionParser()
	parser.add_option('-i', '--infile', dest='infile', default='None', type='string',
	                  help='Input configuration file. e.g. CAENCalibration.ini')
	parser.add_option('-v', '--verbose', dest='verb', default=False, help='Toggles verbose', action='store_true')
	parser.add_option('-a', '--automatic', dest='auto', default=False, help='Toggles automatic conversion and analysis afterwards', action='store_true')
	parser.add_option('-c', '--calibration_run', dest='calrun', default=False, action='store_true', help='Used for calibration runs with pulser')
	parser.add_option('-t', '--time', dest='stabilizationtime', type='int', default=180, help='Time in seconds to wait before taking data due to HV stabilization. Default: 180. For calibration runs, only 10 seconds is used')
	parser.add_option('--hvtest', dest='hvtest', default=False, action='store_true', help='option for HV test')

	(options, args) = parser.parse_args()
	infile = str(options.infile)
	auto = bool(options.auto)
	verb = bool(options.verb)
	iscal = bool(options.calrun)
	time_wait = int(options.stabilizationtime) if not iscal else 10
	ishvtest = bool(options.hvtest)
	ccd = CCD_Caen(infile, verb, iscal=iscal)

	if auto or iscal:
		if not iscal:
			ccd.StartHVControl()
			if not ishvtest:
				ccd.AdjustBaseLines()
		ccd.SavePickles()
		time.sleep(time_wait)
		written_events = ccd.GetData()
		if ccd.stop_run: print 'Run stopped because current is too high'
		ccd.settings.num_events = written_events
		ccd.SavePickles()  # update pickles with the real amount of written events
		ccd.settings.MoveBinaryFiles()
		ccd.settings.RenameDigitiserSettings()
		if not iscal:
			ccd.CloseHVClient()
		if not ccd.settings.simultaneous_conversion:
			ccd.CreateRootFile(files_moved=True)
			while ccd.pconv.poll() is None:
				time.sleep(3)
			ccd.CloseSubprocess('converter', stdin=False, stdout=False)

	print 'Finished :)'
	sys.stdout.write('\a\a\a')
	sys.stdout.flush()
	return ccd

if __name__ == '__main__':
	ccd = main()
