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
from copy import deepcopy
import glob
from Utils import *


# from DataAcquisition import DataAcquisition

class HV_Control:
	def __init__(self, settings):
		self.settings = settings
		self.hv_supply, self.ch, self.bias, self.current_limit, self.hot_start, self.address = settings.hv_supply, settings.hv_ch, settings.bias, settings.current_limit, settings.hot_start, settings.hv_address
		self.filename, self.dut = settings.filename, settings.dut
		self.Pics_folder_path = settings.pics_folder_path
		self.doControlHV = False if self.hv_supply == '' else True
		self.logs_dir = '{f}/{d}_CH{ch}'.format(f=self.filename, d=self.hv_supply, ch=self.ch)
		self.log_file = None
		self.last_line = {'event': 0, 'seconds': 0, 'nanoseconds': 0, 'voltage': 0, 'current': 0}
		self.ramp = settings.hv_ramp
		self.supply_number = 0
		self.time_update = 2.0
		self.r_passive = self.settings.r_passive if 'r_passive' in list(self.settings.__dict__.keys()) else 230e6
		self.abort_percentage_drop = 0.1
		self.stop_run = False
		self.out_file = None
		self.hv_struct = self.settings.hv_struct_fmt
		self.hv_struct_len = self.settings.hv_struct_len
		self.hv_struct_pack = None
		self.out_file_name = 'hvfile_{f}.dat'.format(f=self.settings.filename)
		self.time0 = None
		if self.Pics_folder_path == '':
			print('Cannot control voltage because Pics folder (Micha) was not found XD')
			self.doControlHV = False
		self.process = None
		if self.doControlHV:
			self.LinkPicsFolder()
			self.CreateConfigFiles()
			if self.hot_start:
				self.process = subp.Popen(['HVClient.py', '-H'], bufsize=-1, stdin=subp.PIPE, stdout=subp.PIPE, close_fds=True)
			else:
				self.process = subp.Popen(['HVClient.py'], bufsize=-1, stdin=subp.PIPE, stdout=subp.PIPE, close_fds=True)
				# print 'Only hot start has been implemented :P'
			self.time0 = time.time()
			time.sleep(3)
			self.process.stdin.write('yes\n')
			self.process.stdin.flush()

	def LinkPicsFolder(self):
		if not os.path.isdir('Pics'):
			os.symlink(self.Pics_folder_path, 'Pics')

	def CreateConfigFiles(self):
		if not os.path.isdir('config'):
			os.mkdir('config')
		self.CreateHVClientConfig()

	def CreateHVClientConfig(self):
		num_supplies = 9
		supplies_num = list(range(1, num_supplies + 1))
		supplies_ids = {1: 'Keithley1', 2: 'Keithley2410', 3: 'Keithley237', 4: 'Keithley6517', 5: '', 6: 'Keithley2657A', 7: 'ISEG-NHS-6220x', 8: 'ISEG-NHS-6220n', 9: 'Keithley6517B'}
		conf_file = open('config/hv_{f}.cfg'.format(f=self.filename), 'w')

		conf_file.write('[Main]\n')
		conf_file.write('devices = [0]\n'.format(n=num_supplies))
		conf_file.write('testbeam_name: {f}\n'.format(f=self.filename))

		conf_file.write('\n[Names]\n')
		conf_file.write('HV{s}: {d}\n'.format(s=self.supply_number, d=self.dut))

		if self.hv_supply == 'ISEG-NHS-6220x' or self.hv_supply == 'ISEG-NHS-6220n':
			conf_file.write('\n[HV{s}]\n'.format(s=self.supply_number))
			conf_file.write('name: {n}\n'.format(n=self.hv_supply))
			conf_file.write('model: ' + '-'.join(self.hv_supply.split('-')[1:]) + '\n')
			conf_file.write('module_name: ISEG\n')
			conf_file.write('nChannels: 6\n')
			conf_file.write('active_channels: [{ch}]\n'.format(ch=self.ch))
			conf_file.write('address: {a}\n'.format(a=self.address))
			conf_file.write('# in V/s\n')
			conf_file.write('ramp: {r}\n'.format(r=self.ramp))
			conf_file.write('config_file: iseg.cfg\n')
		elif self.hv_supply == 'Keithley2410':
			conf_file.write('\n[HV{s}]\n'.format(s=self.supply_number))
			conf_file.write('name: {n}\n'.format(n=self.hv_supply))
			conf_file.write('model: 2410\n')
			conf_file.write('address: {a}\n'.format(a=self.address))
			conf_file.write('compliance: {c} nA\n'.format(c=self.current_limit*1e9))
			conf_file.write('ramp: {r}\n'.format(r=self.ramp))
			conf_file.write('max_step: 10\n')
			conf_file.write('bias: 0\n')
			conf_file.write('min_bias: -1000\n')
			conf_file.write('max_bias: 1000\n')
			conf_file.write('baudrate: 57600\n')
			conf_file.write('output: front\n')
			# TODO: write the other cases for the other supplies
		elif self.hv_supply == 'Keithley6517B':
			conf_file.write('\n[HV{s}]\n'.format(s=self.supply_number))
			conf_file.write('name: {n}\n'.format(n=self.hv_supply))
			conf_file.write('model: 6517\n')
			conf_file.write('address: {a}\n'.format(a=self.address))
			conf_file.write('compliance: {c} nA\n'.format(c=self.current_limit*1e9))
			conf_file.write('ramp: {r}\n'.format(r=self.ramp))
			conf_file.write('max_step: 10\n')
			conf_file.write('bias: 0\n')
			conf_file.write('min_bias: -1000\n')
			conf_file.write('max_bias: 1000\n')
			conf_file.write('baudrate: 57600\n')
			conf_file.write('output: rear\n')
			# TODO: write the other cases for the other supplies

		conf_file.close()
		self.UnlinkConfigFile('keithley.cfg')
		os.symlink('hv_{f}.cfg'.format(f=self.filename), 'config/keithley.cfg')

		if self.hv_supply == 'ISEG-NHS-6220x' or self.hv_supply == 'ISEG-NHS-6220n':
			self.CreateIsegConfigFile()

	def CreateIsegConfigFile(self):
		iseg_chs = 6
		conf_file = open('config/iseg_{f}.cfg'.format(f=self.filename), 'w')
		conf_file.write('[Names]\n')
		for ch in range(iseg_chs):
			if ch == self.ch:
				conf_file.write('CH{ch}: {d}\n'.format(ch=ch, d=self.dut))
			else:
				conf_file.write('CH{ch}: None\n'.format(ch=ch))
		for ch in range(iseg_chs):
			conf_file.write('\n[CH{ch}]\n'.format(ch=ch))
			conf_file.write('name: CH{ch}\n'.format(ch=ch))
			compliance = self.current_limit if ch == self.ch else 250e-9
			conf_file.write('compliance: {c}\n'.format(c=compliance))
			conf_file.write('measure_range: 10e-6\n')
			conf_file.write('bias: 0\n')
			min_bias = -1 if self.bias >= 0 else min(self.bias * 2, -1110)
			max_bias = max(self.bias * 2, 1110) if self.bias >= 0 else 1
			if ch == self.ch:
				conf_file.write('min_bias: {m}\n'.format(m=min_bias))
				conf_file.write('max_bias: {m}\n'.format(m=max_bias))
			else:
				conf_file.write('min_bias: -1\n')
				conf_file.write('max_bias: 1\n')
		conf_file.close()
		self.UnlinkConfigFile('iseg.cfg')
		os.symlink('iseg_{f}.cfg'.format(f=self.filename), 'config/iseg.cfg')

	def UnlinkConfigFile(self, name):
		if os.path.isfile('config/{n}'.format(n=name)) or os.path.islink('config/{n}'.format(n=name)):
			if os.path.islink('config/{n}'.format(n=name)):
				os.unlink('config/{n}'.format(n=name))
			else:
				os.remove('config/{n}'.format(n=name))

	def CheckVoltage(self):
		max_tries = 2
		self.ReadLastLine()
		delta_voltage = abs(self.last_line['voltage'] - self.bias)
		do_ramp = False
		if delta_voltage >= 1:
			do_ramp = True
			print('Ramping voltage... ', end=' ') ; sys.stdout.flush()
		while delta_voltage >= 1 and max_tries != 0:
			self.CorrectBias(delta_voltage)
			self.ReadLastLine()
			delta_voltage = abs(self.last_line['voltage'] - self.bias)
			max_tries -= 1
			if max_tries == 0:
				print('\nCould not set the desired voltage. Taking data with {v}V\n'.format(v=self.last_line['voltage']))
		if do_ramp: print('Done')

	def GetLastLogFilePath(self):
		list_logs = glob.glob('{d}/*.log'.format(d=self.logs_dir))
		if not list_logs:
			return
		self.log_file = max(list_logs, key=os.path.getmtime)
		del list_logs

	def FindLogFilePath(self, timesec, timens):
		list_logs = glob.glob('{d}/*.log'.format(d=self.logs_dir))
		if not list_logs:
			return
		list_logs.sort(key=lambda x: os.path.getmtime(x))
		position = 0
		for it, filet in enumerate(list_logs):
			if os.path.getmtime(filet) >= timesec + timens * 1e-9:
				position = it
				break
		position = -1 if position == 0 else position
		self.log_file = list_logs[position]

	def FindLineInLog(self, timesec, timens):
		lines = []
		if self.log_file:
			current_log = open('{f}'.format(f=self.log_file), 'r')
			lines = current_log.readlines()
			lines = [line.split() for line in lines if len(line.split()) >= 3 and IsFloat(line.split()[1]) and IsFloat(line.split()[2])]
			current_log.close()
		tempTime = ro.TTimeStamp()
		tempTime.Set(1970, 1, 1, 0, 0, timesec, timens, True, 0)
		tempYear, tempMonth, tempDay = np.zeros(1, 'int32'), np.zeros(1, 'int32'), np.zeros(1, 'int32')
		tempHour, tempMinute, tempSecond = np.zeros(1, 'int32'), np.zeros(1, 'int32'), np.zeros(1, 'int32')
		tempTime.GetDate(False, 0, tempYear, tempMonth, tempDay)
		tempTime.GetTime(False, 0, tempHour, tempMinute, tempSecond)
		if len(lines) > 0:
			lines2 = np.array([ro.TTimeStamp(int(tempYear), int(tempMonth), int(tempDay), int(line[0].split(':')[0]), int(line[0].split(':')[1]), int(line[0].split(':')[2]), 0, False, 0).AsDouble() for line in lines], 'f8')
			pos = abs(lines2 - tempTime.AsDouble()).argmin()
			tempdic = {'time_s': lines2[pos], 'voltage': float(lines[pos][1]), 'current': float(lines[pos][2])}
			return tempdic
		return {'time_s': timesec + 1e-9 * timens, 'voltage': 0, 'current': 0}

	def ReadLastLine(self, nEvent=0):
		self.GetLastLogFilePath()
		if self.log_file:
			current_log = open('{f}'.format(d=self.logs_dir, f=self.log_file), 'r')
			lines = current_log.readlines()
			current_log.close()
		else:
			lines = None
		temp_time = time.time()
		self.last_line['event'] = nEvent
		self.last_line['seconds'] = int(temp_time)
		self.last_line['nanoseconds'] = int(1e9 * abs(temp_time - self.last_line['seconds']))
		if not lines:
			return
		if len(lines) >= 1:
			temp_line = lines[-1].split()
			if len(temp_line) >= 3:
				if IsFloat(temp_line[1]) and IsFloat(temp_line[2]):
					self.last_line['voltage'] = float(temp_line[1])
					self.last_line['current'] = float(temp_line[2]) if abs(self.last_line['current']) < 100e-6 else 0
					if self.last_line['voltage'] != 0:
						if abs(self.last_line['current'] * self.r_passive) > abs(self.abort_percentage_drop * self.last_line['voltage']):
							print('due to high current, the voltage in the diamond is below 90% of the intended value. Sending termination signal to the run')
							self.stop_run = True
		return

	def CorrectBias(self, delta_volts):
		if self.hv_supply.startswith('ISE'):
			self.process.stdin.write('ON HV{s} CH{c}\n'.format(s=self.supply_number, c=self.ch))
			self.process.stdin.write('BIAS HV{s} CH{c} {v}\n'.format(s=self.supply_number, c=self.ch, v=self.bias))
		else:
			self.process.stdin.write('ON HV{s}\n'.format(s=self.supply_number))
			self.process.stdin.write('BIAS HV{s} {v}\n'.format(s=self.supply_number, v=self.bias))
		self.process.stdin.flush()
		wait_time = delta_volts/float(self.ramp) + 5
		time.sleep(wait_time)

	def UpdateHVFile(self, nEvent=0):
		self.stop_run = False
		if time.time() - self.time0 >= self.time_update:
			self.ReadLastLine(nEvent)
			self.time0 = time.time()
			self.WriteHVFile()
		return self.stop_run

	def WriteHVFile(self):
		temp_array = [self.last_line['event'], self.last_line['seconds'], self.last_line['nanoseconds'], self.last_line['voltage'], self.last_line['current']]
		self.hv_struct_pack = struct.pack(self.hv_struct, *temp_array)
		with open(self.out_file_name, 'ab') as self.out_file:
			self.out_file.write(self.hv_struct_pack)
		del self.out_file, temp_array
		self.out_file = None

	def SearchForData(self, time_sec, time_ns):
		# todo
		pass

	def CloseClient(self):
		if self.process:
			self.process.stdin.write('exit\n')
			self.process.stdin.flush()
			time.sleep(1)
			self.MoveLogsAndConfig()
			# os.remove(self.out_file_name)

	def MoveLogsAndConfig(self):
		path_dir = '{d}/Runs/{f}/HV_{f}'.format(d=self.settings.outdir, f=self.filename)
		if os.path.isdir(path_dir) or os.path.isfile(path_dir):
			print('HV directory already exists. Recreating')
			if os.path.isdir(path_dir):
				shutil.rmtree(path_dir)
			else:
				os.remove(path_dir)
		shutil.move(self.filename, path_dir)
		'Moved hv logs to output folder'
		if os.path.isfile('config/hv_{f}.cfg'.format(f=self.filename)):
			shutil.move('config/hv_{f}.cfg'.format(f=self.filename), path_dir)
			print('Moved hv config file to output folder')
			if self.hv_supply.lower().startswith('iseg'):
				shutil.move('config/iseg_{f}.cfg'.format(f=self.filename), path_dir)
			shutil.move(self.out_file_name, path_dir)
			print('Moved hv data structure to output folder')
		del path_dir
		if os.path.islink('config/keithley.cfg'):
			os.unlink('config/keithley.cfg')
			print('Unlinked config/keithley.cfg')

if __name__ == '__main__':
	print('blaaaa')


