#!/usr/bin/env python
import os, glob
import shutil
import struct
import subprocess as subp
import sys
import time
from ConfigParser import ConfigParser
from optparse import OptionParser

from collections import OrderedDict

import ROOT as ro
import numpy as np
import cPickle as pickle

from Channel_Caen import Channel_Caen
from Settings_Caen import Settings_Caen
from HV_Control import HV_Control
from Utils import *
from Langaus import LanGaus
from Crystalball import Crystalball
from GausTail import GausTail
from VcalToElectrons import VcalToElectrons
# from memory_profiler import profile

trig_rand_time = 0.2
wait_time_hv = 7
BRANCHES1DTOTAL = ['event', 'vetoedEvent', 'badShape', 'badPedestal', 'satEvent', 'voltageDia', 'voltageHV', 'currentHV',
				   'timeHV', 'peakPosition', 'peakPositionCF', 'pedestal', 'pedestalCF', 'pedestalSigma', 'pedestalSigmaCF',
				   'signalAndPedestal', 'signalAndPedestalCF', 'signalAndPedestalSigma', 'signalAndPedestalSigmaCF', 'signal',
				   'signalCF', 'deltaTimeCF']
BRANCHES1DTYPE = {'event': 'uint32', 'vetoedEvent': 'bool', 'badShape': 'int8', 'badPedestal': 'bool', 'satEvent': 'bool',
				  'voltageDia': 'float32', 'voltageHV': 'float32', 'currentHV': 'float32', 'timeHV.AsDouble()': 'float64',
				  'timeHV.Convert()': 'uint32', 'peakPosition': 'float32', 'peakPositionCF': 'float32', 'pedestal': 'float32',
				  'pedestalCF': 'float32', 'pedestalSigma': 'float32', 'pedestalSigmaCF': 'float32', 'signalAndPedestal': 'float32',
				  'signalAndPedestalCF': 'float32','signalAndPedestalSigma': 'float32','signalAndPedestalSigmaCF': 'float32',
				  'signal': 'float32', 'signalCF': 'float32', 'deltaTimeCF': 'float64'}
BRANCHESWAVESTOTAL = ['time', 'voltageSignal', 'voltageTrigger', 'voltageVeto']
BRANCHESWAVESTYPE = {'time': 'float64', 'voltageSignal': 'float64', 'voltageTrigger': 'float64', 'voltageVeto': 'float64'}
BRANCHES1DLOAD = ['event', 'voltageDia', 'voltageHV','currentHV', 'timeHV.Convert()', 'timeHV.AsDouble()', 'peakPosition', 'peakPositionCF', 'deltaTimeCF']
BRANCHESWAVESLOAD = ['time', 'voltageSignal']
ANALYSISSCALARBRANCHES = ['pedestal', 'pedestalCF', 'pedestalSigma', 'pedestalSigmaCF', 'signalAndPedestal', 'signalAndPedestalCF', 'signalAndPedestalSigma', 'signalAndPedestalSigmaCF', 'signal', 'signalCF']
fit_method = ('Minuit2', 'Migrad', )
LOAD_IGNORE_NAMES = ['analysis', 'pedestal', 'waveform', 'voltage', 'signal', 'dist', 'currents', 'ph', 'veto', 'trigger', 'peak', 'blag']

ro.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(100000)
ro.Math.MinimizerOptions.SetDefaultTolerance(0.000001)
ro.Math.MinimizerOptions.SetDefaultMinimizer('Minuit2', 'Migrad')

class AnalysisCaenCCD:
	def __init__(self, directory='.', config='CAENAnalysisConfig.cfg', infile='', bias=0.0, overw=False, verbose=False, doDebug=False):
		print 'Starting CCD Analysis ...'
		self.config = config
		self.verb = verbose
		self.overw = overw
		self.doDebug = doDebug
		self.inputFile = ''
		self.inDir = directory
		self.in_root_file = None
		self.in_tree_name = ''
		self.in_root_tree = None
		self.settings = None
		self.signal_ch = None
		self.trigger_ch = None
		self.veto_ch = None
		self.max_events = 0
		self.is_cal_run = False
		self.ch_caen_signal = 3
		self.cal_run_type = ''
		self.bias = bias
		self.emptyAnalysis = False
		self.pedestalIntegrationTime = 0.4e-6
		self.pedestalTEndPos = -20e-9
		self.peakTime = 2.14e-6
		self.peakTimeCF = 2.14e-6
		self.doPeakPos = True
		self.peakForward = self.pedestalIntegrationTime / 2.0
		self.peakBackward = self.pedestalIntegrationTime / 2.0
		self.doBadPedestalCut = True
		self.badShapeCut = 2
		self.doVetoedEventCut = True
		self.doSatCut = True
		self.peakTimeCut = 50e-9
		self.currentCut = 10e-9
		self.voltageDiaMaxOffset = 10  # value in
		self.voltageDiaMaxSigmas = 3  # value in sigmas
		self.voltageDiaMean = self.bias
		self.voltageDiaSpread = 0
		self.fit_min = -10000000
		self.fit_max = -10000000
		self.pedestal_sigma = 0
		self.pedestal_sigma_err = 0
		self.pedestal_vcal_sigma = 0
		self.pedestal_charge_sigma = 0
		self.pedestal_sigmaCF = 0
		self.pedestal_sigma_errCF = 0
		self.pedestal_vcal_sigmaCF = 0
		self.pedestal_charge_sigmaCF = 0

		self.analysisTreeExisted = False

		self.ptsWave, self.event, self.events, self.max_events = 0, np.zeros(1, 'I'), 0, 0
		self.eventVect = np.empty(0, 'f8')
		self.timeVect, self.signalWaveVect, self.triggerWaveVect, self.vetoWaveVect = np.empty(0, 'f8'), np.empty(0, 'f8'), np.empty(0, 'f8'), np.empty(0, 'f8')

		self.pedVect, self.pedSigmaVect, self.sigAndPedVect, self.sigAndPedSigmaVect, self.sigVect = None, None, None, None, None
		self.pedVectCF, self.pedSigmaVectCF, self.sigAndPedVectCF, self.sigAndPedSigmaVectCF, self.sigVectCF = None, None, None, None, None
		self.ped, self.pedSigma, self.sigAndPed, self.sigAndPedSigma, self.sig = np.zeros(1, 'f'), np.zeros(1, 'f'), np.zeros(1, 'f'), np.zeros(1, 'f'), np.zeros(1, 'f')
		self.pedCF, self.pedSigmaCF, self.sigAndPedCF, self.sigAndPedSigmaCF, self.sigCF = np.zeros(1, 'f'), np.zeros(1, 'f'), np.zeros(1, 'f'), np.zeros(1, 'f'), np.zeros(1, 'f')
		self.vetoedEvent, self.badShape, self.badPedestal = np.empty(0, '?'), np.empty(0, np.dtype('int8')), np.empty(0, '?')
		self.voltageHV, self.currentHV = np.empty(0, 'f8'), np.empty(0, 'f8')
		self.timeHV = np.empty(0, 'f8')

		self.signalWaveMeanVect, self.signalWaveSigmaVect = None, None

		self.cut0 = ro.TCut('cut0', '')
		self.cut0CF = ro.TCut('cut0CF', '')

		self.branches1DTotal = BRANCHES1DTOTAL[:]
		self.branches1DType = BRANCHES1DTYPE.copy()
		self.branchesWavesTotal = BRANCHESWAVESTOTAL[:]
		self.branchesWavesType = BRANCHESWAVESTYPE.copy()
		self.branchesAll = self.branches1DTotal + self.branchesWavesTotal

		self.branches1DLoad = BRANCHES1DLOAD[:]
		self.branchesWavesLoad = BRANCHESWAVESLOAD[:]
		self.dic1DVectLoaded = {}
		self.dicWavesVectLoaded = {}

		self.analysisScalarsBranches = ANALYSISSCALARBRANCHES[:]

		self.dicBraVect1D = OrderedDict()
		self.dicBraVectWaves = OrderedDict()

		self.hasBranch = {}

		self.pedestalTimeIndices = None  # has the indices for the pedestal for each event
		self.pedestalTimeIndicesCF = None  # has the indices for the pedestal for each event
		self.peak_positions = None  # has the peak position in time for the peak of the signal for each event
		self.peak_positionsCF = None  # has the peak position in time for the peak of the signal for each event
		self.peak_position = np.zeros(1, 'f')
		self.peak_positionCF = np.zeros(1, 'f')
		self.signalTimeIndices = None  # has the indices for the integral of the signal for each event
		self.signalTimeIndicesCF = None  # has the indices for the integral of the signal for each event

		self.cut0 = ro.TCut('cut0', '')

		self.utils = Utils()

		self.outDir = ''

		self.canvas = {}
		self.profile = {}
		self.histo = {}
		self.langaus = {}
		self.graph = {}
		self.line = {}
		self.random = None
		self.toy_histos = []

		self.signal_cal_folder = ''
		self.signal_cal_fit_params = {'p0': 0, 'p1': 0, 'p0_error': 0, 'p1_error': 0, 'prob': 0, 'chi2': 0, 'ndf': 0}
		self.cal_circuit_settings = ''
		self.vcal_to_q = VcalToElectrons()

		if infile == '' and directory != '.':
			print 'Is analysis of data after 06/18...'
			self.LoadInputTree()
			self.LoadPickles()
			self.outDir = self.inDir
			self.inputFile = self.in_tree_name + '.root'
			self.SetFromSettingsFile()

		elif infile != '' and directory == '.':
			print 'Is analysis of data before 06/18...'
			self.outDir, self.inputFile = '/'.join(infile.split('/')[:-1]), infile.split('/')[-1]
			self.in_tree_name = '.'.join(self.inputFile.split('.')[:-1])
			self.inDir = self.outDir
			self.LoadInputTree()
			self.LoadPickles()
		else:
			ExitMessage('I don\'t know what to do. If you want to run the analysis for old files, give inputFile, configFile, bias. If you want to run the analysis for new files, just give the directory where the pickles, data files and root files are.', os.EX_CONFIG)

		self.analysisFile = None
		self.analysisTree = None
		self.analysisTreeNameStem = self.in_tree_name + '.analysis'
		self.analysisTreeName = self.analysisTreeNameStem

		self.load_entries = 100000  # needed for sparing memory while loading vectors.
		self.start_entry = 0
		self.loaded_entries = 0
		self.suffix = None

		# self.delta_v_signal = self.signal_ch.adc_to_volts_cal['p1'] * 100 * 1000  # in mV
		self.delta_v_signal = 5

		self.Load_Config_File()
		#
		self.timeCFDelay = 1800 # delay time for CF in ns
		self.attenCF = 75 # percentage for attenuation in CF calculation
		self.time_CF_deltas = np.empty(1, 'f8')
		self.timeCFDeltasPeak = 0
		self.timeCFDeltasCut = 500e-9
		# self.timeCF = np.empty(self.settings.points, 'f8')
		# self.timeCFVect = np.empty(self.settings.points, 'f8')
		#

	def __del__(self):
		self.cut0.Delete()
		self.cut0CF.Delete()
		self.CloseInputROOTFiles()
		self.CloseAnalysisROOTFile()
		for key in self.line.keys():
			del self.line[key]
		for key in self.graph.keys():
			del self.graph[key]
		for key in self.langaus.keys():
			del self.langaus[key]
		for key in self.histo.keys():
			del self.histo[key]
		for key in self.profile.keys():
			del self.profile[key]
		for key in self.canvas.keys():
			self.canvas[key].Close()
			del self.canvas[key]
		print 'Deleted analysis object\n'

	def SetRandomGenerator(self, seed=0):
		seedi = seed if seed != 0 else int(time.time() % 1000000)
		print 'Using seed', seedi, 'for random generator TRandom3'
		self.random = ro.TRandom3(seedi)
		ro.gRandom = self.random

	def Load_Config_File(self):
		parser = ConfigParser()
		if os.path.isfile(self.config):
			print 'Reading configuration file:', self.config, '...', ; sys.stdout.flush()
			parser.read(self.config)

			if parser.has_section('ANALYSIS'):
				if parser.has_option('ANALYSIS', 'bias'):
					self.bias = parser.getfloat('ANALYSIS', 'bias')
				if parser.has_option('ANALYSIS', 'input_file') and self.inputFile == '':
					self.inputFile = parser.get('ANALYSIS', 'input_file')
				if parser.has_option('ANALYSIS', 'out_dir') and self.outDir == '':
					self.outDir = parser.get('ANALYSIS', 'out_dir')
				if parser.has_option('ANALYSIS', 'max_events'):
					self.max_events = parser.getint('ANALYSIS', 'max_events')
				if parser.has_option('ANALYSIS', 'peak_time'):
					self.peakTime = parser.getfloat('ANALYSIS', 'peak_time') * 1e-6
				if parser.has_option('ANALYSIS', 'integration_time'):
					self.pedestalIntegrationTime = parser.getfloat('ANALYSIS', 'integration_time') * 1e-6
				if parser.has_option('ANALYSIS', 'do_peak_positioning'):
					self.doPeakPos = parser.getboolean('ANALYSIS', 'do_peak_positioning')
				if parser.has_option('ANALYSIS', 'transition_time'):
					self.pedestalTEndPos = parser.getfloat('ANALYSIS', 'transition_time') * 1e-9
				if parser.has_option('ANALYSIS', 'forward_bakward_ratio'):
					self.peakForward = self.pedestalIntegrationTime / (1.0 + 1.0 / parser.getfloat('ANALYSIS', 'forward_bakward_ratio'))
					self.peakBackward = self.pedestalIntegrationTime / (1.0 + parser.getfloat('ANALYSIS', 'forward_bakward_ratio'))
				if parser.has_option('ANALYSIS', 'fit_min'):
					self.fit_min = parser.getfloat('ANALYSIS', 'fit_min')
				if parser.has_option('ANALYSIS', 'fit_max'):
					self.fit_max = parser.getfloat('ANALYSIS', 'fit_max')
				if parser.has_option('ANALYSIS', 'delta_v'):
					self.delta_v_signal = parser.getfloat('ANALYSIS', 'delta_v')

			if parser.has_section('CALIBRATION'):
				if parser.has_option('CALIBRATION', 'signal_cal_folder'):
					self.signal_cal_folder = Correct_Path(parser.get('CALIBRATION', 'signal_cal_folder'))
				if parser.has_option('CALIBRATION', 'cal_circuit_settings'):
					self.cal_circuit_settings = Correct_Path(parser.get('CALIBRATION', 'cal_circuit_settings'))

			if parser.has_section('CUTS'):
				if parser.has_option('CUTS', 'bad_pedestal'):
					self.doBadPedestalCut = parser.getboolean('CUTS', 'bad_pedestal')
				if parser.has_option('CUTS', 'bad_shape'):
					self.badShapeCut = parser.getint('CUTS', 'bad_shape')
				if parser.has_option('CUTS', 'vetoed_events'):
					self.doVetoedEventCut = parser.getboolean('CUTS', 'vetoed_events')
				if parser.has_option('CUTS', 'peak_position'):
					self.peakTimeCut = parser.getfloat('CUTS', 'peak_position') * 1e-6
				if parser.has_option('CUTS', 'current_cut'):
					self.currentCut = parser.getfloat('CUTS', 'current_cut') * 1e-9
				if parser.has_option('CUTS', 'sat_events'):
					self.doSatCut = parser.getboolean('CUTS', 'sat_events')
				if parser.has_option('CUTS', 'dia_voltage_offset'):
					self.voltageDiaMaxOffset = parser.getfloat('CUTS', 'dia_voltage_offset')
				if parser.has_option('CUTS', 'dia_voltage_spread'):
					self.voltageDiaMaxSigmas = parser.getfloat('CUTS', 'dia_voltage_spread')
			print 'Done'

	def SetFromSettingsFile(self):
		if self.settings:
			self.is_cal_run = self.settings.is_cal_run if 'is_cal_run' in self.settings.__dict__.keys() else False
			self.bias = self.settings.bias if not self.is_cal_run else self.settings.pulser_amplitude
			self.ch_caen_signal = self.settings.sigCh
			if self.is_cal_run:
				self.cal_run_type = 'in' if 'in' in self.in_tree_name else 'out'

	def LoadInputTree(self):
		if not os.path.isdir(self.inDir):
			ExitMessage('The directory {d} does not exist! Exiting...'.format(d=self.inDir), os.EX_DATAERR)
		if self.in_tree_name == '':
			root_files = glob.glob('{d}/*.root'.format(d=self.inDir))
			for name in LOAD_IGNORE_NAMES:
				root_files = [filei for filei in root_files if name not in filei.split('/')[-1].lower()]
			if len(root_files) == 1:
				self.in_tree_name = root_files[0].split('/')[-1].split('.root')[0]
			elif len(root_files) == 0:
				ExitMessage('There is no root file inside the directory {d}. Exiting...'.format(d=self.inDir), os.EX_DATAERR)
			else:
				print 'The following files were encountered:'
				for it in root_files:
					print it
				file_to_open = raw_input('Copy and paste the one that you want to open (should be input converted root file): ')
				if file_to_open in root_files:
					self.in_tree_name = file_to_open.split('/')[-1].split('.root')[0]
				else:
					ExitMessage('The given file: {f} is not part of the given list. Exiting...'.format(f=file_to_open), os.EX_DATAERR)
		if os.path.isfile('{d}/{f}.root'.format(d=self.inDir, f=self.in_tree_name)):
			self.in_root_file = ro.TFile('{d}/{f}.root'.format(d=self.inDir, f=self.in_tree_name))
			if self.in_root_file.GetListOfKeys().GetSize() == 0:
				print 'The file {d}/{f}.root is empty. Run won\'t be analysed.'.format(d=self.inDir, f=self.in_tree_name)
			else:
				self.in_root_tree = self.in_root_file.Get(self.in_tree_name)
				print 'Loaded raw tree'
		else:
			ExitMessage('The file {f} is not inside {d}. Exiting...'.format(d=self.inDir, f=self.in_tree_name), os.EX_DATAERR)

	def LoadPickles(self):
		if self.in_tree_name != '':
			if os.path.isdir(self.inDir):
				if os.path.isfile('{d}/{f}.settings'.format(d=self.inDir, f=self.in_tree_name)):
					with open('{d}/{f}.settings'.format(d=self.inDir, f=self.in_tree_name), 'rb') as pklsets:
						self.settings = pickle.load(pklsets)
				if os.path.isfile('{d}/{f}.signal_ch'.format(d=self.inDir, f=self.in_tree_name)):
					with open('{d}/{f}.signal_ch'.format(d=self.inDir, f=self.in_tree_name)) as pklsign:
						self.signal_ch = pickle.load(pklsign)
				if os.path.isfile('{d}/{f}.trigger_ch'.format(d=self.inDir, f=self.in_tree_name)):
					with open('{d}/{f}.trigger_ch'.format(d=self.inDir, f=self.in_tree_name)) as pkltrig:
						self.trigger_ch = pickle.load(pkltrig)
				if os.path.isfile('{d}/{f}.veto_ch'.format(d=self.inDir, f=self.in_tree_name)):
					with open('{d}/{f}.veto_ch'.format(d=self.inDir, f=self.in_tree_name)) as pklveto:
						self.veto_ch = pickle.load(pklveto)
				elif os.path.isfile('{d}/{f}.veto'.format(d=self.inDir, f=self.in_tree_name)):
					with open('{d}/{f}.veto'.format(d=self.inDir, f=self.in_tree_name)) as pklveto:
						self.veto_ch = pickle.load(pklveto)

	def CheckLoadedVectors(self):
		if self.eventVect.size == 0:
			print 'There are no waves that pass the cuts:', self.cut0.GetTitle(), '. Exiting'
			self.CloseAnalysisROOTFile()
			self.CloseInputROOTFiles()
			self.emptyAnalysis = True

	def AnalysisWaves(self, doCuts0=True):
		optiont = 'RECREATE' if self.overw else 'UPDATE'
		self.OpenAnalysisROOTFile(optiont)
		if doCuts0: self.CreateCut0()
		if not self.hasBranch['peakPosition'] or not np.array([self.hasBranch[key0] for key0 in self.analysisScalarsBranches]).all():
			if self.in_root_tree:
				self.suffix = None if self.in_root_tree.GetEntries() <= self.load_entries else self.start_entry
				while self.loaded_entries < self.in_root_tree.GetEntries():
					self.CloseAnalysisROOTFile()
					self.OpenAnalysisROOTFile('RECREATE')
					self.LoadVectorsFromTree()
					if self.doDebug: ipdb.set_trace()
					self.ExplicitVectorsFromDictionary()
					if self.doDebug: ipdb.set_trace()
					self.CheckLoadedVectors()
					if self.emptyAnalysis:
						break
					self.DoConstantFraction()
					if self.doDebug: ipdb.set_trace()
					if self.doPeakPos:
						self.FindRealPeakPosition()
						if self.doDebug: ipdb.set_trace()
					else:
						self.peak_positions = np.full(self.events, self.peakTime)
						self.peak_positionsCF = np.full(self.events, self.peakTimeCF)
					if self.doDebug: ipdb.set_trace()
					self.FindPedestalPosition()
					if self.doDebug: ipdb.set_trace()
					self.FindSignalPositions(self.peakBackward, self.peakForward)
					if self.doDebug: ipdb.set_trace()
					self.CalculatePedestalsAndSignals()
					if self.doDebug: ipdb.set_trace()
					self.FillAnalysisBranches()
					self.CloseAnalysisROOTFile()
					self.CloseInputROOTFiles()
					self.Reset_Braches_Lists_And_Dictionaries()
					self.LoadInputTree()
					self.OpenAnalysisROOTFile('READ')
					if self.suffix >= 0:
						self.suffix += 1
				if self.suffix >= 0:
					self.CloseAnalysisROOTFile()
					if self.suffix >=0:
						print 'Merging analysis files...'
						if os.path.isfile('{d}/{f}.root'.format(d=self.outDir, f=self.analysisTreeNameStem)):
							os.remove('{d}/{f}.root'.format(d=self.outDir, f=self.analysisTreeNameStem))
						mergep = subp.Popen(['hadd', '-k', '-O', '-n', '0', '{d}/{f}.root'.format(d=self.outDir, f=self.analysisTreeNameStem)] + ['{d}/{f}{s}.root'.format(d=self.outDir, f=self.analysisTreeNameStem, s=si) for si in xrange(self.suffix)],  bufsize=-1, stdin=subp.PIPE, stdout=subp.PIPE, close_fds=True)
						while mergep.poll() is None:
							time.sleep(2)
						pid = mergep.pid
						mergep.stdin.close()
						mergep.stdout.close()
						if mergep.wait() is None:
							mergep.kill()
						try:
							os.kill(pid, 15)
						except OSError:
							pass
						del mergep
						print 'Finished merging analysis files'
						for si in xrange(self.suffix):
							os.remove('{d}/{f}{s}.root'.format(d=self.outDir, f=self.analysisTreeNameStem, s=si))
						print 'Merging CF files...'
						if os.path.isfile('{d}/{f}.timeCF.{de}.{att}.root'.format(d=self.outDir, f=self.analysisTreeNameStem, de=self.timeCFDelay, att=self.attenCF)):
							os.remove('{d}/{f}.timeCF.{de}.{att}.root'.format(d=self.outDir, f=self.analysisTreeNameStem, de=self.timeCFDelay, att=self.attenCF))
						mergep = subp.Popen(['hadd', '-k', '-O', '-n', '0', '{d}/{f}.timeCF.{de}.{att}.root'.format(d=self.outDir, f=self.analysisTreeNameStem, de=self.timeCFDelay, att=self.attenCF)] + ['{d}/{f}{s}.timeCF.{de}.{att}.root'.format(d=self.outDir, f=self.analysisTreeNameStem, s=si, de=self.timeCFDelay, att=self.attenCF) for si in xrange(self.suffix)],  bufsize=-1, stdin=subp.PIPE, stdout=subp.PIPE, close_fds=True)
						while mergep.poll() is None:
							time.sleep(2)
						pid = mergep.pid
						mergep.stdin.close()
						mergep.stdout.close()
						if mergep.wait() is None:
							mergep.kill()
						try:
							os.kill(pid, 15)
						except OSError:
							pass
						del mergep
						print 'Finished merging CF files'
						for si in xrange(self.suffix):
							os.remove('{d}/{f}{s}.timeCF.{de}.{att}.root'.format(d=self.outDir, f=self.analysisTreeNameStem, s=si, de=self.timeCFDelay, att=self.attenCF))
					self.suffix = None
					self.OpenAnalysisROOTFile('READ')
		if self.emptyAnalysis: return
		if not self.is_cal_run or self.cal_run_type == 'out':
			if self.doDebug: ipdb.set_trace()
			if self.hasBranch['deltaTimeCF']:
				self.PlotCFDeltaDistribution()
				self.AddCFShiftCut()
			if self.hasBranch['peakPosition']:
				self.PlotPeakPositionDistributions('peakPosDist', isCF=False)
				self.AddPeakPositionCut()
			if self.hasBranch['peakPositionCF']:
				self.PlotPeakPositionDistributions('peakPosDistCF', isCF=True)
				self.AddPeakPositionCutCF()
		if self.hasBranch['voltageDia']:
			if self.doDebug: ipdb.set_trace()
			self.PlotHVDiaDistribution(cut=self.cut0.GetTitle())
			self.AddDiamondVoltageCut()

	def OpenAnalysisROOTFile(self, mode='READ'):
		self.analysisTreeName = self.analysisTreeNameStem + str(self.suffix) if self.suffix >= 0 else self.analysisTreeNameStem
		if not os.path.isdir(self.outDir):
			ExitMessage('The directory {d} does not exist. Exiting...'.format(d=self.outDir), os.EX_DATAERR)

		if self.analysisFile:
			if self.analysisFile.IsOpen():
				if self.analysisFile.GetOption().lower() != mode.lower():
					if self.analysisFile.ReOpen(mode) == -1:
						ExitMessage('Could not reopen file {f}.root in mode: {m}. Exiting...'.format(f=self.analysisTreeName, m=mode))
				return
			else:
				self.analysisFile = None
				self.analysisFile = ro.TFile('{o}/{f}.root'.format(o=self.outDir, f=self.analysisTreeName), mode)
				self.LoadAnalysisTree()
		else:
			if mode.lower() in ['new', 'create', 'update', 'recreate']:
				mode2 = mode if os.path.isfile('{d}/{f}.root'.format(d=self.outDir, f=self.analysisTreeName)) else 'RECREATE'
				self.analysisFile = ro.TFile('{o}/{f}.root'.format(o=self.outDir, f=self.analysisTreeName), mode2)
				self.LoadAnalysisTree()
			else:
				if os.path.isfile('{d}/{f}.root'.format(d=self.outDir, f=self.analysisTreeName)):
					self.analysisFile = ro.TFile('{o}/{f}.root'.format(o=self.outDir, f=self.analysisTreeName), mode)
					self.LoadAnalysisTree()
				else:
					ExitMessage('Can\'t open the file {f}.root in {m} mode because it does not exist!!! Exiting...'.format(f=self.analysisTreeName, m=mode), os.EX_DATAERR)

	def LoadAnalysisTree(self, mode='UPDATE'):
		self.analysisTreeExisted = True
		if self.analysisFile:
			if not self.analysisFile.IsOpen():
				print 'Analysis file is closed. Opening it in {m} mode...'.format(m=mode.lower())
				self.OpenAnalysisROOTFile(mode)
		else:
			print 'The file has not been opened. Call method OpenAnalysisROOTFile first'
			return
		self.analysisFile.cd()
		self.analysisTree = self.analysisFile.Get(self.analysisTreeName)
		if not self.analysisTree:
			self.analysisTreeExisted = False
			self.analysisTree = ro.TTree(self.analysisTreeNameStem, self.analysisTreeNameStem)
		# else:
		# 	if not self.analysisTree.GetFriend(self.in_tree_name):
		# 		self.analysisTree.AddFriend(self.in_tree_name, '{d}/{f}.root'.format(d=self.inDir, f=self.in_tree_name))
		self.AddExistingFriends()
		self.hasBranch = {branch: self.TreeHasBranch(branch) for branch in self.branchesAll}
		self.IsTimeHVaTimeStamp()
		self.UpdateBranchesLists()

	def AddExistingFriends(self):
		if self.analysisTree:
			if not self.analysisTree.GetFriend('vcalTree'):
				if os.path.isfile('{d}/{f}.vcal.root'.format(d=self.outDir, f=self.analysisTreeName)) and not self.overw:
					self.AddVcalFriend(False, False)
			if not self.analysisTree.GetFriend('vcalTreeCF'):
				if os.path.isfile('{d}/{f}.vcalCF.root'.format(d=self.outDir, f=self.analysisTreeName)) and not self.overw:
					self.AddVcalFriend(False, True)
			if not self.analysisTree.GetFriend('chargeTree'):
				if os.path.isfile('{d}/{f}.charge.root'.format(d=self.outDir, f=self.analysisTreeName)) and not self.overw:
					self.AddChargeFriend(False, False)
			if not self.analysisTree.GetFriend('chargeTreeCF'):
				if os.path.isfile('{d}/{f}.chargeCF.root'.format(d=self.outDir, f=self.analysisTreeName)) and not self.overw:
					self.AddChargeFriend(False, True)
			if not self.analysisTree.GetFriend('timeCFTree'):
				if os.path.isfile('{d}/{f}.timeCF.{de}.{att}.root'.format(d=self.outDir, f=self.analysisTreeName, de=self.timeCFDelay, att=self.attenCF)) and not self.overw:
					self.AddTimeCFFriend(self.timeCFDelay, self.attenCF, False)

	def TreeHasBranch(self, branch):
		if self.in_root_tree:
			if self.in_root_tree.GetBranch(branch) or self.analysisTree.GetBranch(branch):
				return True
		return False

	def IsTimeHVaTimeStamp(self):
		if self.hasBranch['timeHV']:
			if self.in_root_tree.GetLeaf('timeHV').GetTypeName() != 'TDatime':
				# self.branches1DLoad = ['timeHV.AsDouble()' if branch == 'timeHV' else branch for branch in self.branches1DLoad]
				if 'timeHV.Convert()' in self.branches1DLoad: self.branches1DLoad.remove('timeHV.Convert()')
				# del self.branches1DType['timeHV']
				# self.branches1DType['timeHV.AsDouble()'] = 'float64'
				if self.branches1DType.has_key('timeHV.Convert()'): del self.branches1DType['timeHV.Convert()']
			else:
				# self.branches1DLoad = ['timeHV.Convert()' if branch == 'timeHV' else branch for branch in self.branches1DLoad]
				if 'timeHV.AsDouble()' in self.branches1DLoad: self.branches1DLoad.remove('timeHV.AsDouble()')
				# del self.branches1DType['timeHV']
				# self.branches1DType['timeHV.Convert()'] = 'uint32'
				if self.branches1DType.has_key('timeHV.AsDouble()'): del self.branches1DType['timeHV.AsDouble()']

	def UpdateBranchesLists(self):
		for branch in self.branches1DTotal[:]:
			if not self.hasBranch[branch]:
				self.branches1DTotal.remove(branch)
		for branch in self.branches1DLoad[:]:
			if branch.startswith('timeHV'):
				if not self.hasBranch['timeHV']:
					self.branches1DLoad.remove(branch)
				else:
					if self.dicBraVect1D.has_key(branch):
						self.dic1DVectLoaded[branch] = True if self.dicBraVect1D[branch] else False
					else:
						self.dic1DVectLoaded[branch] = False
			elif not self.hasBranch[branch]:
				self.branches1DLoad.remove(branch)
			else:
				if self.dicBraVect1D.has_key(branch):
					self.dic1DVectLoaded[branch] = True if self.dicBraVect1D[branch] else False
				else:
					self.dic1DVectLoaded[branch] = False

		for branch in self.branchesWavesTotal[:]:
			if not self.hasBranch[branch]:
				self.branchesWavesTotal.remove(branch)
		for branch in self.branchesWavesLoad[:]:
			if not self.hasBranch[branch]:
				self.branchesWavesLoad.remove(branch)
			else:
				if self.dicBraVectWaves.has_key(branch):
					self.dicWavesVectLoaded[branch] = True if self.dicBraVectWaves[branch] else False
				else:
					self.dicWavesVectLoaded[branch] = False

	def CreateCut0(self):
		if self.cut0.GetTitle() != '':
			self.cut0.SetTitle('')
		tempCut = self.ReturnBasicCut0()
		self.cut0 += tempCut
		if self.cut0CF.GetTitle() != '':
			self.cut0CF.SetTitle('')
		self.cut0CF += tempCut

	def ReturnBasicCut0(self):
		tempCut = ro.TCut('tempCut', '')
		if self.doBadPedestalCut and 'badPedestal' in self.branches1DTotal:
			tempCut += ro.TCut('badPedCut', 'badPedestal==0')
		if self.doVetoedEventCut and 'vetoedEvent' in self.branches1DTotal:
			tempCut += ro.TCut('vetoedEventCut', 'vetoedEvent==0')
		if self.badShapeCut == 1 and 'badShape' in self.branches1DTotal:
			tempCut += ro.TCut('badShapeCut', 'badShape!=1')
		elif self.badShapeCut == 2 and 'badShape' in self.branches1DTotal:
			tempCut += ro.TCut('badShapeCut', 'badShape==0')
		if self.doSatCut and 'satEvent' in self.branches1DTotal:
			tempCut += ro.TCut('satEventCut', 'satEvent==0')
		if 'currentHV' in self.branches1DTotal:
			tempCut += ro.TCut('currentCut', 'abs(currentHV)<{cc}'.format(cc=self.currentCut))
		return tempCut

	def ResetCut0(self):
		self.cut0.Clear()
		self.cut0 = ro.TCut('cut0', '')
		self.cut0CF.Clear()
		self.cut0CF = ro.TCut('cut0CF', '')

	def LoadVectorsFromTree(self):
		self.max_events = self.in_root_tree.GetEntries() if self.max_events == 0 else self.max_events
		self.max_events = self.load_entries if self.max_events > self.load_entries else self.max_events
		self.ptsWave = self.in_root_tree.GetLeaf('time').GetLen()
		working_tree = self.analysisTree if self.analysisTreeExisted else self.in_root_tree
		branches_to_load_1D = [branch for branch in self.branches1DLoad if not self.dic1DVectLoaded[branch]]
		options = 'goff' if len(branches_to_load_1D) == 1 else 'goff para'
		if len(branches_to_load_1D) > 0:
			leng = working_tree.Draw(':'.join(branches_to_load_1D), self.cut0, options, self.max_events, self.start_entry)
			if leng == -1:
				print 'Error, could not load the branches: {b}. Try again :('.format(b=':'.join(branches_to_load_1D))
				return
			while leng > working_tree.GetEstimate():
				working_tree.SetEstimate(leng)
				leng = working_tree.Draw(':'.join(branches_to_load_1D), self.cut0, options, self.max_events, self.start_entry)
			self.events = leng
			for pos, branch in enumerate(branches_to_load_1D):
				if self.verb: print 'Vectorising branch:', branch, '...', ; sys.stdout.flush()
				temp = working_tree.GetVal(pos)
				self.dicBraVect1D[branch] = np.array([temp[ev] for ev in xrange(self.events)], dtype=np.dtype(self.branches1DType[branch]))
				self.dic1DVectLoaded[branch] = True
				if self.verb: print 'Done'
				del temp

		branches_to_load_waves = [branch for branch in self.branchesWavesLoad if not self.dicWavesVectLoaded[branch]]
		if len(branches_to_load_waves) > 0:
			leng = working_tree.Draw(':'.join(branches_to_load_waves), self.cut0, 'goff para', self.max_events, self.start_entry)
			if leng == -1:
				print 'Error, could not load the branches {b}. Try again :('.format(b=':'.join(branches_to_load_waves))
				return
			while leng > working_tree.GetEstimate():
				working_tree.SetEstimate(leng)
				leng = working_tree.Draw(':'.join(branches_to_load_waves), self.cut0, 'goff para', self.max_events, self.start_entry)
			for pos, branch in enumerate(branches_to_load_waves):
				if self.verb: print 'Vectorising branch:', branch, '...', ; sys.stdout.flush()
				temp = working_tree.GetVal(pos)
				self.dicBraVectWaves[branch] = np.array([[temp[ev * self.ptsWave + pt] for pt in xrange(self.ptsWave)] for ev in xrange(self.events)], dtype=np.dtype(self.branchesWavesType[branch]))
				self.dicWavesVectLoaded[branch] = True
				if self.verb: print 'Done'
				del temp

	def LoadSignalScalars(self):
		try:
			temp = self.eventVect[0] + self.sigVect[0] + self.pedVect[0]
			temp2 = self.eventVectCF[0] + self.sigVectCF[0] + self.pedVectCF[0]
		except Exception:
			branchesLoad = ['event', 'signal', 'pedestal', 'signalCF', 'pedestalCF']
			print 'Loading {b} branches...'.format(b=':'.join(branchesLoad)), ; sys.stdout.flush()
			tempCut = ro.TCut('cutScalars', '')
			tempCut += self.ReturnBasicCut0()
			leng = self.analysisTree.Draw(':'.join(branchesLoad), tempCut, 'goff para')
			if leng == -1:
				ExitMessage('Error, could not load the branches: {b}. Try again :('.format(b=':'.join(branchesLoad)))
			while leng > self.analysisTree.GetEstimate():
				self.analysisTree.SetEstimate(leng)
				leng = self.analysisTree.Draw(':'.join(branchesLoad), tempCut, 'goff para')
			events = leng
			dicBranches = {}
			for pos, branch in enumerate(branchesLoad):
				temp = self.analysisTree.GetVal(pos)
				dicBranches[branch] = np.array([temp[ev] for ev in xrange(events)], 'f4')
			self.eventVect = dicBranches['event'].astype('uint32')
			self.sigVect = dicBranches['signal'].astype('f4')
			self.sigVectCF = dicBranches['signalCF'].astype('f4')
			self.pedVect = dicBranches['pedestal'].astype('f4')
			self.pedVectCF = dicBranches['pedestalCF'].astype('f4')
			print 'Done'

	def ExplicitVectorsFromDictionary(self):
		if self.hasBranch['voltageSignal']:
			if self.dicWavesVectLoaded['voltageSignal']:
				self.signalWaveVect = self.dicBraVectWaves['voltageSignal']
		if self.hasBranch['time']:
			if self.dicWavesVectLoaded['time']:
				self.timeVect = self.dicBraVectWaves['time']
		if self.hasBranch['event']:
			if self.dic1DVectLoaded['event']:
				self.eventVect = self.dicBraVect1D['event']
		if self.hasBranch['voltageHV']:
			if self.dic1DVectLoaded['voltageHV']:
				self.voltageHV = self.dicBraVect1D['voltageHV']
		if self.hasBranch['currentHV']:
			if self.dic1DVectLoaded['currentHV']:
				self.currentHV = self.dicBraVect1D['currentHV']
		if self.hasBranch['timeHV']:
			key = 'timeHV.Convert()' if 'timeHV.Convert()' in self.branches1DLoad else 'timeHV.AsDouble()'
			if self.dic1DVectLoaded[key]:
				self.timeHV = self.dicBraVect1D[key]
		if self.hasBranch['peakPosition']:
			if self.dic1DVectLoaded['peakPosition']:
				self.peak_positions = self.dicBraVect1D['peakPosition']
		if self.hasBranch['peakPositionCF']:
			if self.dic1DVectLoaded['peakPositionCF']:
				self.peak_positionsCF = self.dicBraVect1D['peakPositionCF']

	def DoConstantFraction(self):
		self.AddTimeCFFriend(self.timeCFDelay, self.attenCF, self.overw)

	def FindRealPeakPosition(self):
		def GetConcavity(x, fit_pol4):
			return 12. * fit_pol4.GetParameter(4) * x ** 2 + 6 * fit_pol4.GetParameter(3) * x + 2 * fit_pol4.GetParameter(2)

		print 'Getting real peak positions...'
		# mpos = self.signalWaveVect.argmin(axis=1) if self.bias >= 0 else self.signalWaveVect.argmax(axis=1)
		# time_mpos = self.timeVect[:, mpos].diagonal()
		# time_mpos = np.array([self.timeVect[it[0], pos] for it, pos in np.ndenumerate(mpos)])
		# xmin, xmax = time_mpos - self.pedestalIntegrationTime, time_mpos + self.pedestalIntegrationTime
		# par0lim = {'low': 1e-30, 'up': 1e30} if self.bias >= 0 else {'low': -1e30, 'up': -1e-30}
		# par0ini = 3.14 if self.bias >= 0 else -3.14
		ro.Math.MinimizerOptions.SetDefaultMinimizer(*fit_method)
		ro.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(10000)
		# ro.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(1000000)
		ro.Math.MinimizerOptions.SetDefaultTolerance(0.1)
		ro.gErrorIgnoreLevel = ro.kFatal
		self.peak_positions = []
		self.peak_positionsCF = []
		print 'Calculating peak positions...'
		self.utils.CreateProgressBar(len(self.timeVect))
		self.utils.bar.start()
		# fit_fcn = ro.TF1('fit_peak', '[0]*exp(-((x-[1])/[2])^2)+pol1(3)', self.peakTime - 2e-6 if not self.settings.is_cal_run else self.peakTime - 1e-6, self.peakTime + 2e-6 if not self.settings.is_cal_run else self.peakTime + 1e-6)
		fit_fcn = ro.TF1('fit_peak', 'pol4', self.peakTime - 1.5e-6, self.peakTime + 1.5e-6)
		fit_fcn.SetNpx(1000)
		for it, timei in enumerate(self.timeVect):
			# 	par0 = 1 if self.bias < 0 else -1
			par1 = self.peakTime
			# 	par2 = 1e-6
			# 	par3 = -1 if self.bias < 0 else 1
			# 	par4 = 0.1e-6 if self.bias < 0 else -0.1e-6
			# 	fit_fcn.SetParameter(0, -1)
			# 	fit_fcn.SetParLimits(0, -100, 100)
			# 	if not self.is_cal_run:
			# 		fit_fcn.SetParameter(0, par0)
			# 		fit_fcn.SetParLimits(0, 0 if self.bias < 0 else -100, 100 if self.bias < 0 else 0)
			par1Min, par1Max = par1 - 1.5e-6, par1 + 1.5e-6
			# 	fit_fcn.SetParameter(1, par1)
			# 	fit_fcn.SetParLimits(1, par1Min, par1Max)
			# 	fit_fcn.SetParameter(2, par2)
			# 	fit_fcn.SetParLimits(2, 0, 50e-6)
			# 	fit_fcn.SetParameter(3, 1)
			# 	if not self.is_cal_run:
			# 		fit_fcn.SetParameter(3, par3)
			# 	fit_fcn.SetParLimits(3, -100, 100)
			# 	fit_fcn.SetParameter(4, -0.1e-6)
			# 	fit_fcn.SetParLimits(4, -1e6, 1e6)
			# 	if not self.is_cal_run:
			# 		fit_fcn.SetParameter(4, par4)
			# 		fit_fcn.SetParLimits(4, 0 if self.bias < 0 else -1e6, 1e6 if self.bias < 0 else 0)
			xmin, xmax = par1Min, par1Max
			# fit_fcn = ro.TF1('fit_{it}'.format(it=it), '[0]*(x-[1])^2+[2]', xmin[it], xmax[it])
			# par2limFact = {'low': -10.0, 'up': -0.1} if self.bias >= 0 else {'low': 0.1, 'up': 10.0}
			# bin1, bin2 = int(mpos[it] - np.round(self.pedestalIntegrationTime / self.settings.time_res)), int(mpos[it] + np.round(self.pedestalIntegrationTime / self.settings.time_res))
			# y1, y2, y0 = self.signalWaveVect[it, bin1], self.signalWaveVect[it, bin2], self.signalWaveVect[it, mpos[it]]
			xpeak = self.peakTime
			fit_res = None
			grafit = ro.TGraph(len(timei), np.add(timei, self.time_CF_deltas[it], dtype='f8'), self.signalWaveVect[it])
			grafit.GetXaxis().SetRangeUser(xmin, xmax)
			ro.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(100000)
			ro.Math.MinimizerOptions.SetDefaultTolerance(0.000001)
			ro.Math.MinimizerOptions.SetDefaultMinimizer('Minuit2', 'Migrad')
			xmin, xmax = par1 - 1e-6, par1 + 1e-6
			fit_res = grafit.Fit('fit_peak', 'QBNM0S', '', xmin, xmax)
			fit_res = grafit.Fit('fit_peak', 'QBNM0S', '', xmin, xmax)
			for it2 in xrange(5):
				fit_res = grafit.Fit('fit_peak', 'QBNM0S', '', xmin, xmax)
				fit_res = grafit.Fit('fit_peak', 'QBNM0S', '', xmin, xmax)
				# y1, y2, y3 = fit_fcn.Eval(xmin), fit_fcn.Eval(xpeak), fit_fcn.Eval(xmax)
				# concavity = 1 if y2 < y1 and y2 < y3 else -1 if y2 > y1 and y2 > y3 else 0
				concavity = GetConcavity(xpeak, fit_fcn)
				ymax, ymin = fit_fcn.GetMaximumX(xmin, xmax), fit_fcn.GetMinimumX(xmin, xmax)
				# xpeak = fit_fcn.GetMaximumX(xmin, xmax) if self.bias < 0 else fit_fcn.GetMinimumX(xmin, xmax)
				# concavity1 = GetConcavity(max(xmin, xpeak - 0.5e-6), fit_fcn)
				# concavity3 = GetConcavity(max(xmin, xpeak + 0.5e-6), fit_fcn)
				xpeak = ymax if concavity < 0 else ymin if concavity > 0 else ymax if concavity == 0 and abs(xpeak - ymax) < abs(xpeak - ymin) else ymin
				xwindow = 900e-9 if it2 == 0 else 800e-9 if it2 == 1 else 700e-9 if it2 == 2 else 600e-9 if it2 == 3 else 500e-9
				if xpeak + xwindow > par1Max:
					xmin, xmax = par1Max - 2 * xwindow, par1Max
				elif xpeak - xwindow < par1Min:
					xmin, xmax = par1Min, par1Min + 2 * xwindow
				else:
					xmin, xmax = xpeak - xwindow, xpeak + xwindow
			if fit_res.IsValid() and (par1Min < xpeak < par1Max):
				self.peak_positionsCF.append(xpeak)
				self.peak_positions.append(xpeak - self.time_CF_deltas[it])
			else:
				self.peak_positionsCF.append(0)
				self.peak_positions.append(0)
			# print 'Could not find peak for event', it, '; setting peak position to 0 for this event'
			# fit_fcn.Delete()
			# fit.Delete()
			self.utils.bar.update(it + 1)
		ro.gErrorIgnoreLevel = ro.kInfo
		# fit = [ro.TGraph(len(timei), timei, self.signalWaveVect[it]).Fit('pol2', 'QMN0FS', '', xmin[it], xmax[it]) for it, timei in enumerate(self.timeVect)]
		# fit = [ro.TGraph(len(timei), timei, self.signalWaveVect[it]).Fit('fit_{it}'.format(it=it), 'QBMN0FS', '', xmin[it], xmax[it]) for it, timei in enumerate(self.timeVect)]
		# for it, timei in enumerate(self.timeVect):
		# 	pass
		# b, a = np.array([fiti.Parameter(1) for fiti in fit]), np.array([fiti.Parameter(2) for fiti in fit])
		# self.peak_positions = np.divide(-b, 2 * a)
		# self.peak_positions = np.array([fiti.Parameter(1) for fiti in fit])
		self.peak_positionsCF = np.array(self.peak_positionsCF)
		self.peak_positions = np.array(self.peak_positions)
		self.utils.bar.finish()
		print 'Done getting real peak positions'

	def FindRealPeakPosition2(self):
		print 'Getting real peak positions...'
		mpos = self.signalWaveVect.argmin(axis=1) if self.bias >= 0 else self.signalWaveVect.argmax(axis=1)
		# time_mpos = self.timeVect[:, mpos].diagonal()
		time_mpos = np.array([self.timeVect[it[0], pos] for it, pos in np.ndenumerate(mpos)])
		xmin, xmax = time_mpos - self.pedestalIntegrationTime, time_mpos + self.pedestalIntegrationTime
		par0lim = {'low': 1e-30, 'up': 1e30} if self.bias >= 0 else {'low': -1e30, 'up': -1e-30}
		par0ini = 3.14 if self.bias >= 0 else -3.14
		ro.Math.MinimizerOptions.SetDefaultMinimizer('Minuit2', 'Migrad')
		self.peak_positions = []
		print 'Calculating peak positions...'
		self.utils.CreateProgressBar(len(self.timeVect))
		self.utils.bar.start()
		for it, timei in enumerate(self.timeVect):
			# fit_fcn = ro.TF1('fit_{it}'.format(it=it), '[0]*(x-[1])^2+[2]', xmin[it], xmax[it])
			fit_fcn = ro.TF1('fit_{it}'.format(it=it), '[0]*exp(-((x-[1])/[2])^2+pol1(3)', self.peakTime - 1e-6, self.peakTime - 1e-6)
			# par2limFact = {'low': -10.0, 'up': -0.1} if self.bias >= 0 else {'low': 0.1, 'up': 10.0}
			bin1, bin2 = int(mpos[it] - np.round(self.pedestalIntegrationTime / self.settings.time_res)), int(mpos[it] + np.round(self.pedestalIntegrationTime / self.settings.time_res))
			y1, y2, y0 = self.signalWaveVect[it, bin1], self.signalWaveVect[it, bin2], self.signalWaveVect[it, mpos[it]]
			par0 = (y1 + y2 - 2 * y0) / (2 * self.pedestalIntegrationTime ** 2)
			par1 = time_mpos[it] + (self.pedestalIntegrationTime * (y2 - y1)) / (2 * (2 * y0 - y1 - y2))
			par2 = ((y1 - 4 * y0) ** 2 - 2 * (4 * y0 + y1) * y2 + y2 ** 2) / (8 * (2 * y0 - y1 - y2))
			fit_fcn.SetParameter(0, par0)
			fit_fcn.SetParLimits(0, par0/10.0, par0 * 10)
			fit_fcn.SetParameter(1, par1)
			fit_fcn.SetParLimits(1, time_mpos[it] - 2 * self.pedestalIntegrationTime, time_mpos[it] + 2 * self.pedestalIntegrationTime)
			fit_fcn.SetParameter(2, par2)
			fit_fcn.SetParLimits(2, -2.0 * abs(par2), 2.0 * abs(par2))
			fit = ro.TGraph(len(timei), timei, self.signalWaveVect[it]).Fit('fit_{it}'.format(it=it), 'QBMN0FS', '', xmin[it], xmax[it])
			self.peak_positions.append(fit.Parameter(1))
			# fit_fcn.Delete()
			# fit.Delete()
			self.utils.bar.update(it + 1)
		# fit = [ro.TGraph(len(timei), timei, self.signalWaveVect[it]).Fit('pol2', 'QMN0FS', '', xmin[it], xmax[it]) for it, timei in enumerate(self.timeVect)]
		# fit = [ro.TGraph(len(timei), timei, self.signalWaveVect[it]).Fit('fit_{it}'.format(it=it), 'QBMN0FS', '', xmin[it], xmax[it]) for it, timei in enumerate(self.timeVect)]
		# for it, timei in enumerate(self.timeVect):
		# 	pass
		# b, a = np.array([fiti.Parameter(1) for fiti in fit]), np.array([fiti.Parameter(2) for fiti in fit])
		# self.peak_positions = np.divide(-b, 2 * a)
		# self.peak_positions = np.array([fiti.Parameter(1) for fiti in fit])
		self.peak_positions = np.array(self.peak_positions)
		self.utils.bar.finish()
		print 'Done getting real peak positions'

	def FillTreePeakPositions(self):
		print 'Filling tree with peak positions...'
		peakPosBra = self.analysisTree.Branch('peakPosition', self.peak_position, 'peakPosition/F')
		entries = self.in_root_tree.GetEntries()
		self.CloseInputROOTFiles()
		self.utils.CreateProgressBar(entries)
		self.utils.bar.start()
		self.analysisFile.cd()
		for ev in xrange(entries):
			# self.in_root_tree.GetEntry(ev)
			if ev in self.eventVect:
				try:
					self.peak_position.itemset(self.peak_positions[np.argwhere(ev == self.eventVect).flatten()])
				except ValueError:
					ExitMessage('Could not fill event {ev}; it should have a peak position of: {v}. Exiting'.format(ev=ev, v=self.peak_positions[np.argwhere(ev == self.eventVect).flatten()]), os.EX_DATAERR)
			else:
				self.peak_position.itemset(0)
			# peakPosBra.Fill()
			self.analysisTree.Fill()
			self.utils.bar.update(ev + 1)
		self.analysisFile.cd()
		if not self.analysisTree.GetFriend(self.in_tree_name):
			self.analysisTree.AddFriend(self.in_tree_name, '{d}/{f}.root'.format(d=os.path.abspath(self.inDir), f=self.in_tree_name))
		self.analysisTree.Write()
		self.utils.bar.finish()

	def CloseAnalysisROOTFile(self):
		if self.analysisFile:
			if self.analysisFile.IsOpen():
				self.analysisFile.Close()
			if self.analysisTree:
				del self.analysisTree
		self.analysisTree = None
		self.analysisFile = None

	def CloseInputROOTFiles(self):
		if self.in_root_file:
			if self.in_root_file.IsOpen():
				self.in_root_file.Close()
			if self.in_root_tree:
				del self.in_root_tree
		self.in_root_file = None
		self.in_root_tree = None

	def Reset_Braches_Lists_And_Dictionaries(self):
		self.branches1DTotal = BRANCHES1DTOTAL[:]
		self.branches1DType = BRANCHES1DTYPE.copy()
		self.branchesWavesTotal = BRANCHESWAVESTOTAL[:]
		self.branchesWavesType = BRANCHESWAVESTYPE.copy()
		self.branchesAll = self.branches1DTotal + self.branchesWavesTotal
		self.branches1DLoad = BRANCHES1DLOAD[:]
		self.branchesWavesLoad = BRANCHESWAVESLOAD[:]
		self.analysisScalarsBranches = ANALYSISSCALARBRANCHES[:]
		self.dicBraVect1D = OrderedDict()
		self.dicBraVectWaves = OrderedDict()
		self.hasBranch = {}

	def AddCFShiftCut(self):
		self.cut0CF += ro.TCut('CFShiftCut', 'abs(deltaTimeCF-{pp})<={ppc}'.format(pp=self.timeCFDeltasPeak, ppc=self.timeCFDeltasCut))

	def AddPeakPositionCut(self):
		self.cut0 += ro.TCut('peakTimeCut', 'abs(peakPosition-{pp})<={ppc}'.format(pp=self.peakTime, ppc=self.peakTimeCut))

	def AddPeakPositionCutCF(self):
		self.cut0CF += ro.TCut('peakTimeCutCF', 'abs(peakPositionCF-{pp})<={ppc}'.format(pp=self.peakTimeCF, ppc=self.peakTimeCut))

	def AddDiamondVoltageCut(self):
		if 'voltageDia' in self.branches1DTotal:
			if self.voltageDiaMaxSigmas != 0 and self.voltageDiaSpread != 0:
				self.cut0 += ro.TCut('voltageDiaSpreadCut', 'abs(voltageDia-{vm})<{s}*{vs}'.format(vm=self.voltageDiaMean, s=self.voltageDiaMaxSigmas, vs=self.voltageDiaSpread))
			self.cut0 += ro.TCut('voltageDiaOffsetCut', 'abs(voltageDia-{b})<{v}'.format(b=self.bias, v=self.voltageDiaMaxOffset))

	def FindPedestalPosition(self):
		print 'Calculating position of pedestals...', ;sys.stdout.flush()
		self.pedestalTimeIndices = [np.argwhere(np.bitwise_and(self.pedestalTEndPos - self.pedestalIntegrationTime <= timeVectEvi, timeVectEvi <= self.pedestalTEndPos)).flatten() for it, timeVectEvi in enumerate(self.timeVect)]
		self.pedestalTimeIndicesCF = [np.argwhere(np.bitwise_and(self.pedestalTEndPos - self.pedestalIntegrationTime <= np.add(timeVectEvi, self.time_CF_deltas[it]), np.add(timeVectEvi, self.time_CF_deltas[it]) <= self.pedestalTEndPos)).flatten() for it, timeVectEvi in enumerate(self.timeVect)]
		print 'Done'

	def FindSignalPositions(self, backward, forward):
		print 'Calculating position of signals...', ;sys.stdout.flush()
		self.signalTimeIndices = [np.argwhere(abs(timeVectEvi - self.peak_positions[it] - (forward - backward)/2.0) <= (forward + backward)/2.0).flatten() for it, timeVectEvi in enumerate(self.timeVect)]
		self.signalTimeIndicesCF = [np.argwhere(abs(np.add(timeVectEvi, self.time_CF_deltas[it]) - self.peak_positionsCF[it] - (forward - backward)/2.0) <= (forward + backward)/2.0).flatten() for it, timeVectEvi in enumerate(self.timeVect)]
		print 'Done'

	def CalculatePedestalsAndSignals(self):
		print 'Calculating pedestals and signals...', ;sys.stdout.flush()
		self.pedVect = np.array([self.signalWaveVect[ev, pedTimeIndxs].mean() if pedTimeIndxs.size > 0 else -10 for ev, pedTimeIndxs in enumerate(self.pedestalTimeIndices)])
		self.pedVectCF = np.array([self.signalWaveVect[ev, pedTimeIndxs].mean() if pedTimeIndxs.size > 0 else -10 for ev, pedTimeIndxs in enumerate(self.pedestalTimeIndicesCF)])
		self.pedSigmaVect = np.array([self.signalWaveVect[ev, pedTimeIndxs].std() if pedTimeIndxs.size > 1 else -10 for ev, pedTimeIndxs in enumerate(self.pedestalTimeIndices)])
		self.pedSigmaVectCF = np.array([self.signalWaveVect[ev, pedTimeIndxs].std() if pedTimeIndxs.size > 1 else -10 for ev, pedTimeIndxs in enumerate(self.pedestalTimeIndicesCF)])
		self.sigAndPedVect = np.array([self.signalWaveVect[ev, sigTimeIndxs].mean() if sigTimeIndxs.size > 0 else -10 for ev, sigTimeIndxs in enumerate(self.signalTimeIndices)])
		self.sigAndPedVectCF = np.array([self.signalWaveVect[ev, sigTimeIndxs].mean() if sigTimeIndxs.size > 0 else -10 for ev, sigTimeIndxs in enumerate(self.signalTimeIndicesCF)])
		self.sigAndPedSigmaVect = np.array([self.signalWaveVect[ev, sigTimeIndxs].std() if sigTimeIndxs.size > 1 else -10 for ev, sigTimeIndxs in enumerate(self.signalTimeIndices)])
		self.sigAndPedSigmaVectCF = np.array([self.signalWaveVect[ev, sigTimeIndxs].std() if sigTimeIndxs.size > 1 else -10 for ev, sigTimeIndxs in enumerate(self.signalTimeIndicesCF)])
		self.sigVect = np.subtract(self.sigAndPedVect, self.pedVect)
		self.sigVectCF = np.subtract(self.sigAndPedVectCF, self.pedVectCF)
		print 'Done'

	def FillPedestalsAndSignals(self):
		print 'Filling tree with scalars...'
		pedBra = self.analysisTree.Branch('pedestal', self.ped, 'pedestal/F')
		pedSigmaBra = self.analysisTree.Branch('pedestalSigma', self.pedSigma, 'pedestalSigma/F')
		pedSignalBra = self.analysisTree.Branch('signalAndPedestal', self.sigAndPed, 'signalAndPedestal/F')
		pedSignalSigmaBra = self.analysisTree.Branch('signalAndPedestalSigma', self.sigAndPed, 'signalAndPedestalSigma/F')
		sigBra = self.analysisTree.Branch('signal', self.sig, 'signal/F')
		entries = self.in_root_tree.GetEntries()
		self.CloseInputROOTFiles()
		self.analysisFile.cd()
		self.utils.CreateProgressBar(self.max_events)
		self.utils.bar.start()
		for ev in xrange(entries):
			# self.in_root_tree.GetEntry(ev)
			self.analysisTree.GetEntry(ev)
			if ev in self.eventVect:
				try:
					argum = np.argwhere(ev == self.eventVect).flatten()
					self.ped.itemset(self.pedVect[argum])
					self.pedSigma.itemset(self.pedSigmaVect[argum])
					self.sigAndPed.itemset(self.sigAndPedVect[argum])
					self.sigAndPedSigma.itemset(self.sigAndPedSigmaVect[argum])
					self.sig.itemset(self.sigVect[argum])
				except ValueError:
					ExitMessage('Could not fill event {ev}. Exiting...'.format(ev=ev), os.EX_DATAERR)
			else:
				self.ped.itemset(0)
				self.pedSigma.itemset(0)
				self.sigAndPed.itemset(0)
				self.sigAndPedSigma.itemset(0)
				self.sig.itemset(0)
			self.analysisTree.Fill()
			# pedBra.Fill()
			# pedSigmaBra.Fill()
			# pedSignalBra.Fill()
			# pedSignalSigmaBra.Fill()
			# sigBra.Fill()
			self.utils.bar.update(ev + 1)
		if not self.analysisTree.GetFriend(self.in_tree_name):
			self.analysisTree.AddFriend(self.in_tree_name, '{d}/{f}.root'.format(d=os.path.abspath(self.inDir), f=self.in_tree_name))
		# self.analysisTree.Write('', ro.TObject.kOverwrite)
		self.analysisTree.Write('', ro.TObject.kWriteDelete)
		self.utils.bar.finish()

	def FillAnalysisBranches(self):
		print 'Filling analysis tree ...'
		peakPosBra = self.analysisTree.Branch('peakPosition', self.peak_position, 'peakPosition/F')
		peakPosBra = self.analysisTree.Branch('peakPositionCF', self.peak_positionCF, 'peakPositionCF/F')
		pedBra = self.analysisTree.Branch('pedestal', self.ped, 'pedestal/F')
		pedBra = self.analysisTree.Branch('pedestalCF', self.pedCF, 'pedestalCF/F')
		pedSigmaBra = self.analysisTree.Branch('pedestalSigma', self.pedSigma, 'pedestalSigma/F')
		pedSigmaBra = self.analysisTree.Branch('pedestalSigmaCF', self.pedSigmaCF, 'pedestalSigmaCF/F')
		pedSignalBra = self.analysisTree.Branch('signalAndPedestal', self.sigAndPed, 'signalAndPedestal/F')
		pedSignalBra = self.analysisTree.Branch('signalAndPedestalCF', self.sigAndPedCF, 'signalAndPedestalCF/F')
		pedSignalSigmaBra = self.analysisTree.Branch('signalAndPedestalSigma', self.sigAndPedSigma, 'signalAndPedestalSigma/F')
		pedSignalSigmaBra = self.analysisTree.Branch('signalAndPedestalSigmaCF', self.sigAndPedSigmaCF, 'signalAndPedestalSigmaCF/F')
		sigBra = self.analysisTree.Branch('signal', self.sig, 'signal/F')
		sigBra = self.analysisTree.Branch('signalCF', self.sigCF, 'signalCF/F')
		entries = self.in_root_tree.GetEntries() - self.loaded_entries
		entries = entries if entries < self.load_entries else self.load_entries
		self.CloseInputROOTFiles()
		self.analysisFile.cd()
		self.utils.CreateProgressBar(entries)
		self.utils.bar.start()
		for ev in xrange(self.start_entry, self.start_entry + entries):
			# self.in_root_tree.GetEntry(ev)
			if ev in self.eventVect:
				try:
					argum = np.argwhere(ev == self.eventVect).flatten()
					self.peak_position.itemset(self.peak_positions[argum])
					self.peak_positionCF.itemset(self.peak_positionsCF[argum])
					self.ped.itemset(self.pedVect[argum])
					self.pedCF.itemset(self.pedVectCF[argum])
					self.pedSigma.itemset(self.pedSigmaVect[argum])
					self.pedSigmaCF.itemset(self.pedSigmaVectCF[argum])
					self.sigAndPed.itemset(self.sigAndPedVect[argum])
					self.sigAndPedCF.itemset(self.sigAndPedVectCF[argum])
					self.sigAndPedSigma.itemset(self.sigAndPedSigmaVect[argum])
					self.sigAndPedSigmaCF.itemset(self.sigAndPedSigmaVectCF[argum])
					self.sig.itemset(self.sigVect[argum])
					self.sigCF.itemset(self.sigVectCF[argum])
				except ValueError:
					ExitMessage('Could not fill event {ev}. Exiting...'.format(ev=ev), os.EX_DATAERR)
			else:
				self.peak_position.itemset(-100000)
				self.peak_positionCF.itemset(-100000)
				self.ped.itemset(-100000)
				self.pedCF.itemset(-100000)
				self.pedSigma.itemset(-100000)
				self.pedSigmaCF.itemset(-100000)
				self.sigAndPed.itemset(-100000)
				self.sigAndPedCF.itemset(-100000)
				self.sigAndPedSigma.itemset(-100000)
				self.sigAndPedSigmaCF.itemset(-100000)
				self.sig.itemset(-100000)
				self.sigCF.itemset(-100000)
			self.analysisTree.Fill()
			self.utils.bar.update(ev - self.start_entry + 1)
		self.loaded_entries += entries
		self.start_entry = self.loaded_entries
		if not self.analysisTree.GetFriend(self.in_tree_name):
			self.analysisTree.AddFriend(self.in_tree_name, '{d}/{f}.root'.format(d=os.path.abspath(self.inDir), f=self.in_tree_name))
		# self.analysisTree.Write('', ro.TObject.kOverwrite)
		self.analysisTree.Write()
		self.utils.bar.finish()

	def ExtractMeanOfWaveforms(self):
		self.signalWaveMeanVect = self.signalWaveVect.mean(axis=0)
		self.signalWaveSigmaVect = self.signalWaveVect.std(axis=0)

	def DrawHisto(self, name, xmin, xmax, deltax, var, varname, cuts='', option='e'):
		if not IsFloat(xmin) or not IsFloat(xmax) or not IsFloat(deltax):
			print 'Won\'t create histogram as the limits are not well defined (xmin, xmax, deltax): {mi}, {ma}, {dx}'.format(mi=xmin, ma=xmax, dx=deltax)
			return
		elif deltax <= 0 or xmin >= xmax:
			print 'Won\'t create histogram as the limits are not well defined (xmin, xmax, deltax): {mi}, {ma}, {dx}'.format(mi=xmin, ma=xmax, dx=deltax)
			return
		ro.TFormula.SetMaxima(100000)
		if self.histo.has_key(name):
			if self.histo[name]:
				self.histo[name].Delete()
			del self.histo[name]
		self.histo[name] = ro.TH1F('h_' + name, 'h_' + name, RoundInt(float(xmax - xmin) / deltax), xmin, xmax)
		self.histo[name].GetXaxis().SetTitle(varname)
		self.histo[name].GetYaxis().SetTitle('entries')
		if 'goff' not in option:
			if self.canvas.has_key(name):
				if self.canvas[name]:
					self.canvas[name].Close()
				del self.canvas[name]
			self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
			self.canvas[name].cd()
		cuts0 = self.cut0.GetTitle() if cuts == '' else cuts
		self.analysisTree.Draw('{v}>>h_{n}'.format(v=var, n=name), cuts0, option)
		if 'goff' not in option:
			self.canvas[name].SetGridx()
			self.canvas[name].SetGridy()
			self.canvas[name].SetTicky()
			ro.gPad.Update()
			SetDefault1DStats(self.histo[name])
			ro.gPad.Update()
		ro.TFormula.SetMaxima(1000)

	def DrawProfile(self, name, varx, xmin, xmax, deltax, xname, vary, ymin, ymax, yname, cuts='', options='e hist'):
		if not IsFloat(xmin) or not IsFloat(xmax) or not IsFloat(deltax) or not IsFloat(ymin) or not IsFloat(ymax):
			print 'Won\'t create profile as the limits are not well defined (xmin, xmax, deltax, ymin, ymax): {mi}, {ma}, {dx}, {ym}, {yma}'.format(mi=xmin, ma=xmax, dx=deltax, ym=ymin, yma=ymax)
			return
		elif deltax <= 0 or xmin >= xmax or ymin >= ymax:
			print 'Won\'t create profile as the limits are not well defined (xmin, xmax, deltax, ymin, ymax): {mi}, {ma}, {dx}, {ym}, {yma}'.format(mi=xmin, ma=xmax, dx=deltax, ym=ymin, yma=ymax)
			return
		ro.TFormula.SetMaxima(100000)
		if self.profile.has_key(name):
			if self.profile[name]:
				self.profile[name].Delete()
			del self.profile[name]
		self.profile[name] = ro.TProfile('h_' + name, 'h_' + name, int(RoundInt(float(xmax - xmin) / deltax) + 2), xmin - deltax, xmax + deltax, ymin, ymax)
		self.profile[name].GetXaxis().SetTitle(xname)
		self.profile[name].GetYaxis().SetTitle(yname)
		if 'goff' not in options:
			if self.canvas.has_key(name):
				if self.canvas[name]:
					self.canvas[name].Close()
				del self.canvas[name]
			self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
			self.canvas[name].cd()
		cuts0 = cuts
		self.analysisTree.Draw('{vy}:{vx}>>h_{n}'.format(vy=vary, vx=varx, n=name), cuts0, options)
		if 'goff' not in options:
			self.canvas[name].SetGridx()
			self.canvas[name].SetGridy()
			self.canvas[name].SetTicky()
			ro.gPad.Update()
			# SetDefault1DStats(self.profile[name])
			ro.gPad.Update()
		ro.TFormula.SetMaxima(1000)

	def DrawHisto2D(self, name, varx, xmin, xmax, deltax, xname, vary, ymin, ymax, deltay, yname, cuts='', option='colz', num_evts=1000000000, start_ev=0):
		if not IsFloat(xmin) or not IsFloat(xmax) or not IsFloat(deltax) or not IsFloat(ymin) or not IsFloat(ymax) or not IsFloat(deltay):
			print 'Won\'t create histogram as the limits are not well defined (xmin, xmax, deltax, ymin, ymax, deltay): {mi}, {ma}, {dx}, {ym}, {yma}, {dy}'.format(mi=xmin, ma=xmax, dx=deltax, ym=ymin, yma=ymax, dy=deltay)
			return
		elif deltax <= 0 or xmin >= xmax or deltay <= 0 or ymin >= ymax:
			print 'Won\'t create histogram as the limits are not well defined (xmin, xmax, deltax, ymin, ymax, deltay): {mi}, {ma}, {dx}, {ym}, {yma}, {dy}'.format(mi=xmin, ma=xmax, dx=deltax, ym=ymin, yma=ymax, dy=deltay)
			return
		ro.TFormula.SetMaxima(100000)
		if self.histo.has_key(name):
			if self.histo[name]:
				self.histo[name].Delete()
			del self.histo[name]
		self.histo[name] = ro.TH2F('h_' + name, 'h_' + name, int(RoundInt(float(xmax - xmin) / deltax) + 2), xmin - deltax, xmax + deltax, int(RoundInt(float(ymax - ymin) / deltay) + 2), ymin - deltay, ymax + deltay)
		self.histo[name].GetXaxis().SetTitle(xname)
		self.histo[name].GetYaxis().SetTitle(yname)
		self.histo[name].GetZaxis().SetTitle('entries')
		if 'goff' not in option:
			if self.canvas.has_key(name):
				if self.canvas[name]:
					self.canvas[name].Close()
				del self.canvas[name]
			self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
			self.canvas[name].cd()
		cuts0 = self.cut0.GetTitle() if cuts == '' else cuts
		self.analysisTree.Draw('{y}:{x}>>h_{n}'.format(y=vary, x=varx, n=name), cuts0, option, num_evts, start_ev)
		if 'goff' not in option:
			self.canvas[name].SetGridx()
			self.canvas[name].SetGridy()
			self.canvas[name].SetTickx()
			self.canvas[name].SetTicky()
			ro.gPad.Update()
			SetDefault2DStats(self.histo[name])
		ro.TFormula.SetMaxima(1000)

	def DrawGraphFromTree(self, name, varx, xname, vary, yname, cuts='', option='AP', num_evts=100, start_ev=0):
		ro.TFormula.SetMaxima(100000)
		if self.graph.has_key(name):
			if self.graph[name]:
				self.graph[name].Delete()
			del self.graph[name]
		option = option.replace('goff', '')  # does not work if goff is present. Root has to draw the graph!
		CreateCanvasInDic(self.canvas, name)
		self.canvas[name].cd()
		self.analysisTree.Draw('{y}:{x}'.format(y=vary, x=varx), cuts, option, num_evts, start_ev)
		self.graph[name] = ro.gPad.GetPrimitive('Graph')
		self.graph[name].SetTitle(name)
		self.graph[name].GetXaxis().SetTitle(xname)
		self.graph[name].GetYaxis().SetTitle(yname)

	def PlotCFDeltaDistribution(self, name='CFDeltaDist', low_t0=-1, up_t0=1, nbins=100, cut=''):
		nbins2 = nbins
		funcArgs = (name, low_t0, up_t0, float(up_t0 - low_t0) / nbins2, 'deltaTimeCF*1000000', 'CF Time Shift [us]', cut, 'e')
		self.DrawHisto(*funcArgs)
		if not name in self.histo.keys():
			print 'There was a problem creating the histogram {n}'.format(n=name)
			return
		if not self.histo[name] or IsHistogramEmpty(self.histo[name]):
			print 'There was a problem with the created histogram {n}'.format(n=name)
			return
		tempbin0, tempbinf = self.histo[name].FindFirstBinAbove(), self.histo[name].FindLastBinAbove()
		self.histo[name].SetBinContent(tempbin0, 0)
		self.histo[name].SetBinContent(tempbinf, 0)
		self.histo[name].GetXaxis().SetRangeUser(-0.5, 0.5)
		peakbin = self.histo[name].GetMaximumBin()
		self.timeCFDeltasPeak = self.histo[name].GetBinCenter(peakbin) * 1e-6
		self.histo[name].GetXaxis().SetRangeUser(low_t0, up_t0)
		deltax = CheckBinningForFit(self, name, self.DrawHisto, funcArgs, 3, 12) if not self.settings.is_cal_run else (up_t0 - low_t0)/100.
		tempbin0, tempbinf = self.histo[name].FindFirstBinAbove(), self.histo[name].FindLastBinAbove()
		self.histo[name].SetBinContent(tempbin0, 0)
		self.histo[name].SetBinContent(tempbinf, 0)
		self.histo[name].GetXaxis().SetRangeUser(max(low_t0, self.histo[name].GetMean() - 5 * self.histo[name].GetRMS()), min(up_t0, self.histo[name].GetMean() + 5 * self.histo[name].GetRMS()))
		fb, lb = self.histo[name].FindFirstBinAbove(0), self.histo[name].FindLastBinAbove(0)
		fb = fb + 1 if self.histo[name].GetBinContent(fb) > 3 * self.histo[name].GetBinContent(fb + 1) else fb
		lb = lb - 1 if self.histo[name].GetBinContent(lb) > 3 * self.histo[name].GetBinContent(lb - 1) else lb
		low_t, up_t = max(low_t0, self.histo[name].GetBinLowEdge(fb)), min(up_t0, self.histo[name].GetBinLowEdge(lb + 1))
		if self.histo[name].Integral(fb, lb) > 9:
			func = ro.TF1('fit_' + name, 'gaus(0)+gaus(3)', low_t, up_t)
			func.SetNpx(1000)
			params = np.array([self.histo[name].GetBinContent(peakbin), self.timeCFDeltasPeak, self.histo[name].GetRMS() / 2., self.histo[name].GetBinContent(peakbin) / 10., self.histo[name].GetMean(), self.histo[name].GetRMS() * 2], 'f8')
			func.SetParameters(params)
			func.SetParLimits(0, 1, 10 * self.histo[name].GetBinContent(peakbin))
			func.SetParLimits(1, self.timeCFDeltasPeak, up_t)
			func.SetParLimits(2, self.histo[name].GetRMS() / 10., 2 * self.histo[name].GetRMS())
			func.SetParLimits(3, 0, 2 * self.histo[name].GetBinContent(peakbin))
			func.SetParLimits(4, -10 * abs(low_t), 10 * abs(up_t))
			func.SetParLimits(5, 0.5 * self.histo[name].GetRMS(), 10 * self.histo[name].GetRMS())
			fit = self.histo[name].Fit('fit_' + name, 'QEMSB', '', low_t, up_t)
			xpeak = func.GetMaximumX(low_t, up_t)
			fit = self.histo[name].Fit('fit_' + name, 'QEMSB', '', max(low_t, xpeak - 2 * self.histo[name].GetRMS()), min(up_t, xpeak + 2 * self.histo[name].GetRMS()))
			# params = np.array([fit.Parameter(i) for i in xrange(6)], 'f8')
			SetDefaultFitStats(self.histo[name], func)
			xpeak = func.GetMaximumX(max(low_t, xpeak - 2 * self.histo[name].GetRMS()), min(up_t, xpeak + 2 * self.histo[name].GetRMS()))
			self.timeCFDeltasPeak = np.divide(xpeak, 1e6, dtype='f8')
		self.canvas[name].Modified()
		ro.gPad.Update()

	def PlotPeakPositionDistributions(self, name='peakPosDist', low_t=1, up_t=4, nbins=200, cut='', isCF=False):
		appendCF = '' if not isCF else 'CF'
		nbins2 = nbins
		funcArgs = (name, low_t, up_t, float(up_t - low_t) / nbins2, 'peakPosition{a}*1000000'.format(a=appendCF), 'Peak Position [us]', cut, 'e')
		self.DrawHisto(*funcArgs)
		if not name in self.histo.keys():
			print 'There was a problem creating the histogram {n}'.format(n=name)
			return
		if not self.histo[name] or IsHistogramEmpty(self.histo[name]):
			print 'There was a problem with the created histogram {n}'.format(n=name)
			return
		deltax = CheckBinningForFit(self, name, self.DrawHisto, funcArgs, 3, 12)
		self.histo[name].GetXaxis().SetRangeUser(max(low_t, self.histo[name].GetMean() - 5 * self.histo[name].GetRMS()), min(up_t, self.histo[name].GetMean() + 5 * self.histo[name].GetRMS()))
		self.peakTime = self.histo[name].GetBinCenter(self.histo[name].GetMaximumBin())
		# peakTimeMP = self.histo[name].GetBinCenter(self.histo[name].GetMaximumBin())
		# self.peakTime = self.histo[name].GetMean()
		if self.histo[name].Integral() > 9:
			func = ro.TF1('fit_' + name, 'gaus(0)+gaus(3)', low_t, up_t)
			func.SetNpx(1000)
			skew = self.histo[name].GetSkewness()
			const_p = [self.histo[name].GetMaximum() * i for i in [1.5, 0.5]]
			mean_p_i = [-0.5, 0.5] if skew >= 0 else [0.5, -0.5]
			mean_p = [self.histo[name].GetMean() + i * self.histo[name].GetRMS() for i in mean_p_i]
			rms_p_i = [0.5, 1.5] if skew >= 0 else [1.5, 0.5]
			rms_p = [self.histo[name].GetRMS() * i for i in rms_p_i]
			params = np.array((const_p[0], mean_p[0], rms_p[0], const_p[1], mean_p[1], rms_p[1]), 'float64')
			func.SetParameters(params)
			func.SetParLimits(0, 0, self.histo[name].GetMaximum() * 2)
			func.SetParLimits(2, max(0.002, 0.5 * self.histo[name].GetRMS()) if self.is_cal_run else 0, self.histo[name].GetRMS())
			func.SetParLimits(3, 0, self.histo[name].GetMaximum() * 2)
			func.SetParLimits(5, max(0.002, self.histo[name].GetRMS()), self.histo[name].GetRMS() * 10)
			if skew >= 0:
				func.SetParLimits(1, self.histo[name].GetMean() - 5 * self.histo[name].GetRMS(), self.histo[name].GetMean() + self.histo[name].GetRMS())
				func.SetParLimits(4, self.histo[name].GetMean() - self.histo[name].GetRMS(), self.histo[name].GetMean() + 5 * self.histo[name].GetRMS())
			else:
				func.SetParLimits(1, self.histo[name].GetMean() - self.histo[name].GetRMS(), self.histo[name].GetMean() + 5 * self.histo[name].GetRMS())
				func.SetParLimits(4, self.histo[name].GetMean() - 5 * self.histo[name].GetRMS(), self.histo[name].GetMean() + self.histo[name].GetRMS())
			fit = self.histo[name].Fit('fit_' + name, 'QEMSB', '', max(low_t, self.histo[name].GetMean() - 4 * self.histo[name].GetRMS()), min(up_t, self.histo[name].GetMean() + 4 * self.histo[name].GetRMS()))
			params = np.array([fit.Parameter(i) for i in xrange(6)], 'f8')
			xpeak = func.GetMaximumX(max(low_t, self.histo[name].GetMean() - 4 * self.histo[name].GetRMS()), min(up_t, self.histo[name].GetMean() + 4 * self.histo[name].GetRMS()))
			fit2 = self.histo[name].Fit('fit_' + name, 'QEMSB', '', xpeak - 3 * self.histo[name].GetRMS(), xpeak + 3* self.histo[name].GetRMS())
			params = np.array([fit2.Parameter(i) for i in xrange(6)], 'f8')
			# fit = self.histo[name].Fit('fit_' + name, 'QEMSB', '', params[1] - 2 * params[2], params[1] + 2 * params[2])
			SetDefaultFitStats(self.histo[name], func)
			#TODO check if fit is good enough to update peakTime
			if not self.settings.is_cal_run:
				xpeak = func.GetMaximumX(max(low_t, xpeak - 3 * self.histo[name].GetRMS()), min(up_t, xpeak + 3 * self.histo[name].GetRMS()))
			else:
				xpeak = self.histo[name].GetMean()
			if not isCF:
				self.peakTime = np.divide(xpeak, 1e6, dtype='f8')
			else:
				self.peakTimeCF = np.divide(xpeak, 1e6, dtype='f8')
		self.canvas[name].Modified()
		ro.gPad.Update()

	def PlotHVDiaDistribution(self, name='DUTVoltageHisto', low_v=0, up_v=0, nbins=500, cut=''):
		resVoltage = 0.025
		minv = low_v if low_v != 0 and up_v != 0 else GetMinimumBranch(self.analysisTree, 'voltageDia', 'voltageDia!=0')
		maxv = up_v if low_v != 0 and up_v != 0 else GetMaximumBranch(self.analysisTree, 'voltageDia', 'voltageDia!=0')
		deltav = max(float(maxv - minv) / nbins, resVoltage)
		funcArgs = (name, minv, maxv + resVoltage, resVoltage, 'voltageDia', 'Voltage on Diamond [V]', cut, 'e')
		self.DrawHisto(*funcArgs)
		if not name in self.histo.keys():
			print 'There was a problem creating the histogram' + name
			return
		if not self.histo[name] or IsHistogramEmpty(self.histo[name]):
			print 'There was a problem with the created histogram' + name
			return
		deltax = CheckBinningForFit(self, name, self.DrawHisto, funcArgs, 3, 6, resVoltage, True)
		self.voltageDiaMean = self.histo[name].GetMean()
		self.voltageDiaSpread = self.histo[name].GetRMS()

	def PlotSignal(self, name='signal', bins=0, cuts='', option='e', minx=-50, maxx=950, branch='signal', optimizeBinning=True):
		if not self.hasBranch['signal']:
			return
		if 'vcal' in branch.lower():
			if self.bias >= 0:
				plotvar = '1000*' + branch if not self.is_cal_run or self.cal_run_type == 'out' else '-1000*' + branch
				# vmax, vmin, deltav = -self.analysisTree.GetMinimum('signal'), -self.analysisTree.GetMaximum('signal'), self.signal_ch.adc_to_volts_cal['p1'] * 100
				plotVarName = branch + ' [mV]' if not self.is_cal_run or self.cal_run_type == 'out' else '-' + branch + ' [mV]'
			else:
				plotvar = '-1000*' + branch if not self.is_cal_run or self.cal_run_type == 'out' else '1000*' + branch
				# vmax, vmin, deltav = -self.analysisTree.GetMinimum('signal'), -self.analysisTree.GetMaximum('signal'), self.signal_ch.adc_to_volts_cal['p1'] * 100
				plotVarName = '-' + branch + ' [mV]' if not self.is_cal_run or self.cal_run_type == 'out' else branch + ' [mV]'

		elif 'charge' in branch.lower():
			if self.bias >= 0:
				plotvar = branch if not self.is_cal_run or self.cal_run_type == 'out' else '-' + branch
				# vmax, vmin, deltav = -self.analysisTree.GetMinimum('signal'), -self.analysisTree.GetMaximum('signal'), self.signal_ch.adc_to_volts_cal['p1'] * 100
				plotVarName = branch + ' [e]' if not self.is_cal_run or self.cal_run_type == 'out' else '-' + branch + ' [e]'
			else:
				plotvar = '-' + branch if not self.is_cal_run or self.cal_run_type == 'out' else branch
				# vmax, vmin, deltav = -self.analysisTree.GetMinimum('signal'), -self.analysisTree.GetMaximum('signal'), self.signal_ch.adc_to_volts_cal['p1'] * 100
				plotVarName = '-' + branch + ' [e]' if not self.is_cal_run or self.cal_run_type == 'out' else branch + ' [e]'

		else:
			if self.bias >= 0:
				plotvar = '-1000*' + branch if not self.is_cal_run or self.cal_run_type == 'out' else '1000*' + branch
				# vmax, vmin, deltav = -self.analysisTree.GetMinimum('signal'), -self.analysisTree.GetMaximum('signal'), self.signal_ch.adc_to_volts_cal['p1'] * 100
				plotVarName = '-' + branch + ' [mV]' if not self.is_cal_run or self.cal_run_type == 'out' else branch + ' [mV]'
			else:
				plotvar = '1000*' + branch if not self.is_cal_run or self.cal_run_type == 'out' else '-1000*' + branch
				# vmin, vmax, deltav = self.analysisTree.GetMinimum('signal'), self.analysisTree.GetMaximum('signal'), self.signal_ch.adc_to_volts_cal['p1'] * 100
				plotVarName = branch + ' [mV]' if not self.is_cal_run or self.cal_run_type == 'out' else '-' + branch + ' [mV]'
		(vmin, vmax) = (minx, maxx) if 'charge' not in branch.lower() else (TruncateFloat(self.vcal_to_q.Q_in_e_from_mV(minx).nominal_value, 2000.), TruncateFloat(self.vcal_to_q.Q_in_e_from_mV(maxx).nominal_value, 2000.))
		deltav = self.delta_v_signal
		deltav = deltav if bins == 0 else (vmax - vmin) / float(bins)
		deltav = deltav if 'charge' not in branch.lower() else TruncateFloat(deltav if bins!=0 else deltav / self.vcal_to_q.vgain.nominal_value, 10.)
		deltav = deltav if deltav != 0 else 10.
		self.DrawHisto(name, vmin - deltav/2.0, vmax + deltav/2.0, deltav, plotvar, plotVarName, cuts, option)
		if optimizeBinning:
			funcArgs = (name, vmin - deltav/2.0, vmax + deltav/2.0, deltav, plotvar, plotVarName, cuts, option)
			truncateResol = 0.1 if 'charge' not in branch.lower() else 10.
			temp = CheckBinningForFit(self, name, self.DrawHisto, funcArgs, 3, 7, truncateResol)
			# self.delta_v_signal = deltav if temp == 0 else temp

	def PlotPedestal(self, name='pedestal', bins=0, cuts='', option='e', minx=0, maxx=0, branch='pedestal'):
		if not self.hasBranch['pedestal']:
			return
		vmin = minx if minx != 0 else GetMinimumBranch(self.analysisTree, branch, cuts) * 1000 if 'charge' not in branch.lower() else GetMinimumBranch(self.analysisTree, branch, cuts)
		vmax = maxx if maxx != 0 else GetMaximumBranch(self.analysisTree, branch, cuts) * 1000 if 'charge' not in branch.lower() else GetMaximumBranch(self.analysisTree, branch, cuts)
		deltax = (vmax - vmin) / 100.0 if bins == 0 else (vmax - vmin) / float(bins)
		self.DrawHisto(name, vmin - deltax/2.0, vmax + deltax/2.0, deltax, '{f}*{b}'.format(f=1000 if 'charge' not in branch.lower() else 1, b=branch), '{b} {u}'.format(b=branch, u='[e]' if 'charge' in branch.lower() else '[mV]'), cuts, option)
		if not name in self.histo.keys():
			print 'The dictionary histogram does not have the key ' + name
			return
		if not self.histo[name] or IsHistogramEmpty(self.histo[name]):
			print 'The created histogram {n} is empty'.format(n=name)
			return
		if self.histo[name].GetMean() == 0 and self.histo[name].GetRMS() == 0:
			print 'The created histogram {n} is empty'.format(n=name)
			return
		lowbin, highbin = self.histo[name].FindFirstBinAbove(0), self.histo[name].FindLastBinAbove(0)
		vmin, vmax = self.histo[name].GetBinLowEdge(lowbin), self.histo[name].GetBinLowEdge(highbin + 1)
		funcArgs = (name, vmin, vmax, deltax, '{f}*{b}'.format(f=1000 if 'charge' not in branch.lower() else 1, b=branch), '{b} {u}'.format(b=branch, u='[e]' if 'charge' in branch.lower() else '[mV]'), cuts, option)
		truncateResol = 0.1 if 'charge' not in branch.lower() else 10.
		temp = CheckBinningForFit(self, name, self.DrawHisto, funcArgs, 3, 6, truncateResol)
		deltax = deltax if temp == 0 else temp
		func = ro.TF1('fit_' + name, 'gaus', vmin, vmax)
		func.SetNpx(1000)
		mean_p, sigma_p = self.histo[name].GetMean(), self.histo[name].GetRMS()
		vmin = min(vmin, mean_p - 4 * sigma_p)
		vmax = max(vmax, mean_p + 4 * sigma_p)
		width = max(abs(vmin), abs(vmax))
		vmin, vmax = -width, width
		fit = self.histo[name].Fit('fit_' + name, 'QEMS', '', mean_p - 2 * sigma_p, mean_p + 2 * sigma_p)
		params = np.array((fit.Parameter(0), fit.Parameter(1), fit.Parameter(2)), 'float64')
		func.SetParameters(params)
		self.DrawHisto(name, vmin - deltax / 2.0, vmax + deltax / 2.0, deltax, '{f}*{b}'.format(f=1000 if 'charge' not in branch.lower() else 1, b=branch), '{b} {u}'.format(b=branch, u='[e]' if 'charge' in branch.lower() else '[mV]'), cuts, option)
		fit = self.histo[name].Fit('fit_' + name, 'QEMS', '', params[1] - 2 * params[2], params[1] + 2 * params[2])
		SetDefaultFitStats(self.histo[name], func)
		# self.pedestal_sigma = func.GetParameter(2) if not self.is_cal_run or self.cal_run_type == 'out' else self.histo[name].GetRMS()
		if 'charge' in branch.lower():
			if not 'cf' in branch.lower():
				self.pedestal_charge_sigma = self.histo[name].GetRMS()
			else:
				self.pedestal_charge_sigmaCF = self.histo[name].GetRMS()
		elif 'vcal' in branch.lower():
			if not 'cf' in branch.lower():
				self.pedestal_vcal_sigma = self.histo[name].GetRMS()
			else:
				self.pedestal_vcal_sigmaCF = self.histo[name].GetRMS()
		else:
			if not 'cf' in branch.lower():
				self.pedestal_sigma = self.histo[name].GetRMS()
				self.pedestal_sigma_err = self.histo[name].GetRMSError()
			else:
				self.pedestal_sigmaCF = self.histo[name].GetRMS()
				self.pedestal_sigma_errCF = self.histo[name].GetRMSError()

	def PlotWaveforms(self, name='SignalWaveform', sigtype='signal', vbins=0, cuts='', option='colz', start_ev=0, num_evs=0, do_logz=False, do_cf=False):
		if do_cf:
			if not self.hasBranch['deltaTimeCF']:
				return
		if ('signal' in sigtype.lower() and not self.hasBranch['voltageSignal']) or ('trig' in sigtype.lower() and not self.hasBranch['voltageTrigger']) or ('veto' in sigtype.lower() and not self.hasBranch['voltageVeto']):
			return
		if not do_cf:
			var = '1000*(voltageSignal-pedestal)' if 'signal_ped_cor' in sigtype.lower() else '1000*voltageSignal' if 'signal' in sigtype.lower() else '1000*voltageTrigger' if 'trig' in sigtype.lower() else '1000*voltageVeto' if 'veto' in sigtype.lower() else ''
		else:
			var = '1000*(voltageSignal-pedestalCF)' if 'signal_ped_cor' in sigtype.lower() else '1000*voltageSignal' if 'signal' in sigtype.lower() else '1000*voltageTrigger' if 'trig' in sigtype.lower() else '1000*voltageVeto' if 'veto' in sigtype.lower() else ''
		if var == '':
			print 'sigtype should be "signal", "signal_ped_corrected", "trigger" or "veto"'
			return
		# cuts0 = self.cut0.GetTitle() if cuts == '' else cuts
		vname = ('signal - ped' if var == '1000*(voltageSignal-pedestal)' else 'signal' if var == '1000*voltageSignal' else 'trigger' if var == '1000*voltageTrigger' else 'veto' if var == '1000*voltageVeto' else '') + ' [mV]'
		num_events = self.analysisTree.GetEntries() if num_evs == 0 else num_evs
		# tmin, tmax, deltat = 1e6*self.analysisTree.GetMinimum('time'), 1e6*self.analysisTree.GetMaximum('time'), 1e6*self.settings.time_res
		tmin, tmax, deltat = 1e6*GetMinimumBranch(self.analysisTree, 'time', cuts), 1e6*GetMaximumBranch(self.analysisTree, 'time', cuts), 1e6*self.settings.time_res
		varxname = '1000000*time' if not do_cf else '1000000*(time+deltaTimeCF)'
		varxAxisName = 'time [us]' if not do_cf else 'CF time [us]'
		if do_cf:
			tmin += 1e6 * GetMinimumBranch(self.analysisTree, 'deltaTimeCF', cuts)
			tmax += 1e6 * GetMaximumBranch(self.analysisTree, 'deltaTimeCF', cuts)
		vmin, vmax, deltav = -1000, 1000, (1000*self.signal_ch.adc_to_volts_cal['p1'] * 10 if 'voltageSignal' in var else 1000*self.settings.sigRes * 10)
		miny, maxy = vmin, vmax
		if 'signal' in var.lower():
			extra_y = self.pedestal_sigma if self.pedestal_sigma != 0 else 10
			extra_y = extra_y if not self.is_cal_run or self.cal_run_type == 'out' else 10
			miny, maxy = self.analysisTree.GetMinimum('voltageSignal') * 1000 - 3 * extra_y, self.analysisTree.GetMaximum('voltageSignal') * 1000 + 3 * extra_y
			vmin, vmax = TruncateFloat(miny, deltav) - deltav / 2.0, TruncateFloat(maxy, deltav) + deltav / 2.0
		if vbins == 0:
			self.DrawHisto2D(name, varxname, tmin - deltat/2.0, tmax + deltat/2.0, deltat, varxAxisName, var, vmin, vmax, deltav, vname, cuts, option, num_events, start_ev)
		else:
			self.DrawHisto2D(name, varxname, tmin - deltat/2.0, tmax + deltat/2.0, deltat, varxAxisName, var, vmin, vmax, (vmax-vmin)/vbins, vname, cuts, option, num_events, start_ev)
		if do_logz and 'goff' not in option and name in self.histo.keys():
			self.canvas[name].SetLogz()
		if name in self.histo.keys(): self.histo[name].GetYaxis().SetRangeUser(miny, maxy)

	def PeakPositionStudy(self, xmin, xmax, deltax, pos_or_cut='pos', do_fit=True):
		if not (pos_or_cut.lower().startswith('pos') or pos_or_cut.lower().startswith('cut')):
			return 'pos_or_cut must be "pos" or "cut" to specify the type of analysis you want to do. Try again'
		t0_vec = np.append(np.arange(xmin, xmax, deltax, dtype='float64'), xmax)
		store_peakTimeCut = self.peakTimeCut
		store_peakTime = 0
		if self.peakTime:
			store_peakTime = self.peakTime
		else:
			print 'Find first the peak of the distribution with PlotPeakDistributions'
			return
		self.utils.CreateProgressBar(len(t0_vec))
		self.utils.bar.start()
		for index, t0 in np.ndenumerate(t0_vec):
			self.ResetCut0()
			self.CreateCut0()
			if pos_or_cut.lower().startswith('pos'):
				self.peakTime = t0
			else:
				self.peakTimeCut = t0
			self.AddPeakPositionCut()
			peakT, timeC = self.peakTime, self.peakTimeCut
			self.PlotWaveforms('SigWFs_peak_{p:.3f}us_cut_{c:.3f}ns'.format(p=peakT * 1e6, c=timeC * 1e9), 'signal', cuts=self.cut0.GetTitle(), do_logz=True)
			self.PlotSignal('PH_peak_{p:.3f}us_cut_{c:.3f}ns'.format(p=peakT * 1e6, c=timeC * 1e9), cuts=self.cut0.GetTitle())
			if do_fit:
				self.FitLanGaus('PH_peak_{p:.3f}us_cut_{c:.3f}ns'.format(p=peakT * 1e6, c=timeC * 1e9))
			self.utils.bar.update(index[0] + 1)
		self.utils.bar.finish()
		self.peakTimeCut = store_peakTimeCut
		self.peakTime = store_peakTime

	def ResetTreeToOriginal(self, keepBranches=['event','time','voltageSignal','voltageTrigger','voltageVeto','vetoedEvent','badShape','badPedestal','voltageHV','currentHV','timeHV']):
		print 'Restoring tree with the following branches:', keepBranches, '...'
		raw_input('Press a key and Enter to continue: ')
		self.OpenAnalysisROOTFile('READ')
		self.LoadAnalysisTree()
		self.in_root_tree.SetBranchStatus('*', 0)
		for branch in keepBranches:
			if self.TreeHasBranch(branch):
				self.in_root_tree.SetBranchStatus(branch, 1)
		newFile = ro.TFile('{o}/temp.root'.format(o=self.outDir), 'recreate')
		newTree = self.in_root_tree.CloneTree()
		newTree.Print()
		newFile.Write()
		del self.in_root_file
		del newFile
		self.in_root_file = None
		checkFile = ro.TFile('{o}/temp.root'.format(o=self.outDir), 'READ')
		checkTree = checkFile.Get(self.in_tree_name)
		doMoveFile = True
		if checkTree:
			for branch in keepBranches:
				if not checkTree.GetLeaf(branch):
					doMoveFile = False
					break
			if doMoveFile:
				print 'The file was cloned successfully :)'
				checkFile.Close()
				del checkFile
				shutil.move('{o}/temp.root'.format(o=self.outDir), '{o}/{f}'.format(o=self.outDir, f=self.inputFile))
				return
		print 'The file was not cloned successfully :S. Check original tree and "temp.root"'

	def PrintPlotLimits(self, ti=-5.12e-7, tf=4.606e-6, vmin=-0.7, vmax=0.05):
		print np.double([(tf-ti)/float(self.settings.time_res) +1, ti-self.settings.time_res/2.0,
						 tf+self.settings.time_res/2.0, (vmax-vmin)/self.settings.sigRes, vmin, vmax])

	def PlotHVCurrents(self, name='HVCurrents', cuts='', deltat=30., options='e hist'):
		if not self.in_root_tree:
			return
		print 'Plotting HV current'
		if self.in_root_tree.GetLeaf('timeHV').GetTypeName() == 'TDatime':
			print 'Time format is TDataime. Not implemented yet :P'
		else:
			leng = self.analysisTree.Draw('timeHV.AsDouble()', cuts, 'goff')
			timehv = self.analysisTree.GetVal(0)
			timehv = np.array([timehv[i] for i in xrange(leng)])
			xmin, xmax, deltax = np.floor(timehv.min()), np.ceil(timehv.max()), deltat
			ymin, ymax = self.analysisTree.GetMinimum('currentHV'), self.analysisTree.GetMaximum('currentHV')
			if ymin == ymax:
				ymin -= 1e-10
				ymax += 1e-10
			self.DrawProfile(name, 'timeHV.AsDouble()', xmin, xmax, deltax, 'time', 'currentHV', ymin, ymax, 'Current', cuts, options)
			self.profile[name].SetStats(0)
			# def FitLanGaus(self, name, conv_steps=100, color=ro.kRed, xmin=-10000000, xmax=-10000000):
			self.line[name+'_p'] = ro.TLine(self.profile[name].GetXaxis().GetXmin(), self.currentCut, self.profile[name].GetXaxis().GetXmax(), self.currentCut)
			self.line[name+'_p'].SetLineStyle(2)
			self.line[name+'_p'].SetLineColor(ro.kRed)
			self.line[name+'_p'].Draw('same')
			self.line[name+'_n'] = ro.TLine(self.profile[name].GetXaxis().GetXmin(), -self.currentCut, self.profile[name].GetXaxis().GetXmax(), -self.currentCut)
			self.line[name+'_n'].SetLineStyle(2)
			self.line[name+'_n'].SetLineColor(ro.kRed)
			self.line[name+'_n'].Draw('same')
			self.graph[name] = ro.TGraphErrors(1, np.array([float(xmax + xmin) / 2.], 'f8'), np.array([0], 'f8'), np.array([float(xmax - xmin) / 2.], 'f8'), np.array([self.currentCut], 'f8'))
			self.graph[name].SetNameTitle('g_'+ name, 'g_'+ name)
			self.graph[name].SetFillColor(ro.kRed)
			self.graph[name].SetFillStyle(3003)
			self.graph[name].Draw('same f 2')
			self.profile[name].GetXaxis().SetTimeDisplay(1)

	def PlotDiaVoltage(self, name='DUTVoltage', cuts='', deltat=30, options='e hist'):
		if not self.in_root_tree:
			return
		print 'Plotting DUT voltage'
		if self.in_root_tree.GetLeaf('timeHV').GetTypeName() == 'TDatime':
			print 'Time format is TDataime. Not implemented yet :P'
		else:
			leng = self.analysisTree.Draw('timeHV.AsDouble()', cuts, 'goff')
			timehv = self.analysisTree.GetVal(0)
			timehv = np.array([timehv[i] for i in xrange(leng)])
			xmin, xmax, deltax = np.floor(timehv.min()), np.ceil(timehv.max()), deltat
			ymin, ymax = self.analysisTree.GetMinimum('voltageDia'), self.analysisTree.GetMaximum('voltageDia')
			if ymin == ymax:
				ymin -= 1
				ymax += 1
			self.DrawProfile(name, 'timeHV.AsDouble()', xmin, xmax, deltax, 'time', 'voltageDia', ymin, ymax, 'DUT Voltage [V]', cuts, options)
			self.profile[name].SetStats(0)
			self.line[name + '1_p'] = ro.TLine(self.profile[name].GetXaxis().GetXmin(), self.bias + self.voltageDiaMaxOffset, self.profile[name].GetXaxis().GetXmax(), self.bias + self.voltageDiaMaxOffset)
			self.line[name + '1_p'].SetLineStyle(2)
			self.line[name + '1_p'].SetLineColor(ro.kRed)
			self.line[name + '1_p'].Draw('same')
			self.line[name + '1_n'] = ro.TLine(self.profile[name].GetXaxis().GetXmin(), self.bias - self.voltageDiaMaxOffset, self.profile[name].GetXaxis().GetXmax(), self.bias - self.voltageDiaMaxOffset)
			self.line[name + '1_n'].SetLineStyle(2)
			self.line[name + '1_n'].SetLineColor(ro.kRed)
			self.line[name + '1_n'].Draw('same')
			self.graph[name + '1'] = ro.TGraphErrors(1, np.array([float(xmax + xmin) / 2.], 'f8'), np.array([self.bias], 'f8'), np.array([float(xmax - xmin) / 2.], 'f8'), np.array([self.voltageDiaMaxOffset], 'f8'))
			self.graph[name + '1'].SetNameTitle('g_'+ name + '1', 'g_'+ name + '1')
			self.graph[name + '1'].SetFillColor(ro.kRed)
			self.graph[name + '1'].SetFillStyle(3003)
			self.graph[name + '1'].Draw('same f 2')
			self.profile[name].GetYaxis().SetRangeUser(self.bias - self.voltageDiaMaxOffset, self.bias + self.voltageDiaMaxOffset)
			if self.voltageDiaMaxSigmas != 0:
				yming = self.bias - self.voltageDiaMaxOffset if self.bias - self.voltageDiaMaxOffset < self.voltageDiaMean - self.voltageDiaMaxSigmas * self.voltageDiaSpread else self.voltageDiaMean - self.voltageDiaMaxSigmas * self.voltageDiaSpread
				ymaxg = self.bias + self.voltageDiaMaxOffset if self.bias + self.voltageDiaMaxOffset > self.voltageDiaMean + self.voltageDiaMaxSigmas * self.voltageDiaSpread else self.voltageDiaMean + self.voltageDiaMaxSigmas * self.voltageDiaSpread
				self.profile[name].GetYaxis().SetRangeUser(yming, ymaxg)
				self.line[name + '2_p'] = ro.TLine(self.profile[name].GetXaxis().GetXmin(), self.voltageDiaMean + self.voltageDiaMaxSigmas * self.voltageDiaSpread, self.profile[name].GetXaxis().GetXmax(), self.voltageDiaMean + self.voltageDiaMaxSigmas * self.voltageDiaSpread)
				self.line[name + '2_p'].SetLineStyle(2)
				self.line[name + '2_p'].SetLineColor(ro.kBlue)
				self.line[name + '2_p'].Draw('same')
				self.line[name + '2_n'] = ro.TLine(self.profile[name].GetXaxis().GetXmin(), self.voltageDiaMean - self.voltageDiaMaxSigmas * self.voltageDiaSpread, self.profile[name].GetXaxis().GetXmax(), self.voltageDiaMean - self.voltageDiaMaxSigmas * self.voltageDiaSpread)
				self.line[name + '2_n'].SetLineStyle(2)
				self.line[name + '2_n'].SetLineColor(ro.kBlue)
				self.line[name + '2_n'].Draw('same')
				self.graph[name + '2'] = ro.TGraphErrors(1, np.array([float(xmax + xmin) / 2.], 'f8'), np.array([self.voltageDiaMean], 'f8'), np.array([float(xmax - xmin) / 2.], 'f8'), np.array([self.voltageDiaMaxSigmas * self.voltageDiaSpread], 'f8'))
				self.graph[name + '2'].SetNameTitle('g_' + name + '2', 'g_' + name + '2')
				self.graph[name + '2'].SetFillColor(ro.kBlue)
				self.graph[name + '2'].SetFillStyle(3003)
				self.graph[name + '2'].Draw('same f 2')
			self.profile[name].GetXaxis().SetTimeDisplay(1)

	def FitLanGaus(self, name, conv_steps=100, color=ro.kRed, doToyStats=False):
		if not name in self.canvas.keys() or not name in self.histo.keys():
			print 'Can\'t do Langaus fit as there is a problem with the histogram'
			return
		if not self.canvas[name] or not self.histo[name]:
			print 'Can\'t do Langaus fit as there is a problem with the histogram'
			return
		if self.histo[name].GetEntries() == 0 or self.histo[name].GetMean()== 0 and self.histo[name].GetRMS() == 0:
			print 'Can\'t do Langaus fit as there is a problem with the histogram'
			return
		if self.histo[name].Integral(1, self.histo[name].GetNbinsX()) < 2:
			print 'Can\'t do Langaus fit with only 2 filled bins'
			return
		self.canvas[name].cd()
		self.langaus[name] = LanGaus(self.histo[name])
		self.langaus[name].LanGausFit(conv_steps, xmin=self.fit_min, xmax=self.fit_max)
		xlow, xhigh = self.langaus[name].fit_range['min'], self.langaus[name].fit_range['max']
		self.line[name] = ro.TLine(xlow, 0, xhigh, 0)
		self.line[name].SetLineColor(ro.kViolet + 1)
		self.line[name].SetLineWidth(4)
		fitmean = self.langaus[name].fit.Mean(xlow, xhigh)
		self.histo[name].FindObject('stats').SetOptFit(1)
		ro.gPad.Update()
		self.langaus[name].fit.Draw('same')
		self.langaus[name].fit.SetLineColor(color)
		self.line[name].Draw('same')
		ro.gPad.Update()
		print u'Histo {n}: <PH> = {f} \u00B1 {f2}'.format(n=name, f=self.histo[name].GetMean(), f2=self.histo[name].GetMeanError())
		print 'Fit {n}: <PH> = {f}'.format(n=name, f=fitmean)
		print 'Fit {n}: MP_fit = {f}'.format(n=name, f=self.langaus[name].fit_mp)
		lineToAddInStats = ['Mean_{Fit}', 'MP_{Fit}']
		valuesLineToAddInStats = [fitmean, self.langaus[name].fit_mp]
		if doToyStats:
			self.toy_histos = [self.histo[name].Clone('h_' + name + '_toy_' + str(it)) for it in xrange(100)]
			self.utils.CreateProgressBar(len(self.toy_histos))
			self.utils.bar.start()
			for itj, toyhisto in enumerate(self.toy_histos):
				toyhisto.SetTitle(toyhisto.GetName())
				toyhisto.Reset()
				for it in xrange(self.langaus[name].entries_under_curve):
					toyhisto.Fill(self.langaus[name].fit.GetRandom(self.langaus[name].fit.GetXmin(), self.langaus[name].fit.GetXmax()))
				self.utils.bar.update(itj + 1)
			self.utils.bar.finish()
			toys_means = np.array([toyhisto.GetMean() for toyhisto in self.toy_histos], 'f8')
			toys_means_errors = np.array([toyhisto.GetMeanError() for toyhisto in self.toy_histos], 'f8')
			toys_means_weights = np.divide(1., np.power(toys_means_errors, 2., dtype='f8'), dtype='f8')
			mean_toys = np.average(toys_means, weights=toys_means_weights)
			lineToAddInStats.append('Mean_{100toys}')
			valuesLineToAddInStats.append(mean_toys)
			print '100 Toys Fit {n}: <PH> = {f}'.format(n=name, f=mean_toys)
		self.histo[name].FindObject('stats').SetX1NDC(0.55)
		self.histo[name].FindObject('stats').SetX2NDC(0.9)
		self.histo[name].FindObject('stats').SetY1NDC(0.5)
		self.histo[name].FindObject('stats').SetY2NDC(0.9)
		AddLineToStats(self.canvas[name], lineToAddInStats, valuesLineToAddInStats)
		self.histo[name].SetStats(0)
		self.canvas[name].Modified()
		ro.gPad.Update()

	def FitConvolutedGaussians(self, name, conv_steps=100, color=ro.kRed):
		if not name in self.canvas.keys() or not name in self.histo.keys():
			print 'Can\'t do convoluted gaussians fit as there is aproblem with the histogram'
			return
		if not self.canvas[name] or not self.histo[name]:
			print 'Can\'t do convoluted gaussians fit as there is a problem with the histogram'
			return
		if self.histo[name].GetEntries() == 0 or self.histo[name].GetMean()== 0 and self.histo[name].GetRMS() == 0:
			print 'Can\'t do convoluted gaussians fit as there is a problem with the histogram'
			return
		self.canvas[name].cd()
		# self.langaus[name] = LanGaus(self.histo[name])
		# self.langaus[name].LanGausFit(conv_steps, xmin=self.fit_min, xmax=self.fit_max)
		# xlow, xhigh = self.langaus[name].fit_range['min'], self.langaus[name].fit_range['max']
		xlow, xhigh = max(self.histo[name].GetBinLowEdge(self.histo[name].FindFirstBinAbove(0)), self.histo[name].GetMean() - 3 * self.histo[name].GetRMS()), min(self.histo[name].GetBinLowEdge(self.histo[name].FindLastBinAbove(0) + 1), self.histo[name].GetMean() + 3 * self.histo[name].GetRMS())
		funccg = ro.TF1('fit_' + name, '[0]*exp(-0.5*((x-[1])^2/([2]^2+[3]^2)))', xlow, xhigh)
		funccg.SetParameter(0, self.histo[name].GetMaximum())
		funccg.SetParLimits(0, 0, 100 * self.histo[name].GetMaximum())
		funccg.SetParameter(1, self.histo[name].GetBinCenter(self.histo[name].GetMaximumBin()))
		funccg.SetParLimits(1, self.histo[name].GetMean() - self.histo[name].GetMeanError(), self.histo[name].GetMean() + self.histo[name].GetMeanError())
		funccg.SetParameter(2, (abs(self.histo[name].GetRMS() ** 2 - self.pedestal_sigma ** 2) ** 0.5) if self.histo[name].GetRMS() > self.pedestal_sigma else self.histo[name].GetRMS())
		funccg.SetParLimits(2, 0, 10 * self.histo[name].GetRMS())
		funccg.SetParameter(3, self.pedestal_sigma)
		funccg.SetParLimits(3, self.pedestal_sigma - self.pedestal_sigma_err if self.histo[name].GetRMS() > self.pedestal_sigma else 0., self.pedestal_sigma + self.pedestal_sigma_err if self.histo[name].GetRMS() > self.pedestal_sigma else self.pedestal_sigma / 1000.)
		self.histo[name].Fit('fit_' + name, 'IQMR')
		mean = funccg.GetParameter(1)
		sigm = np.sqrt(funccg.GetParameter(2) ** 2 + funccg.GetParameter(3) ** 2, dtype='f8')
		self.histo[name].Fit('fit_' + name, 'IQM', '', mean - 3 * sigm, mean + 3 * sigm)
		self.langaus[name] = funccg
		self.line[name] = ro.TLine(xlow, 0, xhigh, 0)
		self.line[name].SetLineColor(ro.kViolet + 1)
		self.line[name].SetLineWidth(4)
		fitmean = self.langaus[name].Mean(xlow, xhigh)
		self.histo[name].FindObject('stats').SetOptFit(1)
		ro.gPad.Update()
		self.langaus[name].Draw('same')
		self.langaus[name].SetLineColor(color)
		self.line[name].Draw('same')
		ro.gPad.Update()
		print u'Histo {n}: <PH> = {f} \u00B1 {f2}'.format(n=name, f=self.histo[name].GetMean(), f2=self.histo[name].GetMeanError())
		print 'Fit {n}: <PH> = {f}'.format(n=name, f=fitmean)
		print 'Fit {n}: signal_sigma = {f}'.format(n=name, f=self.langaus[name].GetParameter(2))
		print 'Fit {n}: pedestal_sigma = {f}'.format(n=name, f=self.langaus[name].GetParameter(3))
		xlow, xhigh = self.langaus[name].GetMaximumX(xlow, xhigh) - 5 * self.histo[name].GetRMS(), self.langaus[name].GetMaximumX(xlow, xhigh) + 5 * self.histo[name].GetRMS()
		self.histo[name].GetXaxis().SetRangeUser(xlow, xhigh)
		self.histo[name].FindObject('stats').SetX1NDC(0.55)
		self.histo[name].FindObject('stats').SetX2NDC(0.9)
		self.histo[name].FindObject('stats').SetY1NDC(0.5)
		self.histo[name].FindObject('stats').SetY2NDC(0.9)
		AddLineToStats(self.canvas[name], ['Mean_{Fit}', 'signal_sigma_{Fit}', 'pedestal_sigma_{Fit}'], [fitmean, self.langaus[name].GetParameter(2), self.langaus[name].GetParameter(3)])
		self.histo[name].SetStats(0)
		self.canvas[name].Modified()
		ro.gPad.Update()

	def SaveCanvasInlist(self, lista):
		if not os.path.isdir('{d}'.format(d=self.inDir)):
			ExitMessage('The directory does not exist!!!!', os.EX_UNAVAILABLE)
		for canv in lista:
			if self.canvas.has_key(canv):
				if self.canvas[canv]:
					self.canvas[canv].SaveAs('{d}/{c}.png'.format(d=self.inDir, c=canv))
					self.canvas[canv].SaveAs('{d}/{c}.root'.format(d=self.inDir, c=canv))

	def SaveAllCanvas(self):
		self.SaveCanvasInlist(self.canvas.keys())

	def AddVcalFriend(self, overwr=False, isCF=False):
		appendCF = '' if not isCF else 'CF'
		if not self.analysisTree.GetFriend('vcalTree' + appendCF):
			if os.path.isfile('{d}/{f}.vcal{a}.root'.format(d=self.outDir, f=self.analysisTreeName, a=appendCF)) and not overwr:
				self.analysisTree.AddFriend('vcalTree' + appendCF, '{d}/{f}.vcal{a}.root'.format(d=self.outDir, f=self.analysisTreeName, a=appendCF))
				self.hasBranch['signalVcal' + appendCF] = self.TreeHasBranch('signalVcal' + appendCF)
				self.hasBranch['pedestalVcal' + appendCF] = self.TreeHasBranch('pedestalVcal' + appendCF)
			else:
				if self.hasBranch['signal' + appendCF]:
					self.CreateVcalFriend(isCF)
					self.AddVcalFriend(False, isCF)

	def CreateVcalFriend(self, isCF=False):
		appendCF = '' if not isCF else 'CF'
		if self.hasBranch['signal' + appendCF]:
			if self.signal_cal_folder != '':
				if os.path.isdir(self.signal_cal_folder):
					root_files = glob.glob('{d}/Signal_vs_CalStep_ch*.root'.format(d=self.signal_cal_folder))
					if len(root_files) > 0:
						tempf = ro.TFile(root_files[0], 'READ')
						cal_name = tempf.GetName().strip('.root').split('/')[-1]
						if tempf.FindKey('c_' + cal_name):
							tempc = tempf.Get('c_' + cal_name)
							if tempc.FindObject('fit_' + cal_name):
								tempfit = tempc.GetPrimitive('fit_' + cal_name)
								self.signal_cal_fit_params['p0'] = tempfit.GetParameter(0) / 1000.  # fit in mV
								self.signal_cal_fit_params['p1'] = tempfit.GetParameter(1)
								self.signal_cal_fit_params['p0_error'] = tempfit.GetParError(0) / 1000.  # fit in mV
								self.signal_cal_fit_params['p1_error'] = tempfit.GetParError(1)
								self.signal_cal_fit_params['prob'] = tempfit.GetProb()
								self.signal_cal_fit_params['chi2'] = tempfit.GetChisquare()
								self.signal_cal_fit_params['ndf'] = tempfit.GetNDF()
						tempf.Close()
						if self.signal_cal_fit_params['p1'] != 0:
							self.FillVcalFriend(isCF)
				else:
					ExitMessage('signal cal folder does not exist!')

	def SignalToVcal(self, vsignal, typenp='f4'):
		if self.signal_cal_fit_params['p1'] != 0:
			return np.divide(np.subtract(vsignal, self.signal_cal_fit_params['p0'], dtype='f8'), self.signal_cal_fit_params['p1'], dtype='f8').astype(typenp)
		return 0

	def FillVcalFriend(self, isCF=False):
		appendCF = '' if not isCF else 'CF'
		self.LoadSignalScalars()
		print 'Creating vcal{a} root tree:'.format(a=appendCF)
		vcalfile = ro.TFile('{d}/{f}.vcal{a}.root'.format(d=self.outDir, f=self.analysisTreeName, a=appendCF), 'RECREATE')
		vcaltree = ro.TTree('vcalTree' + appendCF, 'vcalTree' + appendCF)
		vcalSignalEv = np.zeros(1, 'f4')
		vcalPedEv = np.zeros(1, 'f4')
		vcaltree.Branch('signalVcal' + appendCF, vcalSignalEv, 'signalVcal{a}/F'.format(a=appendCF))
		vcaltree.Branch('pedestalVcal' + appendCF, vcalPedEv, 'pedestalVcal{a}/F'.format(a=appendCF))
		totEvents = self.settings.num_events
		self.utils.CreateProgressBar(totEvents + 1)
		self.utils.bar.start()
		for ev in xrange(totEvents):
			if ev in self.eventVect:
				try:
					argum = np.argwhere(ev == self.eventVect).flatten()
					sigVcal = self.SignalToVcal(self.sigVect[argum]) if not isCF else self.SignalToVcal(self.sigVectCF[argum])
					pedVcal = self.SignalToVcal(self.pedVect[argum]) if not isCF else self.SignalToVcal(self.pedVectCF[argum])
					vcalSignalEv.itemset(sigVcal)
					vcalPedEv.itemset(pedVcal)
				except ValueError:
					ExitMessage('Could not fill event {ev}. :S'.format(ev=ev))
			else:
				vcalSignalEv.itemset(0)
				vcalPedEv.itemset(0)
			vcaltree.Fill()
			self.utils.bar.update(ev + 1)
		vcalfile.Write()
		vcalfile.Close()
		self.utils.bar.finish()
		print 'Finished creating vcalTree{a} in {d}/{f}.vcal{a}.root'.format(d=self.outDir, f=self.analysisTreeName, a=appendCF)

	def AddChargeFriend(self, overwr=False, isCF=False):
		appendCF = '' if not isCF else 'CF'
		if not self.analysisTree.GetFriend('chargeTree' + appendCF):
			if os.path.isfile('{d}/{f}.charge{a}.root'.format(d=self.outDir, f=self.analysisTreeName, a=appendCF)) and not overwr:
				self.analysisTree.AddFriend('chargeTree' + appendCF, '{d}/{f}.charge{a}.root'.format(d=self.outDir, f=self.analysisTreeName, a=appendCF))
				self.hasBranch['signalCharge' + appendCF] = self.TreeHasBranch('signalCharge' + appendCF)
				self.hasBranch['pedestalCharge' + appendCF] = self.TreeHasBranch('pedestalCharge' + appendCF)
			else:
				if self.hasBranch['signal' + appendCF]:
					self.CreateChargeFriend(isCF)
					self.AddChargeFriend(False, isCF)

	def CreateChargeFriend(self, isCF=False):
		if self.hasBranch['signal']:
			if self.cal_circuit_settings != '':
				if os.path.isfile(self.cal_circuit_settings):
					self.vcal_to_q = VcalToElectrons(self.cal_circuit_settings)
			self.FillChargeFriend(isCF)

	def FillChargeFriend(self, isCF=False):
		appendCF = '' if not isCF else 'CF'
		self.LoadSignalScalars()
		print 'Creating charge{a} root tree:'.format(a=appendCF)
		chargefile = ro.TFile('{d}/{f}.charge{a}.root'.format(d=self.outDir, f=self.analysisTreeName, a=appendCF), 'RECREATE')
		chargetree = ro.TTree('chargeTree' + appendCF, 'chargeTree' + appendCF)
		chargeSignalEv = np.zeros(1, 'f4')
		chargePedEv = np.zeros(1, 'f4')
		chargetree.Branch('signalCharge' + appendCF, chargeSignalEv, 'signalCharge{a}/F'.format(a=appendCF))
		chargetree.Branch('pedestalCharge' + appendCF, chargePedEv, 'pedestalCharge{a}/F'.format(a=appendCF))
		totevents = self.in_root_tree.GetEntries()
		self.utils.CreateProgressBar(totevents + 1)
		self.utils.bar.start()
		for ev in xrange(totevents):
			if ev in self.eventVect:
				try:
					argum = np.argwhere(ev == self.eventVect).flatten()
					sigVcal = self.SignalToVcal(self.sigVect[argum] if not isCF else self.sigVectCF[argum]) * 1000.  # in mV
					pedVcal = self.SignalToVcal(self.pedVect[argum] if not isCF else self.pedVectCF[argum]) * 1000.  # in mV
					chargeSignalEv.itemset(self.vcal_to_q.Q_in_e_from_mV(sigVcal).nominal_value)
					chargePedEv.itemset(self.vcal_to_q.Q_in_e_from_mV(pedVcal).nominal_value)
				except ValueError:
					ExitMessage('Could not fill event {ev}. :S'.format(ev=ev))
			else:
				chargeSignalEv.itemset(0)
				chargePedEv.itemset(0)
			chargetree.Fill()
			self.utils.bar.update(ev + 1)
		chargefile.Write()
		chargefile.Close()
		self.utils.bar.finish()
		print 'Finished creating chargeTree{a} in {d}/{f}.charge{a}.root'.format(d=self.outDir, f=self.analysisTreeName, a=appendCF)

	def AddTimeCFFriend(self, delay=1800, attPercent=75, overwr=False):
		if not self.analysisTree.GetFriend('timeCFTree'):
			if os.path.isfile('{d}/{f}.timeCF.{de}.{att}.root'.format(d=self.outDir, f=self.analysisTreeName, de=delay, att=attPercent)) and not overwr:
				self.analysisTree.AddFriend('timeCFTree', '{d}/{f}.timeCF.{de}.{att}.root'.format(d=self.outDir, f=self.analysisTreeName, de=delay, att=attPercent))
				self.hasBranch['deltaTimeCF'] = self.TreeHasBranch('deltaTimeCF')
			else:
				if self.hasBranch['time']:
					self.CreateTimeCFFriend(delay, attPercent)
					self.AddTimeCFFriend(delay, attPercent, False)

	def CreateTimeCFFriend(self, delay=1800, attPercent=75):
		entries = self.in_root_tree.GetEntries()
		print 'Creating shifted constant fraction time tree:'
		if self.settings.is_cal_run:
			print 'Constant Fraction shifts will be 0 because it is a calibration run'
		crossPosition0 = 1.2195e-6
		ndelay = np.floor(np.multiply(delay, 1e-9, dtype='f8') / self.settings.time_res + 0.5).astype('int32')
		func = ro.TF1('fit_CF', '[0]*cos([1]*x+pi/(2e-6*[1])-1+[2])+[3]', 0.5e-6, 2e-6)
		func.SetNpx(10000)
		timeCFfile = ro.TFile('{d}/{f}.timeCF.{de}.{att}.root'.format(d=self.outDir, f=self.analysisTreeName, de=delay, att=attPercent), 'RECREATE')
		timeCFtree = ro.TTree('timeCFTree', 'timeCFTree')
		# shiftTime = np.zeros(1, 'f8')
		# timeCF = np.zeros(self.settings.points, 'f8')
		deltaTimeCF = np.zeros(1, 'f8')
		# timeCFtree.Branch('timeCF', timeCF, 'timeCF[{p}]/D'.format(p=self.settings.points))
		timeCFtree.Branch('deltaTimeCF', deltaTimeCF, 'deltaTimeCF/D')
		self.utils.CreateProgressBar(entries)
		self.utils.bar.start()
		# self.timeCFVect = []
		self.time_CF_deltas = []
		# for it, timei in enumerate(self.timeVect):
		for ev in xrange(entries):
			deltaTimeCF.fill(-100e-6)
			if ev in self.eventVect:
				if not self.settings.is_cal_run:
					it = np.argwhere(ev == self.eventVect).flatten()[0]
					timei = self.timeVect[it]
					sigi = self.signalWaveVect[it]
					sig0 = np.multiply(-0.01 * attPercent, np.roll(sigi, -ndelay), dtype='f8')
					sig1 = np.add(sigi, sig0, dtype='f8')
					grafi = ro.TGraph(timei.size, timei, sig1)
					# name = 'CFSignal' + str(it)
					# grafi.SetName('g_' + name)
					grafi.GetXaxis().SetRangeUser(0, 5e-6)
					p0 = 0.100 if self.bias > 0 else -0.100
					p0l = [0, 1] if self.bias > 0 else [-1, 0]
					p1 = 1.7e6
					p1l = [1e6, 3e6]
					p2 = -0.25
					p2l = [-1.5, 1]
					p3 = 0
					p3l = [-1, 1]
					func.SetParameters(np.array([p0, p1, p2, p3], 'f8'))
					func.SetParLimits(0, p0l[0], p0l[1])
					func.SetParLimits(1, p1l[0], p1l[1])
					func.SetParLimits(2, p2l[0], p2l[1])
					func.SetParLimits(3, p3l[0], p3l[1])
					tfit = grafi.Fit('fit_CF', 'QM0RSB', '')
					if tfit.Prob() < 0.9:
						tfit = grafi.Fit('fit_CF', 'QM0RSB', '')
					if tfit.Prob() < 0.9:
						tfit = grafi.Fit('fit_CF', 'QM0RSB', '')
					# shiftTime.fill(0)
					if tfit.Prob() >= 0.01:
						cfCrossing = func.GetX(0, 0.5e-6, 2e-6)
						if (func.Derivative(cfCrossing) < 0 < self.bias) or (func.Derivative(cfCrossing) > 0 > self.bias):
							# shiftTime.fill(np.subtract(crossPosition0, cfCrossing, dtype='f8'))
							deltaTimeCF.fill(np.subtract(crossPosition0, cfCrossing, dtype='f8'))
					# np.putmask(timeCF, np.ones(self.settings.points, '?'), np.add(timei, shiftTime, dtype='f8'))
					# self.timeCFVect.append(timeCF)
				else:
					deltaTimeCF.fill(0)
				self.time_CF_deltas.append(deltaTimeCF[0])
			timeCFtree.Fill()
			self.utils.bar.update(ev + 1)
		timeCFfile.Write()
		timeCFfile.Close()
		self.utils.bar.finish()
		# self.timeCFVect = np.array(self.timeCFVect, 'f8')
		self.time_CF_deltas = np.array(self.time_CF_deltas, 'f8')
		print 'Finished creating timeCFTree in {d}/{f}.timeCF.{de}.{att}.root'.format(d=self.outDir, f=self.analysisTreeName, de=delay, att=attPercent)

	def RemovePedestalFromSignal(self, namePlot0, namePlotZoom, pedSigma, xmin=10000, xmax=-10000, namePlotNewSuffix='NoPedestal'):
		namec = namePlot0+namePlotNewSuffix
		if not self.histo.has_key(namePlotZoom):
			return namec
		elif self.histo[namePlotZoom].Integral(0, self.histo[namePlotZoom].GetNbinsX()) < 1:
			return namec
		if self.histo.has_key(namePlot0) and self.histo.has_key(namePlotZoom):
			if self.histo.has_key(namec):
				if self.histo[namec]:
					self.histo[namec].Delete()
				del self.histo[namec]
			self.histo[namec] = self.histo[namePlot0].Clone('h_' + namec)
			histo0 = self.histo[namePlot0]
			histoz = self.histo[namePlotZoom]
			histoc = self.histo[namec]
			xlow = histoz.GetBinLowEdge(1) if xmin == 10000 else xmin
			xhigh = histoz.GetBinLowEdge(histoz.GetNbinsX()) if xmax == -10000 else xmax
			funcg = ro.TF1('fit_' + namePlotZoom, 'gaus', xlow, xhigh)
			funcg.SetNpx(1000)
			funcg.SetParameter(0, 10)
			funcg.SetParLimits(0, 0, 1e6)
			par1l, par1h = max(xlow, histoz.GetBinLowEdge(histoz.FindFirstBinAbove(0))), min(xhigh, 1000 if 'charge' in namec.lower() else 10)
			funcg.SetParameter(1, 0.5 * (par1l + par1h))
			funcg.SetParLimits(1, xlow, par1h)
			funcg.SetParameter(2, pedSigma * 2)
			funcg.SetParLimits(2, 0.5 * pedSigma, 2.1 * pedSigma)
			fit = histoz.Fit('fit_' + namePlotZoom, 'QEMBS', '', par1l, par1h)
			params = np.array([fit.Parameter(0), fit.Parameter(1), fit.Parameter(2)], 'f8')
			fit = histoz.Fit('fit_' + namePlotZoom, 'QEMBS', '', par1l, min(params[1] + 5 * params[2], par1h))
			SetDefaultFitStats(histoz, funcg)
			params = np.array([fit.Parameter(0), fit.Parameter(1), fit.Parameter(2)], 'f8')
			histozBinW = histoz.GetBinWidth(1)
			n0, mu, sigma = params[0], params[1], params[2]
			maxc = max(histo0.GetBinLowEdge(histo0.GetMaximumBin()), self.langaus[namePlot0].fit.GetMaximumX()) / 4.0
			xEndCompensation = min(params[1] + 3 * params[2], maxc)
			bini = 1
			xlow, xhigh = histo0.GetBinLowEdge(bini), histo0.GetBinLowEdge(bini + 1)
			doBin = xhigh < xEndCompensation
			while doBin:
				compBin = RoundInt(n0 * np.sqrt(np.pi / 2.) * sigma * (np.math.erf((xhigh - mu)/ (sigma * np.sqrt(2))) - np.math.erf((xlow - mu)/ (sigma * np.sqrt(2)))) / histozBinW)
				# compBin = -1 * min(compBin, histo0.GetBinContent(bini))
				# histoc.Fill(histo0.GetBinCenter(bini), compBin)
				histoc.SetBinContent(bini, int(max(0, histo0.GetBinContent(bini) - compBin)))
				bini += 1
				xlow, xhigh = histo0.GetBinLowEdge(bini), histo0.GetBinLowEdge(bini + 1)
				doBin = xhigh < xEndCompensation
			xhigh = min(xhigh, xEndCompensation)
			compBin = RoundInt(n0 * np.sqrt(np.pi / 2.) * sigma * (np.math.erf((xhigh - mu) / (sigma * np.sqrt(2))) - np.math.erf((xlow - mu) / (sigma * np.sqrt(2)))) / histozBinW)
			# compBin = -1 * min(compBin, histo0.GetBinContent(bini))
			# histoc.Fill(histo0.GetBinCenter(bini), compBin)
			histoc.SetBinContent(bini, max(0, histo0.GetBinContent(bini) - compBin))
			histoc.ResetStats()
			if self.canvas.has_key(namec):
				if self.canvas[namec]:
					self.canvas[namec].Close()
				del self.canvas[namec]
			self.canvas[namec] = ro.TCanvas('c_' + namec, 'c_' + namec, 1)
			self.canvas[namec].cd()
			histoc.Draw('e')
			if histoc.FindObject('mystats'):
				histoc.FindObject('mystats').Delete()
			histoc.SetStats(False)
			histoc.SetStats(True)
			self.canvas[namec].SetGridx()
			self.canvas[namec].SetGridy()
			self.canvas[namec].SetTicky()
			ro.gPad.Update()
			SetDefault1DStats(self.histo[namec])
			ro.gPad.Update()
			return namec

	def PlotConstantFractionSignal(self, event=0, delay=1.8e-6, att=0.75):
		ndelay = np.floor(delay / self.settings.time_res + 0.5).astype('int32')
		leng = self.analysisTree.Draw('voltageSignal:time', 'event==' + str(event), 'goff')
		sig = self.analysisTree.GetVal(0)
		sig = np.array([sig[i] * 1e3 for i in xrange(leng)], 'f8')
		tim = self.analysisTree.GetVal(1)
		tim = np.array([tim[i] * 1e6 for i in xrange(leng)], 'f8')
		sig0 = np.multiply(-att, np.roll(sig, -ndelay), dtype='f8')
		sig1 = np.add(sig, sig0, dtype='f8')
		graf = ro.TGraph(len(tim), tim, sig1)
		name = 'CFSignal' + str(event)
		graf.SetName('g_' + name)
		graf.SetTitle('g_' + name)
		graf.GetXaxis().SetTitle('time [us]')
		graf.GetYaxis().SetTitle('CF signal [mV]')
		graf.GetXaxis().SetRangeUser(0, 5)
		# self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
		# self.canvas[name].SetGridx()
		# self.canvas[name].SetGridy()
		self.graph[name] = graf
		# self.graph[name].Draw('AP')
		func = ro.TF1('fit_' + name,'[0]*cos([1]*x+pi/(2*[1])-1+[2])+[3]', 0.5, 2)
		p0 = 100 if self.bias > 0 else -100
		p0l = [0, 1000] if self.bias > 0 else [-1000, 0]
		p1 = 1.7
		p1l = [1, 3]
		p2 = -0.25
		p2l = [-1.5, 1]
		p3 = 0
		p3l = [-1000, 1000]
		func.SetNpx(10000)
		func.SetLineColor(ro.kRed)
		func.SetParameters(np.array([p0, p1, p2, p3], 'f8'))
		func.SetParLimits(0, p0l[0], p0l[1])
		func.SetParLimits(1, p1l[0], p1l[1])
		func.SetParLimits(2, p2l[0], p2l[1])
		func.SetParLimits(3, p3l[0], p3l[1])
		# func = ro.TF1('fit_' + name,'[0]*cos([1]*x+[2])+[3]', 0.5, 2)
		# p0 = 100 if self.bias > 0 else -100
		# p0l = [0, 1000] if self.bias > 0 else [-1000, 0]
		# p1 = 1.7
		# p1l = [1, 3]
		# p2 = -0.5
		# p2l = [-1, 0.5]
		# p3 = 0
		# p3l = [-1000, 1000]
		# func.SetNpx(10000)
		# func.SetLineColor(ro.kRed)
		# func.SetParameters(np.array([p0, p1, p2, p3], 'f8'))
		# func.SetParLimits(0, p0l[0], p0l[1])
		# func.SetParLimits(1, p1l[0], p1l[1])
		# func.SetParLimits(2, p2l[0], p2l[1])
		# func.SetParLimits(3, p3l[0], p3l[1])
		tfit = self.graph[name].Fit('fit_' + name, 'QM0RS', '')
		if tfit.Prob() < 0.9:
			tfit = self.graph[name].Fit('fit_' + name, 'QM0RS', '')
		if tfit.Prob() < 0.9:
			tfit = self.graph[name].Fit('fit_' + name, 'QM0RS', '')
		# print 'fit prob', tfit.Prob() * 100, '%'
		# print 'Parameters', tfit.Parameter(0), tfit.Parameter(1), tfit.Parameter(2), tfit.Parameter(3)
		self.langaus[name] = func.Clone()
		# self.canvas[name].cd()
		# self.langaus[name].Draw('same')
		# self.canvas[name].Modified()
		# ro.gPad.Update()
		peakbla = -100000
		peakbla2 = -100000
		peakposition2 = -100000
		if tfit.Prob() >= 0.01:
			peakbla2 = self.langaus[name].GetX(0, 0.5, 2)
			if (self.langaus[name].Derivative(peakbla2) < 0 < self.bias) or (self.langaus[name].Derivative(peakbla2) > 0 > self.bias):
				# peak offset 0.023 behind for calibration run for +250V in 20200916 V2
				# crossPosition0 = 1.19654 + 0.023 = 1.21954
				crossPosition0 = 1.2195 # using the pol0 fit for 20200316 cal 0 which gives a peak position at 2.124 us. This is the offset needed tested with some measurements.
				deltaCrossPosition = crossPosition0 - peakbla2
				peakbla = peakbla2
				# peakbla += 0.90584 # for a peak position of ~2.10159 us
				peakbla += 0.928 # for a peak position of ~2.124 us
				peakposition2 = 1 + 9.16383e-01 * peakbla2
		# print 'Peak position set to', peakbla
		return peakbla2

	def HistoBla(self):
		blah0 = ro.TH1F('blah0', 'blah0', 500, 1, 3)
		blah1 = ro.TH1F('blah1', 'blah1', 500, -1, 1)
		blah2 = ro.TH1F('blah2', 'blah2', 1000, -1, 3)
		blah3 = ro.TH2F('blah3', 'blah3', 1000, -1, 3, 500, 1, 3)
		for i in xrange(int(self.analysisTree.GetEntries())):
			if i % 1000 == 0: print i
			self.analysisTree.GetEntry(i)
			pi = self.analysisTree.peakPosition * 1e6
			if 1 < pi < 3:
				blah0.Fill(pi)
				blah1.Fill(2.124 - pi)
			pi2 = self.PlotConstantFractionSignal(i)
			if -1 < pi2 < 3:
				blah2.Fill(pi2)
			if 1 < pi < 3 and -1 < pi2 < 3:
				blah3.Fill(pi2, pi)
		blac0 = ro.TCanvas('blac0', 'blac0', 1)
		blah0.Draw('e')
		blac1 = ro.TCanvas('blac1', 'blac1', 1)
		blah1.Draw('e')
		blac2 = ro.TCanvas('blac2', 'blac2', 1)
		blah2.Draw('e')
		blac3 = ro.TCanvas('blac3', 'blac3', 1)
		blah3.Draw('colz')
		ipdb.set_trace()

	def DoAll(self, doDebug=False, saveCanvases=True):
		ddbg = doDebug or self.doDebug
		overw = self.overw
		if self.in_root_tree:
			self.AnalysisWaves()
			if self.emptyAnalysis: return
			if ddbg: ipdb.set_trace()
			self.AddVcalFriend(overw, False)
			self.AddVcalFriend(overw, True)
			if ddbg: ipdb.set_trace()
			self.AddChargeFriend(overw, False)
			self.AddChargeFriend(overw, True)
			if ddbg: ipdb.set_trace()
			self.PlotPedestal('Pedestal', cuts=self.cut0.GetTitle(), branch='pedestal')
			self.PlotPedestal('PedestalCF', cuts=self.cut0CF.GetTitle(), branch='pedestalCF')
			if ddbg: ipdb.set_trace()
			if not self.is_cal_run:
				self.PlotPedestal('PedestalVcal', cuts=self.cut0.GetTitle(), branch='pedestalVcal')
				self.PlotPedestal('PedestalVcalCF', cuts=self.cut0CF.GetTitle(), branch='pedestalVcalCF')
				if ddbg: ipdb.set_trace()
				self.PlotPedestal('PedestalCharge', cuts=self.cut0.GetTitle(), branch='pedestalCharge')
				self.PlotPedestal('PedestalChargeCF', cuts=self.cut0CF.GetTitle(), branch='pedestalChargeCF')
				if ddbg: ipdb.set_trace()
			self.PlotWaveforms('SelectedWaveforms', 'signal', cuts=self.cut0.GetTitle(), do_logz=True, do_cf=False)
			self.PlotWaveforms('SelectedWaveformsCF', 'signal', cuts=self.cut0CF.GetTitle(), do_logz=True, do_cf=True)
			if ddbg: ipdb.set_trace()
			# if 'SelectedWaveforms' in self.canvas.keys():
			# 	self.canvas['SelectedWaveforms'].SetLogz()
			self.PlotWaveforms('SelectedWaveformsPedCor', 'signal_ped_corrected', cuts=self.cut0.GetTitle(), do_logz=True, do_cf=False)
			self.PlotWaveforms('SelectedWaveformsPedCorCF', 'signal_ped_corrected', cuts=self.cut0CF.GetTitle(), do_logz=True, do_cf=True)
			if ddbg: ipdb.set_trace()
			# if 'SelectedWaveformsPedCor' in self.canvas.keys():
			# 	self.canvas['SelectedWaveformsPedCor'].SetLogz()
			self.PlotSignal('PH', cuts=self.cut0.GetTitle(), branch='signal')
			self.PlotSignal('PH_CF', cuts=self.cut0CF.GetTitle(), branch='signalCF')
			if ddbg: ipdb.set_trace()
			if not self.is_cal_run:
				self.FitLanGaus('PH')
				self.FitLanGaus('PH_CF')
				self.PlotSignal('PHvcal', cuts=self.cut0.GetTitle(), branch='signalVcal')
				self.PlotSignal('PHvcal_CF', cuts=self.cut0CF.GetTitle(), branch='signalVcalCF')
				if ddbg: ipdb.set_trace()
				self.FitLanGaus('PHvcal')
				self.FitLanGaus('PHvcal_CF')
				self.PlotSignal('PHcharge', cuts=self.cut0.GetTitle(), branch='signalCharge')
				self.PlotSignal('PHcharge_CF', cuts=self.cut0CF.GetTitle(), branch='signalChargeCF')
				if ddbg: ipdb.set_trace()
				self.FitLanGaus('PHcharge')
				self.FitLanGaus('PHcharge_CF')

				phPpos, phGsigma = max(self.histo['PH'].GetBinLowEdge(self.histo['PH'].GetMaximumBin()), self.langaus['PH'].fit.GetMaximumX()), self.langaus['PH'].fit.GetParameter(3)
				if phPpos > 10:
					minx, maxx = -10, max(max(TruncateFloat(phPpos / 4., 10), 40), phPpos - 3.5 * phGsigma)
					binsx = RoundInt((maxx - minx) / 2.)
					self.PlotSignal('PHPedestal', binsx, cuts=self.cut0.GetTitle(), branch='signal', minx=minx, maxx=maxx, optimizeBinning=False)
					histocName = self.RemovePedestalFromSignal('PH', 'PHPedestal', self.pedestal_sigma, xmax=max(TruncateFloat(phPpos / 4., 10), TruncateFloat(phPpos - 3.5 * phGsigma, 10)))
					self.FitLanGaus(histocName)
				phvcalPpos, phvcalGsigma = max(self.histo['PHvcal'].GetBinLowEdge(self.histo['PHvcal'].GetMaximumBin()), self.langaus['PHvcal'].fit.GetMaximumX()), self.langaus['PHvcal'].fit.GetParameter(3)
				if phvcalPpos > 10:
					minx, maxx = -10, max(max(TruncateFloat(phvcalPpos / 4., 10), 40), phvcalPpos - 3.5 * phvcalGsigma)
					binsx = RoundInt((maxx - minx) / 2.)
					self.PlotSignal('PHvcalPedestal', binsx, cuts=self.cut0.GetTitle(), branch='signalVcal', minx=minx, maxx=maxx, optimizeBinning=False)
					histocName = self.RemovePedestalFromSignal('PHvcal', 'PHvcalPedestal', self.pedestal_vcal_sigma, xmax=max(TruncateFloat(phvcalPpos / 4., 10), TruncateFloat(phvcalPpos - 3.5 * phvcalGsigma, 10)))
					self.FitLanGaus(histocName)
				phchPpos, phchGsigma = max(self.histo['PHcharge'].GetBinLowEdge(self.histo['PHcharge'].GetMaximumBin()), self.langaus['PHcharge'].fit.GetMaximumX()), self.langaus['PHcharge'].fit.GetParameter(3)
				if phchPpos > 1200:
					minx, maxx = -10, max(max(TruncateFloat(phvcalPpos / 4., 10), 40), phvcalPpos - 3.5 * phvcalGsigma)
					binsx = RoundInt((maxx - minx) / 2.)
					self.PlotSignal('PHchargePedestal', binsx, cuts=self.cut0.GetTitle(), branch='signalCharge', minx=minx, maxx=maxx, optimizeBinning=False)
					histocName = self.RemovePedestalFromSignal('PHcharge', 'PHchargePedestal', self.pedestal_charge_sigma, xmax=max(TruncateFloat(phchPpos / 4., 1000), TruncateFloat(phchPpos - 3.5 * phchGsigma, 1000)))
					self.FitLanGaus(histocName)

				phPpos, phGsigma = max(self.histo['PH_CF'].GetBinLowEdge(self.histo['PH_CF'].GetMaximumBin()), self.langaus['PH_CF'].fit.GetMaximumX()), self.langaus['PH_CF'].fit.GetParameter(3)
				if phPpos > 10:
					minx, maxx = -10, max(max(TruncateFloat(phPpos / 4., 10), 40), phPpos - 3.5 * phGsigma)
					binsx = RoundInt((maxx - minx) / 2.)
					self.PlotSignal('PHPedestal_CF', binsx, cuts=self.cut0CF.GetTitle(), branch='signalCF', minx=minx, maxx=maxx, optimizeBinning=False)
					histocName = self.RemovePedestalFromSignal('PH_CF', 'PHPedestal_CF', self.pedestal_sigmaCF, xmax=max(TruncateFloat(phPpos / 4., 10), TruncateFloat(phPpos - 3.5 * phGsigma, 10)))
					self.FitLanGaus(histocName)
				phvcalPpos, phvcalGsigma = max(self.histo['PHvcal_CF'].GetBinLowEdge(self.histo['PHvcal_CF'].GetMaximumBin()), self.langaus['PHvcal_CF'].fit.GetMaximumX()), self.langaus['PHvcal_CF'].fit.GetParameter(3)
				if phvcalPpos > 10:
					minx, maxx = -10, max(max(TruncateFloat(phvcalPpos / 4., 10), 40), phvcalPpos - 3.5 * phvcalGsigma)
					binsx = RoundInt((maxx - minx) / 2.)
					self.PlotSignal('PHvcalPedestal_CF', binsx, cuts=self.cut0CF.GetTitle(), branch='signalVcalCF', minx=minx, maxx=maxx, optimizeBinning=False)
					histocName = self.RemovePedestalFromSignal('PHvcal_CF', 'PHvcalPedestal_CF', self.pedestal_vcal_sigmaCF, xmax=max(TruncateFloat(phvcalPpos / 4., 10), TruncateFloat(phvcalPpos - 3.5 * phvcalGsigma, 10)))
					self.FitLanGaus(histocName)
				phchPpos, phchGsigma = max(self.histo['PHcharge_CF'].GetBinLowEdge(self.histo['PHcharge_CF'].GetMaximumBin()), self.langaus['PHcharge_CF'].fit.GetMaximumX()), self.langaus['PHcharge_CF'].fit.GetParameter(3)
				if phchPpos > 1200:
					minx, maxx = -10, max(max(TruncateFloat(phvcalPpos / 4., 10), 40), phvcalPpos - 3.5 * phvcalGsigma)
					binsx = RoundInt((maxx - minx) / 2.)
					self.PlotSignal('PHchargePedestal_CF', binsx, cuts=self.cut0CF.GetTitle(), branch='signalChargeCF', minx=minx, maxx=maxx, optimizeBinning=False)
					histocName = self.RemovePedestalFromSignal('PHcharge_CF', 'PHchargePedestal_CF', self.pedestal_charge_sigmaCF, xmax=max(TruncateFloat(phchPpos / 4., 1000), TruncateFloat(phchPpos - 3.5 * phchGsigma, 1000)))
					self.FitLanGaus(histocName)
				self.PlotHVCurrents('HVCurrents', '', 5)
				if ddbg: ipdb.set_trace()
				self.PlotDiaVoltage('DUTVoltage', '', 5)
				if ddbg: ipdb.set_trace()

			else:
				self.FitConvolutedGaussians('PH')
				self.FitConvolutedGaussians('PH_CF')
			if ddbg: ipdb.set_trace()
			if saveCanvases:
				self.SaveAllCanvas()


# def main():
# 	parser = OptionParser()
# 	parser.add_option('-d', '--inDir', dest='inDir', default='.', type='string', help='Directory containing the run files')
# 	parser.add_option('-c', '--configFile', dest='config', default='CAENAnalysisConfig.cfg', type='string', help='Path to file containing Analysis configuration file')
# 	parser.add_option('-i', '--inputFile', dest='infile', default='', type='string', help='Path to root file to be analysed')
# 	parser.add_option('-b', '--bias', dest='bias', default=0, type='float', help='Bias voltage used')
# 	parser.add_option('-v', '--verbose', dest='verb', default=False, help='Toggles verbose', action='store_true')
# 	parser.add_option('-a', '--automatic', dest='auto', default=False, help='Toggles automatic basic analysis', action='store_true')
#
# 	(options, args) = parser.parse_args()
# 	directory = str(options.inDir)
# 	config = str(options.config)
# 	infile = str(options.infile)
# 	bias = float(options.bias)
# 	verb = bool(options.verb)
# 	autom = bool(options.auto)
#
# 	ana = AnalysisCaenCCD(directory, config, infile, bias, verb)
#
# 	# ana.LoadAnalysisTree()
# 	# ana.LoadPickles()
# 	if autom:
# 		ana.AnalysisWaves()
# 		ana.PlotWaveforms('SelectedWaveforms', 'signal', cuts=ana.cut0.GetTitle())
# 		ana.canvas['SelectedWaveforms'].SetLogz()
# 		ana.PlotSignal('PH', cuts=ana.cut0.GetTitle())
# 		ana.FitLanGaus('PH')
# 		ana.PlotPedestal('Pedestal', cuts=ana.cut0.GetTitle())
# 		ana.SaveAllCanvas()
# 	return ana
#
# 	# if auto:
# 	# 	ccd.StartHVControl()
# 	# 	ccd.AdjustBaseLines()
# 	# 	ccd.SavePickles()
# 	# 	written_events = ccd.GetData()
# 	# 	ccd.settings.num_events = written_events
# 	# 	ccd.SavePickles()  # update pickles with the real amount of written events
# 	# 	ccd.settings.MoveBinaryFiles()
# 	# 	ccd.settings.RenameDigitiserSettings()
# 	# 	ccd.CloseHVClient()
# 	# 	if not ccd.settings.simultaneous_conversion:
# 	# 		ccd.CreateRootFile(files_moved=True)
# 	# 		while ccd.pconv.poll() is None:
# 	# 			continue
# 	# 		ccd.CloseSubprocess('converter', stdin=False, stdout=False)
# 	#
# 	# print 'Finished :)'
# 	# sys.stdout.write('\a\a\a')
# 	# sys.stdout.flush()

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-d', '--inDir', dest='inDir', default='.', type='string', help='Directory containing the run files')
	parser.add_option('-c', '--configFile', dest='config', default='CAENAnalysisConfig.cfg', type='string', help='Path to file containing Analysis configuration file')
	parser.add_option('-i', '--inputFile', dest='infile', default='', type='string', help='Path to root file to be analysed')
	parser.add_option('-b', '--bias', dest='bias', default=0, type='float', help='Bias voltage used')
	parser.add_option('-v', '--verbose', dest='verb', default=False, help='Toggles verbose', action='store_true')
	parser.add_option('-a', '--automatic', dest='auto', default=False, help='Toggles automatic basic analysis', action='store_true')
	parser.add_option('-o', '--overwrite', dest='overwrite', default=False, help='Toggles overwriting of the analysis tree', action='store_true')
	parser.add_option('--batch', dest='dobatch', default=False, help='Toggles batch mode with no plotting', action='store_true')
	parser.add_option('--debug', dest='dodebug', default=False, help='Toggles debug mode', action='store_true')

	(options, args) = parser.parse_args()
	directory = str(options.inDir)
	config = str(options.config)
	infile = str(options.infile)
	bias = float(options.bias)
	verb = bool(options.verb)
	verb = True or verb
	autom = bool(options.auto)
	overw = bool(options.overwrite)
	doBatch = bool(options.dobatch)
	doDebug = bool(options.dodebug)

	ana = AnalysisCaenCCD(directory, config, infile, bias, overw, verb, doDebug)

	ro.gROOT.SetBatch(0)
	if doBatch:
		print 'Running in Batch mode!'
		ro.gROOT.SetBatch(1)

	# ana.LoadAnalysisTree()
	# ana.LoadPickles()
	if autom:
		ana.DoAll(doDebug)
	# return ana
	# ana = main()
