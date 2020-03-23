#!/usr/bin/env python
import numpy as np
from struct import unpack
import time, os, sys
from optparse import OptionParser
import ipdb
import cPickle as pickle
from Utils import *
from Settings_Caen import Settings_Caen
from Converter_Caen import Converter_Caen
import subprocess as subp
from AnalysisCaenCCD import AnalysisCaenCCD
import ROOT as ro

fit_method = ('Minuit2', 'Migrad', )

class AnalysisAllRunsInFolder:
	def __init__(self, runsdir='', config='', configinput='', overwrite=False):
		self.time0 = time.time()
		self.runsdir = Correct_Path(runsdir)
		self.config = Correct_Path(config)
		self.configInput = Correct_Path(configinput) if configinput != '' else ''
		self.overwrite = overwrite
		self.are_cal_runs = False if self.configInput == '' else True
		runstemp = glob.glob('{d}/*'.format(d=self.runsdir))
		if len(runstemp) < 1:
			ExitMessage('The directory does not have any runs', os.EX_USAGE)
		# self.num_cores = 1
		self.runs = runstemp
		self.num_runs = len(self.runs)
		if self.num_runs < 1: ExitMessage('There is not even the required data to convert one run', os.EX_DATAERR)
		self.voltages = []
		self.signalOut = {}
		self.signalOutSigma = {}
		self.signalIn = {}
		self.signalInSigma = {}
		self.caen_ch = 3
		self.graph = ro.TGraphErrors()
		self.canvas = None
		self.fit = None
		self.cal_pickle = None

	def DoAll(self):
		self.time0 = time.time()
		self.LoopRuns()
		self.PlotSignals()
		if self.are_cal_runs:
			self.FitLine()
			self.FillPickle()
			self.SavePickle()
		self.SaveCanvas()

	def LoopRuns(self):
		for run in self.runs:
			if os.path.isdir(run):
				configfile = self.config if not self.are_cal_runs or '_out_' in run else self.configInput
				ro.gROOT.SetBatch(1)
				print '\nAnalysing in batch mode run:', run
				print 'Using config file:', configfile
				print 'Overwriting analysis tree if it exists' if self.overwrite else 'Not overwriting analysis tree if it exists'
				anaRun = AnalysisCaenCCD(run, configfile, overw=self.overwrite)
				anaRun.AnalysisWaves()
				self.voltages.append(anaRun.bias)
				self.caen_ch = anaRun.ch_caen_signal if self.caen_ch != anaRun.ch_caen_signal else self.caen_ch
				anaRun.PlotPedestal('Pedestal', cuts=anaRun.cut0.GetTitle())
				anaRun.PlotWaveforms('SelectedWaveforms', 'signal', cuts=anaRun.cut0.GetTitle())
				anaRun.canvas['SelectedWaveforms'].SetLogz()
				anaRun.PlotWaveforms('SelectedWaveformsPedCor', 'signal_ped_corrected', cuts=anaRun.cut0.GetTitle())
				anaRun.canvas['SelectedWaveformsPedCor'].SetLogz()
				anaRun.PlotSignal('PH', cuts=anaRun.cut0.GetTitle())
				if not anaRun.is_cal_run:
					anaRun.FitLanGaus('PH')
					anaRun.PlotHVCurrents('HVCurrents', '', 5)
				else:
					anaRun.FitConvolutedGaussians('PH')
				signalRun = np.double(anaRun.histo['PH'].GetMean())
				signalSigmaRun = np.sqrt(np.power(anaRun.histo['PH'].GetRMS(), 2, dtype='f8') + np.power(anaRun.pedestal_sigma, 2, dtype='f8'), dtype='f8')
				if not self.are_cal_runs or '_out_' in run:
					self.signalOut[anaRun.bias] = signalRun if anaRun.bias < 0 else -signalRun
					self.signalOutSigma[anaRun.bias] = signalSigmaRun
				else:
					self.signalIn[anaRun.bias] = -signalRun if anaRun.bias < 0 else signalRun
					self.signalInSigma[anaRun.bias] = signalSigmaRun
		self.voltages = sorted(set(self.voltages))

	def PlotSignals(self):
		if len(self.voltages) > 0:
			if self.are_cal_runs:
				voltages2 = []
				for volt in self.voltages:
					if volt in self.signalOut.keys() and volt in self.signalIn.keys():
						voltages2.append(volt)
				self.voltages = voltages2
				stepsIn = np.array([self.signalIn[volt] for volt in self.voltages], dtype='f8')
				stepsInErrs = np.array([self.signalInSigma[volt] for volt in self.voltages], dtype='f8')
				signalO = np.array([self.signalOut[volt] for volt in self.voltages], dtype='f8')
				signalOutErrs = np.array([self.signalOutSigma[volt] for volt in self.voltages], dtype='f8')
				self.graph = ro.TGraphErrors(len(self.voltages), stepsIn, signalO, stepsInErrs, signalOutErrs)
				self.graph.SetNameTitle('Signal_vs_CalStep_ch_' + str(self.caen_ch), 'Signal_vs_CalStep_ch_' + str(self.caen_ch))
				self.graph.GetXaxis().SetTitle('vcal step [mV]')
				self.graph.GetYaxis().SetTitle('signal [mV]')
			else:
				signalO = np.array([self.signalOut[volt] for volt in self.voltages], dtype='f8')
				signalOutErrs = np.array([self.signalOutSigma[volt] for volt in self.voltages], dtype='f8')
				self.graph = ro.TGraphErrors(len(self.voltages), np.array(self.voltages, 'f8'), signalO, np.zeros(len(self.voltages), 'f8'), signalOutErrs)
				self.graph.SetNameTitle('Signal_vs_HV', 'Signal_vs_HV')
				self.graph.GetXaxis().SetTitle('HV [V]')
				self.graph.GetYaxis().SetTitle('signal [mV]')
			self.graph.SetMarkerStyle(7)
			self.graph.SetMarkerColor(ro.kBlack)
			self.graph.SetLineColor(ro.kBlack)

			ro.gROOT.SetBatch(0)

			self.canvas = ro.TCanvas('c_' + self.graph.GetName(), 'c_' + self.graph.GetName(), 1)
			self.graph.Draw('AP')

	def FitLine(self):
		ro.Math.MinimizerOptions.SetDefaultMinimizer(*fit_method)
		ro.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(1000000)
		ro.Math.MinimizerOptions.SetDefaultTolerance(0.00001)
		ro.gStyle.SetOptFit(1111)
		func = ro.TF1('fit_' + self.graph.GetName(), 'pol1', -1000, 1000)
		func.SetLineColor(ro.kRed)
		func.SetNpx(int(10 * (2 ** 14 - 1)))
		self.graph.Fit('fit_' + self.graph.GetName(), 'QM0', '', -1000, 1000)
		if func.GetProb() < 0.9:
			self.graph.Fit('fit_' + self.graph.GetName(), 'QM0', '', -1000, 1000)
		if func.GetProb() < 0.9:
			self.graph.Fit('fit_' + self.graph.GetName(), 'QM0', '', -1000, 1000)
		self.fit = func
		self.fit.Draw('same')

	def FillPickle(self):
		self.cal_pickle = {'voltages': self.voltages,
						   'signal_in': self.signalIn,
						   'signal_in_sigma': self.signalInSigma,
						   'signal_out': self.signalOut,
						   'signal_out_sigma': self.signalOutSigma,
						   'caen_ch': self.caen_ch,
						   'fit_p0': self.fit.GetParameter(0),
						   'fit_p0_error': self.fit.GetParError(0),
						   'fit_p1': self.fit.GetParameter(1),
						   'fit_p1_error': self.fit.GetParError(1),
						   'fit_prob': self.fit.GetProb(),
						   'fit_chi2': self.fit.GetChisquare(),
						   'fit_ndf': self.fit.GetNDF()
						   }

	def SavePickle(self):
		pickleName = 'signal_cal_{c}.cal'.format(c=self.caen_ch)
		if not self.cal_pickle:
			self.FillPickle()
		if os.path.isfile('{d}/{f}'.format(d=self.runsdir, f=pickleName)):
			if not self.overwrite:
				print 'The file', pickleName, 'already exists in', self.runsdir
				return
		with open('{d}/{f}'.format(d=self.runsdir, f=pickleName), 'wb') as fpickle:
			pickle.dump(self.cal_pickle, fpickle, pickle.HIGHEST_PROTOCOL)
		print 'Saved calibration pickle', pickleName, 'in', self.runsdir

	def SaveCanvas(self):
		self.canvas.SaveAs('{d}/{n}.png'.format(d=self.runsdir, n=self.graph.GetName()))
		self.canvas.SaveAs('{d}/{n}.root'.format(d=self.runsdir, n=self.graph.GetName()))

def main():
	parser = OptionParser()
	parser.add_option('-d', '--runsdir', dest='runsdir', type='str', default='', help='path to folder containing all the run folders to modify')
	parser.add_option('-c', '--config', dest='config', type='str', default='', help='path to analysis config file used for all runs inside the folder')
	parser.add_option('-o', '--overwrite', dest='overwrite', default=False, action='store_true', help='Sets overwrite of analysis tree for the run if it exists')
	parser.add_option('--configinput', dest='configinput', type='str', default='', help='path to analysis config file used for all runs inside the folder that are of calibration type "in". Only necessary when doing signal calibration')
	(options, args) = parser.parse_args()
	runsdir = str(options.runsdir)
	config = str(options.config)
	configinput = str(options.configinput)
	overwrite = bool(options.overwrite)
	arif = AnalysisAllRunsInFolder(runsdir, config, configinput, overwrite)
	return arif

if __name__ == '__main__':
	arif = main()
