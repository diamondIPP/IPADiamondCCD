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
		self.are_cal_runs = False if self.configInput == '' else True # this is only true for signal calibration runs
		runstemp = glob.glob('{d}/*'.format(d=self.runsdir))
		if len(runstemp) < 1:
			ExitMessage('The directory does not have any runs', os.EX_USAGE)
		# self.num_cores = 1
		self.runs = [runi for runi in runstemp if os.path.isdir(runi)]
		if self.are_cal_runs:
			self.runs.sort(key=lambda x: float(x.split('_')[-1].split('mV')[0].split('V')[0]))
		else:
			self.runs.sort(key=lambda x: float(x.split('_Pos_')[-1].split('V')[0]) if 'Pos' in x else -1*float(x.split('_Neg_')[-1].split('V')[0]))
		self.num_runs = len(self.runs)
		if self.num_runs < 1: ExitMessage('There is not even the required data to convert one run', os.EX_DATAERR)
		self.voltages = []
		self.diaVoltages = {}
		self.diaVoltagesSigma = {}
		self.signalOut = {}
		self.signalOutSigma = {}
		self.signalOutVcal = {}
		self.signalOutVcalSigma = {}
		self.signalOutCharge = {}
		self.signalOutChargeSigma = {}
		self.signalIn = {}
		self.signalInSigma = {}
		self.caen_ch = 3
		self.graph = ro.TGraphErrors()
		self.graphVcal = ro.TGraphErrors()
		self.graphCharge = ro.TGraphErrors()
		self.canvas = None
		self.canvasVcal = None
		self.canvasCharge = None
		self.fit = None
		self.cal_pickle = None
		self.LoadPickle()

	def LoadPickle(self):
		cal_files = glob.glob('{d}/signal_*.cal'.format(d=self.runsdir))
		if len(cal_files) > 0:
			self.cal_pickle = pickle.load(open(cal_files[0], 'rb'))
			self.voltages = self.cal_pickle['voltages']
			self.signalIn = self.cal_pickle['signal_in']
			self.signalInSigma = self.cal_pickle['signal_in_sigma']
			self.signalOut = self.cal_pickle['signal_out']
			self.signalOutSigma = self.cal_pickle['signal_out_sigma']
			self.caen_ch = self.cal_pickle['caen_ch']
			nameVcal = ''
			nameCharge = ''
			if 'signal_out_vcal' in self.cal_pickle.keys():
				self.signalOutVcal = self.cal_pickle['signal_out_vcal']
				self.signalOutVcalSigma = self.cal_pickle['signal_out_vcal_sigma']
				nameVcal = '' if self.are_cal_runs else 'SignalVcal_vs_HV'
			if 'signal_out_charge' in self.cal_pickle.keys():
				self.signalOutCharge = self.cal_pickle['signal_out_charge']
				self.signalOutChargeSigma = self.cal_pickle['signal_out_charge_sigma']
				nameCharge = '' if self.are_cal_runs else 'SignalCharge_vs_HV'
			name = 'Signal_vs_CalStep_ch_' + str(self.caen_ch) if self.are_cal_runs else 'Signal_vs_HV'
			xpoints = np.array([self.signalIn[volt] for volt in self.voltages], 'f8') if len(self.signalIn.keys()) >= 1 else np.array(self.voltages, 'f8')
			xpointserrs = np.array([self.signalInSigma[volt] for volt in self.voltages], 'f8') if len(self.signalInSigma.keys()) >= 1 else np.zeros(len(self.voltages), 'f8')
			self.graph = ro.TGraphErrors(len(self.voltages), xpoints, np.array([self.signalOut[volt] for volt in self.voltages], 'f8'), xpointserrs, np.array([self.signalOutSigma[volt] for volt in self.voltages], 'f8'))
			self.graph.SetNameTitle(name, name)
			self.graph.GetXaxis().SetTitle('vcal step [mV]' if self.are_cal_runs else 'HV [V]')
			self.graph.GetYaxis().SetTitle('signal [mV]')
			self.graph.SetMarkerStyle(7)
			self.graph.SetMarkerColor(ro.kBlack)
			self.graph.SetLineColor(ro.kBlack)
			if self.are_cal_runs:
				self.fit = ro.TF1('fit_' + name, 'pol1', -1000, 1000)
				self.fit.SetLineColor(ro.kRed)
				self.fit.SetNpx(10000)
				self.fit.SetParameters(np.array([self.cal_pickle['fit_p0'], self.cal_pickle['fit_p1']], 'f8'))
				self.fit.SetParErrors(np.array([self.cal_pickle['fit_p0_error'], self.cal_pickle['fit_p1_error']], 'f8'))
				self.fit.SetNDF(self.cal_pickle['fit_ndf'])
				self.fit.SetChisquare(self.cal_pickle['fit_chi2'])
			else:
				if nameVcal != '':
					self.graphVcal = ro.TGraphErrors(len(self.voltages), xpoints, np.array([self.signalOutVcal[volt] for volt in self.voltages], 'f8'), xpointserrs, np.array([self.signalOutVcalSigma[volt] for volt in self.voltages], 'f8'))
					self.graphVcal.SetNameTitle(nameVcal, nameVcal)
					self.graphVcal.GetXaxis().SetTitle('HV [V]')
					self.graphVcal.GetYaxis().SetTitle('signalVcal [mV]')
					self.graphVcal.SetMarkerStyle(7)
					self.graphVcal.SetMarkerColor(ro.kBlack)
					self.graphVcal.SetLineColor(ro.kBlack)
				if nameCharge != '':
					self.graphCharge = ro.TGraphErrors(len(self.voltages), xpoints, np.array([self.signalOutCharge[volt] for volt in self.voltages], 'f8'), xpointserrs, np.array([self.signalOutChargeSigma[volt] for volt in self.voltages], 'f8'))
					self.graphCharge.SetNameTitle(nameVcal, nameVcal)
					self.graphCharge.GetXaxis().SetTitle('HV [V]')
					self.graphCharge.GetYaxis().SetTitle('signalCharge [e]')
					self.graphCharge.SetMarkerStyle(7)
					self.graphCharge.SetMarkerColor(ro.kBlack)
					self.graphCharge.SetLineColor(ro.kBlack)
				self.fit = None
			print 'Loaded pickle', cal_files[0]
			return
		print 'There is no pickle to load yet (or it is not a calibration run)'

	def PlotFromPickle(self):
		if self.cal_pickle:
			if self.graph:
				self.canvas = ro.TCanvas('c_' + self.graph.GetName(), 'c_' + self.graph.GetName(), 1)
				self.graph.Draw('AP')
			if self.fit:
				self.fit.Draw('same')
			if self.graphVcal:
				self.canvasVcal = ro.TCanvas('c_' + self.graphVcal.GetName(), 'c_' + self.graphVcal.GetName(), 1)
				self.graphVcal.Draw('AP')
			if self.graphCharge:
				self.canvasCharge = ro.TCanvas('c_' + self.graphCharge.GetName(), 'c_' + self.graphCharge.GetName(), 1)
				self.graphCharge.Draw('AP')

	def DoAll(self, updateHistos=False):
		self.time0 = time.time()
		if not self.cal_pickle or self.overwrite or updateHistos:
			self.LoopRuns()
		self.PlotSignals()
		if self.are_cal_runs:
			self.FitLine()
			self.FillPickle()
			self.SavePickle()
		self.SaveCanvas()
		print 'Finished in', time.time() - self.time0, 'seconds'

	def LoopRuns(self):
		# ro.gROOT.SetBatch(True)
		for run in self.runs:
			print 'Analysing run:', run
			if os.path.isdir(run):
				configfile = self.config if not self.are_cal_runs or '_out_' in run else self.configInput
				print 'Using config file:', configfile
				print 'Overwriting analysis tree if it exists' if self.overwrite else 'Not overwriting analysis tree if it exists'
				anaRun = AnalysisCaenCCD(run, configfile, overw=self.overwrite)
				anaRun.AnalysisWaves()
				anaRun.AddVcalFriend()
				anaRun.AddChargeFriend()
				self.voltages.append(anaRun.bias)
				self.caen_ch = anaRun.ch_caen_signal if self.caen_ch != anaRun.ch_caen_signal else self.caen_ch
				anaRun.PlotPedestal('Pedestal', cuts=anaRun.cut0.GetTitle(), branch='pedestal')
				if not anaRun.is_cal_run:
					anaRun.PlotPedestal('PedestalVcal', cuts=anaRun.cut0.GetTitle(), branch='pedestalVcal')
					anaRun.PlotPedestal('PedestalCharge', cuts=anaRun.cut0.GetTitle(), branch='pedestalCharge')
				anaRun.PlotWaveforms('SelectedWaveforms', 'signal', cuts=anaRun.cut0.GetTitle())
				if 'SelectedWaveforms' in anaRun.canvas.keys():
					anaRun.canvas['SelectedWaveforms'].SetLogz()
				anaRun.PlotWaveforms('SelectedWaveformsPedCor', 'signal_ped_corrected', cuts=anaRun.cut0.GetTitle())
				if 'SelectedWaveformsPedCor' in anaRun.canvas.keys():
					anaRun.canvas['SelectedWaveformsPedCor'].SetLogz()
				anaRun.PlotSignal('PH', cuts=anaRun.cut0.GetTitle())
				if not anaRun.is_cal_run:
					anaRun.FitLanGaus('PH')
					anaRun.PlotSignal('PHvcal', cuts=anaRun.cut0.GetTitle(), branch='signalVcal')
					anaRun.FitLanGaus('PHvcal')
					anaRun.PlotSignal('PHcharge', cuts=anaRun.cut0.GetTitle(), branch='signalCharge')
					anaRun.FitLanGaus('PHcharge')
					anaRun.PlotHVCurrents('HVCurrents', '', 5)
					anaRun.PlotDiaVoltage('DUTVoltage', '', 5)
				else:
					anaRun.FitConvolutedGaussians('PH')
				anaRun.SaveAllCanvas()

				signalRun = 0
				signalSigmaRun = 0
				if 'PH' in anaRun.histo.keys():
					signalRun = np.double(anaRun.histo['PH'].GetMean())
					signalSigmaRun = np.sqrt(np.power(anaRun.histo['PH'].GetRMS(), 2, dtype='f8') - np.power(anaRun.pedestal_sigma, 2, dtype='f8'), dtype='f8') if anaRun.histo['PH'].GetRMS() > anaRun.pedestal_sigma else anaRun.histo['PH'].GetRMS()
				self.diaVoltages[anaRun.bias] = anaRun.voltageDiaMean
				self.diaVoltagesSigma[anaRun.bias] = anaRun.voltageDiaSpread
				if not anaRun.is_cal_run:
					signalRunVcal = 0
					signalSigmaRunVcal = 0
					signalRunCharge = 0
					signalSigmaRunCharge = 0
					if 'PHvcal' in anaRun.histo.keys():
						signalRunVcal = np.double(anaRun.histo['PHvcal'].GetMean())
						signalSigmaRunVcal = np.sqrt(np.power(anaRun.histo['PHvcal'].GetRMS(), 2, dtype='f8') - np.power(anaRun.pedestal_vcal_sigma, 2, dtype='f8'), dtype='f8') if anaRun.histo['PHvcal'].GetRMS() > anaRun.pedestal_vcal_sigma else anaRun.histo['PHvcal'].GetRMS()
					if 'PHcharge' in anaRun.histo.keys():
						signalRunCharge = np.double(anaRun.histo['PHcharge'].GetMean())
						signalSigmaRunCharge = np.sqrt(np.power(anaRun.histo['PHcharge'].GetRMS(), 2, dtype='f8') - np.power(anaRun.pedestal_charge_sigma, 2, dtype='f8'), dtype='f8')  if anaRun.histo['PHcharge'].GetRMS() > anaRun.pedestal_charge_sigma else anaRun.histo['PHcharge'].GetRMS()
				if not self.are_cal_runs or '_out_' in run:
					self.signalOut[anaRun.bias] = signalRun if anaRun.bias < 0 else -signalRun
					self.signalOutSigma[anaRun.bias] = signalSigmaRun
					if not self.are_cal_runs:
						self.signalOutVcal[anaRun.bias] = -signalRunVcal if anaRun.bias < 0 else signalRunVcal
						self.signalOutVcalSigma[anaRun.bias] = signalSigmaRunVcal
						self.signalOutCharge[anaRun.bias] = -signalRunCharge if anaRun.bias < 0 else signalRunCharge
						self.signalOutChargeSigma[anaRun.bias] = signalSigmaRunCharge
				else:
					self.signalIn[anaRun.bias] = -signalRun if anaRun.bias < 0 else signalRun
					self.signalInSigma[anaRun.bias] = signalSigmaRun
				del anaRun
		self.voltages = sorted(set(self.voltages))
		# ro.gROOT.SetBatch(False)

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
				signalOVcal = np.array([self.signalOutVcal[volt] for volt in self.voltages], dtype='f8')
				signalOutVcalErrs = np.array([self.signalOutVcalSigma[volt] for volt in self.voltages], dtype='f8')
				signalOCharge = np.array([self.signalOutCharge[volt] for volt in self.voltages], dtype='f8')
				signalOutChargeErrs = np.array([self.signalOutChargeSigma[volt] for volt in self.voltages], dtype='f8')
				diaVoltages = np.array([self.diaVoltages[volt] for volt in self.voltages], dtype='f8')
				diaVoltagesSigma = np.array([self.diaVoltagesSigma[volt] for volt in self.voltages], dtype='f8')
				self.graph = ro.TGraphErrors(len(self.voltages), diaVoltages, signalO, diaVoltagesSigma, signalOutErrs)
				self.graph.SetNameTitle('Signal_vs_HV', 'Signal_vs_HV')
				self.graph.GetXaxis().SetTitle('HV [V]')
				self.graph.GetYaxis().SetTitle('signal [mV]')

				self.graphVcal = ro.TGraphErrors(len(self.voltages), diaVoltages, signalOVcal, diaVoltagesSigma, signalOutVcalErrs)
				self.graphVcal.SetNameTitle('SignalVcal_vs_HV', 'SignalVcal_vs_HV')
				self.graphVcal.GetXaxis().SetTitle('HV [V]')
				self.graphVcal.GetYaxis().SetTitle('signalVcal [mV]')

				self.graphCharge = ro.TGraphErrors(len(self.voltages), diaVoltages, signalOCharge, diaVoltagesSigma, signalOutChargeErrs)
				self.graphCharge.SetNameTitle('SignalCharge_vs_HV', 'SignalCharge_vs_HV')
				self.graphCharge.GetXaxis().SetTitle('HV [V]')
				self.graphCharge.GetYaxis().SetTitle('signalCharge [e]')

			self.graph.SetMarkerStyle(7)
			self.graph.SetMarkerColor(ro.kBlack)
			self.graph.SetLineColor(ro.kBlack)

			self.canvas = ro.TCanvas('c_' + self.graph.GetName(), 'c_' + self.graph.GetName(), 1)
			self.graph.Draw('AP')
			self.canvas.SetGridx()
			self.canvas.SetGridy()

			if not self.are_cal_runs:
				self.graphVcal.SetMarkerStyle(7)
				self.graphVcal.SetMarkerColor(ro.kBlack)
				self.graphVcal.SetLineColor(ro.kBlack)

				self.canvasVcal = ro.TCanvas('c_' + self.graphVcal.GetName(), 'c_' + self.graphVcal.GetName(), 1)
				self.graphVcal.Draw('AP')
				self.canvasVcal.SetGridx()
				self.canvasVcal.SetGridy()

				self.graphCharge.SetMarkerStyle(7)
				self.graphCharge.SetMarkerColor(ro.kBlack)
				self.graphCharge.SetLineColor(ro.kBlack)

				self.canvasCharge = ro.TCanvas('c_' + self.graphCharge.GetName(), 'c_' + self.graphCharge.GetName(), 1)
				self.graphCharge.Draw('AP')
				self.canvasCharge.SetGridx()
				self.canvasCharge.SetGridy()

	def FitLine(self):
		ro.Math.MinimizerOptions.SetDefaultMinimizer(*fit_method)
		ro.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(1000000)
		ro.Math.MinimizerOptions.SetDefaultTolerance(0.00001)
		ro.gStyle.SetOptFit(1111)
		func = ro.TF1('fit_' + self.graph.GetName(), 'pol1', -1000, 1000)
		func.SetLineColor(ro.kRed)
		func.SetNpx(10000)
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
						   'signal_out_vcal': self.signalOutVcal,
						   'signal_out_vcal_sigma': self.signalOutVcalSigma,
						   'signal_out_charge': self.signalOutCharge,
						   'signal_out_charge_sigma': self.signalOutChargeSigma,
						   'caen_ch': self.caen_ch,
						   'fit_p0': self.fit.GetParameter(0) if self.fit else 0,
						   'fit_p0_error': self.fit.GetParError(0) if self.fit else 0,
						   'fit_p1': self.fit.GetParameter(1) if self.fit else 0,
						   'fit_p1_error': self.fit.GetParError(1) if self.fit else 0,
						   'fit_prob': self.fit.GetProb() if self.fit else 0,
						   'fit_chi2': self.fit.GetChisquare() if self.fit else 0,
						   'fit_ndf': self.fit.GetNDF() if self.fit else 0
						   }

	def SavePickle(self):
		pickleName = 'signal_cal_{c}.cal'.format(c=self.caen_ch) if self.are_cal_runs else 'signal_{c}.cal'.format(c=self.caen_ch)
		if not self.cal_pickle:
			self.FillPickle()
		if os.path.isfile('{d}/{f}'.format(d=self.runsdir, f=pickleName)):
			if not self.overwrite:
				print 'The file', pickleName, 'already exists in', self.runsdir
				return
		with open('{d}/{f}'.format(d=self.runsdir, f=pickleName), 'wb') as fpickle:
			pickle.dump(self.cal_pickle, fpickle, pickle.HIGHEST_PROTOCOL)
		print 'Saved pickle', pickleName, 'in', self.runsdir

	def SaveCanvas(self):
		self.canvas.SaveAs('{d}/{n}.png'.format(d=self.runsdir, n=self.graph.GetName()))
		self.canvas.SaveAs('{d}/{n}.root'.format(d=self.runsdir, n=self.graph.GetName()))
		if self.canvasVcal:
			self.canvasVcal.SaveAs('{d}/{n}.png'.format(d=self.runsdir, n=self.graphVcal.GetName()))
			self.canvasVcal.SaveAs('{d}/{n}.root'.format(d=self.runsdir, n=self.graphVcal.GetName()))
		if self.canvasCharge:
			self.canvasCharge.SaveAs('{d}/{n}.png'.format(d=self.runsdir, n=self.graphCharge.GetName()))
			self.canvasCharge.SaveAs('{d}/{n}.root'.format(d=self.runsdir, n=self.graphCharge.GetName()))


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
