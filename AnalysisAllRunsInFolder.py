#!/usr/bin/env python
import numpy as np
from struct import unpack
import time, os, sys
from optparse import OptionParser
import ipdb
import pickle as pickle
from Utils import *
from Settings_Caen import Settings_Caen
from Converter_Caen import Converter_Caen
import subprocess as subp
from AnalysisCaenCCD import AnalysisCaenCCD
import ROOT as ro

fit_method = ('Minuit2', 'Migrad', )

class AnalysisAllRunsInFolder:
	def __init__(self, runsdir='', config='', configinput='', overwrite=False, doDebug=False):
		self.time0 = time.time()
		self.runsdir = Correct_Path(runsdir)
		self.config = Correct_Path(config)
		self.configInput = Correct_Path(configinput) if configinput != '' else ''
		self.overwrite = overwrite
		self.doDebug = doDebug
		self.are_cal_runs = False if self.configInput == '' else True # this is only true for signal calibration runs
		runstemp = glob.glob('{d}/*'.format(d=self.runsdir))
		if len(runstemp) < 1:
			ExitMessage('The directory does not have any runs', os.EX_USAGE)
		# self.num_cores = 1
		self.runs = [runi for runi in runstemp if os.path.isdir(runi)]
		if self.are_cal_runs:
			self.runs.sort(key=lambda x: float(x.split('_')[-1].split('mV')[0].split('V')[0]))
		else:
			self.runs.sort(key=lambda x: float(x.split('_Pos_')[-1].split('V')[0]) if '_Pos_' in x else -1*float(x.split('_Neg_')[-1].split('V')[0]))
		self.num_runs = len(self.runs)
		if self.num_runs < 1: ExitMessage('There is not even the required data to convert one run', os.EX_DATAERR)

		self.caen_ch = 3

		self.voltages = []

		self.diaVoltages = {}
		self.diaVoltagesSigma = {}

		self.signalOut = {}
		self.signalOutSigma = {}
		self.signalOutVcal = {}
		self.signalOutVcalSigma = {}
		self.signalOutCharge = {}
		self.signalOutChargeSigma = {}

		self.signalOutCF = {}
		self.signalOutSigmaCF = {}
		self.signalOutVcalCF = {}
		self.signalOutVcalSigmaCF = {}
		self.signalOutChargeCF = {}
		self.signalOutChargeSigmaCF = {}

		self.signalOutNP = {}
		self.signalOutSigmaNP = {}
		self.signalOutVcalNP = {}
		self.signalOutVcalSigmaNP = {}
		self.signalOutChargeNP = {}
		self.signalOutChargeSigmaNP = {}

		self.signalOutCFNP = {}
		self.signalOutSigmaCFNP = {}
		self.signalOutVcalCFNP = {}
		self.signalOutVcalSigmaCFNP = {}
		self.signalOutChargeCFNP = {}
		self.signalOutChargeSigmaCFNP = {}

		self.signalIn = {}
		self.signalInSigma = {}

		self.signalPeakTime = {}
		self.signalPeakTimeSigma = {}
		self.signalPeakTimeCF = {}
		self.signalPeakTimeSigmaCF = {}

		self.graphPeakTime = ro.TGraphErrors()
		self.graphPeakTimeCF = ro.TGraphErrors()

		self.graphPH = ro.TGraphErrors()
		self.graphVcal = ro.TGraphErrors()
		self.graphCharge = ro.TGraphErrors()
		self.graphPHCF = ro.TGraphErrors()
		self.graphVcalCF = ro.TGraphErrors()
		self.graphChargeCF = ro.TGraphErrors()

		self.graphPHNP = ro.TGraphErrors()
		self.graphVcalNP = ro.TGraphErrors()
		self.graphChargeNP = ro.TGraphErrors()
		self.graphPHCFNP = ro.TGraphErrors()
		self.graphVcalCFNP = ro.TGraphErrors()
		self.graphChargeCFNP = ro.TGraphErrors()

		self.canvasPeakTime = None
		self.canvasPeakTimeCF = None

		self.canvasPH = None
		self.canvasVcal = None
		self.canvasCharge = None
		self.canvasPHCF = None
		self.canvasVcalCF = None
		self.canvasChargeCF = None

		self.canvasPHNP = None
		self.canvasVcalNP = None
		self.canvasChargeNP = None
		self.canvasPHCFNP = None
		self.canvasVcalCFNP = None
		self.canvasChargeCFNP = None

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
			namePeakTime = ''
			if 'signal_out_vcal' in list(self.cal_pickle.keys()):
				self.signalOutVcal = self.cal_pickle['signal_out_vcal']
				self.signalOutVcalSigma = self.cal_pickle['signal_out_vcal_sigma']
				nameVcal = '' if self.are_cal_runs else 'SignalVcal_vs_HV'
			if 'signal_out_charge' in list(self.cal_pickle.keys()):
				self.signalOutCharge = self.cal_pickle['signal_out_charge']
				self.signalOutChargeSigma = self.cal_pickle['signal_out_charge_sigma']
				nameCharge = '' if self.are_cal_runs else 'SignalCharge_vs_HV'
			if 'signal_peak_time' in list(self.cal_pickle.keys()):
				self.signalPeakTime = self.cal_pickle['signal_peak_time']
				self.signalPeakTimeSigma = self.cal_pickle['signal_peak_time_sigma']
				namePeakTime = 'SignalPeakTime_vs_Signal_ch_' + str(self.caen_ch) if self.are_cal_runs else ''
			name = 'Signal_vs_CalStep_ch_' + str(self.caen_ch) if self.are_cal_runs else 'Signal_vs_HV'
			xpoints = np.array([self.signalIn[volt] for volt in self.voltages], 'f8') if len(list(self.signalIn.keys())) >= 1 else np.array(self.voltages, 'f8')
			xpointserrs = np.array([self.signalInSigma[volt] for volt in self.voltages], 'f8') if len(list(self.signalInSigma.keys())) >= 1 else np.zeros(len(self.voltages), 'f8')
			self.graphPH = ro.TGraphErrors(len(self.voltages), xpoints, np.array([self.signalOut[volt] for volt in self.voltages], 'f8'), xpointserrs, np.array([self.signalOutSigma[volt] for volt in self.voltages], 'f8'))
			self.graphPH.SetNameTitle(name, name)
			self.graphPH.GetXaxis().SetTitle('vcal step [mV]' if self.are_cal_runs else 'HV [V]')
			self.graphPH.GetYaxis().SetTitle('signal [mV]')
			self.graphPH.SetMarkerStyle(7)
			self.graphPH.SetMarkerColor(ro.kBlack)
			self.graphPH.SetLineColor(ro.kBlack)
			if self.are_cal_runs:
				self.fit = ro.TF1('fit_' + name, 'pol1', -1000, 1000)
				self.fit.SetLineColor(ro.kRed)
				self.fit.SetNpx(10000)
				self.fit.SetParameters(np.array([self.cal_pickle['fit_p0'], self.cal_pickle['fit_p1']], 'f8'))
				self.fit.SetParErrors(np.array([self.cal_pickle['fit_p0_error'], self.cal_pickle['fit_p1_error']], 'f8'))
				self.fit.SetNDF(self.cal_pickle['fit_ndf'])
				self.fit.SetChisquare(self.cal_pickle['fit_chi2'])
				if namePeakTime != '':
					self.graphPeakTime = ro.TGraphErrors(len(list(self.signalPeakTime.values())), np.array([self.signalOut[volt] for volt in self.voltages], 'f8'), np.array([self.signalPeakTime[volt] for volt in self.voltages], 'f8'), np.array([self.signalOutSigma[volt] for volt in self.voltages], 'f8'), np.array([self.signalPeakTimeSigma[volt] for volt in self.voltages], 'f8'))
					self.graphPeakTime.SetNameTitle(namePeakTime, namePeakTime)
					self.graphPeakTime.GetXaxis().SetTitle('signal [mV]')
					self.graphPeakTime.GetYaxis().SetTitle('Peak Time [us]')
					self.graphPeakTime.SetMarkerStyle(7)
					self.graphPeakTime.SetMarkerColor(ro.kBlack)
					self.graphPeakTime.SetLineColor(ro.kBlack)
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
			print('Loaded pickle', cal_files[0])
			return
		print('There is no pickle to load yet (or it is not a calibration run)')

	def PlotFromPickle(self):
		if self.cal_pickle:
			if self.graphPH:
				self.canvasPH = ro.TCanvas('c_' + self.graphPH.GetName(), 'c_' + self.graphPH.GetName(), 1)
				self.graphPH.Draw('AP')
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
			self.PlotPeakTimes()
		self.SaveAllCanvas()
		print('Finished in', time.time() - self.time0, 'seconds')

	def LoopRuns(self):
		# ro.gROOT.SetBatch(True)
		for run in self.runs:
			print('\nAnalysing run:', run)
			if os.path.isdir(run):
				root_files = glob.glob('{d}/*.root'.format(d=run))
				if len(root_files) > 0:
					configfile = self.config if not self.are_cal_runs or '_out_' in run else self.configInput
					print('Using config file:', configfile)
					print('Overwriting analysis tree if it exists' if self.overwrite else 'Not overwriting analysis tree if it exists')
					anaRun = AnalysisCaenCCD(run, configfile, overw=self.overwrite, doDebug=self.doDebug)
					anaRun.DoAll()
					self.voltages.append(anaRun.bias)
					self.caen_ch = anaRun.ch_caen_signal if self.caen_ch != anaRun.ch_caen_signal else self.caen_ch

					peakTime = 0
					peakTimeSigma = 0
					peakTimeCF = 0
					peakTimeSigmaCF = 0
					if self.are_cal_runs and '_out_' in run:
						if 'peakPosDist' in list(anaRun.histo.keys()):
							peakTime = np.double(anaRun.peakTime * 1e6)
							peakTimeSigma = np.double(anaRun.histo['peakPosDist'].GetRMS())
						if 'peakPosDistCF' in list(anaRun.histo.keys()):
							peakTimeCF = np.double(anaRun.peakTimeCF * 1e6)
							peakTimeSigmaCF = np.double(anaRun.histo['peakPosDistCF'].GetRMS())

					self.diaVoltages[anaRun.bias] = anaRun.voltageDiaMean
					self.diaVoltagesSigma[anaRun.bias] = anaRun.voltageDiaSpread

					signalRun = 0
					signalSigmaRun = 0
					signalRunCF = 0
					signalSigmaRunCF = 0

					signalRunNP = 0
					signalSigmaRunNP = 0
					signalRunCFNP = 0
					signalSigmaRunCFNP = 0

					if 'PH' in list(anaRun.histo.keys()):
						signalRun = np.double(anaRun.langaus['PH'].GetParameter(1)) if self.are_cal_runs else np.double(anaRun.histo['PH'].GetMean())
						signalSigmaRun = np.double(anaRun.langaus['PH'] .GetParameter(2)) if self.are_cal_runs else np.sqrt(np.power(anaRun.histo['PH'].GetRMS(), 2, dtype='f8') - np.power(anaRun.pedestal_sigma, 2, dtype='f8'), dtype='f8') if anaRun.histo['PH'].GetRMS() > anaRun.pedestal_sigma else anaRun.histo['PH'].GetRMS()
						signalRunCF = np.double(anaRun.langaus['PH_CF'].GetParameter(1)) if self.are_cal_runs else np.double(anaRun.histo['PH_CF'].GetMean())
						signalSigmaRunCF = np.double(anaRun.langaus['PH_CF'] .GetParameter(2)) if self.are_cal_runs else np.sqrt(np.power(anaRun.histo['PH_CF'].GetRMS(), 2, dtype='f8') - np.power(anaRun.pedestal_sigma, 2, dtype='f8'), dtype='f8') if anaRun.histo['PH_CF'].GetRMS() > anaRun.pedestal_sigma else anaRun.histo['PH_CF'].GetRMS()

						signalRunNP = signalRun
						signalSigmaRunNP = signalSigmaRun
						signalRunCFNP = signalRunCF
						signalSigmaRunCFNP = signalSigmaRunCF

						if 'PHNoPedestal' in anaRun.langaus or 'PHNoPedestal' in anaRun.histo:
							signalRunNP = np.double(anaRun.langaus['PHNoPedestal'].GetParameter(1)) if self.are_cal_runs else np.double(anaRun.histo['PHNoPedestal'].GetMean())
							signalSigmaRunNP = np.double(anaRun.langaus['PHNoPedestal'] .GetParameter(2)) if self.are_cal_runs else np.sqrt(np.power(anaRun.histo['PHNoPedestal'].GetRMS(), 2, dtype='f8') - np.power(anaRun.pedestal_sigma, 2, dtype='f8'), dtype='f8') if anaRun.histo['PHNoPedestal'].GetRMS() > anaRun.pedestal_sigma else anaRun.histo['PHNoPedestal'].GetRMS()
						if 'PH_CFNoPedestal' in anaRun.langaus or 'PH_CFNoPedestal' in anaRun.histo:
							signalRunCFNP = np.double(anaRun.langaus['PH_CFNoPedestal'].GetParameter(1)) if self.are_cal_runs else np.double(anaRun.histo['PH_CFNoPedestal'].GetMean())
							signalSigmaRunCFNP = np.double(anaRun.langaus['PH_CFNoPedestal'] .GetParameter(2)) if self.are_cal_runs else np.sqrt(np.power(anaRun.histo['PH_CFNoPedestal'].GetRMS(), 2, dtype='f8') - np.power(anaRun.pedestal_sigma, 2, dtype='f8'), dtype='f8') if anaRun.histo['PH_CFNoPedestal'].GetRMS() > anaRun.pedestal_sigma else anaRun.histo['PH_CFNoPedestal'].GetRMS()

					if not anaRun.is_cal_run:
						signalRunVcal = 0
						signalSigmaRunVcal = 0
						signalRunCharge = 0
						signalSigmaRunCharge = 0
						signalRunVcalCF = 0
						signalSigmaRunVcalCF = 0
						signalRunChargeCF = 0
						signalSigmaRunChargeCF = 0

						if 'PHvcal' in list(anaRun.histo.keys()):
							signalRunVcal = np.double(anaRun.histo['PHvcal'].GetMean())
							signalSigmaRunVcal = np.sqrt(np.power(anaRun.histo['PHvcal'].GetRMS(), 2, dtype='f8') - np.power(anaRun.pedestal_vcal_sigma, 2, dtype='f8'), dtype='f8') if anaRun.histo['PHvcal'].GetRMS() > anaRun.pedestal_vcal_sigma else anaRun.histo['PHvcal'].GetRMS()
						if 'PHcharge' in list(anaRun.histo.keys()):
							signalRunCharge = np.double(anaRun.histo['PHcharge'].GetMean())
							signalSigmaRunCharge = np.sqrt(np.power(anaRun.histo['PHcharge'].GetRMS(), 2, dtype='f8') - np.power(anaRun.pedestal_charge_sigma, 2, dtype='f8'), dtype='f8')  if anaRun.histo['PHcharge'].GetRMS() > anaRun.pedestal_charge_sigma else anaRun.histo['PHcharge'].GetRMS()
						if 'PHvcal_CF' in list(anaRun.histo.keys()):
							signalRunVcalCF = np.double(anaRun.histo['PHvcal_CF'].GetMean())
							signalSigmaRunVcalCF = np.sqrt(np.power(anaRun.histo['PHvcal_CF'].GetRMS(), 2, dtype='f8') - np.power(anaRun.pedestal_vcal_sigma, 2, dtype='f8'), dtype='f8') if anaRun.histo['PHvcal_CF'].GetRMS() > anaRun.pedestal_vcal_sigma else anaRun.histo['PHvcal_CF'].GetRMS()
						if 'PHcharge_CF' in list(anaRun.histo.keys()):
							signalRunChargeCF = np.double(anaRun.histo['PHcharge_CF'].GetMean())
							signalSigmaRunChargeCF = np.sqrt(np.power(anaRun.histo['PHcharge_CF'].GetRMS(), 2, dtype='f8') - np.power(anaRun.pedestal_charge_sigma, 2, dtype='f8'), dtype='f8')  if anaRun.histo['PHcharge_CF'].GetRMS() > anaRun.pedestal_charge_sigma else anaRun.histo['PHcharge_CF'].GetRMS()

						signalRunVcalNP = signalRunVcal
						signalSigmaRunVcalNP = signalSigmaRunVcal
						signalRunChargeNP = signalRunCharge
						signalSigmaRunChargeNP = signalSigmaRunCharge
						signalRunVcalCFNP = signalRunVcalCF
						signalSigmaRunVcalCFNP = signalSigmaRunVcalCF
						signalRunChargeCFNP = signalRunChargeCF
						signalSigmaRunChargeCFNP = signalSigmaRunChargeCF

						if 'PHvcalNoPedestal' in list(anaRun.histo.keys()):
							signalRunVcalNP = np.double(anaRun.histo['PHvcalNoPedestal'].GetMean())
							signalSigmaRunVcalNP = np.sqrt(np.power(anaRun.histo['PHvcalNoPedestal'].GetRMS(), 2, dtype='f8') - np.power(anaRun.pedestal_vcal_sigma, 2, dtype='f8'), dtype='f8') if anaRun.histo['PHvcalNoPedestal'].GetRMS() > anaRun.pedestal_vcal_sigma else anaRun.histo['PHvcalNoPedestal'].GetRMS()
						if 'PHchargeNoPedestal' in list(anaRun.histo.keys()):
							signalRunChargeNP = np.double(anaRun.histo['PHchargeNoPedestal'].GetMean())
							signalSigmaRunChargeNP = np.sqrt(np.power(anaRun.histo['PHchargeNoPedestal'].GetRMS(), 2, dtype='f8') - np.power(anaRun.pedestal_charge_sigma, 2, dtype='f8'), dtype='f8')  if anaRun.histo['PHchargeNoPedestal'].GetRMS() > anaRun.pedestal_charge_sigma else anaRun.histo['PHchargeNoPedestal'].GetRMS()
						if 'PHvcal_CFNoPedestal' in list(anaRun.histo.keys()):
							signalRunVcalCFNP = np.double(anaRun.histo['PHvcal_CFNoPedestal'].GetMean())
							signalSigmaRunVcalCFNP = np.sqrt(np.power(anaRun.histo['PHvcal_CFNoPedestal'].GetRMS(), 2, dtype='f8') - np.power(anaRun.pedestal_vcal_sigma, 2, dtype='f8'), dtype='f8') if anaRun.histo['PHvcal_CFNoPedestal'].GetRMS() > anaRun.pedestal_vcal_sigma else anaRun.histo['PHvcal_CFNoPedestal'].GetRMS()
						if 'PHcharge_CFNoPedestal' in list(anaRun.histo.keys()):
							signalRunChargeCFNP = np.double(anaRun.histo['PHcharge_CFNoPedestal'].GetMean())
							signalSigmaRunChargeCFNP = np.sqrt(np.power(anaRun.histo['PHcharge_CFNoPedestal'].GetRMS(), 2, dtype='f8') - np.power(anaRun.pedestal_charge_sigma, 2, dtype='f8'), dtype='f8')  if anaRun.histo['PHcharge_CFNoPedestal'].GetRMS() > anaRun.pedestal_charge_sigma else anaRun.histo['PHcharge_CFNoPedestal'].GetRMS()

					if not self.are_cal_runs or '_out_' in run:
						self.signalOut[anaRun.bias] = signalRun if anaRun.bias < 0 else -signalRun
						self.signalOutSigma[anaRun.bias] = signalSigmaRun
						self.signalOutCF[anaRun.bias] = signalRunCF if anaRun.bias < 0 else -signalRunCF
						self.signalOutSigmaCF[anaRun.bias] = signalSigmaRunCF

						self.signalOutNP[anaRun.bias] = signalRunNP if anaRun.bias < 0 else -signalRunNP
						self.signalOutSigmaNP[anaRun.bias] = signalSigmaRunNP
						self.signalOutCFNP[anaRun.bias] = signalRunCFNP if anaRun.bias < 0 else -signalRunCFNP
						self.signalOutSigmaCFNP[anaRun.bias] = signalSigmaRunCFNP

						if not self.are_cal_runs:
							self.signalOutVcal[anaRun.bias] = -signalRunVcal if anaRun.bias < 0 else signalRunVcal
							self.signalOutVcalSigma[anaRun.bias] = signalSigmaRunVcal
							self.signalOutCharge[anaRun.bias] = -signalRunCharge if anaRun.bias < 0 else signalRunCharge
							self.signalOutChargeSigma[anaRun.bias] = signalSigmaRunCharge
							self.signalOutVcalCF[anaRun.bias] = -signalRunVcalCF if anaRun.bias < 0 else signalRunVcalCF
							self.signalOutVcalSigmaCF[anaRun.bias] = signalSigmaRunVcalCF
							self.signalOutChargeCF[anaRun.bias] = -signalRunChargeCF if anaRun.bias < 0 else signalRunChargeCF
							self.signalOutChargeSigmaCF[anaRun.bias] = signalSigmaRunChargeCF

							self.signalOutVcalNP[anaRun.bias] = -signalRunVcalNP if anaRun.bias < 0 else signalRunVcalNP
							self.signalOutVcalSigmaNP[anaRun.bias] = signalSigmaRunVcalNP
							self.signalOutChargeNP[anaRun.bias] = -signalRunChargeNP if anaRun.bias < 0 else signalRunChargeNP
							self.signalOutChargeSigmaNP[anaRun.bias] = signalSigmaRunChargeNP
							self.signalOutVcalCFNP[anaRun.bias] = -signalRunVcalCFNP if anaRun.bias < 0 else signalRunVcalCFNP
							self.signalOutVcalSigmaCFNP[anaRun.bias] = signalSigmaRunVcalCFNP
							self.signalOutChargeCFNP[anaRun.bias] = -signalRunChargeCFNP if anaRun.bias < 0 else signalRunChargeCFNP
							self.signalOutChargeSigmaCFNP[anaRun.bias] = signalSigmaRunChargeCFNP
						if '_out_' in run:
							self.signalPeakTime[anaRun.bias] = peakTime
							self.signalPeakTimeSigma[anaRun.bias] = peakTimeSigma
							self.signalPeakTimeCF[anaRun.bias] = peakTimeCF
							self.signalPeakTimeSigmaCF[anaRun.bias] = peakTimeSigmaCF
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
					if volt in list(self.signalOut.keys()) and volt in list(self.signalIn.keys()):
						voltages2.append(volt)
				self.voltages = voltages2
				stepsIn = np.array([self.signalIn[volt] for volt in self.voltages], dtype='f8')
				stepsInErrs = np.array([self.signalInSigma[volt] for volt in self.voltages], dtype='f8')
				signalO = np.array([self.signalOut[volt] for volt in self.voltages], dtype='f8')
				signalOutErrs = np.array([self.signalOutSigma[volt] for volt in self.voltages], dtype='f8')
				signalOCF = np.array([self.signalOutCF[volt] for volt in self.voltages], dtype='f8')
				signalOutErrsCF = np.array([self.signalOutSigmaCF[volt] for volt in self.voltages], dtype='f8')
				self.graphPH = ro.TGraphErrors(len(self.voltages), stepsIn, signalO, stepsInErrs, signalOutErrs)
				self.graphPH.SetNameTitle('Signal_vs_CalStep_ch_' + str(self.caen_ch), 'Signal_vs_CalStep_ch_' + str(self.caen_ch))
				self.graphPH.GetXaxis().SetTitle('vcal step [mV]')
				self.graphPH.GetYaxis().SetTitle('signal [mV]')
				self.graphPHCF = ro.TGraphErrors(len(self.voltages), stepsIn, signalOCF, stepsInErrs, signalOutErrsCF)
				self.graphPHCF.SetNameTitle('Signal_CF_vs_CalStep_ch_' + str(self.caen_ch), 'Signal_CF_vs_CalStep_ch_' + str(self.caen_ch))
				self.graphPHCF.GetXaxis().SetTitle('vcal step [mV]')
				self.graphPHCF.GetYaxis().SetTitle('signal_CF [mV]')

				signalONP = np.array([self.signalOutNP[volt] for volt in self.voltages], dtype='f8')
				signalOutErrsNP = np.array([self.signalOutSigmaNP[volt] for volt in self.voltages], dtype='f8')
				signalOCFNP = np.array([self.signalOutCFNP[volt] for volt in self.voltages], dtype='f8')
				signalOutErrsCFNP = np.array([self.signalOutSigmaCFNP[volt] for volt in self.voltages], dtype='f8')
				self.graphPHNP = ro.TGraphErrors(len(self.voltages), stepsIn, signalONP, stepsInErrs, signalOutErrsNP)
				self.graphPHNP.SetNameTitle('Signal_NoPedestal_vs_CalStep_ch_' + str(self.caen_ch), 'Signal_NoPedestal_vs_CalStep_ch_' + str(self.caen_ch))
				self.graphPHNP.GetXaxis().SetTitle('vcal step [mV]')
				self.graphPHNP.GetYaxis().SetTitle('signal [mV]')
				self.graphPHCFNP = ro.TGraphErrors(len(self.voltages), stepsIn, signalOCFNP, stepsInErrs, signalOutErrsCFNP)
				self.graphPHCFNP.SetNameTitle('Signal_NoPedestal_CF_vs_CalStep_ch_' + str(self.caen_ch), 'Signal_NoPedestal_CF_vs_CalStep_ch_' + str(self.caen_ch))
				self.graphPHCFNP.GetXaxis().SetTitle('vcal step [mV]')
				self.graphPHCFNP.GetYaxis().SetTitle('signal_CF [mV]')

			else:
				diaVoltages = np.array([self.diaVoltages[volt] for volt in self.voltages], dtype='f8')
				diaVoltagesSigma = np.array([self.diaVoltagesSigma[volt] for volt in self.voltages], dtype='f8')

				signalO = np.array([self.signalOut[volt] for volt in self.voltages], dtype='f8')
				signalOutErrs = np.array([self.signalOutSigma[volt] for volt in self.voltages], dtype='f8')
				signalOVcal = np.array([self.signalOutVcal[volt] for volt in self.voltages], dtype='f8')
				signalOutVcalErrs = np.array([self.signalOutVcalSigma[volt] for volt in self.voltages], dtype='f8')
				signalOCharge = np.array([self.signalOutCharge[volt] for volt in self.voltages], dtype='f8')
				signalOutChargeErrs = np.array([self.signalOutChargeSigma[volt] for volt in self.voltages], dtype='f8')
				signalOCF = np.array([self.signalOutCF[volt] for volt in self.voltages], dtype='f8')
				signalOutErrsCF = np.array([self.signalOutSigmaCF[volt] for volt in self.voltages], dtype='f8')
				signalOVcalCF = np.array([self.signalOutVcalCF[volt] for volt in self.voltages], dtype='f8')
				signalOutVcalErrsCF = np.array([self.signalOutVcalSigmaCF[volt] for volt in self.voltages], dtype='f8')
				signalOChargeCF = np.array([self.signalOutChargeCF[volt] for volt in self.voltages], dtype='f8')
				signalOutChargeErrsCF = np.array([self.signalOutChargeSigmaCF[volt] for volt in self.voltages], dtype='f8')

				signalONP = np.array([self.signalOutNP[volt] for volt in self.voltages], dtype='f8')
				signalOutErrsNP = np.array([self.signalOutSigmaNP[volt] for volt in self.voltages], dtype='f8')
				signalOVcalNP = np.array([self.signalOutVcalNP[volt] for volt in self.voltages], dtype='f8')
				signalOutVcalErrsNP = np.array([self.signalOutVcalSigmaNP[volt] for volt in self.voltages], dtype='f8')
				signalOChargeNP = np.array([self.signalOutChargeNP[volt] for volt in self.voltages], dtype='f8')
				signalOutChargeErrsNP = np.array([self.signalOutChargeSigmaNP[volt] for volt in self.voltages], dtype='f8')
				signalOCFNP = np.array([self.signalOutCFNP[volt] for volt in self.voltages], dtype='f8')
				signalOutErrsCFNP = np.array([self.signalOutSigmaCFNP[volt] for volt in self.voltages], dtype='f8')
				signalOVcalCFNP = np.array([self.signalOutVcalCFNP[volt] for volt in self.voltages], dtype='f8')
				signalOutVcalErrsCFNP = np.array([self.signalOutVcalSigmaCFNP[volt] for volt in self.voltages], dtype='f8')
				signalOChargeCFNP = np.array([self.signalOutChargeCFNP[volt] for volt in self.voltages], dtype='f8')
				signalOutChargeErrsCFNP = np.array([self.signalOutChargeSigmaCFNP[volt] for volt in self.voltages], dtype='f8')

				self.graphPH = ro.TGraphErrors(len(self.voltages), diaVoltages, signalO, diaVoltagesSigma, signalOutErrs)
				self.graphPH.SetNameTitle('Signal_vs_HV', 'Signal_vs_HV')
				self.graphPH.GetXaxis().SetTitle('HV [V]')
				self.graphPH.GetYaxis().SetTitle('signal [mV]')
				self.graphPHCF = ro.TGraphErrors(len(self.voltages), diaVoltages, signalOCF, diaVoltagesSigma, signalOutErrsCF)
				self.graphPHCF.SetNameTitle('Signal_CF_vs_HV', 'Signal_CF_vs_HV')
				self.graphPHCF.GetXaxis().SetTitle('HV [V]')
				self.graphPHCF.GetYaxis().SetTitle('signal_CF [mV]')

				self.graphVcal = ro.TGraphErrors(len(self.voltages), diaVoltages, signalOVcal, diaVoltagesSigma, signalOutVcalErrs)
				self.graphVcal.SetNameTitle('SignalVcal_vs_HV', 'SignalVcal_vs_HV')
				self.graphVcal.GetXaxis().SetTitle('HV [V]')
				self.graphVcal.GetYaxis().SetTitle('signalVcal [mV]')
				self.graphVcalCF = ro.TGraphErrors(len(self.voltages), diaVoltages, signalOVcalCF, diaVoltagesSigma, signalOutVcalErrsCF)
				self.graphVcalCF.SetNameTitle('SignalVcal_CF_vs_HV', 'SignalVcal_CF_vs_HV')
				self.graphVcalCF.GetXaxis().SetTitle('HV [V]')
				self.graphVcalCF.GetYaxis().SetTitle('signalVcal_CF [mV]')

				self.graphCharge = ro.TGraphErrors(len(self.voltages), diaVoltages, signalOCharge, diaVoltagesSigma, signalOutChargeErrs)
				self.graphCharge.SetNameTitle('SignalCharge_vs_HV', 'SignalCharge_vs_HV')
				self.graphCharge.GetXaxis().SetTitle('HV [V]')
				self.graphCharge.GetYaxis().SetTitle('signalCharge [e]')
				self.graphChargeCF = ro.TGraphErrors(len(self.voltages), diaVoltages, signalOChargeCF, diaVoltagesSigma, signalOutChargeErrsCF)
				self.graphChargeCF.SetNameTitle('SignalCharge_CF_vs_HV', 'SignalCharge_CF_vs_HV')
				self.graphChargeCF.GetXaxis().SetTitle('HV [V]')
				self.graphChargeCF.GetYaxis().SetTitle('signalCharge_CF [e]')

				self.graphPHNP = ro.TGraphErrors(len(self.voltages), diaVoltages, signalONP, diaVoltagesSigma, signalOutErrsNP)
				self.graphPHNP.SetNameTitle('Signal_NoPedestal_vs_HV', 'Signal_NoPedestal_vs_HV')
				self.graphPHNP.GetXaxis().SetTitle('HV [V]')
				self.graphPHNP.GetYaxis().SetTitle('signal [mV]')
				self.graphPHCFNP = ro.TGraphErrors(len(self.voltages), diaVoltages, signalOCFNP, diaVoltagesSigma, signalOutErrsCFNP)
				self.graphPHCFNP.SetNameTitle('Signal_NoPedestal_CF_vs_HV', 'Signal_NoPedestal_CF_vs_HV')
				self.graphPHCFNP.GetXaxis().SetTitle('HV [V]')
				self.graphPHCFNP.GetYaxis().SetTitle('signal_CF [mV]')

				self.graphVcalNP = ro.TGraphErrors(len(self.voltages), diaVoltages, signalOVcalNP, diaVoltagesSigma, signalOutVcalErrsNP)
				self.graphVcalNP.SetNameTitle('SignalVcal_NoPedestal_vs_HV', 'SignalVcal_NoPedestal_vs_HV')
				self.graphVcalNP.GetXaxis().SetTitle('HV [V]')
				self.graphVcalNP.GetYaxis().SetTitle('signalVcal [mV]')
				self.graphVcalCFNP = ro.TGraphErrors(len(self.voltages), diaVoltages, signalOVcalCFNP, diaVoltagesSigma, signalOutVcalErrsCFNP)
				self.graphVcalCFNP.SetNameTitle('SignalVcal_NoPedestal_CF_vs_HV', 'SignalVcal_NoPedestal_CF_vs_HV')
				self.graphVcalCFNP.GetXaxis().SetTitle('HV [V]')
				self.graphVcalCFNP.GetYaxis().SetTitle('signalVcal_CF [mV]')

				self.graphChargeNP = ro.TGraphErrors(len(self.voltages), diaVoltages, signalOChargeNP, diaVoltagesSigma, signalOutChargeErrsNP)
				self.graphChargeNP.SetNameTitle('SignalCharge_NoPedestal_vs_HV', 'SignalCharge_NoPedestal_vs_HV')
				self.graphChargeNP.GetXaxis().SetTitle('HV [V]')
				self.graphChargeNP.GetYaxis().SetTitle('signalCharge [e]')
				self.graphChargeCFNP = ro.TGraphErrors(len(self.voltages), diaVoltages, signalOChargeCFNP, diaVoltagesSigma, signalOutChargeErrsCFNP)
				self.graphChargeCFNP.SetNameTitle('SignalCharge_NoPedestal_CF_vs_HV', 'SignalCharge_NoPedestal_CF_vs_HV')
				self.graphChargeCFNP.GetXaxis().SetTitle('HV [V]')
				self.graphChargeCFNP.GetYaxis().SetTitle('signalCharge_CF [e]')

			self.graphPH.SetMarkerStyle(7)
			self.graphPH.SetMarkerColor(ro.kBlack)
			self.graphPH.SetLineColor(ro.kBlack)
			self.graphPHCF.SetMarkerStyle(7)
			self.graphPHCF.SetMarkerColor(ro.kBlack)
			self.graphPHCF.SetLineColor(ro.kBlack)

			self.canvasPH = ro.TCanvas('c_' + self.graphPH.GetName(), 'c_' + self.graphPH.GetName(), 1)
			self.graphPH.Draw('AP')
			self.canvasPH.SetGridx()
			self.canvasPH.SetGridy()
			self.canvasPHCF = ro.TCanvas('c_' + self.graphPHCF.GetName(), 'c_' + self.graphPHCF.GetName(), 1)
			self.graphPHCF.Draw('AP')
			self.canvasPHCF.SetGridx()
			self.canvasPHCF.SetGridy()

			self.graphPHNP.SetMarkerStyle(7)
			self.graphPHNP.SetMarkerColor(ro.kBlack)
			self.graphPHNP.SetLineColor(ro.kBlack)
			self.graphPHCFNP.SetMarkerStyle(7)
			self.graphPHCFNP.SetMarkerColor(ro.kBlack)
			self.graphPHCFNP.SetLineColor(ro.kBlack)

			self.canvasPHNP = ro.TCanvas('c_' + self.graphPHNP.GetName(), 'c_' + self.graphPHNP.GetName(), 1)
			self.graphPHNP.Draw('AP')
			self.canvasPHNP.SetGridx()
			self.canvasPHNP.SetGridy()
			self.canvasPHCFNP = ro.TCanvas('c_' + self.graphPHCFNP.GetName(), 'c_' + self.graphPHCFNP.GetName(), 1)
			self.graphPHCFNP.Draw('AP')
			self.canvasPHCFNP.SetGridx()
			self.canvasPHCFNP.SetGridy()

			if not self.are_cal_runs:
				self.graphVcal.SetMarkerStyle(7)
				self.graphVcal.SetMarkerColor(ro.kBlack)
				self.graphVcal.SetLineColor(ro.kBlack)
				self.graphVcalCF.SetMarkerStyle(7)
				self.graphVcalCF.SetMarkerColor(ro.kBlack)
				self.graphVcalCF.SetLineColor(ro.kBlack)

				self.canvasVcal = ro.TCanvas('c_' + self.graphVcal.GetName(), 'c_' + self.graphVcal.GetName(), 1)
				self.graphVcal.Draw('AP')
				self.canvasVcal.SetGridx()
				self.canvasVcal.SetGridy()
				self.canvasVcalCF = ro.TCanvas('c_' + self.graphVcalCF.GetName(), 'c_' + self.graphVcalCF.GetName(), 1)
				self.graphVcalCF.Draw('AP')
				self.canvasVcalCF.SetGridx()
				self.canvasVcalCF.SetGridy()

				self.graphCharge.SetMarkerStyle(7)
				self.graphCharge.SetMarkerColor(ro.kBlack)
				self.graphCharge.SetLineColor(ro.kBlack)
				self.graphChargeCF.SetMarkerStyle(7)
				self.graphChargeCF.SetMarkerColor(ro.kBlack)
				self.graphChargeCF.SetLineColor(ro.kBlack)

				self.canvasCharge = ro.TCanvas('c_' + self.graphCharge.GetName(), 'c_' + self.graphCharge.GetName(), 1)
				self.graphCharge.Draw('AP')
				self.canvasCharge.SetGridx()
				self.canvasCharge.SetGridy()
				self.canvasChargeCF = ro.TCanvas('c_' + self.graphChargeCF.GetName(), 'c_' + self.graphChargeCF.GetName(), 1)
				self.graphChargeCF.Draw('AP')
				self.canvasChargeCF.SetGridx()
				self.canvasChargeCF.SetGridy()

				self.graphVcalNP.SetMarkerStyle(7)
				self.graphVcalNP.SetMarkerColor(ro.kBlack)
				self.graphVcalNP.SetLineColor(ro.kBlack)
				self.graphVcalCFNP.SetMarkerStyle(7)
				self.graphVcalCFNP.SetMarkerColor(ro.kBlack)
				self.graphVcalCFNP.SetLineColor(ro.kBlack)

				self.canvasVcalNP = ro.TCanvas('c_' + self.graphVcalNP.GetName(), 'c_' + self.graphVcalNP.GetName(), 1)
				self.graphVcalNP.Draw('AP')
				self.canvasVcalNP.SetGridx()
				self.canvasVcalNP.SetGridy()
				self.canvasVcalCFNP = ro.TCanvas('c_' + self.graphVcalCFNP.GetName(), 'c_' + self.graphVcalCFNP.GetName(), 1)
				self.graphVcalCFNP.Draw('AP')
				self.canvasVcalCFNP.SetGridx()
				self.canvasVcalCFNP.SetGridy()

				self.graphChargeNP.SetMarkerStyle(7)
				self.graphChargeNP.SetMarkerColor(ro.kBlack)
				self.graphChargeNP.SetLineColor(ro.kBlack)
				self.graphChargeCFNP.SetMarkerStyle(7)
				self.graphChargeCFNP.SetMarkerColor(ro.kBlack)
				self.graphChargeCFNP.SetLineColor(ro.kBlack)

				self.canvasChargeNP = ro.TCanvas('c_' + self.graphChargeNP.GetName(), 'c_' + self.graphChargeNP.GetName(), 1)
				self.graphChargeNP.Draw('AP')
				self.canvasChargeNP.SetGridx()
				self.canvasChargeNP.SetGridy()
				self.canvasChargeCFNP = ro.TCanvas('c_' + self.graphChargeCFNP.GetName(), 'c_' + self.graphChargeCFNP.GetName(), 1)
				self.graphChargeCFNP.Draw('AP')
				self.canvasChargeCFNP.SetGridx()
				self.canvasChargeCFNP.SetGridy()

	def PlotPeakTimes(self):
		if len(self.voltages) > 0:
			if self.are_cal_runs:
				if len(list(self.signalPeakTime.values())) < 1:
					for run in self.runs:
						print('Analysing run:', run)
						if os.path.isdir(run):
							root_files = glob.glob('{d}/*.root'.format(d=run))
							if len(root_files) > 0:
								if '_out_' in run:
									configfile = self.config
									print('Using config file:', configfile)
									print('Loading analysis tree if it exists')
									anaRun = AnalysisCaenCCD(run, configfile, overw=False)
									anaRun.AnalysisWaves()
									peakTime = 0
									peakTimeSigma = 0
									peakTimeCF = 0
									peakTimeSigmaCF = 0
									if 'peakPosDist' in list(anaRun.histo.keys()):
										peakTime = np.double(anaRun.peakTime * 1e6)
										peakTimeSigma = np.double(anaRun.histo['peakPosDist'].GetRMS())
									if 'peakPosDistCF' in list(anaRun.histo.keys()):
										peakTimeCF = np.double(anaRun.peakTime * 1e6)
										peakTimeSigmaCF = np.double(anaRun.histo['peakPosDistCF'].GetRMS())
									self.signalPeakTime[anaRun.bias] = peakTime
									self.signalPeakTimeSigma[anaRun.bias] = peakTimeSigma
									self.signalPeakTimeCF[anaRun.bias] = peakTimeCF
									self.signalPeakTimeSigmaCF[anaRun.bias] = peakTimeSigmaCF

				signalO = np.array([self.signalOut[volt] for volt in self.voltages], dtype='f8')
				signalOutErrs = np.array([self.signalOutSigma[volt] for volt in self.voltages], dtype='f8')
				signalPeaks = np.array([self.signalPeakTime[volt] for volt in self.voltages], 'f8')
				signalPeaksSigma = np.array([self.signalPeakTimeSigma[volt] for volt in self.voltages], 'f8')

				signalOCF = np.array([self.signalOutCF[volt] for volt in self.voltages], dtype='f8')
				signalOutErrsCF = np.array([self.signalOutSigmaCF[volt] for volt in self.voltages], dtype='f8')
				signalPeaksCF = np.array([self.signalPeakTimeCF[volt] for volt in self.voltages], 'f8')
				signalPeaksSigmaCF = np.array([self.signalPeakTimeSigmaCF[volt] for volt in self.voltages], 'f8')

				self.graphPeakTime = ro.TGraphErrors(len(self.voltages), signalO, signalPeaks, signalOutErrs, signalPeaksSigma)
				self.graphPeakTime.SetNameTitle('SignalPeakTime_vs_Signal_ch_' + str(self.caen_ch), 'SignalPeakTime_vs_Signal_ch_' + str(self.caen_ch))
				self.graphPeakTime.GetXaxis().SetTitle('signal [mV]')
				self.graphPeakTime.GetYaxis().SetTitle('Peak Time [us]')
				self.graphPeakTime.SetMarkerStyle(7)
				self.graphPeakTime.SetMarkerColor(ro.kBlack)
				self.graphPeakTime.SetLineColor(ro.kBlack)

				self.canvasPeakTime = ro.TCanvas('c_' + self.graphPeakTime.GetName(), 'c_' + self.graphPeakTime.GetName(), 1)
				self.graphPeakTime.Draw('AP')
				self.canvasPeakTime.SetGridx()
				self.canvasPeakTime.SetGridy()

				self.graphPeakTimeCF = ro.TGraphErrors(len(self.voltages), signalOCF, signalPeaksCF, signalOutErrsCF, signalPeaksSigmaCF)
				self.graphPeakTimeCF.SetNameTitle('SignalPeakTime_CF_vs_Signal_CF_ch_' + str(self.caen_ch), 'SignalPeakTime_CF_vs_Signal_CF_ch_' + str(self.caen_ch))
				self.graphPeakTimeCF.GetXaxis().SetTitle('signal [mV]')
				self.graphPeakTimeCF.GetYaxis().SetTitle('Peak Time [us]')
				self.graphPeakTimeCF.SetMarkerStyle(7)
				self.graphPeakTimeCF.SetMarkerColor(ro.kBlack)
				self.graphPeakTimeCF.SetLineColor(ro.kBlack)

				self.canvasPeakTimeCF = ro.TCanvas('c_' + self.graphPeakTimeCF.GetName(), 'c_' + self.graphPeakTimeCF.GetName(), 1)
				self.graphPeakTimeCF.Draw('AP')
				self.canvasPeakTimeCF.SetGridx()
				self.canvasPeakTimeCF.SetGridy()

	def FitLine(self):
		ro.Math.MinimizerOptions.SetDefaultMinimizer(*fit_method)
		ro.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(1000000)
		ro.Math.MinimizerOptions.SetDefaultTolerance(0.00001)
		ro.gStyle.SetOptFit(1111)
		func = ro.TF1('fit_' + self.graphPH.GetName(), 'pol1', -1000, 1000)
		func.SetLineColor(ro.kRed)
		func.SetNpx(10000)
		tfit = self.graphPH.Fit('fit_' + self.graphPH.GetName(), 'QM0S', '', -1000, 1000)
		if tfit.Prob() < 0.9:
			tfit = self.graphPH.Fit('fit_' + self.graphPH.GetName(), 'QM0S', '', -1000, 1000)
		if tfit.Prob() < 0.9:
			tfit = self.graphPH.Fit('fit_' + self.graphPH.GetName(), 'QM0S', '', -1000, 1000)
		# func.SetParameters(np.array([tfit.Parameter(0), tfit.Parameter(1)], 'f8'))
		# func.SetChisquare(tfit.Chi2())
		# func.SetNDF(tfit.Ndf())
		func = tfit.FittedFunction().GetFunction()
		self.fit = func.Clone()
		self.canvasPH.cd()
		self.fit.Draw('same')
		self.canvasPH.Modified()
		ro.gPad.Update()

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
		                   'signal_peak_time': self.signalPeakTime,
		                   'signal_peak_time_sigma': self.signalPeakTimeSigma,
						   'caen_ch': self.caen_ch,
						   'fit_p0': self.fit.GetParameter(0) if self.fit else 0,
						   'fit_p0_error': self.fit.GetParError(0) if self.fit else 0,
						   'fit_p1': self.fit.GetParameter(1) if self.fit else 0,
						   'fit_p1_error': self.fit.GetParError(1) if self.fit else 0,
						   'fit_prob': self.fit.GetProb() if self.fit else 0,
						   'fit_chi2': self.fit.GetChisquare() if self.fit else 0,
						   'fit_ndf': self.fit.GetNDF() if self.fit else 0
						   }

	def SavePickle(self, overWritePickle=False):
		pickleName = 'signal_cal_{c}.cal'.format(c=self.caen_ch) if self.are_cal_runs else 'signal_{c}.cal'.format(c=self.caen_ch)
		if not self.cal_pickle:
			self.FillPickle()
		if os.path.isfile('{d}/{f}'.format(d=self.runsdir, f=pickleName)):
			if not self.overwrite and not overWritePickle:
				print('The file', pickleName, 'already exists in', self.runsdir, '. Not saving!')
				return
		with open('{d}/{f}'.format(d=self.runsdir, f=pickleName), 'wb') as fpickle:
			pickle.dump(self.cal_pickle, fpickle, pickle.HIGHEST_PROTOCOL)
		print('Saved pickle', pickleName, 'in', self.runsdir)

	def SaveCanvas(self, canvas, graph):
		if canvas and graph:
			canvas.SaveAs('{d}/{n}.png'.format(d=self.runsdir, n=graph.GetName()))
			canvas.SaveAs('{d}/{n}.root'.format(d=self.runsdir, n=graph.GetName()))


	def SaveAllCanvas(self):
		self.SaveCanvas(self.canvasPH, self.graphPH)
		self.SaveCanvas(self.canvasPHCF, self.graphPHCF)
		self.SaveCanvas(self.canvasPHNP, self.graphPHNP)
		self.SaveCanvas(self.canvasPHCFNP, self.graphPHCFNP)

		self.SaveCanvas(self.canvasVcal, self.graphVcal)
		self.SaveCanvas(self.canvasVcalCF, self.graphVcalCF)
		self.SaveCanvas(self.canvasVcalNP, self.graphVcalNP)
		self.SaveCanvas(self.canvasVcalCFNP, self.graphVcalCFNP)

		self.SaveCanvas(self.canvasCharge, self.graphCharge)
		self.SaveCanvas(self.canvasChargeCF, self.graphChargeCF)
		self.SaveCanvas(self.canvasChargeNP, self.graphChargeNP)
		self.SaveCanvas(self.canvasChargeCFNP, self.graphChargeCFNP)

		self.SaveCanvas(self.canvasPeakTime, self.graphPeakTime)
		self.SaveCanvas(self.canvasPeakTimeCF, self.graphPeakTimeCF)

def main():
	parser = OptionParser()
	parser.add_option('-d', '--runsdir', dest='runsdir', type='str', default='', help='path to folder containing all the run folders to modify')
	parser.add_option('-c', '--config', dest='config', type='str', default='', help='path to analysis config file used for all runs inside the folder')
	parser.add_option('-o', '--overwrite', dest='overwrite', default=False, action='store_true', help='Sets overwrite of analysis tree for the run if it exists')
	parser.add_option('--configinput', dest='configinput', type='str', default='', help='path to analysis config file used for all runs inside the folder that are of calibration type "in". Only necessary when doing signal calibration')
	parser.add_option('--debug', dest='dodebug', default=False, help='Toggles debug mode', action='store_true')
	(options, args) = parser.parse_args()
	runsdir = str(options.runsdir)
	config = str(options.config)
	configinput = str(options.configinput)
	overwrite = bool(options.overwrite)
	dodebug = bool(options.dodebug)
	arif = AnalysisAllRunsInFolder(runsdir, config, configinput, overwrite, dodebug)
	return arif

if __name__ == '__main__':
	arif = main()
