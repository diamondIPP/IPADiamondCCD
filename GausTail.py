#!/usr/bin/env python

import ROOT as ro
import numpy as np
import ipdb

class GausTail:
	def __init__(self, graph, bias=1, peakTime=2e-6):
		# (0)amp, (1)peakShift, (2)sigma, (3)p0, (4)p1
		self.bias = bias
		self.peakTime = peakTime
		self.graph = graph
		yvect, npoints = self.graph.GetY(), self.graph.GetN()
		yvect = [yvect[i] for i in xrange(npoints)]
		par0 = max(yvect) if bias < 0 else min(yvect)
		self.params = np.array([par0, peakTime, 0.5e-6, -1, 0.1e-6], 'f8') if bias < 0 else np.array([-1, peakTime, 0.5e-6, 1, -0.1e-6], 'f8')
		self.paramErrors = np.zeros(5, 'f8')
		self.paramsLimitsLow = np.zeros(5, 'f8')
		self.paramsLimitsHigh = np.zeros(5, 'f8')
		self.paramsFitErrors = np.zeros(5, 'f8')
		self.chi2 = 0
		self.ndf = 0
		# self.fit_range = [max(self.histo.GetMean()*0.3, 0), min(self.histo.GetXaxis().GetXmax(), self.histo.GetMean()*3)]
		self.fit_range = {'min': 0.5e-6, 'max': 5e-6}
		# self.fit_range = [0, min(self.histo.GetXaxis().GetXmax(), self.histo.GetMean()*5)]
		self.fit = None
		self.resFit = None
		self.ParametersLimits()
		self.fit_mp = 0
		self.fit_minx_1 = 0
		self.histo_bin_1 = 1
		self.entries_under_curve = 0

	def ParametersLimits(self):
		(self.paramsLimitsLow[0], self.paramsLimitsHigh[0]) = (0, 100) if self.bias < 0 else (-100, 0)
		(self.paramsLimitsLow[1], self.paramsLimitsHigh[1]) = (self.peakTime -1.5e-6, self.peakTime + 1.5e-6)
		(self.paramsLimitsLow[2], self.paramsLimitsHigh[2]) = (0, 10e-6)
		(self.paramsLimitsLow[3], self.paramsLimitsHigh[3]) = (-100, 100)
		(self.paramsLimitsLow[4], self.paramsLimitsHigh[4]) = (0, 1e6) if self.bias < 0 else (-1e6, 0)

	def FitGausTail(self, xmin=0.5e-6, xmax=5e-6):
		fit_name = 'fit_' + self.graph.GetTitle()
		fit_old = ro.gROOT.GetListOfFunctions().FindObject(fit_name)
		if fit_old:
			fit_old.Delete()
			del fit_old
		self.fit = ro.TF1(fit_name, '[0]*exp(-((x-[1])/(1.4142135623730951*[2]))^2)+pol1(3)', self.fit_range['min'], self.fit_range['max'])
		self.fit.SetNpx(1000)
		self.fit.SetParameters(self.params)
		self.fit.SetParNames('N', 'peakPosShift', 'sigma', 'peakPos', 'Sigma')
		options = 'QBN0S'
		ro.Math.MinimizerOptions.SetDefaultMinimizer('Minuit2', 'Migrad')
		for i in xrange(len(self.params)):
			self.fit.SetParLimits(i, self.paramsLimitsLow[i], self.paramsLimitsHigh[i])
		self.resFit = self.graph.Fit(fit_name, options, '', xmin, xmax)
		self.resFit = self.graph.Fit(fit_name, options, '', xmin, xmax)
		self.fit.GetParameters(self.params)
		for i in xrange(len(self.params)):
			self.paramsFitErrors[i] = self.fit.GetParError(i)
		self.chi2 = self.fit.GetChisquare()
		self.ndf = self.fit.GetNDF()
		self.peakTime = self.fit.GetMaximumX(xmin, xmax) if self.bias < 0 else self.fit.GetMinimumX(xmin, xmax)


if __name__ == '__main__':
	print 'Bla'