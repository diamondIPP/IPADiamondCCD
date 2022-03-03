#!/usr/bin/env python

import ROOT as ro
import numpy as np
import ipdb

class Crystalball:
	def __init__(self, graph, bias=1, alpha=0.6, n=4.5, peakpos=2e-6, sigma=1e-6):
		# (0)amp, (1)alpha, (2)n, (3)peak, (4)sigma, (5)offset
		self.graph = graph
		self.bias = bias
		yvect, npoints = self.graph.GetY(), self.graph.GetN()
		yvect = [yvect[i] for i in range(npoints)]
		par0 = max(yvect) if bias < 0 else min(yvect)
		self.params = np.array([par0, alpha, n, peakpos, sigma, 0], 'f8')
		self.paramErrors = np.zeros(6, 'f8')
		self.paramsLimitsLow = np.zeros(6, 'f8')
		self.paramsLimitsHigh = np.zeros(6, 'f8')
		self.paramsFitErrors = np.zeros(6, 'f8')
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
		self.peakTime = peakpos

	def ParametersLimits(self):
		self.paramsLimitsLow[0], self.paramsLimitsHigh[0] = -10,1000
		# self.paramsLimitsLow[0], self.paramsLimitsHigh[0] = -10, 10
		self.paramsLimitsLow[1], self.paramsLimitsHigh[1] = 0, 5
		# self.paramsLimitsLow[1], self.paramsLimitsHigh[1] = 0, 100
		self.paramsLimitsLow[2], self.paramsLimitsHigh[2] = 0, 20
		# self.paramsLimitsLow[2], self.paramsLimitsHigh[2] = 0, 100
		self.paramsLimitsLow[3], self.paramsLimitsHigh[3] = 70, 100
		# self.paramsLimitsLow[3], self.paramsLimitsHigh[3] = self.fit_range['min'], self.fit_range['max']
		self.paramsLimitsLow[4], self.paramsLimitsHigh[4] = 0, 20
		# self.paramsLimitsLow[4], self.paramsLimitsHigh[4] = 0, 10e-6
		self.paramsLimitsLow[5], self.paramsLimitsHigh[5] = -1000, 1000
		# self.paramsLimitsLow[5], self.paramsLimitsHigh[5] = -2, 2

	def CrystalballFunc(self, x, params):
		if x[0] - params[3] < params[4] * params[1]:
			# return params[0] * np.exp(-np.power((x[0] - params[3]) / (np.sqrt(2) * params[4]), 2, dtype='f8')) + params[5]
			return params[0] * ro.TMath.Exp(-(((x[0] - params[3]) / (1.4142135623730951 * params[4])) ** 2.)) + params[5]
		else:
			a = (params[2] / params[1]) ** params[2] * ro.TMath.Exp(-(params[1] ** 2.) / 2.)
			# a = np.power(params[2] / params[1], params[2], dtype='f8') * np.exp(-np.power(params[1], 2, dtype='f8')/2.)
			b = params[2] / float(params[1]) - params[1]
			# b = np.divide(params[2], params[1], dtype='f8') - params[1]
			return params[0] * a * (b + float(x[0] - params[3]) / params[4]) ** (-params[2]) + params[5]
			# return params[0] * a * np.power(b + float(x[0] - params[3]) / params[4], -params[2], dtype='f8') + params[5]

	def FitCrystalball(self, xmin=0.5e-6, xmax=5e-6):
		fit_name = 'fit_' + self.graph.GetTitle()
		fit_old = ro.gROOT.GetListOfFunctions().FindObject(fit_name)
		if fit_old:
			fit_old.Delete()
			del fit_old
		self.fit = ro.TF1(fit_name, self.CrystalballFunc, self.fit_range['min'], self.fit_range['max'], 6)
		self.fit.SetNpx(1000)
		self.fit.SetParameters(self.params)
		self.fit.SetParNames('N', 'Alpha', 'n', 'peakPos', 'Sigma')
		options = 'QBN0S'
		ro.Math.MinimizerOptions.SetDefaultMinimizer('Minuit2', 'Migrad')
		for i in range(len(self.params)):
			self.fit.SetParLimits(i, self.paramsLimitsLow[i], self.paramsLimitsHigh[i])
		self.resFit = self.graph.Fit(fit_name, options, '', xmin, xmax)
		self.resFit = self.graph.Fit(fit_name, options, '', xmin, xmax)
		self.fit.GetParameters(self.params)
		for i in range(len(self.params)):
			self.paramsFitErrors[i] = self.fit.GetParError(i)
		self.chi2 = self.fit.GetChisquare()
		self.ndf = self.fit.GetNDF()
		self.peakTime = self.params[3]


if __name__ == '__main__':
	print('Bla')