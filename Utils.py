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
from ConfigParser import ConfigParser
import subprocess as subp
import struct
import ROOT as ro
import shutil
from copy import deepcopy
import glob


# from DataAcquisition import DataAcquisition

class Utils:
	def __init__(self):
		self.bar = None

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

def FindRedundantEvents(settings):
	channels = [settings.sigCh, settings.trigCh, settings.acCh] if settings.ac_enable else [settings.sigCh, settings.trigCh]
	types = {settings.sigCh: 'signal_ch', settings.trigCh: 'trigger_ch', settings.acCh: 'veto'} if settings.ac_enable else {settings.sigCh: 'signal_ch', settings.trigCh: 'trigger_ch'}
	for ch in channels:
		print '\n', types[ch], ':\n'
		filename = '{d}/Runs/{f}/{f}_{t}.dat'.format(d=settings.outdir, f=settings.filename, t=types[ch])
		f0 = open(filename, 'rb')
		events = int(np.floor(os.path.getsize(filename) / settings.struct_len))
		for ev0 in xrange(events):
			f0.seek(ev0 * settings.struct_len)
			dataev0 = f0.read(settings.struct_len)
			dataev0 = struct.Struct(settings.struct_fmt).unpack_from(dataev0)
			dataev0 = np.array(dataev0, 'H')
			for ev1 in xrange(ev0+1, events):
				f0.seek(ev1 * settings.struct_len)
				dataev1 = f0.read(settings.struct_len)
				dataev1 = struct.Struct(settings.struct_fmt).unpack_from(dataev1)
				dataev1 = np.array(dataev1, 'H')

				if (dataev0 == dataev1).all():
					print '', ev0, '\t', ev1
		f0.close()

def IsInt(i):
	try:
		int(i)
		return True
	except ValueError:
		return False

def IsFloat(f):
	try:
		float(f)
		return True
	except ValueError:
		return False

def PlotHisto2DLimits(xmin=-0.48e-6, xmax=4.68e-6, ymin=-0.5, ymax=0.1, xres=2e-9, yres=2.15/(2.0**14-1)):
	xbins = int(round((xmax - xmin)/xres))
	ybins = int(round((ymax - ymin)/yres))
	print '({nx},{minx},{maxx},{ny},{miny},{maxy})'.format(nx=xbins, minx=xmin, maxx=xmax, ny=ybins, miny=ymin,maxy=ymax)

def ExitMessage(msg, code=os.EX_SOFTWARE):
	print '\n##########'
	print msg
	print '##########'
	sys.exit(code)

def SetDefault2DStats(histo):
	if histo.FindObject('stats'):
		histo.FindObject('stats').SetOptStat(11)
		histo.FindObject('stats').SetX1NDC(0.7)
		histo.FindObject('stats').SetX2NDC(0.9)
		histo.FindObject('stats').SetY1NDC(0.9)
		histo.FindObject('stats').SetY2NDC(0.975)
		ro.gPad.Update()
	if histo.FindObject('palette'):
		histo.FindObject('palette').SetX1NDC(0.87)
		histo.FindObject('palette').SetX2NDC(0.92)
		ro.gPad.Update()

def SetDefault1DStats(histo):
	if histo.FindObject('stats'):
		histo.FindObject('stats').SetOptStat(112211)
		histo.FindObject('stats').SetX1NDC(0.6)
		histo.FindObject('stats').SetX2NDC(0.9)
		histo.FindObject('stats').SetY1NDC(0.6)
		histo.FindObject('stats').SetY2NDC(0.9)
		ro.gPad.Update()

def SetDefaultFitStats(histo, fit, color=ro.kRed):
	fit.SetLineColor(color)
	histo.FindObject('stats').SetOptFit(1)
	histo.FindObject('stats').SetX1NDC(0.6)
	histo.FindObject('stats').SetX2NDC(0.9)
	histo.FindObject('stats').SetY1NDC(0.6)
	histo.FindObject('stats').SetY2NDC(0.9)

def SetDefault1DCanvasSettings(canvas):
	canvas.SetGridx()
	canvas.SetGridy()
	canvas.SetTicky()
	ro.gPad.Update()


def AddLineToStats(canvas, keys=['Mean_{Fit}'], values=[0], samplelinekey='Mean'):
	if canvas:
		ps = canvas.GetPrimitive('stats')
		ps.SetName('mystats')
		lol = ps.GetListOfLines()
		sampleline = ps.GetLineWith(samplelinekey)
		for it, key in enumerate(keys):
			line = ro.TLatex(0, 0, '{k} = {v:.4f}'.format(k=key, v=values[it]))
			line.SetTextAlign(sampleline.GetTextAlign())
			line.SetTextAngle(sampleline.GetTextAngle())
			line.SetTextColor(sampleline.GetTextColor())
			line.SetTextFont(sampleline.GetTextFont())
			line.SetTextSize(sampleline.GetTextSize())
			lol.Add(line)
		canvas.Modified()

def GetMinimumBranch(tree, bra, cut=''):
	tree.Draw('>>list{b}'.format(b=bra), cut)
	event_list = ro.gDirectory.Get('list{v}'.format(v=bra))
	tree.SetEventList(event_list)
	val = tree.GetMinimum(bra)
	tree.SetEventList(0)
	event_list.Delete()
	return val

def GetMaximumBranch(tree, bra, cut=''):
	tree.Draw('>>list{b}'.format(b=bra), cut)
	event_list = ro.gDirectory.Get('list{v}'.format(v=bra))
	tree.SetEventList(event_list)
	val = tree.GetMaximum(bra)
	tree.SetEventList(0)
	event_list.Delete()
	return val

def Correct_Path(path, times=2):
	abs_path = ''
	if path[0] == '~':
		abs_path += os.path.expanduser('~')
		abs_path += path[1:]
	elif os.path.isabs(path):
		abs_path += path
	else:
		abs_path += os.path.abspath(path)
	if times != 1:
		return Correct_Path(abs_path, 1)
	return abs_path

def RoundInt(n, nptype='int32'):
	val = np.floor(np.add(n, 0.5, dtype='f8'), dtype='f8').astype(nptype)
	if nptype.lower().startswith('i') or nptype.lower().startswith('ui'):
		return int(val)
	elif nptype.lower().startswith('f'):
		return float(val)
	return val

def TruncateFloat(num_float, resol=0):
	'''
	Truncates a float to the specified number of digits and decimals
	:param num_float: float to be truncated
	:param resol: gives the resolution. if 0, it means it does not matter.
	:return: truncated value to the specified resolution
	'''
	temp_float = num_float if resol == 0 else float(np.multiply(resol, RoundInt(np.divide(num_float, resol, dtype='f8')), dtype='f8'))
	return temp_float

def CheckEmptyBinsHisto(histo):
	emptybins = 0
	nbinsx = histo.GetNbinsX()
	nbinsy = histo.GetNbinsY()
	nbinsz = histo.GetNbinsZ()
	nbinsTot = nbinsx * nbinsy * nbinsz
	for bin in xrange(1, nbinsTot + 1):
		emptybins += 1 if histo.GetBinContent(bin) < 1 else 0
	return emptybins

def IsHistogramEmpty(histo):
	if histo.GetEntries() > 0:
		return False
	return True

def CheckFilledBinsHisto(histo):
	emptyBins = CheckEmptyBinsHisto(histo)
	nbinsx = histo.GetNbinsX()
	nbinsy = histo.GetNbinsY()
	nbinsz = histo.GetNbinsZ()
	nbinsTot = nbinsx * nbinsy * nbinsz
	return nbinsTot - emptyBins

def CheckBinningForFit(anaObject, histoKey, drawFunc, funcArgs, deltaPosInfuncArgs, fitParamsMin=6, truncateResol=0):
	good_binning = False
	while not good_binning:
		drawFunc(*funcArgs)
		if histoKey in anaObject.histo.keys():
			if anaObject.histo[histoKey] and not IsHistogramEmpty(anaObject.histo[histoKey]):
				filledBins = CheckFilledBinsHisto(anaObject.histo[histoKey])
				entries = anaObject.histo[histoKey].GetEntries()
				minBins = RoundInt(max(fitParamsMin, entries ** (2/5.)))
				maxBins = RoundInt(entries ** (10/17.))
				if minBins <= maxBins:
					delta = funcArgs[deltaPosInfuncArgs]
					good_binning = False
					if filledBins < minBins:
						delta /= 2.
						funcArgs = tuple([funcArgs[i] if i != deltaPosInfuncArgs else delta for i in xrange(len(funcArgs))])
					elif filledBins > maxBins:
						delta *= 1.7179869184
						funcArgs = tuple([funcArgs[i] if i != deltaPosInfuncArgs else delta for i in xrange(len(funcArgs))])
					else:
						good_binning = True
					delta = TruncateFloat(delta, truncateResol) if truncateResol != 0 else delta
				else:
					# only take into account the necessary filled Parameters for the fit
					delta = funcArgs[deltaPosInfuncArgs]
					good_binning = False
					if filledBins > minBins:
						delta *= 1.25
						funcArgs = tuple([funcArgs[i] if i != deltaPosInfuncArgs else delta for i in xrange(len(funcArgs))])
					elif filledBins <= minBins - 1:
						delta /= 1.5
						funcArgs = tuple([funcArgs[i] if i != deltaPosInfuncArgs else delta for i in xrange(len(funcArgs))])
					else:
						good_binning = True
					delta = TruncateFloat(delta, truncateResol) if truncateResol != 0 else delta

	return funcArgs[deltaPosInfuncArgs]

if __name__ == '__main__':
	print 'blaaaa'


