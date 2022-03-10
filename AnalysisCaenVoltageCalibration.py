#!/usr/bin/env python
import glob
from optparse import OptionParser
from argparse import ArgumentParser
import ROOT as ro
import pickle as pickle
from Utils import *
import dill

# accuracy and resolution of reference multimeter UT803. Change accordingly if using another reference device. All values are given in mV
reference_ranges = [600, 6000, 60000, 600000, 1000000]
reference_accuracy = {600: {'percent': 0.6, 'digits': 0.2}, 6000: {'percent': 0.3, 'digits': 2}, 60000: {'percent': 0.3, 'digits': 20}, 600000: {'percent': 0.3, 'digits': 200}, 6000000: {'percent': 0.5, 'digits': 3000}}
fit_method = ('Minuit2', 'Migrad', )

class AnalysisCaenVoltageCalibration:
	def __init__(self, directory='.'):
		print('Starting Caen Voltage Calibration Analysis ...')
		self.utils = Utils()
		self.graphs = {}
		self.canvas = {}
		self.fits = {}
		self.cal_pickle = None
		self.runs = []
		self.vcals = []
		self.vcals_valid = []
		self.vcals_uncertainty = []
		self.adcs_mean = []
		self.adcs_std = []
		self.cal_pickle_name = ''
		self.inDir = Correct_Path(directory)
		if not os.path.isdir(self.inDir):
			ExitMessage('The given directory does not exist. Exiting!', os.EX_DATAERR)
		self.cal_pickles = self.LookForCalPickles()
		if len(self.cal_pickles) > 0:
			self.LoadPickle()

		if not self.cal_pickle:
			self.DoAsFirstTime()

		self.working_dir = ''
		self.working_vcal = 0
		self.working_vcal_un = 0
		self.working_settings_path = ''
		self.working_settings = None
		self.working_channel_path = ''
		self.working_channel = None
		self.working_signal_path = ''
		self.working_signal = None
		self.working_struct_fmt = ''
		self.working_struct_len = 0
		self.working_bits_adc = 14
		self.working_caen_ch = 3
		self.working_caen_ch_dc_off_percent = 0
		self.working_num_events = 0
		self.working_adcs = []
		self.working_total_points = 0

	def DoAsFirstTime(self):
		self.cal_pickle = None
		self.runs = []
		self.vcals = []
		self.vcals_valid = []
		self.vcals_uncertainty = []
		self.runs = glob.glob(self.inDir + '/*mV')
		if len(self.runs) < 2: ExitMessage('Can\'t make calibration with only ' + str(len(self.runs)) + ' runs. Exiting!', os.EX_USAGE)
		self.runs.sort(key=lambda x: float(x.split('_')[-1].split('mV')[0]))
		self.vcals = [float(x.split('_')[-1].split('mV')[0]) for x in self.runs]
		if len(self.vcals) != len(self.runs): ExitMessage(f'There was an error. Check the runs inside {self.inDir}', os.EX_DATAERR)
		self.EstimateVcalsUncertainty()

	def LookForCalPickles(self):
		pick_files = glob.glob(self.inDir + '/*.cal')
		return pick_files

	def LoadPickle(self, pos=0):
		dill._dill._reverse_typemap["ObjectType"] = object
		if pos < len(self.cal_pickles):
			self.cal_pickle = pickle.load(open(self.cal_pickles[pos], 'rb'), encoding="latin1")
			self.cal_pickle_name = self.cal_pickle['file_name']
			self.vcals = self.cal_pickle['vcals']
			self.vcals_uncertainty = self.cal_pickle['vcals_sigma']
			self.vcals_valid = self.cal_pickle['vcals_valid']
			self.adcs_mean = self.cal_pickle['adcs_mean']
			self.adcs_std = self.cal_pickle['adcs_std']
			self.adcs_std = self.cal_pickle['adcs_std']
			self.fits['ADC_Voltage_cal'] = ro.TF1('fit_' + 'ADC_Voltage_cal', 'pol1', 0, 2**14 -1)
			self.fits['ADC_Voltage_cal'].SetParameter(0, self.cal_pickle['fit_p0'])
			self.fits['ADC_Voltage_cal'].SetParError(0, self.cal_pickle['fit_p0_error'])
			self.fits['ADC_Voltage_cal'].SetParameter(1, self.cal_pickle['fit_p1'])
			self.fits['ADC_Voltage_cal'].SetParError(1, self.cal_pickle['fit_p1_error'])
			self.fits['ADC_Voltage_cal'].SetChisquare(self.cal_pickle['fit_chi2'])
			self.fits['ADC_Voltage_cal'].SetNDF(self.cal_pickle['fit_ndf'])
			print('Loaded pickle:', self.cal_pickles[pos], '. Run LoadPickle again with another argument if you want to load another pickle')

	def EstimateSystematicUncertainty(self, reading, reading_percent=0.6, reading_fixed=0.2):
		resol = reading_fixed / float(str(reading_fixed).split('.')[-1].strip('0')) if reading_fixed != 0 else ExitMessage('reading_fixed cannot be 0!. Exiting', os.EX_PROTOCOL)
		system_un = TruncateFloat(abs(reading) * reading_percent / 100. + reading_fixed, resol)
		return system_un

	def EstimateVcalsUncertainty(self):
		for vcal in self.vcals:
			pos = np.less_equal(abs(vcal), reference_ranges).argmax()
			pos = -1 if pos == 0 and np.all(np.greater_equal(abs(vcal), reference_ranges)) else pos
			range_used = reference_ranges[pos]
			ref_acc = reference_accuracy[range_used]
			self.vcals_uncertainty.append(self.EstimateSystematicUncertainty(vcal, ref_acc['percent'], ref_acc['digits']))

	def LoopRuns(self):
		print('Looping over all vcals:')
		self.utils.CreateProgressBar(len(self.vcals))
		self.utils.bar.start()
		for pos in range(len(self.vcals)):
			self.working_adcs = []
			self.LoadRun(pos)
			self.LoadEvents()
			self.working_adcs = [val for sublist in self.working_adcs for val in sublist]
			self.working_total_points = len(self.working_adcs)
			# check if the adcs were saturated more than 0.5%
			valid_vcal = (np.equal(self.working_adcs, 2 ** self.working_bits_adc - 1).sum()+ np.equal(self.working_adcs, 0).sum()) / float(self.working_total_points) < 0.005
			self.vcals_valid.append(valid_vcal)
			self.adcs_mean.append(np.mean(self.working_adcs, dtype='f8') if valid_vcal else 0)
			self.adcs_std.append(np.std(self.working_adcs, dtype='f8') if valid_vcal else 0)
			self.utils.bar.update(pos + 1)
		self.utils.bar.finish()
		print('Finished with all vcals :)')

	def LoadEvents(self):
		unpack_fmt = '@' + str(self.working_num_events * self.working_settings.points) + 'H'
		self.working_signal.seek(0, 0)
		print('Types:', type(self.working_struct_len), type(self.working_num_events))
		tempdata = self.working_signal.read(self.working_struct_len * self.working_num_events)
		tempstruct = struct.Struct(unpack_fmt).unpack_from(tempdata)
		self.working_adcs.append(tempstruct)

	def LoadRun(self, pos):
		self.CloseBinary()
		self.working_dir = self.runs[pos]
		self.working_vcal = self.vcals[pos]
		self.working_vcal_un = self.vcals_uncertainty[pos]

		temp_list = glob.glob(self.working_dir + '/*.settings')
		if len(temp_list) != 1: ExitMessage(f'There should be one and only one settings file pickle in {self.working_dir}', os.EX_DATAERR)
		self.working_settings_path = temp_list[0]
		# convert py2 to py3
		dill._dill._reverse_typemap["ObjectType"] = object
		self.working_settings = pickle.load(open(self.working_settings_path, 'rb'), encoding="latin1")
		print('Working Settings:', self.working_settings.struct_fmt, self.working_settings.struct_len, self.working_settings.dig_bits)
		self.working_struct_fmt = self.working_settings.struct_fmt
		self.working_struct_len = self.working_settings.struct_len
		self.working_bits_adc = self.working_settings.dig_bits

		temp_list = glob.glob(self.working_dir + '/*.signal_ch')
		if len(temp_list) != 1:
			ExitMessage(f'There should be one and only one signal_ch file pickle in {self.working_dir}', os.EX_DATAERR)
		self.working_channel_path = temp_list[0]
		self.working_channel = pickle.load(open(self.working_channel_path, 'rb'), encoding="latin1")
		self.working_caen_ch = self.working_channel.ch
		self.working_caen_ch_dc_off_percent = self.working_channel.dc_offset_percent

		temp_list = glob.glob(self.working_dir + '/*signal.dat')
		if len(temp_list) != 1:
			ExitMessage(f'There should be one and only one signal.dat binary file in {self.working_dir}', os.EX_DATAERR)
		self.working_signal_path = temp_list[0]
		self.LoadBinary()

	def LoadBinary(self):
		if os.path.isfile(self.working_signal_path):
			self.working_signal = open(self.working_signal_path, 'rb')
			self.working_num_events = os.path.getsize(self.working_signal_path) // self.working_settings.struct_len

	def CloseBinary(self):
		if self.working_signal:
			if not self.working_signal.closed:
				self.working_signal.close()
				self.working_signal = None

	def CreateResultsGraph(self, name='ADC_Voltage_cal'):
		ypoints = np.multiply(np.extract(self.vcals_valid, self.vcals).astype('f8'), 0.001, dtype='f8')
		ypoints_errs = np.multiply(np.extract(self.vcals_valid, self.vcals_uncertainty).astype('f8'), 0.001, dtype='f8')
		xpoints = np.extract(self.vcals_valid, self.adcs_mean).astype('f8')
		xpoints_errs = np.extract(self.vcals_valid, self.adcs_std).astype('f8')
		npoints = int(np.sum(self.vcals_valid))
		self.CheckExistingGraph(name)
		self.graphs[name] = ro.TGraphErrors(npoints, xpoints, ypoints, xpoints_errs, ypoints_errs)
		self.graphs[name].SetNameTitle('g_' + name, 'g_' + name)
		self.graphs[name].SetMarkerStyle(7)
		self.graphs[name].SetMarkerColor(ro.kBlack)
		self.graphs[name].GetXaxis().SetTitle('adc')
		self.graphs[name].SetLineColor(ro.kBlack)
		self.graphs[name].GetXaxis().SetTitle('ADC')
		self.graphs[name].GetYaxis().SetTitle('Voltage [V]')

	def ExcludeVcal(self, vcal):
		pos = self.vcals.index(vcal)
		self.vcals_valid[pos] = False

	def IncludeExcludedVcal(self, vcal):
		pos = self.vcals.index(vcal)
		self.vcals_valid[pos] = True

	def DrawGraph(self, name='ADC_Voltage_cal'):
		self.CreateResultsGraph(name)
		self.CheckExistingCanvas(name)
		self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
		self.graphs[name].Draw('AP')
		SetDefault1DCanvasSettings(self.canvas[name])
		self.FitGraph(name)
		self.canvas[name].SaveAs(f'{self.inDir}/{name}_{self.working_caen_ch}_{self.working_caen_ch_dc_off_percent}.png')
		self.canvas[name].SaveAs(f'{self.inDir}/{name}_{self.working_caen_ch}_{self.working_caen_ch_dc_off_percent}.root')

	def FitGraph(self, name='ADC_Voltage_cal'):
		if name in list(self.graphs.keys()):
			if self.graphs[name]:
				ro.Math.MinimizerOptions.SetDefaultMinimizer(*fit_method)
				ro.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(1000000)
				ro.Math.MinimizerOptions.SetDefaultTolerance(0.00001)
				ro.gStyle.SetOptFit(1111)
				func = ro.TF1('fit_' + name, 'pol1', 0, 2**14 -1)
				func.SetLineColor(ro.kRed)
				func.SetNpx(int(10 * (2**14 - 1)))
				self.graphs[name].Fit('fit_' + name, 'Q0', '', 0, 2**14 -1)
				if func.GetProb() < 0.9:
					self.graphs[name].Fit('fit_' + name, 'Q0', '', 0, 2**14 -1)
				if func.GetProb() < 0.9:
					self.graphs[name].Fit('fit_' + name, 'Q0', '', 0, 2**14 -1)
				self.fits[name] = func
				# self.fits[name].Draw('same')

	def CheckExistingCanvas(self, name='ADC_Voltage_cal'):
		if name in list(self.canvas.keys()):
			if self.canvas[name]:
				self.canvas[name].Close()
				del self.canvas[name]
				
	def CheckExistingGraph(self, name='ADC_Voltage_cal'):
		if name in list(self.graphs.keys()):
			if self.graphs[name]:
				del self.graphs[name]

	def CheckExistingFits(self, name='ADC_Voltage_cal'):
		if name in list(self.fits.keys()):
			if self.fits[name]:
				del self.fits[name]

	def FillPickle(self, name='ADC_Voltage_cal'):
		self.cal_pickle = {'file_name': self.cal_pickle_name,
		                   'vcals': self.vcals,
		                   'vcals_sigma': self.vcals_uncertainty,
		                   'vcals_valid': self.vcals_valid,
		                   'adcs_mean': self.adcs_mean,
		                   'adcs_std': self.adcs_std,
		                   'fit_p0': self.fits[name].GetParameter(0),
		                   'fit_p0_error': self.fits[name].GetParError(0),
		                   'fit_p1': self.fits[name].GetParameter(1),
		                   'fit_p1_error': self.fits[name].GetParError(1),
		                   'fit_prob': self.fits[name].GetProb(),
		                   'fit_chi2': self.fits[name].GetChisquare(),
		                   'fit_ndf': self.fits[name].GetNDF()
		                   }

	def SavePickle(self, name='ADC_Voltage_cal', overwrite=False):
		self.cal_pickle_name = f'adc_cal_{self.working_caen_ch}_{self.working_caen_ch_dc_off_percent}.cal'
		if not self.cal_pickle:
			self.FillPickle(name)
		if os.path.isfile(f'{self.inDir}/{self.cal_pickle_name}'):
			if not overwrite:
				print('The file', self.cal_pickle_name, 'already exists in', self.inDir)
				return
		with open(f'{self.inDir}/{self.cal_pickle_name}', 'wb') as fpickle:
			pickle.dump(self.cal_pickle, fpickle, pickle.HIGHEST_PROTOCOL)
		print('Saved calibration pickle', self.cal_pickle_name, 'in', self.inDir)


if __name__ == '__main__':
	parser = ArgumentParser()
	parser.add_argument('-d', '--inDir', dest='inDir', default='.', type=str,
						help='Directory containing the subdirectories with different voltages run files')
	parser.add_argument('-a', '--automatic', dest='auto', default=False,
						help='Toggles automatic basic analysis', action='store_true')
	parser.add_argument('-o', '--overwrite', dest='overwrite', default=False,
						help='Toggles overwriting of the analysis tree', action='store_true')

	args = parser.parse_args()
	directory = str(args.inDir)
	auto = bool(args.auto)
	overwrite = bool(args.overwrite)

	ana = AnalysisCaenVoltageCalibration(directory)

	if auto:
		if not ana.cal_pickle:
			ana.LoopRuns()
			ana.CreateResultsGraph()
			ana.FitGraph()
			ana.SavePickle()
	# return ana
	# ana = main()
