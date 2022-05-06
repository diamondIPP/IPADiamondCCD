#!/usr/bin/env python
import time
from optparse import OptionParser
from argparse import ArgumentParser
from Settings_Caen import Settings_Caen
from Utils import *
from CCD_Caen import CCD_Caen

trig_rand_time = 0.2
wait_time_hv = 7


class VoltageScan:
    def __init__(self, infile='None', vini=5, vend=10, vstep=5, lista=[], listn=[], timebla=1, verbose=False):
        print('Starting CCD program ...')
        self.infile = infile
        self.verb = verbose
        self.vini = vini
        self.vend = vend
        self.vstep = vstep
        self.voltages = np.linspace(self.vini, self.vend,
                                    int(round(float(self.vend - self.vini) / float(self.vstep)) + 1),
                                    dtype='int32') if lista == [] else np.array(lista, 'int32')
        self.time_sleep = timebla
        self.settings = Settings_Caen(self.infile, self.verb)
        self.settings.read_input_file()
        self.settings.Get_Calibration_Constants()
        self.settings.SetOutputFiles()
        self.wd = os.getcwd()

        self.p = {volt: None for volt in self.voltages}
        self.time0 = time.time()
        self.num_events = self.settings.num_events * np.ones(len(self.voltages), 'int32') if listn == [] else np.array(
            listn, 'int32')
        print(self.num_events)

    def DoVoltageScan(self):

        for it, volt in enumerate(self.voltages):
            self.ResetSettings()
            self.settings.bias = float(volt)
            self.settings.num_events = int(self.num_events[it])
            self.settings.Get_Calibration_Constants()
            self.settings.SetOutputFiles()

            self.p[volt] = CCD_Caen(settings=self.settings)
            self.p[volt].StartHVControl()
            self.p[volt].AdjustBaseLines()
            self.p[volt].SavePickles()
            print(f'Waiting {self.time_sleep} seconds before starting run...', end=' ')
            sys.stdout.flush()
            time.sleep(self.time_sleep)
            print('Done')
            written_events = self.p[volt].GetData()
            if self.p[volt].stop_run:
                print('Run stopped because current is too high')
            self.settings.num_events = written_events
            self.p[volt].SavePickles()
            self.p[volt].settings.MoveBinaryFiles()
            self.p[volt].settings.RenameDigitiserSettings()
            self.p[volt].CloseHVClient()
            if not self.p[volt].settings.simultaneous_conversion:
                self.p[volt].CreateRootFile(files_moved=True)
                while self.p[volt].pconv.poll() is None:
                    time.sleep(3)
                self.p[volt].CloseSubprocess('converter', stdin=False, stdout=False)

            print(f'\nFinished voltage scan for {volt} Volts :)\n')
            sys.stdout.write('\a\a\a')
            sys.stdout.flush()

            self.p[volt] = None

        print('Finished all voltage scan :D')
        sys.stdout.write('\a\a\a')
        sys.stdout.flush()

    def ResetSettings(self):
        self.settings = Settings_Caen(self.infile, self.verb)
        self.settings.read_input_file()
        self.settings.Get_Calibration_Constants()
        self.settings.SetOutputFiles()


def main():
    parser = ArgumentParser()
    parser.add_argument('-i', '--infile', type=str, help='Input configuration file. e.g. CAENRunConfig.cfg')
    parser.add_argument('--start', type=int, help='Starting scan voltage', default=0)
    parser.add_argument('--stop', type=int, help='Stopping scan voltage', default=0)
    parser.add_argument('--step', type=int, help='Voltage step between scans', default=0)
    parser.add_argument('-v', '--verbose', dest='verb', help='Toggle verbosity', action='store_true')
    parser.add_argument('-a', '--auto', help='Toggles between automatic conversion and analysis afterwards',
                        action='store_true')
    parser.add_argument('-l', '--list', dest='voltages', type=str,
                        help='List of bias voltages to analyse. It overrides "--step, --start, --stop option".')
    parser.add_argument('-n', '--numevts', dest='events', type=str,
                        help='List of number of events to collect at each bias voltage.')
    parser.add_argument('-t', '--time', dest='timebetween', type=int, help='time in between measurements in seconds.',
                        default=120)

    args = parser.parse_args()
    infile = str(args.infile)
    auto = args.auto
    verb = args.verb
    vini = args.start
    vend = args.stop
    vstep = args.step
    voltages = [float(v) for v in args.voltages.split()]
    n_events = [int(n) for n in args.events.split()]
    print(voltages)
    timebla = int(args.timebetween)
    print(n_events)
    if n_events and voltages:
        if len(n_events) != len(voltages):
            ExitMessage('The lists of voltages and event number have different dimensions.')
    vs = VoltageScan(infile, vini, vend, vstep, voltages, n_events, timebla, verb)
    if auto:
        vs.DoVoltageScan()
    return vs


if __name__ == '__main__':
    vs = main()
