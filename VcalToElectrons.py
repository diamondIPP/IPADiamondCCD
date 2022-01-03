#!/usr/bin/env python
import os
import shutil
import sys
import time
from ConfigParser import ConfigParser
from optparse import OptionParser

import ROOT as ro
import numpy as np
import scipy.constants as spc
from uncertainties import ufloat
import cPickle as pickle

from Settings_Caen import Settings_Caen
from Utils import *

#
# vcal------------r2--------------calCermic--------JFET
#            |          |     |_____capLin_____|
#            r1         r3
#            |          |
#            |          |
#           GND        GND

class VcalToElectrons:
    def __init__(self, configfile=''):
        self.configPath = Correct_Path(configfile) if configfile != '' else ''
        # in Ohms
        self.r2 = ufloat(467, 4.036)
        self.r3 = ufloat(4.505, 0.032)
        self.vgain = ufloat(0, 0)
        self.UpdateVoltageGain()
        # in pF
        self.cap_ceramic = ufloat(1.846, 0.002)
        self.cap_line = ufloat(2.26, 0.05)
        self.cap_in = ufloat(0, 0)
        self.UpdateInputCapacitance()

        self.LoadFromConfigFile()

    def UpdateVoltageGain(self):
        self.vgain = np.divide(self.r3, np.add(self.r2, self.r3))

    def UpdateInputCapacitance(self):
        self.cap_in = np.add(self.cap_line, self.cap_ceramic)

    def LoadFromConfigFile(self):
        if self.configPath != '':
            if os.path.isfile(self.configPath):
                confparser = ConfigParser()
                confparser.read(self.configPath)
                if confparser.has_section('COMPONENTS'):
                    if confparser.has_option('COMPONENTS', 'r2'):
                        r2_temp = confparser.getfloat('COMPONENTS', 'r2')
                        r2_sigma_temp = confparser.getfloat('COMPONENTS', 'r2_sigma') if confparser.has_option('COMPONENTS', 'r2_sigma') else np.multiply(r2_temp, 0.01)
                        self.r2 = ufloat(r2_temp, r2_sigma_temp)
                    if confparser.has_option('COMPONENTS', 'r3'):
                        r3_temp = confparser.getfloat('COMPONENTS', 'r3')
                        r3_sigma_temp = confparser.getfloat('COMPONENTS', 'r3_sigma') if confparser.has_option('COMPONENTS', 'r3_sigma') else np.multiply(r3_temp, 0.01)
                        self.r3 = ufloat(r3_temp, r3_sigma_temp)
                    if confparser.has_option('COMPONENTS', 'cap_ceramic'):
                        c_temp = confparser.getfloat('COMPONENTS', 'cap_ceramic')
                        c_sigma_temp = confparser.getfloat('COMPONENTS', 'cap_ceramic_sigma') if confparser.has_option('COMPONENTS', 'cap_ceramic_sigma') else np.multiply(c_temp, 0.1)
                        self.cap_ceramic = ufloat(c_temp, c_sigma_temp)
                    if confparser.has_option('COMPONENTS', 'cap_line'):
                        cl_temp = confparser.getfloat('COMPONENTS', 'cap_line')
                        cl_sigma_temp = confparser.getfloat('COMPONENTS', 'cap_line_sigma') if confparser.has_option('COMPONENTS', 'cap_line_sigma') else np.multiply(cl_temp, 0.1)
                        self.cap_line = ufloat(cl_temp, cl_sigma_temp)

        self.UpdateVoltageGain()
        self.UpdateInputCapacitance()

    def Q_in_e_from_mV(self, vcal=0, vcal_sigma=0):
        vtemp = np.multiply(self.vgain, np.divide(ufloat(vcal, vcal_sigma), 1000))
        return np.divide(np.multiply(vtemp, np.multiply(self.cap_in, 1e-12)), spc.elementary_charge)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c', '--config', dest='config', default='', type='string', help='configuration file with the values for the calibration circuit. Default is \'\'')
    (options, args) = parser.parse_args()
    configfile = str(options.config)

    v2e = VcalToElectrons(configfile)
