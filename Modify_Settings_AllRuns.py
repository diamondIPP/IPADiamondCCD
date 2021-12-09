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


from Modify_Settings_Caen import Modify_Pickles_Caen

class Modify_Settings_AllRuns:
    def __init__(self, runsdir='', Triggerval = 0):
        self.runsdir = Correct_Path(runsdir)
        runstemp = glob.glob('{d}/*'.format(d=self.runsdir))
        self.runs = [runi for runi in runstemp if os.path.isdir(runi)]
        self.TriggerVal = Triggerval



    def LoopRunsTrigger(self):
        for run in self.runs:
            print '\nAnalysing run:', run
            Modif = Modify_Pickles_Caen(run)
            Modif.ModifyTriggerThreshold(self.TriggerVal)
            Modif.SavePickles()

def main():
    parser = OptionParser()
    parser.add_option('-d', '--runsdir', dest='runsdir', type='str', default='', help='path to folder containing all the run folders to modify')
    parser.add_option('-v', '--value', dest='value', type='int', default=0, help='Value to set Trigger Threshold')

    (options, args) = parser.parse_args()
    runsdir = str(options.runsdir)
    value =  int(options.value)
    ms = Modify_Settings_AllRuns(runsdir, value)
    return ms

if __name__ == '__main__':
	ms = main()






