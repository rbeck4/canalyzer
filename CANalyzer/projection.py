import numpy as np
import pandas as pd
from CANalyzer.load import Load
import sys

class Projection():
    def __init__(self, logfile, fchkfile1, fchkfile2, filename, displaywidth):
        self.filename = filename
        self.displaywidth = displaywidth
        self.fchk1 = Load(logfile, fchkfile1, None, None, None)
        self.fchk2 = Load(logfile, fchkfile2, None, None, None)
        self.overlap = None
        self.mo1a = None
        self.mo1b = None
        self.mo2a = None
        self.mo2b = None
        self.ctsc = None # fchk1 if projected onto fchk2


    def project(self):
        self.ctsc = np.real(np.conj(self.mo1a.T) @ self.overlap @ self.mo2a)
        if self.mo1b and self.mo2b:
            ctscb = np.real(np.conj(self.mo1b.T) @ self.overlap @ self.mo2b)
            self.ctsc += ctscb


    def print_project(self):
        with open(self.filename, 'w') as sys.stdout, pd.option_context('display.max_rows', None, 'display.max_columns', None):
            print("Projection of fchk1 MOs onto fchk2 MOs \n")
            nmo1, nmo2 = self.ctsc.shape
            for i in range(nmo1):
                for j in range(nmo2):
                    contrib = self.ctsc[j, i]
                    if np.abs(contrib) > 0.01:
                        print(f"Contribution of fchk1 MO {i+1} from fchk2 MO {i+1}: {contrib}\n")


    def start(self):
        self.fchk1.parse_constants()
        self.fchk1.parse_log()
        self.overlap = self.fchk1.read_overlap()
        self.mo1a, self.mo1b = self.fchk1.read_mo()

        self.fchk2.parse_constants()
        self.fchk2.parse_log()
        self.mo2a, self.mo2b = self.fchk2.read_mo()

        if self.fchk1.xhf[0] != self.fchk2.xhf[0]:
            raise Exception("Both jobs must be the same component.")




