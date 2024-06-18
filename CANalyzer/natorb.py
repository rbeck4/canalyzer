import os
import numpy as np
import numpy.linalg as npl
import subprocess
import sys
import pandas as pd
from CANalyzer.load import Load
from CANalyzer.utilities import write_fchk

class NaturalOrbitals(Load):
    def __init__(self, logfile, fchkfile, filename, states, displaywidth):
        Load.__init__(self, logfile, fchkfile, filename, None, displaywidth)
        self.MO = None
        self.niorb = None
        self.naorb = None
        self.nstates = None
        self.states = states


    def start(self):
        self.parse_constants()
        self.parse_log()

        # getting number of active orbitals
        line = str(subprocess.check_output("grep 'NAOrb=' " + self.logfile, shell=True)).split(" ")
        refine_line = []
        for l in line:
            try:
                refine_line.append(int(l))
            except:
                pass
        self.niorb = refine_line[1]
        self.naorb = refine_line[2]

        # defining states of interest
        if self.states:
            soi = self.states.split(',')
            refine_soi = []
            for i in soi:
                try:
                    refine_soi.append(int(i))
                except:
                    bounds = i.split('-')
                    for j in range(int(bounds[0]), int(bounds[1])+1):
                        refine_soi.append(j)
            refine_soi.sort()
            self.states = refine_soi
        else:
            self.states = [1]

        # load orbitals
        self.MO, void = self.read_mo(separate=False)


    def compute_natorb(self):
        pdm = None
        for istate in self.states:
            real_pdm = self.readlog_matrix("1PDM Matrix (real):", self.naorb, self.naorb, instance=istate)
            imag_pdm = self.readlog_matrix("1PDM Matrix (imag):", self.naorb, self.naorb, instance=istate)
            pdm = real_pdm + 1j*imag_pdm
            noon, no_tranform = npl.eig(pdm)

            acMO = self.MO[:, self.niorb:self.niorb+self.naorb]
            acnatorb = acMO @ no_tranform
            natorb = np.copy(self.MO)
            natorb[:, self.niorb:self.niorb+self.naorb] = acnatorb

            newfchk = f"natorb-state{istate}.fchk"
            matsize = self.nbasis * self.nbsuse * self.ncomp * self.ncomp
            write_fchk("Alpha MO coefficients", self.nri, natorb, matsize, self.fchkfile, newfchk)

            with open(self.filename, 'w') as sys.stdout, pd.option_context('display.max_rows', None,
                                                                           'display.max_columns', None):
                print(f"Natural Orbital Occupation Numbers for State {istate}")
                for i in range(len(noon)):
                    print(f'Natural Orbital {i}: {np.real(noon[i])}')
                print(" ")










