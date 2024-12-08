import os
import numpy as np
#import numpy.linalg as npl
import subprocess
import sys
import pandas as pd
from CANalyzer.load import Load
from CANalyzer.utilities import write_fchk, eig

class NaturalOrbitals(Load):
    def __init__(self, logfile, fchkfile, filename, states, displaywidth):
        super().__init__(logfile, fchkfile, filename, None, displaywidth)
        self.MO = None
        self.niorb = None
        self.naorb = None
        self.nstates = None
        self.states = states


    def start(self):
        super().start()

        if self.software == "GDV":
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



        elif self.software == "CQ":
            self.niorb = int(str(subprocess.check_output("grep 'Number of Inactive Core Orbitals:' " + self.logfile, shell=True)).split()[-1].split("\\")[0])
            self.naorb = int(str(subprocess.check_output("grep 'Number of Correlated Orbitals:' " + self.logfile, shell=True)).split()[-1].split("\\")[0])
            print(self.naorb, self.niorb)
            exit()

        # defining states of interest
        if self.states:
            soi = self.states.split(',')
            refine_soi = []
            for i in soi:
                try:
                    refine_soi.append(int(i))
                except:
                    bounds = i.split('-')
                    for j in range(int(bounds[0]), int(bounds[1]) + 1):
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
            if self.software == "GDV":
                real_pdm = self.readlog_matrix("1PDM Matrix (real):", self.naorb, self.naorb, instance=istate)
                imag_pdm = self.readlog_matrix("1PDM Matrix (imag):", self.naorb, self.naorb, instance=istate)
                pdm = real_pdm + 1j*imag_pdm
            elif self.software == "CQ":

            noon, no_transform = eig(pdm)

            acMO = self.MO[:, self.niorb:self.niorb+self.naorb]
            acnatorb = acMO @ no_transform
            natorb = np.copy(self.MO)
            natorb[:, self.niorb:self.niorb+self.naorb] = acnatorb

            no_transform_abs = np.abs(no_transform)
            np.savetxt(f'natorb-state{istate}-ascontrib.csv', no_transform_abs, delimiter=",")

            newfchk = f"natorb-state{istate}.fchk"
            matsize = self.nbasis * self.nbsuse * self.ncomp * self.ncomp
            write_fchk("Alpha MO coefficients", self.nri, natorb, matsize, self.fchkfile, newfchk)

            with open(f"{self.filename}-state{istate}.txt", 'a') as sys.stdout, pd.option_context('display.max_rows', None,
                                                                           'display.max_columns', None):
                print(f"Natural Orbital Occupation Numbers for State {istate}")
                for i in range(len(noon)):
                    nn = round(np.real(noon[i]), 6)
                    print(f'Natural Orbital {i}: {nn}')
                print(" ")










