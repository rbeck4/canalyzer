"""
This is meant to be called as a library from a Jupyter notebook.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
import linecache
from CANalyzer.mo import MO

class Plot(MO):
    def __init__(self, logfile, fchkfile, groups, orbital_decomposition, separate_ml=False):
        super().__init__(logfile, fchkfile, None, groups, None, None)
        self.orbital_decomposition = orbital_decomposition
        self.oscstr = None
        self.energy = None
        self.jobtype = None
        self.separate_ml = separate_ml
        self.nstates = None


    def start(self):
        super().start()

        while self.jobtype not in ["TDDFT", "CI", "EOMCC"]:
            self.jobtype = input("TDDFT, CI, or EOMCC?")

        if self.software == "GDV":
            if self.jobtype == "TDDFT":
                startline = int(
                    str(subprocess.check_output("grep -n 'Excited State' " + self.logfile, shell=True)).split()[
                        0].split("'")[-1][:-1])
                self.nstates = int(str(subprocess.check_output("grep -n 'Excited State' " + self.logfile, shell=True)).split()[-8][:-1])

                statecounter = 1
                offswitch = True
                while offswitch:
                    parseline = linecache.getline(self.logfile, startline)
                    print(parseline)
                    if "Excited State" in parseline:
                        statecounter += 1
                    if statecounter > self.nstates:
                        if "->" not in parseline and "<-" not in parseline and "Excited State" not in parseline:
                            offswitch = False
                    startline += 1

        elif self.software == "CQ":
            raise Exception("NYI")
