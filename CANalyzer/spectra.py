"""
This is meant to be called as a library from a Jupyter notebook.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
import linecache
from itertools import chain

import CANalyzer.utilities
from CANalyzer.mo import MO

class Spectra(MO):
    def __init__(self, logfile, fchkfile, groups, orbital_decomposition, separate_ml=False):
        super().__init__(logfile, fchkfile, None, groups, None, None)
        self.orbital_decomposition = orbital_decomposition
        self.oscstr = []
        """
            TDDFT: self.oscstr is stored in excited state ordering
            CI: self.oscstr is stored as a list ordered in a compound index that represents (from_state, to_state) 
                where the from_state is the fast-running index; compound index can be obtained using self.ci_os_index() 
        """
        self.energy = []
        """
            TDDFT: self.energy stores excitation energies in eV
            CI: self.energy stores raw energies in Hartree
        """
        self.jobtype = None
        self.separate_ml = separate_ml
        self.nstates = None
        """
            TDDFT: self.nstates is the number of excited states
            CI: self.nstates is the number of total states
        """
        self.orbital_contributions = None
        """
            TDDFT: self.orbital_contributions is [[(from_mo, to_mo, contribution)] for each state]
            CI: self.orbital_contributions is nstates x naorb matrix where elements are orbital occupation numbers
        """
        self.decomp_os = None
        #self.spin = None
        #self.plots = []
        #"""
        #    List of different plots using this data
        #"""


    def start(self):
        super().start()

        while self.jobtype not in ["TDDFT", "CI", "EOMCC"]:
            self.jobtype = input("TDDFT, CI, or EOMCC? ")

            if self.jobtype == "TDDFT":
                self.parse_tddft()
            elif self.jobtype == "CI":
                self.parse_ci()


    def parse_tddft(self):
        if self.software == "GDV":
            startline = int(
                str(subprocess.check_output("grep -n 'Excited State' " + self.logfile, shell=True)).split()[
                    0].split("'")[-1][:-1])
            self.nstates = int(
                str(subprocess.check_output("grep -n 'Excited State' " + self.logfile, shell=True)).split()[-8][:-1])
            raw_orbital_contributions = [[] for i in range(self.nstates)]
            self.orbital_contributions = [[] for i in range(self.nstates)]

            statecounter = 1
            current_state = 0
            offswitch = True
            while offswitch:
                parseline = linecache.getline(self.logfile, startline)
                if "Excited State" in parseline:
                    splitline = parseline.split()
                    current_state = int(splitline[2][:-1])
                    self.energy.append(float(splitline[4]))
                    self.oscstr.append(float(splitline[8].split("=")[1]))
                    # state_spin = splitline[9].split("=")[1]
                    statecounter += 1
                elif "->" in parseline and self.orbital_decomposition:
                    splitline = parseline.split()
                    from_orbital = int(splitline[0])
                    to_orbital = int(splitline[2])
                    contribution = float(splitline[3])
                    raw_orbital_contributions[current_state - 1].append([from_orbital, to_orbital, contribution])
                elif '<-' in parseline and self.orbital_decomposition:
                    splitline = parseline.split()
                    from_orbital = int(splitline[0])
                    to_orbital = int(splitline[2])
                    contribution = float(splitline[3])
                    raw_orbital_contributions[current_state - 1].append([from_orbital, to_orbital, contribution])

                if statecounter > self.nstates:
                    if "->" not in parseline and "<-" not in parseline and "Excited State" not in parseline:
                        offswitch = False
                startline += 1

            if self.orbital_decomposition:
                for i in range(len(raw_orbital_contributions)):
                    state_contribution_list = raw_orbital_contributions[i]
                    done_pairs = []
                    unnormalized_contributions = []
                    for pair_index in range(len(state_contribution_list)):
                        current_pair = (state_contribution_list[pair_index][0], state_contribution_list[pair_index][1])
                        contribution = np.abs(state_contribution_list[pair_index][2])

                        if current_pair not in done_pairs:
                            unnormalized_contributions.append(contribution)
                            done_pairs.append(current_pair)

                        elif current_pair in done_pairs:
                            done_pair_index = done_pairs.index(current_pair)
                            contribution += unnormalized_contributions[done_pair_index]
                            unnormalized_contributions[done_pair_index] = contribution

                    unnormalized_total_contribution = sum(unnormalized_contributions)
                    normalized_contributions = [x / unnormalized_total_contribution for x in unnormalized_contributions]
                    self.orbital_contributions[i] = [(done_pairs[j][0], done_pairs[j][1], normalized_contributions[j]) for j
                                                     in range(len(done_pairs))]
        elif self.software == "CQ":
            raise Exception("NYI")


    def parse_ci(self):
        if self.software == "GDV":
            energyline = str(subprocess.check_output("grep -n '(Hartree)' " + self.logfile, shell=True))
            energyline_split = energyline.split("\\n")
            self.energy = [float(x.split("(Hartree):")[-1]) for x in energyline_split[:-1]]
            self.nstates = int(energyline_split[-2].split(":")[2].split()[0])
            pdmdiag_list = []  # just 1PDM diagonals in MO basis
            for i in range(self.nstates):
                pdmdiag = self.readlog_matrix("For Simplicity", self.nstates, 1, instance=i+1).flatten()
                pdmdiag_list.append(pdmdiag)
            self.orbital_contributions = np.array(pdmdiag_list)

            self.oscstr = []
            oslines = str(subprocess.check_output("grep 'Oscillator Strength For States' " + self.logfile, shell=True)).split("\\n")[:-1]
            oslines[0] = oslines[0].split("'")[-1]
            for i in range(len(oslines)):
                oscstr = float(oslines[i].split("f=")[-1])
                self.oscstr.append(oscstr)

        elif self.software == "CQ":
            raise Exception("NYI")


    def ci_os_index(self, from_state, to_state):
        states_per_layers_before = [self.nstates - i for i in range(1, from_state)]
        skip = sum(states_per_layers_before)
        return skip + to_state - from_state - 1


    def purge(self, threshold=1E-7):
        dark_states = []
        for i in range(self.nstates):
            if self.oscstr[i] < threshold:
                dark_states.append(i)
        for i in sorted(dark_states, reverse=True):
            del self.oscstr[i]
            del self.energy[i]
            if self.jobtype == "TDDFT":
                del self.orbital_contributions[i]
            elif self.jobtype == "CI":
                self.orbital_contributions[i] = np.delete(self.orbital_contributions, i, axis=0)
        ndark_state = len(dark_states)
        self.nstates = self.nstates - ndark_state


    def make_spectra(self, xstart, xend, os, broadening, lineshape=CANalyzer.utilities.lorentzian, npoints=5000):
        line_list = []
        xspace = np.linspace(xstart, xend, npoints)
        for i in range(self.nstates):
            root = self.energy[i]
            oscstr = os[i]
            if root >= xstart and root <= xend:
                line_list.append(lineshape(root, broadening, oscstr, xspace))
        spectrum = np.zeros(npoints)
        for i in line_list:
            spectrum += i
        return xspace, spectrum


    def decompose_os(self, occupied_groups, virtual_groups):
        """
        Decomposes oscillator strengths based on contribution from excitations between defined groups.
        :param occupied_groups: Input groups of occupied orbitals as dict {"Name1" : range(start, end + 1)}
        :param virtual_groups:
        """
        occupied_groups_keys = occupied_groups.keys()
        virtual_groups_keys = virtual_groups.keys()

        occupied_superrange = chain.from_iterable((occupied_groups[x] for x in occupied_groups_keys))
        virtual_superrange = chain.from_iterable((virtual_groups[x] for x in virtual_groups_keys))

        min_specified = min(list(occupied_superrange))
        max_specified = max(list(virtual_superrange))

        occupied_groups["Other Occ"] = range(1, min_specified + 1)
        virtual_groups["Other Virt"] = range(max_specified + 1, self.nbsuse)

        occupied_groups_keys = occupied_groups.keys()
        virtual_groups_keys = virtual_groups.keys()

        self.ex_categories_keys = [(x, y) for x in occupied_groups_keys for y in virtual_groups_keys]
        self.ncat = len(self.ex_categories_keys)
        self.decomp_os = np.zeros((self.nstates, self.ncat))
        for state in range(self.nstates):
            oscstr = self.oscstr[state]
            excitations = self.orbital_contributions[state]
            for pair in excitations:
                from_mo, to_mo, contribution = pair
                for icat in range(self.ncat):
                    cat = self.ex_categories_keys[icat]
                    from_range = occupied_groups[cat[0]]
                    to_range = virtual_groups[cat[1]]
                    if from_mo in from_range and to_mo in to_range:
                        scaled_os = oscstr * contribution
                        self.decomp_os[state, icat] = self.decomp_os[state, icat] + scaled_os

