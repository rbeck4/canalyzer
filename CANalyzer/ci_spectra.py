import numpy as np
import subprocess
import linecache

from CANalyzer.spectra import Spectra
from CANalyzer.mo import MO
import matplotlib.pyplot as plt
import matplotlib as mpl

class CI_spectra(Spectra, MO):
    def __init__(self, logfile, fchkfile, groups=None, separate_ml=False):
        MO.__init__(self, logfile, fchkfile, None, groups, None, separate_ml)
        Spectra.__init__(self)
        self.occnum = None
        # approximate occupation numbers using 1PDM diagoanls
        # nstates x nactive matrix where element (i, j) is the occupation number for orbital j in state i
        # index offset by 1 due to difference between Python and state indexing

        self.nactive = None
        # number of active orbitals

        self.nactelec = None
        # number of electrons in active space

        self.nocc = None
        # number of occupied active space orbitals

        self.unocc = None
        # number of unoccupied active space orbitals

        self.decomp_byspaces = None
        # numfromstates x nstates x nspaces tensor where element (i, j, k)
        # is the state i -> j oscillator strength contribution from the change in electron population in space k
        # state index offset by 1 due to difference between Python and state indexing
        # elements between two leaving states are not populated and remains zero

        self.decomp_byorbital = None
        # numfromstates x nstates x ngroups tensor where element (i, j, k)
        # is the state i -> j oscillator strength contribution from the change 
        #     in electron population in group k
        # state index offset by 1 due to difference between Python and state indexing
        # elements between two leaving states are not populated and remains zero

        self.nstates = None
        # number of total states

        self.numfromstates = None
        # number of leaving states

        self.ndet = None
        # Total number of determinants

        self.nspaces = None
        # number of orbital partitions

        self.spaces = None
        # orbital partitioning as a dictionary of a list of tuples with start and end indices (inclusive, AS index)
        # {"Name" : (start, end)}

        self.groups = groups
        #MO decomp groups

        self.nroots = None
        # number of non-zero excitation energies to be plotted

        self.eigenVals = None
        # numstates is the energy of each state

        self.energy = None
        # numfromstates x nstates matrix where element (i, j) is the excitation energy between states i and j
        # energies between from states are not populated and remain zero
        # indexing is offset by 1 because output file state indexing starts at 1 and Python starts at 0

        self.oscstr = None
        # numfromstates x nstates matrix where element (i, j) is the oscillator strength between states i and j
        # oscillator strengths between from states are not populated and remain zero
        # indexing is offset by 1 because output file state indexing starts at 1 and Python starts at 0

        self.num_ex_electrons = None
        # number of electrons excited compared between from and to state

        self.parse_ci()


    def parse_ci(self):
        super().start()
        if self.software == "GDV":
            """
            GDV parse not done yet
            """
            energyline = str(subprocess.check_output("grep -n '(Hartree)' " + self.logfile, shell=True))
            energyline_split = energyline.split("\\n")
            self.energy = [float(x.split("(Hartree):")[-1]) for x in energyline_split[:-1]]
            self.nstates = int(energyline_split[-2].split(":")[2].split()[0])
            self.nactive = int(subprocess.check_output("grep 'Number of CAS Orbitals' " + self.fchkfile, shell=True).split()[-1])
            pdmdiag_list = []  # just 1PDM diagonals in MO basis
            stateEnergies = []
            
            enlines = str(subprocess.check_output("grep 'Energy (Hartree):' " + self.logfile, shell=True)).split("\\n")[:-1]
            for i in range(self.nstates):
                print('\r', "READING PDM FOR STATE: %i / %i" %(i, self.nstates), end='')
                stateEnergies.append(float(enlines[i].split()[-1])) 
                pdmdiag = self.readlog_matrix("For Simplicity", self.nactive, 1, instance=i+1).flatten()
                pdmdiag_list.append(pdmdiag)
            enlines = None
            self.occnum = np.array(pdmdiag_list)

            self.oscstr = []
            oslines = str(subprocess.check_output("grep 'Oscillator Strength For States' " + self.logfile, shell=True)).split("\\n")[:-1]
            oslines[0] = oslines[0].split("'")[-1]
            self.numfromstates = int(oslines[-1].split()[-5])
            
            self.oscstr = np.zeros([self.numfromstates,self.nstates])
            self.energy = np.zeros([self.numfromstates,self.nstates])
           
            for i in range(len(oslines)):
                oscstr = float(oslines[i].split("f=")[-1])
                sts    = oslines[i].split(":")
                tmpFrom = int(sts[0].split()[-1])-1
                tmpArrv = int(sts[1].split()[0])-1
                self.oscstr[tmpFrom,tmpArrv] = oscstr
                self.energy[tmpFrom,tmpArrv] = stateEnergies[tmpArrv] - stateEnergies[tmpFrom]

        elif self.software == "CQ":
            self.nactive = int(str(subprocess.check_output("grep 'Number of Correlated Orbitals:' " + self.logfile, shell=True)).split()[-1].split("\\n")[0])
            self.nactelec = int(str(subprocess.check_output("grep 'Number of Correlated Electrons:' " + self.logfile, shell=True)).split()[-1].split("\\n")[0])
            # unocc and nocc differs for 1-component. currently, cq does not have 1-component
            self.unocc = self.nactive - self.nactelec
            self.nocc = self.nactive - self.unocc
            grepline = str(subprocess.check_output("grep -n 'Excited State:' " + self.logfile, shell=True))

            #In case job doesn't completely finish all property evals
            #NOTE this relies on the job summary printing at the beginning of the job:
            self.ndet = int(str(subprocess.check_output("grep -n 'Total Number of Determinants' " + self.logfile, shell=True)).split('=')[-1].split('\\n')[0])
            self.nstates = int(str(subprocess.check_output("grep -n 'NROOTS:' " + self.logfile, shell=True)).split(':')[-1].split('\\n')[0])
            #In case number of requested roots > number det:
            self.nstates = min(self.ndet, self.nstates)
            self.numfromstates = int(str(subprocess.check_output("grep -n 'OSCISTREN:' " + self.logfile, shell=True)).split(':')[-1].split('\\n')[0])

            self.oscstr = np.zeros((self.numfromstates, self.nstates))
            self.energy = np.zeros((self.numfromstates, self.nstates))
            self.occnum = np.zeros((self.nstates, self.nactive))
            self.nroots = self.numfromstates * self.nstates

            # start parsing for excitation energies and oscillator strength
            startline_os = int(grepline.split("'")[1].split(":")[0])
            statecounter = self.numfromstates
            fromstatecounter = 1
            for x in range(self.nroots*2):
                #X2 as CQ double spaces them...
                parseline = linecache.getline(self.logfile, startline_os)
                try:
                    #
                    #CQ printing:
                    #Ex State: endSt to state: fromSt En F
                    #
                    #Energy and Osc. storage in 2D array (len(fromStates) by 
                    #len(nStates)).
                    endSt = int(parseline.split()[2])-1
                    fromSt = int(parseline.split()[5].split(':')[0])-1
                    self.energy[fromSt, endSt] = float(parseline.split("=")[1].split("f")[0])
                    self.oscstr[fromSt, endSt] = float(parseline.split("=")[-1].split("f")[0])
                    statecounter += 1
                except:
                    pass
                startline_os += 1
                if statecounter >= self.nstates:
                    statecounter = self.numfromstates
                    fromstatecounter += 1

            # start parsing for PDM diagonals
            grepline = str(subprocess.check_output("grep -n 'diagonal elements of real 1-RDM' " + self.logfile, shell=True))
            startline_pdm = int(grepline.split("'")[1].split(":")[0]) + 2
            statecounter = 0
            lastcycle = False

            while statecounter <= self.nstates:
                parseline = linecache.getline(self.logfile, startline_pdm).split()
                try:
                    if "State" in parseline[0]:
                        statecounter += 1
                        occupation_numbers = [(int(x.split("(")[0]), float(x.split("(")[1][:-1])) for x in parseline[2:]]
                        for n in occupation_numbers:
                            self.occnum[statecounter - 1, n[0] - 1] = n[1]
                        if statecounter == self.nstates:
                            lastcycle = True
                    else:
                        #Final round, presumes '--' is used as closer...
                        if lastcycle and "-------------------" in parseline[0]:
                            statecounter += 1
                        occupation_numbers = [(int(x.split("(")[0]), float(x.split("(")[1][:-1])) for x in parseline]
                        for n in occupation_numbers:
                            self.occnum[statecounter - 1, n[0] - 1] = n[1]
                except:
                    pass
                startline_pdm += 1
            
            # start parsing for EigenValues
            self.eigenVals = np.asarray([0. for x in range(self.nstates)])
            grepline = str(subprocess.check_output("grep -m%i 'State:' " %self.nstates + self.logfile, shell=True)).split('\\n')
            for i in range(self.nstates):
                self.eigenVals[i] = float(grepline[i].split(':')[-1])


    def decompose_byspaces(self):
        if not self.spaces:
            raise Exception("self.spaces attributes must be defined")

        spaces_names = list(self.spaces.keys())
        self.nspaces = len(spaces_names)

        self.decomp_byspaces = np.zeros((self.numfromstates, self.nstates, self.nspaces))

        electron_pop_change = np.zeros((self.numfromstates, self.nstates, self.nactive))
        for i in range(self.numfromstates):
            for j in range(self.nstates):
                electron_pop_change[i, j, :] = self.occnum[j, :] - self.occnum[i, :]

        self.num_ex_electrons = np.zeros((self.numfromstates, self.nstates))
        for i in range(self.numfromstates):
            for j in range(self.nstates):
                self.num_ex_electrons[i, j] = np.sum(np.abs(electron_pop_change[i, j, :])) / 2

        for n in range(self.numfromstates):
            for nspace in range(self.nspaces):
                space_name = spaces_names[nspace]
                space_start = self.spaces[space_name][0] - 1
                space_end = self.spaces[space_name][1]
                for state in range(self.numfromstates, self.nstates):
                    if self.num_ex_electrons[n, state] < 0.6:
                        print(f"States {n+1} and {state+1} has nearly the same electron occupation. ")
                        if self.num_ex_electrons[n, state] < 0.01:
                            print(f"Omitting contributions from this excitation. Beware of results.")
                            continue
                    pop_change = np.sum(electron_pop_change[n, state, space_start:space_end])
                    scaling = pop_change / self.num_ex_electrons[n, state]
                    self.decomp_byspaces[n, state, nspace] = scaling * self.oscstr[n, state]

    def decompose_byorbital(self):
        MO.start()
        MO.mulliken_analysis()
        if not self.groups:
            self.groups = MO.reduced_groups

        groups_names = list(self.groups.keys())
        self.nspaces = len(groups_names)

        self.decomp_byorbital = np.zeros((self.numfromstates, self.nstates, self.nspaces))

        electron_pop_change = np.zeros((self.numfromstates, self.nstates, self.nactive))
        for i in range(self.numfromstates):
            for j in range(self.nstates):
                electron_pop_change[i, j, :] = self.occnum[j, :] - self.occnum[i, :]

        self.num_ex_electrons = np.zeros((self.numfromstates, self.nstates))
        for i in range(self.numfromstates):
            for j in range(self.nstates):
                self.num_ex_electrons[i, j] = np.sum(np.abs(electron_pop_change[i, j, :])) / 2

        for n in range(self.numfromstates):
            for nspace in range(self.nspaces):
                space_name = groups_names[nspace]
                space_start = self.groups[space_name][0] - 1
                space_end = self.groups[space_name][1]
                for state in range(self.numfromstates, self.nstates):
                    if self.num_ex_electrons[n, state] < 0.6:
                        print(f"States {n+1} and {state+1} has nearly the same electron occupation. ")
                        if self.num_ex_electrons[n, state] < 0.01:
                            print(f"Omitting contributions from this excitation. Beware of results.")
                            continue
                    pop_change = np.sum(electron_pop_change[n, state, space_start:space_end])
                    scaling = pop_change / self.num_ex_electrons[n, state]
                    self.decomp_byspaces[n, state, nspace] = scaling * self.oscstr[n, state]
