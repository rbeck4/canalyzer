import numpy as np
import subprocess
import linecache
from itertools import chain
import matplotlib.pyplot as plt

from CANalyzer.spectra import Spectra
from CANalyzer.mo import MO

class TDDFT_spectra(Spectra, MO):
    def __init__(self, logfile, fchkfile, groups, orbital_decomposition=True, separate_ml=False):
        MO.__init__(self, logfile, fchkfile, None, groups, None, separate_ml)
        Spectra.__init__(self)
        self.categories = None
        # categories of excitations between partitions of orbitals

        self.occupied_spaces = None
        # partitioning of occupied orbitals
        # orbital partitioning as a dictionary of a list of tuples with start and end indices (inclusive, AS index)
        # {"Name" : (start, end)}

        self.virtual_spaces = None
        # partitioning of virtual orbitals
        # orbital partitioning as a dictionary of a list of tuples with start and end indices (inclusive, AS index)
        # {"Name" : (start, end)}

        self.orbital_decomposition = orbital_decomposition
        # True to decompose spectra into excitation categories

        self.decomp_byexcat = None
        # nroots x ncat matrix where element (i,j) is the oscillator strength of root i for excitation category j

        self.nstates = None
        # number of excited states from TDDFT calculation

        self.nroots = None
        # number of non-zero oscillator strengths

        self.orbital_contributions = None
        # (i, j, k) where i is the leaving orbital, j is the arriving orbital, and k is the oscillator strength contributed by i -> j

        self.ncat = None
        # number of excitation categories

        self.oscstr = []
        # original osc str

        self.energy = []
        # excitation energies

        super().start()
        self.parse_tddft()


    def parse_tddft(self):
        if self.software == "GDV":
            startline = int(
                str(subprocess.check_output("grep -n 'Excited State' " + self.logfile, shell=True)).split()[
                    0].split("'")[-1][:-1])
            self.nstates = int(
                str(subprocess.check_output("grep -n 'Excited State' " + self.logfile, shell=True)).split()[-8][:-1])
            self.nroots = self.nstates
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
            self.energy = np.array(self.energy)
            self.oscstr = np.array(self.oscstr)

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


    def purge(self, threshold=1E-7):
        dark_states = []
        self.energy = list(self.energy)
        self.oscstr = list(self.oscstr)
        for i in range(self.nstates):
            if self.oscstr[i] < threshold:
                dark_states.append(i)
        for i in sorted(dark_states, reverse=True):
            del self.oscstr[i]
            del self.energy[i]
            del self.orbital_contributions[i]
        ndark_state = len(dark_states)
        self.nroots = self.nstates - ndark_state
        self.energy = np.array(self.energy)
        self.oscstr = np.array(self.oscstr)


    def decompose_byexcat(self):
        occupied_spaces_keys = self.occupied_spaces.keys()
        virtual_spaces_keys = self.virtual_spaces.keys()

        occupied_ranges = [range(self.occupied_spaces[x][0], self.occupied_spaces[x][1]+1) for x in
                           occupied_spaces_keys]
        virtual_ranges = [range(self.virtual_spaces[x][0], self.virtual_spaces[x][1]+1) for x in
                           virtual_spaces_keys]

        occupied_superrange = chain.from_iterable(occupied_ranges)
        virtual_superrange = chain.from_iterable(virtual_ranges)

        min_specified = min(list(occupied_superrange))
        max_specified = max(list(virtual_superrange))

        self.occupied_spaces["Other Occ"] = (0, min_specified)
        self.virtual_spaces["Other Virt"] = (max_specified + 1, self.nbsuse)

        occupied_ranges.append(range(0, min_specified))
        virtual_ranges.append(range(max_specified + 1, self.nbsuse))

        occupied_spaces_keys = self.occupied_spaces.keys()
        virtual_spaces_keys = self.virtual_spaces.keys()

        occupied_ranges_dict = dict(zip(occupied_spaces_keys, occupied_ranges))
        virtual_ranges_dict = dict(zip(virtual_spaces_keys, virtual_ranges))

        self.categories = [(x, y) for x in occupied_spaces_keys for y in virtual_spaces_keys]
        self.ncat = len(self.categories)
        self.decomp_byexcat = np.zeros((self.nroots, self.ncat))

        for root in range(self.nroots):
            oscstr = self.oscstr[root]
            excitations = self.orbital_contributions[root]
            for pair in excitations:
                from_mo, to_mo, contribution = pair
                for icat in range(self.ncat):
                    cat = self.categories[icat]
                    from_range = occupied_ranges_dict[cat[0]]
                    to_range = virtual_ranges_dict[cat[1]]
                    if from_mo in from_range and to_mo in to_range:
                        scaled_os = oscstr * contribution
                        self.decomp_byexcat[root, icat] = self.decomp_byexcat[root, icat] + scaled_os


    def plot(self, spectra_names, xstart, xend, plotname=None, vline=None, ifsticks=True, stick_color="black", sticks=None):
        fig, ax = super().plot(spectra_names, xstart, xend, None, False, vline)
        if ifsticks:
            if sticks is None:
                sticks = self.energy

            for n in range(self.nroots):
                height = self.oscstr[n]*self.sticks_yscale
                root = sticks[n] - self.redshift
                ax.plot((root, root), (0, height), c=stick_color, linewidth=self.linewidth, alpha=self.alpha)

        if plotname:
            plt.savefig(plotname, dpi=200, format="png", bbox_inches='tight')

        plt.show()




