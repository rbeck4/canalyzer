import numpy as np
import pandas as pd
from CANalyzer.load import Load
import CANalyzer.utilities as util
import sys

class MO(Load):
    def __init__(self, logfile, fchkfile, filename, groups, displaywidth, separate_ml=False, grouptotal=False, renormalize_negatives=False):
        super().__init__(logfile, fchkfile, filename, groups, displaywidth)
        self.moalpha = None
        self.mobeta = None
        self.overlap = None
        self.alpha_pop = None
        self.beta_pop = None
        self.sorted_alphapop = None
        self.sorted_betapop = None
        self.reduced_groups = None
        self.alpha_orbital_energy = None
        self.beta_orbital_energy = None
        self.separate_ml = separate_ml
        self.grouptotal = grouptotal
        self.renormalize_negatives = renormalize_negatives

    def mulliken_analysis(self):
        self.overlap = self.read_overlap()
        self.moalpha, self.mobeta = self.read_mo()
        cs = np.matmul(self.overlap, self.moalpha)
        if self.nri == 1:
            pop = np.multiply(self.moalpha, cs)
        else:
            conj_moalpha = np.conj(np.copy(self.moalpha))
            pop = np.real(np.multiply(conj_moalpha, cs))

        if self.xhf in ['GHF', 'GCAS', 'UHF', 'DHF']:
            self.alpha_pop = pop
            csbeta = np.matmul(self.overlap, self.mobeta)
            if self.nri == 1:
                self.beta_pop = np.multiply(self.mobeta, csbeta)
            elif self.nri == 2:
                conj_mobeta = np.conj(np.copy(self.mobeta))
                self.beta_pop = np.real(np.multiply(conj_mobeta, csbeta))
        else:
            self.alpha_pop = pop

        if self.xhf == 'DHF':
            moalpha_sp, mobeta_sp = self.read_mo_sp()

            kinetic_energy = self.read_kinetic()
            csalpha = kinetic_energy @ moalpha_sp
            csbeta = kinetic_energy @ mobeta_sp

            conj_moalpha_sp = np.conj(np.copy(moalpha_sp))
            conj_mobeta_sp = np.conj(np.copy(mobeta_sp))

            alpha_pop_sp = (1/(2*util.c*util.c)) * np.multiply(conj_moalpha_sp, csalpha)
            beta_pop_sp = (1/(2*util.c*util.c)) * np.multiply(conj_mobeta_sp, csbeta)

            self.alpha_pop += np.real(alpha_pop_sp)
            self.beta_pop += np.real(beta_pop_sp)

        if self.separate_ml:
            get_ml = 4
        else:
            get_ml = 3

        if not self.groups:
            self.reduced_groups = []
            [self.reduced_groups.append(x[1:get_ml]) for x in self.subshell if x[1:get_ml] not in self.reduced_groups]
            self.sorted_alphapop = np.zeros((len(self.reduced_groups), self.nbsuse*self.ncomp),dtype="complex128")
            self.sorted_betapop = np.zeros((len(self.reduced_groups), self.nbsuse*self.ncomp),dtype="complex128")
            for j in range(self.nbsuse*self.ncomp):
                for i in range(self.nbasis):
                    atomshellpair = self.subshell[i][1:get_ml]
                    index = self.reduced_groups.index(atomshellpair)
                    self.sorted_alphapop[index, j] += self.alpha_pop[i, j]
                    if self.xhf in ['UHF', 'GHF', 'GCAS', 'DHF']:
                        self.sorted_betapop[index, j] += self.beta_pop[i, j]
        else:
            self.separate_ml = False
            self.reduced_groups = [(g, l) for g in self.groupnames for l in range(self.maxL + 1)]
            self.sorted_alphapop = np.zeros((len(self.reduced_groups), self.nbsuse*self.ncomp))
            self.sorted_betapop = np.zeros((len(self.reduced_groups), self.nbsuse*self.ncomp))
            for j in range(self.nbsuse*self.ncomp):
                for i in range(self.nbasis):
                    atomshellpair = self.subshell[i]
                    atom_num = atomshellpair[0]
                    l = util.OAM[atomshellpair[2]]
                    group = self.identify_group(atom_num)
                    index = self.reduced_groups.index((group, l))
                    self.sorted_alphapop[index, j] += self.alpha_pop[i, j]
                    if self.xhf in ['UHF', 'GHF', 'GCAS', 'DHF']:
                        self.sorted_betapop[index, j] += self.beta_pop[i, j]
        self.sorted_alphapop = np.real(self.sorted_alphapop)
        self.sorted_betapop = np.real(self.sorted_betapop)

        if self.renormalize_negatives or self.grouptotal:
            self.renorm_neg()

        if self.grouptotal:
            # must account for case with no beta matrix
            try:
                num_groups = len(self.groupnames)
            except:
                num_groups = len(set(self.atoms))
            sum_over_groups_alphapop = np.zeros((num_groups, self.nbsuse))
            if self.xhf in ["UHF", "GHF", "GCAS"]:
                sum_over_groups_betapop = np.zeros((num_groups, self.nbsuse))
            auxlist = []
            for redgroup in self.reduced_groups:
                auxlist.append(redgroup[0])
            val, index, count = np.unique(auxlist, return_counts=True, return_index=True)
            val = val[np.argsort(index)]
            count = count[np.argsort(index)]
            self.reduced_groups = [("", i) for i in val]

            global_count = 0
            for i in range(num_groups):
                num_subgroup = count[i]
                for l in range(num_subgroup):
                    sum_over_groups_alphapop[i,:] +=  self.sorted_alphapop[global_count,:]
                    if self.xhf in ["UHF", "GHF", "GCAS"]:
                        sum_over_groups_betapop[i,:] +=  self.sorted_betapop[global_count,:]
                    global_count += 1
            self.sorted_alphapop = sum_over_groups_alphapop
            if self.xhf in ["UHF", "GHF", "GCAS"]:
                self.sorted_betapop = sum_over_groups_betapop


    def renorm_neg(self):
        # treats all negative populations as positive and renormalize
        try:
            num_groups = len(self.groupnames)
        except:
            num_groups = len(set(self.atoms))
        num_subgroups = len(self.reduced_groups)
        unnorm_pos_alpha = np.abs(self.sorted_alphapop)
        if self.xhf in ["UHF", "GHF", "GCAS"]:
            unnorm_pos_beta = np.abs(self.sorted_betapop)
        for i in range(self.nbsuse):
            unnorm_pos_alpha[:,i] = unnorm_pos_alpha[:,i]*(1/np.sum(unnorm_pos_alpha[:,i]))
            if self.xhf in ["UHF", "GHF", "GCAS"]:
                unnorm_pos_beta[:,i] = unnorm_pos_beta[:,i]*(1/np.sum(unnorm_pos_beta[:,i]))
        self.sorted_alphapop = unnorm_pos_alpha
        if self.xhf in ["UHF", "GHF", "GCAS"]:
            self.sorted_betapop = unnorm_pos_beta


    def identify_group(self, atom_number):
        group='undef'
        for name in self.groupnames:
            grouprange = self.groups[name]
            for r in grouprange:
                start = r[0]
                end = r[1] + 1
                if atom_number in range(start, end):
                    group = name
                    break
        return group


    def print_mulliken(self):
        self.alpha_orbital_energy, self.beta_orbital_energy = self.read_orbitalenergy()
        if self.separate_ml and not self.grouptotal:
            header = ('Orbital Energy', '(Hartree)', 'ml')
        else:
            header = ('Orbital Energy', '(Hartree)')
        if self.xhf in ['GHF', 'GCAS', 'DHF']:
            results_alpha = dict(zip([header] + self.reduced_groups,
                                     np.append(self.alpha_orbital_energy.round(5).T,
                                               (self.sorted_alphapop + self.sorted_betapop).round(3), axis=0)))
            remark_alpha = f"\n{self.xhf} Orbitals\n"
            self.alpha_pop = np.real(self.alpha_pop.round(3))
            self.beta_pop = np.real(self.beta_pop.round(3))
            alpha_spin = [np.sum(self.alpha_pop[:, i]) for i in range(self.nbsuse*self.ncomp)]
            beta_spin = [np.sum(self.beta_pop[:, i]) for i in range(self.nbsuse*self.ncomp)]
            self.spin = dict(zip(['Orbital Energy (Hartree)', 'Alpha', 'Beta'],
                                 [self.alpha_orbital_energy.flatten().round(5), alpha_spin, beta_spin]))
        elif self.xhf in ['RHF', 'ROHF', 'RCAS']:
            results_alpha = dict(zip([header] + self.reduced_groups,
                                     np.append(self.alpha_orbital_energy.round(5).T,self.sorted_alphapop.round(3), axis=0)))
            remark_alpha = f"\n{self.xhf} Orbitals\n"
        else:
            results_alpha = dict(zip([header] + self.reduced_groups,
                                     np.append(self.alpha_orbital_energy.round(5).T,self.sorted_alphapop.round(3), axis=0)))
            remark_alpha = f"\n{self.xhf} Alpha Orbitals\n"
            results_beta = dict(zip([header] + self.reduced_groups,
                                    np.append(self.beta_orbital_energy.round(5).T,self.sorted_betapop.round(3), axis=0)))
            remark_beta = f"\n{self.xhf} Beta Orbitals\n"

        with open(self.filename, 'w') as sys.stdout, pd.option_context('display.max_rows', None, 'display.max_columns', None):
            if self.software == 'CQ':
                print(f'Number of Total Electrons: {self.nae}')
            else:
                print(f'Number of Alpha Electrons: {self.nae}   Beta Electrons: {self.nbe}')
            print(f'Number of AOs: {self.nbasis*self.ncomp}    MOs: {self.nbsuse*self.ncomp}')
            pd.set_option('display.width', self.displaywidth)
            results_alpha_df = pd.DataFrame(results_alpha)
            results_alpha_df.index += 1
            print(remark_alpha)
            print(results_alpha_df)
            if self.xhf in ['UHF']:
                results_beta_df = pd.DataFrame(results_beta)
                results_beta_df.index += 1
                print(remark_beta)
                print(results_beta_df)
            elif self.xhf in ['GHF', 'GCAS', 'DHF']:
                spin_df = pd.DataFrame(self.spin)
                spin_df.index += 1
                print('\nSpin Contribution\n')
                print(spin_df)




















