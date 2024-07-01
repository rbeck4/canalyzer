import numpy as np
import subprocess
import linecache
import CANalyzer.utilities as util
import ast

class Load:
    def __init__(self, logfile, fchkfile, filename, groups, displaywidth):
        self.logfile = logfile
        self.fchkfile = fchkfile
        self.filename = filename
        self.nbasis = None
        self.nbsuse = None
        self.maxL = None
        self.atoms = []
        self.groups = groups
        self.groupnames = None
        self.natoms = None
        self.nae = None
        self.nbe = None
        self.subshell = None
        self.xhf = 'RHF'
        self.ri = 'Real'
        self.pureD = 0
        self.pureF = 0
        self.ncomp = 1
        self.nri = 1
        self.thresh = 1e-10
        self.displaywidth = displaywidth


    def parse_constants(self):
        # getting atoms
        self.natoms = int(
            str(subprocess.check_output("grep 'Atomic numbers ' " + self.fchkfile, shell=True)).split(" ")[
                -1].split("\\")[0])
        atomsline = int(str(subprocess.check_output("grep -n 'Atomic numbers' " + self.fchkfile, shell=True)).split(" ")
                        [0].split("'")[1].split(":")[0]) + 1
        atomsrows = int(np.ceil(self.natoms / 6))
        for i in range(atomsrows):
            atoms = linecache.getline(self.fchkfile, atomsline + i).split(" ")
            self.atoms += [util.ptable[x] for x in util.onlytype(atoms,int)]

        # converting groups to dict
        if self.groups:
            rawgroups = ast.literal_eval(self.groups)
            self.groupnames = rawgroups.keys()
            grouprange = []
            for key in self.groupnames:
                rawrange = rawgroups[key]
                rawrange = rawrange.replace('[', '(')
                rawrange = rawrange.replace(']', ')')
                rawrange = rawrange.replace(',', '),(')
                rawrange = rawrange.replace(':', ',')
                rawrange = ast.literal_eval(f"[{rawrange}]")
                grouprange.append(rawrange)
            self.groups = dict(zip(self.groupnames, grouprange))

        # getting parameters
        self.nae = int(
            str(subprocess.check_output("grep 'Number of alpha electrons ' " + self.fchkfile, shell=True)).split(" ")[
                -1].split("\\")[0])
        self.nbe = int(
            str(subprocess.check_output("grep 'Number of beta electrons ' " + self.fchkfile, shell=True)).split(" ")[
                -1].split("\\")[0])
        self.nbasis = int(
            str(subprocess.check_output("grep 'Number of basis functions ' " + self.fchkfile, shell=True)).split(" ")[
                -1].split("\\")[0])
        self.nbsuse = int(str(subprocess.check_output("grep 'Number of independent functions ' " + self.fchkfile,
                shell=True)).split(" ")[-1].split("\\")[0])
        self.pureD = int(
            str(subprocess.check_output("grep 'Pure/Cartesian d shells' " + self.fchkfile, shell=True)).split(
                " ")[-1].split("\\")[0])
        self.pureF = int(
            str(subprocess.check_output("grep 'Pure/Cartesian f shells' " + self.fchkfile, shell=True)).split(
                " ")[-1].split("\\")[0])
        self.maxL = int(
            str(subprocess.check_output("grep 'Highest angular momentum' " + self.fchkfile, shell=True)).split(
                " ")[-1].split("\\")[0])


    def parse_log(self):
        # what kind of HF/KS
        overlay3 = self.overlay_route(3)
        set_iops_3 = overlay3.keys()
        if 116 in set_iops_3:
            iop116 = overlay3[116]
            if iop116 == 1:
                self.xhf = "RHF"
                self.ri = 'Real'
            elif iop116 == 2:
                self.xhf = "UHF"
                self.ri = 'Real'
            elif iop116 == 3:
                self.xhf = "RHF"
                self.ri = "Complex"
            elif iop116 == 4:
                self.xhf = "UHF"
                self.ri = "Complex"
            elif iop116 == 7:
                self.xhf = "GHF"
                self.ri = "Complex"
            elif iop116 == 101:
                while self.xhf not in ["ROHF", "RCAS"]:
                    self.xhf = input("Is this ROHF or RCAS?\n")
                self.ri = "Real"
            elif iop116 == 103:
                while self.xhf not in ["ROHF", "RCAS"]:
                    self.xhf = input("Is this ROHF or RCAS?\n")
                self.ri = "Complex"
            elif iop116 == 107:
                self.xhf = "GCAS"
                self.ri = "Complex"
            else:
                while self.xhf not in ["ROHF", "RHF", "UHF", "GHF", "GCAS", "RCAS"]:
                    self.xhf = input("Is this RHF, UHF, ROHF, GHF, RCAS, GCAS?\n")
                while self.ri not in ['Real', 'Complex']:
                    self.ri = input("Is this Real or Complex?\n")
                if self.xhf == 'GHF' and self.ri == 'Real':
                    raise Exception("Real GHF NYI")
        else:
            while self.xhf not in ["ROHF", "RHF", "UHF", "GHF", "GCAS", "RCAS"]:
                self.xhf = input("Is this RHF, UHF, ROHF, GHF, RCAS, GCAS?\n")
            while self.ri not in ['Real', 'Complex']:
                self.ri = input("Is this Real or Complex?\n")
                if self.xhf == 'GHF' and self.ri == 'Real':
                    raise Exception("Real GHF NYI")

        if self.xhf in ['GHF', 'GCAS']:
            self.ncomp = 2
        else:
            self.ncomp = 1

        if self.ri == 'Real':
            self.nri = 1
        elif self.ri == 'Complex':
            self.nri = 2

        # basis info
        basis = [[] for i in range(self.natoms)]
        shells = list(util.OAM.keys())
        shells.append('SP')
        basisstart = int(str(subprocess.check_output(
            "grep -n 'AO basis set in the form of general basis input (Overlap normalization):' " + self.logfile,
            shell=True)).split(" ")[0].split("'")[1].split(":")[0]) + 2
        atomcount = 0
        while atomcount < self.natoms:
            line = linecache.getline(self.logfile, basisstart).split(" ")
            if line[1] in shells:
                basis[atomcount].append(line[1])
            if "****\n" in line:
                atomcount += 1
            basisstart += 1

        new_basis = []
        for a in basis:
            try:
                spindex = a.index("SP",0, len(a))
                spcount = 0
                while 'SP' in a:
                    a.remove("SP")
                    spcount += 1
                for p in range(0, spcount):
                    a.insert(spindex + p, 'P')
                for p in range(0, spcount):
                    a.insert(spindex + p, 'S')
            except:
                pass

            new_a = []
            for ss in range(len(a)):
                l = int(util.OAM[a[ss]])
                multiplicity = 2 * l + 1
                if l == 2:
                    if self.pureD == 1:
                        multiplicity = 6
                    elif self.pureD == 0:
                        multiplicity = 5
                else:
                    if self.pureF == 1:
                        multiplicity = int((l + 1) * (l + 2) / 2)
                    elif self.pureF == 0:
                        multiplicity = int(2 * l + 1)

                for i in range(multiplicity):
                    new_a.append(a[ss])

            new_basis.append(new_a)

        self.subshell = []
        for i in range(self.natoms):
            for j in new_basis[i]:
                self.subshell.append((i, self.atoms[i], j))

    def read_overlap(self):
        return self.readlog_matrix(r"\*\*\* Overlap \*\*\*", self.nbasis, self.nbasis, True,False)


    def read_mo(self, separate=True):
        moalpha = None
        mobeta = None
        if self.xhf in ['RHF', 'ROHF', 'RCAS']:
            moalpha = self.readfchk_matrix("Alpha MO coefficients", self.nbasis*self.ncomp,
                                                self.nbsuse*self.ncomp, False, False)
        elif self.xhf == 'UHF':
            moalpha = self.readfchk_matrix("Alpha MO coefficients", self.nbasis, self.nbsuse, False,
                                                False)
            mobeta = self.readfchk_matrix("Beta MO coefficients", self.nbasis, self.nbsuse, False,
                                               False)
        else:
            moalpha = self.readfchk_matrix("Alpha MO coefficients", self.nbasis * self.ncomp,
                                           self.nbsuse * self.ncomp, False, False)
            if separate:
                mobeta = moalpha[1::2, :]
                moalpha = moalpha[::2, :]

        return moalpha, mobeta


    def read_orbitalenergy(self):
        beta_energy = None
        alpha_energy = self.readfchk_matrix("Alpha Orbital Energies",  self.nbsuse*self.ncomp, 1,
                                            False, False, True)
        if self.xhf == 'UHF':
            beta_energy = self.readfchk_matrix("Beta Orbital Energies",  self.nbsuse * self.ncomp, 1,
                                                False, False, True)
        return alpha_energy, beta_energy


    def readlog_matrix(self, startstr, nrows, ncol, ifltt=False, ifantisymm=False, instance=1):
        # for files with multiple matrices with the same startstr, this pick which instance to read
        if instance == 1:
            startline = int(str(subprocess.check_output(f"grep -n '{startstr}' {self.logfile}", shell=True)).split(" ")
                        [0].split("'")[1].split(":")[0]) + 2
        elif instance > 1:
            grep_results = str(
                subprocess.check_output(f"grep -n '{startstr}' {self.logfile}", shell=True)).split("\\n")
            linenumbers = []
            for line in grep_results:
                words = line.split(':')
                try:
                    linenumbers.append(int(words[0]))
                except:
                    pass
            startline = linenumbers[instance-2] + 2

        else:
            raise Exception('Instance must be greater than 0')

        # initialize matrix
        matrix = np.zeros((nrows, ncol))

        # reading matrix
        #if ifltt:
        blocks = int(np.ceil(ncol / 5))
        col_offset = 0
        for b in range(blocks):
            for line in range(nrows):
                linear2= []
                matrixline = linecache.getline(self.logfile, startline)
                if 'D' in matrixline:
                    matrixline = matrixline.replace("D", "e")
                if 'd' in matrixline:
                    matrixline = matrixline.replace('d', 'e')

                matrixlist = matrixline.split()
                for x in matrixlist:
                    try:
                        xx = float(x)
                        if abs(xx) < self.thresh:
                            xx = 0.0
                        linear2.append(xx)
                    except:
                        pass

                try:
                    dump = linear2.pop(0)
                    dump = int(dump - 1)
                except:
                    pass

                startline += 1

                for x in range(len(linear2)):
                    matrix[dump, x + col_offset] = linear2[x]
            if ifltt:
                nrows -= 5
            col_offset += 5
            startline += 1

        if ifltt:
            if ifantisymm:
                matrix = np.tril(matrix, -1).T - matrix
            elif not ifantisymm:
                matrix = np.tril(matrix, -1).T + matrix
        else:
            matrix = matrix.reshape((nrows, ncol))

        return matrix


    def readfchk_matrix(self, startstr, nrows, ncol, ifltt=False, ifantisymm=False, realonly=False):
        startline = int(str(subprocess.check_output(f"grep -n '{startstr}' {self.fchkfile}", shell=True)).split(" ")[
            0].split("'")[1].split(":")[0]) + 1
        rawmatrix = []
        matrixsize = nrows * ncol
        if self.nri == 2 and not realonly:
            matrixsize = matrixsize * self.nri
        num_lines = int(np.ceil(matrixsize / 5))
        for i in range(num_lines):
            matrixline = linecache.getline(self.fchkfile, startline).split(" ")
            for x in matrixline:
                try:
                    xx = float(x)
                    if abs(xx) < self.thresh:
                        xx = 0.0
                    rawmatrix.append(xx)
                except:
                    pass
            startline += 1

        rawmatrix = np.array(rawmatrix)
        if self.nri == 2 and not realonly:
            real_matrix = rawmatrix[::2].reshape((ncol, nrows)).T
            imag_matrix = rawmatrix[1::2].reshape((ncol, nrows)).T
            matrix = real_matrix + 1j * imag_matrix
        else:
            matrix = rawmatrix.reshape((ncol, nrows)).T

        return matrix

    def overlay_route(self, overlay):
        line1 = str(subprocess.check_output("grep ' 3/' " + self.logfile, shell=True)).split(" ")[-1].split(",")
        line2 = []
        first = True
        for x in line1:
            if "=" in x:
                if "/" in x:
                    id = x.index("/")
                    if first:
                        x = x[id + 1:]
                        first = False
                    elif not first:
                        x = x[:id]
                x = x.replace("=", ":")
                line2.append(x)
        line2[0] = "{" + line2[0]
        line2[-1] = line2[-1] + "}"
        final = ast.literal_eval(",".join(line2))
        return final

    def start(self):
        self.parse_constants()
        self.parse_log()





