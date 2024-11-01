import numpy as np
import subprocess
import linecache
import CANalyzer.utilities as util
import ast
import h5py

class Load:
    def __init__(self, logfile, fchkfile, filename, groups, displaywidth):
        self.software = None
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

        if self.fchkfile[-4:] == ".bin" and self.logfile[-4:] == ".out":
            self.software = 'CQ'
        else:
            self.software = 'GDV'

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


    def parse_constants_gdv(self):
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
            line = linecache.getline(self.logfile, basisstart).split()
            if line[0] in shells:
                basis[atomcount].append(line[0])
            if "****" in line:
                atomcount += 1
            basisstart += 1

        new_basis = []
        for a in range(len(basis)):
            for ss in basis[a]:
                if len(ss) == 1:
                    l = int(util.OAM[ss])
                    multiplicity = 2 * l + 1
                    ml = []
                    for i in range(l+1):
                        if i == 0:
                            ml.append(i)
                        else:
                            ml.append(i)
                            ml.append(-i)
                    for i in range(multiplicity):
                        new_basis.append((a+1, self.atoms[a], ss, ml[i]))
                else:
                    ss_split = list(ss)
                    for sss in ss_split:
                        l = int(util.OAM[sss])
                        multiplicity = 2 * l + 1
                        if sss != "P":
                            ml = []
                            for i in range(l+1):
                                if i == 0:
                                    ml.append(i)
                                else:
                                    ml.append(i)
                                    ml.append(-i)
                            for i in range(multiplicity):
                                new_basis.append((a+1, self.atoms[a], sss, ml[i]))
                        elif sss == "P":
                            new_basis.append((a+1, self.atoms[a], sss, 1))
                            new_basis.append((a+1, self.atoms[a], sss, -1))
                            new_basis.append((a+1, self.atoms[a], sss, 0))

        self.subshell = new_basis


    def read_overlap(self):
        if self.software == 'CQ':
            overlap = self.readbin_matrix("/INTS/OVERLAP")
            pass
        else:
            overlap = self.readlog_matrix(r"\*\*\* Overlap \*\*\*", self.nbasis, self.nbasis, True,False)
        return overlap


    def read_mo(self, separate=True):
        moalpha = None
        mobeta = None
        if self.software == 'CQ':
            if self.xhf in ['RHF', 'ROHF', 'RCAS']:
                moalpha = self.readbin_matrix("/SCF/MO1").T
            elif self.xhf == "UHF":
                moalpha = self.readbin_matrix("/SCF/MO1").T
                mobeta = self.readbin_matrix("/SCF/MO2").T
            elif self.xhf == 'DHF':
                moalpha = self.read_mo_full4c()
                moalpha = moalpha[:self.nbasis * 2, :self.nbasis * 2]
                if separate:
                    mobeta = moalpha[self.nbasis:, :]
                    moalpha = moalpha[:self.nbasis, :]
            else:
                moalpha = self.readbin_matrix("/SCF/MO1").T
                if separate:
                    mobeta = moalpha[self.nbasis:, :]
                    moalpha = moalpha[:self.nbasis, :]
        else:
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


    def read_mo_sp(self):
        moalpha = self.read_mo_full4c()
        moalpha = moalpha[self.nbasis*2:, :self.nbasis*2]
        mobeta = moalpha[self.nbasis:, :]
        moalpha = moalpha[:self.nbasis, :]
        return moalpha, mobeta


    def read_mo_full4c(self):
        if self.xhf != 'DHF':
            raise Exception("read_mo_full4c only for DHF")

        C = np.zeros((4 * self.nbasis, 4 * self.nbasis), dtype="complex128")
        preC = np.zeros((4 * self.nbasis, 4 * self.nbasis), dtype="complex128")
        rawC = self.readbin_matrix("/SCF/MO1").T

        preC[:, :2 * self.nbasis] = rawC[:, 2 * self.nbasis:]
        preC[:, 2 * self.nbasis:] = rawC[:, :2 * self.nbasis]

        C[:self.nbasis, :] = preC[:self.nbasis, :]
        C[3 * self.nbasis:, :] = preC[3 * self.nbasis:, :]
        C[self.nbasis:2 * self.nbasis, :] = preC[2 * self.nbasis:3 * self.nbasis, :]
        C[2 * self.nbasis:3 * self.nbasis, :] = preC[self.nbasis:2 * self.nbasis, :]

        return C


    def read_kinetic(self):
        if self.software == 'CQ':
            return self.readbin_matrix("/INTS/KINETIC")


    def read_orbitalenergy(self):
        beta_energy = None
        if self.software == 'CQ':
            if self.xhf != 'UHF':
                startline_occ = int(str(subprocess.check_output(f"grep -n 'Orbital Eigenenergies ' {self.logfile}", shell=True)).split(":")[0].split("'")[1]) + 3
                alpha_energy = []
                while True:
                    try:
                        energies = [float(x) for x in linecache.getline(self.logfile, startline_occ).split()]
                        alpha_energy += energies
                        startline_occ += 1
                    except:
                        break
                startline_virt = startline_occ + 1
                while True:
                    try:
                        energies = [float(x) for x in linecache.getline(self.logfile, startline_virt).split()]
                        alpha_energy += energies
                        startline_virt += 1
                    except:
                        break
                alpha_energy = np.array([alpha_energy]).T

            else:
                raise Exception("UHF NYI for CQ")

        else:
            alpha_energy = self.readfchk_matrix("Alpha Orbital Energies",  self.nbsuse*self.ncomp, 1,
                                                False, False, True)
            if self.xhf == 'UHF':
                beta_energy = self.readfchk_matrix("Beta Orbital Energies",  self.nbsuse * self.ncomp, 1,
                                                    False, False, True)
        return alpha_energy, beta_energy


    def readlog_matrix(self, startstr, nrows, ncol, ifltt=False, ifantisymm=False, instance=1):
        # for files with multiple matrices with the same startstr, this pick which instance to read
        if instance == 1:
            startline = int(str(subprocess.check_output(f"grep -n '{startstr}' {self.logfile}", shell=True)).split()[0].split("'")[1][:-1]) + 2
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
                        linear2.append(x)

                dump = linear2.pop(0)
                try:
                    dump = int(dump - 1)
                except:
# catches the case where matrix rows are labeled by a/b
# may be source of failure in the future if a/b is capitalized of L/S component is added
                    fact = 2 - (ord(dump[-1]) - 96)
                    dump = int(dump[:-1]) * 2 - fact - 1

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

        if ifltt:
            raise Exception("Reading LTT from fchk NYI")

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


    def readbin_matrix(self, h5path, ndim=None, realonly=False):
        h5bin = h5py.File(self.fchkfile, 'r')
        matrix = h5bin[h5path][:]
        if ndim:
            matrix = np.reshape(matrix, (ndim, ndim))
        return matrix


    def overlay_route(self, overlay):
        line1 = str(subprocess.check_output(r"grep '\ " + f"{overlay}/' " + self.logfile, shell=True)).split()[-1].split(",")
        line2 = [x.replace("=",":") for x in line1 if "=" in x]
        first = "{" + line2[0].split("/")[1]
        last = line2[-1].split("/")[0] + "}"
        line2[0] = first
        line2[-1] = last
        final = ast.literal_eval(",".join(line2))
        return final


    def parse_constants_cq(self):
        with h5py.File(self.fchkfile, "r") as h5f:
            # loading ncomp
            self.ncomp = int(h5f['REF']['NCOMP'][0])

            # loading reference type
            reftype = int(h5f['REF']['REFTYPE'][0])
            if reftype == 1:
                self.xhf = "RHF"
            elif reftype == 2:
                self.xhf = "UHF"
            elif reftype == 3:
                self.xhf = "ROHF"
            elif reftype == 4:
                self.xhf = "GHF"
            elif reftype == 5:
                self.xhf = "DHF"
                # we only look at LL block in DHF MO analysis
                self.ncomp = 2
                self.nri = 2
            else:
                raise Exception("Invalid reference type detected")

        # loading NRI
        self.ri = str(subprocess.check_output(f"grep 'Reference:   ' {self.logfile}", shell=True)).split()[2]
        if self.ri == 'Real':
            self.nri = 1
        elif self.ri == 'Complex':
            self.nri = 2
        else:
            while self.ri not in ['Real', 'Complex']:
                self.ri = input("Is this calculation Real or Complex?   ")
                if self.ri == 'Real':
                    self.nri = 1
                elif self.ri == 'Complex':
                    self.nri = 2

        # loading pureD, pureF
        try:
            isForceCart = str(subprocess.check_output(f"grep 'NAtoms' {self.logfile}", shell=True)).split()
            if "true" in isForceCart or "True" in isForceCart:
                self.pureD = 1
                self.pureF = 1
        except:
            pass

        # loading atoms
        self.natoms = int(str(subprocess.check_output(f"grep 'NAtoms' {self.logfile}", shell=True)).split(" ")[
            -1][:-3])
        atomstart = int(str(subprocess.check_output(
            f"grep -n 'geom:' {self.logfile}",
            shell=True)).split(" ")[0].split("'")[1].split(":")[0])
        atomnumber = 1
        while atomnumber <= self.natoms:
            atom = linecache.getline(self.logfile, atomstart + atomnumber).split()[0]
            self.atoms.append(atom)
            atomnumber += 1

        # loading nbasis
        self.nbasis = int(str(subprocess.check_output(f"grep 'NBasis' {self.logfile}", shell=True)).split(" ")[
                              -1][:-3])
        self.nbsuse = self.nbasis

        # loading number of electrons
        self.nae = int(str(subprocess.check_output(f"grep 'Total Electrons' {self.logfile}", shell=True)).split(" ")[
                              -1][:-3])
        if self.xhf in ['GHF', 'DHF']:
            self.nbe = 0
        elif self.xhf == 'RHF':
            self.nbe = self.nae / 2
            self.nae = self.nbe
        else:
            totale = self.nae
            self.nae = input("How many alpha electrons?")
            self.nbe = input("How many beta electrons?")
            if self.nae + self.nbe != totale:
                raise Exception("Invalid number of alpha and beta electrons")


        # loading basis
        self.maxL = int(str(subprocess.check_output(f"grep '  Max L  ' {self.logfile}", shell=True)).split(" ")[
                           -1][:-3])

        basisstart = int(str(subprocess.check_output(f"grep -n 'EigV --' {self.logfile}", shell=True)).split(" ")[0].split("'")[-1][:-1]) + 2
        basiscount = 1
        linenumber = 0
        atomcount = 0
        basis = [[] for i in range(self.natoms)]
        while basiscount <= self.nbasis:
            line = linecache.getline(self.logfile, basisstart + linenumber)
            if line == "\n":
                atomcount += 1
            else:
                oam = line[18]
                str_ml = line[19] + line[20]
                try:
                    ml = int(str_ml)
                except:
                    if str_ml == "X":
                        ml = 1
                    elif str_ml == "Y":
                        ml = -1
                    elif str_ml == "Z":
                        ml = 0
                    else:
                        ml = 0
                basis[atomcount].append((oam, ml))
                basiscount += 1
            linenumber += 1

        self.subshell = []
        for i in range(self.natoms):
            for j in basis[i]:
                self.subshell.append((i+1, self.atoms[i], j[0], j[1]))


    def start(self):
        if self.software == 'CQ':
            self.parse_constants_cq()
        else:
            self.parse_constants_gdv()
            self.parse_log()





