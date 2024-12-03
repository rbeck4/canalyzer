import numpy as np
from numpy.linalg import svd, norm
import pandas as pd
from CANalyzer.load import Load
import CANalyzer.utilities as util
import subprocess

class Harmonics(Load):
    def __init__(self, logfile, fchkfile):
        super().__init__(logfile, fchkfile, None, None, None)
        #Load.__init__(self, logfile, fchkfile, None, None, None)
        self.normalmodes = None # rank-3 tensor with indices (mode, atom, xyz)
        self.geometry = None
        self.nmodes = None


    def start(self):
        super().start()

        self.geometry = self.readfchk_matrix(r"Current cartesian coordinates", 3, self.natoms).T / 1.8897259886
        self.nmodes = int(str(subprocess.check_output(f"grep 'Number of Normal Modes' {self.fchkfile}", shell=True)).split("\\")[-2].split()[-1])

        self.normalmodes = self.readfchk_matrix(r"Vib-Modes", self.natoms*3*self.nmodes, 1)
        self.normalmodes = np.reshape(self.normalmodes, (self.nmodes, self.natoms, 3))


    def distort(self, mode, amount):
        return self.geometry + amount * self.normalmodes[mode, :, :]


    """def pca_normalmodes(self, weight):
        if len(weight) != self.nmodes:
            raise Exception("Number of weights do not match number of modes.")
        weighted_modes = np.reshape(self.normalmodes, (self.nmodes, self.natoms*3))
        sum_weights = np.sum(weight)
        for w in range(self.nmodes):
            weighted_modes[w, :] = weighted_modes[w, :] * weight[w]
        weighted_modes = weighted_modes / sum_weights
        try:
            U, S, Vt = svd(weighted_modes, full_matrices=False)
        except:
            print("SVD did not converged. Eliminating zero columns.")
            deleted_modes = []
            for i in range(self.nmodes):
                n = norm(weighted_modes[:, i])
                print(n)
                if n < 1e-6:
                    weighted_modes = np.delete(weighted_modes, i, axis=1)
                    deleted_modes.append(n)
            print("Deleted modes: ", deleted_modes)
            print("New dimensions: ", weighted_modes.shape)
            U, S, Vt = svd(weighted_modes, full_matrices=False)

        return U, S, Vt"""


