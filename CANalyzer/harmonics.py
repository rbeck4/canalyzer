import numpy as np
from numpy.linalg import svd, norm
import pandas as pd
from CANalyzer.load import Load
import CANalyzer.utilities as util
import subprocess
from math import factorial
from numpy.polynomial.hermite import Hermite

class Harmonics(Load):
    def __init__(self, logfile, fchkfile):
        Load.__init__(self, logfile, fchkfile, None, None, None)
        #Load.__init__(self, logfile, fchkfile, None, None, None)
        self.normalmodes = None # rank-3 tensor with indices (mode, atom, xyz)
        self.geometry = None
        self.nmodes = None
        self.harmonic_data = None
        # nmodes x 4: col 0 frequency (au), col 1 reduced mass (au), col 2 force constant (au), col 3 IR Intensity (km/mol)

        self.start()


    def start(self):
        super().start()

        self.geometry = self.readfchk_matrix(r"Current cartesian coordinates", 3, self.natoms).T / 1.88972612457
        self.nmodes = int(str(subprocess.check_output(f"grep 'Number of Normal Modes' {self.fchkfile}", shell=True)).split("\\")[-2].split()[-1])

        self.normalmodes = self.readfchk_matrix(r"Vib-Modes", self.natoms*3*self.nmodes, 1)
        self.normalmodes = np.reshape(self.normalmodes, (self.nmodes, self.natoms, 3))

        self.harmonic_data = self.readfchk_matrix(r"Vib-E2",  self.nmodes, 4)
        self.harmonic_data[:, 0] = self.harmonic_data[:, 0] / 219474.63
        self.harmonic_data[:, 1] = self.harmonic_data[:, 1] * 1822.8885
        self.harmonic_data[:, 2] = self.harmonic_data[:, 2] * 0.0642199266


    def distort(self, mode, amount, geometry=None):
        # bug here not allowing ints as argument
        if type(mode) != type(amount) or len(mode) != len(amount):
            raise Exception("Inconsistent mode and amount length")

        if geometry:
            new_geometry = geometry
        else:
            new_geometry = self.geometry.copy()

        if type(mode) == int:
            new_geometry += amount * self.normalmodes[mode, :, :]
        elif type(mode) == list:
            for i in range(len(mode)):
                m = mode[i]
                new_geometry += amount[i] * self.normalmodes[m, :, :]

        return new_geometry


    def harmonic_wavefunction(self, mode, v=0):
        # Returned wave function input must be in Bohr
        redmass = self.harmonic_data[mode, 1]
        freq = self.harmonic_data[mode, 0]
        if freq <= 0:
            raise Exception("Frequency must be positive")
        factor = np.sqrt((1 / (factorial(v)*(2**v)))) * (redmass*freq/np.pi)**(0.25)

        hermite_weights = np.zeros(v+1)
        hermite_weights[v] = 1
        hermite = Hermite(hermite_weights)

        def wavefunction(x):
            return factor*np.exp(-redmass*freq*(x**2)/2)*hermite(np.sqrt(redmass*freq)*x)

        return wavefunction


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


