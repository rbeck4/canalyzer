from CANalyzer.load import Load
from argparse import ArgumentParser, Namespace
import numpy as np

mat_string1 = "X2C Decoupling Derivative --- Imag --- Component:                         3"
mat_string2 = "X(r) matrix derivatives for IMat=    3 (imaginary):"
nrows = 44
ncol = 44

parser = ArgumentParser()
parser.add_argument('--log', help='Full directory path to log file or out file', type=str)
parser.add_argument('--fchk', help='Full directory path to fchk or bin file', type=str)
parser.add_argument('--log2', help='Full directory path to log file or out file', type=str)
parser.add_argument('--fchk2', help='Full directory path to fchk file', type=str)
args: Namespace = parser.parse_args()

job1 = Load(args.log, args.fchk, None, None, None)
job2 = Load(args.log2, args.fchk2, None, None, None)

job1.start()
job2.start()

mat1 = job1.readlog_matrix(mat_string1, nrows, ncol)
mat2 = job2.readlog_matrix(mat_string2, nrows, ncol)

np.savetxt("mat1.csv", mat1, delimiter=",")
np.savetxt("mat2.csv", mat2, delimiter=",")

matsub = np.abs(mat1 - mat2)
rms = np.sqrt(np.mean(np.multiply(matsub, matsub)))
maximum = np.max(matsub)

print("RMS Difference: ", rms)
print("Max Difference: ", maximum, "Index: ", np.unravel_index(np.argmax(matsub, axis=None), matsub.shape))




