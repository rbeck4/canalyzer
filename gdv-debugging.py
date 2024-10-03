from CANalyzer.load import Load
from argparse import ArgumentParser, Namespace
import numpy as np

mat_string1 = "X2C Renormalization 1st Derivative (Real) Component:    1"
mat_string2 = "R(r) matrix derivatives for IMat=    1 (real):"
"""
Numerical vs analytical derivative indexing:
4 : 1    : xx
5 : 2, 4 : xy, yx
6 : 5    : yy
7 : 3, 7 : xz, zx
8 : 6, 8 : yz, zy
9 : 9    : zz
"""

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

matsub = np.abs(mat1 - mat2)
rms = np.sqrt(np.mean(np.multiply(matsub, matsub)))
maximum = np.max(matsub)

print("RMS Difference: ", rms)
print("Max Difference: ", maximum, "Index: ", np.unravel_index(np.argmax(matsub, axis=None), matsub.shape))
print("RMS Element Analytical: ", np.sqrt(np.mean(np.multiply(mat1, mat1))))
print("Max Element Analytical: ", np.max(mat1))
print("RMS Element Numerical: ", np.sqrt(np.mean(np.multiply(mat2, mat2))))
print("Max Element Numerical: ", np.max(mat2))



