from argparse import ArgumentParser, Namespace
import os
from pathlib import Path
import CANalyzer.mo

parser = ArgumentParser()

parser.add_argument('filepath', help='Full directory path to log file', type=str)
parser.add_argument('fchk', help='Full directory path to fchk file', type=str)
parser.add_argument('-f', '--filename', help='Custom filename to save created files', type=str)
parser.add_argument("--groups", help='Dictionary of custom atom groupings', type=str)
args: Namespace = parser.parse_args()

directory = os.getcwd()
filename = directory + "/"
if args.filename:
    filename += args.filename
else:
    filename += Path(args.filepath).stem

Analyze = CANalyzer.mo.MO(args.filepath, args.fchk, args.filename, args.groups)
Analyze.start()
Analyze.mulliken_analysis()
Analyze.print_mulliken()






