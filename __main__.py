from argparse import ArgumentParser, Namespace
import os
import CANalyzer.mo
import CANalyzer.projection

parser = ArgumentParser()
parser.add_argument('jobtype', help='Type of job to perform', type=str)
"""
    jobtype: MOAnalyzer - provides AO contribution to MOs partitioned by angular momentum and atom (groups of atoms)
             Projection - Provides projection.py of MOs in fchk onto MOs in fchk2
"""
parser.add_argument('logfile', help='Full directory path to log file', type=str)
parser.add_argument('fchk', help='Full directory path to fchk file', type=str)
parser.add_argument('fchk2', help='Full directory path to fchk file', type=str)
parser.add_argument('-f', '--filename', help='Custom filename to save created files', type=str)
parser.add_argument("--groups", help='Dictionary of custom atom groupings', type=str)
parser.add_argument("--displaywidth", help='Display width of output file before skipping line', type=int)
args: Namespace = parser.parse_args()

# setting defaults
directory = os.getcwd()
filename = directory + "/"
if args.filename:
    filename += args.filename

if args.displaywidth:
    displaywidth = args.displaywidth
else:
    displaywidth = 100000

# running jobs
if args.jobtype == "MOAnalyzer":
    if not args.filename:
        filename += 'orbitals.txt'
    MOAnalyzer = CANalyzer.mo.MO(args.logfile, args.fchk, filename, args.groups, displaywidth)
    MOAnalyzer.start()
    MOAnalyzer.mulliken_analysis()
    MOAnalyzer.print_mulliken()
elif args.jobtype == 'Projection':
    if not args.filename:
        filename += 'projections.txt'
    Projection = CANalyzer.projection.Projection(args.logfile, args.fchk, args.fchk2, filename, displaywidth)
    Projection.start()
    Projection.project()
    Projection.print_project()
else:
    print("Invalid jobtype")






