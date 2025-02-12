from argparse import ArgumentParser, Namespace
import os
import ast
import CANalyzer.mo
import CANalyzer.projection
import CANalyzer.natorb
import CANalyzer.transfer

parser = ArgumentParser()

valid_jobtypes = ["MOANALYZER", "PROJECTION", "NATORB", "SWAPMO", "MOVEMO"]

jobtype_str = """
    jobtype: MOAnalyzer - provides AO contribution to MOs partitioned by angular momentum and atom (groups of atoms)
             Projection - Provides projection of MOs in fchk onto MOs in fchk2 (Only for Gaussian)
             NatOrb - Computes natural orbitals and stores it into a new fchk (Only for Gaussian)
             SwapMO - Moves MO ordering
             MoveMO - Move MOs stored in CQ/GDV bin/fchk to GDV/CQ fchk/bin
                    - CQ -> GDV creates a new fchk
                    - GDV -> CQ overwrites input bin
                    - make sure that basis, geometry, complex/real, R/U/GHF is consistent between jobs
                    - NYI for UHF 
"""
parser.add_argument('jobtype', choices = valid_jobtypes, help=jobtype_str, type=str.upper)
parser.add_argument('--log', help='Full directory path to log file or out file', type=str)
parser.add_argument('--fchk', help='Full directory path to fchk or bin file', type=str)
parser.add_argument('--filename', help='Custom filename to save created files', type=str)
parser.add_argument('--log2', help='Full directory path to log file or out file', type=str)
parser.add_argument('--fchk2', help='Full directory path to fchk file', type=str)
parser.add_argument('--ml', help='True: MOAnalyzer separates based on ml number as well', type=bool)
parser.add_argument('--grouptotal', help='True: MOAnalyzer returns total AO contribution from group')
parser.add_argument('--renormalize_negatives', help="Treat negative populations as positive and renormalize total population to 1")
parser.add_argument("--groups", help='Dictionary of custom atom groupings or range of states for NatOrb', type=str)
"""
    groups: MOAnalyzer - Make groups of atoms you want to group together in output. If you have atoms (listed in order)
                         C, F, Cl, Br, C, H, H, C, O, O, H
                         and you want to group the Halogens, Carbons, and Others.
                         --groups="{'Halogens':'[2:4]', 'Carbons':'[1:1,5:5,8:8]', 'Others':'[6:7,9:11]'}"
                         
            NatOrb - CI state which natural orbitals are generated from
"""
parser.add_argument("--displaywidth", help='Display width of output file before skipping line', type=int)
parser.add_argument("--swaps", help='List of tuples of MO indices to swap', type=str)
args: Namespace = parser.parse_args()
args.jobtype = args.jobtype.upper()

# setting defaults
directory = os.getcwd()
filename = directory + "/"
if args.filename:
    if args.filename[0] == "/":
      filename = args.filename
    else:
      filename += args.filename

if args.displaywidth:
    displaywidth = args.displaywidth
else:
    displaywidth = None

# running jobs
if args.jobtype == "MOANALYZER":
    if not args.filename:
        filename += 'orbitals.txt'
    if not args.ml:
        args.ml = False
    MOAnalyzer = CANalyzer.mo.MO(args.log, args.fchk, filename, args.groups, displaywidth, args.ml, args.grouptotal, args.renormalize_negatives)
    MOAnalyzer.start()
    MOAnalyzer.mulliken_analysis()
    MOAnalyzer.print_mulliken()
elif args.jobtype == 'PROJECTION':
    if not args.filename:
        filename += 'projections.txt'
    Projection = CANalyzer.projection.Projection(args.log, args.fchk, args.fchk2, filename, displaywidth)
    Projection.start()
    Projection.project()
    Projection.print_project()
elif args.jobtype == 'NATORB':
    if not args.filename:
        filename += 'natorb'
    NatOrb = CANalyzer.natorb.NaturalOrbitals(args.log, args.fchk, filename, args.groups, displaywidth)
    NatOrb.start()
    NatOrb.compute_natorb()
elif args.jobtype == 'SWAPMO':
    swaps = ast.literal_eval(args.swaps)
    MoveMO = CANalyzer.transfer.Tranfer(args.log, args.fchk, args.log2, args.fchk2)
    MoveMO.swapmo(swaps)
elif args.jobtype == 'MOVEMO':
    MoveMO = CANalyzer.transfer.Tranfer(args.log, args.fchk, args.log2, args.fchk2)
    MoveMO.movemo()
else:
    print(args.jobtype + " is an Invalid jobtype")






