MO analysis from Gaussian SCF calculations. Need Gaussian log and fchk.\
\
How to use:\
python \_\_main\_\_.py \<logfile\> \<fchkfile\>
\
Flags:\
-f output file name. Default is orbitals.txt\
--groups groups atoms and returns MO analysis with respect to groups instead of atom type\
Example: \
--group "{'methyl':'[0:3]', 'methene':'[4:6,7:9]', 'O':'[10:10]','Acidic H':'[11:11]'}"\
defines atoms 0-3 as "methyl", 4-6 and 7-9 as 'methene', atom 10 as "O", and atom 11 as "Acidic H"
