# PDB2TRENT 
rotational-translational entropy calculation from thermodynamic ensembles 

The program computes the rotational-translational entropy of one molecule
with respect to the other. The input is a PDB with multiple models for the 
molecular complex.
First a superposition is performed to put the complex in a common frame 
where the first molecule is superposed on the first reference frame,
using the atoms specified on the command-line. 
Then rotation-translation distances based on the superposition of the atoms
specified on the command-line are computed, and finally the entropy is computed 
using the nearest-neighbor method.
The theory and details of the method are reported in:

F. Fogolari, C. J. Dongmo Foumthuim, S. Fortuna, M. A. Soler, A. Corazza, G. Esposito
Accurate estimation of the entropy of rotation-translation probability distributions
J. Chem. Theory. Comput., 12, 1-8, 2016.    

============================================================================

In the built-in superposition tool, routines are suitable modifications of
those written by D.L. Theobald, therefore, if you use these routines, please 
cite:
Theobald, D. L. (2005). Rapid calculation of rmsds using a quaternion-based
characteristic polynomial. Acta. Crystallogr. A, 61, 478–480.

============================================================================

INSTALLATION:

The program is compiled with: 

- if OpenMP is installed: 
cc pdb2trent.c -o pdb2trent -lm -fopenmp

- otherwise:
cc pdb2trent.c -o pdb2trent -lm 

RUNNING TRENT

./pdb2trent without arguments will print options available

./pdb2trent expects at least the following arguments (in this order):
 - the name of the input pdb file 
 - the name of the torsion-adjacency definition file
 - a string for indication of atoms used for frame 
 - a string for indication of atoms used to compute rotational-translational entropy 

The string takes the form of multiple ranges separated by spaces. 

Each range is specified as "chain:residue numbers:atom", e.g.: 

"A:12-18,23-24:N,CA,C" for a single range or 

"A:12-18,23,24:N,CA,C A:104-111,115-123:N,CA,C" for two ranges 

residue numbers are ranges as 12-18 or single numbers separated by commas,
chains are a string of letters, e.g AB, and atom names are separated by commas,
with no space. There is only minimal check on the string, so care must be taken 
to follow strictly the expected format.

Other options are listed hereafter

Usage:

./pdb2trent pdb_infile def_infile outfile [Options]

Options:

-n (max k neighbours for listing entropies (20 default))

-nt X (number of threads to be used, if less than 1, e.g. with -nt 0, the program finds the number of threads 
available)

-wp pdb_file (write superimposed structures in pdb_file)

-v (verbose mode)

USAGE EXAMPLES:

--- compute rotational-translational entropy. For frame superposition use chain A, residues 
12-18,23-24,29-35,41-48,54-55,67-73,75-81,88-97,104-111,115-123 and atoms N,CA,C 

To calculate rotation-translation distances use chain B, 
residues 12-18,23-24,29-35,41-48,54-55,67-73,75-81,88-97,104-111,115-123 and atoms N,CA,C

Use 8 threads for parallel computation:

./pdb2trent sample.pdb sample.out "A:12-18,23-24,29-35,41-48,54-55,67-73,75-81,88-97,104-111,115-123:N,CA,C" "B:12-18,23-24,29-35,41-48,54-55,67-73,75-81,88-97,104-111,115-123:N,CA,C" -nt 8 

--- compute rotational-translational entropy using the same atoms as above. Use the maximum number of 
threads available for parallel computation and write superposed structures in file sup.pdb:

./pdb2trent sample.pdb sample.out "A:12-18,23-24,29-35,41-48,54-55,67-73,75-81,88-97,104-111,115-123:N,CA,C" "B:12-18,23-24,29-35,41-48,54-55,67-73,75-81,88-97,104-111,115-123:N,CA,C" -nt -1 -wp sup.pdb

OUTPUT

1) The output lists the rotational-translational entropy in entropic units and then the translational and rotational entropies, computed isolatedly. 

The reference state is: 1 M concentration, random orientation.

The output lists:
- the k^th nearest neighbour (k = 1..20 by default)
- the entropy value (in R units)
- the average distance to the k^th nearest neighbour
- the average log(distance to the k^th nearest neighbour)

sample.pdb is provided here only for demonstrative purposes, and to reduce the 
computational time. Many more conformational samples (in the range of thousands) 
are needed for accurate estimations of entropy.

