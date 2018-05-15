# pdb2trent
Translational-rotational entropy from pdb configurational ensembles 

The program computes the rotational-translational entropy of one molecule
with respect to the other. The input is a PDB file with multiple models for the 
molecular complex.  
First a superposition is performed to put the complex in a common frame,
in other words the global rotation and translation of the complex is removed.
The first molecule (i.e. specific atoms of the first molecule) is superposed 
on the first reference frame, using the atoms specified on the command-line. 
Then rotation-translation distances among all conformational samples based on 
the superposition of the second molecule atoms specified on the command-line, 
and finally the translational-rotational entropy is 
computed using the nearest-neighbor method.

CONTACT:  

Federico Fogolari  
Dipartimento di Scienze Matematiche, Informatiche e Fisiche  
Universita' di Udine  
Via delle Scienze 206  
33100 Udine - Italy  
Tel ++39 0432 494320  
E-mail federico.fogolari@uniud.it  

REFERENCE:  

Please cite:  
F. Fogolari, O. Maloku, C.J. Dongmo Foumthuim, A. Corazza, G. Esposito  
PDB2ENTROPY and PDB2TRENT: conformational and translational-rotational entropy from molecular ensembles  
J. Chem. Inf. Model. (submitted)  

The theory and details of the method are reported in:  
F. Fogolari, C. J. Dongmo Foumthuim, S. Fortuna, M. A. Soler, A. Corazza, G. Esposito  
Accurate estimation of the entropy of rotation-translation probability distributions  
J. Chem. Theory. Comput., 12, 1-8, 2016.    

============================================================================

In the built-in superposition tool, routines are suitable modifications of
those written by D.L. Theobald, therefore, if you use these routines, please 
cite:  
Theobald, D. L.   
Rapid calculation of rmsds using a quaternion-based characteristic polynomial.  
Acta. Crystallogr. A, 61, 478–480, 2005.  

============================================================================

COMPILATION:

The program is compiled with: 

- if OpenMP is installed:  
cc pdb2trent.c -o pdb2trent -lm -fopenmp

- otherwise:  
cc pdb2trent.c -o pdb2trent -lm 

RUNNING PBB2TRENT

./pdb2trent without arguments will print options available

./pdb2trent expects at least the following arguments (in this order):  
 - the name of the input pdb_file with configurational samples between the MODEL and ENDMDL lines  
 - the name of the output file  
 - two strings for indication of atoms used for frame superposition  
   and for entropy calculation.  
The string takes the form of multiple ranges separated by spaces.   
Each range is specified as "chain:residue_numbers:atom", e.g.:  
"A:12-18,23-24:N,CA,C" for a single range and
"A:12-18,23,24:N,CA,C A:104-111,115-123:N,CA,C" for two ranges

Usage:  
./pdb2trent pdb_infile outfile "string1" "string2" [Options]  
string1 (atoms for superposition) format: "chains:r1,r2-r3,r4....:ATNAM1,ATNAM2... chains:r1,r2-r3,r4....:ATNAM1,ATNAM2..."  
string2 (atoms for rot-trans entropy calc.) format: "chain:r1:ATNAM1,ATNAM2,ATNAM3..."  

Options:  
-n (max k neighbours for listing entropies (20 default))  
-b X (length to mix translational and rotational degrees of freedom)  
-s k (use only one snapshot every k snapshots)  
-nt X (number of threads to be used, if less than 1 the program finds the number of threads available)  
-wp pdb_file (write superimposed structures in pdb_file)  
-v (verbose mode)  

USAGE EXAMPLES:  

--- compute rotational-translational entropy. For frame superposition use chain A, residues 
12-18,23-24,29-35,41-48,54-55,67-73,75-81,88-97,104-111,115-123 and atoms N,CA,C  
To calculate rotation-translation distances use chain B, 
residues 12-18,23-24,29-35,41-48,54-55,67-73,75-81,88-97,104-111,115-123 and atoms N,CA,C  
Use 8 threads for parallel computation:  

./pdb2trent sample.pdb sample.out "A:12-18,23-24,29-35,41-48,54-55,67-73,75-81,88-97,104-111,115-123:N,CA,C" \   "B:12-18,23-24,29-35,41-48,54-55,67-73,75-81,88-97,104-111,115-123:N,CA,C" -nt 8 

--- compute rotational-translational entropy using the same atoms as above. Use the maximum number of 
threads available for parallel computation and write superposed structures in file sup.pdb:  

./pdb2trent sample.pdb sample.out "A:12-18,23-24,29-35,41-48,54-55,67-73,75-81,88-97,104-111,115-123:N,CA,C" \   "B:12-18,23-24,29-35,41-48,54-55,67-73,75-81,88-97,104-111,115-123:N,CA,C" -nt -1 -wp sup.pdb

OUTPUT  

1) The output lists the rotational-translational entropy in entropic units
and then the translational and rotational entropies, computed isolatedly.
The reference state is: 1 M concentration, random orientation.  
The output lists:
- the k^th nearest neighbour (k = 1..20 by default) 
- the entropy value (in R units)  
- the average distance to the k^th nearest neighbour  
- the average log(distance to the k^th nearest neighbour)  

sample.pdb is provided here only for demonstrative purposes, and to reduce the 
computational time. Many more conformational samples (in the range of thousands) 
are needed for accurate estimations of entropy.
More example files are available in the Download menu at biophysics.uniud.it
