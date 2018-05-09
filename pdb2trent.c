/***** includes *********************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/***** defines *********************/
#define HUGE_INT -999999
#define MAX_N_SUP_AT 5000
#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_set_num_threads(num_threads) 0
#endif


/****** structs ********/

/* chain, residue, atom for superposition, i is the atom index */
struct Cra {
char c;
int r;
char a[5];
int i;
} Cra;

/* structs for the molecular system */
struct System { 
int n_atoms;
struct Atom *atoms;
int n_residues;
struct Residue *residues;
int n_chains;
struct Chain *chains;
int n_segments;
struct Segment *segments;
int n_models;
struct Model *models;
}  System;


struct Atom {
char at_name[5]; 
char alt_loc;
char res_name[4];
char chain; 
char element[3]; 
int model;
int at_n;
int res_n;
char res_ins; 
double coor[3]; 
double occ,temp;
char segid[5];
char pdb_chrg[3]; 
} Atom;

struct Residue { 
char res_name[4];
char chain;
char res_ins;
int res_n;
int model;
char seg_name[5];
int beg, end ;
int n_alt_loc;  
int prev, next;
} Residue;

struct Chain { char ch_name;
int beg, end;
int model;
char seg_name[5];
} Chain;

struct Segment { char seg_name[5];
int beg, end;
int model;
}  Segment;

struct Model { int model_n;
int beg, end;
}  Model;
/* struct for program parameters and options */
struct Flag_par {  
	       char file_in_pdb[120];
               char file_out[120];
               char pdb_outfile[120];
               double bond;
               char s1[1024],s2[1024];
               int n, flag_n;
               int skip;
               int nt;
               int wp;
               int verbose;
	     } Flag_par;

/* function prototypes */
int omp_thread_count();
void read_PDB_atoms(FILE *fp1, int *n_atoms, struct Atom *atoms, int skip);
void read_atom_pdb(char *buf, struct Atom *atom);
double distv(double *r1, double *r2);
void init_flag_par(struct Flag_par *flag_par);
void check_cmd_line(int argc, char *argv[], struct Flag_par *flag_par);
void print_info_flag_par(struct Flag_par flag_par);
int comp(const void * elem1, const void * elem2);
void copy_atom(struct Atom *atom_to, struct Atom *atom_from);
void make_system(struct System *system);
void eliminate_alt_loc(struct System *system);
double torsionv(double *x1, double *x2, double *x3, double *x4);
void hpsort(struct Atom *ra, int n);
void suppos(double **ref, double **mob, int len, double *t, double **R, double *rmsd, double **moball, int lenall);
void quat_to_rot(double *q, double **R);
void CalcQuarticCoeffs(double **coords1, double **coords2, int len, double *coeff);
void CalcQuarticCoeffs_M(double **coords1, double **coords2, int len, double *coeff, double **M);
double CoordsInnerProd(double **coords, int len);
double QCProot(double *coeff, double guess, double delta);
double eval_horn_NR_corrxn(double *c, double x);
double eval_horn_quart(double *c, double x);
double eval_horn_quart_deriv(double *c, double x);
void parse(struct Cra *cra, int *ncra, char *s);
FILE *file_open(char *fname,char *acc);
void diffv(double *v, double *r2, double *r1);
double modv(double *x1);


int main(int argc, char *argv[]) {
      FILE *fp_in_1, *fp_in_2, *fp_out_1;
      char buf[1024], tmpres[5];
      struct Flag_par flag_par;
      struct System system;
      struct Cra *cra;
      int ncra;
      int i, j, k, l, m, n_tors,n_atoms_per_model,n_res_per_model,K, num_threads,nsupat;
      double c,L;
      double *dt, *dr, *dtr, V; 
      double *ent_k_tr, *d_mean_tr, *ld_mean_tr;
      double *ent_k_t, *d_mean_t, *ld_mean_t;
      double *ent_k_r, *d_mean_r, *ld_mean_r;
      double *ent_k, *d, *d_mean, *ld_mean, *ent_k_tot, **h1, **h2,logdk;
      double **R, *t, **coords1,**coords2, rmsd, **coords, ***Ra, **ta, tt;
     
/* Initialize options in struct flag_par */ 
        init_flag_par(&flag_par);      
/* Check command line for options */
        check_cmd_line(argc, argv, &flag_par);      
        if(flag_par.nt < 1)
        num_threads = omp_thread_count();
        else num_threads = flag_par.nt;
        omp_set_num_threads(num_threads);
        print_info_flag_par(flag_par);      
        srand48(1);

/* Count atoms in pdb file, allocate memory and read atoms */
	fp_in_1 = file_open(flag_par.file_in_pdb,"r");
        system.n_atoms = 0;
        while(fgets(buf,120,fp_in_1) != NULL )
            if(!strncmp(buf,"ATOM  ",strlen("ATOM  "))) system.n_atoms++;
        printf("%i atoms found in file %s\n",system.n_atoms,flag_par.file_in_pdb); 
         
        system.atoms=calloc(system.n_atoms, sizeof(struct Atom));
        if(system.atoms != NULL) 
        printf("memory allocated...\n");
        else
        {printf("I could not allocate memory for %i atoms... exiting...\n", system.n_atoms); exit(0);}
        fclose(fp_in_1); 
 
	fp_in_1 = file_open(flag_par.file_in_pdb,"r");
        read_PDB_atoms(fp_in_1, &(system.n_atoms), system.atoms,flag_par.skip);  
        fclose(fp_in_1); 
        printf("atoms stored...\n");
/* Create the structure of the system, by sorting atoms, residues, chains and segments */
        make_system(&system); 
printf("atoms sorted ...\n");

        /* allocate for strings to define atoms for superposition */
        cra = calloc(MAX_N_SUP_AT, sizeof(struct Cra));

/* check that the number of residues and atoms modulo the number of models is zero and compute number of residues and atoms per model */
         if(system.n_residues % system.n_models == 0 && system.n_atoms % system.n_models == 0) 
               {
               n_res_per_model = system.n_residues / system.n_models;
               n_atoms_per_model = system.n_atoms / system.n_models;
               printf("OK: n. atoms and residues modulo n_models == 0\n");
               printf("n. residues per model = %i\n",n_res_per_model);
               }
          else {printf("n. atoms (%i) and residues (%i) modulo n_models (%i) != 0... exiting\n",system.n_atoms, system.n_residues,system.n_models); exit(0);}

/* parse the string defining atoms for superposition */
parse(cra,&ncra,flag_par.s1);
       for(k=0,l=0;k<ncra;k++)
       cra[k].i = -1;

/* prepare for superposition */
    R = malloc(3 * sizeof(double *));
    for(i=0; i <3; i++)
        R[i] = malloc(3 * sizeof(double));
        t = malloc(3 * sizeof(double));
    coords1 = malloc(n_atoms_per_model * sizeof(double *));
    coords2 = malloc(n_atoms_per_model * sizeof(double *));
    coords = malloc(n_atoms_per_model * sizeof(double *));
    for(i=system.models[0].beg, l =0; i<= system.models[0].end; i++)
       for(k=0;k<ncra;k++)
         if(
           (system.atoms[i].chain == cra[k].c)  &&
           (system.atoms[i].res_n == cra[k].r)  &&
           !strcmp(system.atoms[i].at_name,cra[k].a)  
           )
           cra[k].i = i;
       for(k=0,l=0;k<ncra;k++)
       if(cra[k].i != -1)  
       coords1[l++] = system.atoms[cra[k].i].coor;
    nsupat = l;

/* superimpose models to the first one using the first set of atoms to remove global motions */
    for(j=0; j<system.n_models; j++)
      { 
        l = 0;
        for(k=0;k<ncra;k++)
         if(cra[k].i != -1)
         coords2[l++] = system.atoms[cra[k].i + j * n_atoms_per_model].coor;
      for(i=system.models[j].beg; i<= system.models[j].end; i++)
      coords[i-system.models[j].beg] = system.atoms[i].coor;
      suppos(coords1, coords2, nsupat, t, R, &rmsd, coords, n_atoms_per_model);
      }
/* if wanted write out superimposed models */
     if(flag_par.wp)
     {
     fp_out_1=fopen(flag_par.pdb_outfile,"w");
     for(j=0; j<system.n_models; j++)
      {
      fprintf(fp_out_1,"MODEL  %5i\n",j); 
      for(i=system.models[j].beg; i<= system.models[j].end; i++)
        if((int) strlen(system.atoms[i].at_name) == 1)
	fprintf(fp_out_1,"ATOM  %5i  %1s  %c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2lf%6.2lf      %4s%2s%2s\n",
             system.atoms[i].at_n, system.atoms[i].at_name, system.atoms[i].alt_loc, system.atoms[i].res_name, system.atoms[i].chain,
             system.atoms[i].res_n, system.atoms[i].res_ins,system.atoms[i].coor[0], system.atoms[i].coor[1], system.atoms[i].coor[2], 
             system.atoms[i].occ, system.atoms[i].temp, system.atoms[i].segid, system.atoms[i].element, system.atoms[i].pdb_chrg); 
        else if((int) strlen(system.atoms[i].at_name) == 2)
	fprintf(fp_out_1,"ATOM  %5i  %2s %c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2lf%6.2lf      %4s%2s%2s\n",
             system.atoms[i].at_n, system.atoms[i].at_name, system.atoms[i].alt_loc, system.atoms[i].res_name, system.atoms[i].chain,
             system.atoms[i].res_n, system.atoms[i].res_ins,system.atoms[i].coor[0], system.atoms[i].coor[1], system.atoms[i].coor[2], 
             system.atoms[i].occ, system.atoms[i].temp, system.atoms[i].segid, system.atoms[i].element, system.atoms[i].pdb_chrg); 
        else if((int) strlen(system.atoms[i].at_name) == 3)
	fprintf(fp_out_1,"ATOM  %5i  %3s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2lf%6.2lf      %4s%2s%2s\n",
             system.atoms[i].at_n, system.atoms[i].at_name, system.atoms[i].alt_loc, system.atoms[i].res_name, system.atoms[i].chain,
             system.atoms[i].res_n, system.atoms[i].res_ins,system.atoms[i].coor[0], system.atoms[i].coor[1], system.atoms[i].coor[2], 
             system.atoms[i].occ, system.atoms[i].temp, system.atoms[i].segid, system.atoms[i].element, system.atoms[i].pdb_chrg); 
        else
	fprintf(fp_out_1,"ATOM  %5i %4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2lf%6.2lf      %4s%2s%2s\n",
             system.atoms[i].at_n, system.atoms[i].at_name, system.atoms[i].alt_loc, system.atoms[i].res_name, system.atoms[i].chain,
             system.atoms[i].res_n, system.atoms[i].res_ins,system.atoms[i].coor[0], system.atoms[i].coor[1], system.atoms[i].coor[2], 
             system.atoms[i].occ, system.atoms[i].temp, system.atoms[i].segid, system.atoms[i].element, system.atoms[i].pdb_chrg); 
      fprintf(fp_out_1,"ENDMDL\n"); 
      }
     fclose(fp_out_1);
     }

/* prepare for finding relative rotations and translations of the second set of atoms */
Ra = calloc(system.n_models, sizeof(double **));
ta = calloc(system.n_models, sizeof(double *));
for(i=0; i< system.n_models; i++)
{
Ra[i] = calloc(3, sizeof(double *));
ta[i] = calloc(3, sizeof(double));
for(j=0; j< 3; j++)
Ra[i][j] = calloc(3, sizeof(double));
}

    for(i=0; i <3; i++)
        free(R[i]);
    free(R);
    free(t);
ncra = 0;
/* parse second set of atoms */
parse(cra,&ncra,flag_par.s2);

/* find atoms */
       for(k=0;k<ncra;k++)
       cra[k].i = -1;
    for(i=system.models[0].beg, l=0; i<= system.models[0].end; i++)
       for(k=0;k<ncra;k++)
         if(
           (system.atoms[i].chain == cra[k].c)  &&
           (system.atoms[i].res_n == cra[k].r)  &&
           !strcmp(system.atoms[i].at_name,cra[k].a)  
           )
         cra[k].i = i;

       for(k=0,l=0;k<ncra;k++)
       if(cra[k].i != -1)  
       coords1[l++] = system.atoms[cra[k].i].coor;
       nsupat = l;
    for(i=0; i< n_atoms_per_model; i++)
        coords[i] = calloc(3, sizeof(double));
/* find rotation and translation to superimpose atoms on the first model */
    for(j=0; j < system.n_models; j++)
      { 
        l = 0;
        R = Ra[j];
        t = ta[j];
        for(k=0;k<ncra;k++)
         if(cra[k].i != -1)
         coords2[l++] = system.atoms[cra[k].i + j * n_atoms_per_model].coor;

        nsupat = l;
       suppos(coords1, coords2, nsupat, t, R, &rmsd, coords, n_atoms_per_model);       
      }
/*... calculate the distance between rotation-translations 
and calculate entropy based on the nearest neighbor method */
t = NULL;
R = NULL;
t = calloc(3, sizeof(double));
R = calloc(3, sizeof(double *));
    for(i=0; i< 3; i++)
     R[i] = calloc(3, sizeof(double));
    tt = flag_par.bond * flag_par.bond;
K = flag_par.n + 1;
if(K > system.n_models) 
{
K = system.n_models;
flag_par.n = K - 1;
}
ent_k_r = calloc(K-1,sizeof(double));
d_mean_r = calloc(K,sizeof(double)); 
ld_mean_r = calloc(K,sizeof(double));
ent_k_t = calloc(K-1,sizeof(double));
d_mean_t = calloc(K,sizeof(double)); 
ld_mean_t = calloc(K,sizeof(double));
ent_k_tr = calloc(K-1,sizeof(double));
d_mean_tr = calloc(K,sizeof(double)); 
ld_mean_tr = calloc(K,sizeof(double));

#pragma omp parallel for num_threads(num_threads) private(i,j,k,l,dt,dr,dtr,V) shared(ent_k_t, d_mean_t, ld_mean_t,ent_k_r, d_mean_r, ld_mean_r,ent_k_tr, d_mean_tr, ld_mean_tr)    
for(i=0; i< system.n_models; i++)
      {
dt = calloc(system.n_models,sizeof(double));
dr = calloc(system.n_models,sizeof(double));
dtr = calloc(system.n_models,sizeof(double));
      for(j=0,dtr[i] = 0.0,dt[i] = 0.0, dr[i] = 0.0; j< system.n_models; j++)
        if(i != j)
         {
         dt[j] = distv(ta[i],ta[j]); 
            
         for(k=0, dr[j] = 0.0; k<3; k++)
         {
         for(l=0; l<3; l++)
         dr[j] = dr[j] + Ra[i][k][l] * Ra[j][k][l];
         }
         if(fabs(dr[j]) <= 3.0)
         dr[j] = acos(0.5 * (  dr[j] - 1));
         else if(dr[j] > 0) dr[j] = 0.0;
         else if(dr[j] < 0) dr[j] = M_PI;
         dtr[j] = sqrt(dt[j] * dt[j] + tt * dr[j] *dr[j]);
         }
         qsort(dt,system.n_models,sizeof(double), comp);
         qsort(dr,system.n_models,sizeof(double), comp);
         qsort(dtr,system.n_models,sizeof(double), comp);

for(k = 1; k<=flag_par.n; k++)
{
/* apply the approximation to the volume of a ball in the rotation-translation space */
V = M_PI * M_PI * M_PI * pow(dtr[k],3.0) * 
   (  pow(dtr[k]/flag_par.bond,3.0) / 12.0  
    - pow(dtr[k]/flag_par.bond,5.0) / 384.0 
    + pow(dtr[k]/flag_par.bond,7.0) / 23040.0 
    - pow(dtr[k]/flag_par.bond,9.0) / 2211840.0 
    + pow(dtr[k]/flag_par.bond,11.0) / 309657600.0 ) ;
#pragma omp critical
{
ent_k_tr[k-1] = ent_k_tr[k-1] + log(V);
d_mean_tr[k] = d_mean_tr[k] + dtr[k];
ld_mean_tr[k] = ld_mean_tr[k] + log(dtr[k]);
V = 4.0 * M_PI * pow(dt[k],3.0)/3.0;
ent_k_t[k-1] = ent_k_t[k-1] + log(V);
d_mean_t[k] = d_mean_t[k] + dt[k];
ld_mean_t[k] = ld_mean_t[k] + log(dt[k]);
V = 4.0 * M_PI * (dr[k] - sin(dr[k]));
ent_k_r[k-1] = ent_k_r[k-1] + log(V);
d_mean_r[k] = d_mean_r[k] + dr[k];
ld_mean_r[k] = ld_mean_r[k] + log(dr[k]);
}
}
free(dtr);
free(dt);
free(dr);
}
    free(system.atoms); 
    free(system.residues);
    free(system.chains);
    free(system.models);
    free(system.segments);
fp_out_1 = file_open(flag_par.file_out,"w");
for(k = 1; k <= flag_par.n; k++)
{
ent_k_tr[k-1] = ent_k_tr[k-1]/(double) (system.n_models -1);
d_mean_tr[k] = d_mean_tr[k]/(double) (system.n_models -1);
ld_mean_tr[k] = ld_mean_tr[k]/(double) (system.n_models -1);
ent_k_t[k-1] = ent_k_t[k-1]/(double) (system.n_models -1);
d_mean_t[k] = d_mean_t[k]/(double) (system.n_models -1);
ld_mean_t[k] = ld_mean_t[k]/(double) (system.n_models -1);
ent_k_r[k-1] = ent_k_r[k-1]/(double) (system.n_models -1);
d_mean_r[k] = d_mean_r[k]/(double) (system.n_models -1);
ld_mean_r[k] = ld_mean_r[k]/(double) (system.n_models -1);
}
c = 0.5722 + log( (double) (system.n_models -1));
for(k = 1,L=0; k <= flag_par.n; k++)
{
ent_k_tr[k-1] = ent_k_tr[k-1] + c - L;
ent_k_t[k-1] = ent_k_t[k-1] + c - L;
ent_k_r[k-1] = ent_k_r[k-1] + c - L;
L = L + 1.0/(double) k;
}
/* output rotation-translation and rotation-only and translation only entropy */
fprintf(fp_out_1,"Entropy in R units, reference state 1M, random orientation\n");
for(k = 0; k<flag_par.n ; k++)
fprintf(fp_out_1,"ROT. TRANS. (R units): ent_k %3i %e d_mean %e ld_mean %e\n", k+1, ent_k_tr[k] -2.0*log(2*M_PI) -7.414898, d_mean_tr[k+1],ld_mean_tr[k+1]);
for(k = 0; k<flag_par.n ; k++)
fprintf(fp_out_1,"TRANS. (R units)     : ent_k %3i %e d_mean %e ld_mean %e\n", k+1, ent_k_t[k] -7.414898, d_mean_t[k+1],ld_mean_t[k+1]);
for(k = 0; k<flag_par.n ; k++)
fprintf(fp_out_1,"ROT. (R units)       : ent_k %3i %e d_mean %e ld_mean %e\n", k+1, ent_k_r[k] -2.0*log(2*M_PI), d_mean_r[k+1],ld_mean_r[k+1]);

if(system.n_models >= 11) k = 10;
else
k = system.n_models - 1;
fprintf(fp_out_1,"\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n");
fprintf(fp_out_1,"Total entropy estimate on %5i^th nearest neighbours:\n%10.2lf R units\n", k, ent_k_tr[k-1] -2.0*log(2*M_PI) -7.414898);
fprintf(fp_out_1,"\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
printf("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n");
printf("Total entropy estimate on %5i^th nearest neighbours:\n%10.2lf R units\n", k, ent_k_tr[k-1] -2.0*log(2*M_PI) -7.414898);
printf("\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");


fclose(fp_out_1);
}

/* read the command line */
void check_cmd_line(int argc, char *argv[], struct Flag_par *flag_par)
{
	int i;
        char tmp[100];
        char extension[100];

	if(argc < 5) 
	{
	printf("Usage:\n"); 
	printf("./trent pdb_infile outfile \"string1\" \"string2\" [Options]\n"); 
        printf("string1 (atoms for superposition) format: \"chains:r1,r2-r3,r4....:ATNAM1,ATNAM2... chains:r1,r2-r3,r4....:ATNAM1,ATNAM2...\"\n");
        printf("string2 (atoms for rot-trans entropy calc.) format: \"chain:r1:ATNAM1,ATNAM2,ATNAM3...\"\n");

	printf("Options:\n"); 
	printf("-n (max k neighbours for listing entropies (20 default))\n"); 
	printf("-b X (length of bond to mix translational and rotational degrees of freedom)\n"); 
        printf("-s k (keep only one sample every k samples)\n");
	printf("-nt X (number of threads to be used, if less than 1 the program finds the number of threads available)\n"); 
	printf("-wp pdb_file (write superimposed structures in pdb_file)\n"); 
	printf("-v (verbose mode)\n"); 
	printf("\n"); 
	exit(1);
	}

        strcpy((*flag_par).file_in_pdb, argv[1]);
        strcpy((*flag_par).file_out, argv[2]);
        strcpy((*flag_par).s1, argv[3]);
        strcpy((*flag_par).s2, argv[4]);


	for (i = 5; i < argc; i++) {
		if (!strncmp(argv[i],"-v",3)) (*flag_par).verbose = 1;
		else if (!strncmp(argv[i],"-b",3)) (*flag_par).bond = atoi(argv[++i]);
		else if (!strncmp(argv[i],"-nt",4)) (*flag_par).nt = atoi(argv[++i]);
   	        else if (!strncmp(argv[i],"-s",3)) (*flag_par).skip = atoi(argv[++i]);
else if (!strncmp(argv[i],"-wp",4)) 
                 {
                 flag_par->wp = 1;
                 strcpy(flag_par->pdb_outfile,argv[++i]);
                 }
		else if (!strncmp(argv[i],"-n",3)) 
                 {
                 (*flag_par).flag_n = 1;
                 (*flag_par).n = atoi(argv[++i]);
                 }
		else 
                {
                 printf("I don't know option %s\n", argv[i]);
                 exit(2);
                }
}
printf("\n########################################################\n\n"); 

}

/* initialize options */
void init_flag_par(struct Flag_par *flag_par)
{
(*flag_par).flag_n=1;
(*flag_par).nt=1;
(*flag_par).skip=1;
(*flag_par).wp=0;
(*flag_par).bond=1.0;
(*flag_par).n=20;
(*flag_par).verbose=0;
}

/* print options */
void print_info_flag_par(struct Flag_par flag_par)
{
	printf("########################################################\n"); 
	printf("# ACCORDING TO YOUR PARAMETER CHOICES:                 #\n"); 
	printf("########################################################\n\n"); 

        printf("pdb file: %s\n", flag_par.file_in_pdb);
        printf("out file: %s\n", flag_par.file_out);
        if(flag_par.nt > 0)
        printf("I will use %i threads\n", flag_par.nt);
        else
        printf("I will use all threads available\n");
        printf("I will use 1 snapshot every %i snapshots\n", flag_par.skip);
        printf("I will print entropy only for the first %i neighbours\n", flag_par.n);
        printf("I will superimpose all structures on the first one\n");
        if(flag_par.wp)
        printf("I will write superimposed structures in file %s\n", flag_par.pdb_outfile);
	if(flag_par.verbose) 
		printf("I will be verbose\n");
	else
		printf("I will NOT be verbose\n");
	printf("########################################################\n\n"); 
}

/* implicit compare function */
int comp (const void * elem1, const void * elem2) {
    double f1 = *((double *)elem1);
    double f2 = *((double *)elem2);
    if (f1 > f2) return  1;
    if (f1 < f2) return -1;
    return 0;
}

/* wrap of fopen() */
FILE *file_open(char *fname,char *acc) {
    FILE *fp;
    fp =fopen(fname,acc);
    if (fp==NULL) 
	{
        fprintf(stderr,"unable to open file \"%s\"... \n",fname);
        exit(1);
    }
    return(fp);
}

/* read coordinates in PDB file */ 
void read_PDB_atoms(FILE *fp1, int *n_atoms, struct Atom *atoms, int skip)
{
	char buf[120];
	int i=0;
        int mod_id=1;
        int j=1;

	while(fgets(buf,120,fp1) != NULL )
    	if(!strncmp("ATOM",buf,4)) 
            {
	    atoms[i].model=j; 
	    read_atom_pdb(buf, &atoms[i]);
            if((mod_id % skip) == 0)
            i++;
            if(i%100000 == 0) printf("%i atoms read...\n",i); 
	    }
	    else
	      if(!strncmp("ENDMDL",buf,6)) 
              {
              mod_id++; 
              if((mod_id % skip) == 0)
              j++;
              }
	*n_atoms = i;
}

/* read atom line in PDB file */
void read_atom_pdb(char *buf, struct Atom *atom)
{

    char at_rec[5];
    char tok[10];              
 
    strncpy(tok,buf,4);
    tok[4] = '\0';
    sscanf(tok,"%s", at_rec);
    if(strncmp("ATOM",at_rec,4))
    {
    	printf("What is supposed to be an ATOM line does not start with string ATOM... exiting...\n");
        exit(1);
    }

    strncpy(tok,buf + 6,5);
    tok[5] = '\0';
    sscanf(tok,"%i",&(atom->at_n));

    strncpy(tok,buf + 12,4);
    tok[4] = '\0';
    sscanf(tok,"%s", atom->at_name);
 
    strncpy(tok,buf + 16,1);
    tok[1] = '\0';
    if(sscanf(tok,"%c", &(atom->alt_loc)) == -1) atom->alt_loc=' ';

	strncpy(tok,buf + 17,3);
    tok[3] = '\0'; 
    sscanf(tok,"%s", atom->res_name);

	strncpy(tok,buf + 21,1);
    tok[1] = '\0';
    if(sscanf(tok,"%c", &(atom->chain)) == EOF) atom->chain = ' ';

    strncpy(tok,buf + 22,4);
    tok[4] = '\0';
    sscanf(tok,"%i", &(atom->res_n));

	strncpy(tok,buf + 26,1);
    tok[1] = '\0';
    if (sscanf(tok,"%c", &(atom->res_ins)) == EOF) atom->res_ins=' ';

    strncpy(tok,buf + 30,8);
    tok[8] = '\0';
    sscanf(tok,"%lf", &(atom->coor[0]));

	strncpy(tok,buf + 38,8);
    tok[8] = '\0';
    sscanf(tok,"%lf", &(atom->coor[1]));

    strncpy(tok,buf + 46,8);
    tok[8] = '\0';
    sscanf(tok,"%lf", &(atom->coor[2]));

    strncpy(tok,buf + 54,6);
    tok[6] = '\0';
    sscanf(tok,"%lf", &(atom->occ));

	strncpy(tok,buf + 60,6);
    tok[6] = '\0';
    sscanf(tok,"%lf", &(atom->temp));
    if(strlen(buf) > 76)
    { 
      	strncpy(tok,buf + 72,4);
       	tok[4] = '\0';
		if (sscanf(tok,"%s", (atom->segid)) == EOF) 
			strcpy(atom->segid,"    ");
    }
    else strcpy(atom->segid,"    ");


    if(strlen(buf) > 78)
	{
       	strncpy(tok,buf + 76,2);
       	tok[2] = '\0';
       	if (sscanf(tok,"%s", (atom->element)) == EOF) 
			strcpy(atom->element,"UN");
	}

    if(strlen(buf) > 80)
	{
       	strncpy(tok,buf + 78,2);
       	tok[2] = '\0';
       	if (sscanf(tok,"%s", (atom->pdb_chrg)) == EOF) 
			strcpy(atom->pdb_chrg,"  ");
	}

}

/* implicit compare function for sorting atoms */
int cmp_atoms(const void *p1, const void *p2)
{
	struct Atom A_atom, B_atom;
	int check = 0 ; 

	A_atom = *((struct Atom *)p1); 
	B_atom = *((struct Atom *)p2); 



	if( A_atom.model < B_atom.model) check = -1; 
	else if (A_atom.model == B_atom.model) 
    {
        check = 0;
    	if( strcmp(A_atom.segid,B_atom.segid) < 0) check = -1; 
    	else if (!strcmp(A_atom.segid,B_atom.segid)) 
        {
        	check = 0;
           	if( (int) A_atom.chain < (int) B_atom.chain ) check = -1;
           	else if( A_atom.chain ==  B_atom.chain )
            {
            	check = 0;
                if( A_atom.res_n <  B_atom.res_n ) check = -1;
                else if(A_atom.res_n ==  B_atom.res_n )
                {
                	check = 0;
                    if((int) A_atom.res_ins < (int) B_atom.res_ins ) check = -1;
                    else if( A_atom.res_ins ==  B_atom.res_ins )
                    {
                    	check = 0;
                        if( strcmp(A_atom.at_name,B_atom.at_name) < 0) check = -1; 
                        else if (!strcmp(A_atom.at_name,B_atom.at_name)) 
                        {
                        check = 0;
                        if( (int) A_atom.alt_loc < (int) B_atom.alt_loc ) check = -1;
                        else if( A_atom.alt_loc ==  B_atom.alt_loc ) check = 0;
                        else check = 1;
			}
                        else check = 1;
                    }
                        else check = 1;
                }
                        else check = 1;
            }
                        else check = 1;
        }
                        else check = 1;
	}
        else check = 1;
	return check;

}

/* Create the structure of the system, by sorting atoms, residues, chains and segments */
void make_system(struct System *system)
{
	int n_atoms = system->n_atoms;
	int *p_n_residues = &(system->n_residues); 
	int *p_n_chains = &(system->n_chains); 
	int *p_n_segments  = &(system->n_segments); 
	int *p_n_models = &(system->n_models); 

	int imodel=-1, isegment=-1, ichain=-1, iresidue=-1;
	char res_ins = '*';
	char chain = '*';
	char segid[5] = "****";
	int model = -1;
	int res_n = -HUGE_INT;

	int i;

	hpsort(system->atoms, system->n_atoms);

        eliminate_alt_loc(system);

        printf("after eliminate_alt_loc %i atoms left\n", system->n_atoms);  

	n_atoms = system->n_atoms;


             iresidue = ichain = isegment = imodel = 0;
	for(i=0; i<n_atoms; i++)
             {
		if(system->atoms[i].model != model)
		{
			imodel++;
			isegment++;
			ichain++;
			iresidue++;
                }
                else if(strcmp(system->atoms[i].segid,segid))
                        {
                                isegment++;
                                ichain++;
                                iresidue++;
                        }  
                 else if(system->atoms[i].chain != chain)
                        {
                                ichain++;
                                iresidue++;
                        }
	         else if(system->atoms[i].res_n != res_n) 
                                iresidue++;
                 else if(system->atoms[i].res_ins != res_ins) 
			        iresidue++;
 
        	model = system->atoms[i].model; 
        	chain = system->atoms[i].chain; 
        	res_n = system->atoms[i].res_n; 
        	res_ins = system->atoms[i].res_ins; 
            	strcpy(segid, system->atoms[i].segid); 
                }
                 system->residues = calloc(iresidue,sizeof(struct Residue));
        if(system->residues == NULL) 
        {printf("I could not allocate memory for %i residues... exiting...\n", iresidue); exit(0);}
                 system->chains = calloc(ichain,sizeof(struct Chain));
                 system->segments = calloc(isegment,sizeof(struct Segment));
                 system->models = calloc(imodel,sizeof(struct Model));
             iresidue = ichain = isegment = imodel = -1;
	res_ins = '*';
	chain = '*';
	strcpy(segid,"****");
	model = -1;
	res_n = -HUGE_INT;
	for (i=0; i<n_atoms; i++)
	{
		if(system->atoms[i].model != model)
		{
			imodel++;
			isegment++;
			ichain++;
			iresidue++;
 
			strcpy(system->residues[iresidue].res_name, system->atoms[i].res_name);
			system->residues[iresidue].res_n = system->atoms[i].res_n;
			system->residues[iresidue].res_ins = system->atoms[i].res_ins;
			system->residues[iresidue].chain = system->atoms[i].chain;
			system->residues[iresidue].model = system->atoms[i].model;
			strcpy(system->residues[iresidue].seg_name, system->atoms[i].segid);
			system->residues[iresidue].beg = i;
			if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
        	chain = system->atoms[i].chain; 
        	res_n = system->atoms[i].res_n; 
        	res_ins = system->atoms[i].res_ins; 
            	strcpy(segid, system->atoms[i].segid); 
        
			model=system->atoms[i].model;        

			system->chains[ichain].ch_name = chain;
			system->chains[ichain].beg = i;
                 	system->chains[ichain].model=system->atoms[i].model;
			strcpy(system->chains[ichain].seg_name, system->atoms[i].segid);
			if (ichain != 0) system->chains[ichain - 1].end = i-1;

			strcpy(system->segments[isegment].seg_name, segid);
        	system->segments[isegment].model=system->atoms[i].model;
        	system->segments[isegment].beg = i;
 			if (isegment != 0) system->segments[isegment - 1].end = i-1;

        	system->models[imodel].beg=i;
        	system->models[imodel].model_n=model;
        	if (imodel != 0) system->models[imodel - 1].end =i-1;

		}
		else if(strcmp(system->atoms[i].segid,segid))
			{
				isegment++;
  				ichain++;
  				iresidue++;

				strcpy(system->residues[iresidue].res_name, system->atoms[i].res_name);
				system->residues[iresidue].res_n = system->atoms[i].res_n;
				system->residues[iresidue].res_ins = system->atoms[i].res_ins;
				system->residues[iresidue].chain = system->atoms[i].chain;
				system->residues[iresidue].model = system->atoms[i].model;
				strcpy(system->residues[iresidue].seg_name, system->atoms[i].segid);
				system->residues[iresidue].beg = i;
				if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
			    chain = system->atoms[i].chain; 
    		    res_n = system->atoms[i].res_n; 
       		 	res_ins = system->atoms[i].res_ins; 
            		strcpy(segid, system->atoms[i].segid); 
        
				system->chains[ichain].ch_name = chain;
				system->chains[ichain].beg = i;
        		system->chains[ichain].model=system->atoms[i].model;
				strcpy(system->chains[ichain].seg_name, system->atoms[i].segid);
				if (ichain != 0) system->chains[ichain - 1].end = i-1;

	      		strcpy(system->segments[isegment].seg_name, segid);
    	    	system->segments[isegment].model=system->atoms[i].model;
       	 		system->segments[isegment].beg = i;
 				if (isegment != 0) system->segments[isegment - 1].end = i-1;

			}
			else if(system->atoms[i].chain != chain) 
				{
					ichain++;
					iresidue++;

					strcpy(system->residues[iresidue].res_name, system->atoms[i].res_name);
					system->residues[iresidue].res_n = system->atoms[i].res_n;
					system->residues[iresidue].res_ins = system->atoms[i].res_ins;
					system->residues[iresidue].chain = system->atoms[i].chain;
					system->residues[iresidue].model = system->atoms[i].model;
					strcpy(system->residues[iresidue].seg_name, system->atoms[i].segid);
					system->residues[iresidue].beg = i;
					if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
        			chain = system->atoms[i].chain; 
       		 		res_n = system->atoms[i].res_n; 
        			res_ins = system->atoms[i].res_ins; 
					system->chains[ichain].ch_name = chain;
					system->chains[ichain].beg = i;
        			system->chains[ichain].model=system->atoms[i].model;
					strcpy(system->chains[ichain].seg_name, system->atoms[i].segid);
					if (ichain != 0) system->chains[ichain - 1].end = i-1;
				}
				else if(system->atoms[i].res_n != res_n) 
					{
					iresidue++;

					strcpy(system->residues[iresidue].res_name, system->atoms[i].res_name);
					res_n = system->atoms[i].res_n;
					res_ins = system->atoms[i].res_ins;
					system->residues[iresidue].res_n = res_n;
					system->residues[iresidue].res_ins = res_ins;
					system->residues[iresidue].chain = chain;
					system->residues[iresidue].model = system->atoms[i].model;
					strcpy(system->residues[iresidue].seg_name, system->atoms[i].segid);
					system->residues[iresidue].beg = i;
					if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
					}
					else if(system->atoms[i].res_ins != res_ins) 
					{
					iresidue++;
					strcpy(system->residues[iresidue].res_name, system->atoms[i].res_name);
					res_n = system->atoms[i].res_n;
					res_ins = system->atoms[i].res_ins;
					system->residues[iresidue].res_n = res_n;
					system->residues[iresidue].res_ins = res_ins;
					system->residues[iresidue].chain = chain;
					system->residues[iresidue].model = system->atoms[i].model;
					strcpy(system->residues[iresidue].seg_name, system->atoms[i].segid);
					system->residues[iresidue].beg = i;
					if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
					}
	}
if(n_atoms != 0)
    {
    system->segments[isegment].end = i-1;
	system->chains[ichain].end = i-1;
	system->residues[iresidue].end = i-1;
    system->models[imodel].end = i-1;
    }
	*p_n_chains = ichain+1;
	*p_n_segments = isegment+1;
	*p_n_residues = iresidue+1;
	*p_n_models = imodel+1;

	printf("########################################################\n"); 
	printf("# SYSTEM:                                              #\n"); 
	printf("########################################################\n\n"); 
	printf("atoms = %8i, residues = %8i, chains = %8i, segments = %8i, models = %8i\n\n", n_atoms, *p_n_residues, *p_n_chains,*p_n_segments,*p_n_models  );
	printf("########################################################\n\n"); 
}

/* heapsort */
void hpsort(struct Atom *ra, int n)
{  
    int N, i, parent, child;  
    struct Atom rra;  
    N = n;
    i = n/2;
    for (;;) { 
        if (i > 0) { 
            i--;           
            copy_atom(&(rra),&(ra[i]));
        } else {     
            n--;           
            if (n == 0) return; 
            copy_atom(&(rra),&(ra[n]));
            copy_atom(&(ra[n]),&(ra[0]));
        }  
  
        parent = i; 
        child = i*2 + 1; 
  
        while (child < n) {  
            if (child + 1 < n  && (cmp_atoms(&(ra[child + 1]),&(ra[child])) > 0)) {
                child++; 
            }  
            if (cmp_atoms(&(ra[child]),&(rra)) > 0) {   
                copy_atom(&(ra[parent]),&(ra[child]));
                parent = child;  
                child = parent*2+1;
            } else {  
                break; 
            }  
        }  
       copy_atom(&(ra[parent]),&(rra));
    }  
}

/* copy struct Atom */
void copy_atom(struct Atom *atom_to, struct Atom *atom_from)
{
int i;
strcpy((*atom_to).at_name,(*atom_from).at_name);
(*atom_to).alt_loc = (*atom_from).alt_loc;
strcpy((*atom_to).res_name,(*atom_from).res_name);
(*atom_to).chain = (*atom_from).chain;
(*atom_to).model = (*atom_from).model;
(*atom_to).at_n = (*atom_from).at_n;
(*atom_to).res_n = (*atom_from).res_n;
(*atom_to).res_ins = (*atom_from).res_ins;
for(i=0;i<3;i++)
(*atom_to).coor[i] = (*atom_from).coor[i];
(*atom_to).occ = (*atom_from).occ;
(*atom_to).temp = (*atom_from).temp;
strcpy((*atom_to).segid,(*atom_from).segid);
strcpy((*atom_to).pdb_chrg,(*atom_from).pdb_chrg);
}

/* if two sorted atoms are the same, e.g. with different occupancies keep the first one */
void eliminate_alt_loc(struct System *system)
{
int i,j;
char buf[120];
j = 1;
for (i=1; i< (*system).n_atoms; i++)
{
if(
strcmp((*system).atoms[i-1].at_name,(*system).atoms[i].at_name) ||
strcmp((*system).atoms[i-1].res_name,(*system).atoms[i].res_name) ||
strcmp((*system).atoms[i-1].segid,(*system).atoms[i].segid) ||
((*system).atoms[i-1].chain != (*system).atoms[i].chain) ||
((*system).atoms[i-1].model != (*system).atoms[i].model) ||
((*system).atoms[i-1].res_n != (*system).atoms[i].res_n) ||
((*system).atoms[i-1].res_ins != (*system).atoms[i].res_ins)
)
{
copy_atom(&((*system).atoms[j]), &((*system).atoms[i]));
(*system).atoms[j].alt_loc = ' ';
j++;
}
}
if ((*system).n_atoms > 0)
(*system).n_atoms = j;
else
{
printf("n. atoms less than 1.... exiting \n");
exit(0);
}
}

/* vector math utilities */
void scale(double *v, double a)
{
	v[0] = v[0] * a;
	v[1] = v[1] * a;
	v[2] = v[2] * a;
}

void diffv(double *v, double *r2, double *r1)
{
	v[0] = r2[0] - r1[0];
	v[1] = r2[1] - r1[1];
	v[2] = r2[2] - r1[2];
}

void sumv(double *v, double *r2, double *r1)
{
	v[0] = r2[0] + r1[0];
	v[1] = r2[1] + r1[1];
	v[2] = r2[2] + r1[2];
}

double dotv(double *x1, double *x2)
{       
    double d;
    d = x1[0]*x2[0] +  x1[1]*x2[1] + x1[2]*x2[2];
    return d;
}

double modv(double *x1)
{       
    double d;
    d = sqrt(x1[0]*x1[0] +  x1[1]*x1[1] + x1[2]*x1[2]);
    return d;
}

void vecv(double *x, double *x1, double *x2)
{
    x[0] = x1[1]*x2[2]-x1[2]*x2[1];
    x[1] = x1[2]*x2[0]-x1[0]*x2[2];
    x[2] = x1[0]*x2[1]-x1[1]*x2[0];
}
 
double distv(double *r1, double *r2)
{
	double d,v[3];
	diffv(v,r2,r1);
	d = dotv(v, v );
	d =  sqrt(d);
	return d;
}

/* superpose two sets of coordinates using and adapting the routines of Theobald */ 
void suppos(double **ref, double **mob, int len, double *t, double **R, double *rmsd, double **moball, int lenall)
{
    double          innerprod;
    double          lambdamax;
    double          coeff[3];
    double          **x1, **x2, c1[3], c2[3], tmp[3];
    double          **M; 
    double          M_copy[4][4]; 
    double          ev[4], v[4], mod;
    double          max;
    int i, j, k, found;
    int ind_max[4];


    M = malloc(4 * sizeof(double *));
    for(i=0; i <4; i++)
        M[i] = malloc(4 * sizeof(double));

    x1 = malloc(len * sizeof(double *));
    x2 = malloc(len * sizeof(double *));
    for(i=0; i <len; i++)
        {
        x1[i] = malloc(3 * sizeof(double));
        x2[i] = malloc(3 * sizeof(double));
        }

    for(j=0; j <3; j++)
    {
    c1[j] = c2[j] = 0.0;
    for(i=0; i <len; i++)
      {
      c1[j] = c1[j] + ref[i][j];
      c2[j] = c2[j] + mob[i][j];
      }
    c1[j] = c1[j] / (double) len; 
    c2[j] = c2[j] / (double) len; 
    for(i=0; i <len; i++)
      {
      x1[i][j] = ref[i][j] - c1[j];
      x2[i][j] = mob[i][j] - c2[j];
      }
    }
    for(j=0; j <3; j++)
    t[j] = c1[j] - c2[j];
    innerprod = CoordsInnerProd(x1, len) + CoordsInnerProd(x2, len);
    CalcQuarticCoeffs_M(x1, x2, len, coeff, M);
    lambdamax = QCProot(coeff, 0.5 * innerprod, 1e-3);
    mod = (innerprod - (2.0 * lambdamax))/(double) len;
    if(mod >= 0.0) *rmsd = sqrt(mod);
    else if(fabs(mod) < 1e-4) *rmsd = 0.0; 
    for(i=0; i<4; i++)
    for(j=0; j<4; j++)
        M_copy[i][j] = M[i][j];

    for(j=0; j<4; j++)
    M[j][j] = M[j][j] - lambdamax;


    for(j=1; j<4; j++) ind_max[j] = -1;
    for(j=1; j<4; j++)
      {
      max = 0;
      for(i=0; i<4; i++)
      if((fabs(M[i][j])) > max) 
      {
      found = 0;
      for(k=1; k<=j; k++)
      if( ind_max[k] == i)  found = 1;
      if(!found)
      {
      max = fabs(M[i][j]); 
      ind_max[j] = i;
      }
      }

      for( k = 0; k < 4; k++) 
      if(k != ind_max[j])
      {
      max = M[k][j];   
      for(i = 0; i<4; i++)
      M[k][i] = M[k][i] - M[ind_max[j]][i] * max/M[ind_max[j]][j];
      }
      }

      for(j=3; j>0; j--)
        {
        ev[j] = -M[ind_max[j]][0]/M[ind_max[j]][j];
        for(k=1; k<j; k++)
           M[ind_max[k]][0] = M[ind_max[k]][0] - M[ind_max[k]][j] *ev[j];
        }
      ev[0] = 1.0;
      mod = 0.0;
      for(j=0; j<4; j++)
      mod += ev[j]*ev[j];
      mod = sqrt(mod);

      for(j=0; j<4; j++)
        ev[j] = ev[j]/mod;

    quat_to_rot(ev,R);

    for(i=0; i <lenall; i++)
    {
    for(j=0; j < 3; j++)
     tmp[j] = moball[i][j] - c2[j];
    for(j=0; j < 3; j++)
     moball[i][j] = 0.0; 
    for(j=0; j < 3; j++)
    for(k=0; k < 3; k++)
    moball[i][j] = moball[i][j] + R[j][k]*tmp[k];
    for(j=0; j < 3; j++)
    moball[i][j] = moball[i][j]  + t[j] + c2[j];
    }
    for(i=0; i <len; i++)
        {
        free(x1[i]); 
        free(x2[i]); 
        }
    free(x1);
    free(x2);

    for(i=0; i <4; i++)
        free(M[i]); 
    free(M);
}



/*******************************************************************************
 *  -/_|:|_|_\- 
 *
 *  File:           QuatCharPoly.c
 *
 *  Function:       Rapid calculation of RMSD using a quaternion-based
 *                  characteristic polynomial
 *
 *  Author(s):      Douglas Theobald
 *                  Department of Chemistry and Biochemistry
 *                  UCB 215
 *                  University of Colorado at Boulder
 *                  Boulder, CO 80309-0215
 *
 *                  theobal@colorado.edu
 *                  dtheobald@hotmail.com
 *
 *  Copyright:      Copyright (c) 2005 Douglas L. Theobald
 *
 *  QuatCharPoly.c is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published
 *  by the Free Software Foundation; either version 2 of the License,
 *  or (at your option) any later version.
 *
 *  QuatCharPoly.c is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with theseus.c in the file 'COPYING'; if not, write to the:
 *
 *  Free Software Foundation, Inc.,
 *  59 Temple Place, Suite 330,
 *  Boston, MA  02111-1307  USA
 *
 *  Source:         started anew.
 *
 *  Change History:
 *    5/5/05   6:21 PM  Started source
 *    7/18/05 11:22 AM  Fixed errors in calls to CoordsInnerProd()
 *                      Removed undefined mysquare() call in CalcQuarticCoeffs()
 *                      (thanks Chris Pettitt!)
 *  
 ******************************************************************************/
/* 
    This method is highly derived from:

		Horn, B. K. P. (1987).
		"Closed-form solution of absolute orientation using unit quaternions."
		J Opt Soc Am A 4(4):629Ð642.

    Here the superposition problem is solved using a simple and very numerically
    stable Newton-Raphson procedure. The minimizing rotation is ignored, only
    the minimum RMSD is calculated.

    The main function is QuatCharPoly(), which takes as arguments:

    double **coords1      A 3 x N array of structure coordinates.
    double **coords2      The target structure coordinates.
    int len               The length N of the coords (# of atoms)
    double *coeff               A pointer to an array of 3 doubles, to hold
                                the last three quartic coefficients

    QuatCharPoly() returns the sum of squared deviations at the least squares
    minimum for coords1 and coords2. It *does not* return the RMSD.

    For the example provided below, the minimum least-squares RMSD for the two
    7-atom fragments should be 0.719106 A.

    NB #1: QuatCharPoly() returns the sum of squared deviations (sumdev^2)
           between the two structures in coords1 and coords2.
           RMSD = sqrt(sumdev^2 / N).

    NB #2: If you are doing a full superposition (the usual least squares way),
           you MUST center each structure first. That is, you must translate
           each structure so that its centroid is at the origin.
           You can use CenterCoords() for this.

    NB #3: Please note how I store structure coordinates in the double **coords
           arrays. They are 3xN arrays, not Nx3 arrays as is also commonly
           used (where the x, y, z axes are interleaved). The difference is 
           something like this for storage of a structure with 8 atoms:

           Nx3: xyzxyzxyzxyzxyzxyzxyzxyz
           3xN: xxxxxxxxyyyyyyyyzzzzzzzz

           The functions can be easily modified, however, to accomodate any
           data format preference. I chose this format because it is readily
           used in vectorized functions (SIMD, Altivec, MMX, SSE2, etc.).

    If you use this QCP RMSD calculation method in a publication, please
    reference:

        Douglas L. Theobald (2005)
        "Rapid calculation of RMSD using a quaternion-based characteristic
        polynomial."
        Acta Crystallographica A 61(4):478-480.
*/

/* gcc -O3 -Wall -Werror -ansi -o QuatCharPoly QuatCharPoly.c */

double CoordsInnerProd(double **coords, int len)
{
    int             i;
    double          sum, tmpx, tmpy, tmpz;

    sum = 0.0;
    for (i = 0; i < len; ++i)
    {
        tmpx = coords[i][0];
        tmpy = coords[i][1];
        tmpz = coords[i][2];
        sum += (tmpx*tmpx + tmpy*tmpy + tmpz*tmpz);
    }

    return(sum);
}

void quat_to_rot(double *q, double **R)
{
int i;
double t = 0.0;
R[0][0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3]; 
R[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3]; 
R[2][2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3]; 
R[0][1] = 2.0*(q[1]*q[2] + q[0]*q[3]); 
R[1][0] = 2.0*(q[1]*q[2] - q[0]*q[3]); 
R[0][2] = 2.0*(q[1]*q[3] - q[0]*q[2]); 
R[2][0] = 2.0*(q[1]*q[3] + q[0]*q[2]); 
R[1][2] = 2.0*(q[2]*q[3] + q[0]*q[1]); 
R[2][1] = 2.0*(q[2]*q[3] - q[0]*q[1]); 
}

void CalcQuarticCoeffs(double **coords1, double **coords2, int len, double *coeff)
{
    double          Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
    double          Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
                    SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
                    SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
                    SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
    double          x1, x2, y1, y2, z1, z2;
    int             i;

    Sxx = Sxy = Sxz = Syx = Syy = Syz = Szx = Szy = Szz = 0.0;
    for (i = 0; i < len; ++i)
    {
        x1 = coords1[i][0];
        y1 = coords1[i][1];
        z1 = coords1[i][2];
        x2 = coords2[i][0];
        y2 = coords2[i][1];
        z2 = coords2[i][2];
   
        Sxx += (x1 * x2);
        Sxy += (x1 * y2);
        Sxz += (x1 * z2);

        Syx += (y1 * x2);
        Syy += (y1 * y2);
        Syz += (y1 * z2);

        Szx += (z1 * x2);
        Szy += (z1 * y2);
        Szz += (z1 * z2);  
    }

    Sxx2 = Sxx * Sxx;
    Syy2 = Syy * Syy;
    Szz2 = Szz * Szz;

    Sxy2 = Sxy * Sxy;
    Syz2 = Syz * Syz;
    Sxz2 = Sxz * Sxz;

    Syx2 = Syx * Syx;
    Szy2 = Szy * Szy;
    Szx2 = Szx * Szx;

    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

    coeff[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
    coeff[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

    SxzpSzx = Sxz+Szx;
    SyzpSzy = Syz+Szy;
    SxypSyx = Sxy+Syx;
    SyzmSzy = Syz-Szy;
    SxzmSzx = Sxz-Szx;
    SxymSyx = Sxy-Syx;
    SxxpSyy = Sxx+Syy;
    SxxmSyy = Sxx-Syy;
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

    coeff[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
             + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
             + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
             + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));
}

void CalcQuarticCoeffs_M(double **coords1, double **coords2, int len, double *coeff, double **M)
{
    double          Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
    double          Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
                    SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
                    SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
                    SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
    double          x1, x2, y1, y2, z1, z2;
    int             i;
    Sxx = Sxy = Sxz = Syx = Syy = Syz = Szx = Szy = Szz = 0.0;
    for (i = 0; i < len; ++i)
    {
        x1 = coords1[i][0];
        y1 = coords1[i][1];
        z1 = coords1[i][2];
        x2 = coords2[i][0];
        y2 = coords2[i][1];
        z2 = coords2[i][2];
   
        Sxx += (x1 * x2);
        Sxy += (x1 * y2);
        Sxz += (x1 * z2);

        Syx += (y1 * x2);
        Syy += (y1 * y2);
        Syz += (y1 * z2);

        Szx += (z1 * x2);
        Szy += (z1 * y2);
        Szz += (z1 * z2);  
    }
M[0][0] = Sxx + Syy + Szz;
M[1][1] = Sxx - Syy - Szz;
M[2][2] = -Sxx + Syy - Szz;
M[3][3] = -Sxx - Syy + Szz;
M[0][1] = M[1][0] = Syz - Szy;
M[0][2] = M[2][0] = Szx - Sxz;
M[0][3] = M[3][0] = Sxy - Syx;
M[1][2] = M[2][1] = Sxy + Syx;
M[1][3] = M[3][1] = Szx + Sxz;
M[2][3] = M[3][2] = Syz + Szy;
    Sxx2 = Sxx * Sxx;
    Syy2 = Syy * Syy;
    Szz2 = Szz * Szz;

    Sxy2 = Sxy * Sxy;
    Syz2 = Syz * Syz;
    Sxz2 = Sxz * Sxz;

    Syx2 = Syx * Syx;
    Szy2 = Szy * Szy;
    Szx2 = Szx * Szx;

    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

    coeff[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
    coeff[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

    SxzpSzx = Sxz+Szx;
    SyzpSzy = Syz+Szy;
    SxypSyx = Sxy+Syx;
    SyzmSzy = Syz-Szy;
    SxzmSzx = Sxz-Szx;
    SxymSyx = Sxy-Syx;
    SxxpSyy = Sxx+Syy;
    SxxmSyy = Sxx-Syy;
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

    coeff[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
             + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
             + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
             + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));
}

/* Newton-Raphson root finding */
double QCProot(double *coeff, double guess, double delta)
{
    int             i;
    double          oldg;
    for (i = 0; i < 50; ++i)
    {
        oldg = guess;
        /* guess -= (eval_horn_quart(coeff, guess) / eval_horn_quart_deriv(coeff, guess)); */
        guess -= eval_horn_NR_corrxn(coeff, guess);
/* Qui ho modificato: se l'rmsd e' prossimo a zero accetto*/
        if ( (fabs(guess - oldg) < fabs(delta*guess)) || fabs(guess+oldg) < 10*delta)
            return(guess);

    }

    fprintf(stderr,
            "\n\n ERROR21: Newton-Raphson root-finding in \'QCProot()\' did not converge \n");
    fprintf(stderr,
            "last rmsd in iteration: %f %f %f\n", oldg ,guess, delta);
    exit(EXIT_FAILURE);
}


/* Evaluates the Newton-Raphson correction for the Horn quartic.
   only 11 FLOPs */
double eval_horn_NR_corrxn(double *c, double x)
{
    double x2 = x*x;
    double b = (x2 + c[2])*x;
    double a = b + c[1];

    return((a*x + c[0])/(2.0*x2*x + b + a));
}


/* Evaluates the Horn quartic for coefficients c and given x. */
double eval_horn_quart(double *c, double x)
{
    return(((x*x + c[2])*x + c[1])*x + c[0]);
}


/* Evaluates the derivative of the Horn quartic for
   coefficients c and given x. */
double eval_horn_quart_deriv(double *c, double x)
{
    return(2.0*(2.0*x*x + c[2])*x + c[1]);
}

/* parse a string with chain, residue number and atom names */
void parse(struct Cra *cra, int *ncra, char *s)
{
int i,j,k,l,kk,in;
int ntok,r[140], tmpi, nc, nr, na;
char p,punct[50];
char c[256],at[128][5],tok[128][256];

j = k = 0;
while(k == 0 && s[j] == ' ') j++; 
for(; j< (int) strlen(s); j++,k++) 
{
if(s[j] == ' ') 
{
s[k++]=s[j];
while(s[++j] == ' ') ;
}
s[k] = s[j];
}
s[k] = '\0';
while(s[k-1] == ' ' || s[k-1] == '\0') s[--k] = '\0';
kk = 0;
for(j=0, i=0, k=0,l=0; j< (int) strlen(s); j++)
{
if(s[j] != ' ' && s[j] != ':' && s[j] != '-' && s[j] != ',' )
tok[i][k++] = s[j];
else
{
tok[i][k] = '\0';
i++;
k=0;
punct[l++] = s[j];
}
}
tok[i][k] = '\0';
punct[l] = '\0';
ntok = i+1;
k = 0;
for(i = 0; i < ntok - 1;) 
{
l = -1;
strcpy(c,tok[i++]);
nc=strlen(c);
if(nc==0)
{
strcpy(c," ");
nc = 1;
}
if(punct[i-1] != ':') 
 {
 printf("I expected colon here... exiting...\n"); 
 exit(0);
 }
do
{
r[++l] = atoi(tok[i++]);
if(punct[i-1] == '-') 
 {
 tmpi = atoi(tok[i++]);
 while(r[l] < tmpi)
 {
 r[l+1] = r[l] + 1;
 l++;
 }
 }
}
while(punct[i-1] != ':');

nr = l+1;
l = 0;
while(punct[i-1] != ' ' && i < ntok)
strcpy(at[l++],tok[i++]);

na = l;
for(j = 0; j<nc; j++)
for(k = 0 ; k<nr; k++)
for(l = 0 ; l<na; l++)
{
cra[kk].c = c[j];
cra[kk].r = r[k];
strcpy(cra[kk].a,at[l]);
kk++;
}
}
*ncra = kk;
} 

/* count how many threads are available for openMP */
int omp_thread_count() {
    int n = 0;
    #pragma omp parallel reduction(+:n)
    n += 1;
    return n;
}


