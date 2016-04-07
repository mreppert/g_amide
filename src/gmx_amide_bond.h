#ifndef GMX_AMIDE_BOND_H
#define GMX_AMIDE_BOND_H

#include "statutil.h"
#include "copyrite.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "futil.h"
#include "tpxio.h"
#include "physics.h"
#include "macros.h"
#include "gmx_fatal.h"
#include "index.h"
#include "smalloc.h"
#include "vec.h"
#include "xvgr.h"
#include "gstat.h"
#include "matio.h"
#include "string2.h"
#include "pbc.h"

#define MAX_PATH_LENGTH (50)
typedef struct {
        int length;
        char Path[MAX_PATH_LENGTH][6];
} t_atom_path;

typedef struct {
        char resname1[4],resname2[4];		// Residue names ("PRO", "GLY", etc.)
        int resind1,resind2;			// The residue indces in top.atoms.resinfo[]
	int resnum1,resnum2;    	      	// The residue number (e.g. from pdb)
	int isPro;				// Does the amide bond involve an N-CD bond?
	int nexcluded;			// The number of atoms which will be excluded from electrostatic calculations
	int *Excluded;			// Indices for excluded atoms
	int N, H, C, O;
	int CAN, CAC;
	int phiN[4], psiN[4];
	int phiC[4], psiC[4];
	int dih_set_N;
	int dih_set_C;
	int **MapAtoms;
	int nsites;
	int *natoms;
	int **CoupAtoms;
	int ncoupsites; 
	int *ncoupatoms;
} t_protbond;

int set_dihedral_indices_N( t_pbc*, t_topology, rvec*, t_protbond* );
int set_dihedral_indices_C( t_pbc*, t_topology, rvec*, t_protbond* );
int find_covalent( t_pbc*, t_topology, rvec*, int, char*, double);
int find_all_covalent( t_pbc*, t_topology, rvec*, int, int, int** );
int print_bond(t_protbond, t_topology);
int trace_path( t_pbc*, t_topology, rvec*, int, t_atom_path);
int trace_open_path( t_pbc*, t_topology, rvec*, int, t_atom_path, int*, int*);
int find_bonds(t_pbc*, t_topology, rvec*, t_protbond**);
int initialize_bond(t_pbc*, t_topology, rvec*, t_protbond*, int, int, int, int, int);
int free_bond(t_protbond pb);

#endif



