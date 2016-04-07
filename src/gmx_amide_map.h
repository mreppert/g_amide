#ifndef GMX_AMIDE_MAP_H
#define GMX_AMIDE_MAP_H

#include "stdio.h"
#include "string.h"

#include "gmx_amide_bond.h"

#include "statutil.h"

typedef struct {
	int natoms;
	t_atom_path* AtomPaths;
	real shift[10], dipx[10], dipy[10], dipz[10];
} t_mapsite;

typedef struct {
	int natoms; 
	t_atom_path* AtomPaths;
	real q;
	real dq;
	real nux;
	real nuy;
	real nuz;
} t_coupsite;

typedef struct {
	real freq;
	real proshift;
	rvec dip;
	char name[64];
	int nsites;
	t_mapsite* MapSites;
	int dimNNFSN[2], dimNNFSC[2], dimNNC[2], dimDNNC[2];
	real **NNFSN, **NNFSC, **NNC, **DNNC;
	int npaths;
	t_atom_path* ExcludedAtomPaths;
	int elec_used[10];
	char couptype[64];
	int ncoupsites;
	t_coupsite* CoupSites;
	int dipset;
	int tcset;
} t_amide_map;


int read_amide_map( char*, t_amide_map*, int verbose );
int free_amide_map( t_amide_map );

#endif
