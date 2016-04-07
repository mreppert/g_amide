#ifndef GMX_AMIDE_CHARGE_H
#define GMX_AMIDE_CHARGE_H

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

#include "string.h"

typedef struct { 
	char resname[10];
	int nnames;	
	int natoms;
	char **names;
	char **atomnames;
	float *atomcharges;
} t_residue;

int free_residue(t_residue);
int free_resarray(t_residue*, int);
int swap_atomnames(t_topology, char **repnames[2], int, int);
int read_charge_map(char*, t_residue**, int*, char **resnames[2], int*, int);
int identify_residues(t_pbc*, t_topology, rvec*, t_residue *, int, int);

#endif
