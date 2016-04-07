/* gmx_amide_bond.c contains subroutines for operations involving the
 * t_protbond structures used for storing bond inforation by g_amide.  
 *
 * Written by Mike Reppert (2015)
 */


#include "gmx_amide_bond.h"

// Cut-off lengths for covalent bonds
const double rcutDEFAULT = 0.200;
const double rcutOC = 0.2000; 
const double rcutCN = 0.2000; 
const double rcutNH = 0.2000; 
const double rcutNCD = 0.2000; 
const double rcutNCA = 0.2000;
const double rcutCAC = 0.2000;  
const double rcutCH = 0.2000; 
const double rH = 0.05;
const double rX = 0.1;

// Chech whether two atoms are bonded to one another. We distinguish two atomic 
// radii: rH for hydrogen and rX for all others. 
int is_bonded( t_pbc *p_pbc, t_topology top, rvec *x, int at1, int at2) {
	int num1, num2;
	double r1, r2;
	double d;
	rvec dx;
	
	num1 = top.atoms.atom[at1].atomnumber;
	num2 = top.atoms.atom[at2].atomnumber;
	if(num1<3) r1 = rH;
	else r1 = rX;
	if(num2<3) r2 = rH;
	else r2 = rX;

	printf("%d: (%6.5f, %6.5f)\n", num1, r1, top.atomtypes.radius[at1]);
	printf("%d: (%6.5f, %6.5f)\n", num2, r2, top.atomtypes.radius[at2]);

	pbc_dx(p_pbc, x[at1], x[at2], dx);
	d = norm(dx);
	if(d<(r1+r2)) return 1;
	else return 0;
}

// Find covalent bonds from atom at to any other atom with atomname name, using the cutoff
// distance rcut for the bond length. 
int find_covalent( t_pbc *p_pbc, t_topology top, rvec *x, int at, char* name, double rcut ) {
	int at2;
	rvec dx;
	double d;
	for(at2=0; at2<top.atoms.nr; at2++) {
		if(at2!=at) {
			pbc_dx(p_pbc, x[at], x[at2], dx);
			d = norm(dx);
			if(d<rcut) {
				if(strncmp((*(top.atoms.atomname[at2])), name, 6)==0) return at2;
			}
		}
	}
	return -1;
}

// Find covalent bonds from atom at to any other atom with atomic number num using the cutoff
// distance rcut for the bond length. One specific atom may be excluded using the exclude flag 
// (this is used in tracing open ended paths, preventing one from re-doubling along the same path). 
int find_covalent_by_num( t_pbc *p_pbc, t_topology top, rvec *x, int at, int num, double rcut, int exclude) {
	int at2;
	rvec dx;
	double d;
	for(at2=0; at2<top.atoms.nr; at2++) {
		if(at2!=at) {
			pbc_dx(p_pbc, x[at], x[at2], dx);
			d = norm(dx);
			if(d<rcut) {
				if( (top.atoms.atom[at2].atomnumber==num) && (at2!=exclude) ) return at2;
			}
		}
	}
	return -1;
}

int set_dihedral_indices_N( t_pbc *p_pbc, t_topology top, rvec *x, t_protbond *p_pb ) {
	int at;
	int success_N = 1;
	t_protbond pb = *p_pb;

	// pb.phiN[0-3] and pb.psiN[0-3] are the atom indices for the four atoms defining the 
	// N-terminal dihedral angles for the protein bond pb, i.e. the dihedral angles
	// for the alpha carbon on the N-terminal side of the bond. 
	pb.psiN[3] = pb.N;
	pb.psiN[2] = pb.C;
	pb.phiN[3] = pb.C;
	at = find_covalent(p_pbc, top, x, pb.C, "CA", rcutDEFAULT);
	if(at<0) success_N = 0;
	else {
		pb.psiN[1] = at;	// N-terminal C-alpha
		pb.phiN[2] = at;
		at = find_covalent(p_pbc, top, x, at, "N", rcutDEFAULT);
		if(at<0) success_N = 0;
		else {
			pb.psiN[0] = at;
			pb.phiN[1] = at;
			at = find_covalent(p_pbc, top, x, at, "C", rcutDEFAULT);
			if(at<0) success_N = 0;
			else pb.phiN[0] = at;
		}
	}
	*p_pb = pb;
	// Note that we will fail for N-terminal residues. 
	return success_N;
}

// Identify atoms defining C-terminal phi and psi angles. 
int set_dihedral_indices_C( t_pbc *p_pbc, t_topology top, rvec *x, t_protbond *p_pb ) {
	int at;
	int success_C = 1;
	t_protbond pb = *p_pb;

	// pb.phiC[0-3] and pb.psiC[0-3] are the atom indices for the four atoms defining the 
	// C-terminal dihedral angles for the protein bond pb, i.e. the dihedral angles
	// for the alpha carbon on the C-terminal side of the bond. 
	pb.phiC[0] = pb.C;
	pb.phiC[1] = pb.N;
	pb.psiC[0] = pb.N;
	at = find_covalent(p_pbc, top, x, pb.N, "CA", rcutDEFAULT);
	if(at<0) success_C = 0;
	else {
		pb.phiC[2] = at;
		pb.psiC[1] = at;
		at = find_covalent(p_pbc, top, x, at, "C", rcutDEFAULT);
		if(at<0) success_C = 0;
		else {
			pb.phiC[3] = at;
			pb.psiC[2] = at;
			at = find_covalent(p_pbc, top, x, at, "N", rcutDEFAULT);
			if(at<0) success_C = 0;
			else pb.psiC[3] = at;
		}
	}
	*p_pb = pb;
	// Note that we will fail to set indices for C-terminal residues. 
	return success_C;
}

// Trace a specified atom path from start to end. Return the index of the terminating atom
// in top.atoms.atom[] array. 
int trace_path( t_pbc *p_pbc, t_topology top, rvec *x, int at0, t_atom_path aPath ) {
	int i, at;
	// First check whether the start atom is what we believe it to be. 
	if(strncmp(aPath.Path[0], (*(top.atoms.atomname[at0])), 6)!=0) {
		if( !((strncmp(aPath.Path[0], "H",1)==0) && (strncmp((*(top.atoms.atomname[at0])), "CD", 2)==0)) ) return -1;
	} 
	at = at0;
	// Now step through path entries until we come to the end. 
	for(i=1; i<aPath.length; i++) {
		at = find_covalent(p_pbc, top, x, at, aPath.Path[i], rcutDEFAULT);
		if(at<0) break;
	}
	return at;
}

// Print out all bond properties. 
int print_bond(t_protbond pb, t_topology top) {
	int i,j;
	printf("%s%d--%s%d:\n", pb.resname1, pb.resnum1, pb.resname2, pb.resnum2);
	printf("\tBond Atoms: \n");
	printf("\tCaN\tC\tO\tN\tH\tCaC\n");
	printf("\t%d\t%d\t%d\t%d\t%d\t%d\n\n", pb.CAN, pb.C, pb.O, pb.N, pb.H, pb.CAC);
	printf("\tMap Atoms:\n");
	for(i=0; i<pb.nsites; i++) {
		printf("\t(");
		for(j=0; j<pb.natoms[i]; j++) {
			if(j>0) printf(", ");
			printf("%d", pb.MapAtoms[i][j]);
		}
		printf(")  ");
	}
	printf("\n\n");
	printf("\tProline: ");
	if(pb.isPro) printf("yes\n");
	else printf("no\n");
	printf("\n");
	printf("\tExcluding %d atoms from electrostatics: \n", pb.nexcluded);
	for(i=0; i<pb.nexcluded; i++) {
		int atnum = pb.Excluded[i];
		int resndx = top.atoms.atom[atnum].resind;
		int resnum = top.atoms.resinfo[resndx].nr;
		printf("\t%s %d:\t%s (%d)\n", (*(top.atoms.resinfo[resndx].name)), resnum, (*(top.atoms.atomname[atnum])), atnum);
	}
	printf("\n\tCoupling Sites: \n");
	for(j=0; j<pb.ncoupsites; j++) {
		int k;
		printf("\tSite %d: ", j);
		for(k=0; k<pb.ncoupatoms[j]; k++) {
			if(k>0) printf(" & ");
			printf("%d", pb.CoupAtoms[j][k]);
		}
		printf("\n");
	}

	printf("\n\n\tLocated N-terminal dihedrals? ");
	if(pb.dih_set_N) {
		printf("Yes\n");
		printf("\tPhiN:\t%d\t%d\t%d\t%d\n", pb.phiN[0], pb.phiN[1], pb.phiN[2], pb.phiN[3]);
		printf("\tPsiN:\t%d\t%d\t%d\t%d\n", pb.psiN[0], pb.psiN[1], pb.psiN[2], pb.psiN[3]);
	}
	else printf("No\n");
	printf("\n\tLocated C-terminal dihedrals? ");
	if(pb.dih_set_C) {
		printf("Yes\n");
		printf("\tPhiC:\t%d\t%d\t%d\t%d\n", pb.phiC[0], pb.phiC[1], pb.phiC[2], pb.phiC[3]);
		printf("\tPsiC:\t%d\t%d\t%d\t%d\n", pb.psiC[0], pb.psiC[1], pb.psiC[2], pb.psiC[3]);
	}
	else printf("No\n");
	printf("\n\n");
	return 1;
}


// Locate bonds from the topology and structure and store them in the p_pb array. 
int find_bonds(t_pbc *p_pbc, t_topology top, rvec *x, t_protbond** pp_pb) {
	int bonds = 0;			// How many bonds we have located. 
	int nPro = 0;
	t_protbond pb;			// A dummy protein bond to store preliminary data
	t_protbond *p_pb;		// A t_protbond array so we don't have to write *(*pp_pb+i) all the time
	t_protbond *new_p_pb; 		// A dummy pointer for re-allocation purposes
	p_pb = *pp_pb;			// p_pr now points to what will be an array of protein bonds

	// We'll pre-allocate an array for nb_guess bonds in the protein. 
	// If there are more, we'll allocate more memory. If less, we'll 
	// decrease the size of the array. 
	int nb_guess = 1000;
	int nb_allocated = nb_guess;
	p_pb = (t_protbond*) malloc( sizeof(t_protbond)*nb_guess);
	if(p_pb==NULL) {
		printf("Error allocating memory for protein bond array\n");
		return 0;
	}
	
	int N, C, O, H;
	int foundBond;
	int isPro;
	// We'll now loop through atoms, looking for N, H, C, O
	// If these four atoms are located, we consider a bond to have been formed
	// First check for O atom. 
	for(O=0; O<top.atoms.nr; O++) {
		// Located an O atom. 
		if(strncmp((*(top.atoms.atomname[O])), "O", 6)==0) {
			isPro = 0;
			foundBond = 0;
			// Check whether the O atom is bonded to a C atom. 
			C = find_covalent(p_pbc, top, x, O, "C", rcutOC);
			if(C!=-1) {
				// Check whether C is bonded to N. 
				N = find_covalent(p_pbc, top, x, C, "N", rcutCN);
				if(N!=-1) {
					// And check whether N is bonded to H. 
					// Look first for H. If not, check for CD (proline).
					H = find_covalent(p_pbc, top, x, N, "H", rcutNH);
					// If we do locate an "H", we have identified a bond, and it is not proline.
					if(H!=-1) {
						foundBond = 1;
						isPro = 0;
					} else {
						// If we don't find an "H" check whether this may be proline. 
						H = find_covalent(p_pbc, top, x, N, "CD", rcutNCD);
						if(H!=-1) {
							nPro++;
							foundBond = 1;
							isPro = nPro;
						}
					}
					if(foundBond) {
						// Only the atoms C, O, N, and H are necessary to allocate a new bond. 
						initialize_bond(p_pbc, top, x, &pb, N, H, C, O, isPro);
						bonds++;
						// Check whether we need to allocate more memory in our protein bond array. 
						if(bonds==nb_allocated+1) {
							new_p_pb = (t_protbond*) realloc(p_pb, (nb_allocated+nb_guess)*sizeof(t_protbond));
							if(new_p_pb!=NULL) {
								nb_allocated += nb_guess;
								p_pb = new_p_pb;
								printf("Allocated more memory for %d bonds\n", nb_allocated);
							} else {
								printf("Error re-allocating memory for protein bond array\n");
								// Free everything and return 0.
								int i;
								for(i=0; i<bonds-1; i++) {
									free_bond(p_pb[i]);
								}
								free(p_pb);
								return 0;
							}
						}
						// Add new bond to the array. 
						p_pb[bonds-1] = pb;
					}
				}
			}
		}
	}
	// Finally, trim any extra allocated memory for the residue array and return
	new_p_pb = (t_protbond*) realloc(p_pb, (bonds)*sizeof(t_protbond));
	if(new_p_pb!=NULL) {
		p_pb = new_p_pb;
		nb_allocated = bonds;
	} else {
		if(bonds==0) return 0;	// If bonds==0, realloc has just freed the array
		printf("Error trimming extra memory in protein bond array\n");
		int i; 
		for(i=0; i<bonds; i++) {
			free_bond(p_pb[i]);
		}
		free(p_pb);
		return 0;
	}

	*pp_pb = p_pb;
	return bonds;
}

// Initialize new t_protbond variable. Set N, H, C, O, and isPro values. All others go to default. 
int initialize_bond(t_pbc *p_pbc, t_topology top, rvec *x, t_protbond *p_pb, int N, int H, int C, int O, int isPro ) {
	// res1 is the residue which contributes the CO
	// res2 is the residue which contributes the NH
	int res1 = top.atoms.atom[C].resind;
	int res2 = top.atoms.atom[N].resind;
	t_protbond pb;
	
	pb.N = N;
	pb.H = H;
	pb.C = C;
	pb.O = O;
	pb.isPro = isPro;
	pb.resind1 = res1;
	pb.resind2 = res2;
	pb.resnum1 = top.atoms.resinfo[res1].nr;
	pb.resnum2 = top.atoms.resinfo[res2].nr;	
	
	strncpy(pb.resname1, (*(top.atoms.resinfo[res1].name)), 3);
	pb.resname1[3] = '\0';  // Null character to make the char* array a true string
	strncpy(pb.resname2, (*(top.atoms.resinfo[res2].name)), 3);
	pb.resname2[3] = '\0';  // Null character to make the char* array a true string

	pb.CAN = find_covalent_by_num(p_pbc, top, x, pb.C, 6, rcutCAC, -1);
	pb.CAC = find_covalent_by_num(p_pbc, top, x, pb.N, 6, rcutNCA, pb.C);

	pb.dih_set_N = set_dihedral_indices_N( p_pbc, top, x, &pb );
	pb.dih_set_C = set_dihedral_indices_C( p_pbc, top, x, &pb );

	pb.nexcluded = 0;
	pb.nsites = 0;
	
	*p_pb = pb;
	return 1;
}


// Deallocate a protein bond. 
int free_bond(t_protbond pb) {
	int i;
	for(i=0; i<pb.nsites; i++) {
		if(pb.natoms[i]>0) free(pb.MapAtoms[i]);
		pb.natoms[i] = 0;
	}
	if(pb.nsites>0) {
		free(pb.natoms);
		free(pb.MapAtoms);
		pb.nsites = 0;
	}
	if(pb.ncoupsites>0) {
		free(pb.ncoupatoms);
		free(pb.CoupAtoms);
		pb.ncoupsites = 0;
	}
	if(pb.nexcluded>0) free(pb.Excluded);
	return 1;
}

