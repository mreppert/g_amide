/* gmx_amide_map.c contains subroutines used for reading and implementing 
 * the t_amide_map structures used by g_amide. 
 *
 * Written by Mike Reppert (2015)
 */


#include "gmx_amide_map.h"
#include "mem_helper.h"

/* MAP FILE FORMAT
 * Map files may contain the following entries: 
 * 
 * 	(1) Zero-field frequency value (usually ~ 1700 cm-1 for amide I):
 * 		FREQ: 1707.0
 *
 * 	(2) Proline-shift value (usually ~ -35 cm-1). The value may be much 
 * 	larger for maps which reference CD atom electrostatics:
 * 		PROSHIFT: -35.0
 *
 * 	(3) Zero-field dipole moment. X, Y, and Z components. 
 *
 * 		DIP: 0.242141 -0.085961 0.0
 *
 * 	(4) Map site specification. Four fundamental map sites are recognized: 
 * 	N, H, C, and O. Midpoints between sites X and Y can be specified as 
 * 	X & Y. Neighboring atoms can be specified by writing out a hyphen-
 * 	delimited path away from the fundamental site, e.g. the C-terminal 
 * 	alpha carbon is denoted C-CA. Any number of map sites can be specified, 
 * 	but the number of map sites MUST be specified following the SITES entry
 * 	and must match the number of entries in the SHIFT, DIPX, DIPY, and DIPZ 
 * 	entries. NB: the SITES entry must precede these entries in the text file. 
 * 	A few examples follow. 
 *
 * 	-- Two sites, C and O: 
 * 		SITES: 2
 * 		C
 * 		O
 *
 * 	-- One site, midpoint between C and O atoms: 
 * 		SITES: 1
 * 		C & O
 *
 * 	-- Four sites, O, C, N, and C-terminal CA:
 * 		SITES: 4
 * 		O
 * 		C
 * 		N
 * 		C-CA
 *
 * 	(5) Electrostatic-dependent frequency shifts. Each line of this entry 
 * 	specifies a set of electrostatic-shift components corresponding to the 
 * 	map sites specified earlier in the file by the SITES entry. The ten 
 * 	components in each SHIFT row represent (in order) phi, Ex, Ey, Ez, Gxx, 
 * 	Gxy, Gxz, Gyy, Gyz, and Gzz, where phi represents electrical potential,
 * 	Ei the i^th component of the electric field, and Gij the (i,j) component of 
 * 	the gradient of the electric field (i.e. second derivative of the potential). 
 * 	For example, if we specify a map which relies solely on the x-component 
 * 	of the electric field evaluated at the oxygen atom, our SITES and SHIFT
 * 	entries would have the form
 *
 * 		SITES: 1
 * 		O
 *
 * 		SHIFT: 
 * 		0.0, xxx, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
 *
 * 	(6) Electrostatic-dependent dipole shifts. Each line of this entry 
 * 	specifies a set of electrostatic dipole shift components. The format
 * 	is the same as for the SHIFT entry. X, Y, and Z components are independent,
 * 	and one can specify any one with or without the others. 
 *
 * 		DIPX: 
 * 		0.0, xxx, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
 *
 * 		DIPY: 
 * 		0.0, xxx, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
 *
 * 		DIPZ: 
 * 		0.0, xxx, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
 *
 *
 * 	(7) Coupling model specification. Options include tcc (transition
 * 	charge coupling), pdc (point-dipole coupling), edc (extended dipole 
 * 	coupling). 
 *
 * 		COUPLING: tcc
 *
 *	(8) Nearest-neighbor frequency shifts and coupling models. The NNFSN, 
 *	NNFSC, and NNC entries specify nearest-neighbor dihedral maps based
 *	on the phi and psi angles between amide bonds. NNFSN and NNFSC specify
 *	frequency shifts due to the N- and C-terminal amide groups, while the 
 *	NNC entry specifies the coupling between two amide bonds. 
 *
 *	The format for all three entries is the same. The flag (NNFSN, NNFSC,
 *	or NNC) should be followed by a pair of integers specifying the size of 
 *	the matrix to follow. The matrix should then be printed column-delimited,
 *	one row per line. This matrix specifies the corresponding value as a 
 *	function of phi (row) and psi (column) dihedral angles. The first entry
 *	will be interpreted as corresponding to (phi,psi) = (-180,-180), and 
 *	the grid will be assumed to extend from -180 to +180 for both angles. 
 *	The grid spacing will thus be determined by the size of the matrix.
 *	NB: the first and last entries of every line and every column should
 *	be identical since e.g. (-180,xxx) is equivalent to (+180,xxx), and 
 *	both values should be reported in the matrix. For example, the entry
 *
 *		NNC: 13 13
 *
 *	specifies that a 13x13 matrix will follow, specifying a dihedral 
 *	coupling map in 360/(13-1) = 30 degree increments. In spectral
 *	calculations, these values will be linearly interpolated to the 
 *	real (phi,psi) values from the protein structure. 
 *
 *	For example, suppose we have a peptide with the sequence 
 * 	AAA-BBB-CCC-DDD, i.e. four residues linked by three amide bonds. 
 *	The frequency of the AAA-BBB bond will be determined by the 
 *	electrostatic maps and a nearest-neighbor frequency shift from its
 *	C-terminal neighbor CCC; this shift will be specified by the NNFSC 
 *	matrix and the (phi,psi) dihedral angles at the BBB C-alpha atom. 
 *	The coupling between bonds AAA-BBB and BBB-CCC will likewise be 
 *	determined by the same set of dihedral angles. No N-terminal frequency
 *	shift will be used since the bond lacks an N-terminal neighbor. 
 *
 *	The frequency of the CCC-DDD bond will conversely contain a 
 *	component due to the NNFSN map and the dihedral angles of the
 *	CCC C-alpha carbon; these same angles will determine the 
 *	CCC-DDD to BBB-CCC coupling constant. In this case, no C-terminal
 *	frequency shift will be in play. 
 *
 *	Finally, the frequency of the BBB-CCC bond will contain both a 
 *	C-terminal contribution (determined by the CCC C-alpha angles)
 *	and an N-terminal contribution (from the BBB-C-alpha angles). 
 *
 *
 * 	(9) Excluded atoms. This will generally be the longest entry
 * 	in a map file, and is intended to provide flexibility in which 
 * 	atoms are excluded from electrostatic calculations. The entry
 * 	should ALWAYS contain the map sites, along with their bonded 
 * 	neighbors. Additional excluded atoms depend on the map. The number
 * 	of excluded atoms should be specified following the EXCLUDE header.
 *
 * 	For example, to exclude the C, O, N, and H atoms, along with 
 * 	their neighboring CA atoms, the mapfile entry should be
 *
 * 		EXCLUDE: 6
 * 		C
 * 		O
 * 		N
 * 		H
 * 		N-CA
 * 		C-CA
 *
 * 	The last two entries specify the C- and N-terminal alpha
 * 	carbons, respectively (note that N-CA is the C-terminal 
 * 	carbon since the NH group is on the C-terminal end of the 
 * 	amide bond. 
 * 		
 * 	Note that the CD atom in proline will automatically be 
 * 	interpreted as an H atom for map purposes, since it corrresponds
 * 	in position to the H atom in normal amide bonds. Thus, for a map
 * 	parametrizing the N atom, it is not necessary to separately specify 
 * 	the CD atom. In contrast, if the H atom itself is a map site, 
 * 	in the case of proline-containing bonds, the HD atoms will by 
 * 	default be included in the electrostatic calculations, unless
 * 	separate CD-HD1 and CD-HD2 entries are explicitly specified in
 * 	the EXCLUDE section. 
 *
 * */


const int maxchar = 512;
const char* path_delim = "&";
const char* bond_delim = "-";
const char* mat_delim = " ,";
const int nelec = 10;

int count_tokens(const char* line, const char* delim) {
	int count = 0;
	char line_copy[maxchar];
	strcpy(line_copy, line);
	char* p = strtok(line_copy, delim);
	while(p!=NULL) {
		p = strtok(NULL, delim);
		count++;
	}
	return count;
}

int mystrtok( char* remainder, char* token, const char* delim_str) {
	char* p;
	char buffer[maxchar];
	token[0] = '\0';
	p = strstr(remainder, delim_str);
	if(p!=NULL) {
		*p = '\0';
		strcpy(token, remainder);
		strcpy(buffer, p+strlen(delim_str));
		strcpy(remainder, buffer);
	} else {
		strcpy(token, remainder);
		remainder[0] = '\0';
	}
	return 1;
}

int set_atom_path(t_atom_path *p_ap, const char* str) {
	t_atom_path ap = *p_ap;
	char str_copy[maxchar];
	strcpy(str_copy, str);
	char token[maxchar], remainder[maxchar];
	strcpy(remainder, str);
	ap.length = 0;
	while(strlen(remainder)>0) {
		mystrtok(remainder, token, bond_delim);
		if(strlen(token)>5) {
			printf("Error: Atom identifier %s is too long (5 character is pdb maximum).\n", token);
			ap.length = 0;
			*p_ap = ap;
			return 0;
		}
		sscanf(token, "%s", ap.Path[ap.length]);
		ap.length++;
		if(ap.length==MAX_PATH_LENGTH) { 
			printf("Error setting atom path. Path exceeds maximum length (%d bonds): %s\n", MAX_PATH_LENGTH, str);
			return 0;
		}
	}
	*p_ap = ap;
	return 1;
}

int read_sites(FILE *fp, t_amide_map *p_map) {
	int i;
	int success = 0;	// Number of map sites successfully allocated
	char line[maxchar];
	t_amide_map map = *p_map;
	
	map.MapSites = (t_mapsite*) malloc(map.nsites*sizeof(t_mapsite));
	if(map.MapSites==NULL) {
		map.nsites = 0;
		success = 0;
	}
	for(i=0; i<map.nsites; i++) {
		if(fgets(line, maxchar, fp)!=NULL) {
			int natoms = count_tokens(line, path_delim);
			int error = 0;
			if(natoms>0) {
				map.MapSites[i].AtomPaths = (t_atom_path*) malloc(natoms*sizeof(t_atom_path));
				if(map.MapSites[i].AtomPaths!=NULL) {
					map.MapSites[i].natoms = natoms;
					char token[maxchar], remainder[maxchar];
					int path = 0;
					strcpy(remainder, line);
					while(strlen(remainder)!=0) {
						mystrtok(remainder, token, path_delim);
						set_atom_path(&map.MapSites[i].AtomPaths[path], token);
						path++;
					}
				} else error = 1;
			} else error = 1;
			if(error) {
				i = map.nsites;
				break;
			}
			success++;
		} else return 0;
	}
	if(success<map.nsites) {
		for(i=0; i<success; i++) {
			free(map.MapSites[i].AtomPaths);
		}
		free(map.MapSites);
		map.nsites = 0;
	}
	*p_map = map;
	return 1;
}

int read_coupsites(FILE *fp, t_amide_map *p_map) {
	int i;
	int success = 0;	// Number of map sites successfully allocated
	char line[maxchar];
	t_amide_map map = *p_map;
	
	map.CoupSites = (t_coupsite*) malloc(map.ncoupsites*sizeof(t_coupsite));
	if(map.CoupSites==NULL) {
		map.ncoupsites = 0;
		success = 0;
	}
	for(i=0; i<map.ncoupsites; i++) {
		if(fgets(line, maxchar, fp)!=NULL) {
			int natoms = count_tokens(line, path_delim);
			int error = 0;
			if(natoms>0) {
				map.CoupSites[i].AtomPaths = (t_atom_path*) malloc(natoms*sizeof(t_atom_path));
				if(map.CoupSites[i].AtomPaths!=NULL) {
					map.CoupSites[i].natoms = natoms;
					char token[maxchar], remainder[maxchar];
					int path = 0;
					strcpy(remainder, line);
					while(strlen(remainder)!=0) {
						mystrtok(remainder, token, path_delim);
						set_atom_path(&map.CoupSites[i].AtomPaths[path], token);
						path++;
					}
				} else error = 1;
			} else error = 1;
			if(error) {
				i = map.ncoupsites;
				break;
			}
			success++;
		} else return 0;
	}
	if(success<map.ncoupsites) {
		for(i=0; i<success; i++) {
			free(map.CoupSites[i].AtomPaths);
		}
		free(map.CoupSites);
		map.ncoupsites = 0;
	}
	*p_map = map;
	return 1;
}


int read_site_coeff( FILE* fp, t_amide_map *p_map, char* str ) {
	t_amide_map map  = *p_map;
	char line[maxchar];
	real Array[10];
	int i,j;
	if(map.nsites==0) {
		printf("Error parsing site coefficients: no map sites to parametrize.\n");
		printf("Please check that the SITES flag appears first in the map file.\n");
		return 0;
	}
	for(i=0; i<map.nsites; i++) {
		if(fgets(line, maxchar, fp)!=NULL) {
			int j = 0; 
			char* p;
			p = strtok(line, mat_delim);
			do {
				if(sscanf(p, "%f", &Array[j])!=1) break;
				p = strtok(NULL, mat_delim);
				j++;
			} while ( (p!=NULL) && (j<nelec) );
			if( j!=nelec ) {
				printf("Error parsing shift parameters for site %d. Check line %s", i, line);
				return 0;
			}
		}
		if(strncmp(str, "SHIFT", 5)==0) {
			for(j=0; j<nelec; j++) map.MapSites[i].shift[j] = Array[j];
		} else if(strncmp(str, "DIPX", 4)==0) {
			for(j=0; j<nelec; j++) map.MapSites[i].dipx[j] = Array[j];
		} else if(strncmp(str, "DIPY", 4)==0) {
			for(j=0; j<nelec; j++) map.MapSites[i].dipy[j] = Array[j];
		} else if(strncmp(str, "DIPZ", 4)==0) {
			for(j=0; j<nelec; j++) map.MapSites[i].dipz[j] = Array[j];
		}
	}
	*p_map = map;
	return 1;
}


// Read amide map from text file. 
int read_amide_map( char* fnm, t_amide_map *p_map, int verbose ) {
	t_amide_map map;
	FILE* fp = fopen(fnm, "r");
	if(fp==NULL) {
		printf("Error opening map file. Please check input.\n");
		return 0;
	}

	const int nelec = 10;
	char line[maxchar];

	map.name[0] = '\0';
	map.nsites = 0;
	map.ncoupsites = 0;
	map.npaths = 0;
	map.freq = 0.0;
	map.proshift = 0.0;
	map.dimNNFSN[0] = 0;
	map.dimNNFSC[0] = 0;
	map.dimNNC[0] = 0;
	map.dimNNFSN[1] = 0;
	map.dimNNFSC[1] = 0;
	map.dimNNC[1] = 0;
	map.dipset = 0;
	map.tcset = 0;
	int i,j,k;
	for(i=0; i<nelec; i++) {
		map.elec_used[i] = 0;
	}

	while(fgets(line, maxchar, fp)!=NULL) {
		if(!strncmp(line, "SITES:", 6)) {
			if(sscanf(line, "%*s%d", &map.nsites)>0) {
				if(!read_sites(fp, &map)) {
					printf("Error reading map sites. Please check input.\n");
					return 0;
				}
				for(i=0; i<map.nsites; i++) {
					for(j=0; j<nelec; j++) {
						map.MapSites[i].shift[j] = 0.0;
						map.MapSites[i].dipx[j] = 0.0;
						map.MapSites[i].dipy[j] = 0.0;
						map.MapSites[i].dipz[j] = 0.0;
					}
				}
				if(verbose){
					 printf("Number of sites: %d\n", map.nsites);
					for(i=0; i<map.nsites; i++) {
						printf("\tSite %d: ", i);
						for(j=0; j<map.MapSites[i].natoms; j++) {
							for(k=0; k<map.MapSites[i].AtomPaths[j].length; k++) {
								printf("%s", map.MapSites[i].AtomPaths[j].Path[k]);
								if(k<map.MapSites[i].AtomPaths[j].length-1) printf("-");
							}
							if(j<map.MapSites[i].natoms-1) printf(" %s ", path_delim);
						}
						printf("\n");
					}
					printf("\n");
				}
			} else {
				printf("Error parsing number of sites from map file %s.\n", fnm);
				printf("The bad line seems to be: \n\t%s", line);
				return 0;
			}
		} else if(!strncmp(line, "FREQ:", 5)) {
			if(sscanf(line, "%*s%f", &map.freq)>0) {
				if(verbose) printf("Successfully set gas phase frequency at %6.2f cm-1\n\n", map.freq);
			} else {
				printf("Error setting gas phase frequency.\n");
				printf("The bad line seems to be: \n\t%s\n", line);
				return 0;
			}
		} else if(!strncmp(line, "PROSHIFT:", 9)) {
			if(sscanf(line, "%*s%f", &map.proshift)>0) {
				if(verbose) printf("Successfully set proline shift at %6.2f cm-1\n\n", map.proshift);
			} else {
				printf("Error setting proline shift frequency.\n");
				printf("The bad line seems to be: \n\t%s\n", line);
				return 0;
			}
		} else if(!strncmp(line, "DIP:", 4)) {
			if(sscanf(line, "%*s%f%f%f", &map.dip[0], &map.dip[1], &map.dip[2])==3) {
				if(verbose) printf("Successfully set gas phase dipole: \t%6.4f\t%6.4f\t%6.4f\n\n", map.dip[0], map.dip[1], map.dip[2]);
				map.dipset = 1;
			} else {
				printf("Error setting gas phase transition dipole moment.\n");
				printf("The bad line seems to be: \n\t%s\n", line);
				return 0;
			}
		} else if(!strncmp(line, "SHIFT:", 6)) {
			if(!read_site_coeff(fp, &map, line)) {
				printf("Error reading site freq coefficients. Please check input.\n");
				return 0;
			} else {
				if(verbose) {
					printf("Successfully read site freq coefficients: \n");
					for(i=0; i<map.nsites; i++) {
						for(j=0; j<nelec; j++) {
							printf("%6.1f\t", map.MapSites[i].shift[j]);
						}
						printf("\n");
					}
					printf("\n");
				}
			}
		} else if(!strncmp(line, "DIPX:", 5)) {
			if(!read_site_coeff(fp, &map, line)) {
				printf("Error reading site dipx coefficients. Please check input.\n");
			} else {
				if(verbose) {
					printf("Successfully read site dipx coefficients: \n");
					for(i=0; i<map.nsites; i++) {
						for(j=0; j<nelec; j++) {
							printf("%6.3f\t", map.MapSites[i].dipx[j]);
						}
						printf("\n");
					}
					printf("\n");
				}
			}
		} else if(!strncmp(line, "DIPY:", 5)) {
			if(!read_site_coeff(fp, &map, line)) {
				printf("Error reading site dipy coefficients. Please check input.\n");
			} else {
				if(verbose) {
					printf("Successfully read site dipy coefficients: \n");
					for(i=0; i<map.nsites; i++) {
						for(j=0; j<nelec; j++) {
							printf("%6.3f\t", map.MapSites[i].dipy[j]);
						}
						printf("\n");
					}
					printf("\n");
				}
			}
		} else if(!strncmp(line, "DIPZ:", 5)) {
			if(!read_site_coeff(fp, &map, line)) {
				printf("Error reading site dipz coefficients. Please check input.\n");
			} else {
				if(verbose) {
					printf("Successfully read site dipz coefficients: \n");
					for(i=0; i<map.nsites; i++) {
						for(j=0; j<nelec; j++) {
							printf("%6.3f\t", map.MapSites[i].dipz[j]);
						}
						printf("\n");
					}
					printf("\n");
				}
			}
		} else if(!strncmp(line, "NNFSN:", 6)) {
			if(sscanf(line, "%*s%d%d", &map.dimNNFSN[0], &map.dimNNFSN[1])==2) {
				if(map.dimNNFSN[0]*map.dimNNFSN[1]==0) break;
				if(!allocate_2d_array_real(&map.NNFSN, map.dimNNFSN[0], map.dimNNFSN[1])) {
					map.dimNNFSN[0] = 0;
					map.dimNNFSN[1] = 0;
					printf("Error allocating memory for NNFSN map\n");
					return 0;
				}
				for(i=0; i<map.dimNNFSN[0]; i++) {
					j=0;
					if(fgets(line, maxchar, fp)!=NULL) {
						char* p;
						p = strtok(line, mat_delim);
						do {
							if(sscanf(p, "%f", &map.NNFSN[i][j])!=1) break;
							p = strtok(NULL, mat_delim);
							j++;
						} while( (p!=NULL) && j<map.dimNNFSN[1] );
					} 
					if( j==map.dimNNFSN[1] ) {
						if(i==map.dimNNFSN[0]-1) {
							if(verbose) {
								printf("\nSuccessfully read %d-by-%d NNFSN map\n", map.dimNNFSN[0], map.dimNNFSN[1]);
								int k,l;
								for(k=0; k<map.dimNNFSN[0]; k++) {
									for(l=0; l<map.dimNNFSN[1]; l++) printf("%5.2f\t", map.NNFSN[k][l]);
									printf("\n");
								}
							}
						}
					} else {
						printf("Error reading NNFSN parameters (expected %d-by-%d matrix).\n", map.dimNNFSN[0], map.dimNNFSN[1]);
						return 0;
					}
				}
			} else {
				printf("Error reading dimension of N-terminal nearest-neighbor frequency shift map.\n");
				return 0;
			}
		} else if(!strncmp(line, "NNFSC:", 6)) {
			if(sscanf(line, "%*s%d%d", &map.dimNNFSC[0], &map.dimNNFSC[1])==2) {
				if( (map.dimNNFSC[0]*map.dimNNFSC[1])==0 ) break;
				if(!allocate_2d_array_real(&map.NNFSC, map.dimNNFSC[0], map.dimNNFSC[1])) {
					map.dimNNFSC[0] = 0;
					map.dimNNFSC[1] = 0;
					printf("Error allocating memory for NNFSC map\n");
					return 0;
				}
				
				for(i=0; i<map.dimNNFSC[0]; i++) {
					j=0;
					if(fgets(line, maxchar, fp)!=NULL) {
						char* p;
						p = strtok(line, mat_delim);
						do {
							if(sscanf(p, "%f", &map.NNFSC[i][j])!=1) break;
							p = strtok(NULL, mat_delim);
							j++;
						} while( (p!=NULL) && j<map.dimNNFSC[1] );
					} 
					if( j==map.dimNNFSC[1] ) {
						if(i==map.dimNNFSC[0]-1) {
							if(verbose) {
								printf("\nSuccessfully read %d-by-%d NNFSC map\n", map.dimNNFSC[0], map.dimNNFSC[1]);
								int k,l;
								for(k=0; k<map.dimNNFSC[0]; k++) {
									for(l=0; l<map.dimNNFSC[1]; l++) printf("%5.2f\t", map.NNFSC[k][l]);
									printf("\n");
								}
							}
						}
					} else {
						printf("Error reading NNFSC parameters (expected %d-by-%d matrix).\n", map.dimNNFSC[0], map.dimNNFSC[1]);
						return 0;
					}
				}
			} else {
				printf("Error reading dimension of C-terminal nearest-neighbor frequency shift map.\n");
				return 0;
			}
		} else if(!strncmp(line, "COUPLING:", 9)) {
			int error = 0;
			char couptype[64];
			if(sscanf(line, "%*s %s", couptype)) {
				if(!strcmp(couptype, "pdc")) {
					strcpy(map.couptype, "pdc");
				} else if(!strcmp(couptype, "tcc")) {
					strcpy(map.couptype, "tcc");
				} else error = 1;
			} else error = 1;
			if(error) {
				printf("Error reading couling model type. Please enter pdc or tcc\n");
				free_amide_map(map);
				return 0;
			}

		} else if(!strncmp(line, "COUPSITES:", 10)) {
			if(sscanf(line, "%*s%d", &map.ncoupsites)>0) {
				if(!read_coupsites(fp, &map)) {
					printf("Error reading map coupling sites. Please check input.\n");
					return 0;
				}
				for(i=0; i<map.ncoupsites; i++) {
					map.CoupSites[i].nux = 0.0;
					map.CoupSites[i].nuy = 0.0;
					map.CoupSites[i].nuz = 0.0;
					map.CoupSites[i].q = 0.0;
					map.CoupSites[i].dq = 0.0;
					map.tcset = 1;
				}
				if(verbose) {
					printf("\nLooking for %d coupling sites:", map.ncoupsites);
					for(i=0; i<map.ncoupsites; i++) {
						printf("\n\tSite %d: ", i);
						for(j=0; j<map.CoupSites[i].natoms; j++) {
							for(k=0; k<map.CoupSites[i].AtomPaths[j].length; k++) {
								printf("%s", map.CoupSites[i].AtomPaths[j].Path[k]);
								if(k<map.CoupSites[i].AtomPaths[j].length-1) printf("-");
							}
							if(j<map.CoupSites[i].natoms-1) printf(" %s ", path_delim);
						}
					}
					printf("\n");
				}
			} else {
				printf("Error parsing number of coupling sites from map file %s.\n", fnm);
				printf("The bad line seems to be: \n\t%s\n", line);
				return 0;
			}

		} else if(!strncmp(line, "NUX:", 4)) {
			if(map.ncoupsites==0) {
				printf("Error! Found NUX parameter flag with no coupling sites specified.\n");
				printf("Please ensure that COUPSITES flag precedes NUX, NUY, NUZ, Q, and DQ flags.\n");
				return 0;
			}
			for(i=0; i<map.ncoupsites; i++) {
				if( (fgets(line, maxchar, fp)==NULL) || (sscanf(line, "%f", &map.CoupSites[i].nux)!=1) ) {
					printf("Error reading NUX component for coupling site %d.\n", i);
					return 0;
				}
			}
			if(verbose) {
				printf("\nSuccessfully read all NUX components: \n");
				for(i=0; i<map.ncoupsites; i++) printf("Site %d: %6.10f\n", i, map.CoupSites[i].nux);
			}
		} else if(!strncmp(line, "NUY:", 4)) {
			if(map.ncoupsites==0) {
				printf("Error! Found NUY parameter flag with no coupling sites specified.\n");
				printf("Please ensure that COUPSITES flag precedes NUX, NUY, NUZ, Q, and DQ flags.\n");
				return 0;
			}
			for(i=0; i<map.ncoupsites; i++) {
				if( (fgets(line, maxchar, fp)==NULL) || (sscanf(line, "%f", &map.CoupSites[i].nuy)!=1) ) {
					printf("Error reading NUX component for coupling site %d.\n", i);
					return 0;
				}
			}
			if(verbose) {
				printf("\nSuccessfully read all NUY components: \n");
				for(i=0; i<map.ncoupsites; i++) printf("Site %d: %6.10f\n", i, map.CoupSites[i].nuy);
			}
		} else if(!strncmp(line, "NUZ:", 4)) {
			if(map.ncoupsites==0) {
				printf("Error! Found NUZ parameter flag with no coupling sites specified.\n");
				printf("Please ensure that COUPSITES flag precedes NUX, NUY, NUZ, Q, and DQ flags.\n");
				return 0;
			}
			for(i=0; i<map.ncoupsites; i++) {
				if( (fgets(line, maxchar, fp)==NULL) || (sscanf(line, "%f", &map.CoupSites[i].nuz)!=1) ) {
					printf("Error reading NUY component for coupling site %d.\n", i);
					return 0;
				}
			}
			if(verbose) {
				printf("\nSuccessfully read all NUZ components: \n");
				for(i=0; i<map.ncoupsites; i++) printf("Site %d: %6.10f\n", i, map.CoupSites[i].nuz);
			}
		} else if(!strncmp(line, "Q:", 2)) {
			if(map.ncoupsites==0) {
				printf("Error! Found Q parameter flag with no coupling sites specified.\n");
				printf("Please ensure that COUPSITES flag precedes NUX, NUY, NUZ, Q, and DQ flags.\n");
				return 0;
			}
			for(i=0; i<map.ncoupsites; i++) {
				if( (fgets(line, maxchar, fp)==NULL) || (sscanf(line, "%f", &map.CoupSites[i].q)!=1) ) {
					printf("Error reading Q component for coupling site %d.\n", i);
					return 0;
				}
			}
			if(verbose) {
				printf("\nSuccessfully read all Q components: \n");
				for(i=0; i<map.ncoupsites; i++) printf("Site %d: %6.10f\n", i, map.CoupSites[i].q);
			}
		} else if(!strncmp(line, "DQ:", 3)) {
			if(map.ncoupsites==0) {
				printf("Error! Found DQ parameter flag with no coupling sites specified.\n");
				printf("Please ensure that COUPSITES flag precedes NUX, NUY, NUZ, Q, and DQ flags.\n");
				return 0;
			}
			for(i=0; i<map.ncoupsites; i++) {
				if( (fgets(line, maxchar, fp)==NULL) || (sscanf(line, "%f", &map.CoupSites[i].dq)!=1) ) {
					printf("Error reading DQ component for coupling site %d.\n", i);
					return 0;
				}
			}
			if(verbose) {
				printf("\nSuccessfully read all DQ components: \n");
				for(i=0; i<map.ncoupsites; i++) printf("Site %d: %6.10f\n", i, map.CoupSites[i].dq);
			}
		} else if(!strncmp(line, "NNC:", 4)) {
			if(sscanf(line, "%*s%d%d", &map.dimNNC[0], &map.dimNNC[1])==2) {
				if(map.dimNNC[0]==0 || map.dimNNC[1]==0) break;
				map.NNC = (real**) malloc(map.dimNNC[0]*sizeof(real*));
				if(map.NNC==NULL) {
					map.dimNNC[0] = 0;
					map.dimNNC[1] = 0;
					printf("Error allocating memory for NNC map\n");
					return 0;
				}
				for(i=0; i<map.dimNNC[0]; i++) {
					map.NNC[i] = (real*) malloc(map.dimNNC[1]*sizeof(real));
					if(map.NNC[i]==NULL) {
						printf("Error allocating memory for NNC map\n");
						int ix;
						for(ix=0; ix<i; ix++) free(map.NNC[ix]);
						free(map.NNC);
						return 0;
					}
				}
				for(i=0; i<map.dimNNC[0]; i++) {
					j=0;
					if(fgets(line, maxchar, fp)!=NULL) {
						char* p;
						p = strtok(line, mat_delim);
						do {
							if(sscanf(p, "%f", &map.NNC[i][j])!=1) break;
							p = strtok(NULL, mat_delim);
							j++;
						} while( (p!=NULL) && j<map.dimNNC[1] );
					} 
					if( j==map.dimNNC[1] ) {
						if(i==map.dimNNC[0]-1) {
							if(verbose) {
								printf("\nSuccessfully read %d-by-%d NNC map\n", map.dimNNC[0], map.dimNNC[1]);
								int k,l;
								for(k=0; k<map.dimNNC[0]; k++) {
									for(l=0; l<map.dimNNC[1]; l++) printf("%5.2f\t", map.NNC[k][l]);
									printf("\n");
								}
							}
						}
					} else {
						printf("Error reading NNC parameters (expected %d-by-%d matrix).\n", map.dimNNC[0], map.dimNNC[1]);
						return 0;
					}
				}
			} else {
				printf("Error reading dimension of nearest-neighbor coupling map.\n");
				return 0;
			}
		} else if(!strncmp(line, "DNNC:", 5)) {
			if(sscanf(line, "%*s%d%d", &map.dimDNNC[0], &map.dimDNNC[1])==2) {
				if(map.dimDNNC[0]==0 || map.dimDNNC[1]==0) break;
				map.DNNC = (real**) malloc(map.dimDNNC[0]*sizeof(real*));
				if(map.DNNC==NULL) {
					map.dimDNNC[0] = 0;
					map.dimDNNC[1] = 0;
					printf("Error allocating memory for DNNC map\n");
					return 0;
				}
				for(i=0; i<map.dimDNNC[0]; i++) {
					map.DNNC[i] = (real*) malloc(map.dimDNNC[1]*sizeof(real));
					if(map.DNNC[i]==NULL) {
						printf("Error allocating memory for DNNC map\n");
						int ix;
						for(ix=0; ix<i; ix++) free(map.DNNC[ix]);
						free(map.DNNC);
						return 0;
					}
				}
				for(i=0; i<map.dimDNNC[0]; i++) {
					j=0;
					if(fgets(line, maxchar, fp)!=NULL) {
						char* p;
						p = strtok(line, mat_delim);
						do {
							if(sscanf(p, "%f", &map.DNNC[i][j])!=1) break;
							p = strtok(NULL, mat_delim);
							j++;
						} while( (p!=NULL) && j<map.dimDNNC[1] );
					} 
					if( j==map.dimDNNC[1] ) {
						if(i==map.dimDNNC[0]-1) {
							if(verbose) {
								printf("\nSuccessfully read %d-by-%d DNNC map\n", map.dimDNNC[0], map.dimDNNC[1]);
								int k,l;
								for(k=0; k<map.dimDNNC[0]; k++) {
									for(l=0; l<map.dimDNNC[1]; l++) printf("%5.2f\t", map.DNNC[k][l]);
									printf("\n");
								}
							}
						}
					} else {
						printf("Error reading DNNC parameters (expected %d-by-%d matrix).\n", map.dimDNNC[0], map.dimDNNC[1]);
						return 0;
					}
				}
			} else {
				printf("Error reading dimension of nearest-neighbor coupling map.\n");
				return 0;
			}

		} else if(!strncmp(line, "EXCLUDE:", 8)) {
			map.npaths = 0;
			if(sscanf(line, "%*s%d", &map.npaths)==1) {
				map.ExcludedAtomPaths = (t_atom_path*) malloc(map.npaths*sizeof(t_atom_path));
				if(map.ExcludedAtomPaths==NULL) {
					map.npaths = 0;
					printf("Error allocating memory for excluded atom paths.\n");
					return 0;
				}
				for(i=0; i<map.npaths; i++) {
					if(fgets(line, maxchar, fp)!=NULL) {
						if(!set_atom_path(&map.ExcludedAtomPaths[i], line)) {
							printf("Error reading excluded atom #%d from line: \n%s", i+1, line);
							return 0;
						}
					}
				}
				if(verbose) {
					printf("Successfully read %d excluded atom paths:\n", map.npaths);
					for(i=0; i<map.npaths; i++) {
						for(j=0; j<map.ExcludedAtomPaths[i].length; j++) {
							printf("%s", map.ExcludedAtomPaths[i].Path[j]);
							if(j<map.ExcludedAtomPaths[i].length-1) printf("-");
						}
						printf("\n");
					}
				}
			} else {
				printf("Error reading number of excluded atoms.\n");
				return 0;
			}
		}
	}

	double eps = 1e-10;
	for(i=0; i<nelec; i++) {
		map.elec_used[i] = 0;
		for(j=0; j<map.nsites; j++) {
			if(abs(map.MapSites[j].shift[i])>eps) map.elec_used[i] = 1;
			if(abs(map.MapSites[j].dipx[i])>eps) map.elec_used[i] = 1;
			if(abs(map.MapSites[j].dipy[i])>eps) map.elec_used[i] = 1;
		}
	}

	if( (!strcmp(map.couptype, "tcc")) && (map.ncoupsites==0) ) {
		printf("Error! Transition charge coupling requested, but no map sites specified.\n");
		free_amide_map(map);
		return 0;
	} else if( (!strcmp(map.couptype, "pdc")) && (map.ncoupsites!=0) ) {
		printf("Error! Point dipole coupling requested, but transition charge coupling sites seem to be specified. Which model do you want to use?\n");
		free_amide_map(map);
		return 0;
	}
	*p_map = map;
	fclose(fp);
	return 1;
}

int free_amide_map( t_amide_map map ) {
	int i;
	for(i=0; i<map.nsites; i++) {
		free(map.MapSites[i].AtomPaths);
		map.MapSites[i].natoms = 0;
	}
	if(map.nsites>0) free(map.MapSites);
	map.nsites = 0;
	if(map.npaths>0) free(map.ExcludedAtomPaths);
	map.npaths = 0;
	for(i=0; i<map.ncoupsites; i++) {
		free(map.CoupSites[i].AtomPaths);
		map.CoupSites[i].natoms = 0;
	}
	if(map.ncoupsites>0) free(map.CoupSites);
	map.ncoupsites = 0;
	if(map.dimNNFSN[0]>0 && map.dimNNFSN[1]>0) {
		for(i=0; i<map.dimNNFSN[0]; i++) free(map.NNFSN[i]);
		free(map.NNFSN);
	}
	if(map.dimNNFSC[0]>0 && map.dimNNFSC[1]>0) {
		for(i=0; i<map.dimNNFSC[0]; i++) free(map.NNFSC[i]);
		free(map.NNFSC);
	}
	if(map.dimNNC[0]>0 && map.dimNNC[1]>0) {
		for(i=0; i<map.dimNNC[0]; i++) free(map.NNC[i]);
		free(map.NNC);
	}
	if(map.dimDNNC[0]>0 && map.dimDNNC[1]>0) {
		for(i=0; i<map.dimDNNC[0]; i++) free(map.DNNC[i]);
		free(map.DNNC);
	}

	return 1;
}
