/* gmx_amide_charge.c contains subroutines used for manipulating atomic
 * charges and called by g_amide. 
 *
 * Written by Mike Reppert (2015)
 */


#include "gmx_amide_bond.h"
#include "gmx_amide_charge.h"

/* Charge file format. 
 *
 * A charge map file contains two types of entries. The REPLACE entry 
 * specifies global atom name replacements. The most common example is 
 * HN to H. CHARMM27 denotes the amide H atom as HN, while OPLS-AA denotes
 * it as H. To translate between the two all "HN" entries should be 
 * translated to simply "H" before charge assignments are made. Most other
 * common replacements have to do with the naming conventions of terminal 
 * capping groups. A common entry for REPLACE would be 
 * 	REPLACE: 8
 * 		HN	H
 * 		O1	OT
 * 		O2	O
 * 		OT1	OT
 * 		OT2	O
 * 		HT2	HO
 * 		HT1	H1
 * 		HT2	H2
 * 	ENDREPLACE
 *
 * Note that all whitespace should be tabs, not spaces. 
 *
 * 
 * The second type of entry defines residues to be recognized by the program.
 * Ideally, the charge file will contain exactly one entry corresponding to 
 * each residue (both protein and non-protein) in the simulation. Additional
 * entries (beyond what is actually present in the simulation) will be read
 * in from the charge file but will not be used. Every residue definition 
 * MUST consist of 6 entries: 
 *
 * 	(1) Residue start. This entry declares that a new residue is being
 * 	introduced and assigns it a name (to be used internally by g_amide). 
 *
 * 		RESIDUE: HISE
 *
 * 	Note that the name assigned here need not be the name actually used
 * 	in the PDB or GRO structure file or the name assigned by gromacs 
 * 	in the simulation. The user has full discretion over this entry. 
 *
 * 	(2) Number of residue names. This entry specifies how many redundant
 * 	names will be recognized by g_amide() as corresponding to the same 
 * 	residue. Under the NAMES entry below, these names will be explicitly 
 * 	specified. For example, the gromacs implementation of the CHARMM27 
 * 	force field refers to the two neutral histidine isomers as HSD and HSE, 
 * 	while OPLS-AA refers to them as HISD and HISE. Since either name (or 
 * 	simply HIS) should be acceptable as specifying this same residue, 
 * 	NNAMES will be set to 3 in this case. 
 *
 * 		NNAMES: 3
 *
 *	The HISE example will be continued below under the NAMES entry. 
 *
 * 	(3) Number of atoms. This entry specifies how many atoms the residue 
 * 	contains. This number must match the actual number of atoms and charges
 * 	specfied under the ATOMS field. 
 *
 * 		NATOMS: 17
 *
 * 	(4) Residue names. This entry lists all possible acceptable aliases
 * 	which might be encountered in the gromacs input. The same entry may 
 * 	appear in this field for multiple residues; the number and naming of
 * 	atoms will in this case be used to distinguish between residues. For
 * 	example, for the epsilon isomer of neutral histidine, we might have 
 * 	the entry
 *
 * 		NAMES:
 * 			HISE
 * 			HSE
 * 			HIS
 *
 * 	even though the entry HIS also occurs for the charged and HISD forms
 * 	of histidine. Charged is distinguished from neutral HIS by the number
 * 	of atoms, while HISE is distinguished from HISD by the atom names. 
 *
 * 	(5) Atoms names and charges. This entry specifies the names and charges
 * 	of each atom in the residue. Atom names should be inset by one tab from
 * 	the start of the line, and charges should be offset one tab further. 
 * 	The first few entries of the CHARMM27 HISE file should be as follows: 
 *
 * 		N	-0.47
 * 		H	0.31
 * 		CA	0.07
 *
 * 	Note that although the defined name in CHARMM27 for the amide H atom is 
 * 	"HN", it should be specified as "H" in the chargefile, since this is 
 * 	the name used by g_amide() for amide hydrogens. 
 *
 * 	(6) Residue end. This field simply specifies that the residue definition
 * 	is finished. Its primary use is to detect user errors (e.g. declaring 
 * 	17 atoms but listing only 16 or declaring 17 atoms but listing 18). 
 *
 */



int allocate_residue(t_residue *p_res, int nnames, int natoms) {
	int maxchar = 10;
	t_residue res = *p_res;

	if(nnames==0) {
		printf("Error: attempted to allocate residue with no names.\n");
		return -1;
	} else if(natoms==0) { 
		printf("Error: attempted to allocate residue with no atoms.\n");
		return -1;
	}
	
	// Allocate charge array.
	res.atomcharges = (float*) malloc(natoms*sizeof(float));
	if(res.atomcharges==NULL) return -1;
	
	// Alocate residue name array.
	res.names = (char**) malloc(nnames*sizeof(char*));
	if(res.names==NULL) {
		free(res.atomcharges);
		return -1;
	}
	else {
		int i; 
		for(i=0; i<nnames; i++) {
			res.names[i] = (char*) malloc(maxchar*sizeof(char));
			if(res.names[i]==NULL) {
				int j;
				for(j=0; j<i; j++) free(res.names[j]);
				free(res.names);
				free(res.atomcharges);
				return -1;
			} 
		}
		res.nnames = nnames;
	}

	// Allocate atomname array. 
	res.atomnames = (char**) malloc(natoms*sizeof(char*));
	if(res.atomnames==NULL) {
		int j;
		for(j=0; j<res.nnames; j++) free(res.names[j]);
		free(res.names);
		free(res.atomcharges);
		return -1;
	}
	else {
		int i; 
		for(i=0; i<natoms; i++) {
			res.atomnames[i] = (char*) malloc(maxchar*sizeof(char));
			if(res.atomnames[i]==NULL) {
				int j;
				for(j=0; j<i; j++) free(res.atomnames[j]);
				free(res.atomnames);
				for(j=0; j<res.nnames; j++) free(res.names[j]);
				free(res.names);
				free(res.atomcharges);
				return -1;
			}
		}
		res.natoms = natoms;
	}
	*p_res = res;
	return 0;
}

int free_residue(t_residue res) {
	int j;
	for(j=0; j<res.nnames; j++) free(res.names[j]);
	for(j=0; j<res.natoms; j++) free(res.atomnames[j]);
	if(res.names!=NULL) free(res.names);
	if(res.atomnames!=NULL) free(res.atomnames);
	if(res.atomcharges!=NULL) free(res.atomcharges);
	return 0;
}

int free_resarray(t_residue *resarray, int nresidues) {
	if(resarray==NULL) return 0;
	int i;
	for(i=0; i<nresidues; i++) free_residue(resarray[i]);
	free(resarray);
	return 0;
}

int realloc_resarray(t_residue **p_resarray, int oldsize, int newsize) {
	t_residue *newarray;
	newarray = (t_residue*) realloc(*p_resarray, (newsize)*sizeof(t_residue));
	if(newarray==NULL) return -1;
	else {
		*p_resarray = newarray;
	}
	return 0;
}


// Read charge specification input file. See comments above for formatting. 
int read_charge_map(char *fname, t_residue **p_resarray, int* p_nres, char **repnames[2], int *p_nrep, int verbose) {
	int maxchar = 512;
	char line[maxchar];
	char name[64];
	char resname[64];
	int nalloc = 0;
	int allocchunk = 100;
	int nchar = 10;
	int nres = 0;
	int error = 0;

	int nrep = 0;
	repnames[0] = NULL;
	repnames[1] = NULL;
	t_residue *resarray = (t_residue*) malloc(allocchunk*sizeof(t_residue));	
	if(resarray==NULL) {
		printf("Error allocating residue array.\n");
		return -1;
	} else nalloc = allocchunk;
	
	FILE *fp = fopen(fname, "r");
	if(fp==NULL) {
		printf("Error opening charge file. \n");
		error = 1;
	}
	while( (fgets(line, maxchar, fp)!=NULL) && (!error) ) {
		// If we find a REPLACE statement, being adding global atom name replcements. 
		if( (!error) && (!strncmp(line, "REPLACE:", 8))) {
			if(sscanf(line, "%*s%d", &nrep)<1) {
				printf("Error reading charge file. Unable to parse number of global name replacements.\n");
				printf("The offending line was: %s\n", line);
				error = 1;
			}
			if(!error) {
				repnames[0] = (char**) malloc(nrep*sizeof(char*));
				if(repnames[0]==NULL) error = 1;
				else {
					int i; 
					for(i=0; i<nrep; i++) {
						repnames[0][i] = (char*) malloc(nchar*sizeof(char));
						if(repnames[0][i]==NULL) {
							int j;
							for(j=0; j<i; j++) free(repnames[0][j]);
							free(repnames[0]);
							error = 1;
							break;
						}
					}
				}
				if(error) printf("Error allocating memory for %d name replacements in function read_charge_map().\n", nrep);
			}
			if(!error) {
				repnames[1] = (char**) malloc(nrep*sizeof(char*));
				if(repnames[1]==NULL) error = 1;
				else {
					int i; 
					for(i=0; i<nrep; i++) {
						repnames[1][i] = (char*) malloc(nchar*sizeof(char));
						if(repnames[1][i]==NULL) {
							int j;
							for(j=0; j<nrep; j++) free(repnames[0][j]);
							free(repnames[0]);
							for(j=0; j<i; j++) free(repnames[1][j]);
							free(repnames[1]);
							error = 1;
							break;
						}
					}
				}
				if(error) { 
					printf("Error allocating memory for %d name replacements in function read_charge_map().\n", nrep);
					nrep = 0;
				}
			}
			if(!error) {
				char name1[maxchar];
				char name2[maxchar];
				int i;
				for(i=0; i<nrep; i++) {
					if(fgets(line, maxchar, fp)==NULL) error = 1;
					else if(sscanf(line, "%s%s", name1, name2)<2) { 
						printf("Error parsing name replacement at line: %s\n", line);
						error = 1;
					} else {
						strncpy(repnames[0][i], name1, 9);
						repnames[0][i][nchar-1] = '\0';
						strncpy(repnames[1][i], name2, 9);
						repnames[1][i][nchar-1] = '\0';
					}
				}
				if( (fgets(line, maxchar, fp)==NULL) || (strncmp(line, "ENDREPLACE", 10)!=0) ) {
					printf("Error parsing atom name replacements in function read_charge_map().\n");
					printf("Was expecting ENDREPLACE statement after REPLACE: entry.\n");
					printf("Last line read was: %s\n", line);
					error = 1;
				}
			}
		}
		// If we enter the if() statement, we index a new residue. 
		if( (!error) && (!strncmp(line, "RESIDUE:", 8))) {
			int nnames = 0;
			int natoms = 0;
			// Increase resarray size if needed.
			if(nres==nalloc) {
				if(realloc_resarray(&resarray, nalloc, nalloc+allocchunk)!=0) {
					printf("Error re-allocating residue array.\n");
					error = 1;
				} else nalloc+=allocchunk;
			}
			if(sscanf(line, "%*s%s", resname)<1) {
				printf("Error reading charge file. Unable to parse name for RESIDUE entry.\n");
				printf("The offending line is: %s\n", line);
				error = 1;
			}
			if(!error) {
				if(fgets(line, maxchar, fp)==NULL) error = 1;
				else if(strncmp(line, "NNAMES:", 7)!=0) {
					printf("Error reading charge file. Expected NNAMES entry after RESIDUE.\n");
					printf("Was attempting to read residue: %s\n", resname);
					error = 1;
				} else if(sscanf(line, "%*s%d", &nnames)<1) {
					printf("Error reading charge file. Couldn't parse NNAMES value.\n");
					printf("The offensive line is: %s\n", line);
					printf("Was attempting to read residue: %s\n", resname);
					error = 1;
				}
			}
			if(!error) {
				if(fgets(line, maxchar, fp)==NULL) error = 1;
				else if(strncmp(line, "NATOMS:", 7)!=0) {
					printf("Error reading charge file. Expected NATOMS entry after NNAMES.\n");
					printf("Was attempting to read residue: %s\n", resname);
					error = 1;
				} else if(sscanf(line, "%*s%d", &natoms)<1) {
					printf("Error reading charge file. Couldn't parse NATOMS value.\n");
					printf("The offensive line is: %s\n", line);
					printf("Was attempting to read residue: %s\n", resname);
					error = 1;
				}
			}
			if(!error) {
				if(allocate_residue(&resarray[nres], nnames, natoms)!=0) {
					printf("Residue memory allocation error when reading charge file.\n");
					error = 1;
				}
				else {
					strncpy(resarray[nres].resname, resname, 9);
					resarray[nres].resname[9] = '\0';
					nres++;
				}
			}
			if(!error) {
				if(fgets(line, maxchar, fp)==NULL) error = 1;
				else if(strncmp(line, "NAMES:", 6)!=0) {
					printf("Error reading charge file. Expected NAMES entry after NATOMS.\n");
					printf("Was attempting to read residue: %s\n", resname);
					error = 1;
				} else {
					int i;
					for(i=0; i<resarray[nres-1].nnames; i++) {
						if(fgets(line, maxchar, fp)==NULL) error = 1;
						else if(sscanf(line, "%s", name)<1) { 
							printf("Error parsing residue name at line: %s\n", line);
							printf("Was attempting to read residue: %s\n", resname);
							error = 1;
						} else {
							strncpy(resarray[nres-1].names[i], name, 9);
							resarray[nres-1].names[i][9] = '\0';
						}
					}
				}
			}
			if(!error) {
				if(fgets(line, maxchar, fp)==NULL) error = 1;
				else if(strncmp(line, "ATOMS:", 6)!=0) {
					printf("Error reading charge file. Expected ATOMS entry after NAMES.\n");
					printf("Was attempting to read residue: %s\n", resname);
					error = 1;
				} else {
					int i;
					float charge;
					for(i=0; i<resarray[nres-1].natoms; i++) {
						if(fgets(line, maxchar, fp)==NULL) error = 1;
						else if(sscanf(line, "%s%f", name, &charge)<2) { 
							printf("Error parsing atom at line: %s\n", line);
							printf("Was attempting to read residue: %s\n", resname);
							error = 1;
						} else {
							strncpy(resarray[nres-1].atomnames[i], name, 9);
							resarray[nres-1].atomnames[i][9] = '\0';
							resarray[nres-1].atomcharges[i] = charge;
						}
					}
					if(!error) {
						if(fgets(line, maxchar, fp)==NULL) error = 1;
						else if(strncmp(line, "ENDRESIDUE", 10)!=0) error = 1;
						if(error) {
							printf("Error reading charge file. Expected ENDRESIDUE marker after ATOMS entries.\n");
							printf("Was attempting to read residue %s.\n", resname);
						}
					}
				}
			}

		}
	}
	fclose(fp);

	if( (!error) && (realloc_resarray(&resarray, nalloc, nres)!=0)) {
		printf("Memory allocation error when trimming residue array.\n");
		error = 1;
	}

	if(error) {
		printf("Error reading charge file.\n");
		free_resarray(resarray, nres);
		if(nrep>0) {
			int i;
			for(i=0; i<nrep; i++) {
				if(repnames[0][i]!=NULL) free(repnames[0][i]);
				if(repnames[1][i]!=NULL) free(repnames[1][i]);
			}
			if(repnames[0]!=NULL) free(repnames[0]);
			if(repnames[1]!=NULL) free(repnames[1]);
		}
		*p_nres = 0;
		*p_nrep = 0;
		return -1;
	} else {
		printf("Successfully loaded charge map from input file %s.\n", fname);
		if(verbose) {
			printf("Located %d atom name replacements. \n", nrep);
			int i; 
			for(i=0; i<nrep; i++) printf("\t%s\t%s\n", repnames[0][i], repnames[1][i]);
			printf("\nLocated %d residues:\n\n", nres);
			for(i=0; i<nres; i++) {
				printf("RESIDUE: %s\n", resarray[i].resname);
				printf("NNAMES: %d\n", resarray[i].nnames);
				printf("NATOMS: %d\n", resarray[i].natoms);
				printf("NAMES: \n");
				int j;
				for(j=0; j<resarray[i].nnames; j++) printf("\t%s\n", resarray[i].names[j]);
				printf("ATOMS:\n");
				for(j=0; j<resarray[i].natoms; j++) printf("\t%s\t%f\n", resarray[i].atomnames[j], resarray[i].atomcharges[j]);
				printf("ENDRESIDUE\n\n");
			}
		}
	}
	
	*p_nrep = nrep;
	*p_nres = nres;
	*p_resarray = resarray;
	return 0;
}


// Parse through top and make any atom name replacements specified by repnames[][].
int swap_atomnames(t_topology top, char **repnames[2], int nrep, int verbose) {
	int at; 
	int i;
	for(at=0; at<top.atoms.nr; at++) {
		for(i=0; i<nrep; i++) {
			if(strcmp(repnames[0][i], *(top.atoms.atomname[at]))==0) {
				if(verbose) printf("Matched atom %d (%s) to replacement name %s\n", at, *(top.atoms.atomname[at]), repnames[1][i]);
				strncpy(*(top.atoms.atomname[at]), repnames[1][i], 5);
				*(*(top.atoms.atomname[at])+5) = '\0';
				break;
			}
		}
	}
	return 0;
}


// Run through all atoms in the system, identify the residues they belong to, and (if necessary) adjust atomic charges. 
int identify_residues(t_pbc *p_pbc, t_topology top, rvec *x, t_residue *resarray, int nres, int verbose) {
	int resind;
	for(resind=0; resind<top.atoms.nres; resind++) {
		int at;
		int natoms = 0;
		for(at=0; at<top.atoms.nr; at++) if(top.atoms.atom[at].resind==resind) natoms++;
		int res;
		int matchfail = 1;
		for(res=0; res<nres; res++) {
			if(natoms==resarray[res].natoms) {
				int i;
				for(i=0; i<resarray[res].nnames; i++) {
					if(strcmp(resarray[res].names[i], *(top.atoms.resinfo[resind].name))==0) {
						int resmatch = 1;
						for(at=0; at<top.atoms.nr; at++) {
							if(top.atoms.atom[at].resind==resind) {
								int atmatch = 0;
								int atx;
								for(atx=0; atx<resarray[res].natoms; atx++) {
									if(strcmp(resarray[res].atomnames[atx], *(top.atoms.atomname[at]))==0) atmatch = 1;
								}
								if(atmatch==0) {
									resmatch = 0;
								}
							}
						}
						if(resmatch==1) {
							if(verbose) printf("Matched residue %s %d to known residue %s.\n", *(top.atoms.resinfo[resind].name), top.atoms.resinfo[resind].nr, resarray[res].resname);
							matchfail = 0;
							for(at=0; at<top.atoms.nr; at++) {
								if(top.atoms.atom[at].resind==resind) {
									int atx;
									for(atx=0; atx<resarray[res].natoms; atx++) {
										if(strcmp(resarray[res].atomnames[atx], *(top.atoms.atomname[at]))==0) {
											if(top.atoms.atom[at].q!=resarray[res].atomcharges[atx]) {
												if(verbose) printf("Changing atom %s %d %s charge from %6.3f to %6.3f\n", *(top.atoms.resinfo[resind].name), resind, *(top.atoms.atomname[at]), top.atoms.atom[at].q, resarray[res].atomcharges[atx]);
											}
											top.atoms.atom[at].q = resarray[res].atomcharges[atx];
										}
									}
								}
							}
						}
					}
				}
			}
		}
		if(matchfail) {
			printf("Error: unable to match residue %s %d to a known residue.\n", *(top.atoms.resinfo[resind].name), top.atoms.resinfo[resind].nr);
			return -1;
		}
	}
	return 0;
}

