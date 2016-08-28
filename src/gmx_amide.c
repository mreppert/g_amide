/* gmx_amide.c is the main source file for the program g_amide, used for estimating
 * site energies, coupling constants, and transition dipole moment parameters 
 * for Amide I vibrations of protein and peptides systems. 
 *
 * Written by Mike Reppert (2015)
 */

#include "gmx_amide_map.h"
#include "gmx_amide_bond.h"
#include "gmx_amide_charge.h"
#include "mem_helper.h"

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
#include "string2.h"
#include "pbc.h"
#include "bondf.h"

#if OMP_PARALLEL
	#include "omp.h"
#endif

/********************************************************************************
* 				Site Assignments				*
********************************************************************************/

int find_map_sites( t_pbc pbc, t_topology top, rvec *x, t_amide_map map, t_protbond *p_pb, int bonds ) {
	int b,i,j,at,at0;
	for(b=0; b<bonds; b++) {
		p_pb[b].MapAtoms = (int**) malloc(map.nsites*sizeof(int*));
		if(p_pb[b].MapAtoms==NULL) {
			printf("Error allocating memory for MapAtoms of bond %d.\n", b+1);
			return 0;
		}
		p_pb[b].nsites = map.nsites;
		p_pb[b].natoms = (int*) malloc(map.nsites*sizeof(int));
		if(p_pb[b].natoms==NULL) {
			printf("Error allocating memory for AtomSite numbers of bond %d.\n", b+1);
			return 0;
		}
		for(i=0; i<map.nsites; i++) p_pb[b].natoms[i] = 0;
		for(i=0; i<map.nsites; i++) {
			p_pb[b].MapAtoms[i] = (int*) malloc(map.MapSites[i].natoms*sizeof(int));
			if(p_pb[b].MapAtoms[i]==NULL) {
				printf("Error allocating memory for MapAtoms of bond %d, map site %d.\n", b+1, i+1);
				return 0;
			}
			p_pb[b].natoms[i] = map.MapSites[i].natoms;
			for(j=0; j<map.MapSites[i].natoms; j++) {
				at0 = -1;
				if(strncmp(map.MapSites[i].AtomPaths[j].Path[0], "C", 6)==0) at0 = p_pb[b].C;
				else if(strncmp(map.MapSites[i].AtomPaths[j].Path[0], "O", 6)==0) at0 = p_pb[b].O;
				else if(strncmp(map.MapSites[i].AtomPaths[j].Path[0], "N", 6)==0) at0 = p_pb[b].N;
				else if(strncmp(map.MapSites[i].AtomPaths[j].Path[0], "H", 6)==0) at0 = p_pb[b].H;
				if(at0==-1) {
					printf("Error finding mapping atoms. Mapping atom path #%d starts at ", i);
					printf("an atom (%s) not found in the peptide bond.\n", map.MapSites[i].AtomPaths[j].Path[0]);
					return 0;
				}
				// Should there be an error-catch here to be sure that at has been assigned? 
				at = trace_path(&pbc, top, x, at0, map.MapSites[i].AtomPaths[j]);
				p_pb[b].MapAtoms[i][j] = at;
			}

		}
	}
	return 1;
}

int find_coup_sites( t_pbc pbc, t_topology top, rvec *x, t_amide_map map, t_protbond *p_pb, int bonds ) {
	int b,i,j,at,at0;
	for(b=0; b<bonds; b++) {
		p_pb[b].CoupAtoms = (int**) malloc(map.ncoupsites*sizeof(int*));
		if(p_pb[b].CoupAtoms==NULL) {
			printf("Error allocating memory for CoupAtoms of bond %d.\n", b+1);
			return 0;
		}
		p_pb[b].ncoupsites = map.ncoupsites;
		p_pb[b].ncoupatoms = (int*) malloc(map.ncoupsites*sizeof(int));
		if(p_pb[b].ncoupatoms==NULL) {
			printf("Error allocating memory for CoupSite numbers of bond %d.\n", b+1);
			return 0;
		}
		for(i=0; i<map.ncoupsites; i++) p_pb[b].ncoupatoms[i] = 0;
		for(i=0; i<map.ncoupsites; i++) {
			p_pb[b].CoupAtoms[i] = (int*) malloc(map.CoupSites[i].natoms*sizeof(int));
			if(p_pb[b].CoupAtoms[i]==NULL) {
				printf("Error allocating memory for CoupAtoms of bond %d, coupling site %d.\n", b+1, i+1);
				return 0;
			}
			p_pb[b].ncoupatoms[i] = map.CoupSites[i].natoms;
			for(j=0; j<map.CoupSites[i].natoms; j++) {
				at0 = -1;
				if(strncmp(map.CoupSites[i].AtomPaths[j].Path[0], "C", 6)==0) at0 = p_pb[b].C;
				else if(strncmp(map.CoupSites[i].AtomPaths[j].Path[0], "O", 6)==0) at0 = p_pb[b].O;
				else if(strncmp(map.CoupSites[i].AtomPaths[j].Path[0], "N", 6)==0) at0 = p_pb[b].N;
				else if(strncmp(map.CoupSites[i].AtomPaths[j].Path[0], "H", 6)==0) at0 = p_pb[b].H;
				if(at0==-1) {
					printf("Error finding coupling site atoms. Coupling atom path #%d starts at ", i);
					printf("an atom (%s) not found in the peptide bond.\n", map.CoupSites[i].AtomPaths[j].Path[0]);
					return 0;
				}
				// Should there be an error-catch here to be sure that at has been assigned? 
				at = trace_path(&pbc, top, x, at0, map.CoupSites[i].AtomPaths[j]);
				p_pb[b].CoupAtoms[i][j] = at;
			}

		}
	}
	return 1;
}

// Identify excluded atoms for each amide bond. 
int find_excluded_atoms( t_pbc pbc, t_topology top, rvec *x, t_amide_map map, t_protbond *p_pb, int bonds) {
	int b,i,at,at0;
	for(b=0; b<bonds; b++) {
		p_pb[b].nexcluded = 0;
		int nfound = 0;
		int nalloc = 0;
		int* ExArray; 

		// Our initial guess is that there are just map.npaths excluded atoms. Could be more or less.
		ExArray = (int*) malloc(map.npaths*sizeof(int));
		nalloc = map.npaths;
		if(ExArray==NULL) {
			printf("Error allocating memory for excluded atom indices of bond %d.\n", b+1);
			return 0;
		}
		for(i=0; i<map.npaths; i++) {
			at0 = -1;
			// First, identify the starting point for the path, i.e. the index of the atom
			// corresponding to the first path entry. 
			if(strncmp(map.ExcludedAtomPaths[i].Path[0], "C", 6)==0) at0 = p_pb[b].C;
			else if(strncmp(map.ExcludedAtomPaths[i].Path[0], "O", 6)==0) at0 = p_pb[b].O;
			else if(strncmp(map.ExcludedAtomPaths[i].Path[0], "N", 6)==0) at0 = p_pb[b].N;
			else if(strncmp(map.ExcludedAtomPaths[i].Path[0], "H", 6)==0) at0 = p_pb[b].H;
			else if(strncmp(map.ExcludedAtomPaths[i].Path[0], "CD", 6)==0) at0 = p_pb[b].H;
			if(at0==-1) {
				printf("Error finding excluded atoms. Excluded atom path #%d starts at ", i);
				printf("an atom (%s) not found in the peptide bond.\n", map.ExcludedAtomPaths[i].Path[0]);
				return 0;
			}
			// Next, trace the atom indices from at0 to the end of the path. 
			// If the path is incomplete (which it will be most of the time)
			// at<0 is returned. 
			at = trace_path(&pbc, top, x, at0, map.ExcludedAtomPaths[i]);
			if(at>=0) {
				int oldat;
				int isnew = 1;
				// If we find a new atom, first check if it's alredy included in ExArray.
				// If not, go ahead and add it. 
				for(oldat=0; oldat<nfound; oldat++) {
					if(ExArray[oldat]==at) isnew = 0;
				}
				if(isnew) {
					ExArray[nfound] = at;
					nfound++;
				}
			}
		}
		// Now trim the allocated memory to keep only the indices for found atoms.
		p_pb[b].Excluded = (int*) realloc(ExArray, nfound*sizeof(int));
		if( (p_pb[b].Excluded==NULL) && (nfound!=0) ) {
			printf("Error trimming excess memory from protein bond %d excluded atoms array.\n", b+1);
			p_pb[b].nexcluded = 0;
			return 0;
		}
		p_pb[b].nexcluded = nfound;
	}
	return 1;
}

// Identify Coupling sites for each amide bond
int find_coupling_sites( t_pbc pbc, t_topology top, rvec *x, t_amide_map map, t_protbond *p_pb, int bonds) {
	int b,i,at,at0;
	for(b=0; b<bonds; b++) {
		p_pb[b].nexcluded = 0;
		int nfound = 0;
		int nalloc = 0;
		int* ExArray; 

		// Our initial guess is that there are just map.npaths excluded atoms. Could be more or less.
		ExArray = (int*) malloc(map.npaths*sizeof(int));
		nalloc = map.npaths;
		if(ExArray==NULL) {
			printf("Error allocating memory for excluded atom indices of bond %d.\n", b+1);
			return 0;
		}
		for(i=0; i<map.npaths; i++) {
			at0 = -1;
			// First, identify the starting point for the path, i.e. the index of the atom
			// corresponding to the first path entry. 
			if(strncmp(map.ExcludedAtomPaths[i].Path[0], "C", 6)==0) at0 = p_pb[b].C;
			else if(strncmp(map.ExcludedAtomPaths[i].Path[0], "O", 6)==0) at0 = p_pb[b].O;
			else if(strncmp(map.ExcludedAtomPaths[i].Path[0], "N", 6)==0) at0 = p_pb[b].N;
			else if(strncmp(map.ExcludedAtomPaths[i].Path[0], "H", 6)==0) at0 = p_pb[b].H;
			else if(strncmp(map.ExcludedAtomPaths[i].Path[0], "CD", 6)==0) at0 = p_pb[b].H;
			if(at0==-1) {
				printf("Error finding excluded atoms. Excluded atom path #%d starts at ", i);
				printf("an atom (%s) not found in the peptide bond.\n", map.ExcludedAtomPaths[i].Path[0]);
				return 0;
			}
			// Next, trace the atom indices from at0 to the end of the path. 
			// If the path is incomplete (which it will be most of the time)
			// at<0 is returned. 
			at = trace_path(&pbc, top, x, at0, map.ExcludedAtomPaths[i]);
			if(at>=0) {
				int oldat;
				int isnew = 1;
				// If we find a new atom, first check if it's alredy included in ExArray.
				// If not, go ahead and add it. 
				for(oldat=0; oldat<nfound; oldat++) {
					if(ExArray[oldat]==at) isnew = 0;
				}
				if(isnew) {
					ExArray[nfound] = at;
					nfound++;
				}
			}
		}
		// Now trim the allocated memory to keep only the indices for found atoms.
		p_pb[b].Excluded = (int*) realloc(ExArray, nfound*sizeof(int));
		if( (p_pb[b].Excluded==NULL) && (nfound!=0) ) {
			printf("Error trimming excess memory from protein bond %d excluded atoms array.\n", b+1);
			p_pb[b].nexcluded = 0;
			return 0;
		}
		p_pb[b].nexcluded = nfound;
	}
	return 1;
}


/********************************************************************************
* 				Memory Allocation				*
********************************************************************************/

int set_arrays(real ****p_elecData, real ****p_ProElecData, int nelec, real ***p_angleData, real ***p_ProAngleData, int nangles, real ***p_freqData, real ***p_ProFreqData, int nfreq, real ****p_coupData, int ncoup, rvec ***p_coordData, rvec ***p_ProCoordData, rvec **p_dipData, rvec **p_centData, matrix **p_RotMat, matrix **p_ProRotMat, int bonds, int nsites, int nPro, int nProSites) {
	int error = 0;
	int i,ix; 

	*p_coordData = (rvec**) malloc(bonds*sizeof(rvec*));
	if(*p_coordData==NULL) error = 1;
	else {
		for(i=0; i<bonds; i++) {
			*(*p_coordData+i) = (rvec*) malloc(nsites*sizeof(rvec));
			if(*(*p_coordData+i)==NULL) { 
				error = 1;
				break;
			}
		}
		if(error) {
			if(i<nsites) {
				for(ix=0; ix<i; ix++) free(*(*p_coordData+ix));
				free(*p_coordData);
			}
			printf("Error allocating memory for Coordinate data array.\n");
		}
	}

	if(!error) {
		*p_dipData = (rvec*) malloc(bonds*sizeof(rvec));
		if(*p_dipData==NULL) {
			error = 2;
			printf("Error allocating memory for Dipole moment data.\n");
		}
	}
	if(!error) {
		*p_centData = (rvec*) malloc(bonds*sizeof(rvec));
		if(*p_centData==NULL) {
			error = 3;
			printf("Error allocating memory for bond Center data.\n");
		}
	}
	if(!error) {
		*p_RotMat = (matrix*) malloc(bonds*sizeof(matrix));
		if(*p_RotMat==NULL) {
			error = 4;
			printf("Error allocating memory for Rotation Matrix data.\n");
		}
	}
	if(!error) {
		if(!allocate_3d_array_real(p_elecData, bonds, nsites, nelec)) {
			printf("Error allocating Electrostatic variable array\n");
			error = 5;
		} else if(!allocate_2d_array_real(p_angleData, bonds, nangles)) {
			printf("Error allocating Angles variable array\n");
			error = 6;
		} else if(!allocate_2d_array_real(p_freqData, bonds, nfreq)) {
			printf("Error allocating Frequency variable array\n");
			error = 7;
		} else if(!allocate_3d_array_real(p_coupData, bonds, bonds, ncoup)) {
			printf("Error allocating Coupling variable array\n");
			error = 8;
		}
	}
	if(!error) {
		if(nProSites!=0) {
			if(!allocate_3d_array_real(p_ProElecData, nPro, nProSites, nelec)) {
				printf("Error allocating Proline Electrostatic variable array\n");
				error = 9;
			} else if(!allocate_2d_array_real(p_ProAngleData, nPro, nangles)) {
				printf("Error allocating Proline Angles variable array\n");
				error = 10;
			} else if(!allocate_2d_array_real(p_ProFreqData, nPro, nfreq)) {
				printf("Error allocating Proline Frequency variable array\n");
				error = 11;
			}
		} 
		*p_ProCoordData = (rvec**) malloc(nPro*sizeof(rvec*));
		if(*p_ProCoordData==NULL) error = 12;
		else {
			for(i=0; i<nPro; i++) {
				*(*p_ProCoordData+i) = (rvec*) malloc(nProSites*sizeof(rvec));
				if(*(*p_ProCoordData+i)==NULL) { 
					error = 13;
					break;
				}
			}
			if(error) {
				if(i<nProSites) {
					for(ix=0; ix<i; ix++) free(*(*p_ProCoordData+ix));
					free(*p_ProCoordData);
				}
				printf("Error allocating memory for Proline Coordinate data array.\n");
			}
		}
		if(!error) {
			*p_ProRotMat = (matrix*) malloc(nPro*sizeof(matrix));
			if(*p_ProRotMat==NULL) {
				error = 14;
				printf("Error allocating memory for Proline Rotation Matrix data.\n");
			}
		}
	}

	if(!error) return 1;
	else {
		if(error>1) {
			for(i=0; i<bonds; i++) free(*(*p_coordData+i));
			free(*p_coordData);
		}
		if(error>2) free(*p_dipData);
		if(error>3) free(*p_centData);
		if(error>4) free(*p_RotMat);
		if(error>5) free_3d_array_real(*p_elecData, bonds, nsites, nelec);
		if(error>6) free_2d_array_real(*p_angleData, bonds, nangles);
		if(error>7) free_2d_array_real(*p_freqData, bonds, nfreq);
		if(error>8) free_3d_array_real(*p_coupData, bonds, bonds, nfreq);
		if(error>9) free_3d_array_real(*p_ProElecData, nPro, nProSites, nelec);
		if(error>10) free_2d_array_real(*p_ProAngleData, nPro, nangles);
		if(error>11) free_2d_array_real(*p_ProFreqData, nPro, nfreq);
		if(error>12) {
			for(i=0; i<nPro; i++) free(*(*p_ProCoordData+i));
			free(*p_ProCoordData);
		}
		if(error>13) free(*p_ProRotMat);
		return 0;
	}
}

int unset_arrays(real ***elecData, real ***ProElecData, int nelec, real **angleData, real **ProAngleData, int nangles, real **freqData, real **ProFreqData, int nfreq, real ***coupData, int ncoup, rvec **coordData, rvec **ProCoordData, rvec *dipData, rvec *centData, matrix *RotMat, matrix *ProRotMat, int bonds, int nsites, int nPro, int nProSites) {
	int i;
	free_3d_array_real(elecData, bonds, nsites, nelec);
	free_2d_array_real(angleData, bonds, nangles);
	free_2d_array_real(freqData, bonds, nfreq);
	free_3d_array_real(coupData, bonds, bonds, ncoup);
	for(i=0; i<bonds; i++) free(coordData[i]);
	free(dipData);
	free(centData);
	free(RotMat);
	if(nProSites!=0) {
		free_2d_array_real(ProFreqData, nPro, nfreq);
		free_3d_array_real(ProElecData, nPro, nProSites, nelec);
		free_2d_array_real(ProAngleData, nPro, nangles);
		free(ProRotMat);
	}
	return 1;
}

/********************************************************************************
* 				Parameters					*
********************************************************************************/

int get_angles(t_pbc* p_pbc, t_topology top, rvec *x, t_amide_map map, t_protbond *p_pb, int bonds, real **angleData, int nangles) {
	const real rad2deg = 57.29578;
	rvec r_ij,r_jk,r_kl,m,n;
	real sign;
	int t1,t2,t3;
	int b;
	t_protbond pb;

	for(b=0; b<bonds; b++) {
		pb = p_pb[b];
		if(pb.dih_set_N) {
			angleData[b][0] = rad2deg*dih_angle(x[pb.phiN[0]], x[pb.phiN[1]], x[pb.phiN[2]], x[pb.phiN[3]], p_pbc, r_ij, r_jk, r_kl, m, n, &sign, &t1, &t2, &t3);
			angleData[b][1] = rad2deg*dih_angle(x[pb.psiN[0]], x[pb.psiN[1]], x[pb.psiN[2]], x[pb.psiN[3]], p_pbc, r_ij, r_jk, r_kl, m, n, &sign, &t1, &t2, &t3);
		} else {
			angleData[b][0] = -1000.0;
			angleData[b][1] = -1000.0;
		}
		if(pb.dih_set_C) {
			angleData[b][2] = rad2deg*dih_angle(x[pb.phiC[0]], x[pb.phiC[1]], x[pb.phiC[2]], x[pb.phiC[3]], p_pbc, r_ij, r_jk, r_kl, m, n, &sign, &t1, &t2, &t3);
			angleData[b][3] = rad2deg*dih_angle(x[pb.psiC[0]], x[pb.psiC[1]], x[pb.psiC[2]], x[pb.psiC[3]], p_pbc, r_ij, r_jk, r_kl, m, n, &sign, &t1, &t2, &t3);
		} else {
			angleData[b][2] = -1000.0;
			angleData[b][3] = -1000.0;
		}
	}
	return 1;
}

real get_NN_val(real phi, real psi, real **data, int dim1, int dim2) {

	real val;
	real f_phi, f_psi;
	real start = -180.0;
	real stop  = 180.0;
	real step1 = (stop - start) / ( (real) (dim1-1) );
	real step2 = (stop - start) / ( (real) (dim2-1) );

	int phi_ndx = (int) ( (phi-start)/step1 );
	int psi_ndx = (int) ( (psi-start)/step2 );
	f_phi = (phi - (start+phi_ndx*step1) ) / step1;
	f_psi = (psi - (start+psi_ndx*step2) ) / step2;
	if(phi_ndx==-1) phi_ndx = dim1-2;
	if(psi_ndx==-1) psi_ndx = dim2-2;
	if(phi_ndx==dim1-1) phi_ndx = 0;
	if(psi_ndx==dim2-1) psi_ndx = 0;
	val = 0.0;
	if( (phi_ndx>=0) && (psi_ndx>=0) && (phi_ndx<dim1) && (psi_ndx<dim2)) {
		val += (1-f_phi)*((1-f_psi)*data[phi_ndx][psi_ndx] + f_psi*data[phi_ndx][psi_ndx+1] );
		val +=   f_phi * ((1-f_psi)*data[phi_ndx+1][psi_ndx]+f_psi*data[phi_ndx+1][psi_ndx+1] );
	}
	return val;
}

int get_coordinates(t_pbc pbc, t_topology top, rvec *x, t_amide_map map, t_protbond *p_pb, int bonds, matrix *RotMat, rvec **coordData) {

	int b,i,j;
	rvec x_ax, y_ax, z_ax, CN_vec;	
	rvec dx; 

	for(b=0; b<bonds; b++) {
		// RotMat holds the molecular coordinate unit vectors
		pbc_dx(&pbc, x[p_pb[b].O], x[p_pb[b].C], x_ax);
		unitv(x_ax, x_ax);
		pbc_dx(&pbc, x[p_pb[b].N], x[p_pb[b].C], CN_vec);
		cprod(x_ax, CN_vec, z_ax);
		unitv(z_ax, z_ax);
		cprod(z_ax, x_ax, y_ax);
		unitv(y_ax, y_ax);

		for(i=0; i<3; i++) {
			RotMat[b][i][0] = x_ax[i];
			RotMat[b][i][1] = y_ax[i];
			RotMat[b][i][2] = z_ax[i];
		}

		clear_rvec(dx);
		for(i=0; i<map.nsites; i++) {
			clear_rvec(coordData[b][i]); // MapAtoms[i][0] is our reference point. 
			copy_rvec(x[p_pb[b].MapAtoms[i][0]], coordData[b][i]);
			for(j=1; j<map.MapSites[i].natoms; j++) {
				pbc_dx(&pbc, x[p_pb[b].MapAtoms[i][j]], x[p_pb[b].MapAtoms[i][0]], dx);
				svmul((1.0)/map.MapSites[i].natoms, dx, dx);
				rvec_inc(coordData[b][i], dx);
			}
			// Incorrect PBC treatment at boundaries. Replaced 08/25/2014 by MER. 
			//clear_rvec(coordData[b][i]);
			//for(j=0; j<map.MapSites[i].natoms; j++) rvec_inc(coordData[b][i], x[p_pb[b].MapAtoms[i][j]]); 
			//svmul((1.0)/map.MapSites[i].natoms, coordData[b][i], coordData[b][i]);
			//
			//  The lines below compare PBC-corrected and uncorrected values. 
			//rvec x0;
			//clear_rvec(x0);
			//for(j=0; j<map.MapSites[i].natoms; j++) rvec_inc(x0, x[p_pb[b].MapAtoms[i][j]]); 
			//svmul((1.0)/map.MapSites[i].natoms, x0, x0);
			//rvec_sub(x0, coordData[b][i], dx);
			//printf("PBC - NoPBC: (%6.10f, %6.10f, %6.10f)\n", dx[0], dx[1], dx[2]);

		}

	}
	return 1;
}

int get_electrostatics(t_pbc pbc, t_topology top, rvec *x, t_amide_map map, t_protbond *p_pb, matrix *RotMat, int bonds, rvec **coordData, real ***elecData, int nelec, int nthreads, real cutoff, int nchunks, int* START, int* STOP) {
	int b;
	int natoms = top.atoms.nr;
	real const nm2bohr = 18.8973;
	real bcutoff = cutoff*nm2bohr;
	int error = 0;
	
	#if OMP_PARALLEL
	omp_set_num_threads(nthreads);
	#pragma omp parallel for \
	shared(natoms, pbc, top, x, map, p_pb, RotMat, bonds, \
		coordData, elecData, nelec, nchunks, START, STOP)
	#endif
	for(b=0; b<bonds; b++) {
		int i,j,at, k;	// k for looping chunks
		real d, invd, invd2, invd3, invd5;
		real q;
		real DX, DY, DZ;
		rvec x_ax, y_ax, z_ax;
		rvec dx;

		for(i=0; i<3; i++) {
			x_ax[i] = RotMat[b][i][0];
			y_ax[i] = RotMat[b][i][1];
			z_ax[i] = RotMat[b][i][2];
		}
		for(i=0; i<nelec; i++) {
			if(map.elec_used[i]) {
				for(j=0; j<map.nsites; j++) elecData[b][j][i] = 0.0;
			}
		}
		for(at=0; at<natoms; at++) {
			int good = 1;

			// cjfeng 08/27/2016
			if(nchunks) {	// Check if we need to include only part of bath
				for(k=0; k<nchunks; k++) {
					if(at<START[k]) {	// less than the starting point
						good = 0;
					}
					else if(at>STOP[k]) {	// more than the ending point
						good = 0;
					}
					else {			// Within one of the region
						good = 1;
						break;
					}
				}
			}
			
			for(i=0; i<p_pb[b].nexcluded; i++) {
				if(at==p_pb[b].Excluded[i]) {
					good = 0;
					break;
				}
			}
			if(good) {
				for(i=0; i<map.nsites; i++) {
					pbc_dx(&pbc, coordData[b][i], x[at], dx);
					//if((b==55)) if(i==0) printf("%6.10f\t%6.10f\t%6.10f\n", dx[0], dx[1], dx[2]);
					svmul(nm2bohr, dx, dx);
					d = norm(dx);
					if( (d!=0.0) && (d<bcutoff) ) {
						DX = iprod(x_ax, dx);
						DY = iprod(y_ax, dx);
						DZ = iprod(z_ax, dx);
						invd = (1.0)/d;
						invd2 = invd*invd;
						invd3 = invd*invd2;
						invd5 = invd2*invd3;
						q = top.atoms.atom[at].q;
						if(map.elec_used[0]) elecData[b][i][0] += q * invd;			// Potential
						if(map.elec_used[1]) elecData[b][i][1] += ( q * DX ) * ( invd3 );	// x-field
						if(map.elec_used[2]) elecData[b][i][2] += ( q * DY ) * ( invd3 );	// y-field
						if(map.elec_used[3]) elecData[b][i][3] += ( q * DZ ) * ( invd3 );	// z-field
						if(map.elec_used[4]) elecData[b][i][4] += ( q*invd3 ) * ( 1 - 3*DX*DX*invd2 );	// xx
						if(map.elec_used[5]) elecData[b][i][5] -= ( 3*q ) * ( DX*DY*invd5 );		// xy
						if(map.elec_used[6]) elecData[b][i][6] -= ( 3*q ) * ( DX*DZ*invd5 );		// xz
						if(map.elec_used[7]) elecData[b][i][7] += ( q*invd3 ) * ( 1 - 3*DY*DY*invd2 );	// yy
						if(map.elec_used[8]) elecData[b][i][8] -= ( 3*q ) * ( DY*DZ*invd5 );		// yz
						if(map.elec_used[9]) elecData[b][i][9] += ( q*invd3 ) * ( 1 - 3*DZ*DZ*invd2 );	// zz
						//if(b==2) if(i==0) printf("Bond %d, atom %d: field %6.10f\n", b, at, elecData[b][i][1]);
						//if(b==55) if(i==0) printf("Bond %d, atom %d: charge %6.10f\n", b, at, q);
						//if(b==2) if(i==0) printf("Bond %d, atom %d: distance %6.10f\n", b, at, d);
					} else if( d==0.0 )  {
						printf("Division by zero when calculating electrostatics between atom %d (bond %d) and site %d.\n", at, b+1, i);
						error = 0;
					}
				}
			}
		}
	}
	return 1;
}


int get_freq(t_amide_map map, t_protbond *p_pb, int bonds, real ***elecData, int nelec, real **angleData, int nangles, real **freqData, int nfreq) {
	int b,i,j;
	real wel, wnn;
	for(b=0; b<bonds; b++) {
		wel = 0.0;
		wnn = 0.0;
		if(map.dimNNFSN[0]*map.dimNNFSN[1]>0) wnn += get_NN_val(angleData[b][0], angleData[b][1], map.NNFSN, map.dimNNFSN[0], map.dimNNFSN[1]);
		if(map.dimNNFSC[0]*map.dimNNFSC[1]>0) wnn += get_NN_val(angleData[b][2], angleData[b][3], map.NNFSC, map.dimNNFSC[0], map.dimNNFSC[1]);
		for(i=0; i<map.nsites; i++) {
			for(j=0; j<nelec; j++) {
				if(map.elec_used[j]) wel += elecData[b][i][j]*map.MapSites[i].shift[j];
			}
		}
		freqData[b][0] = wel + wnn + map.freq;
		if(p_pb[b].isPro) { 
			freqData[b][0] += map.proshift;
		}
		freqData[b][1] = wel;
		freqData[b][2] = wnn;
		//printf("Bond %d: %6.10f\n", b, elecData[b][0][1]);
		//for(i=0; i<3; i++) if(b==2) printf("%6.10f\n", elecData[b][0][i+1]);
	}
	return 1;
}

int get_dip(t_pbc pbc, rvec *x, t_amide_map map, t_protbond *p_pb, int bonds, real ***elecData, int nelec, matrix *RotMat, rvec *dipData, rvec *centData) {
	int b,i,j;
	rvec CN, CO;
	rvec mu;
	real fCO = 0.0665;	// Remember pdb units are nm, not Angstrom
	real fCN = 0.0258;
	for(b=0; b<bonds; b++) {
		pbc_dx(&pbc, x[p_pb[b].N], x[p_pb[b].C], CN);
		unitv(CN, CN);
		svmul(fCN, CN, CN);
		pbc_dx(&pbc, x[p_pb[b].O], x[p_pb[b].C], CO);
		unitv(CO, CO);
		svmul(fCO, CO, CO);
		copy_rvec(x[p_pb[b].C], centData[b]); 	// The C atom is our reference point. 
		rvec_inc(centData[b], CO);		// Augment position with CO and CN vectors. 
		rvec_inc(centData[b], CN);
		copy_rvec(map.dip, mu);
		for(i=0; i<map.nsites; i++) {
			for(j=0; j<nelec; j++) {
				if(map.elec_used[j]) {
					mu[0] += elecData[b][i][j]*map.MapSites[i].dipx[j];
					mu[1] += elecData[b][i][j]*map.MapSites[i].dipy[j];
					mu[2] += elecData[b][i][j]*map.MapSites[i].dipz[j];
				}
			}
		}
		mvmul(RotMat[b], mu, dipData[b]);
	}
	return 1;
}

int get_pdc(t_pbc pbc, t_protbond *p_pb, int bonds, rvec *dipData, rvec *centData, real ***coupData, int ncoup) {
	int b1, b2;
	const real A = 0.1*(848619.0 / 1650.0) / (0.1011*0.1011);
	// Dipole derivatives can be converted to transition dipole moments using the factor 
	// \sqrt( \hbar / 2*w_i ) = 4.1189e-25 m kg^{1/2} =  0.1011 Angstrom * amu^{1/2}
	// for a 1650 cm-1 vibration. A dipole derivative of 2.73 D/(A*amu^1/2) thus corresponds to 
	// a dipole moment matrix element of 2.73*0.1011 = 0.276 D. 
	//
	// The code below assumes that that the transition dipole moment is in units of Debye. 
	real d, invd;
	rvec mu1, mu2;
	rvec r12;
	for(b1=0; b1<bonds; b1++) {
		coupData[b1][b1][0] = 0.0;
		copy_rvec(dipData[b1], mu1);
		for(b2=0; b2<b1; b2++) {
			copy_rvec(dipData[b2], mu2);
			pbc_dx(&pbc, centData[b2], centData[b1], r12);
			d = 10.0*norm(r12);
			invd = 1.0/d;
			unitv(r12, r12);
			coupData[b1][b2][0] = A*invd*invd*invd*( iprod(mu1,mu2) - 3.0*( iprod(mu1, r12)*iprod(mu2,r12) ) );
			coupData[b2][b1][0] = coupData[b1][b2][0];
		}
	}
	return 0;
}


int get_pdc_skinner(t_pbc pbc, t_protbond *p_pb, int bonds, rvec *dipData, rvec *centData, real ***coupData, int ncoup) {
	int b1, b2;
	const real A = 0.1*(848619.0 / 1650.0) * (2.73*2.73);
	real d, invd;
	rvec mu1, mu2;
	rvec r12;
	for(b1=0; b1<bonds; b1++) {
		coupData[b1][b1][0] = 0.0;
		copy_rvec(dipData[b1], mu1);
		unitv(mu1, mu1);
		for(b2=0; b2<b1; b2++) {
			copy_rvec(dipData[b2], mu2);
			unitv(mu2, mu2);
			pbc_dx(&pbc, centData[b2], centData[b1], r12);
			d = 10.0*norm(r12);
			invd = 1.0/d;
			unitv(r12, r12);
			coupData[b1][b2][0] = A*invd*invd*invd*( iprod(mu1,mu2) - 3.0*( iprod(mu1, r12)*iprod(mu2,r12) ) );
			coupData[b2][b1][0] = coupData[b1][b2][0];
		}
	}
	return 0;
}

	/* 

	Electrostatic coupling models such as TCC seek to evaluate matrix elements of the form
	
		<i| H_C |j>

	where H_C is the Coulomb interaction energy between the various atoms of the two coupled
	groups. To develop the TCC model, we expand in normal mode displacements from Q_i = Q_j = 0. 
	By assumption, this starting point corresponds to the energy minimum, so that all first 
	derivatives are zero. Second derivatives dH/dQ_i^2 and dH/DQ_j^2 corresponds to site energy 
	shifts and are assumed to be accounted for by our electrostatic site energy maps. Mixed
	second derivatives d^2H/(dQ_i dQ_j), however, are expected to be non-zero and will form the 
	basis for our TCC model. From this starting point, we need to evaluate coupling elements of 
	the form 

		<i| Q_i Q_j * (d^2 H_C/ dQ_i dQ_j) |j> = (d^2 H_C/ dQ_i dQ_j) <i| Q_i Q_j |j>
			= (d^2 H_C/ dQ_i dQ_j) \hbar/(2*sqrt{w_i w_j})

	where the derivative term is evaluated at Q_i = Q_j = 0 and w_i and w_j are the frequencies
	of the i and j normal modes. 

	We next neglect what should be a full integral over the ground state charge density as a 
	function of adiabatic displacement along the normal mode coordinates and assume instead a 
	simple sum of Coulomb interaction terms where charges are assigned to each site within the 
	amide unit as a function of normal mode displacement:

		H_C = (1 / (4*\pi*\eps_0) ) * \sum_{m,n} ( q_n(Q_i) q_m(Q_j) ) / ( | r_n(Q_i) - r_m(Q_j) | )

	where q_n(Q_i) and r_n(Q_i) are atomic charges and displacements of atom n as a function of the 
	normal mode coordinate Q_i. The atomic displacements are by definition linear functions of the 
	normal modes, i.e. 

		r_n(Q_i) = r_n(0) + (dr_n/dQ_i) * Q_i. 

	This linear relationship need not be true for the atomic charges, but for simplicity (and 
	since we truncate at second order already), we make the same assumption for them: 

		q_n(Q_i) = q_n(0) + (dq_n/dQ_i) * Q_i.

	The zero-displacement value q_n(0) should ideally be the same set of charges as used in 
	electrostatic calculations for site energies, though discrepancies will not likely cause
	serious issues. 

	Within this model, the mixed second derivatives are easily evaluated, with the result

		4*\pi*\eps_0 * ( d^2 H_C / dQ_i dQ_j )  = \sum_{n,m} {
			+ (dq_n*dq_m) / | r_{nm} |
			+ ( q_n*q_m )*( dr_n \cdot dr_m) / | r_{nm} |^3
			-3( q_n*q_m )*( dr_n \cdot r_{nm} )*( dr_m \cdot r_{nm} ) / | r_{nm} |^5
			+ [ ( dq_n*q_m )*( dr_m \cdot r_{nm} ) - ( q_n*dq_m )*( dr_n \cdot r_{nm} ) ] / | r_{nm} |^5

	where for convenience we have denoted q_n(0) as q_n and (dq_n/dQ_i) as simply dq_n. The sums 
	over n and m are to be taken over the atoms contributing to normal modes i and j, respectively. 

	As input for the TCC calculation, we need a set of atomic charges q_n, charge fluxes dq_n, 
	and atomic displacements dr_n (i.e. the displacement of each atom as a function of 
	the normal mode coordinate). 

	As a ratio of two quantities with units length and length*mass^{1/2}, the units on 
	spatial displacements should be amu^{-1/2}. If U_{nx,j} is the unitary matrix
	required to diagonlize the Hessian matrix (in mass-weighted coordinates), these 
	coefficients can be obtained as 

		(dr_{nx}/dQ_j) = U_{nx,j}/sqrt{M_n}

	where r_{nx} is the x-coordinate of atom n, Q_j is the j^{th} normal mode and M_n is 
	the mass of atom n. 

	The charge flux should likewise be specified in units of e_o / ( Angstrom * amu^{1/2} ) 
	since the derivative is taken along the (mass-weighted) normal mode coordinate. Note 
	that the product of the dr_n and dq_n terms then has units of C / (Angstrom*amu). 
	The suspicious mass unit is eliminated by the matrix element <i| Q_i Q_j |j> 
	reported above. 

	As an example, for a simple C=O oscillator, the correct mass-weighted normal mode 
	coordinate is 

		Q_j = ( sqrt(M_O)*r_{Ox} - sqrt(M_C)*r_{Cx} ) / sqrt(2)
		  = (  2.828 * r_{Ox} -  2.450 * r_{Cx} )  

	In this case, U_{Ox,j} = 1/sqrt{2}, while U_{Cx,j} = -1/sqrt{2}. The map input file should
	contain the entries

		COUPSITES: 2
		C
		O	

		NUX:
		-0.20412
		0.17678

		Q:
		0.5
		-0.5
	
	The x-component of the dipole moment matrix element specified by these parameters is given by 

		[mu_i]_x = (d[\mu_x]/dQ_i) 
		   = \sqrt( \hbar / (2*w_i) ) * \sum_n ( q_n(0)*dr_{nx} + dq_n r_{nx}(0) )
	
	For a typical amide I frequency of w_i = 1650 cm-1, the prefactor evaluates to 

		\sqrt( \hbar / 2*w_i ) = 4.1189e-25 m kg^{1/2} =  0.1011 Angstrom * amu^{1/2}

	while for the atomic charges and displacements specified above (letting dq_n = 0)
	the sum over n becomes
	
		  Q_C*( 1/sqrt{2*M_C} ) + Q_O*( 1/sqrt{2*M_O} )
		= 0.5 * ( -0.20412 ) - 0.5 * ( 0.17678 ) 
		= -0.1905 e_o / amu^{1/2}

	giving together a dipole moment magnitude of 

		0.1011 * 0.1905 = 0.0193 e_o*Angstrom = 0.0193*(1.60217657e-19 C/e_o)*(1e-10 m/Angstrom)
				= (3.0922e-31 C m) / (3.33564e-30 C m/Debye) = 0.0927 Debye

	Note that compared to the experimental value of ~0.35 Debye, this model severely
	underestimates the Amide I oscillator strength (see Ackels et al. Vib. Spec. 50 (2009) 2–9). 
	The discrepancy is due to the neglect of charge flow within the bond. During C=O dispacement
	the bonding characteristics of the amide group change, shifting from a C=O double bond
 	resonance structure in which all atoms have formal charge of zero toward a C--O single 
	bond structure in which the N obtains a positive charge and the O a negative one:
  
                              (+)
              H-N           H-N 
                 \   <-->     \\
                  C=O          C--O 
                                  (-)

	To obtain a simple estimate for the atomic charges required to obtain the experimentally
	observed transition dipole moment strength, take as coordinates for a typical amide bond

		r_C = [0.000  0.000  0.000];
		r_O = [1.230 -0.000  0.000];
		r_N = [-0.65  1.170  0.000];
		r_H = [-1.64  1.080  0.100];

	in units of Angstrom. Let us suppose q_n(0) charges of 0.5, -0.5, -0.3, +0.3 for C, O, N, 
	and H, respectively, and that the the normal mode motion involves only the C and O atoms. 
	The contribution to the dipole moment from the static charges is then the same as calculated 
	above, while the contribution from charge motion remains to be assigned. Assuming as a 
	simplest model that charge flows only from N to O, we can characterize the charge flux 
	with a single coefficient dq_N = -dq_Oj = dq > 0. The contribution to the x-component of the 
	dipole moment then becomes
	
		   \sqrt( \hbar / (2*w_i) ) * \sum_n ( q_n(0)*dr_{nx} + dq_n r_{nx}(0) )
		   = 0.1011 * dq * ( r_{Nx} - r_{Ox} )
		   = 0.1901 * {dq} e_o*Angstrom
		   = 3.0457e-30 {dq} C*m
		   = 0.9131 {dq} Debye
	
	(where {dq} represents the numerical value of the charge flux parameter without its units)
	and to the y-component 

		   0.1011 * dq * ( r_{Ny} - r_{Oy} )
		   =  0.1183 * {dq} e_o*Angstrom
		   = 1.8954e-30 {dq} C*m
		   = 0.5682 {dq} Debye. 

	The total oscillator strength is then

		(0.0927+0.9131 dq)^2 + (0.5682 dq)^2 . 

	Setting equal to the experimental value of 0.12 Debye^2, we solve for dq to obtain

		0.0086 + 0.1693 dq + 1.157 dq^2 = 0.12
		1.157 dq^2 + 0.1693 dq - 0.12 = 0

			--> dq = 0.257 e_o / ( Angstrom * amu^{1/2} )

	A minimal (experimentally feasible) TCC model for amide I would thus consist of the mapfile entries

		COUPSITES: 4
		C
		O
		N
		H

		NUX:
		-0.20412
		0.17678
		0.0
		0.0

		Q:
		0.5
		-0.5
		-0.4
		0.4

		DQ: 
		-0.257
		0.0
		0.257
		0.0


 	To get the value of 1/(4*pi*Eps0) used below, start with Eps0 in standard units:
		Eps0 = 8.85418e-12 F/m
	In what follows, curly brackets {} represent numerical values (i.e. a quantity minus its units). 
	Let F = Farads, m = meters, C = coloumbs, J = Joules, cm-1 = wavenumbers, Ang = angstrom, eo = 
	units of elementary charge. We'll need the fundamental constants
		Eps0 = 8.8541878176e-12 	F/m
		h = 6.62606957e-34 		J s
		c = 2.99792458e10 		cm / s
		E0 = 1.602176565e-19		C
	and the conversion relations
		{X} m = {X*1e10} Ang
		{X} F = {X} C^2/J
		{X} J = {X/hc} cm-1
		{X} C = {X/EO} eo
	Applying these in order, we get
		{X} F/m = {X/1e10} F/Ang 
			= {( X )/( 1e10 )} ( C^2 )/( Ang J )
			= {( X h c )/( 1e10 )} ( C^2 )/( Ang cm-1 )
			= {( X h c )/( EO^2 1e10 )} ( eo^2 )/( Ang cm-1 )
	So we have 
		4*pi*Eps0 = {( 4*pi*(8.8541878176e-12)*(6.62606957e-34)*(2.99792458e10) )/( (1e10)*(1.602176565e-19)^2 )} ( eo^2 )/( Ang cm-1 )
			


	Finally, we note that literature normal modes and charge fluxes are frequently reported in "unitless" coordinates, in 
	which the mass-weighted normal mode coefficients dr/dQ have been multiplied by a factor with mass units. These
	values can be converted to standard units (for input in the map file) by (1) multiplying the reported normal mode 
	coefficient by the square root of the atomic mass, (2) normalizing the resultant vector to have unit length (thus
	producing the corresponding column of the unitary transformation matrix which transforms mass-weighted coordinates
	into mass-weighted normal modes), and (3) dividing by the square root of the mass to produce dr/dQ. The charge 
	fluxes must likewise be re-scaled by the same factor. For the Jansen TCC map parameters reported in [J. Chem. Phys. 125,
	044312 (2006)], the conversion can be done in Matlab by entering the commands: 

	dq = [0.01668, -0.02845, -0.01530, 0.01736, 0.00008, 0.00963]';
	nu = 0.028074*[0, 0, 0; -0.831, 0.105, 0.0; 0.517, -0.047, 0.0; 0.074, -0.036, 0.0; 0.073, -0.133, 0.0; 0.0, 0.0, 0.0]; 
	NU = [nu(:,1); nu(:,2); nu(:,3)];
	m = [12, 12, 16, 14, 2, 12]'; M = [m; m; m];
	U = NU.*sqrt(M); U = U/norm(U); nu = (U./sqrt(M)); dq = dq*(nu(2)/NU(2));

	The resultant values nu and dq should be 
	
	dq = [0.165604207569485
	     -0.282460413989918
	     -0.151903140036758
	      0.172355458237785
	      0.000794264784505922
	      0.0956096234349003]

	and 

	nu = [0.0000
	     -0.231622444056777
	      0.144102050032917
	      0.0206258253432028
	      0.0203470979737001
	      0.0
	      0.0
	      0.0292663737977878
	     -0.0131001863666288
	     -0.0100341853020987
	     -0.0370707401438645
	      0.0
	      0.0
	      0.0
	      0.0
	      0.0
	      0.0
	      0.0]
	
	*/

int get_tcc(t_pbc pbc, t_amide_map map, rvec *x, t_protbond *p_pb, int bonds, matrix *RotMat, real ***coupData, int ncoup) {
	int i,j,n,m,k;
	rvec dx, ri, rj, rij, nu, nui, nuj;
	real d, invd, invd3, invd5, J;
	real qn, qm, dqn, dqm;
	real ONEOVER4PIEPS0 = 1.161409733881537e+05; 	// units of ( cm-1 * Ang )/( eo^2 )
	// For coupling purposes, we take wi = wj = 1650 cm-1. 
	real prefac = 0.0102168287; 			// \hbar/(2*sqrt{w_i w_j}) in Angstrom^2 * amu

	for(i=0; i<bonds; i++) {
		for(j=0; j<i; j++) {
			J = 0.0;
			int divbyzero = 0;
			for(n=0; n<map.ncoupsites; n++) {
				nu[0] = map.CoupSites[n].nux;
				nu[1] = map.CoupSites[n].nuy;
				nu[2] = map.CoupSites[n].nuz;
				mvmul(RotMat[i], nu, nui);
				copy_rvec(x[p_pb[i].CoupAtoms[n][0]], ri);
				for(k=1; k<p_pb[i].ncoupatoms[n]; k++) {
					pbc_dx(&pbc, x[p_pb[i].CoupAtoms[n][k]], x[p_pb[i].CoupAtoms[n][0]], dx);
					svmul((1.0)/p_pb[i].ncoupatoms[n], dx, dx);
					rvec_inc(ri, dx);
				}
				// Dipole moment components in e_o * Angstrom
				for(m=0; m<map.ncoupsites; m++) {
					nu[0] = map.CoupSites[m].nux;
					nu[1] = map.CoupSites[m].nuy;
					nu[2] = map.CoupSites[m].nuz;
					mvmul(RotMat[j], nu, nuj);
					copy_rvec(x[p_pb[j].CoupAtoms[m][0]], rj);
					for(k=0; k<p_pb[j].ncoupatoms[m]; k++) {
						pbc_dx(&pbc, x[p_pb[j].CoupAtoms[m][k]], x[p_pb[j].CoupAtoms[m][0]], dx);
						svmul((1.0)/p_pb[j].ncoupatoms[m], dx, dx);
						rvec_inc(rj, dx);
					}
					pbc_dx(&pbc, ri, rj, rij); // points towards i bond
					svmul(10.0, rij, rij);
					d = norm(rij);
					if(d==0) {
						//printf("Warning: Division by zero when calculating transition charge coupling constant for bonds %d and %d.\n", i, j);
						//printf("Setting coupling constant to zero...\n");
						divbyzero = 1;
						break;
					}
					invd = (1.0)/d;
					invd3 = invd*invd*invd;
					invd5 = invd3*invd*invd;
					qn = map.CoupSites[n].q;
					qm = map.CoupSites[m].q;
					dqn = map.CoupSites[n].dq;
					dqm = map.CoupSites[m].dq;
					J += prefac*ONEOVER4PIEPS0*(dqn*dqm)*invd;
					J -= prefac*ONEOVER4PIEPS0*(3*qn*qm)*iprod(nui,rij)*iprod(nuj,rij)*invd5;
					J -= prefac*ONEOVER4PIEPS0*( -dqn*qm*iprod(nuj,rij) + qn*dqm*iprod(nui,rij) - qn*qm*iprod(nui,nuj) )*invd3;
				}
			}
			if(divbyzero) J = 0.0;
			coupData[i][j][0] = J;
			coupData[j][i][0] = J;
			//if((i==2) && (j==0)) printf("TCC (%d, %d): %6.10f\n", i, j, J);
		}
	}
	return 0;
}

int get_tc_dip(t_pbc pbc, rvec *x, t_amide_map map, t_protbond *p_pb, int bonds, matrix *RotMat, rvec *dipData, rvec *centData) {

	int i,n,k;
	rvec dx, ri, nu, nui;
	rvec dip, ref;
	// real pi = 3.14159265359; 
	// real c = 2.99792458e+10; 		// cm/s
	// real hbar = (6.62606957e-34) / (2.0 * pi); 	// J*s
	// real wo = 1650.0*c*2.0*pi; 		// 1/s
	real prefac = 0.485498038384979; // sqrt( hbar / (2.0*wo) ) * ( 1.60217657e-19 ) / ( sqrt( 1.660538921e-27 ) * ( 3.33564e-30 ) );
	// \sqrt( \hbar / 2*w_i ) * ( 1.60217657e-19 C/eo ) / ( sqrt( 1.660538921e-27 amu/kg) * ( 3.33564e-30 C*m/Debye ) )

	for(i=0; i<bonds; i++) {
		// We use CoupAtoms[0][0] as our spatial reference point. 
		// All coordinates MUST be updated relative to this point; 
		// otherwise PBC may give errors. 
		copy_rvec(x[p_pb[i].CoupAtoms[0][0]], ref);
		for(k=0; k<3; k++) dip[k] = 0.0;
		for(n=0; n<map.ncoupsites; n++) {
			nu[0] = map.CoupSites[n].nux;
			nu[1] = map.CoupSites[n].nuy;
			nu[2] = map.CoupSites[n].nuz;
			mvmul(RotMat[i], nu, nui);
			copy_rvec(x[p_pb[i].CoupAtoms[n][0]], ri);
			for(k=1; k<p_pb[i].ncoupatoms[n]; k++) {
				pbc_dx(&pbc, x[p_pb[i].CoupAtoms[n][k]], x[p_pb[i].CoupAtoms[n][0]], dx);
				svmul((1.0)/p_pb[i].ncoupatoms[n], dx, dx);
				rvec_inc(ri, dx);
			}
			// Now treat PBC: Take pbc_dx of ri from ref to get correct pbc-treated spatial displacement.
			pbc_dx(&pbc, ri, ref, dx);
			// Then add (without pbc-correction) this displacement to ref to get the pbc-corrected ri. 
			// This treatment assures that the relative positions of all coupling sites are consistent
			// with one another within a single pbc image. 
			for(k=0; k<3; k++) ri[k] = ref[k] + dx[k];
			// Dipole moment components in e_o * Angstrom
			for(k=0; k<3; k++) dip[k] += (map.CoupSites[n].q*nui[k] + 10.0*map.CoupSites[n].dq*ri[k]);
		}
		for(k=0; k<3; k++) dipData[i][k] = prefac*dip[k];
	}
	return 0;
}


real get_avg_tc_dip(t_pbc pbc, rvec *x, t_amide_map map, t_protbond *p_pb, int bonds, matrix *RotMat) {

	int i,n,k;
	rvec dx, ri, nu, nui;
	rvec dip;
	real prefac = 0.485498038384979; // sqrt( hbar / (2.0*wo) ) * ( 1.60217657e-19 ) / ( sqrt( 1.660538921e-27 ) * ( 3.33564e-30 ) );
	// \sqrt( \hbar / 2*w_i ) * ( 1.60217657e-19 C/eo ) / ( sqrt( 1.660538921e-27 amu/kg) * ( 3.33564e-30 C*m/Debye ) )
	real avg = 0.0;
	rvec ref; 

	for(i=0; i<bonds; i++) {
		// We use CoupAtoms[0][0] as our spatial reference point. 
		// All coordinates MUST be updated relative to this point; 
		// otherwise PBC may give errors. 
		copy_rvec(x[p_pb[i].CoupAtoms[0][0]], ref);
		for(k=0; k<3; k++) dip[k] = 0.0;
		for(n=0; n<map.ncoupsites; n++) {
			nu[0] = map.CoupSites[n].nux;
			nu[1] = map.CoupSites[n].nuy;
			nu[2] = map.CoupSites[n].nuz;
			mvmul(RotMat[i], nu, nui);
			copy_rvec(x[p_pb[i].CoupAtoms[n][0]], ri);
			for(k=1; k<p_pb[i].ncoupatoms[n]; k++) {
				pbc_dx(&pbc, x[p_pb[i].CoupAtoms[n][k]], x[p_pb[i].CoupAtoms[n][0]], dx);
				svmul((1.0)/p_pb[i].ncoupatoms[n], dx, dx);
				rvec_inc(ri, dx);
			}
			// Now treat PBC: Take pbc_dx of ri from ref to get correct pbc-treated spatial displacement.
			pbc_dx(&pbc, ri, ref, dx);
			// Then add (without pbc-correction) this displacement to ref to get the pbc-corrected ri. 
			// This treatment assures that the relative positions of all coupling sites are consistent
			// with one another within a single pbc image. 
			for(k=0; k<3; k++) ri[k] = ref[k] + dx[k];

			// Dipole moment components in e_o * Angstrom
			for(k=0; k<3; k++) dip[k] += (map.CoupSites[n].q*nui[k] + 10.0*map.CoupSites[n].dq*ri[k]);
		}
		real val = 0.0;
		for(k=0; k<3; k++) val += prefac*prefac*dip[k]*dip[k];
		avg += sqrt(val)/( (real) bonds );
	}
	return avg;
}


int get_tcc_jansen(t_pbc pbc, t_amide_map map, rvec *x, t_protbond *p_pb, int bonds, matrix *RotMat, real ***coupData, int ncoup) {
	/* To get the value of 1/(4*pi*Eps0) used below, start with Eps0 in standard units:
		Eps0 = 8.85418e-12 F/m
	In what follows, curly brackets {} represent numerical values (i.e. a quantity minus its units). 
	Let F = Farads, m = meters, C = coloumbs, J = Joules, cm-1 = wavenumbers, Ang = angstrom, eo = 
	units of elementary charge. We'll need the fundamental constants
		Eps0 = 8.8541878176e-12 	F/m
		h = 6.62606957e-34 		J s
		c = 2.99792458e10 		cm / s
		E0 = 1.602176565e-19		C
	and the conversion relations
		{X} m = {X*1e10} Ang
		{X} F = {X} C^2/J
		{X} J = {X/hc} cm-1
		{X} C = {X/EO} eo
	Applying these in order, we get
		{X} F/m = {X/1e10} F/Ang 
			= {( X )/( 1e10 )} ( C^2 )/( Ang J )
			= {( X h c )/( 1e10 )} ( C^2 )/( Ang cm-1 )
			= {( X h c )/( EO^2 1e10 )} ( eo^2 )/( Ang cm-1 )
	So we have 
		4*pi*Eps0 = {( 4*pi*(8.8541878176e-12)*(6.62606957e-34)*(2.99792458e10) )/( (1e10)*(1.602176565e-19)^2 )} ( eo^2 )/( Ang cm-1 )
			
	*/
	int ndxi[6], ndxj[6];
	int i,j,n,m;
	real ONEOVER4PIEPS0 = 1.161409733881537e+05; 	// units of ( cm-1 * Ang )/( eo^2 )
	rvec rij, nui, nuj;
	real d, invd, invd3, invd5, J;

	real amp = 0.028074;			// units of Angstrom/( normal mode amplitude )
	// Vector nu is the normal mode coordinate in units of the normal mode amplitude
	// just defined. nu*amp has units of Angstrom.
	rvec nu[6] = { { 0.0, 0.0, 0.0}, 
			{-0.831, 0.105, 0.0}, 
			{0.517, -0.047, 0.0}, 
			{0.074, -0.036, 0.0}, 
			{0.073, -0.133, 0.0}, 
			{0.0, 0.0, 0.0} };

	// units of elementary charge
	real q[6] = {0.11072, 0.37173, -0.53632, -0.48418, 0.24278, 0.29527};
	real dq[6] = {0.01668, -0.02845, -0.01530, 0.01736, 0.00008, 0.00963};

	for(i=0; i<bonds; i++) {
		ndxi[0] = p_pb[i].CAN;
		ndxi[1] = p_pb[i].C;
		ndxi[2] = p_pb[i].O;
		ndxi[3] = p_pb[i].N;
		ndxi[4] = p_pb[i].H;
		ndxi[5] = p_pb[i].CAC;
		for(j=0; j<i; j++) {
			// Make sure they're not nearest neighbors--otherwise, get division by zero 
			if( (p_pb[i].CAN!=p_pb[j].CAC) && (p_pb[j].CAN!=p_pb[i].CAC) ) {
				ndxj[0] = p_pb[j].CAN;
				ndxj[1] = p_pb[j].C;
				ndxj[2] = p_pb[j].O;
				ndxj[3] = p_pb[j].N;
				ndxj[4] = p_pb[j].H;
				ndxj[5] = p_pb[j].CAC;
				J = 0.0;
				for(n=0; n<6; n++) {
					mvmul(RotMat[i], nu[n], nui);
					for(m=0; m<6; m++) {
						mvmul(RotMat[j], nu[m], nuj);
						pbc_dx(&pbc, x[ndxi[n]], x[ndxj[m]], rij); // points towards i bond
						svmul(10.0, rij, rij);
						d = norm(rij);
						invd = (1.0)/d;
						invd3 = invd*invd*invd;
						invd5 = invd3*invd*invd;
						J += ONEOVER4PIEPS0*(dq[n]*dq[m])*invd;
						J -= ONEOVER4PIEPS0*amp*amp*(3*q[n]*q[m])*iprod(nui,rij)*iprod(nuj,rij)*invd5;
						J -= ONEOVER4PIEPS0*( -dq[n]*q[m]*amp*iprod(nuj,rij) + q[n]*dq[m]*amp*iprod(nui,rij) - q[n]*q[m]*amp*amp*iprod(nui,nuj) )*invd3;
					}
				}
				coupData[i][j][0] = J;
				coupData[j][i][0] = J;
			}
			coupData[i][i][0] = 0.0;
		}
	}
	return 0;
}

int get_nnc(t_amide_map map, t_protbond *p_pb, int bonds, real **angleData, real ***coupData) {
	real Jnn;
	int b1, b2;
	int NN;
	for(b1=0; b1<bonds; b1++) {
		for(b2=0; b2<b1; b2++) {
			NN = 0;
			// Check if the N-terminal CA of b1 is the C-terminal CA of b2
			// Otherwise, check if the N-terminal CA of b2 is the C-terminal CA of b1
			if( ((p_pb[b1].dih_set_N*p_pb[b2].dih_set_C)==1) && (p_pb[b1].psiN[1]==p_pb[b2].psiC[1]) ) {
				Jnn = get_NN_val(angleData[b1][0], angleData[b1][1], map.NNC, map.dimNNC[0], map.dimNNC[1]);
				NN = 1;
			} else if( ((p_pb[b1].dih_set_C*p_pb[b2].dih_set_N)==1) && (p_pb[b2].psiN[1]==p_pb[b1].psiC[1]) ) {
				Jnn = get_NN_val(angleData[b1][2], angleData[b1][3], map.NNC, map.dimNNC[0], map.dimNNC[1]);
				NN = 1;
			}
			if(NN) {
				coupData[b1][b2][0] = Jnn;
				coupData[b2][b1][0] = Jnn;
				coupData[b1][b2][2] = Jnn;
				coupData[b2][b1][2] = Jnn;
			}
		}
	}
	return 1;
}

int get_dnnc(t_amide_map map, t_protbond *p_pb, int bonds, real **angleData, real ***coupData) {
	real Jnn;
	int b1, b2;
	int NN;
	for(b1=0; b1<bonds; b1++) {
		for(b2=0; b2<b1; b2++) {
			NN = 0;
			// Check if the N-terminal CA of b1 is the C-terminal CA of b2
			// Otherwise, check if the N-terminal CA of b2 is the C-terminal CA of b1
			if( ((p_pb[b1].dih_set_N*p_pb[b2].dih_set_C)==1) && (p_pb[b1].psiN[1]==p_pb[b2].psiC[1]) ) {
				Jnn = get_NN_val(angleData[b1][0], angleData[b1][1], map.DNNC, map.dimDNNC[0], map.dimDNNC[1]);
				NN = 1;
			} else if( ((p_pb[b1].dih_set_C*p_pb[b2].dih_set_N)==1) && (p_pb[b2].psiN[1]==p_pb[b1].psiC[1]) ) {
				Jnn = get_NN_val(angleData[b1][2], angleData[b1][3], map.DNNC, map.dimDNNC[0], map.dimDNNC[1]);
				NN = 1;
			}
			if(NN) {
				coupData[b1][b2][0] += Jnn;
				coupData[b2][b1][0] += Jnn;
				coupData[b1][b2][2] += Jnn;
				coupData[b2][b1][2] += Jnn;
			}
		}
	}
	return 1;
}

int open_info_file(FILE **p_fp, char* namebase, const int maxchar) {
	char infonm[maxchar];
	if( (namebase!=NULL && (strlen(namebase)>0) )) {
		strcpy(infonm, namebase);
		if(namebase[strlen(namebase)-1]!='/') {
			strcat(infonm, "_");
		}
	}
	else {
		infonm[0] = '\0';
	}
	strcat(infonm, "info.txt");
	*p_fp = fopen(infonm, "w");
	if(*p_fp==NULL) {
		printf("Error opening output file for writing.\n");
		return 0;
	} else return 1;
}

int open_spec_files(FILE **p_fp, int nfiles, char* namebase, const int maxchar) {
	int i;
	char sitenm[maxchar];
	char hamnm[maxchar];
	char dipxnm[maxchar];
	char dipynm[maxchar];
	char dipznm[maxchar];
	if( (namebase!=NULL && (strlen(namebase)>0) )) {
		strcpy(sitenm, namebase);
		strcpy(hamnm, namebase);
		strcpy(dipxnm, namebase);
		strcpy(dipynm, namebase);
		strcpy(dipznm, namebase);
		if(namebase[strlen(namebase)-1]!='/') {
			strcat(sitenm, "_");
			strcat(hamnm, "_");
			strcat(dipxnm, "_");
			strcat(dipynm, "_");
			strcat(dipznm, "_");
		}
	}
	else {
		sitenm[0] = '\0';
		hamnm[0] = '\0';
		dipxnm[0] = '\0';
		dipynm[0] = '\0';
		dipznm[0] = '\0';
	}
	strcat(sitenm, "sites.txt");
	strcat(hamnm, "ham.txt");
	strcat(dipxnm, "dipx.txt");
	strcat(dipynm, "dipy.txt");
	strcat(dipznm, "dipz.txt");
	p_fp[0] = fopen(sitenm, "w");
	p_fp[1] = fopen(hamnm, "w");
	p_fp[2] = fopen(dipxnm, "w");
	p_fp[3] = fopen(dipynm, "w");
	p_fp[4] = fopen(dipznm, "w");
	for(i=0; i<nfiles; i++) {
		if(p_fp[i]==NULL) {
			printf("Error opening output file for writing.\n");
			return 0;
		}
	}
	return 1;
}

int write_spec_data(FILE **p_fp, real ***coupData, real **freqData, rvec *dipData, int bonds) {
	int b1,b2;
	for(b1=0; b1<bonds; b1++) {
		fprintf(p_fp[0], "%10.6f\t", freqData[b1][0]);
		for(b2=0; b2<bonds; b2++) {
			if(b1==b2) fprintf(p_fp[1], "%10.6f\t", freqData[b1][0]);
			else fprintf(p_fp[1], "%10.6f\t", coupData[b1][b2][0]);
		}
		fprintf(p_fp[2], "%10.6f\t", dipData[b1][0]);
		fprintf(p_fp[3], "%10.6f\t", dipData[b1][1]);
		fprintf(p_fp[4], "%10.6f\t", dipData[b1][2]);
	}
	fprintf(p_fp[0], "\n");
	fprintf(p_fp[1], "\n");
	fprintf(p_fp[2], "\n");
	fprintf(p_fp[3], "\n");
	fprintf(p_fp[4], "\n");
	return 1;
}

int close_spec_files(FILE **p_fp, int nfiles) {
	int i;
	for(i=0; i<nfiles; i++) fclose(p_fp[i]);
	return 1;
}

int open_elec_files(FILE **elecfp[10], char *outname, t_amide_map map, const int maxchar) {
	int i,j,k,l;
	char fname[maxchar];
	char COMP[10][6] = { "P", "Ex", "Ey", "Ez", "Gxx", "Gxy", "Gxz", "Gyy", "Gyz", "Gzz" };
	for(j=0; j<10; j++) {
		if(map.elec_used[j]) {
			elecfp[j] = (FILE**) malloc(map.nsites*sizeof(FILE*));
			if(elecfp[j]==NULL) {
				printf("Error opening electrostatic output files. Please check input.\n");
				int jx;
				for(jx=0; jx<j; jx++) free(elecfp[j]);
				return 0;
			}
		} else elecfp[j] = NULL;
	}
	for(j=0; j<10; j++) {
		for(i=0; i<map.nsites; i++) {
			if(map.elec_used[j]) {
				fname[0] = '\0';
				if( (outname!=NULL) && (strlen(outname)>0) ) {
					strcpy(fname, outname);
					if(outname[strlen(outname)-1]!='/') strcat(fname, "_");
				}
				strcat(fname, COMP[j]);
				for(k=0; k<map.MapSites[i].natoms; k++) {
					if(k>0) strcat(fname, "_AND");
					for(l=0; l<map.MapSites[i].AtomPaths[k].length; l++) {
						strcat(fname, "_");
						strcat(fname, map.MapSites[i].AtomPaths[k].Path[l]);
					}
				}
				strcat(fname, ".txt");
				elecfp[j][i] = fopen(fname, "w");
				if(elecfp[j][i]==NULL) {
					printf("Error opening electrostatic output files. Please check input.\n");
					int ix, jx;
					for(ix=0; ix<i; ix++) fclose(elecfp[j][ix]);
					free(elecfp[j]);
					for(jx=0; jx<j; jx++) {
						for(ix=0; ix<map.nsites; ix++) fclose(elecfp[jx][ix]);
						free(elecfp[jx]);
					}
					return 0;
				}
			}
		}
	}
	return 1;
}

int open_angle_files(FILE *anglefp[10], char *outname, const int maxchar) {
	int i;
	char fname[maxchar];
	char NAME[4][6] = {"phiN", "psiN", "phiC", "psiC"};
	for(i=0; i<4; i++) {
		fname[0] = '\0';
		if( (outname!=NULL) && (strlen(outname)>0) ) {
			strcpy(fname, outname);
			if(outname[strlen(outname)-1]!='/') strcat(fname, "_");
		}
		strcat(fname, NAME[i]);
		strcat(fname, ".txt");
		anglefp[i] = fopen(fname, "w");
		if(anglefp[i]==NULL) {
			printf("Error opening angle output files. Please check input.\n");
			int ix;
			for(ix=0; ix<i; ix++) fclose(anglefp[ix]);
			return 0;
		}
	}
	return 1;
}

int write_elec_dat(FILE **elecfp[10], t_amide_map map, real ***elecData, int bonds, int nelec) {
	int b,i,j;
	for(i=0; i<nelec; i++) {
		if(map.elec_used[i]) {
			for(j=0; j<map.nsites; j++) {
				for(b=0; b<bonds; b++) {
					fprintf(elecfp[i][j], "%6.10f\t", elecData[b][j][i]);
				}
				fprintf(elecfp[i][j], "\n");
			}
		}
	}
	return 1;
}

int write_angle_dat(FILE *anglefp[4], real **angleData, int bonds, int nangles) {
	int b,i;
	for(i=0; i<nangles; i++) {
		for(b=0; b<bonds; b++) 	fprintf(anglefp[i], "%6.10f\t", angleData[b][i]);
		fprintf(anglefp[i], "\n");
	}
	return 1;
}

int close_elec_files(FILE **elecfp[10], t_amide_map map) {
	int i,j;
	for(j=0; j<10; j++) {
		if(elecfp[j]!=NULL) {
			for(i=0; i<map.nsites; i++) if(elecfp[j][i]!=NULL) fclose(elecfp[j][i]);
			free(elecfp[j]);
		}
	}
	return 1;
}

int close_angle_files(FILE *anglefp[4]) {
	int i;
	for(i=0; i<4; i++) if(anglefp[i]!=NULL) fclose(anglefp[i]);
	return 1;
}

/********************************************************************************
* 				Main code.					*
********************************************************************************/

int main ( int argc, char * argv[] ) {
	/************************************************************************
	This section deals with command line arguments. Most of the action 
	happens in the parse_common_args() command, which checks what command
	line arguments were provided and initializes the corresonding variables 
	**************************************************************************/

	// Variables set by command line arguments
	char* mapfile = NULL;
	char* promapfile = NULL;
	char* chargefile = NULL;
	char* outname = NULL;
	// cjfeng 08/27/2016
	// Added chunk command to selectively choose portion of bath
	// to include in electrostatics.
	char* chunkstr = "";
	int nthreads = 1;
	int print_elec = 0;
	int print_angles = 0;
	int nnc, nnfsn, nnfsc, dnnc;
	real cutoff = 1000.0;
	real osc = -1;
	int i;
	int verbose = 0;


	// A list of command line file flags
	t_filenm fnm[] = {
		{ efTRX, "-f",   NULL,     ffREAD  },
		{ efTPX, NULL,   NULL,     ffREAD  }
	};

	// A list of additional command line arguments
	t_pargs pa [] = {
		{"-mapfile", FALSE, etSTR, {&mapfile}, 
		 "Complete map file, specifying all parameters."},
		{"-chargefile", FALSE, etSTR, {&chargefile}, 
		 "Complete charge file, specifying all charges."},
		{"-promapfile", FALSE, etSTR, {&promapfile}, 
		 "Complete map file, specifying all parameters."},
		{"-outname", FALSE, etSTR, {&outname},
		 "File name base for output"},
		{"-nt", FALSE, etINT, {&nthreads},
		 "Number of threads"},
		{"-print_elec", FALSE, etBOOL, {&print_elec},
		 "Print electrostatic values"},
		{"-print_angles", FALSE, etBOOL, {&print_angles},
		"Print dihedral angles"},
		{"-cutoff", FALSE, etREAL, {&cutoff}, 
		 "Cutoff distance (nm)"}, 
		{"-osc", FALSE, etREAL, {&osc}, 
		"Oscillator strength (Debye^2). If positive, used to normalize dipole moments."},
		{"-verbose", FALSE, etBOOL, {&verbose}, 
		"Verbose (print lots of info)"},
		{"-chunk", FALSE, etSTR, {&chunkstr},
		 "Chunk range for selecting part of the system to be included for electrostatic frequency shift, format: [a-b;...;c-d]. Indexing begins at zero."}
		};	

	// The program description
	const char *desc[] = {""};

	// A description of known bugs
	const char *bugs[] = {""};

	// Here we parse the command line arguments and make them available to call within the program
	output_env_t oenv;
	CopyRight(stderr, argv[0]);
	parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_BE_NICE, asize(fnm), fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs, &oenv);

	/**********************************************************
 	* We first check if only part of the bath is included for *
 	* computing electrostatic frequency shift.                *
 	* cjfeng 08/27/2016                                       *
	**********************************************************/
	int START[1000];	// Starting index for computing electrostatics
	int STOP[1000];		// Ending index for computing electrostatics
	int nchunks = 0;	// Number of chunks

	if( strlen(chunkstr)!=0 ) {
		int i, j;
		printf("Parsing chunk request %s...\n", chunkstr);
		if( (chunkstr[0]!='[') || (chunkstr[strlen(chunkstr)-1]!=']') ) {
			printf("Error parsing chunk string. Format: [a-b;c-d;...;e-f]\n");
			return 0;
		} else {
			// Changed string processing code. MER 10/04/2014
			char chunk[1000];
			for(i=0; i<strlen(chunkstr); i++) {
				if( (chunkstr[i]==';') || (chunkstr[i]=='[') ) {
					for(j=i+1; j<strlen(chunkstr); j++) {
						if( (chunkstr[j]==';') || (chunkstr[j]==']') ) {
							strncpy(chunk, chunkstr+i+1, j-i-1);
							chunk[j-i-1] = '\0';
							break;
						}
					}
					int start, stop;
					if( (j>=strlen(chunkstr)) || (sscanf(chunk, "%d-%d", &start, &stop)!=2) ) {
						printf("Error parsing input segment %s\n", chunk);
						return 0;
					} else {
						START[nchunks] = start;
						STOP[nchunks] = stop;
						nchunks++;
					}
				}
			}
		}
	}

	/************************************************************************
	Now we check the input. First we check and set the mapping parameters,
	then the (first frame of the) provided structure file. 
	**************************************************************************/
	
	// First load map parameters
	t_amide_map map;
	t_amide_map promap;
	if( (mapfile==NULL) ) {
		printf("Please specify a map using the -mapfile flag.\n");
		return 0;
	}

	// read_amide_map opens the provided map file and stores all data in the map structure
	if(!read_amide_map(mapfile, &map, verbose)) {
		printf("Error reading specified map file. Please check input.\n");
		return 0;
	}
	printf("\nSuccessfully loaded all mapping parameters from input file %s. \n", mapfile);

	if((promapfile!=NULL)) {
		if(!read_amide_map(promapfile, &promap, verbose)) {
			printf("Error reading specified proline map file. Please check input.\n");
			return 0;
		}
		printf("\nSuccessfully loaded all mapping parameters from proline map file %s. \n", promapfile);
		printf("\nWarning: Currently only site energy shifts are used from the proline map file. Coupling parameters (both through-space and through-bond) are taken from the default Amide I map.\n");
	} else {
		promap.nsites = 0;
	}

	
	// Check whether nearest-neighbor coupling is used.
	if((map.dimNNC[0]*map.dimNNC[1])>0) nnc = 1;
	else nnc = 0;
	if((map.dimDNNC[0]*map.dimDNNC[1])>0) dnnc = 1;
	else dnnc = 0;
	
	// Check whether nearest-neighbor frequency shifts are used. 
	if((map.dimNNFSN[0]*map.dimNNFSN[1])>0) nnfsn = 1;
	else nnfsn = 0;
	if((map.dimNNFSC[0]*map.dimNNFSC[1])>0) nnfsc = 1;
	else nnfsc = 0;
	
	if(verbose) {
		printf("Coupling model: %s\n", map.couptype);
		if(nnc) printf("Using nearest-neighbor coupling map.\n");
		if(dnnc) printf("Using nearest-neighbor coupling difference map.\n");
		if((nnc+dnnc)==0) printf("No nearest-neighbor coupling map in use.\n");
		if(nnfsn) printf("Using N-terminal nearest-neighbor frequency-shift map.\n");
		if(nnfsc) printf("Using C-terminal nearest-neighbor frequency-shift map.\n");
		if((nnc+dnnc)==0) printf("No nearest-neighbor frequency shifts in use.\n");
	}
	if(map.dipset) {
		if(verbose) printf("Using dipole vector provided in map file\n");
	} else if(map.tcset) {
		if(verbose) printf("Using TCC parameters provided in map file to set dipole moment.\n");
	} else {
		printf("Error! No method specified for calculating dipole moments. \n");
		printf("Please provide either a dipole moment (DIP:) entry or TCC parameters.\n");
		free_amide_map(map);
		return 0;
	}

	// Read atomic charges.
	int ffcharge = 0;
	int nres;
	t_residue *resarray;
	int nrep;
	char **reparray[2];
	if( (chargefile!=NULL) ) {
		ffcharge = 0;
		if(read_charge_map(chargefile, &resarray, &nres, reparray, &nrep, verbose)!=0) {
			printf("Error reading charge file. \n");
			free_amide_map(map);
			free_amide_map(promap);
			return 0;
		}
	} else {
		// If no chargefile is specified, we use force-field charges. 
		// We still need to manually enter a replacement name array since
		// some names (particularly HN for the CHARMM27 force field) will 
		// often need to be swapped. 
		int error = 0;
		int j;
		ffcharge = 1;
		int nchars = 20;
		nrep = 8;
		nres = 0;
		resarray = NULL; 
		reparray[0] = NULL;
		reparray[1] = NULL; 
		reparray[0] = (char**) malloc(nrep*sizeof(char*));
		reparray[1] = (char**) malloc(nrep*sizeof(char*));
		if(reparray[0]==NULL) error = 1;
		if(reparray[1]==NULL) error = 1;
		for(i=0; i<nrep; i++) {
			if(!error) { 
				reparray[0][i] = (char*) malloc(nchars*sizeof(char));
				if(reparray[0][i]==NULL) {
					for(j=0; j<i; j++) free(reparray[0][j]);
					error = 2;
				}
			} else if(reparray[0]!=NULL) reparray[0][i] = NULL; 
			if(!error) { 
				reparray[1][i] = (char*) malloc(nchars*sizeof(char));
				if(reparray[1][i]==NULL) {
					for(j=0; j<nrep; j++) free(reparray[1][j]);
					for(j=0; j<i; j++) free(reparray[1][j]);
					error = 2;
				}
			} else if(reparray[1]!=NULL) reparray[1][i] = NULL; 
		}
		if(!error) {
			char *DEFREP1[] = { "HN", "O1", "O2", "OT1", "OT2", "HT2", "HT1", "HT2" };
			char *DEFREP2[] = { "H", "OT", "O", "OT", "O", "HO", "H1", "H2"};
			for(i=0; i<nrep; i++) {
				strcpy(reparray[0][i], DEFREP1[i]);
				strcpy(reparray[1][i], DEFREP2[i]);
			}
		}
		if(error) {
			if(reparray[0]!=NULL) free(reparray[0]);
			if(reparray[1]!=NULL) free(reparray[1]);
			free_amide_map(map);
			free_amide_map(promap);
			return 0;
		}
	}

	// Now structure files. First topology, then first frame of trajectory
	t_topology top;
	t_inputrec ir;
	t_trxstatus *status;
	int trrStatus = 1;
	real t;
	rvec *x;
	matrix box;
	int natoms;
	int ePBC;
	t_pbc pbc;

	// read_tpx_top() reads the topology file and stores the values in the t_topology type top. The t_inputrec type ir stores information
	// on the box and periodic boundary conditions. The matrix box contains the vectors defining the periodic box. 
	ePBC = read_tpx_top(ftp2fn(efTPX,asize(fnm),fnm), &ir, box, &natoms, NULL, NULL, NULL, &top);
	
	// set_pbc initializes a t_pbc type variable pbc to be used later on in calls to, e.g. pbc_dx for calculating pbc-corrected distances.
	set_pbc(&pbc, ePBC, box);

	// read_first_x opens a trajectory file,  allocates memory for one frame of coordinates (the rvec array x),
	// and reads the first frame from the trajectory.
	natoms = read_first_x(oenv, &status, ftp2fn(efTRX,asize(fnm),fnm), &t, &x, box);
	
	 // Check to see if the number of atoms in topology and trajectory files match. 
	if ( natoms > top.atoms.nr )	gmx_fatal(FARGS,"Topology (%d atoms) does not match trajectory (%d atoms)", top.atoms.nr, natoms);



	/************************************************************************
	Now parse the protein structure to assign charges and identify amide bonds. 
	**************************************************************************/

	int error = 0;
	t_protbond *p_pb = NULL;	// Protein bond array
	int b;				// b is a counter, bonds is the number of bonds.
	int bonds = 0;
	t_protbond *p_ppb = NULL;	// Proline protein bond array
	int nPro = 0;

	// First, make any atom name replacements necessary.
	if(!error) swap_atomnames(top, reparray, nrep, verbose);
	
	// Now reassign atomic charges.
	if( (!error) && (!ffcharge) ) if(identify_residues(&pbc, top, x, resarray, nres, verbose)!=0) error = 1;

	// And look for bonds. 
	if(!error) {
		bonds = find_bonds(&pbc, top, x, &p_pb);
		printf("\nLocated %d amide bonds in the peptide chain\n", bonds);
		if(bonds<1) { 
			printf("Error: No amide bonds located.\n");
			error = 1;
		} else {
			for(b=0; b<bonds; b++) {
				if(p_pb[b].isPro) nPro++;
			}
		}
		printf("\nLocated %d proline bonds in the peptide chain\n", nPro);
		if(promap.nsites!=0) {
			p_ppb = (t_protbond*) malloc( sizeof(t_protbond)*nPro );
			if(p_ppb==NULL) {
				printf("Error allocating memory for proline bonds.\n");
				error = 1;
			}	
			int count = 0;
			for(b=0; b<bonds; b++) {
				if(p_pb[b].isPro) {
					initialize_bond(&pbc, top, x, &p_ppb[count], p_pb[b].N, p_pb[b].H, p_pb[b].C, p_pb[b].O, p_pb[b].isPro);
					count++;
				}
			}
		}
	}

	// Next identify mapping sites for each bond
	if(!error) if(!find_map_sites(pbc, top, x, map, p_pb, bonds)) error = 1;
	if((!error) && (promap.nsites!=0)) if(!find_map_sites(pbc, top, x, promap, p_ppb, nPro)) error = 1;
	
	// And identify coupling sites for each bond
	// For now, coupling sites are taken ONLY from the main map file.
	// Coupling sites in the proline map file are ignored. 
	if(!error) if(!find_coup_sites(pbc, top, x, map, p_pb, bonds)) error = 1;
	if( (!error) && (promap.nsites!=0) ) if(!find_coup_sites(pbc, top, x, map, p_ppb, nPro)) error = 1;

	// And identify excluded atoms for electrostatic calculation
	if(!error) if(!find_excluded_atoms(pbc, top, x, map, p_pb, bonds)) error = 1;
	if((!error) && (promap.nsites!=0)) if(!find_excluded_atoms(pbc, top, x, promap, p_ppb, nPro)) error = 1;

	// If anything went wrong, surrender with dignity. 
	if(error) {
		for(b=0; b<bonds; b++) free_bond(p_pb[b]);
		if(p_pb!=NULL) free(p_pb);
		if(p_ppb!=NULL) free(p_ppb);
		if(reparray[0]!=NULL) for(i=0; i<nrep; i++) if(reparray[0][i]!=NULL) free(reparray[0][i]);
		if(reparray[1]!=NULL) for(i=0; i<nrep; i++) if(reparray[1][i]!=NULL) free(reparray[1][i]);
		if(reparray[0]!=NULL) free(reparray[0]);
		if(reparray[1]!=NULL) free(reparray[1]);
		free_resarray(resarray, nres);
		free_amide_map(map);
		free_amide_map(promap);
		return 0;
	}

	// Tell the user what we found. 
	if(verbose) for(b=0; b<bonds; b++) print_bond(p_pb[b], top);
	if(verbose) for(b=0; b<nPro; b++) print_bond(p_ppb[b], top);

	/************************************************************************
	Allocate memory for data arrays and open output files. 
	**************************************************************************/

	int arraysset = 0; 
	int specfilesopen = 0;
	int elecfilesopen = 0;
	int anglefilesopen = 0;
	real ***elecData;
	real ***ProElecData;
	const int nelec = 10;
	real **angleData;
	real **ProAngleData;
	const int nangles = 4;
	real **freqData;
	real **ProFreqData;
	const int nfreq = 3;
	real ***coupData;
	const int ncoup = 3;
	rvec **coordData;
	rvec **ProCoordData;
	rvec *dipData =  NULL;
	rvec *centData = NULL; 
	matrix *RotMat = NULL; 
	matrix *ProRotMat = NULL; 
	if(!set_arrays(&elecData, &ProElecData, nelec, &angleData, &ProAngleData, nangles, &freqData, &ProFreqData, nfreq, &coupData, ncoup, &coordData, &ProCoordData, &dipData, &centData, &RotMat, &ProRotMat, bonds, map.nsites, nPro, promap.nsites)) {
		printf("Error allocating memory for data arrays\n");
		error = 1;
	} else arraysset = 1;

	const int maxchar = 1024;
	FILE *infofp = NULL;
	int infofileopen = 0;
	if(!open_info_file(&infofp,outname,maxchar)) error = 1;
	else infofileopen = 1;
	FILE *specfp[5];
	if(!open_spec_files(specfp,5,outname,maxchar)) error = 1;
	else specfilesopen = 1;
	
	FILE **elecfp[10];
	if( (print_elec) && (!open_elec_files(elecfp, outname, map, maxchar))) error = 1;
	else if(print_elec) elecfilesopen = 1;
	
	FILE *anglefp[4];
	if( (print_angles) && (!open_angle_files(anglefp, outname, maxchar))) error = 1;
	else if(print_angles) anglefilesopen = 1;

	fprintf(infofp, "BONDS: %d\n", bonds);
	for(b=0; b<bonds; b++) fprintf(infofp, "%s %d %s %d\n", p_pb[b].resname1, p_pb[b].resnum1, p_pb[b].resname2, p_pb[b].resnum2);
	
	if( (strcmp(map.couptype, "pdc")!=0) && (strcmp(map.couptype, "tcc")!=0) ) {
		printf("Error: Please specify a coupling model (pdc or tcc).\n");
		error = 1;
	}

	if(error) {
		if(anglefilesopen==1) close_angle_files(anglefp);
		if(elecfilesopen==1) close_elec_files(elecfp, map);
		if(infofileopen==1) fclose(infofp);
		if(specfilesopen==1) close_spec_files(specfp, 5);
		if(arraysset==1) unset_arrays(elecData, ProElecData, nelec, angleData, ProAngleData, nangles, freqData, ProFreqData, nfreq, coupData, ncoup, coordData, ProCoordData, dipData, centData, RotMat, ProRotMat, bonds, map.nsites, nPro, promap.nsites);
		if(reparray[0]!=NULL) for(i=0; i<nrep; i++) if(reparray[0][i]!=NULL) free(reparray[0][i]);
		if(reparray[1]!=NULL) for(i=0; i<nrep; i++) if(reparray[1][i]!=NULL) free(reparray[1][i]);
		if(reparray[0]!=NULL) free(reparray[0]);
		if(reparray[1]!=NULL) free(reparray[1]);
		free_resarray(resarray, nres);
		free_amide_map(map);
		for(b=0; b<bonds; b++) free_bond(p_pb[b]);
		if(p_pb!=NULL) free(p_pb);
		return 0;
	}


	/************************************************************************
	Step through trajectory and calculate parameters. 
	**************************************************************************/
	real dipstrength = 0.0; 
	do {
		// Read coordinates for all amide bonds and set molecular frame axes/rotation matrix. 
		get_coordinates(pbc, top, x, map, p_pb, bonds, RotMat, coordData);
		if(promap.nsites!=0) get_coordinates(pbc, top, x, promap, p_ppb, nPro, ProRotMat, ProCoordData);
		// Calculate electrostatic values around each bond. 
		get_electrostatics(pbc, top, x, map, p_pb, RotMat, bonds, coordData, elecData, nelec, nthreads, cutoff, nchunks, START, STOP);
		if(promap.nsites!=0) get_electrostatics(pbc, top, x, promap, p_ppb, ProRotMat, nPro, ProCoordData, ProElecData, nelec, nthreads, cutoff, nchunks, START, STOP);
		// Calculate dihedral angles for each bond. 
		get_angles(&pbc, top, x, map, p_pb, bonds, angleData, nangles);
		if(promap.nsites!=0) get_angles(&pbc, top, x, promap, p_ppb, nPro, ProAngleData, nangles);
		// Calculate frequencies. 
		get_freq(map, p_pb, bonds, elecData, nelec, angleData, nangles, freqData, nfreq);
		if(promap.nsites!=0) get_freq(promap, p_ppb, nPro, ProElecData, nelec, ProAngleData, nangles, ProFreqData, nfreq);
		// Graft in proline frequencies
		if(promap.nsites!=0)  {
			for(b=0; b<bonds; b++) {
				if(p_pb[b].isPro) {
					freqData[b][0] = ProFreqData[p_pb[b].isPro-1][0];
					freqData[b][1] = ProFreqData[p_pb[b].isPro-1][1];
				}
			}
		}
		if( (osc>0) && (map.dipset) && (dipstrength==0.0) ) {
			dipstrength = sqrt(map.dip[0]*map.dip[0] + map.dip[1]*map.dip[1] + map.dip[2]*map.dip[2]);
			printf("Average dipole strength: %6.4f Debye\n", dipstrength);
			printf("Re-scaling map dipole moment.\n");
			for(i=0; i<3; i++) map.dip[i] *= (sqrt(osc)/dipstrength);
			dipstrength = sqrt(map.dip[0]*map.dip[0] + map.dip[1]*map.dip[1] + map.dip[2]*map.dip[2]);
			printf("Re-scaled average dipole strength: %6.4f Debye\n", dipstrength);
		} else if( (osc<0) && (map.dipset) && (dipstrength==0.0) ) { 
			real defosc = 0.276*0.276;
			dipstrength = sqrt(map.dip[0]*map.dip[0] + map.dip[1]*map.dip[1] + map.dip[2]*map.dip[2]);
			for(i=0; i<3; i++) map.dip[i] *= (sqrt(defosc)/dipstrength);
			dipstrength = sqrt(map.dip[0]*map.dip[0] + map.dip[1]*map.dip[1] + map.dip[2]*map.dip[2]);
			printf("Assigning average dipole strength: %6.4f Debye\n", dipstrength);
		} else if( (osc<0) && (!map.dipset) && (dipstrength==0.0) ) { 
			dipstrength = get_avg_tc_dip(pbc, x, map, p_pb, bonds, RotMat);
			printf("Using TCC-defined dipole strength: %6.4f Debye\n", dipstrength);

		} else if( (osc>0) && (map.tcset) && (dipstrength==0.0) ) { 
			dipstrength = get_avg_tc_dip(pbc, x, map, p_pb, bonds, RotMat);
			printf("Input dipole strength: %6.4f Debye\n", dipstrength);
			printf("Re-scaling map charges and charge fluxes.\n");
			for(i=0; i<map.ncoupsites; i++) map.CoupSites[i].q *= (sqrt(osc)/dipstrength);
			for(i=0; i<map.ncoupsites; i++) map.CoupSites[i].dq *= (sqrt(osc)/dipstrength);
			dipstrength = get_avg_tc_dip(pbc, x, map, p_pb, bonds, RotMat);
			printf("Re-scaled dipole strength: %6.4f Debye\n", dipstrength);
		}
		// If a dipole moment was specified, use it. 
		if(map.dipset) get_dip(pbc, x, map, p_pb, bonds, elecData, nelec, RotMat, dipData, centData);
		// Otherwise, try to calculate from TCC parameters. 
		else if(map.tcset) get_tc_dip(pbc, x, map, p_pb, bonds, RotMat, dipData, centData);
		// If PDC was specified, calculate PDC couplings. 
		if(!strcmp(map.couptype, "pdc")) if(get_pdc(pbc, p_pb, bonds, dipData, centData, coupData, ncoup)) error = 1;
		// If TCC was specified, calculate TCC couplings. 
		if(!strcmp(map.couptype, "tcc")) if(get_tcc(pbc, map, x, p_pb, bonds, RotMat, coupData, ncoup)) error = 1;
		//if(!strcmp(map.couptype, "tcc")) if(get_tcc_jansen(pbc, map, x, p_pb, bonds, RotMat, coupData, ncoup)) error = 1;
		// If using NNC map, get nearest-neighbor couplings. 
		if(nnc) get_nnc(map, p_pb, bonds, angleData, coupData);
		if(dnnc) get_dnnc(map, p_pb, bonds, angleData, coupData);
		// Write data to files. 
		write_spec_data(specfp, coupData, freqData, dipData, bonds);
		// If saving electrostatic data, write it to file now. 
		if(print_elec) write_elec_dat(elecfp, map, elecData, bonds, nelec);
		// If saving angle data, write it to file now. 
		if(print_angles) write_angle_dat(anglefp, angleData, bonds, nangles);
		// If we ran into an error, don't read any more frames. 
		// Otherwise, continue on toward victory!
		if(error) trrStatus = 0; 
		else trrStatus = read_next_x(oenv, status, &t, natoms, x, box);
		set_pbc(&pbc, ePBC, box); // Re-set pbc using new box
	} while(trrStatus);

	if(anglefilesopen==1) close_angle_files(anglefp);
	if(elecfilesopen==1) close_elec_files(elecfp, map);
	if(infofileopen==1) fclose(infofp);
	if(specfilesopen==1) close_spec_files(specfp, 5);
	if(arraysset==1) unset_arrays(elecData, ProElecData, nelec, angleData, ProAngleData, nangles, freqData, ProFreqData, nfreq, coupData, ncoup, coordData, ProCoordData, dipData, centData, RotMat, ProRotMat, bonds, map.nsites, nPro, promap.nsites);
	if(reparray[0]!=NULL) for(i=0; i<nrep; i++) if(reparray[0][i]!=NULL) free(reparray[0][i]);
	if(reparray[1]!=NULL) for(i=0; i<nrep; i++) if(reparray[1][i]!=NULL) free(reparray[1][i]);
	if(reparray[0]!=NULL) free(reparray[0]);
	if(reparray[1]!=NULL) free(reparray[1]);
	free_resarray(resarray, nres);
	free_amide_map(map);
	for(b=0; b<bonds; b++) free_bond(p_pb[b]);
	if(p_pb!=NULL) free(p_pb);
	return 0;
}





