/* mem_helper.c contains various subroutines useful for memory 
 * allocation purposes. 
 *
 * Written by Mike Reppert (2015)
 */

#include "mem_helper.h"

int allocate_2d_array_real(real*** p_Array, int dim1, int dim2) {
	int i;
	if(dim1*dim2==0) return 0;
	real **Array = (real**) malloc(dim1*sizeof(real*));
	if((Array)==NULL && dim1!=0) {
		return 0;
	} else {
		for(i=0; i<dim1; i++) {
			Array[i] = (real*) malloc(dim2*sizeof(real));
			if(Array[i]==NULL && dim2!=0) {
				int k;
				for(k=0; k<i; k++) free(Array[k]);
				free(Array);
				return 0;
			}
		}
	}
	*p_Array = Array;
	return 1;
}

int allocate_3d_array_real(real**** p_Array, int dim1, int dim2, int dim3) {
	int i,j;
	if(dim1*dim2*dim3==0) return 0;
	real ***Array = (real***) malloc(dim1*sizeof(real**));
	if(Array==NULL) {
		return 0;
	} else {
		for(i=0; i<dim1; i++) {
			Array[i] = (real**) malloc(dim2*sizeof(real*));
			if(Array[i]==NULL) {
				int ix;
				for(ix=0; ix<i; ix++) free(Array[ix]);
				free(Array);
				return 0;
			} else {
				for(j=0; j<dim2; j++) {
					Array[i][j] = (real*) malloc(dim3*sizeof(real));
					if(Array[i][j]==NULL) {
						int ix,jx;
						for(jx=0; jx<j; jx++) free(Array[i][jx]);
						free(Array[i]);
						for(ix=0; ix<i; ix++) {
							for(jx=0; jx<dim2; jx++) free(Array[ix][jx]);
							free(Array[ix]);
						}
						free(Array);
						return 0;
					}
				}
			}
		}
	}
	*p_Array = Array;
	return 1;
}

int free_2d_array_real(real** Array, int dim1, int dim2) {
	int i;
	if(dim1*dim2==0) return 1;
	for(i=0; i<dim1; i++) {
		free(Array[i]);
	}
	free(Array);
	return 1;
}

int free_3d_array_real(real ***Array, int dim1, int dim2, int dim3) {
	int i,j;
	if(dim1*dim2*dim3==0) return 1;
	for(i=0; i<dim1; i++) {
		for(j=0; j<dim2; j++) {
			free(Array[i][j]);
		}
		free(Array[i]);
	}
	free(Array);
	return 1;
}

