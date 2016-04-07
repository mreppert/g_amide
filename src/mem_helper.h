#ifndef MEM_HELPER_H
#define MEM_HELPER_H

#include "statutil.h"

int allocate_2d_array_real(real*** p_Array, int dim1, int dim2);
int allocate_3d_array_real(real**** p_Array, int dim1, int dim2, int dim3);
int free_2d_array_real(real** Array, int dim1, int dim2);
int free_3d_array_real(real ***Array, int dim1, int dim2, int dim3);

#endif
