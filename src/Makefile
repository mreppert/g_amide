CC=gcc
OMP_PARALLEL=TRUE
LIB_DIR=/software/gromacs-4.6-el6-x86_64+intel-12.1/lib
INC_DIRS=/software/gromacs-4.6-el6-x86_64+intel-12.1/include/gromacs/


ifeq ($(OMP_PARALLEL), TRUE) 
CFLAGS=-O3 -Wall -fopenmp -DOMP_PARALLEL=TRUE
LIBFLAGS=-lgmx -ldl -lm -pthread -lgomp
else
CFLAGS=-O3 -Wall
LIBFLAGS=-lgmx -ldl -lm 
endif

OBJECTS=gmx_amide.o gmx_amide_map.o gmx_amide_bond.o gmx_amide_charge.o mem_helper.o

all: $(OBJECTS)
	$(CC) $(OBJECTS) -L $(LIB_DIR) $(LIBFLAGS) -o g_amide

gmx_amide.o: gmx_amide.c
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) gmx_amide.c

gmx_amide_map.o: gmx_amide_map.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) gmx_amide_map.c

gmx_amide_bond.o: gmx_amide_bond.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) gmx_amide_bond.c

gmx_amide_charge.o: gmx_amide_charge.c
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) gmx_amide_charge.c

mem_helper.o: mem_helper.c 
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) mem_helper.c

clean:
	rm *.o
