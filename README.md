# g_amide

OVERVIEW

The program g_amide is a post-processing tool for use with the GROMCAS 
(www.gromacs.org) molecular dynamics (MD) simulation package. The program 
allows the user to predict Amide I vibrational frequencies and coupling 
parameters for individual peptide amide groups as a function of time 
along an MD trajectory. In conjunction with the sister program g_spec (in 
preparation for release on github), the program allows the user to predict
Amide I vibrational spectra for different protein structures, allowing the 
user to link MD simulation ensembles with experimental spectroscopic data. 

g_amide uses a set of electrostatic- and dihedral-based maps for frequency
and coupling constant calculations. For an overview of these methods see
[Reppert and Tokmakoff, J Chem Phys. 142, 125104 (2015)]. In works that 
make use of the g_amide source code, we ask that you cite both the relevant
spectroscopic map (see README.txt file in the maps subdirectory) and the 
following review article: 

"Computational Amide I 2D IR Spectroscopy as a Probe of Protein Structure 
and Dynamics". Ann. Rev. Phys. Chem. Vol. 67 (2016), in press. 
DOI: 10.1146/annurev-physchem-040215-112055

The program calls GROMACS libraries for file parsing and so requires a local
gromacs installation to be accessible (see www.gromacs.org). Compiling the code
requires access to the gromacs home installation directory including the
gromacs lib/ and include/gromacs directories. The lib/directory should contain
a libgmx.so file which links to the installed gromacs libraries, while the
include directory should contain a variety of *.h header files. In a default 
gromacs installation, the home directory is just /usr/local/gromacs and
contains the lib/ and include/gromacs/ directories. 


INSTALLATION

To install the program: 

(0) Choose an installation directory.

This may be any directory where you have read/write permissions, but should be
a permanent installation location (e.g. not your user download folder). For
example, if your username is mike, you might wish to install the program in an 
"apps" directory such as /home/mike/apps/. Move the g_amide directory (where
this README.txt file is located) into the chosen installation directory and at
the command prompt cd to that location, e.g. in our example 

	cd /home/mike/apps/g_amide/


(1) Prepare the Makefile. 

Next, modify the Makefile provided in the src directory. E.g. to open the file
using vi, type at the console 

	vi src/Makefile

Modifications should be required for only the first four lines. Default settings 
are pasted below. 

	CC=gcc
	OMP_PARALLEL=FALSE
	LIB_DIR=/usr/local/gromacs/lib
	INC_DIRS=/usr/local/gromacs/include/gromacs/

The first line specifies a C compiler. The gcc compiler should be available on
most linux systems; in this case no change is necessary. 

The second line specifies whether or not parallelization should be enabled.
g_amide has the ability to use OpenMP libraries to parallelize computation
(i.e. parameters are calculated for multiple bonds simultaneously using
multiple cores). To use this capability, a working installation of OpenMP
must be available, and the flag OMP_PARALLEL must be set to TRUE. The default
setting of FALSE does not require the presence of OpenMP libraries, but will
not carry out calculations in parallel.  

The final two lines specify the gromacs lib and include directories already
discussed. In general, these lines should read 

	LIB_DIR=<gromacs home>/lib
	INC_DIRS=<gromacs home>/include/gromacs/

where <gromacs home> specifies the home directory of your local gromacs
installation. 


(2) Compile the code. 

Once the Makefile has been modified, type at the command line (still in the 
g_amide directory)

	make -C src/

to compile the code located in the g_amide/src directory. This should produce 
a binary file called g_amide, which can be called from the command line. 



(3) Add the binary file to the user path. 

To make this program accessible outside the installation directory, it should
be added to the user path. This can be achieved either by copying the file to
a directory already in the path (e.g. /home/<username>/bin) or by adding the
g_amide/src directory to the user path. If modifying the path, you may wish to
do so in a startup file which is executed every time a user logs onto the
system. For example, for the example above, one might append the line

	PATH=$PATH:/home/mike/apps/g_amide/src/

to the file ~/.bashrc, so that the program src directory will be added to the 
path every time the user logs onto the system. 


(4) Run a test calculation!

To run a test calculation, follow the instructions in TUTORIAL.txt in the 
g_amide/test directory. 





