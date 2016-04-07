
ELECTROSTATIC FREQUENCY SHIFTS

Currently seven options are included in this distribution in the folder 
g_amide/maps. File names and references (for the frequency shift map) 
are as follows:

	File Name	Reference	Recommended force field charges
	DC.txt		[1]		CHARMM27
	DO.txt		[1]		OPLS-AA
	JO.txt		[2]		OPLS-AA
	SG.txt		[3]		GROMOS
	4PN-4.txt	[4]		CHARMM27
	4PN-150.txt	[4]		CHARMM27
	DC15.txt	[4]		CHARMM27

The first four are those identified as DC, DO, JO, and SG in [J. Chem. Phys.
142, 125104 (2015)]. 

Note that parameters for the DC15 map are not actually reported in Ref.[4]. 
This reference does, however, describe the data set on which the map is
parameterized. The frequency-shift parameters used here come from the best-
fit line shown in the supporting information of the paper, essentially a 
modification of the DO map of Ref. [1]. 

[1] Reppert and Tokmakoff, J. Chem. Phys. 138, 134116 (2013)
[2] Jansen and Knoester, J. Chem. Phys. 124, 044502 (2006)
[3] Wang et al. J. Phys. Chem. B, 115, 3713 (2011)
[4] Reppert and Tokmakoff, J. Chem. Phys. 143, 061102 (2015)




ELECTROSTATIC COUPLING

Maps DC, DO, JO, 4PN-4, 4PN-150, and DC15  use the transition charge coupling 
map of [J. Chem. Phys. 125, 044312 (2006)]. 

Map SG uses the transition dipole coupling model of Torii and Tasumi in 
[J. Raman. Spectr. 29, 81 (1998)]. 




TRANSITION DIPOLE VECTOR

Map JO uses the electrostatic model of [J. Chem. Phys. 125, 044312 (2006)] 
in which the transition dipole moment vector depends on the local electric
field. 

Maps DC, DO, JO, 4PN-4, 4PN-150, and DC15 use the zero-field vector of the 
JO map, i.e. a fixed (molecular-frame) transition dipole vector equal to 
the value of the JO map transition dipole vector in vacuum. 

Map SG uses the transition dipole moment of Torii and Tasumii defined in 
[J. Raman. Spectr. 29, 81 (1998)]. 



NEAREST-NEIGHBOR EFFECTS

All maps use the nearest-neighbor coupling map of [J. Chem. Phys. 125, 
044312 (2006)]. 

In addition, the JO and SG maps apply a dihedral-dependent nearest-neighbor 
frequency shift defined by the nearest-neighbor shift map of [J. Chem. Phys. 
125, 044312 (2006)]. 










