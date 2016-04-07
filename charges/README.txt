
CHARGEFILES

Chargefiles allow the user to apply different atomic charges during frequency
shift calculations than are actually applied during the MD run. This is
extremely important when comparing data between force fields or, more 
specifically, when analyzing a simulation run with one force field using a 
spectroscopic map developed against another force field. 

For example, when analyzing an OPLS-AA simulation when using the DC or DC15 
maps, use of a CHARMM27 charge file is essential to achieving accurate results.

Even for CHARMM27 simulations, we recommend that the modified CHARMM27 map 
ch27g.txt be applied when using the DC or DC15 maps. The ch27g charge file 
includes modified charges for glycine residues that (as described in [J. Chem. 
Phys. 143, 061102 (2015)]) give better agreement with experiment. 

The charge files included in this distribution are listed below. 

	ch27_spce.txt  
	oplsaa_spce.txt
	ch27g.txt  

The first two specify native charges for the CHARMM27 adn OPLS-AA force fields. 
These charge files are useful when running simulations with a different force 
field and using a frequency shift map based on either CHARMM27 or OPLS-AA 
charges (e.g. the DC map might be used with charmm27_spce.txt, and the DO map 
might be used with oplsaa_spce.txt). In both files the default water charges
are those assigned by the SPC/E water model. These are the charges originally
applied to the dipeptide data set of [J. Chem. Phys. 138, 134116 (2013)]. 

The last file specifies the CHARMM27 charge file with modified glycine charges
as described in [J. Chem. Phys. 143, 061102 (2015)]. Water model charges are
set to those of the TIP3P water model. 



