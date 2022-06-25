all: main.cpp bed_structures.h bed.cpp bed.h core.cpp core.h edit_dist.h edit_dist.cpp edlib/edlib/src/edlib.cpp edlib/edlib/include/edlib.h
	g++ -Wall -o saber main.cpp bed_structures.h bed.cpp bed.h core.cpp core.h edit_dist.cpp edit_dist.h edlib/edlib/src/edlib.cpp -I edlib/edlib/include
	

sim: rearrangement_sim.cpp bed_structures.h bed.cpp bed.h core.cpp core.h edit_dist.cpp edit_dist.h 
	g++ -Wall -o rearrangement_sim rearrangement_sim.cpp bed.cpp bed.h edit_dist.cpp edit_dist.h bed_structures.h core.cpp core.h edlib/edlib/src/edlib.cpp -I edlib/include
	
