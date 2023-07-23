#ifndef PATH_AROUND_RING_H
#define PATH_AROUND_RING_H

class Molecule;
class Set_of_Atoms;

extern int path_around_edge_of_ring_system (Molecule & m,
                                            const int * process_these_atoms,
                                            int flag,
                                            Set_of_Atoms & s);
#endif
