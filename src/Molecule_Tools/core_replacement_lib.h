#ifndef MOLECULE_TOOLS_CORE_REPLACEMENT_LIB_H
#define MOLECULE_TOOLS_CORE_REPLACEMENT_LIB_H

#include "Molecule_Lib/molecule.h"

namespace separated_atoms {
int
IdentifyAtomsToRemove(Molecule& m,
                      const Set_of_Atoms& matched_atoms,
                      int * to_be_removed);
}  // namespace separated_atoms

#endif  // MOLECULE_TOOLS_CORE_REPLACEMENT_LIB_H
