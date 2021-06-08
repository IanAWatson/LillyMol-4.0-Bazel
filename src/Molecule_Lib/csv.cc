// I/O to smiles CSV form.

#include "istream_and_type.h"
#include "molecule.h"
#include "rwmolecule.h"

static char csv_separator = ',';

int
Molecule::write_molecule_csv(std::ostream& output) {
  output << smiles() << csv_separator << name() << '\n';

  return output.good();
}

int
Molecule::read_molecule_csv_ds (iwstring_data_source & input) {
  if (input.eof())
    return 0;

  const_IWSubstring buffer;
  EXTRA_STRING_RECORD(input, buffer, "read mol csv");

  const_IWSubstring token;
  int i = 0;
  if (! buffer.nextword(token, i, csv_separator)) {
    return 0;
  }

  if (! build_from_smiles(token)) {
    cerr << "Cannot parse smiles '" << token << "'\n";
    return 0;
  }

  if (! buffer.nextword(token, i, csv_separator)) {
    return 1;
  }

  set_name(token);

  return 1;
}
