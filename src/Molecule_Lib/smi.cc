#include <iostream>
#include <memory>

// Functions for reading and writing smiles.

#include "Foundational/data_source/iwstring_data_source.h"

#include "molecule.h"
#include "Molecule_Lib/moleculeio.h"
#include "Molecule_Lib/rwmolecule.h"

using std::cerr;

int
Molecule::read_molecule_smi_ds(iwstring_data_source & input)
{
  assert(input.good());

  if (input.eof())
    return 0;

  IWString buffer;

  EXTRA_STRING_RECORD(input, buffer, "read mol smi");

  if (buffer.length() < 1) {
    return 0;
  }

  if (! build_from_smiles(buffer.rawchars(), buffer.length())) {
    cerr << "Molecule::read_molecule_smi_ds: Cannot interpret smiles\n";
    cerr << buffer << '\n';
    return 0;
  }

  if (unconnect_covalently_bonded_non_organics_on_read()) {
    _do_unconnect_covalently_bonded_non_organics();
  }

  return 1;
}

int
Molecule::write_molecule_usmi (std::ostream & os, const IWString & comment)
{
  assert(ok());
  assert(os.good());

  os << unique_smiles();

  if (comment.length())
    os << " " << comment;
  
  os << moleculeio::newline_string();

  if (moleculeio::flush_files_after_writing_each_molecule()) {
    os.flush();
  }

  return os.good();
}

int
Molecule::write_molecule_nausmi (std::ostream & os, const IWString & comment)
{
  assert(ok());
  assert(os.good());

  if (UNIQUE_SMILES_ORDER_TYPE != _smiles_information.smiles_order_type())
    invalidate_smiles();

  os << non_aromatic_unique_smiles();

  if (comment.length())
    os << " " << comment;
  
  os << moleculeio::newline_string();

  if (moleculeio::flush_files_after_writing_each_molecule())
    os.flush();

  return os.good();
}

int
Molecule::write_molecule_smi (std::ostream & os, const IWString & comment)
{
  assert(ok());
  assert(os.good());

  os << smiles();

  if (comment.length())
    os << ' ' << comment;
  
  os << moleculeio::newline_string();

  if (moleculeio::flush_files_after_writing_each_molecule())
    os.flush();

  return os.good();
}

int
Molecule::write_molecule_rsmi(std::ostream & os, const IWString & comment)
{
  assert(ok());
  assert(os.good());

  os << random_smiles();

  if (comment.length())
    os << " " << comment;
  
  os << moleculeio::newline_string();

  if (moleculeio::flush_files_after_writing_each_molecule())
    os.flush();

  return os.good();
}
