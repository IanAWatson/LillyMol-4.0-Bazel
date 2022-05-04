// Reaction enumeration via smiles ring closures

#include <iostream>
#include <limits>
#include <vector>

#include "re2/re2.h"

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/ring_closure_enumeration_lib.h"

namespace ring_closure_enumeration {

using std::cerr;

void
Usage(int rc) {

  cerr << " -r <bond>         specify bond making. <component>:<isotope><bond><ring>, eg 0:3=99\n";
  cerr << " -x                remove chirality\n";
  cerr << " -l                strip to largest fragment\n";
  cerr << " -v                verbose output\n";

  ::exit(rc);
}

class Enumeration {
  private:
    int _verbose;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    FileType _input_type;

    Chemical_Standardisation _chemical_standardisation;

    IWString _component_separator;

    Reagents * _reagents;
    int _number_reagents;

    resizable_array_p<BondFormation> _bond_formation;

    int _min_natoms;
    int _max_natoms;

    resizable_array_p<Substructure_Query> _must_contain;

    int _molecules_formed;
    int _molecules_written;

  // Private functions

    int ReadReagents(const char * fname, Reagents& reagents);
    int ReadReagents(data_source_and_type<Molecule>& input, Reagents& reagents);
    int Enumerate1(IWString_and_File_Descriptor& output);
    int Enumerate2(IWString_and_File_Descriptor& output);
    int Enumerate3(IWString_and_File_Descriptor& output);

    int SetName(Molecule& m,
                     const IWString& name1,
                     const IWString& name2) const;
    int MaybeWrite(Molecule& m, IWString_and_File_Descriptor& output);
    int MatchesMustHaveQuery(Molecule& m);

  public:
    Enumeration();
    ~Enumeration();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int Report(std::ostream& output) const;

    FileType input_type() const {
      return _input_type;
    }

    int Enumerate(IWString_and_File_Descriptor& output);
};

Enumeration::Enumeration() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _input_type = FILE_TYPE_INVALID;
  _number_reagents = 0;
  _reagents = nullptr;
  _component_separator = " + ";
  _min_natoms = 0;
  _max_natoms = std::numeric_limits<int>::max();

  _molecules_formed = 0;
  _molecules_written = 0;
}

Enumeration::~Enumeration() {
  if (_reagents != nullptr) {
    delete [] _reagents;
  }
}

int
Enumeration::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('x')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove chirality from input molecules\n";
    }
  }

  if (cl.option_present('g')) {
    if (! _chemical_standardisation.construct_from_command_line(cl, _verbose, 'g')) {
      cerr << "Cannot initialise chemical standardisation\n";
      return 0;
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce molecules to largest fragment\n";
    }
  }

  if (! cl.option_present('r')) {
    cerr << "Must specify one or more reactions via the -r option\n";
    return 0;
  }

  const_IWSubstring transformation;
  for (int i = 0; cl.value('r', transformation, i); ++i) {
    std::unique_ptr<BondFormation> bond = std::make_unique<BondFormation>();
    if (! bond->Build(transformation)) {
      cerr << "Enumeration::initialise:cannot initialis bond formation directive " << transformation << '\n';
      return 0;
    }
    _bond_formation << bond.release();
  }

  if (cl.option_present('c')) {
    if (! cl.value('c', _min_natoms) || _min_natoms < 1) {
      cerr << "Enumeration::Initialise:invalid min natoms (-c)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will discard products with fewer than " << _min_natoms << " atoms\n";
    }
  }

  if (cl.option_present('C')) {
    if (! cl.value('C', _max_natoms) || _max_natoms < _min_natoms) {
      cerr << "Enumeration::Initialise:invalid max natoms (-C)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will discard products with more than " << _max_natoms << " atoms\n";
    }
  }

  if (cl.option_present('q')) {
    if (! process_queries(cl, _must_contain, _verbose, 'q')) {
      cerr << "Enumeration::Initialise:cannot read must contain queries (-q)\n";
      return 0;
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring smarts;
    for (int i = 0; cl.value('s', smarts, i); ++i) {
      std::unique_ptr<Substructure_Query> qry;
      if (! qry->create_from_smarts(smarts)) {
        cerr << "Enumeration::Initialise:cannot parse smarts '" << smarts << "'\n";
        return 0;
      }
      _must_contain << qry.release();
    }
  }

  if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) { // reading a pipe, assume smiles
    _input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern all file types, use the -i option\n";
    return 0;
  } else if (!process_input_type(cl, _input_type)) {
    return 0;
  }

  _number_reagents = cl.number_elements();
  _reagents = new Reagents[_number_reagents];
  for (int i = 0; i < _number_reagents; ++i) {
    if (! ReadReagents(cl[i], _reagents[i])) {
      cerr << "Enumeration::Initialise:cannot read '" << cl[i] << "'\n";
      return 0;
    }
  }

  for (const BondFormation* bond : _bond_formation) {
    if (bond->component() >= _number_reagents) {
      cerr << "Enumeration::Initialise:invalid component on BondFormation " << bond->component() << '\n';
      return 0;
    }
  }

  for (int i = 0; i < _number_reagents; ++i) {
    resizable_array<BondFormation*> ptrs;
    for (BondFormation* bond : _bond_formation) {
      if (bond->component() == i) {
        ptrs << bond;
      }
    }
    if (ptrs.empty()) {
      cerr << "Enumeration::Initialise:no BondFormation for reagent " << i << '\n';
      return 0;
    }

    if (! _reagents[i].CreateFragments(ptrs)) {
      cerr << "Enumeration::Initialise:reagent " << i << " cannot create fragments, have " <<
               ptrs.size() << " bond formation directives\n";
      return 0;
    }
  }

  return 1;
}

int
Enumeration::ReadReagents(const char * fname,
                          Reagents& reagents) {
  data_source_and_type<Molecule> input(_input_type, fname);
  if (! input.good()) {
    cerr << "Enumeration::ReadReagents:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadReagents(input, reagents);
}

int
Enumeration::ReadReagents(data_source_and_type<Molecule>& input,
                          Reagents& reagents) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    if (! Preprocess(*m)) {
      delete m;
      return 0;
    }

    reagents.Add(m);
  }

  return 1;
}

int
Enumeration::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (m.natoms() == 0) {
    return 0;
  }

  return 1;
}

int
Enumeration::Report(std::ostream& output) const {
  output << _molecules_formed << " molecules formed\n";
  output << _molecules_written << " molecules written\n";

  return 1;
}

int
Enumeration::Enumerate(IWString_and_File_Descriptor& output) {
  if (_number_reagents == 1) {
    return Enumerate1(output);
  }
  if (_number_reagents == 2) {
    return Enumerate2(output);
  }
  if (_number_reagents == 3) {
    return Enumerate3(output);
  }

  cerr << "Enumeration::Enumerate:do not know how to enumerate " << _number_reagents << " components\n";
  return 0;
}

int
Enumeration::Enumerate1(IWString_and_File_Descriptor& output) {
  int rc = 0;
  for (const IWString* smiles : _reagents[0].smiles()) {
    Molecule m;
    if (! m.build_from_smiles(*smiles)) {
      cerr << "Enumeration::Enumerate1:cannot parse " << smiles << '\n';
      return 0;
    }
    MaybeWrite(m, output);
    ++rc;
  }

  return rc;
}

int
Enumeration::SetName(Molecule& m,
                     const IWString& name1,
                     const IWString& name2) const {
  IWString name;
  name << name1 << _component_separator << name2;
  m.set_name(name);
  return 1;
}

int
Enumeration::Enumerate2(IWString_and_File_Descriptor& output) {
  int rc = 0;
  const int n1 = _reagents[0].number_reagents();
  const int n2 = _reagents[1].number_reagents();

  for (int i = 0; i < n1; ++i) {
    const IWString* smi1 = _reagents[0].smiles()[i];
    const IWString& name1 = _reagents[0].MolName(i);
    for (int j = 0; j < n2; ++j) {
      const IWString* smi2 = _reagents[1].smiles()[j];
      IWString smiles;
      smiles << *smi1 << *smi2;
      Molecule m;
      if (! m.build_from_smiles(smiles)) {
        cerr << "Enumeration::Enumerate2:cannot parse " << smiles << '\n';
        return 0;
      }
      if (m.smiles().length() < 1000) {
        continue;
      }
      SetName(m, name1, _reagents[1].MolName(j));

      MaybeWrite(m, output);
      ++rc;
    }
  }

  return rc;
}

int
Enumeration::Enumerate3(IWString_and_File_Descriptor& output) {
  int rc = 0;
  return rc;
}

int
Enumeration::MatchesMustHaveQuery(Molecule& m) {
  if (_must_contain.empty()) {
    return 1;
  }

  Molecule_to_Match target(&m);
  for (Substructure_Query * q : _must_contain) {
    if (q->substructure_search(target)) {
      return 1;
    }
  }

  return 0;
}

int
Enumeration::MaybeWrite(Molecule& m,
                        IWString_and_File_Descriptor& output) {
  ++_molecules_formed;

  const int matoms = m.natoms();
  if (matoms < _min_natoms) {
    return 1;
  } else if (matoms > _max_natoms) {
    return 1;
  }

  if (! MatchesMustHaveQuery(m)) {
    return 1;
  }

  ++_molecules_written;

  output << m.smiles();
  output << ' ';
  output << m.name();
  output << '\n';
  output.write_if_buffer_holds_more_than(32768);
  return 1;
}

int
RingClosureEnumeration(Enumeration& enumeration,
            Molecule& m,
            IWString_and_File_Descriptor& output) {
  return 1;
}

int
RingClosureEnumeration(Enumeration& enumeration,
            data_source_and_type<Molecule>& input,
            IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! enumeration.Preprocess(*m)) {
      return 0;
    }

    if (! RingClosureEnumeration(enumeration, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
RingClosureEnumeration(Enumeration& enumeration,
            const char * fname,
            IWString_and_File_Descriptor& output) {
  FileType input_type = enumeration.input_type();
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return RingClosureEnumeration(enumeration, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:g:lcr:c:C:xq:s:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (! process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process aromaticity options\n";
    return 1;
  }

  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process standard elements options (-E)\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  Enumeration enumeration;
  if (! enumeration.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  enumeration.Enumerate(output);

  if (verbose) {
    enumeration.Report(cerr);
  }

  return 0;
}

}  // namespace ring_closure_enumeration

int
main(int argc, char** argv) {
  int rc = ring_closure_enumeration::Main(argc, argv);

  return rc;
}
