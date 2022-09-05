// Parse the complementary fragment information from dicer
// Input will be records that look like.
//[3cH]1[2cH]cc[1cH]c1 CHEMBL4553639 [1ClH].[2OH2].Fc1ccc(CN[3CH]=O)cc1 COMP CHEMBL4553639 at.1:3502|2:3502|3:3502
// We want to build a set of files.
// At the top level there are directories for the atom separations
// 6,7,10 would be the directory holding 3 connected linkers where the bond
// distances were 6, 7, and 10 bonds.
// Within each directory, we have individual files that group fragments
// totether by atom types. The file 1044,708,421 contains dicer::DicerFragment
// protos with all the examples of linkers that have these three atom types.

#include <sys/stat.h>
#include <sys/types.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <optional>
#include <vector>

#include "google/protobuf/text_format.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/diced_molecule.pb.h"

namespace dicer_to_cores {

using std::cerr;

constexpr char kCloseBrace = ']';
constexpr char kDirSeparator = '/';

void
Usage(int rc) {
  ::exit(rc);
}

// Fragments are stored in a heirarchy.
// The the top level they are stored by distances between the attachment points.
// Within fragments that have the same attachment points, they are stored by
// atom type of the attachments.
// Within that, they are stored by a mapping from unique smiles to SmilesIdCount.
// This class holds a mapping from unique smiles to SmilesIdCount.
class SameAtomType {
  private:
    IW_STL_Hash_Map<IWString, dicer::SmilesIdCount> _usmi_to_data;
  public:
    int Extra(const IWString& smiles,
              const IWString& id);

    // `fname` is non const because we use null_terminated_chars();
    int Write(IWString& fname) const;
    int Write(IWString_and_File_Descriptor& fname) const;
};

int
SameAtomType::Extra(const IWString& smiles,
                    const IWString& id) {
  auto iter = _usmi_to_data.find(smiles);
  if (iter != _usmi_to_data.end()) {
    auto n = iter->second.n();
    iter->second.set_n(n + 1);
    return 1;
  }

  dicer::SmilesIdCount proto;
  proto.set_smiles(smiles.data(), smiles.length());
  proto.set_id(id.data(), id.length());
  proto.set_n(1);
  auto [_, inserted] = _usmi_to_data.emplace(smiles, std::move(proto));
  if (inserted) {
    return 1;
  }

  cerr << "SameAtomType::Extra:cannot insert '" << smiles << "'\n";
  return 0;
}

int
SameAtomType::Write(IWString& fname) const {
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "SameAtomType::Write:cannot open '" << fname << "'\n";
    return 0;
  }

  return Write(output);
}

int
SameAtomType::Write(IWString_and_File_Descriptor& output) const {
  google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  for (auto & iter : _usmi_to_data) {
    output << iter.first << ' ';

    std::string buffer;
    printer.PrintToString(iter.second, &buffer);
    output << buffer;
    output << '\n';

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

// For each bond separation, 6,8,9, this class is a map from
// different atomtype tuples to examplars.
// The examplars are stored as a map from unique smiles
// to protos that will be written.
// This multi=level map is confusing.
class SameBondSeparation {
  private:
    // The key will be a canonical atom type.
    // The value is a map from unique smiles to a proto.
    IW_STL_Hash_Map<IWString, SameAtomType> _atype_to;

  // private functions.
    int AddNewAtype(const IWString& atype,
                const IWString& smiles,
                const IWString& id);

  public:
    int Extra(const IWString& atype, const IWString& smiles,
              const IWString& id);

    // Within a given directory, we create individual files per atom type.
    int Write(const IWString& dirname) const;
};

int
SameBondSeparation::Extra(const IWString& atype,
                const IWString& smiles,
                const IWString& id) {
  auto iter = _atype_to.find(atype);
  if (iter == _atype_to.end()) {
    return AddNewAtype(atype, smiles, id);
  }

  return iter->second.Extra(smiles, id);
}

// A new atom type combination has been encountered.
// Add it to `_seen`.
int
SameBondSeparation::AddNewAtype(const IWString& atype,
                const IWString& smiles,
                const IWString& id) {
// The value for each atom_type is a map from unique smiles to SmilesIdCount.

  IWString key(atype);
  auto result = _atype_to.emplace(atype, SameAtomType());
  if (! result.second) {
    cerr << "SameBondSeparation:AddNewAtype:cannot add '" << atype << "'\n";
    return 0;
  }

  return result.first->second.Extra(smiles, id);
}

int
SameBondSeparation::Write(const IWString& dirname) const {
  for (auto& iter : _atype_to) {
    IWString fname;
    fname << dirname << '/' << iter.first;
    if (! iter.second.Write(fname)) {
      cerr << "SameBondSeparation::Write:cannot write '" << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

class Options {
  private:
    int _verbose;

    // the number of lines in the input that contain COMP that we examine.
    int _liness_examined;

    // After the line has been processed, now many molecules do we process.
    int _molecules_processed;

    // A mapping from distance string, '4,5,10', to a collection of
    // fragments that have the same bond distance strings.
    IW_STL_Hash_Map<IWString, SameBondSeparation> _distances;

    // The directory into which we write directories.
    IWString _file_name_stem;

  // Private functions.
    int Process(const Molecule& parent,
                 Molecule& fragment,
                 const const_IWSubstring& atype_string);

    int Write(const IWString& distance_string, const SameBondSeparation& value) const;

  public:
    Options();

    int Initialise(Command_Line& cl);

    int Process(const Molecule& parent, const const_IWSubstring& buffer);

    int Write();

    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _molecules_processed = 0;
  _liness_examined = 0;
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (! cl.option_present('S')) {
    cerr << "Must specify the output file stem via the -S option\n";
    return 0;
  }

  if (cl.option_present('S')) {
    cl.value('S', _file_name_stem);
    // Bad if it already exists and is not a directory.
    if (dash_f(_file_name_stem.null_terminated_chars()) &&
        ! dash_d(_file_name_stem.null_terminated_chars())) {
      cerr << "Options::Initialise:the -S option must be a directory '" <<  _file_name_stem << "' invalid\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Data written to '" << _file_name_stem << "'\n";
    }
  }

  return 1;
}

// Parse something that looks like 'ndx:value' and set
// destination[ndx] = value
int
GetAtomType(const const_IWSubstring& buffer,
            extending_resizable_array<uint32_t>& destination) {
  static std::unique_ptr<RE2> rx = std::make_unique<RE2>("^([0-9])\.([0-9]+)$");

  int ndx;
  uint32_t value;
  std::string as_string(buffer.data(), buffer.length());

  if (! rx->FullMatch(as_string, *rx, &ndx, &value)) {
    cerr << "GetAtomType:does not match rx '" << buffer << "'\n";
    return 0;
  }
  destination[ndx] = value;
  return 1;
}

// Parse a token that looks like
// at.1:3502,2:3022,3:5207
// into a vector {0, 3502, 3022, 5207}
// Note local copy of `buffer`.
std::optional<extending_resizable_array<uint32_t>>
AtomTypes(const_IWSubstring buffer) {  
  static IWString at_equals = "at,";

  if (! buffer.starts_with(at_equals)) {
    return 0;
  }

  buffer.remove_leading_chars(3);

  extending_resizable_array<uint32_t> result;

  static constexpr char kComma = ',';

  // Now parse the individual i:j values
  int i = 0;
  const_IWSubstring token;
  while (buffer.nextword(token, i, kComma)) {
    if (! GetAtomType(token, result)) {
      cerr << "AtomTypes:invalid atype specification '" << token << "'\n";
      return 0;
    }
  }

  if (result.empty()) {
    cerr << "AtomTypes:no values\n";
    return std::nullopt;
  }

  return result;
}

// clang-format off
// buffer looks like
// [3cH]1[2cH]cc[1cH]c1 CHEMBL4553639 [1ClH].[2OH2].Fc1ccc(CN[3CH]=O)cc1 COMP CHEMBL4553639 at.1:3502|2:3502|3:3502
// clang-format on
int
Options::Process(const Molecule& parent,
                 const const_IWSubstring& buffer) {
  ++_liness_examined;

  static IWString comp("COMP");

  if (buffer.nwords() != 6) {
    cerr << "Options::Process:wrong token count, got " << buffer.nwords() << " expected 6\n";
    return 0;
  }

  IWString fragment_usmi;

  int i = 0;
  const_IWSubstring token;
  buffer.nextword(fragment_usmi, i);

  Molecule fragment;
  if (! fragment.build_from_smiles(fragment_usmi)) {
    cerr << "Options::Process:Invalid fragment smiles\n";
    return 0;
  }

  // We are not interested in terminal fragments.
  if (fragment.number_isotopic_atoms() < 2) {
    return 1;
  }

  // Next is the name.
  buffer.nextword(token, i);
  if (parent.name() != token) {
    cerr << "Options::Process:name mismatch, parent '" << parent.name() << "' fragment '" << token << "'\n";
    return 0;
  }

  // The complementary smiles.
  buffer.nextword(token, i);

  buffer.nextword(token, i);
  if (token != comp) {
    cerr << "Options::Process:not " << comp << " '" << token << "'\n";
    return 0;
  }

  // The repeated name
  buffer.nextword(token, i);
  if (token != parent.name()) {
    cerr << "Options::Process:second name mismatch, parent '" << parent.name() << " got '" << token << "'\n";
    return 0;
  }

  buffer.nextword(token, i);
  std::optional<extending_resizable_array<uint32_t>> atom_types = AtomTypes(token);
  if (! atom_types) {
    cerr << "optional::Process:cannot parse atom types '" << token << "'\n";
    return 0;
  }

  return Process(parent, fragment, token);
}

int
HighestIsotope(const Molecule& m) {
  const int matoms = m.natoms();
  int result= 0;
  for (int i = 0; i < matoms; ++i) {
    const int iso = m.isotope(i);
    if (iso > result) {
      result = iso;
    }
  }

  return result;
}

Set_of_Atoms
IsotopicAtoms(const Molecule& m) {
  Set_of_Atoms result;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.isotope(i) > 0) {
      result << i;
    }
  }

  return result;
}

// The molecule has isotopic labels at the attachment
// points. Return a canonical representation of the distances
// implied by those atoms. Sorted distances.
IWString
CanonicalDistanceString(Molecule& m) {
  const Set_of_Atoms isotopic_atoms = IsotopicAtoms(m);
  const int niso = isotopic_atoms.number_elements();

  resizable_array<int> distances;

  distances.reserve(niso * (niso - 1) / 2);
  for (int i = 0; i < niso; ++i) {
    atom_number_t a1 = isotopic_atoms[i];
    for (int j = i + 1; j < niso; ++j) {
      const int d = m.bonds_between(a1, isotopic_atoms[j]);
      distances << d;
    }
  }

  std::sort(distances.begin(), distances.end(), [](int d1, int d2) {
    return d1 < d2;
  });
#ifdef DEBUG_CANONICAL_DISTANCE_STRING
  cerr << "Have " << distances.size() << " distances\n";
  for (auto d : distances) {
    cerr << " dist " << d << '\n';
  }
#endif

  IWString result;

  bool needs_separator = false;

  for (int d : distances) {
    if (needs_separator) {
      result << '.';
    } else {
      needs_separator = true;
    }
    result << d;
  }

  return result;
}
// `parent` is the parent molecule
// `fragment` is a generated fragment
// `atype_string` is the at.... string, which is canonical.
// `atom_types` maps isotope to atom type in `fragment`.
// We make the decision to NOT use the fragment unique smiles that
// was read above, just in case...
int
Options::Process(const Molecule& parent,
                 Molecule& fragment,
                 const const_IWSubstring& atype_string) {
  ++_molecules_processed;

  IWString key = CanonicalDistanceString(fragment);
  auto iter = _distances.find(key);
  if (iter != _distances.end()) {
    iter->second.Extra(atype_string, fragment.unique_smiles(), parent.name());
    return 1;
  }

  // New distance combination.
  auto [iter2, inserted] = _distances.try_emplace(key, SameBondSeparation());
  if (! inserted) {
    cerr << "Options::Process:cannot add '" << key << "' to bond separation hash\n";
    return 0;
  }

  return iter2->second.Extra(atype_string, fragment.unique_smiles(), parent.name());
}

int
Options::Write() {

  for (const auto& iter : _distances) {
    if (! Write(iter.first, iter.second)) {
      return 0;
    }
  }

  return 1;
}

// Key is a canonical distance string, '4,7,9'.
// `atype_to_exemplars` is a mapping from canonical atom type
// to a mapping from unique smiles to exemplars.
int
Options::Write(const IWString& distance_string,
               const SameBondSeparation& value) const {
  IWString dirname;
  dirname << _file_name_stem;
  dirname << kDirSeparator << distance_string;

  if (! dash_d(dirname.null_terminated_chars())) {
    int rc = ::mkdir(dirname.null_terminated_chars(), S_IRWXU | S_IRWXG | S_IROTH| S_IXOTH);
    if (rc != 0) {
      cerr << "Options::Write:cannot create bond separation dir '" << dirname << " rc " << rc << '\n';
      return 0;
    }
  }

  return value.Write(dirname);
}

int
Options::Report(std::ostream& output) const {
  output << _liness_examined << " lines of dicer output examined\n";
  output << "Processed " << _molecules_processed << " molecules\n";
  return 1;
}

#ifdef NOLONGEUSER
// Turn the array of distances in `proto` into a canonical
// string having the distances sorted.
IWString
CanonicalDistanceString(const dicer::DicerFragment& proto) {
  std::vector<int> distances = proto.dist();
  std::sort(distances.begin(), distances.end(), [](int d1, int d2) {
    return d1 < d2;
  });

  IWString result;

  bool needs_separator = false;

  for (int d : distances) {
    if (needs_separator) {
      result << '.';
    } else {
      needs_separator = true;
    }
    result << d;
  }

  return result;
}
#endif

// Parse a line that looks like.
// C1(=O)C(=CC(=O)C2=C(NN=C12)CCNC(=O)OCC1=CC=CC=C1)OC CHEMBL590067 B=8

int
InstantiateParent(const const_IWSubstring& buffer,
                  Molecule& m) {
  if (! m.build_from_smiles(buffer)) {
    return 0;
  }

  // Remove the B= from name.

  IWString mname = m.name();
  mname.strip_trailing_blanks();
  for (int i = mname.length() - 1; i >= 0; --i) {
    if (mname[i] != ' ') {
      continue;
    }
    mname.iwtruncate(i);
    m.set_name(mname);
    return 1;
  }
  
  cerr << "InstantiateParent:huh this cannot happen\n";
  return 0;
}

int
DicerToCores(iwstring_data_source& input,
             Options& options,
             IWString_and_File_Descriptor& output) {
  std::unique_ptr<RE2> b_equals = std::make_unique<RE2>(" B=[0-9][0-9]*$");

  IWString comp = " COMP ";

  const_IWSubstring buffer;
  Molecule parent;
  while (input.next_record(buffer)) {
    if (iwre2::RE2PartialMatch(buffer, *b_equals)) {
      if (! InstantiateParent(buffer, parent)) {
        cerr << "Cannot parse parent smiles '" << buffer << "'\n";
        return 0;
      }
      continue;
    }

    if (! buffer.contains(comp)) {
      continue;
    }

    if (! options.Process(parent, buffer)) {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
DicerToCores(const char* fname,
             Options& options,
             IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "DicerToCores:cannot open '" << fname << "'\n";
    return 0;
  }

  return DicerToCores(input, options, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vS:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  Options options;
  if (! options.Initialise(cl)) {
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char* fname : cl) {
    if (! DicerToCores(fname, options, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  options.Write();

  return 0;
}
}  // namespace dicer_to_cores

int
main(int argc, char** argv) {
  int rc = dicer_to_cores::Main(argc, argv);

  return rc;
}
