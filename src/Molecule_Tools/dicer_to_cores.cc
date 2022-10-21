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
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "google/protobuf/text_format.h"

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/diced_molecule.pb.h"

namespace dicer_to_cores {

using std::cerr;

constexpr char kCloseBrace = ']';
constexpr char kDirSeparator = '/';

// Being lazy, file scope variables.

int remove_isotopes_from_products = 1;

IWString_and_File_Descriptor stream_for_transformations;

int products_generated = 0;

void
Usage(int rc) {
  ::exit(rc);
}

#ifdef NOT_USED_NOW
// We need to pass various pieces of information
// down to the classes that do the work.
struct Args {
  Molecule* parent;

  // Isotopes can be removed from products.
  int unset_isotopes;

  // It can be helpful to be able to view the parent, the
  // fragments, the replacement and product.
  IWString_and_File_Descriptor stream_for_transformations;

  int products_generated;

  public:
    Args();

    void set_parent(Molecule& m) {
      parent = &m;
    }

    int MaybeWriteDetails(Molecule& fragment, Molecule& complement,
                          const IWString& replacement_smiles, Molecule& replacement);
};

Args::Args() {
  parent = nullptr;
  unset_isotopes = 0;
  products_generated = 0;
}
#endif

// Increment some namespace scope variables, and if the stream
// for detailed information is active, write to it.
int
MaybeWriteDetails(Molecule& parent, Molecule& fragment, Molecule& complement,
                  const IWString& replacement_smiles, Molecule& replacement) {
  ++products_generated;
  if (! stream_for_transformations.active()) {
    return 1;
  }

  constexpr char kSep = ' ';

  stream_for_transformations << parent.smiles() << kSep << parent.name() << kSep << "parent\n";
  stream_for_transformations << fragment.smiles() << kSep << "frag\n";
  stream_for_transformations << complement.smiles() << kSep << "complement\n";
  stream_for_transformations << replacement_smiles << kSep << "replacement\n";
  stream_for_transformations << replacement.smiles() << kSep << "product\n";

  return 1;
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

  // private functions.
    int Add(const const_IWSubstring& buffer);
  public:
    int Extra(const IWString& smiles,
              const IWString& id);

    int NumberFragments() const {
      return _usmi_to_data.size();
    }


    int Read(const char* fname);
    int Read(iwstring_data_source& input);

    // `fname` is non const because we use null_terminated_chars();
    int Write(IWString& fname) const;
    int Write(IWString_and_File_Descriptor& fname) const;

    int ReplaceFragments(Molecule& parent,
                         Molecule& fragment,
                         const const_IWSubstring& complementary_smiles,
                         IWString_and_File_Descriptor& output);
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

int
SameAtomType::Read(const char* fname) {
  iwstring_data_source input;
  if (! input.open(fname)) {
    cerr << "SameAtomType::Read:cannot open '" << fname << "'\n";
    return 0;
  }
  return Read(input);
}

int
SameAtomType::Read(iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! Add(buffer)) {
      cerr << "SameAtomType::Read:cannot process '" << buffer << "'\n";
      return 0;
    }
  }

  return _usmi_to_data.size();
}

// Return a map from atom number to isotope for those
// atoms that have an isotope.
std::unordered_map<int, int>
GetIsotopicAtoms(const Molecule& m) {
  std::unordered_map<int, int> result;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    const int iso = m.isotope(i);
    if (iso) {
      result[iso] = i;
    }
  }

  return result;
}

int
CompatibleIsotopes(const std::unordered_map<int, int>& iso1, const std::unordered_map<int, int>& iso2) {
  if (iso1.size() != iso2.size()) {
    return 0;
  }

  return 1;
}

int
SameAtomType::ReplaceFragments(Molecule& parent,
                         Molecule& fragment,
                         const const_IWSubstring& complementary_smiles,
                         IWString_and_File_Descriptor& output) {
  constexpr char kSep = ' ';

  Molecule complement;
  if (! complement.build_from_smiles(complementary_smiles)) {
    cerr << "SameAtomType::ReplaceFragments:invalid complementary smiles " << complementary_smiles << '\n';
    return 0;
  }

  const std::unordered_map<int, int> complement_isotopes = GetIsotopicAtoms(complement);

  for (const auto [usmi, value] : _usmi_to_data) {
    if (usmi == parent.unique_smiles()) {
      continue;
    }

    Molecule replacement;
    if (! replacement.build_from_smiles(value.smiles())) {
      cerr << "SameAtomType::ReplaceFragments:cannot parse db smiles '" << value.smiles() << "'\n";
     continue;
    }

    const std::unordered_map<int, int> replacement_isotopes = GetIsotopicAtoms(replacement);
    if (! CompatibleIsotopes(replacement_isotopes, complement_isotopes)) {
      cerr << "SameAtomType::ReplaceFragments:isotope count mismatch\n";
      continue;
    }

    const int initial_natoms = replacement.natoms();

    replacement.add_molecule(&complement);

    for (auto [iso, atom] : replacement_isotopes) {
      //cerr << "isotope " << iso << " on atom " << atom << '\n';
      auto iter = complement_isotopes.find(iso);
      if (iter == complement_isotopes.end()) {
        cerr << "SameAtomType::ReplaceFragments:no isotope " << iso << " in complement\n";
        return 0;
      }
      replacement.set_implicit_hydrogens_known(atom, 0);
      int a2 = initial_natoms + iter->second;
      replacement.set_implicit_hydrogens_known(a2, 0);

      replacement.add_bond(atom, a2, SINGLE_BOND);
    }

    MaybeWriteDetails(parent, fragment, complement, value.smiles(), replacement);
    if (remove_isotopes_from_products) {
      replacement.transform_to_non_isotopic_form();
    }

    output << parent.smiles() << kSep << replacement.smiles() <<
              kSep << parent.name() << kSep << "%%" << value.id() << "%%\n";

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

// Parse a record that looks like
// O=[3CH]CN1CC(c2[1cH]cccc2)[2CH2]C1 smiles: "O=[3CH]CN1CC(c2[1cH]cccc2)[2CH2]C1" id: "CHEMBL3484615" n: 1 
// The first token is the key, and the rest is a textproto.
int
SameAtomType::Add(const const_IWSubstring& buffer) {
  int i = buffer.index(' ');
  IWString key;
  key.set(buffer.data(), i);

  std::string txtproto(buffer.data() + i + 1, buffer.length() - i - 2);
  // cerr << "Txtproto '" << txtproto << "'\n";
  dicer::SmilesIdCount proto;
  if (! google::protobuf::TextFormat::ParseFromString(txtproto, &proto)) {
    cerr << "SameAtomType::Add:cannot parse '" << txtproto << "'\n";
    return 0;
  }

  auto [_, inserted] = _usmi_to_data.try_emplace(key, std::move(proto));
  if (! inserted) {
    cerr << "SameAtomType::Add:cannot insert '" << key << "'\n";
    return 0;
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

    // When we are doing core replacements, we can either respect the atom
    // types, or just replace everything with a plausible geometry.
    int _must_match_atom_types;

  // private functions.
    int AddNewAtype(const IWString& atype,
                const IWString& smiles,
                const IWString& id);

  public:
    SameBondSeparation();

    int Extra(const IWString& atype, const IWString& smiles,
              const IWString& id);

    int NumberFragments() const;

    void set_must_match_atom_types(int s) {
      _must_match_atom_types = s;
    }

    // Read all the files in a directory.
    int Read(const char* dirname);

    // Within a given directory, we create individual files per atom type.
    int Write(const IWString& dirname) const;
  
    int ReplaceFragments(Molecule& parent,
                         Molecule& fragment,
                         const const_IWSubstring& complementary_smiles,
                         const const_IWSubstring& atypes,
                         IWString_and_File_Descriptor& output);
};

SameBondSeparation::SameBondSeparation() {
  _must_match_atom_types = 1;
}

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

int
SameBondSeparation::NumberFragments() const {
  int rc = 0;
  for (const auto& [_, value] : _atype_to) {
    rc += value.NumberFragments();
  }
  return rc;
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

// Parse something that starts with 'at,' and return the part after 'at.'
std::optional<IWString>
ExtractAtomTypes(const std::string& fname) {
  IWString result(fname);
  if (! result.starts_with("at.")) {
    cerr << "ExtractAtomTypes:invalid file name '" << fname << "'\n";
    return std::nullopt;
  }

  result.remove_leading_chars(3);
  return result;
}

int
SameBondSeparation::Read(const char* dirname) {
  const std::filesystem::path dir(dirname);
  for (auto const& dir_entry : std::filesystem::directory_iterator{dir})  {
    std::string path_name = dir_entry.path();
    SameAtomType same_atom_type;
    // cerr << "Reading " << dir_entry.path() << '\n';
    if (! same_atom_type.Read(path_name.c_str())) {
      cerr << "SameBondSeparation::Read:cannot read '" << path_name << "'\n";
      return 0;
    }
    std::optional<IWString> atype = ExtractAtomTypes(dir_entry.path().filename());
    if (! atype) {
      continue;
    }
    auto [_, inserted] = _atype_to.try_emplace(*atype, std::move(same_atom_type));
    if (! inserted) {
      cerr << "SameBondSeparation::Read:cannot insert '" << *atype << '\n';
      return 0;
    }
  }
  return _atype_to.size();
}

int
SameBondSeparation::ReplaceFragments(Molecule& parent,
                         Molecule& fragment,
                         const const_IWSubstring& complementary_smiles,
                         const const_IWSubstring& atypes,
                         IWString_and_File_Descriptor& output) {
  if (_must_match_atom_types) {
    IWString key(atypes);
    key.remove_leading_chars(3);
    const auto iter = _atype_to.find(key);
    if (iter == _atype_to.end()) {
      cerr << "No match for atom type '" << key << "\n";
      return 0;
    }

    return iter->second.ReplaceFragments(parent, fragment, complementary_smiles, output);
  }

  for (auto& [_, same_atom_type] : _atype_to) {
    same_atom_type.ReplaceFragments(parent, fragment, complementary_smiles, output);
  }

  return 1;
}

// a common operation is to break a parent molecule into a core
// and one or more substituents.
// The core will have isotopic labels that allow matching the
// substituents.
struct Decomponsition {
  Molecule& parent;
  Molecule core;
  resizable_array_p<Molecule> substituents;
};

class Options {
  private:
    int _verbose;

    // This class can function in one of two ways, both of which involve
    // reading the output from dicer, with complementary fragments.
    // In the first mode, the dicer output is scanned, and the universe of
    // fragments discovered is accumulated and written to files.
    // In the fragment replacement mode, those files are read, and the
    // input is parsed looking for replacement fragments.
    int _doing_fragment_replacement;

    // the number of lines in the input that contain COMP that we examine.
    int _lines_examined;

    // After the line has been processed, now many molecules do we process.
    int _molecules_processed;

    // A mapping from distance string, '4,5,10', to a collection of
    // fragments that have the same bond distance strings.
    IW_STL_Hash_Map<IWString, SameBondSeparation> _distances;

    IW_STL_Hash_Set _distances_to_find;

    // The directory into which we write directories.
    IWString _file_name_stem;

    // During lookups, we have a set of queries that define the
    // bonds to be broken.
    resizable_array_p<Substructure_Query> _bond_break_queries;

    // when doing lookups, the number of molecules for which we cannot find the
    // canonical distance.
    int _distance_mismatch;

    dicer::Query _query;
    // Some query sets that get instantiated based on items in `_query`.
    resizable_array_p<Substructure_Query> _fragment_must_have;
    resizable_array_p<Substructure_Query> _fragment_must_not_have;

  // Private functions.
    int Process(const Molecule& parent,
                 Molecule& fragment,
                 const const_IWSubstring& atype_string);

    int Write(const IWString& distance_string, const SameBondSeparation& value) const;

    // Recursively cound the number of fragments.
    // The accumulatorA
    int NumberFragments() const;

    int LookingForDistances(const std::string& fname) const;
    int AddDistancesToFind(const const_IWSubstring& input);

    // From a previous run, read the fragment data.
    int Read(IWString& dirname);

    int ReplaceFragments(Molecule& parent,
                Molecule& fragment,
                const const_IWSubstring& complementary_smiles,
                const const_IWSubstring& atypes,
                IWString_and_File_Descriptor& output);
    int MakeCoreAndSubstituents(Molecule& m,
                int& fatal,
                Decomponsition& decomposition);

  public:
    Options();

    int Initialise(Command_Line& cl);

    // `parent` is never changed, but we do ask for its unique smiles, so it is not const.
    int Process(Molecule& parent, const const_IWSubstring& buffer, IWString_and_File_Descriptor& output);

    int Write();

    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _doing_fragment_replacement = 1;
  _molecules_processed = 0;
  _lines_examined = 0;
  _distance_mismatch = 0;
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('R') && cl.option_present('S')) {
    cerr << "Options::Initialise:can specify just one of -R or -S options\n";
    return 0;
  }

  if (cl.option_present('d')) {
    const_IWSubstring d;
    for (int i = 0; cl.value('d', d, i); ++i) {
      if (! AddDistancesToFind(d)) {
        cerr << "Options::Initialise:invalid distance '" << d << "'\n";
        return 0;
      }
    }
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
    _doing_fragment_replacement = 0;

    return 1;
  } 

  if (! cl.option_present('R')) {
    cerr << "When doing fragment replacement must specify -R for directory of existing fragments\n";
    return 0;
  }

  if (!cl.option_present('B')) {
    cerr << "When doing fragment replacement must specify -B for proto directives\n";
    return 0;
  }

  if (cl.option_present('R')) {
    IWString dirname;
    for (int i = 0; cl.value('R', dirname, i); ++i) {
      if (! Read(dirname)) {
        cerr << "Options::Initialise:cannot read '" << dirname << "'\n";
        return 0;
      }
    }
    _doing_fragment_replacement = 1;
    if (_verbose) {
      cerr << "Read " << NumberFragments() << " fragments from '" << dirname << "'\n";
    }
  }

  if (cl.option_present('B')) {
    IWString fname = cl.string_value('B');
    std::optional<dicer::Query> maybe_query = iwmisc::ReadTextProto<dicer::Query>(fname);
    if (! maybe_query) {
      cerr << "Cannot read query options (-B) '" << fname << "'\n";
      return 0;
    }
    _query = std::move(*maybe_query);
  }

  return 1;
}

// Take an unordered set of distances and return a canonical string.
// canonical just means sorted by value.
IWString
CanonicalDistanceString(resizable_array<int>& distances) {
  distances.iwqsort_lambda([](int d1, int d2) {
    if (d1 < d2) {
      return -1;
    } else if (d1 == d2) {
      return 0;
    } else {
      return 1;
    }
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

// The caller is requesting processing of a specific set of
// distances. Convert to canonical form and set _distances_to_find.
int
Options::AddDistancesToFind(const const_IWSubstring& input) {
  resizable_array<int> distances;
  int i = 0;
  const_IWSubstring token;
  while(input.nextword(token, i, ',')) {
    int tmp;
    if (! token.numeric_value(tmp) || tmp < 1) {
      cerr << "Options::AddDistancesToFind:invalid distance '" << token << "'\n";
      return 0;
    }
    distances << tmp;
  }

  if (distances.empty()) {
    cerr << "Options::AddDistancesToFind:no values\n";
    return 0;
  }

  if (distances.size() == 1) {
    IWString key(input);
    _distances_to_find.insert(key);
    return 1;
  }

  IWString key = CanonicalDistanceString(distances);
  _distances_to_find.insert(key);

  return 1;
}

// Parse something that looks like 'ndx:value' and set
// destination[ndx] = value
int
GetAtomType(const const_IWSubstring& buffer,
            extending_resizable_array<uint32_t>& destination) {
  static std::unique_ptr<RE2> rx = std::make_unique<RE2>("^([0-9])\\.([0-9]+)$");

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

// Return the atom if there is a single atom attached to `in_center[ndx]`
// that is further away from all the other atoms in `in_center` than
// in_center[ndx]. This is trying to see if there is a single connection
// from the core to a substituent.
std::optional<atom_number_t>
GetSingleConnection(Molecule& m,
                    const Set_of_Atoms& in_center,
                    int ndx) {
  int n = in_center.size();

  atom_number_t zatom = in_center[ndx];
  const Atom* a = m.atomi(zatom);
  atom_number_t result = INVALID_ATOM_NUMBER;
  for (const Bond* b : *a) {
    if (! b->is_single_bond()) {
      continue;
    }
    if (b->nrings()) {
      continue;
    }
    atom_number_t maybe_outside = b->other(zatom);

    for (int j = 0; j < n; ++j) {
      if (j == ndx) {
        continue;
      }

      if (m.bonds_between(zatom, in_center[j]) <
          m.bonds_between(maybe_outside, in_center[j])) {
        if (result == INVALID_ATOM_NUMBER) {
          result = maybe_outside;
        } else {
          return std::nullopt;
        }
      }
    }
  }

  if (result == INVALID_ATOM_NUMBER) {
    return std::nullopt;
  }

  return result;
}

// Given a set of atoms `in_center`, can we identify a unique set
// of atoms, one per atom in `in_center`, that constitute the
// substituents around the center.
std::optional<Set_of_Atoms>
GetAtomsInSubstituents(Molecule& m,
                       const Set_of_Atoms& in_center) {
  Set_of_Atoms result;

  const int n = in_center.size();
  for (int i = 0; i < n; ++i) {
    std::optional<int> conn = GetSingleConnection(m, in_center, i);
    if (! conn) {
      return 0;
    }
    result << *conn;
  }

  return result;
}

// Recursively visit all atoms attached to `zatom` and mark their
// presence in `to_remove`.
void
IdentifyAtomsBeingRemoved(const Molecule& m, atom_number_t zatom, int * to_remove) {
  to_remove[zatom] = 1;
  const Atom* a = m.atomi(zatom);
  for (const Bond* b : *a) {
    atom_number_t j = b->other(zatom);
    if (to_remove[j]) {
      continue;
    }
    IdentifyAtomsBeingRemoved(m, j, to_remove);
  }
}

// Given a molecule and the queries in `_bond_break_queries`,
// break `m` into a core, which is to be replaced, and substituents
// which are to be preserved. The substituents will each have an isotopic
// label. Store these in `decomponsition`.
// If things fail, return 0, and if a fatal error has occurred, set `fatal`.
// Note that `decomposition` already holds a reference to `m`, but we do
// not bother using that.
// atom typing...
int
Options::MakeCoreAndSubstituents(Molecule& m,
                int& fatal,
                Decomponsition& decomposition) {
  // Force sssr determination.
  m.ring_membership();

  fatal = 0;

  Molecule_to_Match target(&m);

  Set_of_Atoms matched_atoms(_bond_break_queries.size());
  for (Substructure_Query* q : _bond_break_queries) {
    Substructure_Results sresults;
    const int nhits = q->substructure_search(target, sresults);
    if (nhits == 1) {
      matched_atoms << sresults.embedding(0)->front();
    } else if (nhits == 0) {
      if (_query.ignore_molecules_not_hitting_query()) {
        return 0;
      }
      cerr << "Options::MakeCoreAndSubstituents:no match to query " << q->comment() << '\n';
      fatal = 1;
      return 0;
    } else {
      //silently take the first of multiple matches
      matched_atoms << sresults.embedding(0)->front();
    }
  }

  assert(matched_atoms.size() == _bond_break_queries.size());

  std::optional<Set_of_Atoms> atoms_in_substituents = GetAtomsInSubstituents(m, matched_atoms);
  if (! atoms_in_substituents) {
    cerr << "Options::MakeCoreAndSubstituents:cannot identify substituent atoms '" << m.name() << '\n';
    fatal = 1;
    return 0;
  }

  // Now that for each matched atom, we know the adjacent atom that will be
  // retained, break those bonds.

  if (atoms_in_substituents.size() != matched_atoms.size()) {
    cerr << "Size mismatch btw atoms and attachments\n";
    return 0;
  }

  Molecule mcopy(m);
  int n = atoms_in_substituents.number_elements();
  for (int i = 0; i < n; ++i) {
    atom_number_t a1 = matched_atoms[i];
    atom_number_t a2 = atoms_in_substituents[i];
    if (! mcopy.remove_bond_between_atoms(a1, a2)) {
      cerr << "Cannot remove bond btw " << a1 << " and " << a2 << " in " << mcopy.smiles() << '\n';
      return 0;
    }
  }

  std::unique_ptr<int> to_remove(new_int(mcopy.natoms()));
  matched_atoms.set_vector(to_remove.get());

  IdentifyAtomsBeingRemoved(mcopy, matched_atoms[0], to_remove.get());

  return 1;
}

// clang-format off
// buffer looks like
// [3cH]1[2cH]cc[1cH]c1 CHEMBL4553639 [1ClH].[2OH2].Fc1ccc(CN[3CH]=O)cc1 COMP CHEMBL4553639 at.1:3502|2:3502|3:3502
// clang-format on
int
Options::Process(Molecule& parent,
                 const const_IWSubstring& buffer,
                 IWString_and_File_Descriptor& output) {
  ++_lines_examined;

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

  const_IWSubstring complementary_smiles;
  buffer.nextword(complementary_smiles, i);

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

  if (_doing_fragment_replacement) {
    return ReplaceFragments(parent, fragment, complementary_smiles, token, output);
  } else {
    return Process(parent, fragment, token);
  }
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

  return CanonicalDistanceString(distances);

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

// Ok file names are of the form 1.2.5.10
int
OkBondSeparation(const std::string& fname) {
  const_IWSubstring buffer;
  const_IWSubstring token;
  int i = 0;
  int prev_value = -1;
  while (buffer.nextword(token, i, '.')) {
    int value;
    if (! token.numeric_value(value) || value < prev_value) {
      return 0;
    }
    prev_value = value;
  }

  return 1;
}

// We are scanning a directory of distances and need to know
// if we should get the data out of `fname`. Check _distances_to_find.
int
Options::LookingForDistances(const std::string& fname) const {
  if (_distances_to_find.empty()) {
    return 1;
  }
  IWString key(fname);
  key.remove_leading_chars(3);
  return _distances_to_find.contains(key);
}

int
Options::Read(IWString& dirname) {
  const std::filesystem::path dir(dirname.null_terminated_chars());
  for (auto const& dir_entry : std::filesystem::directory_iterator{dir})  {
    std::string path_name = dir_entry.path();
    if (! OkBondSeparation(path_name)) {
      continue;
    }
    if (!LookingForDistances(dir_entry.path().filename())) {
      continue;
    }
    SameBondSeparation sbs;
    if (! sbs.Read(path_name.c_str())) {
      cerr << "Options::Read:fatal error processing '" << path_name << "'\n";
      return 0;
    }
    IWString key(dir_entry.path().filename());
    auto [_, inserted] =  _distances.try_emplace(key, std::move(sbs));
    if (! inserted) {
      cerr << "Options::Read:cannot emplace '" << key << "'\n";
      return 0;
    }
  }

  if (_distances.empty()) {
    cerr << "Options::Read:no files in '" << dirname << "'\n";
    return 0;
  }

  return _distances.size();
}

int
Options::ReplaceFragments(Molecule& parent,
                Molecule& fragment,
                const const_IWSubstring& complementary_smiles,
                const const_IWSubstring& atypes,
                IWString_and_File_Descriptor& output) {
  ++_molecules_processed;

  IWString key = CanonicalDistanceString(fragment);
  auto iter = _distances.find(key);
  // Should be unusual to not have a distance combination.
  if (iter == _distances.end()) {
    ++_distance_mismatch;
    if (_verbose > 2) {
      cerr << "Options::ReplaceFragments:no match for distance '" << key << "'\n";
    }
    return 0;
  }

  return iter->second.ReplaceFragments(parent, fragment, complementary_smiles, atypes, output);
}

int
Options::NumberFragments() const {
  int rc = 0;
  for (const auto& [_, value] : _distances) {
    rc += value.NumberFragments();
  }
  return rc;
}

int
Options::Report(std::ostream& output) const {
  output << _lines_examined << " lines of dicer output examined\n";
  output << "Processed " << _molecules_processed << " molecules\n";
  output << NumberFragments() << " fragments\n";
  if (_doing_fragment_replacement) {
    output << _distance_mismatch << " canonical distances not matched\n";
  }
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

    if (! options.Process(parent, buffer, output)) {
      continue;
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
  Command_Line cl(argc, argv, "vS:R:D:d:");

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

  if (cl.option_present('D')) {
    IWString fname = cl.string_value('D');
    if (! stream_for_transformations.open(fname.null_terminated_chars())) {
      cerr << "Cannot open stream for detailed transformation info '" << fname << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Transformation details written to '" << fname << "'\n";
    }
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

  if (cl.option_present('S')) {
    options.Write();
  }

  return 0;
}
}  // namespace dicer_to_cores

int
main(int argc, char** argv) {
  int rc = dicer_to_cores::Main(argc, argv);

  return rc;
}
