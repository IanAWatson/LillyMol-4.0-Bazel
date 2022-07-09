// Use the output from generate_atom_separations to find
// replacement parts of new molecules.

#include <stdlib.h>

#include <cstdint>
#include <iostream>
#include <optional>

#include "google/protobuf/text_format.h"
#include "google/protobuf/io/zero_copy_stream_impl_lite.h"
#include "leveldb/db.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwbits/fixed_bit_vector.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/atom_separations.h"
#include "Molecule_Tools/embedded_fragment.pb.h"

namespace separated_atoms {

using std::cerr;

using fixed_bit_vector::FixedBitVector;
using embedded_fragment::EmbeddedFragment;

void
Usage(int rc) {
  cerr << "Uses a database built by generate_atom_separations to replace fragments\n";
  cerr << " -d <dbname>      name of previously generated database\n";
  cerr << " -s <smarts>      smarts for break points - specify multiple times\n";
  cerr << " -F <query>       fragments must contain <query>\n";
  cerr << " -P <query>       products  must contain <query>\n";
  cerr << " -c               remove chirality\n";
  cerr << " -l               strip to largest fragment\n";
  cerr << " -v               verbose output\n";
  ::exit(rc);
}

struct Options {
  public:
    int verbose = 0;

    int molecules_read = 0;

    int molecules_not_matching_queries = 0;

    FileType input_type = FILE_TYPE_INVALID;

    Chemical_Standardisation chemical_standardisation;

    int reduce_to_largest_fragment = 0;

    int remove_chirality = 0;

    // If reading fragments from a leveldb database.
    std::unique_ptr<leveldb::DB> database;

    // The bonds to be broken are specified by some number of
    // queries. If two sets of queries, then two bond break
    // fragments are examined.
    resizable_array_p<Substructure_Query> queries1;  // Not used (yet).
    resizable_array_p<Substructure_Query> queries2;
    resizable_array_p<Substructure_Query> queries3;

    Hasher hasher;

    // We can require a minimum number of examplars for doing
    // a replacement.
    uint32_t min_examplar_count = 0;
    int too_few_examplars = 0;

    // We can specify what the linker must look like - before it
    // is inserted.
    resizable_array_p<Substructure_Query> fragment_must_have;
    // We can specify what the final product must look like.
    resizable_array_p<Substructure_Query> product_must_have;

    // When doing a replacement, do we need to match the
    // atomic numbers that are stored in the database.
    int match_atomic_numbers = 0;

    IWString fragment_name_separator = " %% ";

    int molecules_generated = 0;
    int molecules_written = 0;

  // private functions

    int OpenLevelDb(IWString& database_name);

    int DoFragmentReplacement2(Molecule& m,
                               const Set_of_Atoms& bonds_to_break,
                               int * tmp,
                               IWString_and_File_Descriptor& output);
    int DoFragmentReplacement3(Molecule& m,
                               const Set_of_Atoms& bonds_to_break,
                               int * tmp,
                               IWString_and_File_Descriptor& output);

    // Do substructure searches to identify 2, 4, 6... atoms that define the
    // 1,2,3... bonds to be broken.
    std::optional<Set_of_Atoms> IdentifyBondsToBreak(Molecule& m);

    std::optional<std::string> ReadDb(uint32_t hash);

    // Return true if proto.count() is consistent with min_examplar_count.
    int OkExemplarCount(const embedded_fragment::EmbeddedFragment& proto);

    // Do a two bond fragment replacement in `m`.
    // a11-a12...a22-a21
    // where the ... atoms are replaced by the fragment in `proto`.
    // Those atoms will be set in `atoms_between`.
    int DoReplacement(Molecule& m,
                       atom_number_t a11,
                       atom_number_t a12,
                       atom_number_t a21,
                       atom_number_t a22,
                       const FixedBitVector& atoms_between,
                       embedded_fragment::EmbeddedFragment& proto,
                       int * storage,
                       IWString_and_File_Descriptor& output);

    int FragmentContainsRequiredSubstructure(Molecule& m);
    int ProductContainsRequiredSubstructure(Molecule& m);

    int Write(Molecule& m,
               const EmbeddedFragment& proto,
               IWString_and_File_Descriptor& output);

  public:
    Options() {
    }

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);
    
    int DoFragmentReplacement(Molecule& m, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

int
Options::Initialise(Command_Line& cl) {
  verbose = cl.option_count('v');

  if (cl.option_present('c')) {
    remove_chirality = 1;
    if (verbose) {
      cerr << "Will remove chirality from input molecules\n";
    }
  }

  if (cl.option_present('g')) {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose, 'g')) {
      cerr << "Cannot initialise chemical standardisation\n";
      return 0;
    }
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;
    if (verbose) {
      cerr << "Will reduce molecules to largest fragment\n";
    }
  }

  if (cl.option_present('F')) {
    if (! process_queries(cl, fragment_must_have, verbose, 'F')) {
      cerr << "Options::Initialiseannot read fragment must have queries (-F)\n";
      return 0;
    }
  }

  if (cl.option_present('P')) {
    if (! process_queries(cl, product_must_have, verbose, 'P')) {
      cerr << "Options::Initialiseannot read product must have queries (-P)\n";
      return 0;
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring s;
    for (int i = 0; cl.value('s', s, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(s)) {
        cerr << "Options::Initialise:cannot parse smarts '" << s << "'\n";
        return 0;
      }
      if (queries1.empty()) {
        queries1 << q.release();
      } else if (queries2.empty()) {
        queries2 << q.release();
      } else if (queries3.empty()) {
        queries3 << q.release();
      } else {
        cerr << "Too many bond breaking queries (-s), cannot handle\n";
        return 0;
      }
    }
  }

  if (! cl.option_present('d')) {
    cerr << "Options::Initialise:must specify name of fragment database via the -d option\n";
    return 0;
  } else {
    IWString dbname = cl.string_value('d');
    if (! OpenLevelDb(dbname)) {
      cerr << "Options::Initialise:cannot open leveldb '" << dbname << "'\n";
      return 0;
    }
  }

  if (cl.number_elements() == 1 && strcmp("-", cl[0]) == 0) { // reading a pipe, assume smiles
    input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern all file types, use the -i option\n";
    return 0;
  } else if (!process_input_type(cl, input_type)) {
    return 0;
  }

  return 1;
}

int
Options::OpenLevelDb(IWString& database_name) {
  const std::string dbname = database_name.AsString();

  leveldb::Options leveldb_options;
  leveldb_options.create_if_missing = false;

  leveldb::DB* db;
  const leveldb::Status status = leveldb::DB::Open(leveldb_options, dbname, &db);
  if (!status.ok()) {
    cerr << "Options::OpenLevelDb:cannot open " << std::quoted(dbname) << " " << status.ToString() << "\n";
    return 0;
  }

  leveldb::Iterator* it = db->NewIterator(leveldb::ReadOptions());
  cerr << "Check iteration" << it << '\n';
  int items_in_db = 0;
  for (it->SeekToFirst(); it->Valid(); it->Next()) {
//  cerr << it->key().ToString() << ": "  << it->value().ToString() << '\n';
    ++items_in_db;
  }
  cerr << it->status().ToString() << '\n';
  cerr << items_in_db << " items in database\n";
  assert(it->status().ok());  // Check for any errors found during the scan
  delete it;

  database.reset(db);
  return 1;
}

int
Options::Report(std::ostream& output) const {
  cerr << "Read " << molecules_read << " molecules\n";
  cerr << molecules_not_matching_queries << " molecules did not match bond breaking queries\n";
  cerr << molecules_generated << " molecules generated " << molecules_written << " molecules written\n";

  return output.good();
}

int
Options::Preprocess(Molecule& m) {
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  if (remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (m.natoms() == 0) {
    return 0;
  }

  return 1;
}

int
FirstMatch(resizable_array_p<Substructure_Query>& queries,
           Molecule_to_Match& target,
           Set_of_Atoms& destination) {
  Substructure_Results sresults;
  for (Substructure_Query* q : queries) {
    const int nhits = q->substructure_search(target, sresults);
    if (nhits == 0) {
      continue;
    }
    for (int i = 0; i < nhits; ++i) {
      const Set_of_Atoms * e = sresults.embedding(i);
      if (e->number_elements() < 2) {  // Should not happen.
        continue;
      }
      destination << e->item(0);
      destination << e->item(1);
      return 1;
    }
  }

  return 0;
}

std::optional<Set_of_Atoms>
Options::IdentifyBondsToBreak(Molecule& m) {
  Molecule_to_Match target(&m);

  Set_of_Atoms result;
  if (queries1.size() > 0) {
    if (! FirstMatch(queries1, target, result)) {
      return std::nullopt;
    }
  }
  if (queries2.size() > 0) {
    if (! FirstMatch(queries2, target, result)) {
      return std::nullopt;
    }
  }
  if (queries3.size() > 0) {
    if (! FirstMatch(queries3, target, result)) {
      return std::nullopt;
    }
  }
  cerr << "Result of substructure matching " << result << '\n';

  return result;
}

int
Options::DoFragmentReplacement(Molecule& m,
                               IWString_and_File_Descriptor& output) {
  std::optional<Set_of_Atoms> bonds_to_break = IdentifyBondsToBreak(m);
  if (! bonds_to_break) {
    ++molecules_not_matching_queries;
    if (verbose > 1) {
      cerr << m.smiles() << ' ' << m.name() << ' ' << " no query matches\n";
      return 0;
    }
  }

  std::unique_ptr<int[]> storage = std::make_unique<int[]>(m.natoms() * 2);

  if (bonds_to_break->size() == 6) {
    DoFragmentReplacement3(m, *bonds_to_break, storage.get(), output);
  }
  if (bonds_to_break->size() == 4) {
    DoFragmentReplacement2(m, *bonds_to_break, storage.get(), output);
  }
  if (bonds_to_break->size() == 2) {
    cerr << "Single frgament replacement not implemented\n";
    return 0;
  }

  return 1;
}


std::optional<std::string>
Options::ReadDb(uint32_t hash) {
  const leveldb::Slice key = Key(hash);

#ifdef ITERATE_DB_HERE
  leveldb::Iterator* it = database->NewIterator(leveldb::ReadOptions());
  cerr << "Check iteration" << it << '\n';
  for (it->SeekToFirst(); it->Valid(); it->Next()) {
    cerr << it->key().ToString() << ": "  << it->value().ToString() << '\n';
  }
  assert(it->status().ok());  // Check for any errors found during the scan
  delete it;
#endif

  leveldb::ReadOptions read_options;
  std::string result;
  const leveldb::Status status = database->Get(read_options, key, &result);
  cerr << "looked up " << hash << " status " << status.ok() << '\n';
  if (status.ok()) {
    return result;
  }

  return std::nullopt;
}
  
int
Options::DoFragmentReplacement2(Molecule& m,
                                const Set_of_Atoms& bonds_to_break,
                                int * tmp,
                                IWString_and_File_Descriptor& output) {
  assert(bonds_to_break.size() == 4);

  const int matoms = m.natoms();

  std::fill_n(tmp, matoms, 0);
  const int c1 = m.identify_side_of_bond(tmp, bonds_to_break[0], 1, bonds_to_break[1]);
  cerr << "c1 " << c1 << '\n';
  if (c1 == 0) {  // Likely a ring.
    return 0;
  }

  FixedBitVector b0, b1;
  b0.ConstructFromArray(tmp, matoms, 1);
  b1.ConstructFromArray(tmp, matoms, 0);

  std::fill_n(tmp, matoms, 0);
  const int c2 = m.identify_side_of_bond(tmp, bonds_to_break[2], 1, bonds_to_break[3]);
  cerr << "c2 " << c2 << '\n';
  if (c2 == 0) {  // Likely a ring.
    return 0;
  }

  FixedBitVector b2, b3;
  b2.ConstructFromArray(tmp, matoms, 1);
  b3.ConstructFromArray(tmp, matoms, 0);
  FixedBitVector atoms_between;

  atom_number_t a11, a12, a21, a22;
  if (b0.BitsInCommon(b2) == 0) {
    a11 = bonds_to_break[0];
    a12 = bonds_to_break[1];
    a21 = bonds_to_break[2];
    a22 = bonds_to_break[3];
    atoms_between = b1;
    atoms_between.iwand(b3);
  } else if (b0.BitsInCommon(b3) == 0) {
    a11 = bonds_to_break[0];
    a12 = bonds_to_break[1];
    a21 = bonds_to_break[3];
    a22 = bonds_to_break[2];
    atoms_between = b1;
    atoms_between.iwand(b2);
  } else if (b1.BitsInCommon(b2) == 0) {
    a11 = bonds_to_break[1];
    a12 = bonds_to_break[0];
    a21 = bonds_to_break[2];
    a22 = bonds_to_break[3];
    atoms_between = b0;
    atoms_between.iwand(b3);
  } else if (b1.BitsInCommon(b3) == 0) {
    a11 = bonds_to_break[1];
    a12 = bonds_to_break[0];
    a21 = bonds_to_break[3];
    a22 = bonds_to_break[2];
    atoms_between = b0;
    atoms_between.iwand(b2);
  } else {
    cerr << "DoFragmentReplacement2:cannot find common atoms\n";
    return 0;
  }

  int d12 = m.bonds_between(a12, a22);

  uint32_t hash = hasher.Value(m.atomic_number(a12), m.atomic_number(a22), d12);
  cerr << "Distance is " << d12 << " hash " << hash << '\n';

  std::optional<std::string> maybe_data = ReadDb(hash);
  if (! maybe_data) {
    cerr << "No data for " << hash << '\n';
    return 0;
  }

  const_IWSubstring buffer(maybe_data->data(), maybe_data->size());
  int i = 0;
  const_IWSubstring line;
  while (buffer.nextword(line, i, '\n')) {
    embedded_fragment::EmbeddedFragment proto;
    google::protobuf::io::ArrayInputStream buffer(line.data(), line.length());
    google::protobuf::TextFormat::Parse(&buffer, &proto);
    cerr << "Looking up based on " << line << '\n';
    DoReplacement(m, a11, a12, a21, a22, atoms_between, proto, tmp, output);
  }

  return 1;
}

// Identify two isotopically labelled atoms in `m`,
int
FetchTwoIsotopicallyLabeledAtoms(const Molecule& m,
                                 atom_number_t& a1,
                                 atom_number_t& a2) {
  a1 = INVALID_ATOM_NUMBER;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.isotope(i) == 0) {
      continue;
    }

    if (a1 == INVALID_ATOM_NUMBER) {
      a1 = i;
    } else {
      a2 = i;
      return 1;
    }
  }

  return 0;
}

// Atoms `to_be_removed` are going to be removed from a molecule
// containing `matoms` atoms.
// Fill in `xref` which is a mapping from old atom numbers to atom
// numbers in what is left.
// Returns the number of items that remain.
int
EstablishXref(const int* to_be_removed,
              int matoms,
              int * xref) {
  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (to_be_removed[i]) {
      xref[i] = -1;
    } else {
      xref[i] = rc;
      ++rc;
    }
  }
  return rc;
}

int
Options::OkExemplarCount(const embedded_fragment::EmbeddedFragment& proto) {
  if (proto.count() < min_examplar_count) {
    ++too_few_examplars;
    return 0;
  }

  return 1;
}

int
Options::DoReplacement(Molecule& m,
                       atom_number_t a11,
                       atom_number_t a12,
                       atom_number_t a21,
                       atom_number_t a22,
                       const FixedBitVector& atoms_between,
                       embedded_fragment::EmbeddedFragment& proto,
                       int * storage,
                       IWString_and_File_Descriptor& output) {
  if (! OkExemplarCount(proto)) {
    return 0;
  }
  cerr << "Exemplar count ok " << proto.count() << '\n';
  cerr << "a11 " << a11 << " a12 " << a12 << " a21 " << a21 << " a22 " << a22 << '\n';
  const int matoms = m.natoms();
  std::fill_n(storage, matoms, 0);
  atoms_between.Scatter(storage, 1);
  Molecule m2(m);
  for (int i = 0; i < matoms; ++i) {
    m2.set_isotope(i, i);
  }
  cerr << "Starting mol " << m2.smiles() << '\n';

  cerr << "Atoms between\n";
  atoms_between.DebugPrint(cerr);
  int * xref = storage + matoms;
  EstablishXref(storage, matoms, xref);
  Molecule mcopy(m);
  mcopy.remove_atoms(storage);
  for (int i = 0; i < matoms; ++i) {
    cerr << i << " " << storage[i] << " xref " << xref[i] << '\n';
  }
  cerr << "After atom removal\n";
  cerr << mcopy.smiles() << '\n';

  Molecule fragment;
  if (! fragment.build_from_smiles(proto.smiles())) {
    cerr << "Options::DoReplacement:cannot parse fragment smiles " << std::quoted(proto.smiles()) << '\n';
    return 0;
  }
  if (! FragmentContainsRequiredSubstructure(fragment)) {
    cerr << "Fragment does not contain required substructure\n";
    return 0;
  }

  // The atom numbers in the fragment that we will join.
  atom_number_t f1 = INVALID_ATOM_NUMBER;
  atom_number_t f2 = INVALID_ATOM_NUMBER;
  if (! FetchTwoIsotopicallyLabeledAtoms(fragment, f1, f2)) {
    cerr << "Options::DoReplacement:cannot find isotopes " << fragment.smiles() << '\n';
    return 0;
  }

  int initial_mcopy_atoms = mcopy.natoms();
  mcopy.add_molecule(&fragment);

  // At this stage, we may generate 1 or 2 different molecules. If there are no
  // constraints on matching atomic numbers, then we generate 2 molecules.
  // If there are constraints on atomic numbers, we may generate either 1 or 2.

  // Do the easy case first.
  cerr << "match_atomic_numbers " << match_atomic_numbers << '\n';
  if (! match_atomic_numbers) {
    mcopy.add_bond(xref[a11], initial_mcopy_atoms + f1, SINGLE_BOND);
    mcopy.add_bond(xref[a21], initial_mcopy_atoms + f2, SINGLE_BOND);
    Write(mcopy, proto, output);
    mcopy.remove_bond_between_atoms(xref[a11], initial_mcopy_atoms + f1);
    mcopy.remove_bond_between_atoms(xref[a21], initial_mcopy_atoms + f2);

    mcopy.add_bond(xref[a11], initial_mcopy_atoms + f2, SINGLE_BOND);
    mcopy.add_bond(xref[a21], initial_mcopy_atoms + f1, SINGLE_BOND);
    return Write(mcopy, proto, output);
  }

  // Now the harder case where isotopes must match.
  if (fragment.isotope(f1) == m.atomic_number(a11) &&
      fragment.isotope(f2) == m.atomic_number(a21)) {
    mcopy.add_bond(xref[a11], initial_mcopy_atoms + f1, SINGLE_BOND);
    mcopy.add_bond(xref[a21], initial_mcopy_atoms + f2, SINGLE_BOND);
    Write(mcopy, proto, output);
    mcopy.remove_bond_between_atoms(xref[a11], initial_mcopy_atoms + f1);
    mcopy.remove_bond_between_atoms(xref[a21], initial_mcopy_atoms + f2);
  }

  if (fragment.isotope(f1) == m.atomic_number(a21) &&
      fragment.isotope(f2) == m.atomic_number(a11)) {
    mcopy.add_bond(xref[a11], initial_mcopy_atoms + f2, SINGLE_BOND);
    mcopy.add_bond(xref[a21], initial_mcopy_atoms + f1, SINGLE_BOND);
    return Write(mcopy, proto, output);
  }

  return 1;  // We may or may not have generated a molecule here.
}

int
MatchAnyQuery(Molecule& m,
              resizable_array_p<Substructure_Query>& queries) {
  Molecule_to_Match target(&m);
  for (Substructure_Query* q : queries) {
    if (q->substructure_search(target)) {
      return 1;
    }
  }

  return 0;
}

int
Options::FragmentContainsRequiredSubstructure(Molecule& m) {
  if (fragment_must_have.empty()) {
    return 1;
  }
  return MatchAnyQuery(m, fragment_must_have);
}

int
Options::ProductContainsRequiredSubstructure(Molecule& m) {
  if (product_must_have.empty()) {
    return 1;
  }
  return MatchAnyQuery(m, product_must_have);
}

// If `m` passes any substructure constraints that might be in effect,
// write it to `output`.
// The name is formed by combining m.name() with the name of the
// exemplar found in `proto`.
int
Options::Write(Molecule& m,
               const EmbeddedFragment& proto,
               IWString_and_File_Descriptor& output) {
  ++molecules_generated;
  if (! ProductContainsRequiredSubstructure(m)) {
    return 0;
  }

  constexpr char kSep = ' ';

  IWString name(m.name());
  name << fragment_name_separator;
  name << proto.exemplar() << proto.count();
  output << m.smiles() << kSep << name << '\n';
  output.write_if_buffer_holds_more_than(8192);

  ++molecules_written;
  return 1;
}

int
Options::DoFragmentReplacement3(Molecule& m,
                                const Set_of_Atoms& bonds_to_break,
                                int * tmp,
                                IWString_and_File_Descriptor& output) {
  return 1;
}

int
FragmentReplacement(Options& options,
                    Molecule& m,
                    IWString_and_File_Descriptor& output) {
  if (!options.Preprocess(m)) {
    return 0;
  }

  return options.DoFragmentReplacement(m, output);
}

int
FragmentReplacement(Options& options,
                    data_source_and_type<Molecule>& input,
                    IWString_and_File_Descriptor& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);
    ++options.molecules_read;
    if (! FragmentReplacement(options, *m, output)) {
      cerr << "Fatal error processing " << m->smiles() << ' ' << m->name() << '\n';
      return 0;
    }
  }

  return 1;
}

int
FragmentReplacement(Options& options,
                    const char * fname,
                    IWString_and_File_Descriptor& output) {
  FileType input_type = options.input_type;
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose > 1) {
    input.set_verbose(1);
  }

  return FragmentReplacement(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:g:i:lcd:F:P:s:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (! process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process aromaticity options\n";
    return 1;
  }

  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process elements\n";
    return 1;
  }

  // Helpful for when making fragments
  set_copy_name_in_molecule_copy_constructor(1);

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  for (const char * fname : cl) {
    if (! FragmentReplacement(options, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace separated_atoms

int
main(int argc, char** argv) {
  int rc = separated_atoms::Main(argc, argv);

  return rc;
}
