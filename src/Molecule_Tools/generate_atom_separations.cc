// Within a molecule, generate smarts corresponding to separated atoms.
// This is a generalisation of ring_extraction/ring_replacement.

#include <algorithm>
#include <iostream>
#include <memory>
#include <unordered_map>

#include "google/protobuf/text_format.h"
#include "leveldb/db.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwbits/fixed_bit_vector.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/atom_separations.h"
#include "Molecule_Tools/embedded_fragment.pb.h"

namespace separated_atoms {

constexpr int kDefaultBitVector = 64;

using fixed_bit_vector::FixedBitVector;

using std::cerr;

using embedded_fragment::EmbeddedFragment;
// For each arrangement of join points, we have a number of possibilities.
class SetOfLinkers {
  private:
    IW_STL_Hash_Map<IWString, EmbeddedFragment> _usmi_to_frag;
  public:
    // Notify of an additional two bond breakage fragment.
    int Extra(Molecule& parent, Molecule& m, int distance);
    // Notify of an additional three bond breakage fragment.
    int Extra(int store_exemplar_smiles, Molecule& parent, Molecule& m, const int d12, int d13, int d23);

    int Write(int store_exemplar_smiles, IWString_and_File_Descriptor& output) const;

    int WriteTextProtos(IWString& output) const;
};

int
SetOfLinkers::Extra(int store_exemplar_smiles,
                    Molecule& parent,
                    Molecule& m,
                    int distance) {
  const IWString& usmi = m.unique_smiles();

  auto iter = _usmi_to_frag.find(usmi);
  if (iter != _usmi_to_frag.end()) {
    iter->second.set_count(iter->second.count() + 1);
    return 1;
  }

  EmbeddedFragment fragment;
  fragment.set_smiles(usmi.data(), usmi.length());
  fragment.set_exemplar(parent.name().data(), parent.name().length());
  fragment.set_count(1);
  if (store_exemplar_smiles) {
    fragment.set_exemplar_smiles(parent.unique_smiles().data(), parent.unique_smiles().length());
  }
  fragment.mutable_dist()->Add(distance);
  _usmi_to_frag.emplace(usmi, std::move(fragment));

  return 1;
}

int
SetOfLinkers::Extra(int store_exemplar_smiles,
                    Molecule& parent, Molecule& m,
                    int d12, int d13, int d23) {
  const IWString& usmi = m.unique_smiles();

  auto iter = _usmi_to_frag.find(usmi);
  if (iter != _usmi_to_frag.end()) {
    iter->second.set_count(iter->second.count() + 1);
    return 1;
  }

  EmbeddedFragment fragment;
  fragment.set_smiles(usmi.data(), usmi.length());
  fragment.set_exemplar(parent.name().data(), parent.name().length());
  fragment.set_count(1);
  if (store_exemplar_smiles) {
    fragment.set_exemplar_smiles(parent.unique_smiles().data(), parent.unique_smiles().length());
  }
  fragment.mutable_dist()->Add(d12);
  fragment.mutable_dist()->Add(d13);
  fragment.mutable_dist()->Add(d23);
  _usmi_to_frag.emplace(usmi, std::move(fragment));

  return 1;
}

int
SetOfLinkers::Write(IWString_and_File_Descriptor& output) const {
  constexpr char kSep = ' ';
  for (const auto& [usmi, proto] : _usmi_to_frag) {
    output << proto.exemplar_smiles() << kSep << proto.exemplar() << '\n';
    output << proto.smiles() << '\n';
#ifdef PRD_VERSION
    output << proto.smiles() << kSep <<
              proto.exemplar() << kSep <<
              proto.count();
    for (uint32_t d : proto.distance()) {
      output << kSep << d;
    }
    output << '\n';
#endif
    output.write_if_buffer_holds_more_than(8192);
  }
  return 1;
}

int
SetOfLinkers::WriteTextProtos(IWString& output) const {

  output.resize(100 * _usmi_to_frag.size());  // Just a guess.

  google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);
  std::string buffer;
  for (const auto& [usmi, proto] : _usmi_to_frag) {
    printer.PrintToString(proto, &buffer);
    output << buffer;
    output << '\n';
  }
  return 1;
}

struct Job {
  int verbose = 0;

  Chemical_Standardisation chemical_standardisation;

  int reduce_to_largest_fragment = 0;

  int remove_chirality = 0;

  FileType input_type = FILE_TYPE_INVALID;

  int break_amide_bonds = 1;

  int fragments_must_have_ring = 0;

  // By default, we only break bonds where 

  resizable_array<int> bond_separations;

  int molecules_read = 0;

  int molecules_too_small = 0;

  // For each number of breakages, the number of fragments generated - before filtering.
  extending_resizable_array<int> fragments_generated;
  // The number of fragments written - per number of broken bonds.
  extending_resizable_array<int> fragments_written;

  // If specified, only break this kind of bond.
  resizable_array_p<Substructure_Query> break_only;
  // If specified, prevent these bonds from breaking.
  resizable_array_p<Substructure_Query> do_not_break;

  // After subsets are generated, we can get rid of those not
  // matching one of these queries.
  resizable_array_p<Substructure_Query> fragments_must_have;
  int fragments_discarded_for_no_query_match = 0;

  // We can impose a limit on the distance between broken bonds
  int max_bond_separation = std::numeric_limits<int>::max();
  // We can impose limits on the number of atoms in
  // the fragments we generate.
  int min_frag_size = 0;
  int max_frag_size = std::numeric_limits<int>::max();

  // Across all molecules, counters for the number of breakable bonds.
  extending_resizable_array<int> breakable_bonds;
  Accumulator<int> acc_breakable_bonds;

  // Some molecules have too many breakable bonds.
  int skip_if_too_many_breakable_bonds = std::numeric_limits<int>::max();
  // A cound of the number of proposed breakages that are discarded
  // for violating skip_if_too_many_breakable_bonds
  int bond_separation_too_long = 0;

  // The number of molecules violating the constraint.
  int molecules_with_too_many_breakable_bonds = 0;

  // A mapping from hash value to fragments found.
  // The suffix is for the number of bond breaks.
  std::unordered_map<uint32_t, SetOfLinkers> fragments;
  std::unordered_map<uint32_t, SetOfLinkers> fragments3;

  // File name stem for output files.
  IWString stem;

  // If writing to a leveldb database.
  IWString database_name;
  std::unique_ptr<leveldb::DB> database;

  Hasher hasher;

  // By default, we do not store the full exemplar smiles
  // in the database - just to save space.
  int store_exemplar_smiles = 0;

  // private functions.
    int WriteLinkers(uint32_t hash,
                     const SetOfLinkers& linkers) const;
    int WriteLinkersDb(uint32_t hash,
                       const SetOfLinkers& linkers) const;

    int IdentifyBreakableBonds(Molecule& m, int * breakable);

    int OpenLevelDb();

  public:
    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int InitialiseSeparationNeeded(const int max_sep, int* needed);

    // If either min_frag_size or max_frag_size are set, see if the number
    // atoms implied by `atoms_between` is consistent with those limits.
    int OkSize(const FixedBitVector& atoms_between) const;

    // Return whether or not `bonds` is consistent with
    // skip_if_too_many_breakable_bonds.
    // Updates the accumulators on numbers of breakable bonds.
    int OkBreakableBondCount(int bonds);

    // Return whether or not `distance` is consistent with max_bond_separation.
    int OkBondSeparation(int distance);

    void AnotherFragment(int nbreak) {
      ++fragments_generated[nbreak];
    }

    int OkAtomicNumber(int z) const {
      return hasher.OkAtomicNumber(z);
    }

    // A pair has two atomic numbers and a distance.
    uint32_t HashValue(int z1, int z2, int d) const {
      return hasher.Value(z1, z2, d);
    }

    // A triple has 3 atoms and 3 distances.
    uint32_t HashValue(int z1, int z2, int z3, int d12, int d13, int d23) const {
      return hasher.Value(z1, z2, z3, d12, d12, d23);
    }

    int MatchesMustHaveQuery(Molecule& m);

    int WriteLinkers() const;

    int Report(std::ostream& output) const;
};

int
Job::Initialise(Command_Line& cl) {
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

  if (cl.option_present('f')) {
    if (! cl.value('f', min_frag_size) || min_frag_size < 0) {
      cerr << "Job::Initialise:invalid min frag size (-f)\n";
      return 0;
    }

    if (verbose) {
      cerr << "Will discard fragments with fewer than " << min_frag_size << " atoms\n";
    }
  }

  if (cl.option_present('F')) {
    if (! cl.value('F', max_frag_size) || max_frag_size < min_frag_size) {
      cerr << "Job::Initialise:invalid max frag size (-F)\n";
      return 0;
    }

    if (verbose) {
      cerr << "Will discard fragments with more than " << max_frag_size << " atoms\n";
    }
  }

  if (cl.option_present('B')) {
    if (! cl.value('B', skip_if_too_many_breakable_bonds) || skip_if_too_many_breakable_bonds < 0) {
      cerr << "Job::Initialise:Invalid too many breakable bonds (-B)\n";
      return 0;
    }

    if (verbose) {
      cerr << "Molecules with more than " << skip_if_too_many_breakable_bonds << " will be ignored\n";
    }
  }

  if (cl.option_present('D')) {
    if (! cl.value('D', max_bond_separation) || max_bond_separation < 1) {
      cerr << "The max bond esparation must be a whole +ve number\n";
      return 0;
    }
    if (verbose) {
      cerr << "Will discard breakages involving bond separation greater than " << max_bond_separation << '\n';
    }
  }

  if (cl.option_present('K')) {
    if (! process_queries(cl, break_only, verbose, 'K')) {
      cerr << "Cannot read break only queries (-K)\n";
      return 0;
    }
  }

  if (cl.option_present('k')) {
    if (! process_queries(cl, do_not_break, verbose, 'k')) {
      cerr << "Cannot read break only queries (-k)\n";
      return 0;
    }
  }

  if (cl.option_present('m')) {
    if (! process_queries(cl, fragments_must_have, verbose, 'm')) {
      cerr << "Cannot read fragments must have queries (-m)\n";
      return 0;
    }
  }

  if (cl.option_present('Y')) {
    const_IWSubstring y;
    for (int i = 0; cl.value('y', y, i); ++i) {
      if (y == "nbamide") {
        break_amide_bonds = 0;
      } else if (y == "fring") {
        fragments_must_have_ring = 1;
      } else if (y == "help") {
      } else {
        cerr << "Job::Initialise:unrecognised -Y directive '" << y << "'\n";
        return 0;
      }
    }
  }

  if (cl.option_present('S')) {
    cl.value('S', stem);
    if (verbose) {
      cerr << "Results written to file stem '" << stem << "'\n";
    }
  } else if (cl.option_present('d')) {
    cl.value('d', database_name);
    if (! OpenLevelDb()) {
      cerr << "Job::Initialise:cannot open database " << database_name << '\n';
      return 0;
    }
    if (verbose) {
      cerr << "Results written to database '" << database_name << "'\n";
    }
  } else {
    cerr << "Must specify output file name stem via the -S option\n";
    return 0;
  }


  if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) { // reading a pipe, assume smiles
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
Job::OpenLevelDb() {
  const std::string dbname = database_name.AsString();

  leveldb::Options leveldb_options;
  leveldb_options.create_if_missing = true;

  leveldb::DB *db;
  const leveldb::Status status = leveldb::DB::Open(leveldb_options, dbname, &db);
  if (!status.ok()) {
    cerr << "Job::OpenLevelDb:cannot open " << std::quoted(dbname) << " " << status.ToString() << "\n";
    return 0;
  }

  database.reset(db);
  return 1;
}

int
Job::Report(std::ostream& output) const {
  output << "Read " << molecules_read << " molecules\n";
  if (fragments_must_have.size() > 0) {
    output << fragments_discarded_for_no_query_match <<
             " fragments discarded for not matching " <<
             fragments_must_have.size() << " queries\n";
  }

  for (int i = 0; i < breakable_bonds.number_elements(); ++i) {
    if (breakable_bonds[i]) {
      output << breakable_bonds[i] << " molecules had " << i << " breakable bonds\n";
    }
  }
  output << static_cast<float>(acc_breakable_bonds.average()) << " mean number of bonds broken\n";

  if (skip_if_too_many_breakable_bonds < std::numeric_limits<int>::max()) {
    output << molecules_with_too_many_breakable_bonds << " molecules had more than " <<
              skip_if_too_many_breakable_bonds << " breakable bonds\n";
  }

  for (int i = 0; i < fragments_generated.number_elements(); ++i) {
    if (fragments_generated[i] > 0) {
      output << fragments_generated[i] << " fragments from " << i << " breakages\n";
    }
  }

  return output.good();
}

int
Job::Preprocess(Molecule& m) {
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
Job::InitialiseSeparationNeeded(const int max_sep, int* needed) {
  for (int d : bond_separations) {
    if (d >= max_sep) {
      ++molecules_too_small;
      return 0;
    }

    needed[d] += 1;
  }

  return 1;
}

// Run a set of substructure queries against `m`.
// Interpret the first two atoms of each match as describing a
// pair of bonded atoms.
// If not bonded, return 0;
// Otherwise update `result` with `flag` for that bond number.
int
SetMatchedAtoms(resizable_array_p<Substructure_Query>& queries,
                Molecule& m,
                int flag,
                int * result) {
  Molecule_to_Match target(&m);
  int rc = 0;
  for (Substructure_Query * q : queries) {
    Substructure_Results sresults;
    const int nhits = q->substructure_search(target, sresults);
    if (nhits == 0) {
      continue;
    }
    for (int i = 0; i < nhits; ++i) {
      const Set_of_Atoms * e = sresults.embedding(i);
      const atom_number_t a0 = e->item(0);
      const atom_number_t a1 = e->item(1);
      const int  bond_number = m.which_bond(a0, a1);
      if (bond_number < 0) {
        cerr << "Matched atoms not bonded " << m.name() << ' ' << a0 << ' ' << a1 << '\n';
        return 0;
      }
      result[bond_number] = flag;
      ++rc;
    }
  }

  return rc;
}
int
Job::IdentifyBreakableBonds(Molecule& m, int * breakable) {
  const int nedges = m.nedges();

  if (! break_amide_bonds) {
    TurnOffAmideBonds(m, breakable);
  }

  if (break_only.empty() && do_not_break.empty()) {
    std::fill_n(breakable, nedges, 1);
    return nedges;
  }

  if (break_only.size() > 0) {
    std::fill_n(breakable, nedges, 0);
    return SetMatchedAtoms(break_only, m, 1, breakable);
  }

  std::fill_n(breakable, nedges, 1);
  return SetMatchedAtoms(do_not_break, m, 0, breakable);
}

int
Job::OkSize(const FixedBitVector& atoms_between) const {
  if (min_frag_size == 0 && max_frag_size == std::numeric_limits<int>::max()) {
    return 1;
  }

  const int natoms = atoms_between.nset();
  if (natoms < min_frag_size || natoms > max_frag_size) {
    return 0;
  }

  return 1;
}

int
Job::OkBreakableBondCount(int bonds) {
  ++breakable_bonds[bonds];
  acc_breakable_bonds.extra(bonds);

  if (bonds < skip_if_too_many_breakable_bonds) {
    return 1;
  }

  ++molecules_with_too_many_breakable_bonds;
  return 0;
}

int
Job::OkBondSeparation(int distance) {
  if (distance <= max_bond_separation) {
    return 1;
  }

  ++bond_separation_too_long;
  return 0;
}

int
Job::MatchesMustHaveQuery(Molecule& m) {
  Molecule_to_Match target(&m);
  for (Substructure_Query * q : fragments_must_have) {
    Substructure_Results sresults;
    const int nhits = q->substructure_search(target, sresults);
    if (nhits) {
      return 1;
    }
  }

  ++fragments_discarded_for_no_query_match;
  return 0;
}

int
Job::WriteLinkers() const {
  for (auto & [hash, linkers] : fragments) {
    WriteLinkers(hash, linkers);
  }

  for (auto & [hash, linkers] : fragments3) {
    WriteLinkers(hash, linkers);
  }

  return 1;
}

int
Job::WriteLinkers(uint32_t hash,
                  const SetOfLinkers& linkers) const {
  if (database) {
    return WriteLinkersDb(hash, linkers);
  }

  IWString fname;
  fname << stem << '.' << hash << ".txt";
  IWString_and_File_Descriptor output;
  if (! output.open(fname.null_terminated_chars())) {
    cerr << "Job::WriteLinkers:Cannot open '" << fname << "'\n";
    return 0;
  }

  return linkers.Write(output);
}

int
Job::WriteLinkersDb(uint32_t hash,
                    const SetOfLinkers& linkers) const {
  IWString buffer;
  linkers.WriteTextProtos(buffer);

  leveldb::Slice key(reinterpret_cast<const char *>(&hash), sizeof(hash));
  leveldb::Slice value(buffer.data(), buffer.length());

  leveldb::WriteOptions write_options;

  const leveldb::Status status = database->Put(write_options, key, value);
  if (status.ok()) {
    return 1;
  }

  cerr << "Options::DoStore:Did not store " << status.ToString() << "\n";

  return 0;
}

// Each breakable bond has info about the atoms on either side of
// the atoms in the bond.
class BondSlicingInformation {
  private:
    // The atoms that define the bond.
    atom_number_t _a1;
    atom_number_t _a2;
    // The number of atoms associdated with the _a1 side of the bond.
    // Not sure I need this, but since it is computed, why not store it
    int _count1;
    FixedBitVector _from_a1;
    FixedBitVector _from_a2;
  public:
    BondSlicingInformation();

    atom_number_t a1() const {
      return _a1;
    }
    atom_number_t a2() const {
      return _a2;
    }

    int Initialise(Molecule& m, const Bond * b, int * storage);

    int Active() const {
      return _a1 != INVALID_ATOM_NUMBER;
    }

    // Return true if we share an atom with `rhs`.
    int IsAdjacent(const BondSlicingInformation& rhs) const;

    // Identify the atoms that are arranged as
    // a11-a12 ... a22-a21
    int IdentifyAtoms(const BondSlicingInformation& rhs,
                        atom_number_t& a11,
                        atom_number_t& a12,
                        atom_number_t& a21,
                        atom_number_t& a22,
                        FixedBitVector& atoms_between) const;

    // Identify the atoms such that the 2 variant points 'in'.
    //               a11
    //               a12
    //                |
    //  a21-a22  .............a32-a31
    int IdentifyAtoms(const BondSlicingInformation& o1,
                      const BondSlicingInformation& o2,
                        atom_number_t& a11,
                        atom_number_t& a12,
                        atom_number_t& a21,
                        atom_number_t& a22,
                        atom_number_t& a31,
                        atom_number_t& a32,
                        FixedBitVector& atoms_between) const;
};

BondSlicingInformation::BondSlicingInformation() {
  _a1 = INVALID_ATOM_NUMBER;
  _a2 = INVALID_ATOM_NUMBER;
  _count1 = 0;
}

int
BondSlicingInformation::Initialise(Molecule& m,
                        const Bond * b,
                        int * storage) {
  const int matoms = m.natoms();

  std::fill_n(storage, matoms, 0);

  _count1 = m.identify_side_of_bond(storage, b->a1(), 1, b->a2());
#ifdef DEBUG_INITIALISE_BOND_SLICINT_INFO
  cerr << " Atoms " << b->a1() << " to " << b->a2();
  for (int i = 0; i < matoms; ++i) {
    if (storage[i] == 1) {
      cerr << ' ' << i;
    }
  }
  cerr << '\n';
#endif

  // this should never happen, we are only doing chain bonds.
  if (_count1 == 0) {
    return 0;
  }

  _a1 = b->a1();
  _a2 = b->a2();

  _from_a1.ConstructFromArray(storage, matoms, 1);
  _from_a2.ConstructFromArray(storage, matoms, 0);
#ifdef DEBUG_INITIALISE_BOND_SLICINT_INFO
  _from_a1.DebugPrint(cerr);
  _from_a2.DebugPrint(cerr);
#endif

  return 1;
}

// Do two BondSlicingInformation objects share any atoms.
int
BondSlicingInformation::IsAdjacent(const BondSlicingInformation& rhs) const {
  if (_a1 == rhs._a1) {
    return 1;
  }
  if (_a1 == rhs._a2) {
    return 1;
  }
  if (_a2 == rhs._a1) {
    return 1;
  }
  if (_a2 == rhs._a2) {
    return 1;
  }

  return 0;
}

int
SeparatedByOneBond(Molecule& m,
                   const BondSlicingInformation& b1,
                   const BondSlicingInformation& b2) {
  if (m.bonds_between(b1.a1(), b2.a1()) > 3) {
    return 0;
  }
  if (m.bonds_between(b1.a1(), b2.a2()) > 3) {
    return 0;
  }
  if (m.bonds_between(b1.a2(), b2.a1()) > 3) {
    return 0;
  }
  if (m.bonds_between(b1.a2(), b2.a2()) > 3) {
    return 0;
  }

  return 1;
}

// We need to work out how these two bonds are arranged in the molecule.
// Think of heptane
// 0 1 2 3 4 5 6
// The first bond is 1-2 and the second is 4-5.
// We need to identify atoms 2, 3 and 4 as the 'central' atoms.
// We establish the first direction by looking for two sides that
// have no atoms in common.
// We arrange the two bonds into pairs so that the first member points 
// outside the core and the second is part of the core.
// In the case of heptane, a11==1 a12==2 a21==5 a22==4
int
BondSlicingInformation::IdentifyAtoms(const BondSlicingInformation& rhs,
                        atom_number_t& a11,
                        atom_number_t& a12,
                        atom_number_t& a21,
                        atom_number_t& a22,
                        FixedBitVector& atoms_between) const {
  // ...1-2___ 2-1...
  if (_from_a1.BitsInCommon(rhs._from_a1) == 0) {
    a11 = _a1;
    a12 = _a2;
    a21 = rhs._a1;
    a22 = rhs._a2;
#ifdef DEBUG_IDENTIFY_ATOMS_BETWEEN
    _from_a2.DebugPrint(cerr);
    rhs._from_a2.DebugPrint(cerr);
#endif
    atoms_between = _from_a2;
    atoms_between.iwand(rhs._from_a2);
  } else if (_from_a1.BitsInCommon(rhs._from_a2) == 0) {
    a11 = _a1;
    a12 = _a2;
    a21 = rhs._a2;
    a22 = rhs._a1;
#ifdef DEBUG_IDENTIFY_ATOMS_BETWEEN
    _from_a2.DebugPrint(cerr);
    rhs._from_a1.DebugPrint(cerr);
#endif
    atoms_between = _from_a2;
    atoms_between.iwand(rhs._from_a1);
  } else if (_from_a2.BitsInCommon(rhs._from_a1) == 0) {
    a11 = _a2;
    a12 = _a1;
    a21 = rhs._a1;
    a22 = rhs._a2;
#ifdef DEBUG_IDENTIFY_ATOMS_BETWEEN
    _from_a1.DebugPrint(cerr);
    rhs._from_a2.DebugPrint(cerr);
#endif
    atoms_between = _from_a1;
    atoms_between.iwand(rhs._from_a2);
  } else if (_from_a2.BitsInCommon(rhs._from_a2) == 0) {
    a11 = _a2;
    a12 = _a1;
    a21 = rhs._a2;
    a22 = rhs._a1;
#ifdef DEBUG_IDENTIFY_ATOMS_BETWEEN
    _from_a1.DebugPrint(cerr);
    rhs._from_a1.DebugPrint(cerr);
#endif
    atoms_between = _from_a1;
    atoms_between.iwand(rhs._from_a1);
  } else {
    cerr << "BondSlicingInformation this should not happen\n";
    return 0;
  }
#ifdef DEBUG_IDENTIFY_ATOMS_BETWEEN
  cerr << "atoms_between ";
  atoms_between.DebugPrint(cerr);
#endif
  return 1;
}

// There are 3 cut points. Again the objective is to order the atoms
// into pairs with the first member of each pair pointint outside the
// common system, and the second one being part of it.
int
BondSlicingInformation::IdentifyAtoms(const BondSlicingInformation& o2,
                        const BondSlicingInformation& o3,
                        atom_number_t& a11,
                        atom_number_t& a12,
                        atom_number_t& a21,
                        atom_number_t& a22,
                        atom_number_t& a31,
                        atom_number_t& a32,
                        FixedBitVector& atoms_between) const {
  IdentifyAtoms(o2, a11, a12, a21, a22, atoms_between);

#ifdef DEBUG_IDENTIFY_ATOMS_3
  cerr << "atoms";
  cerr << " a11 " << a11 << " a12 " << a12 << " a21 " << a21 << " a22 " << a22 << '\n';
  cerr << "Atoms between are\n";
  atoms_between.DebugPrint(cerr);
  cerr << " o3 involves atoms " << o3._a1 << " and " << o3._a2 << '\n';
  o3._from_a1.DebugPrint(cerr);
  o3._from_a2.DebugPrint(cerr);
#endif
  // If o3 is external to the atoms_between that will form a
  // multi fragment subset. Not of interest.
  const int c1 = atoms_between.BitsInCommon(o3._from_a1);
#ifdef DEBUG_IDENTIFY_ATOMS_3
  cerr << "c1 " << c1 << " c2 " << atoms_between.BitsInCommon(o3._from_a2) << '\n';
#endif
  if (c1 == 0) {
    return 0;
  }
  const int c2 = atoms_between.BitsInCommon(o3._from_a2);
  if (c2 == 0) {
    return 0;
  }
  // We need to figure out if o3._a1 or o3._a2 points into the common
  // area defined by o1 and o2.
  // Get a pointer to the outside of the (o1,o2) region.
  const FixedBitVector* fp;
  if (a12 == _a1) {
    fp = &_from_a2;
  } else {
    fp = &_from_a1;
  }

#ifdef DEBUG_IDENTIFY_ATOMS_3
  cerr << "Bic externa; " << fp->BitsInCommon(o3._from_a1) << " a2 " << fp->BitsInCommon(o3._from_a2) << '\n';
#endif
  if (fp->BitsInCommon(o3._from_a1) == 0) {
    a31 = o3._a1;
    a32 = o3._a2;
    atoms_between.iwand(o3._from_a2);
  } else {
    a31 = o3._a2;
    a32 = o3._a1;
    atoms_between.iwand(o3._from_a1);
  }

  return 1;
}

class MoleculeData {
  private:
    // One per Bond in the molecule.
    BondSlicingInformation* _bond_slicing;

    // Used for constructing subsets.
    std::unique_ptr<int[]> _in_subset;
  // private functions

    int FormPair(Job& options, Molecule& m, int ndx1, int ndx2);
    int FormTriple(Job& options, Molecule& m, int ndx1, int ndx2, int ndx3);
  public:
    MoleculeData();
    ~MoleculeData();

    int Initialise(const Job& options, Molecule& m, const int * breakable);

    // Across all the items in `_bond_slicing` how many are active?
    int NumberBreakableBonds(const Molecule& m) const;

    int FormPairs(Job& options, Molecule& m, const int * bond_symmetry);
    int FormTriples(Job& options, Molecule& m, const int * bond_symmetry);
};

MoleculeData::MoleculeData() {
  _bond_slicing = nullptr;
}

MoleculeData::~MoleculeData() {
  if (_bond_slicing != nullptr) {
    delete [] _bond_slicing;
  }
}

int
MoleculeData::Initialise(const Job& options, Molecule& m,
                const int * breakable) {
  m.ring_membership();
  const int matoms = m.natoms();
  const int nedges = m.nedges();
  std::unique_ptr<int[]> storage = std::make_unique<int[]>(matoms);
  _bond_slicing = new BondSlicingInformation[nedges];
  int bonds_to_be_processed = 0;
#ifdef DEBUG_INITIALISE_MOLECULE_DATA
  cerr << "Molecule contains " << matoms << " and " << nedges << " edges\n";
#endif
  for (int i = 0; i < nedges; ++i) {
    if (! breakable[i]) {
      continue;
    }
    const Bond * b = m.bondi(i);
    if (b->nrings()) {
      continue;
    }
    if (! b->is_single_bond()) {
      continue;
    }
    if (! options.OkAtomicNumber(m.atomic_number(b->a1()))) {
      continue;
    }
    if (! options.OkAtomicNumber(m.atomic_number(b->a2()))) {
      continue;
    }
    if (_bond_slicing[i].Initialise(m, b, storage.get())) {
      ++bonds_to_be_processed;
    }
  }

  _in_subset = std::make_unique<int[]>(matoms);

  return bonds_to_be_processed;
}

int
MoleculeData::FormPairs(Job& options,
                        Molecule& m,
                        const int * bond_symmetry) {
  const int nedges = m.nedges();
  int rc = 0;
  for (int i = 0; i < nedges; ++i) {
    if (bond_symmetry[i] == 0) {
      continue;
    }
    if (! _bond_slicing[i].Active()) {
      continue;
    }
    for (int j = i + 1; j < nedges; ++j) {
      if (! _bond_slicing[j].Active()) {
        continue;
      }
      FormPair(options, m, i, j);
      ++rc;
    }
  }

  return rc;
}

int
MoleculeData::FormTriples(Job& options,
                          Molecule& m,
                          const int * bond_symmetry) {
  const int nedges = m.nedges();
  int rc = 0;
  for (int i = 0; i < nedges; ++i) {
    if (bond_symmetry[i] == 0) {
      continue;
    }
    if (! _bond_slicing[i].Active()) {
      continue;
    }
    for (int j = i + 1; j < nedges; ++j) {
      if (! _bond_slicing[j].Active()) {
        continue;
      }
      for (int k = j + 1; k < nedges; ++k) {
        if (! _bond_slicing[k].Active()) {
          continue;
        }
#ifdef DEBUG_FORM_TRIPLE
        cerr << "Forming triple " << i << j << k << '\n';
#endif
        FormTriple(options, m, i, j, k);
        ++rc;
      }
    }
  }

  return rc;
}

int
MoleculeData::NumberBreakableBonds(const Molecule& m) const {
  const int nedges = m.nedges();
  int rc = 0;
  for (int i = 0; i < nedges; ++i) {
    if (! _bond_slicing[i].Active()) {
      ++rc;
    }
  }

  return rc;
}

//#define DEBUG_FORM_PAIR

// Form a two break fragment by breaking bonds `ndx1` and `ndx2` in `m`.
int
MoleculeData::FormPair(Job& options,
                       Molecule& m,
                       int ndx1,
                       int ndx2) {
  const BondSlicingInformation& b1 = _bond_slicing[ndx1];
  const BondSlicingInformation& b2 = _bond_slicing[ndx2];
#ifdef DEBUG_FORM_PAIR
  cerr << "bonds " << *m.bondi(ndx1) << " and " << *m.bondi(ndx2) << " adjacent " << b1.IsAdjacent(b2) << '\n';
#endif
  if (b1.IsAdjacent(b2)) {
    return 0;
  }
  if (SeparatedByOneBond(m, b1, b2)) {
    return 0;
  }
//cerr << "Bonds " << ndx1 <<  " and " << ndx2 << '\n';

  atom_number_t a11, a12, a21, a22;
  FixedBitVector atoms_between;
  if (! b1.IdentifyAtoms(b2, a11, a12, a21, a22, atoms_between)) {
    return 0;
  }
#ifdef DEBUG_FORM_PAIR
  cerr << "Atoms are " << a11 << ' ' << a12 << ' ' << a21 << ' ' << a22 << '\n';
  atoms_between.DebugPrint(cerr);
#endif
  if (! options.OkSize(atoms_between)) {
    return 0;
  }

  const int distance = m.bonds_between(a12, a21);
#ifdef DEBUG_FORM_PAIR
  cerr << "Attachment points separated by " << distance << " bonds, ok? " << options.OkBondSeparation(distance) << '\n';
#endif
  if (! options.OkBondSeparation(distance)) {
    return 0;
  }

#ifdef DEBUG_FORM_PAIR
  cerr << "Setting isotopes on atoms " << a12 << " and " << a21 << '\n';
#endif
  m.set_isotope(a12, m.atomic_number(a11));
  m.set_isotope(a22, m.atomic_number(a21));
  std::fill_n(_in_subset.get(), m.natoms(), 0);
  atoms_between.Scatter(_in_subset.get(), 1);
#ifdef DEBUG_FORM_PAIR
  for (int i = 0; i < m.natoms(); ++i) {
    if (_in_subset[i]) {
      cerr << i << ' ' << m.smarts_equivalent_for_atom(i) << " in subset\n";
    }
  }
#endif
  Molecule subset;
  m.create_subset(subset, _in_subset.get(), 1);
  m.set_isotope(a12, 0);
  m.set_isotope(a22, 0);
  options.AnotherFragment(2);
  if (! options.MatchesMustHaveQuery(subset)) {
    return 0;
  }
#ifdef DEBUG_FORM_PAIR
  cerr << subset.smiles() << '\n';
#endif

  uint32_t hash = options.HashValue(m.atomic_number(a12),
                                    m.atomic_number(a22),
                                    distance);
  auto iter = options.fragments.find(hash);
  if (iter == options.fragments.end()) {
    auto [ins, notused] = options.fragments.emplace(hash, SetOfLinkers());
    iter = ins;
  }

  iter->second.Extra(options.store_exemplar_smiles, m, subset, distance);

  return 1;
}

int
MoleculeData::FormTriple(Job& options,
                         Molecule& m,
                         int ndx1,
                         int ndx2,
                         int ndx3) {
  const BondSlicingInformation& b1 = _bond_slicing[ndx1];
  const BondSlicingInformation& b2 = _bond_slicing[ndx2];
  const BondSlicingInformation& b3 = _bond_slicing[ndx3];
  if (b1.IsAdjacent(b2) || b1.IsAdjacent(b3) || b2.IsAdjacent(b3)) {
    return 0;
  }
#ifdef DEBUG_FORM_TRIPLE
  cerr << "None are adjacent\n";
#endif
  if (SeparatedByOneBond(m, b1, b2) ||
      SeparatedByOneBond(m, b1, b3) ||
      SeparatedByOneBond(m, b2, b3)) {
    return 0;
  }
#ifdef DEBUG_FORM_TRIPLE
  cerr << "Well separated\n";
#endif

  atom_number_t a11, a12, a21, a22, a31, a32;
  FixedBitVector atoms_between;
  if (! b1.IdentifyAtoms(b2, b3, a11, a12, a21, a22, a31, a32, atoms_between)) {
    return 0;
  }
#ifdef DEBUG_FORM_TRIPLE
  cerr << "Identified atoms\n";
  cerr << "a11 " << a11 << " a12 " << a12 << " a21 " << a21 << " a22 " << a22 << " a31 " << a31 << " a32 " << a32 << '\n';
  atoms_between.DebugPrint(cerr);
#endif

  if (!options.OkBondSeparation(m.bonds_between(a12, a22)) ||
      !options.OkBondSeparation(m.bonds_between(a12, a32)) ||
      !options.OkBondSeparation(m.bonds_between(a22, a32))) {
    return 0;
  }
#ifdef DEBUG_FORM_TRIPLE
  cerr << "Bond separation OK, OkSize " << options.OkSize(atoms_between) << '\n';
#endif
  if (! options.OkSize(atoms_between)) {
    return 0;
  }

  m.set_isotope(a12, m.atomic_number(a11));
  m.set_isotope(a22, m.atomic_number(a21));
  m.set_isotope(a32, m.atomic_number(a31));
  std::fill_n(_in_subset.get(), m.natoms(), 0);
  atoms_between.Scatter(_in_subset.get(), 1);
#ifdef DEBUG_FORM_TRIPLE
  for (int i = 0; i < m.natoms(); ++i) {
    if (_in_subset[i]) {
      cerr << i << ' ' << m.smarts_equivalent_for_atom(i) << " in subset\n";
    }
  }
#endif
  Molecule subset;
  m.create_subset(subset, _in_subset.get(), 1);
  m.set_isotope(a12, 0);
  m.set_isotope(a22, 0);
  m.set_isotope(a32, 0);
  options.AnotherFragment(3);
  if (! options.MatchesMustHaveQuery(subset)) {
    return 0;
  }
#ifdef DEBUG_FORM_TRIPLE
  cerr << subset.smiles() << '\n';
#endif

  const int d12 = m.bonds_between(a12, a22);
  const int d13 = m.bonds_between(a12, a32);
  const int d23 = m.bonds_between(a22, a32);
  uint32_t hash = options.HashValue(m.atomic_number(a12),
                                    m.atomic_number(a22), 
                                    m.atomic_number(a32),
                                    d12, d13, d23);

#ifdef DEBUG_FORM_TRIPLE
  cerr << "Hash value " << hash << '\n';
#endif
  auto iter = options.fragments3.find(hash);
  if (iter == options.fragments3.end()) {
    auto [ins, notused] = options.fragments3.emplace(hash, SetOfLinkers());
    iter = ins;
  }

  iter->second.Extra(options.store_exemplar_smiles, m, subset, d12, d13, d23);

  return 1;
}

#ifdef USE_DATA_
int
SeparatedAtomsMatch2(Molecle& m,
                     job& options,
                     const int * separation_needed,
                     IWString_and_File_Descriptor& output) {
  const int matoms = m.natoms();

  std::unique_ptr<int[]> in_subset(new_int(matoms));

  const int * dm = m.distance_matrix_warning_may_change();

  vertices.resize(options.bond_separations.number_elements());
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      const int dij = dm[i * matoms + j];
      if (separation_needed[dij] == 0) {
        continue;
      }
      --separation_needed[dij];
      for (int k = j + 1; k < matoms; ++k) {
        const int dik = dm[i * matoms + k]
        if (separation_needed[dik] == 0) {
          continue;
        }
      }
      const int djk = dm[j * matoms + k];
      if (separation_needed[djk] == 0) {
        continue;
      }
      // At this stage, atoms i,j,k satisfy geometric constraints.
      IdentifyInnerAtoms(m, i, j, k);
      ++separation_needed[dij];
    }
  }

  return 1;
}

int
SeparatedAtomsMatch1(Molecle& m,
                     job& options,
                     const int * separation_needed,
                     IWString_and_File_Descriptor& output) {
  for (int i = 0; i < matoms; ++i) {
    vertices << i;
    for (int j = i + 1; j < matoms; ++j) {
      const int d = dm[i * matoms + j];
      if (separation_needed[d] == 0) {
        continue;
      }
      vertices << j;
      vertices.resize_keep_storage(1);
    }
    vertices.resize_keep_storage(0);
  }

  return 1;
}

int
SeparatedAtomsMatch(Molecule& m,
                    Job& options,
                    IWString_and_File_Descriptor& output) {
  if (! options.Preprocess(m)) {
    return 0;
  }

  const int * dm = m.distance_matrix_warning_may_change();

  int longest_path = m.longest_path();

  const int matoms = m.natoms();

  int * separation_needed = new_int(matoms); std::unique_ptr<int[]> free_separation_needed(separation_needed);

  if (! options.InitialiseSeparationNeeded(longest_path, separation_needed)) {
    return 1;
  }

  const int number_constraints = options.bond_separations.number_elements();
  switch (number_constraints) {
    case 1:
      return SeparatedAtomsMatch1(m, options, separation_needed, output);
    case 2:
      return SeparatedAtomsMatch2(m, options, separation_needed, output);
    case 3:
      return SeparatedAtomsMatch3(m, options, separation_needed, output);
    default:
      return 0;
  }

  return 1;
}
#endif

int
GenerateDatabase(Job& options,
                 Molecule& m,
                 IWString_and_File_Descriptor& output) {
  if (! options.Preprocess(m)) {
    return 0;
  }

  // First see if there are queries governing what gets matched.
  std::unique_ptr<int[]> breakable = std::make_unique<int[]>(m.nedges());
  if (! options.IdentifyBreakableBonds(m, breakable.get())) {
    return 1;
  }

  MoleculeData mdata;
  if (! mdata.Initialise(options, m, breakable.get())) {
    return 1;
  }

  if (! options.OkBreakableBondCount(mdata.NumberBreakableBonds(m))) {
    return 1;
  }

  std::unique_ptr<int[]> bond_symmetry = BondSymmetry(m);

  if (! mdata.FormPairs(options, m, bond_symmetry.get())) {
    if (options.verbose > 1) {
      cerr << "No pairs in " << m.smiles() << ' ' << m.name() << '\n';
    }
    return 1;
  }

  if (! mdata.FormTriples(options, m, bond_symmetry.get())) {
    if (options.verbose > 1) {
      cerr << "No triples in " << m.smiles() << ' ' << m.name() << '\n';
    }
    return 1;
  }

  return 1;
}

int
GenerateDatabase(Job& options,
                 data_source_and_type<Molecule>& input,
                 IWString_and_File_Descriptor& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);
    ++options.molecules_read;
    if (! GenerateDatabase(options, *m, output)) {
      cerr << "Fatal error processing " << m->smiles() << ' ' << m->name() << '\n';
      return 0;
    }
  }

  return 1;
}

int
GenerateDatabase(Job& options,
                 const char * fname,
                 IWString_and_File_Descriptor& output) {
  FileType input_type = options.input_type;
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose > 1) {
    input.set_verbose(1);
  }

  return GenerateDatabase(options, input, output);
}

void
Usage(int rc) {
  cerr << "Generate fragments by breaking one or more bonds\n";
  cerr << " -B <nbonds>  skip any molecules with more than <nbonds> breakable bonds\n";
  cerr << " -D <nbonds>  maximum bond separation to consider\n";
  cerr << " -f <natoms>  minimum fragment size\n";
  cerr << " -F <natoms>  maximum fragment size\n";
  cerr << " -k <query>   queries specifying bonds to NOT break\n";
  cerr << " -K <query>   queries specifying the only bonds to break\n";
  cerr << " -m <query>   queries that must match in the fragments produced\n";
  cerr << " -S <stem>    file name stem for output files\n";
  cerr << " -Y ...       misc options, entery '-Y help' for info\n";
  cerr << " -c           remove chirality\n";
  cerr << " -l           strip to largest fragment\n";
  cerr << " -v           verbose output\n";

  ::exit(rc);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vcf:F:E:A:i:g:lB:d:S:D:Y:k:K:m:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (! process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process aromaticity options\n";
    return 1;
  }

  Job job_options;
  if (! job_options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! GenerateDatabase(job_options, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (job_options.verbose) {
    job_options.Report(cerr);
  }

  if (! job_options.WriteLinkers()) {
    cerr << "Cannot write linkers\n";
    return 1;
  }

  return 0;
}

}  // namespace separated_atoms

int
main(int argc, char** argv) {
  int rc = separated_atoms::Main(argc, argv);

  return rc;
}
