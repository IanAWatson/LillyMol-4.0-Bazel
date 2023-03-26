// Generate a set of ReplacementRing protos from a set of molecules.
// Protos are written as text_format since there are never that many of them.

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <memory>
#include <optional>
#include <vector>

#include "google/protobuf/text_format.h"

#define IWQSORT_FO_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Duplicate/replacement_ring.pb.h"

namespace ring_extraction {

using std::cerr;

constexpr char kOpenSquareBracket = '[';
constexpr char kCloseSquareBracket = ']';

void
Usage(int rc) {
  cerr << "Extracts rings and ring systems creating ReplacementRing protos that can be used by ring_replacement\n";
  cerr << " -S <stem>         create ring data files with <stem>\n";
  cerr << " -R <rsize>        max ring size to process (def 7)\n";
  cerr << " -k                also generate smarts with connectivity not specified\n";
  cerr << " -x                transform within ring aliphatic double bonds to type any\n";
  cerr << " -P <atype>        label atoms by atom type of exocyclic attached atom\n";
  cerr << " -X ...            more options\n";
  cerr << " -c                remove chirality\n";
  cerr << " -l                strip to largest fragment\n";
  cerr << " -v                verbose output\n";

  ::exit(rc);
}

class ExtractRings {
  private:
    int _verbose;

    int _molecules_read;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    // We can optionally mark the attachment point.
    isotope_t _isotope;

    // We can also have the isotope mean the atom type of what used to be
    // attached. Note that is creates ambiguity if there were multiple atom
    // types attached to the ring at a given atom. For now we ignore those rings.
    // TODO:ianwatson revisit this maybe.
    Atom_Typing_Specification _atype;

    // We can ignore rings that are too large.
    uint32_t _max_ring_size;

    // We can ignore ring systems containing too many rings.
    uint32_t _max_ring_system_size;

    int _transform_ring_double_bonds;

    // For every ring/system type we have a mapping from unique
    // smiles to protos.
    IW_STL_Hash_Map<IWString, IW_STL_Hash_Map<IWString, RplRing::ReplacementRing>> _ring;

    // Outputs will be written to a set of files with prefix `_stem`.
    IWString _stem;

    // As a check, once the smarts for a ring is generated, we can do a search
    // in the starting molecule.
    int _substructure_search_starting_molecule;
    int _substructure_search_failures;

    // Optionally we can generate raw rings, with no substituent information.
    int _generate_substitution_not_specified;

    FileType _input_type;

    Chemical_Standardisation _chemical_standardisation;

  // Private functions.
    int GenerateRing(Molecule& parent, Molecule& m, const IWString& label, const int* include_atom, int include_d);
    int GenerateSmarts(Molecule& m, const int* include_atom, int include_d, IWString& result) const;

    int LabelAttachmentPoints(Molecule& parent,
                              Molecule& r,
                              const int* ring_sys,
                              int sys_num,
                              const int* xref,
                              std::unique_ptr<uint32_t[]>& atypes) const;
    isotope_t IsotopeOfExocyclicAtom(Molecule& m,
                atom_number_t zatom,
                const int* ring_sys,
                int sys_num,
                const std::unique_ptr<uint32_t[]>& atypes) const;

    std::optional<IWString> CanonicalRingName(Molecule& m,
                        const int * ring_sys,
                        int sys) const;

    void ChangeRingDoubleBonds(Molecule& m,
                      IWString& smt) const;
    int MaybeCheckSubstructureMatch(Molecule& m, const IWString& smt);

    int WriteRings(IWString& fname,
                         const IW_STL_Hash_Map<IWString, RplRing::ReplacementRing>& rings) const;
    int WriteRings(IWString_and_File_Descriptor& output,
                         const IW_STL_Hash_Map<IWString, RplRing::ReplacementRing>& rings) const;
    isotope_t IsotopeForAtom(Molecule& m, atom_number_t zatom,
                             const int* ring_sys,
                             int sys_num,
                             const std::unique_ptr<uint32_t[]>& atypes) const;


  public:
    ExtractRings();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    // Extract the rings from `m` and accumulate in our internal data structures.
    int Process(Molecule& m);

    // At the completion of processing, report summary information.
    int Report(std::ostream& output) const;

    // Once data is assembled, write protos.
    int WriteRings() const;

    FileType input_type() const {
      return _input_type;
    }
};

ExtractRings::ExtractRings() {
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _isotope = 1;
  _max_ring_size = 7;
  _max_ring_system_size = 3;
  _transform_ring_double_bonds = 0;
  _substructure_search_starting_molecule = 0;
  _substructure_search_failures = 0;
  _generate_substitution_not_specified = 0;
  _input_type = FILE_TYPE_INVALID;
}

void
DisplayDashXOption(std::ostream& output) {
  output << " -X sss     substructure search the starting molecule as a check\n";

  ::exit(0);
}

int
ExtractRings::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove chirality from input molecules\n";
    }
  }

  if (cl.option_present('g')) {
    if (! _chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
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

  if (! cl.option_present('S')) {
    cerr << "Must specify output file name stem via the -S option\n";
    return 0;
  }

  if (cl.option_present('S')) {
    cl.value('S', _stem);
    if (_verbose) {
      cerr << "Output files created with name step '" << _stem << "'\n";
    }
  }

  if (cl.option_present('R')) {
    if (! cl.value('R', _max_ring_size) || _max_ring_size < 3) {
      cerr << "The max ring size (-R) option must be a valid ring size\n";
      return 0;
    }
    if (_verbose) {
      cerr << "WIll skip rings with more than " << _max_ring_size << " atoms\n";
    }
  }

  if (cl.option_present('x')) {
    _transform_ring_double_bonds = 1;
    if (_verbose) {
      cerr << "Will transform within ring aliphatic double bonds to type any\n";
    }
  }

  if (cl.option_present('k')) {
    _generate_substitution_not_specified = 1;
    if (_verbose) {
      cerr << "Will also generate a smarts variant without connectivity\n";
    }
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');
    if (! _atype.build(p)) {
      cerr << "Invalid atom type '" << p << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Atom typing initialised '" << p << "'\n";
    }
  }

  if (cl.option_present('X')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if (x == "sss") {
        _substructure_search_starting_molecule = 1;
        if (_verbose) {
          cerr << "Will perform a substructure search vs the starting molecule\n";
        }
      } else if (x == "help") {
        DisplayDashXOption(cerr);
      } else {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        DisplayDashXOption(cerr);
      }
    }
  }

  if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) { // reading a pipe, assume smiles
    _input_type = FILE_TYPE_SMI;
  } else if (cl.size() == 1 && strncmp(cl[0], "F:", 2) == 2) {
  } else if (process_input_type(cl, _input_type)) {
  } else if (all_files_recognised_by_suffix(cl)) {
  } else {
    cerr << "Cannot discern all file types, use the -i option\n";
    return 0;
  }

  return 1;
}

int
ExtractRings::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (m.empty()) {
    return 0;
  }

  // Do not process phosphorus.
  if (m.natoms(15) > 0) {
    return 0;
  }

  return 1;
}

int
ExtractRings::Report(std::ostream& output) const {
  output << "ExtractRings:read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }

  if (_substructure_search_starting_molecule) {
    output << _substructure_search_failures << " substructure search failures\n";
  }

  return 1;
}

// For each atom, unset implicit hydrogen information.
void
UnsetImplicitHydrogenInformation(Molecule& m) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    m.unset_all_implicit_hydrogen_information(i);
  }
}

// within `m`, `ring_sys` designates a set of atoms comprising a
// ring system. Extend that to any doubly bonded extensions to the
// ring system.
int
ExtendToDoublyBonded(Molecule& m,
                int * ring_sys) {
  // Force sssr if needed.
  m.ring_membership();

  Set_of_Atoms added_here;
  resizable_array<int> sys;

  int rc = 0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (ring_sys[i] == 0) {
      continue;
    }
    if (m.ncon(i) < 3) {
      continue;
    }

    const Atom& a = m.atom(i);
    if (a.nbonds() == a.ncon()) {
      continue;
    }
    for (const Bond* b : a) {
      if (b->is_single_bond()) {
        continue;
      }
      if (b->nrings()) {
        continue;
      }

      atom_number_t o = b->other(i);
      if (ring_sys[o] == ring_sys[i]) {
        continue;
      }

      if (ring_sys[o] > 0) {
        continue;
      }
      added_here << o;
      sys << ring_sys[i];
      ++rc;
    }
  }

  for (int i = 0; i < added_here.number_elements(); ++i) {
    atom_number_t atom = added_here[i];
    ring_sys[atom] = sys[i];
  }

  return rc;
}

int
ExtractRings::Process(Molecule& m) {
  ++_molecules_read;

  if (m.nrings() == 0) {
    return 1;
  }

  const int matoms = m.natoms();
  std::unique_ptr<int[]> ring_sys = std::make_unique<int[]>(matoms);
  std::unique_ptr<int[]> xref = std::make_unique<int[]>(matoms);
  std::unique_ptr<int[]> include_atom(new_int(matoms, 1));

  std::unique_ptr<uint32_t[]> atypes;
  if (_atype.active()) {
    atypes.reset(new uint32_t[matoms]);
    _atype.assign_atom_types(m, atypes.get());
  }

  m.label_atoms_by_ring_system_including_spiro_fused(ring_sys.get());
  ExtendToDoublyBonded(m, ring_sys.get());
#ifdef DEBUG_PROCESS
  write_isotopically_labelled_smiles(m, false, cerr);
  cerr << '\n';
  for (int i = 0; i < m.natoms(); ++i) {
    cerr << " i " << i << " ring_sys " << ring_sys[i] << '\n';
  }
#endif

  for (int sys = 1; ; ++sys) {
    Molecule r;
    if (! m.create_subset(r, ring_sys.get(), sys, xref.get())) {
      break;
    }

    UnsetImplicitHydrogenInformation(r);

    if (! LabelAttachmentPoints(m, r, ring_sys.get(), sys, xref.get(), atypes)) {
      continue;
    }

    std::optional<IWString> label = CanonicalRingName(m, ring_sys.get(), sys);
    if (! label) {
      continue;
    }
    // cerr << "Label '" << label << "' " << r.unique_smiles() << '\n';

    r.set_name(m.name());

    GenerateRing(m, r, *label, include_atom.get(), 1);
    //Maybe also compute the variant with no substition information.
    if (_generate_substitution_not_specified) {
      r.transform_to_non_isotopic_form();
      GenerateRing(m, r, *label, include_atom.get(), 0);
    }
  }

  return 1;
}

// Atom `zatom` is part of a ring system where ring_sys[zatom] == sys_num.
// If it is bonded to an exocyclic atom, return the `atypes` value for
// that attached atom. If there are no attached atoms outside the ring
// system, return 0;
isotope_t
ExtractRings::IsotopeOfExocyclicAtom(Molecule& m,
                atom_number_t zatom,
                const int* ring_sys,
                int sys_num,
                const std::unique_ptr<uint32_t[]>& atypes) const {
  const Atom& atom = m.atom(zatom);
  // If 2 connected, no exocyclic bonds.
  if (atom.ncon() == 2) {
    return 0;
  }

  for (const Bond* b : atom) {
    atom_number_t j = b->other(zatom);
    if (ring_sys[j] == sys_num) {
      continue;
    }

    return atypes[j];
  }

  // All attached atoms are part of the ring system.
  return 0;
}

// We are applying some kind of isotope to an atom. If we have atom types
// use that, otherwise _isotope.
isotope_t
ExtractRings::IsotopeForAtom(Molecule& m, atom_number_t zatom, 
                             const int* ring_sys,
                             int sys_num,
                             const std::unique_ptr<uint32_t[]>& atypes) const {
//cerr << "IsotopeForAtom " << atom_number << " atypes " << atypes.get() << " value " << atypes[atom_number] << '\n';

  if (atypes) {
    return IsotopeOfExocyclicAtom(m, zatom, ring_sys, sys_num, atypes);
  }
  
  return _isotope;
}

// Return 1 if there any atoms where ring_sys[i] == sys_num that correspond
// to more than chain bonds attached to the ring atom.
int
AnyMultiplyConnectedAtoms(Molecule& m,
                const int* ring_sys,
                int sys_num) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (ring_sys[i] != sys_num) {
      continue;
    }

    if (m.is_aromatic(i)) {
      continue;
    }

    const Atom& atom = m.atom(i);
    if (atom.ncon() <= 3) {
      continue;
    }

    // There can be just 1 connection that is not in a ring.
    if (m.ring_bond_count(i) + 1 != atom.ncon()) {
      return 1;
    }
  }

  return 0;  // none detected.
}

int
ExtractRings::LabelAttachmentPoints(Molecule& parent,
                                    Molecule& r,
                                    const int* ring_sys,
                                    int sys_num,
                                    const int* xref,
                                    std::unique_ptr<uint32_t[]>& atypes) const {
  const int matoms = parent.natoms();

  // If we are applying isotopes for the attachment point, make sure we do not have
  // any ring atoms attached to more than 1 exocyclic chain atom.
  if (atypes) {
    if (AnyMultiplyConnectedAtoms(parent, ring_sys, sys_num)) {
      return 0;
    }
  }
  // Return the number of isotopes applied.
  int rc = 0;
  //cerr << "LabelAttachmentPoints::processing ring " << sys_num << '\n';
  for (int i = 0; i < matoms; ++i) {
    //cerr << "parent atom " << i << " ring_sys " << ring_sys[i] << '\n';
    if (ring_sys[i] != sys_num) {
      continue;
    }
    const Atom& a = parent.atom(i);
    //cerr << "  acon " << a.ncon() << '\n';
    if (a.ncon() < 3) {
      continue;
    }
    for (const Bond* b : a) {
      atom_number_t o = b->other(i);
      //cerr << "   attached to " << o << " ring_sys " << ring_sys[o] << '\n';
      if (ring_sys[o] != ring_sys[i]) {
        r.set_isotope(xref[i], IsotopeForAtom(parent, i, ring_sys, sys_num, atypes));
        ++rc;
        //cerr << "ring set isotope on atom " << xref[i] << " value " << _isotope << '\n';
      }
    }
  }

  return rc;
}

// When building a canonical label, information needed about the
// component rings in a ring system.
struct RType {
  int rsize;
  int aromatic;

  RType() {
    rsize = 0;
    aromatic = 0;
  }
  RType(int s, int a) {
    rsize = s;
    aromatic = a;
  }
};

IWString
operator<< (IWString& buffer, const RType& rtype) {
  buffer << rtype.rsize;
  if (rtype.aromatic) {
    buffer << 'a';
  } else {
    buffer << 'A';
  }

  return buffer;
}

// For sorting ring types.
//First on ring size then aromaticity
struct
CompareRType {
  int operator()(const RType& rt1, const RType& rt2) const {
    if (rt1.rsize < rt2.rsize) {
      return -1;
    }
    if (rt1.rsize > rt2.rsize) {
      return 1;
    }
    if (rt1.aromatic > rt2.aromatic) {
      return -1;
    }
    if (rt1.aromatic < rt2.aromatic) {
      return 1;
    }
    return 0;
  }
};

// Maybe return a canonical name for the ring system defined
// by the atoms `ring_sys[i] == sys`.
std::optional<IWString>
ExtractRings::CanonicalRingName(Molecule& m,
                        const int * ring_sys,
                        int sys) const {
  m.compute_aromaticity_if_needed();

  // Gather the rings in the system.
  std::vector<RType> rtype;
  for (const Ring* r : m.sssr_rings()) {
    if (r->count_members_set_in_array(ring_sys, sys) == 0) {
      continue;
    }
    if (r->size() > _max_ring_size) {
      return std::nullopt;
    }
    rtype.emplace_back(RType(r->number_elements(), r->is_aromatic()));
  }

  if (rtype.size() > 1) {
    if (rtype.size() > _max_ring_system_size) {
      return std::nullopt;
    }

    static CompareRType cmp;
    iwqsort(rtype.data(), rtype.size(), cmp);
  }

  IWString result;
  for (const RType& r: rtype) {
    result << r;
  }

  return result;
}

int
GetPosition(const resizable_array<atom_number_t>& order, int n) {
  return order.index(n);
}

// A smarts `smt` has been formed, if we are doing substructure checks, return
// true if `smt` matches `m`.
int
ExtractRings::MaybeCheckSubstructureMatch(Molecule& m, const IWString& smt) {
  if (! _substructure_search_starting_molecule) {
    return 1;
  }

  Substructure_Query query;
  if (! query.create_from_smarts(smt)) {
    cerr << "ExtractRings::MaybeCheckSubstructureMatch:invalid smarts '" << smt << "'\n";
    return 0;
  }

  Molecule_to_Match target(&m);
  Substructure_Results sresults;
  if (query.substructure_search(target, sresults)) {
    return 1;
  }

  cerr << "No substructure match " << m.smiles() << " smt '" << smt << "' matched " <<
      sresults.max_query_atoms_matched_in_search() << " query atoms\n";
  cerr << "Contains " << m.aromatic_atom_count() << " aromatic atoms\n";
  write_isotopically_labelled_smiles(m, false, cerr);
  cerr << '\n';
  ++_substructure_search_failures;

  return 0;
}

// `smt` is a unique smiles for `m`. If there are ring bonds in `smt`
// that are type double, change to type any.
void
ExtractRings::ChangeRingDoubleBonds(Molecule& m,
                      IWString& smt) const {
  if (! _transform_ring_double_bonds) {
    return;
  }

  if (! smt.contains('=')) {
    return;
  }

  const resizable_array<atom_number_t> & atom_order_in_smiles = m.atom_order_in_smiles();

  IWString new_smt;
  new_smt.reserve(smt.size());

  int atom_number = -1;
  for (int i = 0; i < smt.number_elements(); ++i) {
    const char c = smt[i];
    if (c == kOpenSquareBracket) {
      ++atom_number;
      new_smt << c;
      continue;
    }
    if (smt[i] != '=') {
      new_smt << c;
      continue;
    }
    int previous_atom = atom_order_in_smiles.index(atom_number);
    if (m.ncon(previous_atom) == 1) {
      new_smt << '=';
      continue;
    }
    int next_atom = atom_order_in_smiles.index(atom_number + 1);
    // Trailing ring closures will not have a next atom
    if (next_atom < 0) {
      new_smt << '=';
      continue;
    }

    if (m.is_aromatic(previous_atom) && m.is_aromatic(next_atom)) {
      new_smt << ':';
      continue;
    }
    if (m.ncon(next_atom) == 1) {
      new_smt << '=';
      continue;
    }
    if (m.in_same_ring(previous_atom, next_atom)) {
      new_smt << '~';
    }
  }
  // cerr << "Convert " << smt << " to " << new_smt << '\n';
  smt = new_smt;
}

// Given a subset of atoms in `m` indicated by `include_atom`,
// generate a smarts.
// If `include_d` is set, each atomic smarts will include the D directive,
// either as a fixed value or as a D> directive.
int
ExtractRings::GenerateSmarts(Molecule& m,
               const int * include_atom,
               int include_d,
               IWString& result) const {

  const int matoms = m.natoms();

  Smiles_Information smiles_information(matoms);
  smiles_information.allocate_user_specified_atomic_smarts();
  smiles_information.set_smiles_is_smarts(1);

  for (int i = 0; i < matoms; ++i) {
    if (! include_atom[i]) {
      continue;
    }

    IWString smt;
    smt << kOpenSquareBracket;
    // this is not quite correct. If the atom is exocyclic and aliphatic here
    // it prevents matching an aromatic later. Ignore for now.
    if (m.is_aromatic(i)) {
      smt << 'a';
    } else {
      smt << 'A';
    }
    if (m.ring_bond_count(i)) {
      smt << 'x' << m.ring_bond_count(i);
    }
    // This is inefficient, we could precomute the string ring membership(s) for each atom.
    for (const Ring* r : m.sssr_rings()) {
      if (r->contains(i)) {
        smt << 'r' << r->size();
      }
    }

    if (! include_d) {
    } else if (_atype.active() && m.isotope(i) > 0) {
      smt << "D>" << m.ncon(i);
    }else if (_isotope == m.isotope(i)) {
      smt << "D>" << m.ncon(i);
    } else if (m.ncon(i) > 1) {
      smt << "D" << m.ncon(i);
    }
    smt << kCloseSquareBracket;
    smiles_information.set_user_specified_atomic_smarts(i, smt);
  }

  set_write_smiles_aromatic_bonds_as_colons(1);
  result = m.smiles(smiles_information, include_atom);
  set_write_smiles_aromatic_bonds_as_colons(0);

  ChangeRingDoubleBonds(m, result);

  return 1;
}

// `m` is a subset of `parent`, governed by `include_atom`.
// Create smiles and smarts from `m` and add to `_hash`.
// Currently exocyclic double bonds are not handled properly. Ignoring for now.
int
ExtractRings::GenerateRing(Molecule& parent,
                           Molecule& m, const IWString& label, const int* include_atom,
                           int include_d) {

  const IWString smi = m.smiles();
  const IWString& usmi = m.unique_smiles();

  IWString smt;
  GenerateSmarts(m, include_atom, include_d, smt);

  MaybeCheckSubstructureMatch(parent, smt);

  auto iter_label = _ring.find(label);
  if (iter_label == _ring.end()) {
    IW_STL_Hash_Map<IWString, RplRing::ReplacementRing> r;
    auto [iter, _] = _ring.emplace(std::make_pair(label, std::move(r)));
    iter_label = iter;
  }

  auto iter_usmi = iter_label->second.find(usmi);
  if (iter_usmi == iter_label->second.end()) {
    RplRing::ReplacementRing r;
    r.set_smi(smi.data(), smi.length());
    r.set_smt(smt.data(), smt.length());
    r.set_id(m.name().AsString());
    r.set_usmi(usmi.data(), usmi.length());
    r.set_conn(include_d);
    r.set_n(1);

    iter_label->second.emplace(std::make_pair(usmi, std::move(r)));
  } else {
    const auto n = iter_usmi->second.n();
    iter_usmi->second.set_n(n + 1);
  }
  
  return 1;
}

// For each stored ring write the protos to a file
// with prefix `_stem`.
int
ExtractRings::WriteRings() const {
  for (const auto& [label, rings] : _ring) {
    IWString fname;
    fname << _stem << '_' << label << ".smi";
    if (! WriteRings(fname, rings)) {
      cerr << "ExtractRings::WriteRings:cannot write '" << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

// Write a set of rings `rings` to `fname` as text_proto form.
int
ExtractRings::WriteRings(IWString& fname,
                         const IW_STL_Hash_Map<IWString, RplRing::ReplacementRing>& rings) const {
  IWString_and_File_Descriptor output;
  if (! output.open(fname.null_terminated_chars())) {
    cerr << "ExtractRings::WriteRings:cannot open '" << fname << "'\n";
    return 0;
  }

  return WriteRings(output, rings);
}

// Write a set of rings `rings` to `output` as text_proto form.
int
ExtractRings::WriteRings(IWString_and_File_Descriptor& output,
                         const IW_STL_Hash_Map<IWString, RplRing::ReplacementRing>& rings) const {
  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  std::string buffer;

  for (const auto& [usmi, r] : rings) {
    printer.PrintToString(r, &buffer);
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
ReplaceRings(ExtractRings& extract_rings,
            Molecule& m,
            IWString_and_File_Descriptor& output) {
  return extract_rings.Process(m);
}

int
ReplaceRings(ExtractRings& extract_rings,
            data_source_and_type<Molecule>& input,
            IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! extract_rings.Preprocess(*m)) {
      continue;
    }

    if (! ReplaceRings(extract_rings, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
ReplaceRings(ExtractRings& extract_rings,
            const char * fname,
            IWString_and_File_Descriptor& output) {
  FileType input_type = extract_rings.input_type();
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ReplaceRings(extract_rings, input, output);
}


int
ReplaceRings(ExtractRings& extract_rings,
            iwstring_data_source& input,
            IWString_and_File_Descriptor& output) {
  IWString fname;
  while (input.next_record(fname)) {
    if (! ReplaceRings(extract_rings, fname.null_terminated_chars(), output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:g:lcS:xR:P:");
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

  ExtractRings extract_rings;
  if (! extract_rings.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  for (const char* fname : cl) {
    IWString tmp(fname);
    if (tmp.starts_with("F:")) {
      tmp.remove_leading_chars(2);
      iwstring_data_source input;
      if (! input.open(tmp.null_terminated_chars())) {
        cerr << "ReplaceRings:cannot open '" << tmp << "'\n";
        return 0;
      }
      if (! ReplaceRings(extract_rings, input, output)) {
        cerr << "Fatal error processing '" << fname << "'\n";
        return 1;
      }

      continue;
    }

    if (! ReplaceRings(extract_rings, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  extract_rings.WriteRings();

  if (verbose) {
    extract_rings.Report(cerr);
  }

  return 0;
}

}  // namespace ring_extraction

int
main(int argc, char** argv) {
  int rc = ring_extraction::Main(argc, argv);

  return rc;
}
