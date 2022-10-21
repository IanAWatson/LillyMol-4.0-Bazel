#include <iostream>

#include "Molecule_Tools/unique_molecules_api.h"

namespace unique_molecules {

using std::cerr;

UniqueMoleculesImplementation::UniqueMoleculesImplementation() {
  _verbose = 0;
  _molecules_processed = 0;
  _duplicates_found = 0;
  _unique_molecules = 0;
  _reduce_to_largest_fragment = 0;
  _exclude_chiral_info = 0;
  _exclude_cis_trans_bonding_info = 0;
  _ignore_isotopes = 0;

  _use_atom_hash = 0;
  _perform_formula_check = 0;
  _compare_as_graph = 0;
  _only_same_if_structure_and_name_the_same = 0;

  _molecules_changed_by_reactions = 0;

  _accumulate_unique_smiles_counts = 0;

  return;
}

int
UniqueMoleculesImplementation::Report(std::ostream& output) const {
  output << _molecules_processed << " molecules read, ";
  output << _duplicates_found << " duplicates, ";
  output << _unique_molecules << " unique structures\n";

  if (_use_atom_hash) {
    unsigned int s = 0;
    for (auto & h : _atom_hash) {
      s += h->size();
    }
    cerr << "Atom hash contains " << s << " discrete values\n";
  }

  return output.good();
}

// No checking, no NOTHING!
// Beware if reaction queries overlap!!
// Be very careful!
int
UniqueMoleculesImplementation::PerformReactions(Molecule& m) {
  int rc = 0;

  for (IWReaction* rxn : _reaction) {
    Substructure_Results sresults;

    int nhits = rxn->determine_matched_atoms(m, sresults);
    if (nhits == 0) {
      continue;
    }

    if (_verbose > 2) {
      cerr << nhits << " hits for reaction " << rxn->comment() << '\n';
    }

    rxn->in_place_transformations(m, sresults);

    ++rc;
  }

  if (rc && _verbose > 1) {
    cerr << m.name() << ", " << rc << " reactions changed\n";
  }

  if (rc) {
    _molecules_changed_by_reactions++;
  }

  return rc;
  return 1;
}

int
UniqueMoleculesImplementation::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (_elements_to_remove.active()) {
    _elements_to_remove.process(m);
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  if (_elements_to_remove.active()) {
    _elements_to_remove.process(m);
  }

  if (! _reaction.empty()) {
    PerformReactions(m);
  }

  if (_ignore_isotopes) {
    m.transform_to_non_isotopic_form();
  }

  return 1;
}

int
UniqueMoleculesImplementation::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce input molecules to largest fragment before usmi generation\n";
    }
  }

  if (cl.option_present('I')) {
    _ignore_isotopes = 1;
    if (_verbose) {
      cerr << "Will transform molecules to non isotopic form\n";
    }
  }

  if (cl.option_present('z')) {
    _exclude_cis_trans_bonding_info = 1;
    if (_verbose) {
      cerr << "Will exclude cis trans bonding from unique smiles generation\n";
    }
  }

  if (cl.option_present('c')) {
    _exclude_chiral_info = 1;
    if (_verbose) {
      cerr << "Will exclude chirality information from unique smiles generation\n";
    }
  }

  if (cl.option_present('h')) {
    _use_atom_hash = 1;
    if (_verbose) {
      cerr << "Will use an atom hash from quick_bond_hash()\n";
    }
  }

  if (cl.option_present('a')) {
    _compare_as_graph = 1;
    if (_verbose) {
      cerr << "Will compare graph forms\n";
    }
  }

  if (cl.option_present('j')) {
    _only_same_if_structure_and_name_the_same = 1;
    if (_verbose)
      cerr << "Molecules will be considered identical only if both structure and name match\n";
  }

  if (cl.option_present('R')) {
    Sidechain_Match_Conditions sidechain_match_conditions;

    if (! read_reactions(cl, _reaction, sidechain_match_conditions, 'R')) {
      cerr << "Cannot read reaction(s) (-R option)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Defined " << _reaction.size() << " reactions\n";
    }

    for (IWReaction* r : _reaction) {
      r->set_find_unique_embeddings_only(1);
    }
  }

  if (cl.option_present('g')) {
    if (! _chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisations (-g option)\n";
      return 0;
    }
  }

  if (cl.option_present('X')) {
    if (! _elements_to_remove.construct_from_command_line(cl, _verbose, 'X'))
    {
      cerr << "Cannot discern elements to remove, -X\n";
      return 0;
    }
  }

  if (cl.option_present('T')) {
    if (! _element_transformations.construct_from_command_line(cl, _verbose, 'T')) {
      cerr << "Cannot process element transformations (-t option)\n";
      return 0;
    }
  }

  if (cl.option_present('e')) {
    _accumulate_unique_smiles_counts = 1;
    if (_verbose) {
      cerr << "Unique smiles and counts accumulated\n";
    }
  }

  for (int i = 0; i < 200; i++) {
    _smiles_hash.add(new UMHash());

    if (_perform_formula_check) {
      _formula_hash.add(new UMHash());
    }
    if (_use_atom_hash) {
      _atom_hash.add(new std::unordered_set<uint64_t>());
    }
  }

  return 1;
}

int
UniqueMoleculesImplementation::ActivateChemicalStandardization(const std::string& directive) {
  const IWString s(directive.data(), directive.size());
  return _chemical_standardisation.Activate(s, _verbose);
}

// No special provision for the empty molecule.
int
UniqueMoleculesImplementation::IsUnique(Molecule& m) {
  ++_molecules_processed;

  if (IsUniqueInner(m)) {
    ++_unique_molecules;
    return 1;
  } else {
    ++_duplicates_found;
    return 0;
  }
}

int
UniqueMoleculesImplementation::IngestPreviousMolecule(Molecule& m) {
  return IsUniqueInner(m);
}

// Form the key that will be used for the hash lookup.
void
UniqueMoleculesImplementation::FormUniqueSmilesInner(Molecule& m,
                        IWString& usmi) const {
  if (_compare_as_graph) {
    const IWString formula = m.molecular_formula();

    Molecule g(m);
    g.change_to_graph_form();

    usmi = g.unique_smiles();

    usmi << ':' << formula;
  } else {
    usmi = m.unique_smiles();
  }

  if (_only_same_if_structure_and_name_the_same) {
    usmi << m.name();
  }
  //cerr << "Testing unique smiles '" << usmi << "'\n";
}

void
UniqueMoleculesImplementation::FormUniqueSmiles(Molecule& m, IWString& usmi) const {
  const auto save_chiral = include_chiral_info_in_smiles();
  const auto save_cis_trans = include_cis_trans_in_smiles();

  if (_exclude_chiral_info) {
    set_include_chiral_info_in_smiles(0);
  }

  if (_exclude_cis_trans_bonding_info) {
    set_include_cis_trans_in_smiles(0);
  }

  FormUniqueSmilesInner(m, usmi);

  set_include_chiral_info_in_smiles(save_chiral);
  set_include_cis_trans_in_smiles(save_cis_trans);
}

#ifdef FORMULA_CHECK_SLOWS_THINGS_DOWN
// This never worked, seems it is too expensive to compute and check the formula.
// quick_bond_hash works better.
    IWString mformula;
    m.formula_distinguishing_aromatic(mformula);

    //cerr << m.name() << " usmi " << usmi << "'\n";

    IW_STL_Hash_Set* s = _formula_hash[matoms];

    if (s->end() == s->find(mformula)) {
      s->insert(mformula);

      IW_STL_Hash_Map_int* h =
          smiles_hash[matoms];    // very important to update the smiles hash too
      (*h)[usmi] = 1;

      unique_molecules++;
      return 1;
    }
#endif

int
UniqueMoleculesImplementation::IsUniqueInner(Molecule& m) {
  Preprocess(m);

  IWString usmi;  // The key that will be used for hash lookup.
  FormUniqueSmiles(m, usmi);

  if (_compare_as_graph) {
    const IWString formula = m.molecular_formula();

    Molecule g(m);
    g.change_to_graph_form();

    usmi = g.unique_smiles();

    usmi << ':' << formula;
  } else {
    usmi = m.unique_smiles();
  }

  //cerr << "Testing unique smiles '" << usmi << "'\n";

  if (_only_same_if_structure_and_name_the_same) {
    usmi << m.name();
  }

  const int matoms = m.natoms();

  while (matoms >= _smiles_hash.number_elements()) {
    _smiles_hash.add(new UMHash);
    if (_use_atom_hash) {
      _atom_hash.add(new std::unordered_set<uint64_t>());
    }
  }

  //cerr << matoms << " smiles_hash " << smiles_hash.number_elements() << " formula_hash " << formula_hash.number_elements() << '\n';

  if (_use_atom_hash) {
    const uint64_t h = m.quick_bond_hash();

    std::unordered_set<uint64_t>* ha = _atom_hash[matoms];

    const auto f = ha->find(h);

    // Never seen before, molecule is unique.
    if (f == ha->end()) {
      ha->insert(h);

      UMHash* h = _smiles_hash[matoms];
      h->emplace(usmi);

      ++_unique_molecules;
      return 1;
    }
  }

#ifdef FORMULA_CHECK_SLOWS_THINGS_DOWN
  if (_perform_formula_check) {
    PerformFormulaCheck(usmi);  // not implemented.
  }
#endif

  UMHash* h = _smiles_hash[matoms];

  auto f = h->find(usmi);

  MaybeUpdateGlobalUsmiCounters(usmi);

  if (f == h->end()) {
    h->insert(usmi);
    return 1;
  } else {
    return 0;
  }
}

// If we have been requested to accumulate _usmi_count, add `usmi`.
int
UniqueMoleculesImplementation::MaybeUpdateGlobalUsmiCounters(const IWString& usmi) {
  if (! _accumulate_unique_smiles_counts) {
    return 0;
  }

  auto iter = _usmi_count.find(usmi);
  if (iter != _usmi_count.end()) {
    ++iter->second;
    return 1;
  }

  _usmi_count.emplace(usmi, 1);
  return 1;
}

int
UniqueMoleculesImplementation::WriteUsmiHash(const char* fname) const {
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "UniqueMoleculesImplementation::WriteUsmiHash:cannot open '" << fname << "'\n";
    return 0;
  }

  return WriteUsmiHash(output);
}

int
UniqueMoleculesImplementation::WriteUsmiHash(IWString_and_File_Descriptor& output) const {
  for (const auto& [usmi, count] : _usmi_count) {
    output << usmi << ' ' << count << '\n';
    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

// We should have a better way of communicating back a smiles failure.
int
UniqueMoleculesImplementation::SmilesIsUnique(const std::string& smiles) {
  Molecule m;
  if (! m.build_from_smiles(smiles)) {
    cerr << "UniqueMoleculesImplementation::IsUnique:invalid smiles '" << smiles << "'\n";
    return 0;
  }

  return IsUnique(m);
}

}  // namespace unique_molecules
