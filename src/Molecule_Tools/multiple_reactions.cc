//  Apply a set of reactions to input molecules.

#include <iostream>
#include <memory>
#include <optional>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

namespace apply_multiple_reactions {

using std::cerr;

void
Usage(int rc)
{
  cerr << "Apply a set of reactions. Reactions are applied separately\n";
  cerr << " -R <rxn>          specify one or more reactions\n";
  cerr << " -R F:<fname>      file containing multiple reactions\n";
  cerr << " -r <recursion>    number of recursive invocations (def 0)\n";
  cerr << " -z i              ignore molecules not matching any queries\n";
  cerr << " -z w              write molecules not reaction to the output\n";
  cerr << " -S <fname>        file name stem for output\n";
  cerr << " -o <type>         specify output type(s)\n";
  cerr << " -c                remove chirality\n";
  cerr << " -l                strip to largest fragment\n";
  cerr << " -v                verbose output\n";

  ::exit(rc);
}

class MultipleReactions {
 private:

  int _verbose;

  int _molecules_read;

  int _reduce_to_largest_fragment;

  int _remove_chirality;

  resizable_array_p<IWReaction> _rxn;
  //  Keep track of how many times each reaction matches
  extending_resizable_array<int> _rxn_matched;

  //  We keep track of how many of the reactions match each molecule.
  extending_resizable_array<int> _reactions_matching;

  //  Used to keep identify duplicates.
  IW_STL_Hash_Set _seen;

  //  The number of duplicates we identify;
  int _duplicates_discarded;

  int _append_reaction_name_to_product;

  FileType _input_type;

  Chemical_Standardisation _chemical_standardisation;

  //  private functions.

  int ReadReaction(IWString& fname);
  int ReadFileOfReactions(IWString& fname);
  int ReadFileOfReactions(const IWString& fname, iwstring_data_source& input);

  int Process(Molecule& m, Molecule_to_Match& target, IWReaction& rxn,
              resizable_array_p<Molecule>& result);

  int Process(Molecule& m, const Set_of_Atoms& embedding, IWReaction& rxn,
              resizable_array_p<Molecule>& result);

  int IsDuplicate(Molecule& m);
  void MaybeAppendRxnName(const IWString& rxnname, Molecule& product);

 public:

  MultipleReactions();

  int Initialise(Command_Line& cl);

  int Preprocess(Molecule& m);

  int Process(Molecule& m, resizable_array_p<Molecule>& result);

  int Report(std::ostream& output) const;

  FileType input_type() const { return _input_type; }
};

MultipleReactions::MultipleReactions()
{
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _duplicates_discarded = 0;
  _input_type = FILE_TYPE_INVALID;
  _append_reaction_name_to_product = 0;
}

int
MultipleReactions::Initialise(Command_Line& cl)
{
  set_copy_name_in_molecule_copy_constructor(1);

  _verbose = cl.option_count('v');

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove chirality from input molecules\n";
    }
  }

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
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

  if (!cl.option_present('R')) {
    cerr << "Must specify one or more reactions via the -R option\n";
    return 0;
  }

  if (cl.option_present('R')) {
    IWString fname;
    for (int i = 0; cl.value('R', fname, i); ++i) {
      if (fname.starts_with("F:")) {
        fname.remove_leading_chars(2);
        if (!ReadFileOfReactions(fname)) {
          cerr << "Cannot read reactions from '" << fname << "'\n";
          return 0;
        }
      }
      else {
        if (!ReadReaction(fname)) {
          cerr << "Cannot read reaction '" << fname << "'\n";
          return 0;
        }
      }
    }
  }

  if (cl.option_present('a')) {
    _append_reaction_name_to_product = 1;
    if (_verbose) {
      cerr << "Will append the reaction name to the product\n";
    }
  }

  if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) {  //  reading a pipe, assume smiles
    _input_type = FILE_TYPE_SMI;
  }
  else if (!all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern all file types, use the -i option\n";
    return 0;
  }
  else if (!process_input_type(cl, _input_type)) {
    return 0;
  }

  return 1;
}

int
MultipleReactions::ReadFileOfReactions(IWString& fname)
{
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "MultipleReactions::ReadFileOfReactions:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadFileOfReactions(fname, input);
}

int
MultipleReactions::ReadFileOfReactions(const IWString& fname, iwstring_data_source& input)
{
  IWString dirname(fname);
  dirname.truncate_at_last('/');
  dirname << '/';
  IWString buffer;
  while (input.next_record(buffer)) {
    IWString fname;
    fname << dirname << buffer;
    if (!ReadReaction(fname)) {
      cerr << "MultipleReactions::ReadFileOfReactions:Cannot read '" << buffer << "'\n";
      return 0;
    }
  }

  return _rxn.size();
}

int
MultipleReactions::ReadReaction(IWString& fname)
{
  std::optional<ReactionProto::Reaction> maybe_rxn =
      iwmisc::ReadTextProto<ReactionProto::Reaction>(fname);
  if (!maybe_rxn) {
    cerr << "MultipleReactions::ReadTextProto:cannot open '" << fname << "'\n";
    return 0;
  }

  std::unique_ptr<IWReaction> r = std::make_unique<IWReaction>();
  if (!r->ConstructFromProto(*maybe_rxn)) {
    cerr << "MultipleReactions::ReadReaction:cannot parse '" << maybe_rxn->ShortDebugString()
         << '\n';
    return 0;
  }

  _rxn << r.release();

  return 1;
}

int
MultipleReactions::Preprocess(Molecule& m)
{
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
MultipleReactions::Process(Molecule& m, resizable_array_p<Molecule>& results)
{
  ++_molecules_read;

  _seen.insert(m.unique_smiles());

  Molecule_to_Match target(&m);

  int matches = 0;
  const int nrxn = _rxn.number_elements();
  for (int i = 0; i < nrxn; ++i) {
    if (Process(m, target, *_rxn[i], results)) {
      ++_rxn_matched[i];
      ++matches;
    }
  }

  ++_reactions_matching[matches];

  if (matches == 0) {
    if (_verbose > 1) {
      cerr << "MultipleReactions::Process:None of " << _rxn.size() << " reactions reacted with '"
           << m.name() << "'\n";
    }

    return 0;
  }

  return 1;
}

int
MultipleReactions::Process(Molecule& m, Molecule_to_Match& target, IWReaction& rxn,
                           resizable_array_p<Molecule>& result)
{
  Substructure_Results sresults;
  const int nhits = rxn.substructure_search(target, sresults);
  if (_verbose > 2) {
    cerr << nhits << " hits to " << rxn.name() << " in " << m.name() << '\n';
  }

  if (nhits == 0) {
    return 0;
  }

  int rc = 0;
  for (const Set_of_Atoms* e : sresults.embeddings()) {
    if (Process(m, *e, rxn, result)) {
      ++rc;
    }
  }

#ifdef DEBUG_PROCESS
  cerr << "From " << nhits << " rc " << rc << '\n';
#endif
  return rc;
}

int
MultipleReactions::Process(Molecule& m, const Set_of_Atoms& embedding, IWReaction& rxn,
                           resizable_array_p<Molecule>& result)
{
  std::unique_ptr<Molecule> prod = std::make_unique<Molecule>(m);

  if (!rxn.in_place_transformation(*prod.get(), &embedding)) {
    if (_verbose > 2) {
      cerr << "in_place_transformation failed\n";
    }
    return 0;
  }

  if (IsDuplicate(*prod.get())) {
    if (_verbose > 2) {
      cerr << "Duplicat found\n";
    }
    return 0;
  }

  MaybeAppendRxnName(rxn.name(), *prod.get());

  result << prod.release();
  return 1;
}

void
MultipleReactions::MaybeAppendRxnName(const IWString& rxnname, Molecule& product)
{
  if (!_append_reaction_name_to_product) {
    return;
  }

  IWString new_name(product.name());
  new_name << " %% " << rxnname;
  product.set_name(new_name);
}

int
MultipleReactions::IsDuplicate(Molecule& m)
{
  const IWString& usmi = m.unique_smiles();
  const auto iter = _seen.find(usmi);
  //  if we have seen it before, return now.
  if (iter != _seen.end()) {
    ++_duplicates_discarded;
    return 1;
  }

  _seen.insert(usmi);
  //  Have not seen it before.
  return 0;
}

int
MultipleReactions::Report(std::ostream& output) const
{
  output << "MultipleReactions:read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }
  output << "Discarded " << _duplicates_discarded << " duplicates\n";

  for (int i = 0; i < _reactions_matching.number_elements(); ++i) {
    if (_reactions_matching[i]) {
      output << _reactions_matching[i] << " molecules matched " << i << " reactions\n";
    }
  }

  int istop = _rxn_matched.number_elements();
  if (istop > _rxn.number_elements()) {
    istop = _rxn.number_elements();
  }
  for (int i = 0; i < istop; ++i) {
    output << _rxn_matched[i] << " molecules matched " << _rxn[i]->name() << '\n';
  }

  return 1;
}

//  Options for the main program.
struct Options {
  int verbose;

  int write_starting_molecule;

  int ignore_molecules_not_matching_reactions;
  int write_molecules_not_reacting;

  int molecules_read;
  int molecules_not_reacting;

  int recursion_depth;

  extending_resizable_array<int> products_generated;

 public:

  Options();

  int Initialise(Command_Line& cl);

  int Report(std::ostream& output) const;
};

Options::Options()
{
  verbose = 0;
  write_starting_molecule = 0;
  ignore_molecules_not_matching_reactions = 0;
  write_molecules_not_reacting = 0;

  recursion_depth = 0;

  molecules_read = 0;
  molecules_not_reacting = 0;
}

int
Options::Initialise(Command_Line& cl)
{
  verbose = cl.option_present('v');

  if (cl.option_present('p')) {
    write_starting_molecule = 1;
    if (verbose) {
      cerr << "Will write the starting molecule\n";
    }
  }

  if (cl.option_present('z')) {
    IWString z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if (z == 'i') {
        ignore_molecules_not_matching_reactions = 1;
      }
      else if (z == 'w') {
        write_molecules_not_reacting = 1;
      }
      else {
        cerr << "Options::Initialise:unrecognized -z qualifier '" << z << "'\n";
        return 0;
      }
    }
  }

  if (cl.option_present('r')) {
    if (!cl.value('r', recursion_depth) || recursion_depth < 0) {
      cerr << "option_present::Initialise:the recursion depth (-r) option must be a whole number\n";
      return 0;
    }
    if (verbose) {
      cerr << "Will recurse " << recursion_depth << " times\n";
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const
{
  output << "Options:read " << molecules_read << " molecules\n";
  output << molecules_not_reacting << " molecules did not produce results\n";
  for (int i = 0; i < products_generated.number_elements(); ++i) {
    if (products_generated[i]) {
      output << products_generated[i] << " molecules generated " << i << " products\n";
    }
  }

  return 1;
}

void
HandleMoleculesNotMatching(Options& options, Molecule& m,
                           const resizable_array_p<Molecule>& results)
{
  ++options.molecules_not_reacting;
  ++options.products_generated[results.size()];  //  should be zero...
}

int
ApplyMultipleReactions(MultipleReactions& many_reactions, Options& options, Molecule& m,
                       int recursion_depth, resizable_array_p<Molecule>& results)
{
  if (many_reactions.Process(m, results) == 0) {
    HandleMoleculesNotMatching(options, m, results);
    if (recursion_depth > 0) {
      return 1;
    }

    if (!options.ignore_molecules_not_matching_reactions) {
      cerr << "Error processing " << m.name() << '\n';
      return 0;
    }

    return 1;
  }

  ++options.products_generated[results.size()];

  if (recursion_depth < options.recursion_depth) {
    const int existing = results.number_elements();
    for (int i = 0; i < existing; ++i) {
      ApplyMultipleReactions(many_reactions, options, *results[i], recursion_depth + 1, results);
    }
  }

  return 1;
}

//  Write the results of transformations.
int
DoOutput(const Options& options, Molecule& m, resizable_array_p<Molecule>& results,
         Molecule_Output_Object& output)
{
  if (options.write_starting_molecule) {
    output.write(m);
  }

  for (Molecule* r : results) {
    output.write(*r);
  }

  return 1;
}

int
ApplyMultipleReactions(MultipleReactions& many_reactions, Options& options,
                       data_source_and_type<Molecule>& input, Molecule_Output_Object& output)
{
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (!many_reactions.Preprocess(*m)) {
      return 0;
    }

    ++options.molecules_read;

    resizable_array_p<Molecule> results;
    if (!ApplyMultipleReactions(many_reactions, options, *m, 0, results)) {
      return 0;
    }

    DoOutput(options, *m, results, output);
  }

  return 1;
}

int
ApplyMultipleReactions(MultipleReactions& many_reactions, Options& options, const char* fname,
                       Molecule_Output_Object& output)
{
  FileType input_type = many_reactions.input_type();
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ApplyMultipleReactions(many_reactions, options, input, output);
}

int
Main(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vE:A:i:g:lcR:S:az:pr:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process aromaticity options\n";
    return 1;
  }

  if (!process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process standard elements options (-E)\n";
    return 1;
  }

  MultipleReactions many_reactions;
  if (!many_reactions.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  Options options;
  if (!options.Initialise(cl)) {
    Usage(1);
  }

  Molecule_Output_Object output;
  if (cl.option_present('o')) {
    if (!output.determine_output_types(cl, 'o')) {
      cerr << "Cannot determing output type(s) (-o)\n";
      return 1;
    }
  }
  else {
    output.add_output_type(FILE_TYPE_SMI);
  }

  if (!cl.option_present('S')) {
    cerr << "Must specify output file name stem via the -S option\n";
    return 1;
  }

  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    if (output.would_overwrite_input_files(cl, fname)) {
      cerr << "Cannot overwrite input file(s) '" << fname << "'\n";
      return 1;
    }

    if (!output.new_stem(fname)) {
      cerr << "Cannot open output file (-S) '" << fname << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Output written to '" << fname << "'\n";
    }
  }

  for (const char* fname : cl) {
    if (!ApplyMultipleReactions(many_reactions, options, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
    many_reactions.Report(cerr);
    cerr << "Write " << output.molecules_written() << " molecules\n";
  }

  return 0;
}

}  //  namespace apply_multiple_reactions

int
main(int argc, char** argv)
{
  int rc = apply_multiple_reactions::Main(argc, argv);

  return rc;
}
