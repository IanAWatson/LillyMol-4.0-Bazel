/*
  Programme for converting creating a query file from a structure file
*/

#include <iostream>
#include <fstream>
#include <memory>
#include <optional>
#include <tuple>

#include "google/protobuf/text_format.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/ematch.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/mdl_molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

namespace mol2qry {

const char * prog_name = nullptr;

using std::cerr;

int verbose = 0;

int queries_written = 0;

/*
  We may get molecules from slicer. These have isotopic labels that indicate
  where the fragment is attached to the rest of the molecule
*/

int isotopically_labelled_from_slicer = 0;

Chemical_Standardisation chemical_standardisation;

IWString stem_for_output;

/*
  If we are processing multiple molecules, we need to produce a different .qry
  file for each
*/

int next_file_name_to_produce = 0;

std::ofstream stream_for_names_of_query_files;

/*
  We can recognise R atoms as substitution points
*/

int change_R_groups_to_substitutions = 0;

Element_Matcher rgroup;

/*
  People sometimes draw a two atom molecule across a ring in order to
  signify that every atom on that ring could be a point of substitution
*/

//static int kludge_for_ring_substitution = 0;

int all_queries_in_one_file = 0;

std::ofstream stream_for_all_queries;

int remove_chiral_centres = 0;

int add_explicit_hydrogens = 0;

Substructure_Query coordination_point;

int radius_from_coordination_point = 0;

int remove_isotopes_from_input_molecules = 0;

IWString append_to_comment;

int perform_matching_test = 0;

// We can write serialized protos to a TFDataRecord file.
std::unique_ptr<iw_tf_data_record::TFDataWriter> proto_destination;

// Or write textformat protos to a file.
int write_as_text_proto = 0;

// It can be convenient to append the smiles to the name - just to avoid
// having to join it up later.
int append_smiles_to_query_name = 0;

void
usage(int rc = 1)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
  cerr << "Converts a molecule to a query file\n";
  cerr << "  -m             all ncon and nbonds values are written as minima\n";
  cerr << "  -r             all ring bonds become type ANY\n";
  cerr << "  -j             all atoms conserve their ring membership\n";
  cerr << "  -n             all non ring bonds are marked as \"nrings 0\"\n";
  cerr << "  -a             only aromatic atoms will match aromatic atoms\n";
  cerr << "  -d             the saturated/unsaturated status of atoms will be preserved\n";
  cerr << "  -M <smiles>    specify smiles directly\n";
  cerr << "  -s             only allow substitutions at     isotopically labelled atoms\n";
  cerr << "  -w             only allow substitutions at NON isotopically labelled atoms\n";
  cerr << "  -t             not all isotopic atoms need be substituted\n";
  cerr << "  -c             the isotopic number is the number of extra connections at that atom\n";
  cerr << "  -k             use preference values to resolve symmetric atoms\n";
  cerr << "  -u <smarts>    smarts to specify embedding points\n";
  cerr << "  -f ele=smarts  atoms with element type <ele> should match only\n";
  cerr << "  -x <iso>       atoms with isotope <iso> means match any atom type\n";
  cerr << "  -h             condense explicit hydrogens to hcount directives on their anchor atoms\n";
  cerr << "                 atoms matching <smarts>\n";
  cerr << "  -R <rx>        atoms of type <rx> specify substitution points\n";
  cerr << "                 <rx> is a regular expression, e.g. '^R[0-9]*$', or just 'R'\n";
  cerr << "  -o             remove chirality information from molecules\n";
  cerr << "  -L <smarts>    specify atoms that bind to external group\n";
  cerr << "  -l <nbonds>    include all atoms within <nbonds> of the -L atom(s)\n";
  cerr << "  -I             only include isotopically labelled atoms in the query\n";
  cerr << "  -e             query file contains just element type and connectivity info\n";
  cerr << "  -C <fname>     file containing MoleculeToQuery proto for building the query\n";
  cerr << "  -V <file>      file containing environment specification\n";
  cerr << "  -X <file>      file containing environment_no_match specification\n";
  cerr << "  -F <fname>     create a file containing the names of all the query files\n";
  cerr << "  -i <type>      specify input file type\n";
  cerr << "  -S <fname>     specify output file name\n";
  cerr << "  -P <fname>     name for serialized proto output\n";
  cerr << "  -b             put all queries in a single file rather than separate file for each\n";
  cerr << "  -p             write text_format proto queries\n";
  cerr << "  -D ...         create proto query files with GeometricConstraints\n";
  cerr << "  -Y ...         more obscure options, enter '-Y help' for info\n";
  display_standard_aromaticity_options(cerr);
  cerr << "  -g ...         chemical standardisation options, enter '-g help' for info\n";
  cerr << "  -v             verbose operation\n";
  
  exit(rc);
}

struct GeometryConfig {
  // Has anything been specified.
  bool active = false;
  // by how much can existing distances vary and still match.
  float tolerance = 0.10;
  // Are Hydrogen atoms included or not?
  bool include_hydrogen = false;
  // Perhaps a Substructure_Query to specify the atoms to consider.
};

int
BuildGeometricConfig(Command_Line & cl, char flag,
                     GeometryConfig& geometry_config) {
  IWString d;
  for (int i = 0; cl.value(flag, d, i); ++i) {
    d.to_lowercase();
    if (d.starts_with("tol=")) {
      d.remove_leading_chars(4);
      if (! d.numeric_value(geometry_config.tolerance) || 
            geometry_config.tolerance < 0.0) {
        cerr << "BuildGeometricConfig:invalid tolerance " << d << '\n';
        return 0;
      }
    } else if (d == "inch") {
      geometry_config.include_hydrogen = true;
    } else {
      cerr << "Unrecognised -" << flag << " specification '" << d << "'\n";
      return 0;
    }
  }

  geometry_config.active = true;
  return 1;
}

std::tuple<float, float>
DistanceRange(float distance, const GeometryConfig& config) {
  const float delta = config.tolerance * distance;
  const float min_dist = std::max(0.0f, distance - delta);
  const float max_dist = distance + delta;
  return {min_dist, max_dist};
}

int
ToGeometricConstraints(MDL_Molecule& m,
                       const IWString& name_stem,
                       const GeometryConfig& config) {
  IWString smarts;
  const int matoms = m.natoms();
  std::unique_ptr<int[]> ecount(new_int(HIGHEST_ATOMIC_NUMBER + 1));
  std::unique_ptr<int[]> atom_xref(new_int(matoms, -1));

  int ndx = 0;    // A count of the number of active atoms.
  for (int i = 0; i < matoms; ++i) {
    const Atom * a = m.atomi(i);
    const atomic_number_t atnum = a->atomic_number();
    if (atnum == 0 || atnum > HIGHEST_ATOMIC_NUMBER) {
      cerr << "ToGeometricConstraints:invalid atomic number " << atnum << '\n';
      return 0;
    }
    if (atnum == 1 && ! config.include_hydrogen) {
      continue;
    }
    atom_xref[i] = ndx;
    ndx++;
    IWString atomic_smarts;
    atomic_smarts << "[#" << atnum << ']';
    smarts.append_with_spacer(atomic_smarts, '.');
    ecount[atnum]++;
  }

  SubstructureSearch::SubstructureQuery proto;
  SubstructureSearch::SingleSubstructureQuery * qry = proto.add_query();
  std::string tmp(smarts.data(), smarts.length());
  qry->set_smarts(tmp);
  if (append_smiles_to_query_name) {
    IWString name;
    name << m.name() << ' ' << m.smiles();
    tmp.assign(name.data(), name.length());
    qry->set_comment(tmp);
  } else {
    tmp.assign(m.name().data(), m.name().length());
    qry->set_comment(tmp);
  }

  for (int i = 1; i <= HIGHEST_ATOMIC_NUMBER; ++i) {
    if (ecount[i] == 0) {
      continue;
    }
    SubstructureSearch::ElementsNeeded * needed = qry->add_elements_needed();
    needed->set_atomic_number(i);
    needed->add_hits_needed(ecount[i]);
  }

  GeometricConstraints::SetOfConstraints * constraints = qry->add_geometric_constraints();

  for (int i = 0; i < matoms; ++i) {
    if (atom_xref[i] < 0) {
      continue;
    }
    for (int j = i + 1; j < matoms; ++j) {
      if (atom_xref[j] < 0) {
        continue;
      }
      const float d = m.distance_between_atoms(i, j);
      GeometricConstraints::Distance * dconstraint = constraints->add_distances();
      dconstraint->set_a1(atom_xref[i]);
      dconstraint->set_a2(atom_xref[j]);
      GeometricConstraints::Range * range = dconstraint->mutable_range();
      const auto [min_dist, max_dist] = DistanceRange(d, config);
      range->set_min(min_dist);
      range->set_max(max_dist);
    }
  }

  std::string as_string;
  google::protobuf::TextFormat::PrintToString(proto, &as_string);

  IWString fname(name_stem);
  fname << "proto";

  std::ofstream output(fname.null_terminated_chars(), std::ios::out);
  output << as_string;
  output.close();

  return 1;
}

int
expand_isotopes(MDL_Molecule & m,
                atom_number_t zatom,
                int radius,
                isotope_t iso)
{
  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

//cerr << "Expanding isotopes from " << zatom << '\n';

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (iso == m.isotope(j))
      continue;

    m.set_isotope(j, iso);

    if (radius > 0)
      expand_isotopes(m, j, radius - 1, iso);
  }

  return 1;
}

int
identify_coordination_point_and_adjacent_atoms(MDL_Molecule & m)
{
  Substructure_Results sresults;

  Molecule_to_Match target(&m);

  const int nhits = coordination_point.substructure_search(target, sresults);

  if (0 == nhits)
  {
    cerr << "Zero hits to coordination point substructure search\n";
    return 0;
  }

  if (verbose)
    cerr << m.name() << " " << nhits << " hits to coordination point query\n";

  for (const Set_of_Atoms* e : sresults.embeddings())
  {
    const atom_number_t j = e->item(0);

    MDL_Atom_Data * mdlad = m.mdl_atom_data(j);

    mdlad->set_substitution(-2);   // means exactly as specified

    m.set_isotope(j, 973);

    expand_isotopes(m, j, radius_from_coordination_point, 973);
  }

  return nhits;
}

int
mol2qry(MDL_Molecule & m,
        Molecule_to_Query_Specifications & mqs,
        std::ostream & output)
{
  Set_of_Atoms & substitution_points = mqs.externally_specified_substitution_points();

  substitution_points.resize_keep_storage(0);

  if (change_R_groups_to_substitutions)
    m.change_R_groups_to_substitutions(rgroup, 0);

  if (radius_from_coordination_point <= 0)
    ;
  else if (identify_coordination_point_and_adjacent_atoms(m))
    mqs.set_only_include_isotopically_labeled_atoms(1);
  else
  {
    cerr << "Cannot identify coordination points in '" << m.name() << "'\n";
    return 0;
  }

  Substructure_Query query;
  if (! query.create_from_molecule(m, mqs)) {   // it inherits the molecule name
    cerr << "cannot create query from molecule '" << m.name() << "'\n";
    return 1;
  }

  if (perform_matching_test) {
    cerr << "Performing matching test\n";
    if (0 == query.substructure_search(&m))
    {
      cerr << "No match to searching origin '" << m.name() << "'\n";
      return 0;
    }
  }
 
  if (append_to_comment.length()) {
    IWString tmp(m.name());
    tmp.append_with_spacer(append_to_comment);
    query[0]->set_comment(tmp);
  }

  if (append_smiles_to_query_name) {
    IWString tmp(m.name());
    tmp << ' ' << m.smiles();
    query[0]->set_comment(tmp);
  }

  if (proto_destination) {
    SubstructureSearch::SubstructureQuery proto = query.BuildProto();
    std::string serialized;
    proto.SerializeToString(&serialized);
    return proto_destination->Write(serialized.data(), serialized.size());
  }

  if (write_as_text_proto) {
    SubstructureSearch::SubstructureQuery proto = query.BuildProto();
    // A query comment gets in the way of reading these files.
    proto.clear_comment();

    std::string data;
    if (! google::protobuf::TextFormat::PrintToString(proto, &data)) {
      std::cerr << "mol2qry:cannot write\n";
      return 0;
    }
    output << data;
    return 1;
  }

  if (! query.write_msi(output)) {
    return 0;
  }

  queries_written++;

  return output.good();
}

int
mol2qry(MDL_Molecule & m,
        Molecule_to_Query_Specifications & mqs,
        const GeometryConfig& geometry_config,
        const IWString & output_stem)
{
  if (isotopically_labelled_from_slicer && 0 == m.number_isotopic_atoms())
    cerr << "Warning, only substitute at isotopically labelled atoms, but no isotopes '" << m.name() << "'\n";

  if (all_queries_in_one_file)
    return mol2qry(m, mqs, stream_for_all_queries);

  // If writing protos to a binary file, the output stream is not used, so pass anything.
  if (proto_destination) {
    return mol2qry(m, mqs, std::cout);
  }

  IWString output_fname(output_stem);

  output_fname << next_file_name_to_produce;

  next_file_name_to_produce++;

  output_fname += '.';

  if (geometry_config.active) {
    return ToGeometricConstraints(m, output_fname, geometry_config);
  }

  output_fname += suffix_for_file_type(FILE_TYPE_QRY);

  std::ofstream output(output_fname.null_terminated_chars(), std::ios::out);

  if (! output.good())
  {
    cerr << "Cannot open output file '" << output_fname << "'\n";
    return 0;
  }

  if (stream_for_names_of_query_files.rdbuf()->is_open())
    stream_for_names_of_query_files << output_fname << '\n';

  return mol2qry(m, mqs, output);
}

/*
  Note that there is the potential for serious problems if
  the molecule is changed by the chemical standardisation,
  and the MDL_File_Data component of the MDL_Molecule gets
  out of sync...
*/

void
preprocess(MDL_Molecule & m)
{
  if (remove_isotopes_from_input_molecules)
    m.transform_to_non_isotopic_form();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  if (remove_chiral_centres)
    m.remove_all_chiral_centres();

  if (add_explicit_hydrogens)
    m.make_implicit_hydrogens_explicit();

  return;
}

int
mol2qry(data_source_and_type<MDL_Molecule> & input,
        Molecule_to_Query_Specifications & mqs,
        const GeometryConfig& geometry_config,
        IWString & output_fname)
{
  MDL_Molecule * m;

  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<MDL_Molecule> free_m(m);

    preprocess(*m);

    if (! m->arrays_allocated())
      m->build(*m);

    if (! mol2qry(*m, mqs, geometry_config, output_fname))
      return 0;
  }

  return 1;
}

int
mol2qry(const char * input_fname,
        FileType input_type,
        Molecule_to_Query_Specifications & mqs,
        const GeometryConfig& geometry_config,
        IWString & output_fname)
{
  if (FILE_TYPE_INVALID == input_type)
  {
    input_type = discern_file_type_from_name(input_fname);
    assert (FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<MDL_Molecule> input(input_type, input_fname);
  if (! input.ok())
  {
    cerr << prog_name << ": cannot read '" << input_fname << "'\n";
    return 1;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return mol2qry(input, mqs, geometry_config, output_fname);
}

int
mol2qry(const char * ifile,
        const FileType input_type,
        Molecule_to_Query_Specifications & mqs,
        const GeometryConfig& geometry_config)
{
  IWString output_fname;

  if (all_queries_in_one_file) {  // file already opened elsewhere
  } else if (stem_for_output.length()) {
    output_fname = stem_for_output;
  //} else if (proto_destination) {
  } else {
    output_fname = ifile;
    output_fname.remove_suffix();
  }

  return mol2qry(ifile, input_type, mqs, geometry_config, output_fname);
}

int
do_read_environment(const const_IWSubstring & fname,
                    Molecule_to_Query_Specifications & mqs)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return mqs.read_environment_specification(input);
}

int
do_read_environment_no_match(const const_IWSubstring & fname,
                             Molecule_to_Query_Specifications & mqs)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return mqs.read_environment_no_match_specification(input);
}


int
process_smiles_from_command_line(const IWString & smiles,
                                 Molecule_to_Query_Specifications & mqs,
                                 const GeometryConfig& geometry_config)
{
  MDL_Molecule m;
  if (! m.build_from_smiles(smiles))
  {
    cerr << "Cannot parse -M smiles '" << smiles << "'\n";
    return 54;
  }

  if (0 == stem_for_output.length())
    return mol2qry(m, mqs, std::cout);

  return mol2qry(m, mqs, geometry_config, stem_for_output);
}

void
display_dash_y_options(std::ostream & os)
{
  os << " -Y minextra=n  for a match, target must have at least N extra atoms\n";
  os << " -Y maxextra=n  for a match, target must have at most  N extra atoms\n";
  os << " -Y APPC=<s>    append <s> to the comment field of all queries produced\n";
  os << " -Y exph        add explicit hydrogens, but construct query so anything matched\n";
  os << " -Y ablk        aromatic bonds lose their kekule identity\n";
  os << " -Y minfm=<f>   set the min fraction atoms matched to <f>\n";
  os << " -Y maxfm=<f>   set the max fraction atoms matched to <f>\n";
  os << " -Y A2A=<f>     set aromatic atom translation\n";
  os << " -Y A2A=1       aromatic atoms become 'aromatic'\n";
  os << " -Y A2A=2       aromatic heteroatoms must match aromatic heteroatoms\n";
  os << " -Y A2A=3       aromatic rings must preserve the number of heteroatoms\n";
  os << " -Y rmiso       remove all isotope information from input molecules\n";
  os << " -Y ncon=n      matches must have exactly  <n> connections to unmatched atoms\n";
  os << " -Y min_ncon=n  matches must have at least <n> connections to unmatched atoms\n";
  os << " -Y max_ncon=n  matches must have at most  <n> connections to unmatched atoms\n";
  os << " -Y qnamsmi     append the smiles to the query name\n";
  os << " -Y test        for each query formed, do a match against the starting molecule\n";

  exit(1);
}

int
mol2qry(int  argc, char ** argv) {
  Command_Line cl(argc, argv, "aA:S:P:nrmvE:i:M:sV:X:F:f:R:btg:heu:ojK:Y:kl:L:IcdD:x:pC:");

  verbose = cl.option_count('v');

  if (cl.unrecognised_options_encountered()) {
    usage(2);
  }

  if (! process_elements(cl)) {
    usage(3);
  }

  if (cl.option_present('g')) {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('K')) {
    if (! process_standard_smiles_options(cl, verbose, 'K')) {
      cerr << "Cannot initialise standard smiles options (-K)\n";
      return 4;
    }
  }
  
  if (cl.option_present('m') && cl.option_present('R'))
  {
    cerr << "Sorry, the -m and -R options are mutually incompatible, contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)\n";
    return 3;
  }

// Historical quirk. When I wrote this, the -R option meant regular expression.
// Then in May 2005, I needed to allow both regular expressions and element matches.
// The Element_Matcher object can do that, but for it to process a regular expression,
// the string must start with 'RX='

  if (cl.option_present('R'))
  {
    const_IWSubstring r = cl.string_value('R');
    IWString tmp;

    if (r.starts_with("EMATCH:"))
    {
      r.remove_leading_chars(7);
      tmp = r;
    }
    else
      tmp << "RX=" << r;

    if (! rgroup.construct_from_string(tmp))
    {
      cerr << "Invalid R group matching specification '" << tmp << "'\n";
      return 4;
    }

    change_R_groups_to_substitutions = 1;

    if (verbose)
    {
      cerr << "R groups will be changed to substution point specifications\n";
      rgroup.debug_print(cerr);
    }
  }

  if (cl.option_present('b')) {
    if (cl.option_present('F')) {
      cerr << "The -F and -b options don't make sense together\n";
      usage(3);
    }

    if (! cl.option_present('S')) {
      cerr << "Sorry, must specify the -S option with the -b option\n";
      usage(5);
    }

    all_queries_in_one_file = 1;

    if (verbose) {
      cerr << "All queries written to a single file\n";
    }
  }

  if (cl.option_present('S'))
  {
    cl.value('S', stem_for_output);
    if (verbose)
      cerr << "Stem for output is '" << stem_for_output << "'\n";

    if (cl.number_elements() > 1) {
      cerr << "When specifying a stem, only one file can be processed\n";
      return 1;
    }

    if (all_queries_in_one_file) {
      stem_for_output << ".qry";

      stream_for_all_queries.open(stem_for_output.null_terminated_chars(), std::ios::out);
      if (! stream_for_all_queries.good()) {
        cerr << "Cannot open '" << stem_for_output << "'\n";
        return 5;
      }

      if (verbose)
        cerr << "All queries written to '" << stem_for_output << "'\n";
    }
  }

  if (cl.option_present('P')) {
    IWString fname = cl.string_value('P');
    proto_destination = std::make_unique<iw_tf_data_record::TFDataWriter>();
    if (! proto_destination->Open(fname)) {
      cerr << "Cannot open binary serialized proto file " << fname << '\n';
      return 1;
    }

    if (verbose) {
      cerr << "serialized protos written to " << fname << "'\n";
    }
  }

  if (cl.option_present('p')) {
    write_as_text_proto = 1;
    if (verbose) {
      cerr << "Will write queries as textproto form\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('M') && cl.option_present('i')) {
    cerr << "The -M and -i options are mutually exclusive\n";
    usage(11);
  }

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot parse -i directives\n";
      usage(16);
    }
  }

  Molecule_to_Query_Specifications mqs;

  if (cl.option_present('C')) {
    IWString fname = cl.string_value('C');
    std::optional<molecule_to_query::MoleculeToQuery> proto =
                iwmisc::ReadTextProto<molecule_to_query::MoleculeToQuery>(fname);
    if (! proto) {
      cerr << "Cannot read proto '" << fname << "'\n";
      return 1;
    }

    if (! mqs.Build(*proto)) {
      cerr << "Cannot initialise Molecule_to_Query_Specifications from " << proto->ShortDebugString() << '\n';
    }
  }

  if (cl.option_present('m')) {
    mqs.set_make_embedding(1);
  }

  if (cl.option_present('r')) {
    mqs.set_all_ring_bonds_become_undefined(1);
  }

  if (cl.option_present('n')) {
    mqs.set_non_ring_atoms_become_nrings_0(1);
  }
  if (cl.option_present('j')) {
    mqs.set_atoms_conserve_ring_membership(1);
  }

  if (cl.option_present('h')) {
    mqs.set_condense_explicit_hydrogens_to_anchor_atoms(1);

    if (verbose)
      cerr << "Will merge explicit hydrogen information into anchor atom(s)\n";
  }

  if (cl.option_present('e'))
  {
    mqs.set_just_atomic_number_and_connectivity(1);

    if (verbose)
      cerr << "Queries will contain just atomic number and connectivity info\n";
  }

  if (cl.option_present('s') && cl.option_present('w'))
  {
    cerr << "The -s and -w options are mutually incompatible\n";
    usage(3);
  }

  if (cl.option_present('s') || cl.option_present('c') || cl.option_present('t'))
  {
//  mqs.substitutions_only_at().create_from_smarts("[!0*]");
    isotopically_labelled_from_slicer = 1;

    mqs.set_substituents_only_at_isotopic_atoms(1);

    if (cl.option_present('t'))
    {
      mqs.set_must_have_substituent_at_every_isotopic_atom(0);
      if (verbose)
        cerr << "Not all isotopically labelled atoms need substituents\n";
    }

    if (cl.option_present('c')) {
      mqs.set_isotope_count_means_extra_connections(1);
      if (verbose)
        cerr << "Isotopic number indicates number of extra connections\n";
    }
  }
  else if (cl.option_present('w')) {
    mqs.set_substituents_only_at_non_isotopic_atoms(1);
  }
  else if (cl.option_present('u')) {
    const_IWSubstring smarts;
    cl.value('u', smarts);

    if (! mqs.substitutions_only_at().create_from_smarts(smarts))
    {
      cerr << "Invalid smarts for substitution point(s) '" << smarts << "'\n";
      return 3;
    }
  }

  if (cl.option_present('f')) {
    int i = 0;
    const_IWSubstring f;
    while (cl.value('f', f, i++))
    {
      if (! mqs.set_smarts_for_atom(f))
      {
        cerr << "Invalid smarts for atom '" << f << "'\n";
        return 0;
      }
    }
  }

  if (cl.option_present('a')) {
    mqs.set_only_aromatic_atoms_match_aromatic_atoms(1);
    if (verbose)
      cerr << "Only aromatic atoms will match aromatic atoms\n";
  }

  if (cl.option_present('d')) {
    mqs.set_preserve_saturation(1);
    if (verbose)
      cerr << "Atom saturation will be preserved\n";
  }

  if (cl.option_present('V')) {
    const_IWSubstring v = cl.string_value('V');

    if (! do_read_environment(v, mqs))
    {
      cerr << "Cannot read query environment specification from '" << v << "'\n";
      return 8;
    }

    if (verbose)
      cerr << "Read query environment specification from '" << v << "'\n";
  }

  if (cl.option_present('X')) {
    const_IWSubstring x = cl.string_value('X');

    if (! do_read_environment_no_match(x, mqs))
    {
      cerr << "Cannot read query environment rejection specification from '" << x << "'\n";
      return 8;
    }

    if (verbose)
      cerr << "Read query environment rejection specification from '" << x << "'\n";
  }

  if (cl.option_present('k')) {
    mqs.set_use_preference_values_to_distinguish_symmetry(1);

    if (verbose)
      cerr << "Query atom preference values used to differentiate queries\n";
  }

  if (cl.option_present('o')) {
    remove_chiral_centres = 1;
    if (verbose)
      cerr << "Chiral centres will be removed from input molecules\n";
  }

  if (cl.option_present('L')) {
    if (! cl.option_present('l'))
    {
      cerr << "When specifying a coordination point (-L) must also specify bond radius (-l)\n";
      usage(3);
    }

    if (! cl.value('l', radius_from_coordination_point) || radius_from_coordination_point < 1)
    {
      cerr << "The radius from coordination point option (-l) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will include all atoms within " << radius_from_coordination_point << " bonds of coordination point\n";

    const const_IWSubstring l = cl.string_value('L');

    if (! coordination_point.create_from_smarts(l))
    {
      cerr << "Invalid coordination point smarts '" << l << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Coordination points defined by matches to '" << l << "'\n";

    mqs.set_only_include_isotopically_labeled_atoms(1);
  }

  if (cl.option_present('I')) {
    mqs.set_only_include_isotopically_labeled_atoms(1);

    if (verbose)
      cerr << "Will only include isotopically labelled atoms in the query\n";
  }

  if (cl.option_present('x')) {
    int x;
    if (! cl.value('x', x) || x < 0) {
      cerr << "The -x option (any atom matched) must be a whole +ve number\n";
      return 1;
    }
    mqs.set_isotope_means_match_any_atom(x);
    if (verbose) {
      cerr << "Atoms with isotope " << x << " become a signal for match any atom type\n";
    }
  }

  if (cl.option_present('F')) {
    const char * f = cl.option_value('F');

    stream_for_names_of_query_files.open(f, std::ios::out);

    if (! stream_for_names_of_query_files.good())
    {
      cerr << "Cannot open stream for query files '" << f << "'\n";
      return 8;
    }

    if (verbose)
      cerr << "Query file names written to '" << f << "'\n";
  }

  if (cl.option_present('Y')) {
    int i = 0;
    const_IWSubstring y;

    while (cl.value('Y', y, i++))
    {
      if (y.starts_with("minextra="))
      {
        y.remove_leading_chars(9);
        int e;
        if (! y.numeric_value(e) || e < 0)
        {
          cerr << "The min number extra atoms to be matched '-Y minextra=' must be a whole +ve number\n";
          display_dash_y_options(cerr);
        }

        mqs.set_min_extra_atoms_in_target(e);

        if (verbose)
          cerr << "Matches require at least " << e << " extra atoms\n";
      }
      else if (y.starts_with("maxextra="))
      {
        y.remove_leading_chars(9);
        int e;
        if (! y.numeric_value(e) || e < 0)
        {
          cerr << "The max number extra atoms to be matched '-Y minextra=' must be a whole +ve number\n";
          display_dash_y_options(cerr);
        }

        mqs.set_max_extra_atoms_in_target(e);

        if (verbose)
          cerr << "Matches require at most " << e << " extra atoms\n";
      }
      else if (y.starts_with("ncon="))
      {
        y.remove_leading_chars(5);
        int n;
        if (! y.numeric_value(n) || n < 0)
        {
          cerr << "The number of connections to matched atoms '-Y ncon=' must be a whole +ve number\n";
          display_dash_y_options(cerr);
        }

        mqs.set_ncon(n);

        if (verbose)
          cerr << "Matches can have only " << n << " connections to unmatched atoms\n";
      }
      else if (y.starts_with("min_ncon="))
      {
        y.remove_leading_chars(9);
        int n;
        if (! y.numeric_value(n) || n < 0)
        {
          cerr << "The minimum number of connections to matched atoms '-Y min_ncon=' must be a whole +ve number\n";
          display_dash_y_options(cerr);
        }

        mqs.set_min_ncon(n);

        if (verbose)
          cerr << "Matches must have at least " << n << " connections to unmatched atoms\n";
      }
      else if (y.starts_with("max_ncon="))
      {
        y.remove_leading_chars(9);
        int n;
        if (! y.numeric_value(n) || n < 0)
        {
          cerr << "The maximum number of connections to matched atoms '-Y max_ncon=' must be a whole +ve number\n";
          display_dash_y_options(cerr);
        }

        mqs.set_max_ncon(n);

        if (verbose)
          cerr << "Matches must have at least " << n << " connections to unmatched atoms\n";
      }
      else if ("exph" == y)
      {
        add_explicit_hydrogens = 1;
        if (verbose)
          cerr << "Explicit Hydrogens will be added to the molecules\n";

        mqs.set_convert_explicit_hydrogens_to_match_any_atom(1);
      }
      else if ("ablk" == y)
      {
        set_aromatic_bonds_lose_kekule_identity(1);
        if (verbose)
          cerr << "Aromatic bonds will lose their Kekule identity\n";
      }
      else if (y.starts_with("minfm="))
      {
        y.remove_leading_chars(6);
        float f;
        if (! y.numeric_value(f) || f < 0.0 || f > 1.0)
        {
          cerr << "The min fraction atoms matched directive (minfm=) must be a valid fraction\n";
          return 2;
        }

        mqs.set_min_fraction_atoms_matched(f);
        if (verbose)
          cerr <<  "Matches will require a min fraction atom matched of " << f << '\n';
      }
      else if (y.starts_with("maxfm="))
      {
        y.remove_leading_chars(6);
        float f;
        if (! y.numeric_value(f) || f < 0.0 || f > 1.0)
        {
          cerr << "The max fraction atoms matched directive (maxfm=) must be a valid fraction\n";
          return 2;
        }

        mqs.set_max_fraction_atoms_matched(f);
        if (verbose)
          cerr <<  "Matches will require a max fraction atom matched of " << f << '\n';
      }
      else if (y.starts_with("A2A="))
      {
        y.remove_leading_chars(4);
        int a;
        if (! y.numeric_value(a) || a < 1 || a > 3)
        {
          cerr << "The A2A= qualifier must be an int between 1 and 3\n";
          return 0;
        }

        mqs.set_convert_all_aromatic_atoms_to_generic_aromatic(a);
        if (verbose)
          cerr << "Convert aromatic atoms to generic aromatic directive " << a << '\n';
      }
      else if ("rmiso" == y)
      {
        remove_isotopes_from_input_molecules = 1;

        if (verbose)
          cerr << "Will immediately remove isotopes from molecules being read\n";
      }
      else if (y.starts_with("APPC="))
      {
        append_to_comment = y;
        append_to_comment.remove_leading_chars(5);

        if (verbose)
          cerr << "Will append '" << append_to_comment << "' to each query name\n";
      }
      else if (y == "qnamsmi") {
        append_smiles_to_query_name = 1;
        if (verbose) {
          cerr << "Will append the smiles to the query name\n";
        }
      }
      else if ("test" == y)
      {
        perform_matching_test = 1;
        if (verbose)
          cerr << "Will try a match into the originating molecule for each query formed\n";
      }
      else if ("help" == y)
      {
        display_dash_y_options (cerr);
      }
      else
      {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        display_dash_y_options (cerr);
      }
    }
  }

  if (! cl.option_present('A')) {
    set_global_aromaticity_type (Daylight);
  } else if (! process_standard_aromaticity_options(cl, verbose)) {
    usage(4);
  }

  GeometryConfig geometry_config;
  if (cl.option_present('D')) {
    if (! BuildGeometricConfig(cl, 'D', geometry_config)) {
      cerr << "Cannot determine geometric constraints specifications (-D)\n";
      return 1;
    }
    if (verbose)
      cerr << "Will write geometric constraint query protos\n";
  }

  int rc = 0;

  if (cl.option_present('M')) {
    if (cl.number_elements())
    {
      cerr << "Can specify either the -M option or files on the command line\n";
      usage(29);
    }

    if (1 != cl.option_count('M'))
    {
      cerr << "Sorry, only one -M option allowed\n";
      usage(18);
    }

    IWString smiles;

    cl.value('M', smiles);

    rc = process_smiles_from_command_line(smiles, mqs, geometry_config);
  }
  else if (0 == cl.number_elements())
    usage(1);
  else if (0 == input_type && ! all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern input type(s) of command line files\n";
    return 8;
  }

  if (! cl.option_present('M')) {
    for (int i = 0; i < cl.number_elements(); i++) {
      if (! mol2qry(cl[i], input_type, mqs, geometry_config)) {
        rc = i + 1;
        break;
      }
    }
  }

  if (verbose) {
    cerr << queries_written << " queries written\n";
  }

  return rc;
}

}  // namespace mol2qry


int
main(int argc, char ** argv)
{
  mol2qry::prog_name = argv[0];

  int rc = mol2qry::mol2qry(argc, argv);

  return rc;
}
