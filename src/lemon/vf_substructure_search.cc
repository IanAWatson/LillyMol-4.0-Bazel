// Do VF2 based substructure searches

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Lib/llymol_lemon.h"

namespace vf_substructure_search {

using std::cerr;

using llymol_lemon::LemonQuery;
void
Usage(int rc) {

  cerr << " -s <smarts>         queries as smarts\n";
  cerr << " -q ...              query specification\n";
  cerr << " -M <fname>          molecules to be used as queries\n";
  cerr << " -m <fname>          write     matching molecules to <fname>\n";
  cerr << " -n <fname>          write non matching molecules to <fname>\n";
  cerr << " -o <type>           specify output file type(s), default SMI\n";
  cerr << " -c                  remove chirality\n";
  cerr << " -l                  strip to largest fragment\n";
  cerr << " -g ...              chemical standardisation\n";
  cerr << " -i <type>           input type\n";
  cerr << " -v                  verbose output\n";

  ::exit(rc);
}

class VfSubstructureSearch {
  private:
    int _verbose;

    int _molecules_read;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    Chemical_Standardisation _chemical_standardisation;

    resizable_array_p<Molecule> _query_mol;
    resizable_array_p<lemon::ListGraph> _query_graph;
    resizable_array_p<lemon::ListGraph::NodeMap<int>> _query_nodemap;

    resizable_array_p<Substructure_Query> _query;
    resizable_array_p<LemonQuery> _lemon_query;

    int _molecules_matching;

    extending_resizable_array<int> _queries_matched;
    extending_resizable_array<int> _number_matches;

    Molecule_Output_Object _stream_matches;
    Molecule_Output_Object _stream_non_matches;

  // private functions
    int OpenOutputStream(const Command_Line& cl, const_IWSubstring& fname,
                                        Molecule_Output_Object& destination);

    int Process(Molecule& m, lemon::ListGraph& g2);

  public:
    VfSubstructureSearch();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int Process(Molecule& m);

    int Report(std::ostream& output) const;
};

VfSubstructureSearch::VfSubstructureSearch() {
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;

  _molecules_matching = 0;
}

void
MakeIntNodeMap(Molecule& m,
               const lemon::ListGraph& g,
               lemon::ListGraph::NodeMap<int>& destination) {
  int ndx = 0;
  for (lemon::ListGraph::NodeIt nit(g); nit != lemon::INVALID; ++nit) {
    destination[nit] = 2 * m.atomic_number(ndx);
    if (m.is_aromatic(ndx)) {
      destination[nit] += 1;
    } 
    cerr << "Atom " << ndx << ' ' << m.smarts_equivalent_for_atom(ndx) << " value " << destination[nit] << '\n';
    ++ndx;
  }
}

int
Label(Molecule& m,
      atom_number_t zatom) {
  int result = 2 * m.atomic_number(zatom);
  if (m.is_aromatic(zatom)) {
    ++result;
  }

  return result;
}

void
AssignLabels(Molecule& m1,
             Molecule& m2,
             const lemon::ListGraph& g1,
             const lemon::ListGraph& g2,
             lemon::ListGraph::NodeMap<int>& map1,
             lemon::ListGraph::NodeMap<int>& map2) {

  int next_to_assign = 0;
  extending_resizable_array<int> assigned(-1);

  for ( lemon::ListGraph::NodeIt nit(g1); nit != lemon::INVALID; ++nit) {
    int id = g1.id(nit);
    const int label = Label(m1, id);
    if (assigned[label] == -1) {
      assigned[label] = next_to_assign;
      map1[nit] = next_to_assign;
      ++next_to_assign;
    } else {
      map1[nit] = assigned[label];
    }
  }

  for (lemon::ListGraph::NodeIt nit(g2); nit != lemon::INVALID; ++nit) {
    const int id = g2.id(nit);
    const int label = Label(m2, id);
    if (assigned[label] < 0) {
      assigned[label] = next_to_assign;
      map2[nit] = next_to_assign;
      ++next_to_assign;
    } else {
      map2[nit] = assigned[label];
    }
  }
}

void
AssignLabelsOld(Molecule& m1,
             Molecule& m2,
             const lemon::ListGraph& g1,
             const lemon::ListGraph& g2,
             lemon::ListGraph::NodeMap<int>& map1,
             lemon::ListGraph::NodeMap<int>& map2) {
  std::unordered_map<int, int> atype_to_index;

  for ( lemon::ListGraph::NodeIt nit(g1); nit != lemon::INVALID; ++nit) {
    int id = g1.id(nit);
    const int label = Label(m1, id);
    auto iter = atype_to_index.find(label);
    if (iter == atype_to_index.end()) {
      auto s = atype_to_index.size();
      atype_to_index[label] = s;
      map1[nit] = s;
      // cerr << "Type " << label << " assigned " << s << '\n';
    } else {
      map1[nit] = iter->second;
    }
  }

  for (lemon::ListGraph::NodeIt nit(g2); nit != lemon::INVALID; ++nit) {
    const int id = g2.id(nit);
    const int label = Label(m2, id);
    auto iter = atype_to_index.find(label);
    if (iter == atype_to_index.end()) {
      auto s = atype_to_index.size();
      atype_to_index[label] = s;
      map2[nit] = s;
      // cerr << "Type " << label << " assigned " << s << '\n';
    } else {
      map2[nit] = iter->second;
    }
  }
}

int
VfSubstructureSearch::Initialise(Command_Line& cl) {
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

#ifdef PROCESS_QUERIES_TOTO
  if (cl.option_present('s')) {
    const_IWSubstring smarts;
    for (int i = 0; cl.value('s', smarts, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(smarts)) {
        cerr << "VfSubstructureSearch::initialise:invalid smarts '" << smarts << "'\n";
        return 0;
      }
      _query << q.release();
    }
  }

  if (cl.option_present('q')) {
    if (! process_queries(cl, _query, _verbose, 'q')) {
      cerr << "VfSubstructureSearch::initialise:cannot process queries (-q)\n";
      return 0;
    }
  }

  for (Substructure_Query* q : _query) {
    std::unique_ptr<llymol_lemon::LemonQuery> lq = std::make_unique<llymol_lemon::LemonQuery>();
    if (! lq->Build(*q)) {
      cerr << "VfSubstructureSearch::Initialise:cannot convert query to lemon form\n";
      q->debug_print(cerr);
      return 0;
    }
    _lemon_query << lq.release();
  }

  if (_query.empty()) {
    cerr << "VfSubstructureSearch::Initialise:no queries\n";
    return 0;
  }
#endif

  if (cl.option_present('M')) {
    const_IWSubstring smiles;
    for (int i = 0; cl.value('M', smiles, i); ++i) {
      std::unique_ptr<Molecule> m = std::make_unique<Molecule>();
      if (! m->build_from_smiles(smiles)) {
        cerr << "Cannot interpret smiles '" << smiles << "'\n";
        return 0;
      }
      if (m->empty()) {
        cerr << "Ignoring empty molecule\n";
        continue;
      }
      _query_mol << m.release();
    }
  }

  for (Molecule* m : _query_mol) {
    std::unique_ptr<lemon::ListGraph> g = std::make_unique<lemon::ListGraph>();
    if (! llymol_lemon::BuildLemonGraph(*m, *g)) {
      cerr << "Cannot build lemon graph from '" << m->smiles() << "'\n";
      return 0;
    }

    _query_graph << g.release();
  }

  int nq = _query_graph.size();
  for (int i = 0; i < nq; ++i) {
    std::unique_ptr<lemon::ListGraph::NodeMap<int>> mp = std::make_unique<lemon::ListGraph::NodeMap<int>>(*_query_graph[i]);
    //MakeIntNodeMap(*_query_mol[i], *_query_graph[i], *mp);
    _query_nodemap << mp.release();
  }

  if (_query_graph.empty()) {
    cerr << "No queries defined\n";
  }

  if (_verbose) {
    cerr << "Defined " << _query_graph.size() << " lemon queries\n";
  }

  if (cl.option_present('m')) {
    const_IWSubstring fname = cl.string_value('m');
    if (! OpenOutputStream(cl, fname, _stream_matches)) {
      cerr << "VfSubstructureSearch::Initialise:cannot open stream for matches '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Molecules matching written to '" << fname << "'\n";
    }
  }

  if (cl.option_present('n')) {
    const_IWSubstring fname = cl.string_value('n');
    if (! OpenOutputStream(cl, fname, _stream_non_matches)) {
      cerr << "VfSubstructureSearch::Initialise:cannot open stream for non matches '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Molecules not matching written to '" << fname << "'\n";
    }
  }

  return 1;
}

int
VfSubstructureSearch::OpenOutputStream(const Command_Line& cl,
                                        const_IWSubstring& fname,
                                        Molecule_Output_Object& destination) {
  if (cl.option_present('o')) {
    if (! destination.determine_output_types(cl, 'o')) {
      cerr << "VfSubstructureSearch::OpenOutputStream:cannot determine output type(s)\n";
      return 0;
    }
  } else {
    destination.add_output_type(FILE_TYPE_SMI);
  }

  if (destination.would_overwrite_input_files(cl, fname)) {
    cerr << "VfSubstructureSearch::OpenOutputStream:cannot overwrite input file(s) '" << fname << "'\n";
    return 0;
  }

  if (! destination.new_stem(fname)) {
    cerr << "VfSubstructureSearch::OpenOutputStream:cannot open '" << fname << "'\n";
    return 0;
  }

  return 1;
}

int
VfSubstructureSearch::Preprocess(Molecule& m) {
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

  return 1;
}

int
VfSubstructureSearch::Report(std::ostream& output) const {
  output << "VfSubstructureSearch:read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }

  const float f = iwmisc::Fraction<float>(_molecules_matching, _molecules_read);
  output  << _molecules_matching << " molecules matched " << f << '\n';

  for (int i = 0; i < _queries_matched.number_elements(); ++i) {
    if (_queries_matched[i]) {
      output << _queries_matched[i] << " molecules had " << i << " hits\n";
    }
  }

  for (int i = 0; i < _number_matches.number_elements(); ++i) {
    if (_number_matches[i]) {
      output << _number_matches[i] << " molecules had " << i << " embeddings\n";
    }
  }

  return 1;
}

int
VfSubstructureSearch::Process(Molecule& m) {
  ++_molecules_read;

  lemon::ListGraph graph;
  if (! llymol_lemon::BuildLemonGraph(m, graph)) {
    cerr << "Cannot build graph " << m.smiles() << '\n';
    return 0;
  }

  return Process(m, graph);
}

int
VfSubstructureSearch::Process(Molecule& m,
                              lemon::ListGraph& g2) {
  // cerr << "Checking " << _query_graph.size() << " lemon queries vs " << m.smiles() << '\n';
  int queries_matching = 0;
  int nq = _query_graph.size();

  for (int i = 0; i < nq; ++i) {
    lemon::ListGraph& g1 = *_query_graph[i];

    lemon::ListGraph::NodeMap<int> label1(g1);
    lemon::ListGraph::NodeMap<int> label2(g2);
    AssignLabels(*_query_mol[i], m, *_query_graph[i], g2, label1, label2);

    lemon::ListGraph::NodeMap<lemon::ListGraph::Node> lmap(g1);
    lemon::Vf2pp vf2pp(g1, g2, lmap, label1, label2);
    vf2pp.mappingType(lemon::SUBGRAPH);
    int matches = 0;
    while (vf2pp.find()) {
      ++matches;
      // cerr << "Matches " << matches << '\n';
      for (lemon::ListGraph::NodeIt nit(g1); nit != lemon::INVALID; ++nit) {
        //int i = g1[lmap[nit]];
        //cerr << i << '\n';
      }
    }

    if (matches == 0) {
      //cerr << "No matches to query\n";
    } else {
      ++_number_matches[matches];
      ++queries_matching;
    }
  }

  // cerr << queries_matching <<  " of " << _lemon_query.size() << " queries matched\n";

  ++_queries_matched[queries_matching];
  if (queries_matching) {
    ++_molecules_matching;
    if (_stream_matches.active()) {
      _stream_matches.write(m);
    }
  } else if (_stream_non_matches.active()) {
    _stream_non_matches.write(m);
  }

  return 1;
}

int
ReplaceCore(VfSubstructureSearch& vf_sss,
            Molecule& m,
            IWString_and_File_Descriptor& output) {
  return vf_sss.Process(m);
}

int
ReplaceCore(VfSubstructureSearch& vf_sss,
            data_source_and_type<Molecule>& input,
            IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! vf_sss.Preprocess(*m)) {
      return 0;
    }

    if (! ReplaceCore(vf_sss, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
ReplaceCore(VfSubstructureSearch& vf_sss,
            const char * fname,
            FileType input_type,
            IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ReplaceCore(vf_sss, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:g:lcs:q:m:n:o:M:");
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

  VfSubstructureSearch vf_sss;
  if (! vf_sss.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (all_files_recognised_by_suffix(cl)) {
  } else if (1 == cl.size() && 0 == strcmp("-", cl[0])) {
    input_type = FILE_TYPE_SMI;
  } else {
    cerr << "Cannot discern all file types, use the -i option\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);
  for (const char* fname : cl) {
    if (! ReplaceCore(vf_sss, fname, input_type, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    vf_sss.Report(cerr);
  }

  return 0;
}

}  // namespace vf_substructure_search

int
main(int argc, char** argv) {
  int rc = vf_substructure_search::Main(argc, argv);

  return rc;
}
