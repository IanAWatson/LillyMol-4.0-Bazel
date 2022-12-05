// Read a query and write it.
// Most useful for converting legacy query format to proto.

#include <fstream>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/substructure.pb.h"

namespace echoqry {

using std::cerr;

int verbose = 0;

int write_binary_proto = 0;

void
Usage(int rc) {
  exit(rc);
}

int
WriteBinaryProto(Substructure_Query& query,
                 IWString& fname) {
  SubstructureSearch::SubstructureQuery proto = query.BuildProto();
  std::fstream output(fname.null_terminated_chars(), std::ios::out | std::ios::trunc | std::ios::binary);
  if (! proto.SerializeToOstream(&output)) {
    cerr << "WriteBinaryProto:writing to '" << fname << "' failed\n";
    return 0;
  }

  return 1;
}

int
EchoQry(Substructure_Query& query,
        int outer_ndx,
        int inner_ndx,
        const IWString& output_stem) {
  IWString fname;
  if (write_binary_proto) {
    fname << output_stem << '.' << outer_ndx << '.' << inner_ndx << "_qry.dat";
    return WriteBinaryProto(query, fname);
  } else {
    fname << output_stem << '.' << outer_ndx << '.' << inner_ndx << "_qry.txtproto";
    return query.WriteProto(fname);
  }
}

int
EchoQry(const resizable_array_p<Substructure_Query>& queries,
        int ndx,
        const IWString& output_stem) {
  for (int i = 0; i < queries.number_elements(); ++i) {
    extending_resizable_array<int> uid;
    queries[i]->assign_unique_id_from_atom_number_if_set(uid);
    if (! EchoQry(*queries[i], ndx, i, output_stem)) {
      cerr << "Fatal error processing '" << queries[i]->comment() << "'\n";
      return 0;
    }
  }

  return 1;
}

int
EchoQry(const char * token,
        int ndx,
        const IWString& output_stem) {
  resizable_array_p<Substructure_Query> queries;
  if (! process_cmdline_token('*', token, queries, verbose)) {
    cerr << "Cannot instantiate query '" << token << "'\n";
    return 0;
  }

  return EchoQry(queries, ndx, output_stem);
}

int
EchoQry(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:S:b");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('b')) {
    write_binary_proto = 1;
    if (verbose) {
      cerr << "Will write binary proto files\n";
    }
  }

  if (cl.empty()) {
    cerr << "Must specify name(s) of query file(s)\n";
    Usage(1);
  }

  IWString output_stem = "echoqry";
  if (cl.option_present('S')) {
    cl.value('S', output_stem);
  }

  for (int i = 0; i < cl.number_elements(); ++i) {
    if (! EchoQry(cl[i], i, output_stem)) {
      cerr << "Fatal error processing " << cl[i] << '\n';
      return i + 1;
    }
  }

  return 0;
}
}  // namespace echoqry

int
main(int argc, char ** argv) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  return echoqry::EchoQry(argc, argv);
}

