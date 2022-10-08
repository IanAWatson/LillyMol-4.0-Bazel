// Convert a .gfp file to proto form.

#include <stdlib.h>

#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "gfp.h"
#include "Utilities/GFP_Tools/gfp.pb.h"

namespace gfp_to_proto {

using std::cerr;

void
Usage(int rc) {
  cerr << "Tool to ...\n";
  cerr << " -v                   verbose output\n";
  ::exit(rc);
}

class Options {
  private:
    int _verbose;

  public:
    Options();

    int Initialise(Command_Line& cl);
};

Options::Options() {
  _verbose = 0;
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  return 1;
}

int
GfpToProto(iwstring_data_source& input,
            Options& options,
            IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
  }

  return 1;
}

int
GfpToProto(const char* fname,
            Options& options,
            IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "GfpToProto::cannot open '" << fname << "'\n";
    return 0;
  }

  return GfpToProto(input, options, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vF:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  Options options;
  if (! options.Initialise(cl)) {
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "INsufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! GfpToProto(fname, options, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  return 0;
}

}  // namespace gfp_to_proto

int
main (int argc, char ** argv)
{
  int rc = gfp_to_proto::Main(argc, argv);

  return rc;
}
