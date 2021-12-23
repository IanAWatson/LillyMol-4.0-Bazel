// Tests for TFDataRecord

#include "Foundational/cmdline/cmdline.h"

#include "tfdatarecord.h"

namespace test_tfdata {

const char * prog_name = nullptr;

using std::cerr;

void
Usage(int rc) {
  exit(rc);
}

int
TestTfDataRecord(const char * fname,
                 IWString_and_File_Descriptor& output) {
  iw_tf_data_record::TFDataReader reader(fname);
  const_IWSubstring data;
  int records_read = 0;

  for ( ; ! reader.eof() ; ++records_read) {
    std::optional<const_IWSubstring> maybe_data = reader.Next();
    if (! maybe_data) {
      break;
    }
  }
  cerr << "Read " << records_read << " records from " << fname << '\n';

  return 1;
}

int
TestTfDataRecord(int argc, char** argv) {
  Command_Line cl(argc, argv, "v");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  for (const char * fname : cl) {
    if (! TestTfDataRecord(fname, output)) {
      cerr << "Fatal error processing " << fname << '\n';
      return 1;
    }
  }
  return 0;
}
}  // namespace test_tfdata

int
main(int argc, char ** argv) {
  test_tfdata::prog_name = argv[0];

  int rc = test_tfdata::TestTfDataRecord(argc, argv);

  return rc;
}
