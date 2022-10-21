// Convert a .gfp file to proto form.

#include <stdlib.h>

#include <iostream>

#include "snappy.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Utilities/GFP_Tools/gfp.pb.h"
#include "gfp.h"

namespace gfp_to_proto {

using std::cerr;

const IWString smiles_tag("$SMI");
const IWString identifier_tag("PCN");
const IWString mpr("MPR");

void
Usage(int rc) {
  cerr << "Tool to convert a tdt gfp file to a serialized proto form\n";
  cerr << " -S <fname>           name of the output file\n";
  cerr << " -s                   compress with snappy\n";
  cerr << " -b                   write binary fingerprints\n";
  cerr << " -v                   verbose output\n";
  ::exit(rc);
}

class Options {
  private:
    int _verbose;

    // Compress with snappy.
    int _compress;

    int _binary;

  public:
    Options();

    int Initialise(Command_Line& cl);

    int compress() const {
      return _compress;
    }
    int binary() const {
      return _binary;
    }
};

Options::Options() {
  _verbose = 0;
  _compress = 0;
  _binary = 0;
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('s')) {
    _compress = 1;
    if (_verbose) {
      cerr << "Will compress with snappy\n";
    }
  }

  if (cl.option_present('b')) {
    _binary = 1;
    if (_verbose) {
      cerr << "Will write binary fingerprints\n";
    }
  }


  return 1;
}

int
AsBinary(IW_General_Fingerprint& fp,
         const IWString& smiles,
         Options& options,
         iw_tf_data_record::TFDataWriter& output) {
  gfp::Fingerprint proto;
  const IWString& id = fp.id();
  proto.set_id(id.data(), id.length());

  if (! smiles.empty()) {
    proto.set_smiles(smiles.data(), smiles.length());
  }

  const Molecular_Properties_Integer& mprop_integer = fp.molecular_properties_integer();
  if (mprop_integer.active()) {
    std::string mpr("MPR");
    const int* data = mprop_integer.rawdata();
    std::string value;
    value.assign(reinterpret_cast<const char*>(data), mprop_integer.nproperties() * sizeof(int));
    (*proto.mutable_fp())[mpr] = value;
  }

  return output.WriteSerializedProto<gfp::Fingerprint>(proto);
}

int
AsBinary(IW_TDT& tdt,
            Options& options,
            iw_tf_data_record::TFDataWriter& output) {
  int fatal;
  IW_General_Fingerprint fp;
  if (! fp.construct_from_tdt(tdt, fatal)) {
    cerr << "Cannot build fingerprint from '" << tdt << "'\n";
    return 0;
  }

  IWString smiles;
  tdt.dataitem_value(smiles_tag, smiles);

  return AsBinary(fp, smiles, options, output);
}

int
GfpToProto(IW_TDT& tdt,
            Options& options,
            iw_tf_data_record::TFDataWriter& output) {
  if (options.binary()) {
    return AsBinary(tdt, options, output);
  }

  gfp::Fingerprint proto;

  int ptr = 0;
  const_IWSubstring tag, value;
  while (tdt.next_dataitem_value(tag, value, ptr)) {
    if (tag == smiles_tag) {
      proto.set_smiles(value.data(), value.length());
    } else if (tag == identifier_tag) {
      proto.set_id(value.data(), value.length());
    } else if (tag == mpr) {
      (*proto.mutable_fp())[tag.AsString()] = value.AsString();
    } else if (options.compress()) {
      std::string compressed;
      snappy::Compress(value.data(), value.length(), &compressed);
      (*proto.mutable_fp())[tag.AsString()] = compressed;
    } else {
      (*proto.mutable_fp())[tag.AsString()] = value.AsString();
    }
  }

  return output.WriteSerializedProto<gfp::Fingerprint>(proto);
}

int
GfpToProto(iwstring_data_source& input,
            Options& options,
            iw_tf_data_record::TFDataWriter& output) {
  IW_TDT tdt;
  while (tdt.next(input)) {
    if (! GfpToProto(tdt, options, output)) {
      return 0;
    }
  }

  return 1;
}

int
GfpToProto(const char* fname,
            Options& options,
            iw_tf_data_record::TFDataWriter& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "GfpToProto::cannot open '" << fname << "'\n";
    return 0;
  }

  return GfpToProto(input, options, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vS:sb");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  Options options;
  if (! options.Initialise(cl)) {
    Usage(1);
  }

  if (! cl.option_present('S')) {
    cerr << "Must specify output file name via the -S option\n";
    Usage(1);
  }

  iw_tf_data_record::TFDataWriter output;
  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    if (! output.Open(fname)) {
      cerr << "GfpToProto:cannot open output file '" << fname << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Output written tto '" << fname << "'\n";
    }
  }

  if (cl.empty()) {
    cerr << "INsufficient arguments\n";
    Usage(1);
  }


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
