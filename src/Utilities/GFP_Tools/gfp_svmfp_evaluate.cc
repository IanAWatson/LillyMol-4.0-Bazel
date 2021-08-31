// Score an svmfp model.
// A model is described by a set of weighted fingerprints
// and a bias term.

#include <fstream>

#include "google/protobuf/text_format.h"
#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"

#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Utilities/GFP_Tools/gfp.h"
#include "Utilities/GFP_Tools/gfp_to_svm_lite.pb.h"
#include "Utilities/GFP_Tools/gfp_to_svm_lite.pb.h"
#include "Utilities/GFP_Tools/svmfp_model.pb.h"
#include "Utilities/General/class_label_translation.pb.h"

namespace gfp_svmfp_evaluate {
using std::cerr;

int verbose = 0;

IWString weight_tag = "WEIGHT<";

char output_separator = ' ';

// If we have a classification problem, what do we write?
// By default, we write just the numeric score, but we can
// write other things.
bool write_label = true;
bool write_score = false;

Fraction_as_String fraction_as_string;

void
Usage(int rc) {
  cerr << "Evaluates svmfp model(s)\n";
  cerr << " -mdir <dir>     one or more model directories\n";
  cerr << " -cwrite ...     what to write for classification models, '-cwrite help' for info\n";
  cerr << " -sas <n>        write the score as string with <n> digits of accuracy\n";
  cerr << " -v              verbose output\n";

  exit(rc);
}

void
DisplayCWriteOptions(std::ostream& output) {
  output << "The -cwrite option controls output from classification models\n";
  output << " -cwrite label      write the label\n";
  output << " -cwrite score      write the score\n";

  exit(0);
}

// Read a textproto from `fname` into `proto`.
// The iwstring_data_source is just used to get a file descriptor.
template <typename Proto>
int
ReadProto(IWString& fname,
          Proto& proto) {
  iwstring_data_source input(fname.null_terminated_chars());
  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  google::protobuf::io::FileInputStream file_input_stream(input.fd());

  if (! google::protobuf::TextFormat::Parse(&file_input_stream, &proto)) {
    cerr << "GfpToFeature::FromProto:cannot parse '" << fname << "'\n";
    return 0;
  }

  return 1;
}

// Read a binary encoded proto from `fname` into `proto`.
template <typename Proto>
int
ReadBinaryProto(IWString& fname,
                Proto& proto) {
  using std::ios;
  std::fstream input(fname.null_terminated_chars(), ios::in | ios::binary);
  if (! input) {
    cerr << "ReadBinaryProto:cannot open '" << fname << "'\n";
    return 0;
  }
  if (! proto.ParseFromIstream(&input)) {
    cerr << "ReadBinaryProto::cannot decode '" << fname << "'\n";
    return 0;
  }

  return 1;
}

IWString
PathName(const IWString& dirname, const IWString& fname) {
  IWString result;
  result << dirname << '/' << fname;
  return result;
}

// Read a textproto of type `Proto` from `dir/fname`.
template <typename Proto>
std::optional<Proto>
ReadProto(const IWString& dirname, const std::string& fname) {
  IWString path_name = PathName(dirname, fname);
  Proto result;
  if (! ReadProto(path_name, result)) {
    return std::nullopt;
  }
  return result;
}

// Read a binary encoded proto of type `Proto` from `dirname/fname`.
template <typename Proto>
std::optional<Proto>
ReadBinaryProto(const IWString& dirname, const std::string& fname) {
  IWString path_name = PathName(dirname, fname);
  Proto result;
  if (! ReadBinaryProto(path_name, result)) {
    return std::nullopt;
  }
  return result;
}

// Fingerprint and associated weight.
class WeightedFingerprint : public IW_General_Fingerprint {
  private:
    double _weight;

  public:
    WeightedFingerprint();

    double weight() const { return _weight;}

    int Build(IW_TDT& tdt);

    double Tanimoto(IW_General_Fingerprint& fp);
};

WeightedFingerprint::WeightedFingerprint() {
  _weight = 0.0;
}

int
WeightedFingerprint::Build(IW_TDT& tdt) {
  int fatal = 0;
  if (! IW_General_Fingerprint::construct_from_tdt(tdt, fatal)) {
    cerr << "WeightedFingerprint:Build:cannot built gfp\n";
    return 0;
  }

  if (! tdt.dataitem_value(weight_tag, _weight) || _weight < -1.0 || _weight > 1.0) {
    cerr << "WeightedFingerprint::Built:invalid weight\n";
    return 0;
  }

  return 1;
}

double
WeightedFingerprint::Tanimoto(IW_General_Fingerprint& fp) {
  IW_General_Fingerprint* me = this;
  return me->tanimoto(fp) * _weight;
}

// Model class populated from a svmfp_model proto.
class SvmModel {
  private:
    int _nfingerprints;
    WeightedFingerprint* _gfp;

    double _threshold_b;

    float _weight;

    // Read from the proto.
    IWString _response_name;

    // A quick way to check if this is a classification problem.
    bool _is_regression;

    // If we have a classification problem, the class labels to apply.
    IWString _class_labels[2];

  // private functions

    int Initialise(const SvmfpModel::SvmfpModel& model_proto,
                const GfpBitToFeatureMap::GfpBitXref& bit_xref);
    int ReadSupportVectors(const IWString& dir, const IWString& fname);
    int ReadSupportVectors(iwstring_data_source& input);
    int ReadClassLabelTranslation(const IWString& dir, const std::string& fname);

  public:
    SvmModel();
    ~SvmModel();

    // Initialise an SvmModel based on a SvmfpModel proto.
    // Argument `directory` is needed so that the files
    // specified in `model_proto` can be found.
    int Initialise(const SvmfpModel::SvmfpModel& model_proto,
                   const IWString& directory);

    bool is_regression() const { return _is_regression;}

    const IWString& response_name() const { return _response_name;}

    // Return a model evaluation in various forms.
    double Score(IW_General_Fingerprint& fp) const;
    IWString StringScore(IW_General_Fingerprint& fp) const;

    // Given a numeric score, translate to the appropriate label.
    const IWString& ClassLabelForScore(double score) const;
};

SvmModel::SvmModel() {
  _nfingerprints = 0;
  _gfp = nullptr;
  _threshold_b = 0.0;
  _weight = 1.0;
  _is_regression = true;
}

SvmModel::~SvmModel() {
  if (_gfp != nullptr) {
    delete [] _gfp;
  }
  _nfingerprints = 0;
  _gfp = nullptr;
  _weight = 0.0;
}

// Read information from `model_proto` to build state.
int
SvmModel::Initialise(const SvmfpModel::SvmfpModel& model_proto,
                      const IWString& dir) {
  if (! model_proto.has_bit_xref()) {
    cerr << "SvmModel::Initialise:missing bit_xref\n";
    return 0;
  }

  if (! model_proto.has_support_vectors()) {
    cerr << "SvmModel::Initialise:missing support vectors";
    return 0;
  }

  if (! model_proto.has_threshold_b()) {
    cerr << "SvmModel::Initialise:no threshold_b\n";
    return 0;
  }

  _threshold_b = model_proto.threshold_b();

  _response_name = model_proto.response_name();

  std::string fname = model_proto.bit_xref();
  std::optional<GfpBitToFeatureMap::GfpBitXref> bit_xref = ReadProto<GfpBitToFeatureMap::GfpBitXref>(dir, fname);
  if (! bit_xref) {
    cerr << "Cannot read bit cross reference proto '" << fname << "'\n";
    return 0;
  }

  fname = model_proto.support_vectors();
  if (! ReadSupportVectors(dir, fname)) {
    cerr << "SvmModel::Initialise:cannot read support vectors '" << fname << "'\n";
  }

  if (model_proto.has_class_label_translation()) {
    const std::string fname = model_proto.class_label_translation();
    if (!ReadClassLabelTranslation(dir, fname)) {
      cerr << "SvmModel:cannot read class_label_translation " << fname << '\n';
      return 0;
    }
    _is_regression = false;
  }

  return Initialise(model_proto, *bit_xref);
}

// `fname` contains the name of the file containing the weighted support vectors.
// Read that file and populate _gfp.
int
SvmModel::ReadSupportVectors(const IWString& dir, const IWString& fname) {
  IWString path_name = PathName(dir, fname);
  iwstring_data_source input(path_name.null_terminated_chars());
  if (! input.good()) {
    cerr << "SvmModel::ReadSupportVectors:cannot open '" << path_name << "'\n";
    return 0;
  }

  std::unique_ptr<RE2> vbar = std::make_unique<RE2>("^\\|$");
  _nfingerprints = input.grep(*vbar);
  if (_nfingerprints == 0) {
    cerr << "SvmModel::ReadSupportVectors:no fingerprints in '" << path_name << "'\n";
    return 0;
  }

  _gfp = new WeightedFingerprint[_nfingerprints];

  if (verbose) {
    cerr << "Reading " << _nfingerprints << " support vectors\n";
  }

  return ReadSupportVectors(input);
}

// Populate the already allocated _gfp array with the contents of `input`.
int
SvmModel::ReadSupportVectors(iwstring_data_source& input) {
  IW_TDT tdt;
  for (int i = 0; tdt.next(input); ++i) {
    if (! _gfp[i].Build(tdt)) {
      cerr << "SvmModel::ReadSupportVectors:cannot read " << tdt << '\n';
      return 0;
    }
  }

  if (verbose) {
    cerr << "Support vectors initialized\n";
  }

  return 1;
}

// The file `dir/fname` is a ClassLabelTranslation proto. Read that file
// and use it to establish the _class_labels array.
int
SvmModel::ReadClassLabelTranslation(const IWString& dir, const std::string& fname) {
  std::optional<ClassLabelTranslation::ClassLabelTranslation> mapping =
    ReadProto<ClassLabelTranslation::ClassLabelTranslation>(dir, fname);
  if (! mapping) {
    cerr << "SvmModel::ReadClassLabelTranslation:cannot read '" << fname << "'\n";
    return 0;
  }

  for (auto [key, value] : mapping->to_numeric()) {
    if (value == -1) {
      _class_labels[0] = key;
    } else if (value == 1) {
      _class_labels[1] = key;
    } else {
      cerr << "SvmModel::ReadClassLabelTranslation:unrecognised pair '" << key << "' value " << value << '\n';
      return 0;
    }
  }

  if (_class_labels[0].empty() || _class_labels[1].empty()) {
    cerr << "SvmModel::ReadClassLabelTranslation:incomplete specification " << mapping->ShortDebugString() << '\n';
    return 0;
  }

  return 1;
}

// Score a single fingerprint. The weighted similarity to our support vectors.
double
SvmModel::Score(IW_General_Fingerprint& fp) const {
  double result = _threshold_b;
  for (int i = 0; i < _nfingerprints; ++i) {
    result += _gfp[i].Tanimoto(fp) ;
  }
  return result;
}

#ifdef NOT_NEEDED_ANY_MORE
IWString
SvmModel::StringScore(IW_General_Fingerprint& fp) const {
  double score = Score(fp);
  if (_is_regression) {
    IWString result;
    if (fraction_as_string.active()) {
      return fraction_as_string.string_for_fraction(score);
    } else {
      result << static_cast<float>(score);
    }
    return result;
  }
  return ClassLabelForScore(score);
}
#endif

const IWString&
SvmModel::ClassLabelForScore(double score) const {
  // Classification problem, translate to class label.
  if (score < 0.0) {
    return _class_labels[0];
  } else {
    return _class_labels[1];
  }
}

// Append a string representation of `score` to `output`.
// If fraction_as_string is active, use it.
void
AppendScore(float score,
            IWString_and_File_Descriptor& output) {
  if (fraction_as_string.active()) {
    fraction_as_string.append_number(output, score);
  } else {
    output << output_separator << score;
  }
}

// Seemingly nothing left to do. Maybe remove...
int
SvmModel::Initialise(const SvmfpModel::SvmfpModel& model_proto,
                const GfpBitToFeatureMap::GfpBitXref& bit_xref) {
  return 1;
}

int
GfpSvmfpEvaluate(IW_General_Fingerprint& fp,
                 SvmModel* models,
                 int nmodels,
                 IWString_and_File_Descriptor& output) { 
  append_first_token_of_name(fp.id(), output);

  for (int i = 0; i < nmodels; ++i) {
    if (models[i].is_regression()) {
      AppendScore(models[i].Score(fp), output);
    } else {
      const auto score = models[i].Score(fp);
      if (write_label) {
        output << output_separator << models[i].ClassLabelForScore(score);
      }
      if (write_score) {
        AppendScore(score, output);
      }
    }
  }

  output << '\n';
  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

int
GfpSvmfpEvaluate(IW_TDT& tdt,
                 SvmModel* models,
                 int nmodels,
                 IWString_and_File_Descriptor& output) { 
  IW_General_Fingerprint fp;
  int fatal;
  if (! fp.construct_from_tdt(tdt, fatal)) {
    cerr << "Cannot decode TDT\n";
    return 0;
  }

  return GfpSvmfpEvaluate(fp, models, nmodels, output);
}

int
GfpSvmfpEvaluate(iwstring_data_source& input,
                 SvmModel* models,
                 int nmodels,
                 IWString_and_File_Descriptor& output) { 
  IW_TDT tdt;
  while (tdt.next(input)) {
    if (! GfpSvmfpEvaluate(tdt, models, nmodels, output)) {
      cerr << "Error processing " << tdt << '\n';
      return 0;
    }
  }

  return 1;
}

IWString
IwDirname(const IWString& fname) {
  IWString result(fname);
  int ndx = fname.rindex('/');
  if (ndx < 0) {
    return ".";
  }
  result.iwtruncate(ndx);
  return result;
}

int
GfpSvmfpEvaluate(const char * fname,
                 SvmModel* models,
                 int nmodels,
                 IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return GfpSvmfpEvaluate(input, models, nmodels, output);
}

int
GfpSvmfpEvaluate(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv, "-v-M=sfile-cwrite=s-sas=ipos");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('M')) {
    cerr << "Must specify model proto file (-M)\n";
    Usage(1);
  }

  const int nmodels = cl.option_count("M");
  std::unique_ptr<SvmModel[]> models(new SvmModel[nmodels]);
  for (int i = 0; i < nmodels; ++i) {
    IWString fname = cl.string_value('M', i);
    const IWString dir_name = IwDirname(fname);
    SvmfpModel::SvmfpModel model;
    cerr << "fname " << fname << " dir_name '" << dir_name << "'\n";
    if (! ReadBinaryProto(fname, model)) {
      cerr << "Cannot read model proto file '" << fname << "'\n"; 
      return 1;
    }
    if (! models[i].Initialise(model, dir_name)) {
      cerr << "Cannot initialise model " << model.ShortDebugString() << '\n';
      return 1;
    }
  }

  if (verbose && nmodels > 1) {
    cerr << "Read " << nmodels << " model protos\n";
  }

  // Either csv or separate options.
  if (cl.option_present("cwrite")) {
    write_label = false;
    write_score = false;
    IWString cwrite;
    for (int i = 0; cl.value("cwrite", cwrite, i); ++i) {
      int j = 0;
      const_IWSubstring token;
      while (cwrite.nextword(token, j, ',')) {
        if (token == "label") {
          write_label = true;
        } else if (token == "score") {
          write_score = true;
        } else if (token == "help") {
          DisplayCWriteOptions(cerr);
        } else {
          cerr << "Unrecognised -cwrite qualifier '" << cwrite << "\n";
          Usage(1);
        }
      }
    }
  }

  if (cl.option_present("sas")) {
    int npoints;
    cl.value("sas", npoints);
    fraction_as_string.set_leading_string(IWString(output_separator));
    fraction_as_string.initialise(-1.2, 1.2, npoints);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  output << "ID";
  for (int i = 0; i < nmodels; ++i) {
    if (models[i].is_regression()) {
      output << output_separator << models[i].response_name();
      continue;
    }
    output << output_separator << models[i].response_name() << "_label";
    if (write_label && write_score) {
      output << output_separator << models[i].response_name() << "_score";
    }
  }
  output << '\n';
  for (const char * fname : cl) {
    if (! GfpSvmfpEvaluate(fname, models.get(), nmodels, output)) {
      cerr << "GfpSvmfpEvaluate:cannot process '" << fname << "'\n";
      return 1;
    }
  }

  return 0;
}

}  // namespace gfp_svmfp_evaluate

int
main(int argc, char ** argv) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  return gfp_svmfp_evaluate::GfpSvmfpEvaluate(argc, argv);
}
