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

void
Usage(int rc) {
  exit(rc);
}

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
    int ReadClassLabelTranslation(const IWString& dir, const IWString& fname);

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
    IWString fname = model_proto.class_label_translation();
    if (!ReadClassLabelTranslation(dir, fname)) {
      cerr << "SvmModel:cannot read class_label_translation " << fname << '\n';
      return 0;
    }
    _is_regression = false;
  }

  return Initialise(model_proto, *bit_xref);
}

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

  return ReadSupportVectors(input);
}

int
SvmModel::ReadSupportVectors(iwstring_data_source& input) {
  IW_TDT tdt;
  for (int i = 0; tdt.next(input); ++i) {
    if (! _gfp[i].Build(tdt)) {
      cerr << "SvmModel::ReadSupportVectors:cannot read " << tdt << '\n';
      return 0;
    }
  }

  return 1;
}

int
SvmModel::ReadClassLabelTranslation(const IWString& dir, const IWString& fname) {
  IWString path_name = PathName(dir, fname);
  ClassLabelTranslation::ClassLabelTranslation mapping;
  if (! ReadProto(path_name, mapping)) {
    cerr << "SvmModel::ReadClassLabelTranslation:cannot read '" << fname << "'\n";
    return 0;
  }

  for (auto [key, value] : mapping.to_numeric()) {
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
    cerr << "SvmModel::ReadClassLabelTranslation:incomplete specification " << mapping.ShortDebugString() << '\n';
    return 0;
  }

  return 1;
}

double
SvmModel::Score(IW_General_Fingerprint& fp) const {
  double result = _threshold_b;
  for (int i = 0; i < _nfingerprints; ++i) {
    result += _gfp[i].Tanimoto(fp) ;
  }
  return result;
}

IWString
SvmModel::StringScore(IW_General_Fingerprint& fp) const {
  double score = Score(fp);
  if (_is_regression) {
    IWString result;
    result << static_cast<float>(score);
    return result;
  }
  return ClassLabelForScore(score);
}

const IWString&
SvmModel::ClassLabelForScore(double score) const {
  // Classification problem, translate to class label.
  if (score < 0.0) {
    return _class_labels[0];
  } else {
    return _class_labels[1];
  }
}

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
      output << output_separator;
      output << models[i].StringScore(fp);
    } else {
      const auto score = models[i].Score(fp);
      if (write_label) {
        output << output_separator << models[i].ClassLabelForScore(score);
      }
      if (write_score) {
        output << output_separator << static_cast<float>(score);
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
  Command_Line_v2 cl(argc, argv, "-v-M=sfile-cwrite=s");
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

  if (cl.option_present("cwrite")) {
    write_label = false;
    write_score = false;
    const IWString cwrite = cl.string_value("cwrite");
    int i = 0;
    const_IWSubstring token;
    while (cwrite.nextword(token, i, ',')) {
      if (token == "label") {
        write_label = true;
      } else if (token == "score") {
        write_score = true;
      } else {
        cerr << "Unrecognised -cwrite qualifier '" << cwrite << "\n";
        Usage(1);
      }
    }
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
