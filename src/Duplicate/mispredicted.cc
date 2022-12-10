// Across a set of predicted values identify those molecules that
// are consistently predicted with high errors.

#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <optional>
#include <random>

#define IWQSORT_FO_IMPLEMENTATION 1

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

namespace mispredicted {

using std::cerr;

void
Usage(int rc) {
  cerr << "Examines multiple predicted values and identifies molecules with the worst prediced values\n";
  cerr << "Output is ordered by degree of mismatch\n";
  cerr << " -A <fname>    experimental data file name - correct values\n";
  cerr << " -p <col>      column in -A file containing data values\n";
  cerr << " -C            problem is classificatoin type\n";
  cerr << " -S <fname>    optional smiles file. If present output will be a csv smiles file\n";
  cerr << " -v            verbose output\n";

  ::exit(rc);
}

class Predictions {
  private:
    int _verbose;

    // The number of predicted files we read.
    int _files_read;

    bool _classification;

    // The name of the response being examined.
    IWString _response;

    // The columns in which the identifier and the data are found.
    // Note that we leave the door open for the columns being different
    // in the -A file and in the predicted files.
    int _identifier_column;
    int _expt_column_data;
    int _pred_column;

    // Input token separator.
    char _input_token_sep;

    // Data gets populated depending on whether we are activity or
    // classification data. Thought of making this a template, but
    // decided against it.
    IW_STL_Hash_Map_float _id_to_activity;
    IW_STL_Hash_Map_String _id_to_class;

    // Statistics about the experimental data.
    IW_STL_Hash_Map_int _class_count;
    Accumulator<double> _acc_activity;

    // The predicted values.
    // for continuous data, an array of predicted values.
    IW_STL_Hash_Map<IWString, Accumulator<double>> _id_to_pred_continuous;
    // For categorical data, a mapping from class to count.
    IW_STL_Hash_Map<IWString, IW_STL_Hash_Map_int> _id_to_pred_class;

    // For convenience, if we are given a smiles file, we write a csv smiles.
    IW_STL_Hash_Map_String _id_to_smiles;

  // private functions.

    int ReadExperimentalData(IWString& fname);
    int ReadExperimentalData(iwstring_data_source& input);
    int ReadExperimentalDataRecord(const const_IWSubstring& buffer);

    int DetermineResponse(const const_IWSubstring& header);

    int Accumulate(iwstring_data_source& input);
    int AccumulateLine(const const_IWSubstring& buffer);

    int ReadSmiles(IWString& fname);
    int ReadSmiles(iwstring_data_source& input);
    int ReadSmilesLine(const const_IWSubstring& buffer);

    int AddClassification(const IWString& id, const IWString& expt);
    int AddContinuous(const IWString& id, const IWString& expt);

    int ProcessClassification(IWString_and_File_Descriptor& output);
    int ProcessContinuous(IWString_and_File_Descriptor& output);
  
    int MaybeWriteSmiles(IWString& id,
                              char sep,
                              IWString_and_File_Descriptor& output) const;

  public:
    Predictions();

    int Initialise(Command_Line& cl);

    int Accumulate(const char* fname);

    int Process(IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

Predictions::Predictions() {
  _verbose = 0;
  _classification = false;

  _files_read = 0;

  _input_token_sep = ' ';

  _identifier_column = 0;
  _expt_column_data = 1;
  _pred_column = 1;
}

int
Predictions::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (! cl.option_present('A')) {
    cerr << "Must specify experimental data via the -A option\n";
    return 0;
  }

  if (cl.option_present('C')) {
    _classification = true;
    if (_verbose) {
      cerr << "Will treat as classification data\n";
    }
  }

  if (cl.option_present('p')) {
    if (! cl.value('p', _pred_column) || _pred_column < 1) {
      cerr << "The predicted data column (-p) must be a valid column number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will fetch predicted values from column " << _pred_column << '\n';
    }
    --_pred_column;
  }

  if (cl.option_present('A')) {
    IWString fname = cl.string_value('A');
    if (! ReadExperimentalData(fname)) {
      cerr << "Cannot read experimental data '" << fname << "'\n";
      return 0;
    }

    if (_verbose == 0) {
    } else if (_classification) {
      for (const auto [cls, count] : _class_count) {
        cerr << count << " instances of class '" << cls << "'\n";
      }
    } else {
      cerr << _acc_activity.n() << " activity values btw " << _acc_activity.minval() <<
              " and " << _acc_activity.maxval() << " ave " << _acc_activity.average() << '\n';
    }
  }

  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    if (! ReadSmiles(fname)) {
      cerr << "Cannot read smiles from '" << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Predictions::ReadExperimentalData(IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Predictions::ReadExperimentalData:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadExperimentalData(input);
}

int
Predictions::ReadExperimentalData(iwstring_data_source& input) {
  const_IWSubstring buffer;
  if (!input.next_record(buffer)) {
    cerr << "Predictions::ReadExperimentalData:cannot read header\n";
    return 0;
  }

  if (! DetermineResponse(buffer)) {
    return 0;
  }

  while (input.next_record(buffer)) {
    if (! ReadExperimentalDataRecord(buffer)) {
      cerr << "Predictions::ReadExperimentalData:invalid '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Predictions::ReadExperimentalDataRecord(const const_IWSubstring& buffer) {

  IWString id, expt;
  const_IWSubstring token;
  int i = 0;
  for (int col = 0; buffer.nextword(token, i, _input_token_sep); ++col) {
    if (col == _identifier_column) {
      id = token;
    } else if (col == _expt_column_data) {
      expt = token;
    }
  }

  if (id.empty() || expt.empty()) {
    cerr << "Predictions::ReadExperimentalDataRecord:empty id '" << id << "' or data '" << expt << "'\n";
    return 0;
  }

  if (_classification) {
    _id_to_class.emplace(id, expt);
    _class_count[expt] += 1;
    return 1;
  }
  
  float act;
  if (! expt.numeric_value(act)) {
    cerr << "Predictions::ReadExperimentalDataRecord:invalid numeric '" << expt << "'\n";
    return 0;
  }

  _acc_activity.extra(act);

  _id_to_activity.emplace(id, act);

  return 1;
}

int
Predictions::DetermineResponse(const const_IWSubstring& header) {
  int i = 0;
  const_IWSubstring token;

  for (int col = 0; header.nextword(token, i, _input_token_sep); ++col) {
    if (col == _expt_column_data) {
      _response = token;
      return 1;
    }
  }

  cerr << "Predictions::DetermineResponse:cannot find response in '" << header << "'\n";
  return 0;
}

int
Predictions::ReadSmiles(IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Predictions::ReadSmiles:Cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadSmiles(input);
}

int
Predictions::ReadSmiles(iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! ReadSmilesLine(buffer)) {
      cerr << "Predictions::ReadSmiles:invalid smiles '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Predictions::ReadSmilesLine(const const_IWSubstring& buffer) {
  static constexpr char kSep = ' ';
  int i = 0;
  const_IWSubstring token;
  IWString smiles, id;
  for (int col = 0; buffer.nextword(token, i, kSep); ++col) {
    if (col == 0) {
      smiles = token;
    } else if (col == 1) {
      id = token;
    }
  }

  if (smiles.empty() || id.empty()) {
    return 0;
  }
  
  _id_to_smiles[id] = smiles;

  return 1;
}

// Final output is being done. If we have the smiles for `id`
// write the smiles and `id`, otherwise just write `id`.
int
Predictions::MaybeWriteSmiles(IWString& id,
                              char sep,
                              IWString_and_File_Descriptor& output) const {
  if (_id_to_smiles.empty()) {
    output << id;
    return 1;
  }

  const auto iter = _id_to_smiles.find(id);
  if (iter == _id_to_smiles.end()) {
    cerr << "Predictions::MaybeWriteSmiles:no smiles for '" << id << "'\n";
    return 0;
  }

  output << iter->second << sep << id;

  return 1;
}

int
Predictions::Accumulate(const char* fname) {
  ++_files_read;

  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Predictions::Accumulate:cannot open '" << fname << "'\n";
    return 0;
  }

  return Accumulate(input);
}

int
Predictions::Accumulate(iwstring_data_source& input) {
  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "Predictions::Accumulate:cannot read header record\n";
    return 0;
  }

  while (input.next_record(buffer)) {
    if (! AccumulateLine(buffer)) {
      cerr << "Predictions::Accumulate:invalid data '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Predictions::AccumulateLine(const const_IWSubstring& buffer) {

  IWString id, expt;
  const_IWSubstring token;
  int i = 0;
  for (int col = 0; buffer.nextword(token, i, _input_token_sep); ++col) {
    if (col == _identifier_column) {
      id = token;
    } else if (col == _pred_column) {
      expt = token;
    }
  }

  if (id.empty() || expt.empty()) {
    cerr << "Predictions::AccumulateLine:id '" << id << "' or data '" << expt << "' empty\n";
    return 0;
  }

  if (_classification) {
    return AddClassification(id, expt);
  } else {
    return AddContinuous(id, expt);
  }
}

int
Predictions::AddClassification(const IWString& id,
                               const IWString& expt) {
  auto iter = _id_to_pred_class.find(id);
  if (iter == _id_to_pred_class.end()) {
    auto [iter, _] = _id_to_pred_class.emplace(id, IW_STL_Hash_Map_int());
    iter->second[expt] = 1;
  } else {
    iter->second[expt] += 1;
  }

  return 1;
}

int
Predictions::AddContinuous(const IWString& id,
                           const IWString& expt) {
  float act;
  if (! expt.numeric_value(act)) {
    cerr << "Predictions::AddContinuous:invalid numeric '" << expt << "'\n";
    return 0;
  }

  auto iter = _id_to_pred_continuous.find(id);
  if (iter == _id_to_pred_continuous.end()) {
    auto [iter, _] = _id_to_pred_continuous.emplace(id, Accumulator<double>());
    iter->second << act;
  } else {
    iter->second << act;
  }

  return 1;
}

int
Predictions::Process(IWString_and_File_Descriptor& output) {
  if (_classification) {
    return ProcessClassification(output);
  } else {
    return ProcessContinuous(output);
  }
}

struct IdDiff {
  public:
    IWString _id;

    // Average absolute diff from experimental.
    float _diff;

  public:
    IdDiff();

};

IdDiff::IdDiff() {
  _diff = 0.0f;
}

class ContinuousSorter {
  private:
  public:
    int operator() (const IdDiff& idd1, const IdDiff& idd2) {
      if (idd1._diff < idd2._diff) {
        return 1;
      } else if (idd1._diff > idd2._diff) {
        return -1;
      } else {
        return 0;
      }
    };
};


int
Predictions::ProcessContinuous(IWString_and_File_Descriptor& output) {
  const int n = _id_to_pred_continuous.size();
  std::unique_ptr<IdDiff[]>data = std::make_unique<IdDiff[]>(n);

  auto ndx = 0;
  for (const auto& [id, values] : _id_to_pred_continuous) {
    auto iter = _id_to_activity.find(id);
    if (iter == _id_to_activity.end()) {
      cerr << "Predictions::ProcessContinuous:no data for '" << id << "'\n";
      return 0;
    }
    const float expt = iter->second;

    // For each predicted value, load the abs diff into the accumulator.
    data[ndx]._id = id;
    data[ndx]._diff = abs(values.average() - expt);
    ++ndx;
  }

  ContinuousSorter sorter;
  iwqsort(data.get(), n, sorter);

  char osep = ' ';

  if (! _id_to_smiles.empty()) {
    osep = ',';
    output << "Smiles" << osep;
  }
  output << "Id"
         << osep << _response 
         << osep << "Predictions" 
         << osep << "Minval" 
         << osep << "Maxval" 
         << osep << "Ave" 
         << osep << "AveDiff" 
         << '\n';

  for (int i = 0; i < n; ++i) {
    MaybeWriteSmiles(data[i]._id, osep, output);
    output << osep << _id_to_activity[data[i]._id];

    // First echo stats on the predictions.
    const auto iter = _id_to_pred_continuous.find(data[i]._id);
    const Accumulator<double> acc = iter->second;

    output << osep << acc.n()
           << osep << acc.minval()
           << osep << acc.maxval()
           << osep << acc.average()
           << osep << data[i]._diff
           << '\n';
    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

struct ClassificationSummary {
  public:
    IWString _id;
    int _number_predictions;
    int _correct;
    float _fraction_correct;

  public:
    ClassificationSummary();
};

ClassificationSummary::ClassificationSummary() {
  _number_predictions = 0;
  _correct = 0;
  _fraction_correct = 0.0f;
}

class ClassificationSorter {
  private:
  public:
    int operator() (const ClassificationSummary& cs1, const ClassificationSummary& cs2) {
      if (cs1._fraction_correct < cs2._fraction_correct) {
        return 1;
      } else if (cs1._fraction_correct > cs2._fraction_correct) {
        return -1;
      } else {
        return 0;
      }
    }
};

int
Predictions::ProcessClassification(IWString_and_File_Descriptor& output) {
  const int n = _id_to_pred_class.size();

  std::unique_ptr<ClassificationSummary[]>data = std::make_unique<ClassificationSummary[]>(n);

  int ndx = 0;
  for (const auto& [id, values] : _id_to_pred_class) {
    auto iter = _id_to_class.find(id);
    if (iter == _id_to_class.end()) {
      cerr << "Predictions::ProcessClassification:no data for '" << id << "'\n";
      return 0;
    }

    const IWString& correct = iter->second;

    // Across predicted values get number of times correct.
    data[ndx]._id = id;
    data[ndx]._number_predictions = values.size();
    for (const auto [cls, assigned] : _id_to_pred_continuous) {
      if (cls == correct) {
        ++data[ndx]._correct;
      }
    }
    data[ndx]._fraction_correct = iwmisc::Fraction<float>(data[ndx]._correct, static_cast<int>(values.size()));

    ++ndx;
  }

  ClassificationSorter sorter;
  iwqsort(data.get(), n, sorter);

  char osep = ' ';

  if (! _id_to_smiles.empty()) {
    osep = ',';
    output << "Smiles" << osep;
  }
  output << "Id"
         << osep << _response 
         << osep << "Predictions" 
         << osep << "FractionCorrect" 
         << '\n';

  for (int i = 0; i < n; ++i) {
    MaybeWriteSmiles(data[i]._id, osep, output);
    output << osep << _id_to_class[data[i]._id]
           << osep << data[i]._number_predictions
           << osep << data[i]._fraction_correct
           << '\n';
    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

int
Predictions::Report(std::ostream& output) const {
  output << "Data from " << _files_read << " files\n";
  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:p:CS:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  Predictions predictions;
  if (! predictions.Initialise(cl)) {
    cerr << "Cannot initialise Predictions\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  for (const char* fname : cl) {
    if (! predictions.Accumulate(fname)) {
      cerr << "Cannot read predicted values from '" << fname << "'\n";
      return 1;
    }
  }

  IWString_and_File_Descriptor output(1);
  if (! predictions.Process(output)) {
    cerr << "Error running mismatch determinations\n";
    return 1;
  }
  output.flush();

  if (verbose) {
    predictions.Report(cerr);
  }

  return 0;
}

}  // namespace mispredicted

int
main(int argc, char** argv) {
  int rc = mispredicted::Main(argc, argv);

  return rc;
}
