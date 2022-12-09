// Create stratified samples for either continuous or
// class data.

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

namespace stratified_sample {

using std::cerr;

void
Usage(int rc) {
  cerr << "Generates stratified train/test splits as smiles files\n";
  cerr << "Single argument is a smiles file\n";
  cerr << " -A <fname>        activity file name (id activity)  must be present\n";
  cerr << " -i <sep>          token separator in -A file (def space)\n";
  cerr << " -c <col>          activity data is in column <col> of the -A file (def 2)\n";
  cerr << " -n <nsplits>      the number of splits to create\n";
  cerr << " -t <fraction>     fraction for the training set (or percent)\n";
  cerr << " -R <fname>        file name stem for training set files to create\n";
  cerr << " -E <fname>        file name stem for training set files to create\n";
  cerr << " -C                data is class data\n";
  cerr << " -b <nbuckets>     number of buckets for stratification (def 10)\n";
  cerr << "                   smaller number might be needed for smaller datasets\n";
  cerr << "                   The risk is that a bucket with few items may always get everything selected\n";
  cerr << " -a                create split Activity files in addition to split smiles files\n";
  cerr << " -v                verbose output\n";

  ::exit(rc);
}

// We need something that can be sorted and shuffled.
template <typename T>
struct IdValue {
  // The string must remain in scope.
  const IWString* id;
  // Class name or number.
  T value;
  int bucket;

  public:
    IdValue();
};

template <typename T>
IdValue<T>::IdValue() {
  id = nullptr;
  value = {};
  bucket = -1;
}

template <typename T>
class ByValueSorter {
  private:
  public:
    int operator() (const IdValue<T>& idv1, const IdValue<T>& idv2) const {
      if (idv1.value < idv2.value) {
        return -1;
      } else if (idv1.value > idv2.value) {
        return 1;
      } else {
        return 0;
      }
    }
};

template <typename T>
class ByBucketSorter {
  private:
  public:
    int operator() (const IdValue<T>& idv1, const IdValue<T>& idv2) const {
      if (idv1.bucket < idv2.bucket) {
        return -1;
      } else if (idv1.bucket > idv2.bucket) {
        return 1;
      } else {
        return 0;
      }
    }
};


class StratifiedSample {
  private:
    int _verbose;

    int _lines_read;

    int _nsplits;

    // By default, we only create split smiles files, but if requested
    // we can also create train and test activity files.
    bool _write_split_activity_files;

    // The response name read from the -A file.
    IWString _response;

    float _train_fraction;

    char _input_file_separator;

    int _activity_column;

    // There are two ways of assigning buckets. Based on the activity
    // value relative to its position in the range of activities, or
    // within a sorted list of activity values, its index;
    int _bucketise_by_value;

    // When processing continuous data, bucketize the activity
    // data and sample from those buckets.
    int _nbuckets;

    bool _is_classification;

    std::mt19937_64 _rng;

    // If we have a continuous response.
    IW_STL_Hash_Map_float _id_to_activity;
    // If we have class data.
    IW_STL_Hash_Map_String _id_to_class;

    // We translate from class name to a number.
    IW_STL_Hash_Map_int _class_name_to_number;
    // Translation from class number back to starting name.
    std::unique_ptr<IWString[]> _class_number_to_name;

    // If we have smiles, store them here.
    IW_STL_Hash_Map_String _id_to_smiles;

    IWString _train_file_name_stem;
    IWString _test_file_name_stem;

    // For each ID, how many times is it in the training set.
    IW_STL_Hash_Map_int _times_in_train;

  // private functions.

    int ReadActivityDataRecord(const const_IWSubstring& buffer);
    int ReadActivityData(IWString& fname);
    int ReadActivityData(iwstring_data_source& input);

    int ReadSmiles(const char* fname);
    int ReadSmiles(iwstring_data_source& input);
    int ReadSmilesRecord(const const_IWSubstring& line);

    int ProcessContinuous();
    int ProcessClassification();
    int ProcessClassification(const IW_STL_Hash_Map_int& class_count);
    int ProcessContinuous(const Accumulator<double>& acc_activity);
    int ProcessClassification(IdValue<int>* id_value);
    int ProcessContinuous(IdValue<double>* id_value);

    int BucketiseByValue(const Accumulator<double>& acc_activity,
                        IdValue<double>* id_value);
    int BucketiseByPosition(IdValue<double>* id_value);

    template <typename T>
    int CreateSplit(IdValue<T>* id_value, int ndx);

    template <typename T> void ShuffleWithinBuckets(IdValue<T>* id_value,
                                const int n);

    template <typename T>
    int WriteSplit(const IdValue<T>* id_value, int ndx, const int* train);
    template <typename T>
    int WriteValuesContinuous(const IdValue<T>* id_value,
                              const int* train,
                              int flag,
                              IWString& fname) const;
    template <typename T>
    int WriteValuesClassification(const IdValue<T>* id_value,
                              const int* train,
                              int flag,
                              IWString& fname) const;
    template <typename T>
    int WriteValuesClassification(const IdValue<T>* id_value,
                              const int* train,
                              int flag,
                              IWString_and_File_Descriptor& output) const;
    template <typename T>
    int WriteValuesContinuous(const IdValue<T>* id_value,
                              const int* train,
                              int flag,
                              IWString_and_File_Descriptor& output) const;
    template <typename T> int WriteSmiles(const IdValue<T>* id_value,
                        const int* train,
                        int flag,
                        IWString& fname);
    template <typename T> int WriteSmiles(const IdValue<T>* id_value,
                        const int* train,
                        int flag,
                        IWString_and_File_Descriptor& output);

  public:
    StratifiedSample();

    int Initialise(Command_Line& cl);

    int Report(std::ostream& output) const;

    int Process();

};

StratifiedSample::StratifiedSample() {
  _verbose = 0;
  _is_classification = false;
  _activity_column = 1;
  _input_file_separator = ' ';
  _train_file_name_stem = "train";
  _test_file_name_stem = "test";
  _write_split_activity_files = false;

  _train_fraction = 0.8;
  _nsplits = 10;

  _nbuckets = 10;
  _bucketise_by_value = false;

  std::random_device rd;
  _rng.seed(rd());
}

int
StratifiedSample::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (! cl.option_present('A')) {
    cerr << "Must specify activity file via the -A option\n";
    Usage(1);
  }

  if (cl.option_present('C')) {
    _is_classification = true;
    if (_verbose) {
      cerr << "Data is class data\n";
    }
  }

  if (cl.option_present('c')) {
    if (! cl.value('c', _activity_column) || _activity_column < 1) {
      cerr << "StratifiedSample::Initialise:invalid activity column (-c)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Activity values in column " << _activity_column << '\n';
    }
    --_activity_column;
  }

  if (cl.option_present('A')) {
    IWString fname = cl.string_value('A');
    if (! ReadActivityData(fname)) {
      cerr << "Cannot read activity data '" << fname << "'\n";
      return 0;
    }
  }

  if (cl.option_present('b')) {
    _bucketise_by_value = true;
    if (_verbose) {
      cerr << "Will bucketise by value\n";
    }
  }

  if (cl.option_present('f')) {
    if (! cl.value('f', _train_fraction) || _train_fraction <= 0.0f || _train_fraction >= 1.0f) {
      cerr << "The training fraction (-f) must be a valid fraction\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Training fraction " << _train_fraction << '\n';
    }
  }

  if (cl.option_present('n')) {
    if (! cl.value('n', _nsplits) || _nsplits < 1) {
      cerr << "The number of splits (-n) must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will create " << _nsplits << " splits\n";
    }
  }

  if (cl.empty()) {
    cerr << "insufficient arguments\n";
    return 0;
  }

  if (cl.size() > 1) {
    cerr << "Only takes one command line argument\n";
    return 0;

  }

  const char* fname = cl[0];
  if (! ReadSmiles(fname)) {
    if (! ReadSmiles(fname)) {
      cerr << "Cannot read smiles from '" << fname << "'\n";
      return 0;
    }

    // Make sure every smiles id has an associated activity.
    int no_smiles = 0;
    for (const auto& [id, _] : _id_to_activity) {
      const auto iter = _id_to_smiles.find(id);
      if (iter == _id_to_smiles.end()) {
        cerr << "StratifiedSample::Initialise:no smiles for '" << id << "'\n";
        ++no_smiles;
      }
    }
    if (no_smiles) {
      return 0;
    }
  }

  if (cl.option_present('R')) {
    cl.value('R', _train_file_name_stem);
    if (_verbose) {
      cerr << "Train splits prefix " << _train_file_name_stem << '\n';
    }
  }

  if (cl.option_present('E')) {
    cl.value('E', _test_file_name_stem);
    if (_verbose) {
      cerr << "Test splits prefix " << _test_file_name_stem << '\n';
    }
  }

  if (cl.option_present('a')) {
    _write_split_activity_files = true;
    if (_verbose) {
      cerr << "Will write split activity files\n";
    }
  }

  return 1;
}

int
StratifiedSample::Process() {
  if (_is_classification) {
    return ProcessClassification();
  } else {
    return ProcessContinuous();
  }
}

void
WriteClassCounts(const IW_STL_Hash_Map_int& class_count,
                 std::ostream& output) {
  for (auto& [cl, co] : class_count) {
    output << co << " instances of '" << cl << "'\n";
  }
}

int
StratifiedSample::ProcessClassification() {
  if (_id_to_class.empty()) {
    cerr << "StratifiedSample::ProcessClassification:no data\n";
    return 0;
  }

  IW_STL_Hash_Map_int class_count;
  for (const auto& [id, c] : _id_to_class) {
    auto iter = class_count.find(c);
    if (iter == class_count.end()) {
      class_count.emplace(c, 1);
    } else {
      ++iter->second;
    }
  }

  if (class_count.size() != 2) {
    cerr << "StratifiedSample::ProcessClassification:incorrect class count " << class_count.size() << '\n';
    WriteClassCounts(class_count, cerr);
    return 0;
  }

  if (_verbose) {
    cerr << "StratifiedSample::ProcessClassification:class counts\n";
    WriteClassCounts(class_count, cerr);
  }

  return ProcessClassification(class_count);
}

int
StratifiedSample::ProcessClassification(const IW_STL_Hash_Map_int& class_count) {
  assert(class_count.size() == 2);

  // Since there are just two classes, do this the quick and dirty way.
  IWString majority_class_name, minority_class_name;
  int c1 = 0;
  int c2 = 0;
  int ndx = 0;
  for (const auto& [cname, count] : class_count) {
    if (ndx == 0) {
      majority_class_name = cname;
      c1 = count;
    } else {
      minority_class_name = cname;
      c2 = count;
    }
    ++ndx;
  }

  // This absolutely should not happen.
  if (c1 == 0 || c2 == 0) {
    cerr << "StratifiedSample::ProcessClassification:huh, zero class count???\n";
    return 0;
  }

  if (c1 < c2) {
    std::swap(majority_class_name, minority_class_name);
    std::swap(c1, c2);
  }

  _class_number_to_name = std::make_unique<IWString[]>(2);
  _class_number_to_name[0] = majority_class_name;
  _class_number_to_name[1] = minority_class_name;
  _class_name_to_number[majority_class_name] = 0;
  _class_name_to_number[minority_class_name] = 0;

  const uint32_t nvalues = _id_to_activity.size();

  std::unique_ptr<IdValue<int>[]> id_value = std::make_unique<IdValue<int>[]>(nvalues);
  ndx = 0;
  for (const auto& [id, cname] : _id_to_class) {
    id_value[ndx].id = &id;
    if (cname == majority_class_name) {
      id_value[ndx].value = 0;
    } else if (cname == minority_class_name) {
      id_value[ndx].value = 1;
    } else {
      cerr << "StratifiedSample::ProcessClassification:unrecognised class name '" << cname << "'\n";
      return 0;
    }
  }

  ByValueSorter<int> sorter;
  iwqsort(id_value.get(), nvalues, sorter);

  for (uint32_t i = 0; i < nvalues; ++i) {
    if (id_value[i].value == 1) {
      id_value[i].bucket = 1;
    } else {
      id_value[i].bucket = 0;
    }
  }

  return ProcessClassification(id_value.get());
}

int
StratifiedSample::ProcessClassification(IdValue<int>* id_value) {
  for (int i = 0; i < _nsplits; ++i) {
    CreateSplit(id_value, i);
  }

  return 1;
}

template <typename T>
int
StratifiedSample::CreateSplit(IdValue<T>* id_value, int ndx) {
  const int n = _id_to_activity.size();

  std::unique_ptr<int[]> in_train(new_int(n));

  ShuffleWithinBuckets(id_value, n);

  const int needed = static_cast<int>(_train_fraction * n);

  int needed_per_bucket = needed / _nbuckets;
  if (needed_per_bucket == 0) {
    needed_per_bucket = 1;
  }

  if (_verbose && ndx == 0) {
    cerr << n << " points, in train " << needed << " per bucket " << needed_per_bucket << '\n';
  }

  in_train[0] = 1;
  int got_this_bucket = 1;
  for (int i = 1; i < n; ++i) {
    if (id_value[i].bucket == id_value[i-1].bucket) {
      if (got_this_bucket < needed_per_bucket) {
        ++got_this_bucket;
        in_train[i] = 1;
      }
    } else {
      got_this_bucket = 1;
      in_train[i] = 1;
    }
  }

  for (int i = 0; i < n; ++i) {
    if (in_train[i] == 0) {
      continue;
    }
    auto iter = _times_in_train.find(*(id_value[i].id));
    if (iter == _times_in_train.end()) {
      _times_in_train.emplace(*(id_value[i].id), 1);
    } else {
      ++iter->second;
    }
  }

  return WriteSplit(id_value, ndx, in_train.get());
}

// Sort `id_value` by bucket, and then shuffle the values within
// each bucket.
template <typename T>
void
StratifiedSample::ShuffleWithinBuckets(IdValue<T>* id_value, const int n) {
  static ByBucketSorter<T> sorter;

  iwqsort(id_value, n, sorter);

  IdValue<T>* bstart = id_value;
  for (int i = 1; i < n; ++i) {
    if (id_value[i].bucket == id_value[i-1].bucket) {
      continue;
    }
    std::shuffle(bstart, id_value + i, _rng);
    bstart = id_value + i;
  }

  std::shuffle(bstart, id_value + n, _rng);
}

template <typename T>
int
StratifiedSample::WriteSmiles(const IdValue<T>* id_value,
                        const int* train,
                        int flag,
                        IWString& fname) {
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "StratifiedSample::WriteSmiles:cannot open '" << fname << "'\n";
    return 0;
  }

  if (_verbose > 1) {
    cerr << "Writing smiles '" << fname << "'\n";
  }

  return WriteSmiles(id_value, train, flag, output);
}

// For those values in `id_value` for which train[i] == flag,
// fetch the corresponding smiles from _id_to_smiles and 
// write a smiles file.
template <typename T>
int
StratifiedSample::WriteSmiles(const IdValue<T>* id_value,
                        const int* train,
                        int flag,
                        IWString_and_File_Descriptor& output) {
  static constexpr char kSep = ' ';

  const uint32_t n = _id_to_activity.size();
  for (uint32_t i = 0; i < n; ++i) {
    if (train[i] != flag) {
      continue;
    }

    const auto iter = _id_to_smiles.find(*id_value[i].id);
    if (iter == _id_to_smiles.end()) {
      cerr << "StratifiedSample::WriteSmiles:huh, where is smiles for '" << (*id_value[i].id) << "'\n";
      return 0;
    }
    output << iter->second << kSep << *id_value[i].id << '\n';
    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

template <typename T>
int
StratifiedSample::WriteValuesContinuous(const IdValue<T>* id_value,
                        const int* train,
                        int flag,
                        IWString_and_File_Descriptor& output) const {
  constexpr char kSep = ' ';

  output << "ID" << kSep << _response << '\n';

  const uint32_t n = _id_to_activity.size();
  for (uint32_t i = 0; i < n; ++i) {
    if (train[i] != flag) {
      continue;
    }
    output << *id_value[i].id << kSep << static_cast<float>(id_value[i].value) << '\n';
  }

  return 1;
}

template <typename T>
int
StratifiedSample::WriteValuesClassification(const IdValue<T>* id_value,
                        const int* train,
                        int flag,
                        IWString_and_File_Descriptor& output) const {
  constexpr char kSep = ' ';

  output << "ID" << kSep << _response << '\n';

  const uint32_t n = _id_to_activity.size();
  for (uint32_t i = 0; i < n; ++i) {
    if (train[i] != flag) {
      continue;
    }

    output << *id_value[i].id << kSep << _class_number_to_name[id_value[i].value] << '\n';
  }

  return 1;
}

template <typename T>
int
StratifiedSample::WriteValuesContinuous(const IdValue<T>* id_value,
                              const int* train,
                              int flag,
                              IWString& fname) const {
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "StratifiedSample::WriteValues:cannot open '" << fname << "'\n";
    return 0;
  }

  if (_verbose > 1) {
    cerr << "StratifiedSample::WriteValuesClassification:writing '" << fname << "'\n";
  }

  return WriteValuesContinuous(id_value, train, flag, output);
}

template <typename T>
int
StratifiedSample::WriteValuesClassification(const IdValue<T>* id_value,
                              const int* train,
                              int flag,
                              IWString& fname) const {
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "StratifiedSample::WriteValues:cannot open '" << fname << "'\n";
    return 0;
  }

  if (_verbose > 1) {
    cerr << "StratifiedSample::WriteValuesClassification:writing '" << fname << "'\n";
  }

  return WriteValuesClassification(id_value, train, flag, output);
}

// Return $(stem)$(ndx).$(suffix)
IWString
FileName(const IWString& stem,
         int ndx,
         const char* suffix) {
  IWString result(stem);
  result << ndx << '.' << suffix;
  return result;
}

template <typename T>
int
StratifiedSample::WriteSplit(const IdValue<T>* id_value, int ndx,
                const int* train) {
  static constexpr int kTrain = 1;
  static constexpr int kTest = 0;

  if (! _id_to_smiles.empty()) {
    IWString fname = FileName(_train_file_name_stem, ndx, "smi");
    if (! WriteSmiles(id_value, train, kTrain, fname)) {
      cerr << "StratifiedSample::WriteSplit:cannot write train split '" << fname << "'\n";
      return 0;
    }

    fname = FileName(_test_file_name_stem, ndx, "smi");
    if (! WriteSmiles(id_value, train, kTest, fname)) {
      cerr << "Cannot write test split '" << fname << "'\n";
      return 0;
    }
  }

  // Nothing more to do unless we are also writing activity files.
  if (! _write_split_activity_files) {
    return 1;
  }

  IWString train_fname = FileName(_train_file_name_stem, ndx, "activity");

  int rc = 0;
  if (_is_classification) {
    rc = WriteValuesClassification(id_value, train, kTrain, train_fname);
  } else {
    rc = WriteValuesContinuous(id_value, train, kTrain, train_fname);
  }

  if (rc == 0) {
    cerr << "StratifiedSample::WriteSplit:cannot write '" << train_fname << "'\n";
    return 0;
  }

  IWString test_fname = FileName(_test_file_name_stem, ndx, "activity");

  if (_is_classification) {
    rc = WriteValuesClassification(id_value, train, kTest, train_fname);
  } else {
    rc = WriteValuesContinuous(id_value, train, kTest, train_fname);
  }

  if (rc == 0) {
    cerr << "StratifiedSample::WriteSplit:cannot write '" << test_fname << "'\n";
    return 0;
  }

  return 1;
}

int
StratifiedSample::ProcessContinuous() {
  if (_id_to_activity.empty()) {
    cerr << "StratifiedSample::ProcessContinuous:no data\n";
    return 0;
  }

  Accumulator<double> acc_activity;
  for (const auto& [id, a] : _id_to_activity) {
    acc_activity.extra(a);
  }

  if (_verbose) {
    cerr << "StratifiedSample::ProcessContinuous: " << _response << ' ' << _id_to_activity.size() << 
            " values btw " << acc_activity.minval() << " and " << acc_activity.maxval() << " ave " << acc_activity.average() << '\n';
  }

  const double range = acc_activity.maxval() - acc_activity.minval();
  if (range <= 0.0) {
    cerr << "StratifiedSample::ProcessContinuous:range is empty\n";
    return 0;
  }

  return ProcessContinuous(acc_activity);
}

int
StratifiedSample::ProcessContinuous(const Accumulator<double>& acc_activity) {
  std::unique_ptr<IdValue<double>[]> id_value = std::make_unique<IdValue<double>[]>(_id_to_activity.size());

  if (_bucketise_by_value) {
    BucketiseByValue(acc_activity, id_value.get());
  } else {
    BucketiseByPosition(id_value.get());
  }

  // Report the number of items per bucket. Reported if verbose output
  // or of there are empty buckets.
  extending_resizable_array<int> items_in_bucket;
  const int n = _id_to_activity.size();
  for (int i = 0; i < n; ++i) {
    ++items_in_bucket[id_value[i].bucket];
  }

  for (int i = 0; i < items_in_bucket; ++i) {
    if (items_in_bucket[i]) {
      if (items_in_bucket == 0 || _verbose) {
        cerr << items_in_bucket[i] << " items in bucket " << i << '\n';
      }
    }
  }

  return ProcessContinuous(id_value.get());
}

// In this bucketisation scheme, the range is bucketised, and each bucket
// value is assigned based on where in the range the activity falls.
// Note that in the common case of a small number of active molecules
// we will get buckets that have hardly any molecules in that bucket.
int
StratifiedSample::BucketiseByValue(const Accumulator<double>& acc_activity,
                        IdValue<double>* id_value) {
  const double range = acc_activity.maxval() - acc_activity.minval();
  const double dx = range / static_cast<double>(_nbuckets);
  const double min_activity = acc_activity.minval();

  int ndx = 0;
  for (const auto& [id, a] : _id_to_activity) {
    id_value[ndx].id = &id;
    id_value[ndx].value = a;
    id_value[ndx].bucket = static_cast<int>((a - min_activity) / dx);
    ++ndx;
  }

  ByValueSorter<double> sorter;
  iwqsort(id_value, ndx, sorter);

  return 1;
}

// The `id_value` array is filled and then sorted by value.
// Buckets are assigned based on array index, so there will
// always be roughly the same number of items in each bucket.
int
StratifiedSample::BucketiseByPosition(IdValue<double>* id_value) {
  int ndx = 0;
  for (const auto& [id, a] : _id_to_activity) {
    id_value[ndx].id = &id;
    id_value[ndx].value = a;
    ++ndx;
  }

  ByValueSorter<double> sorter;
  iwqsort(id_value, ndx, sorter);

  int items_per_bucket = ndx / _nbuckets;
  if (items_per_bucket == 0) {
    items_per_bucket = 1;
  }

  for (int i = 0; i < ndx; ++i) {
    id_value[i].bucket = i / items_per_bucket;
    if (id_value[i].bucket == _nbuckets) {
      id_value[i].bucket = _nbuckets - 1;
    }
  }

  return 1;
}

int
StratifiedSample::ProcessContinuous(IdValue<double>* id_value) {
  for (int i = 0; i < _nsplits; ++i) {
    CreateSplit(id_value, i);
  }

  return 1;
}

int
StratifiedSample::ReadActivityData(IWString& fname) {
  iwstring_data_source input;
  if (! input.open(fname)) {
    cerr << "StratifiedSample::ReadActivityData:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadActivityData(input);
}

int
StratifiedSample::ReadActivityData(iwstring_data_source& input) {
  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "StratifiedSample::ReadActivityData:cannot read header record\n";
    return 0;
  }

  const_IWSubstring token;
  int i = 0;
  for (int col = 0; buffer.nextword(token, i, _input_file_separator); ++col) {
    if (col == _activity_column) {
      _response = token;
      break;
    }
  }

  if (_response.empty()) {
    cerr << "StratifiedSample::ReadActivityData:did not find column " << (_activity_column + 1) << 
            " in " << buffer << '\n';
    return 0;
  }

  while (input.next_record(buffer)) {
    if (! ReadActivityDataRecord(buffer)) {
      cerr << "StratifiedSample::ReadActivityData:invalid data '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
StratifiedSample::ReadActivityDataRecord(const const_IWSubstring& buffer) {
  const_IWSubstring token;
  int i = 0;
  // We must find non empty values for each of these in `buffer`.
  IWString id;
  IWString activity;
  for (int col = 0; buffer.nextword(token, i, _input_file_separator); ++col) {
    if (col == 0) {
      id = token;
    } else if (col == _activity_column) {
      activity = token;
    }
  }

  if (id.empty() || activity.empty()) {
    cerr << "StratifiedSample::ReadActivityDataRecord:did not find id '" << id << "' and/or activity '" << activity << "' in '" << buffer << "'\n";
    return 0;
  }

  // No checking for duplicates, the last entry wins.
  if (_is_classification) {
    _id_to_class[id] = activity;
    return 1;
  } else {
    float a;
    if (! activity.numeric_value(a)) {
      cerr << "StratifiedSample::ReadActivityDataRecord:invalid numeric '" << activity << "'\n";
      return 0;
    }
    _id_to_activity[id] = a;
  }

  return 1;
}

int
StratifiedSample::Report(std::ostream& output) const {
  output << "Generated " << _nsplits << " splits\n";

  extending_resizable_array<int> in_train;
  extending_resizable_array<int> in_train_highest_bucket;
  for (const auto& [_, count] : _times_in_train) {
    ++in_train[count];
  }

  for (int i = 0; i < in_train.number_elements(); ++i) {
    if (in_train[i] == 0) {
      continue;
    }
    output << in_train[i] << " items in train " << i << " times\n";
  }

  return 1;
}

int
StratifiedSample::ReadSmiles(const char* fname) {
  iwstring_data_source input;
  if (! input.open(fname)) {
    cerr << "StratifiedSample::ReadSmiles:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadSmiles(input);
}

int
StratifiedSample::ReadSmiles(iwstring_data_source& input) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (! ReadSmilesRecord(buffer)) {
      cerr << "StratifiedSample::ReadSmilesRecord:invalid smiles '" << buffer << "'\n";
      return 0;
    }
  }

  return _id_to_smiles.size();
}

int
StratifiedSample::ReadSmilesRecord(const const_IWSubstring& line) {
  IWString smiles, id;
  int i = 0;
  if (! line.nextword(smiles, i, ' ') ||
      ! line.nextword(id, i, ' ') ||
      smiles.empty() || id.empty()) {
    cerr << "StratifiedSample::ReadSmilesRecord:empty smiles and/or id\n";
    return 0;
  }

  _id_to_smiles.emplace(id, smiles);

  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vc:n:f:A:R:E:Cab");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  StratifiedSample sampler;
  if (! sampler.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (! sampler.Process()) {
    cerr << "Sampling failed\n";
    return 1;
  }

  if (verbose) {
    sampler.Report(cerr);
  }

  return 0;
}

}  // namespace stratified_sample

int
main(int argc, char** argv) {
  int rc = stratified_sample::Main(argc, argv);

  return rc;
}
