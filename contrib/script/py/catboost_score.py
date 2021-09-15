import os
import sys
from typing import List

#import catboost as cb
from catboost import CatBoostClassifier, Pool
import numpy as np
from scipy.sparse import csr_matrix

from absl import flags
from absl import app

from lib import class_label_translation_pb2
from lib import gfp_to_svm_lite_pb2
from lib import gfp_model_pb2

FLAGS = flags.FLAGS

flags.DEFINE_string("mdir", None, "Model directory")
flags.DEFINE_string("csr", None, "Stem for csr files _data.csv, _indices.csv, _indptr.csv")
flags.DEFINE_string("ids", None, "File containing the identifiers - samw row order")

def write_classification_result(ids: List[str],
                                response: str,
                                class_number_to_string: List[str],
                                pred:np.array):
  """Write the results of a classification model.
  Args:
    ids: a list if the identifiers
    response: name of the response
    class_number_to_string: mapping from a class number [0,1] to a label.
    pred: a numpy array containing cartboost output.
  Returns:
  """
  sep = ' '
  print(f"ID{sep}pred_{response}{sep}score")
  for i, id in enumerate(ids):
    score = pred[i][1].round(4)
    id = id.rstrip()
    if score <= 0.5:
      print(f"{id}{sep}{class_number_to_string[0]}{sep}{score}")
    else:
      print(f"{id}{sep}{class_number_to_string[1]}{sep}{score}")

def write_regression_result(ids: List[str],
                            response: str,
                            pred: np.array):
  """Write the results of a regression model.
  Args:
    ids: a list if the identifiers
    response: name of the response
    pred: a numpy array containing cartboost output.
  """
  sep = ' '
  print(f"ID{sep}pred_{response}")
  for i, id in enumerate(ids):
    id = id.rstrip()
    print(f"{id}{sep}{pred[i]}")

def get_class_label_translation(mdir: str, \
                                model: gfp_model_pb2.CatboostModel) -> \
                                class_label_translation_pb2.ClassLabelTranslation:
  """Return the ClassLabelTranslation proto in `mdir`.
  Args:
    mdir: model directory
    model: proto which will have the name of the cross reference file
  """
  fname = f"{mdir}/{model.metadata.class_label_translation}"
  with open(fname, "rb") as inp:
    result = class_label_translation_pb2.ClassLabelTranslation()
    if not result.ParseFromString(inp.read()):
      raise ValueError(f"Cannot parse class label translation proto {fname}")

  return result

def get_bit_xref(mdir, proto: gfp_model_pb2.CatboostModel) -> gfp_to_svm_lite_pb2.GfpBitToFeature:
  """Read the bit_xref attribute from `proto` and return the GfpBitToFeature proto.
  Args:
    mdir: model directory
    proto: the CatboostModel that is assumed to be in `mdir`.
  Returns:
  """
  fname = os.path.join(mdir, proto.bit_xref)
  result = gfp_to_svm_lite_pb2.GfpBitToFeature()
  with open(fname, "rb") as inp:
    if not result.ParseFromString(inp.read()):
      raise ValueError(f"Cannot parse bit xref proto {fname}")

  return result

def get_model_proto(mdir) -> gfp_model_pb2.CatboostModel:
  """Return the GfpModel that is in mdir/model.dat.
  Args:
    mdir: model directory
  """
  fname = os.path.join(mdir, "model.dat")
  result = gfp_model_pb2.CatboostModel()
  with open(fname, "rb") as inp:
    if not result.ParseFromString(inp.read()):
      raise ValueError(f"Cannot read model file {fname}")

  return result

def hash_to_list(class_label_translation: class_label_translation_pb2.ClassLabelTranslation) -> List[str]:
  """Convert `class_label_translation` to a list.
  Args:
   class_label_translation:
  Returns:
   List where the to_numeric attribute of each class is the index
  """

  result = ["0", "1"]
  for (key, value) in class_label_translation.to_numeric.items():
    result[value] = key

  return result
    

def catboost_score(argv):
  """Evaluates a catboost model"""

  mdir = FLAGS.mdir
  csr_stem = FLAGS.csr
  id_file = FLAGS.ids

  model = get_model_proto(mdir)
  bit_xref = get_bit_xref(mdir, model)
  nfeatures = bit_xref.params.highest_feature_number
  print(f"Highest feature {nfeatures}", file=sys.stderr)

  fname = f"{csr_stem}_data.csv"
  data = np.loadtxt(fname)
  fname = f"{csr_stem}_indptr.csv"
  indptr = np.loadtxt(fname, dtype=np.int32)
  fname = f"{csr_stem}_indices.csv"
  indices = np.loadtxt(fname, dtype=np.int32)
  print(f"Sizes {data.size} {indptr.size} {indices.size}", file=sys.stderr)
  if data.size != indices.size:
    raise ValueError(f"Data size mismatch {data.size} vs {indices.size}")
  with open(id_file, "r") as inp:
    ids = inp.readlines()

# if len(ids) != indptr.size:
#   raise ValueError(f"Size mismatch btw ids {len(ids)} and data {indptr.size}")

  csr = csr_matrix(((data, indices, indptr)), shape=(indptr.size - 1, nfeatures + 1))

  data = Pool(csr)

  classifier = CatBoostClassifier()

  model_file_name = f"{mdir}/Catboost.model.bin"
  classifier.load_model(model_file_name)

  response = model.metadata.response_name

  if model.metadata.class_label_translation:
    class_label_translation = get_class_label_translation(mdir, model)
    class_number_to_string = hash_to_list(class_label_translation)
    pred = classifier.predict(data, prediction_type='Probability')
    write_classification_result(ids, response, class_number_to_string, pred)
  else:
    pred = classifier.predict(data, prediction_type='Probability')
    write_regression_result(ids, response, pred)


if __name__ == '__main__':
  flags.mark_flag_as_required("mdir")
  flags.mark_flag_as_required("csr")
  flags.mark_flag_as_required("ids")
  app.run(catboost_score)
