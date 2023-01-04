import os
import sys
from typing import List

#import catboost as cb
from catboost import CatBoostClassifier, Pool
import numpy as np
from scipy.sparse import csr_matrix

from absl import app
from absl import flags
from absl import logging

from lib import class_label_translation_pb2
from lib import feature_scaling_pb2
from lib import gfp_to_svm_lite_pb2
from lib import gfp_model_pb2

FLAGS = flags.FLAGS

flags.DEFINE_string("mdir", None, "Model directory")
flags.DEFINE_string("svml", None, "Name of libsvm input file to be scored")
flags.DEFINE_string("ids", None, "File containing the identifiers - same row order")

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
  # Might be a good idea to prefix the predicted value, but that can break other things.
  # print(f"ID{sep}pred_{response}{sep}score")
  print(f"ID{sep}{response}{sep}score")
  for i, id in enumerate(ids):
    score = pred[i][1].round(4)
    id = id.rstrip()
    if score <= 0.5:
      print(f"{id}{sep}{class_number_to_string[0]}{sep}{score}")
    else:
      print(f"{id}{sep}{class_number_to_string[1]}{sep}{score}")

def write_regression_result(ids: List[str],
                            response: str,
                            feature_scaling: feature_scaling_pb2.FeatureScaling,
                            pred: np.array):
  """Write the results of a regression model.
  Args:
    ids: a list if the identifiers
    response: name of the response
    feature_scaling: unscale numeric results back to original rangeback to original range
    pred: a numpy array containing cartboost output.
  """
  sep = ' '
# print(feature_scaling)
  pred = feature_scaling.min + pred * (feature_scaling.max - feature_scaling.min)
  # print(f"ID{sep}pred_{response}")
  print(f"ID{sep}{response}")
  for i, id in enumerate(ids):
    id = id.rstrip()
    print(f"{id}{sep}{pred[i]}")

def get_response_scaling(mdir: str, model: gfp_model_pb2.CatboostModel) -> feature_scaling_pb2.FeatureScaling:
  """Return the FeatureScaling prpto in `mdir`. 
  Args:
    mdir: model directory
    model: proto which will have the name of the feature scaling file
  Returns:
    FeatureScaling proto.
  """
  fname = os.path.join(mdir, model.metadata.response_scaling)
  with open(fname, "rb") as inp:
    result = feature_scaling_pb2.FeatureScaling()
    if not result.ParseFromString(inp.read()):
      raise ValueError(f"Canot parse response scaling proto {fname}")

    return result

def get_class_label_translation(mdir: str, \
                                model: gfp_model_pb2.CatboostModel) -> \
                                class_label_translation_pb2.ClassLabelTranslation:
  """Return the ClassLabelTranslation proto in `mdir`.
  Args:
    mdir: model directory
    model: proto which will have the name of the cross reference file
  Returns:
    ClassLabelTranslation proto.
  """
  fname = os.path.join(mdir, model.metadata.class_label_translation)
  with open(fname, "rb") as inp:
    result = class_label_translation_pb2.ClassLabelTranslation()
    if not result.ParseFromString(inp.read()):
      raise ValueError(f"Cannot parse class label translation proto {fname}")

    return result

def get_model_proto(mdir) -> gfp_model_pb2.CatboostModel:
  """Return the GfpModel that is in mdir/model.dat.
  Args:
    mdir: model directory
  """
  fname = os.path.join(mdir, "model.dat")
  with open(fname, "rb") as inp:
    result = gfp_model_pb2.CatboostModel()
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
  id_file = FLAGS.ids
  svml_file = FLAGS.svml

  model = get_model_proto(mdir)

  with open(id_file, "r") as inp:
    ids = inp.readlines()
  logging.vlog(1, f"Read {len(ids)} identifiers from {id_file}")

  data = Pool(f"libsvm://{svml_file}")

  classifier = CatBoostClassifier()

  model_file_name = os.path.join(mdir, "Catboost.model.bin")
  classifier.load_model(model_file_name)

  response = model.metadata.response_name

  if model.metadata.class_label_translation:
    class_label_translation = get_class_label_translation(mdir, model)
    class_number_to_string = hash_to_list(class_label_translation)
    pred = classifier.predict(data, prediction_type='Probability')
    write_classification_result(ids, response, class_number_to_string, pred)
  else:
    scaling = get_response_scaling(mdir, model)
    pred = classifier.predict(data, prediction_type="RawFormulaVal")
    write_regression_result(ids, response, scaling, pred)


if __name__ == '__main__':
  flags.mark_flag_as_required("mdir")
  flags.mark_flag_as_required("svml")
  flags.mark_flag_as_required("ids")
  app.run(catboost_score)
