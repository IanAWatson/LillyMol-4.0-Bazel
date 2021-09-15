"""Build a catboost model"""

import os
import pathlib
import sys
from typing import Callable, Dict, List

from absl import flags
from absl import app

from catboost import Pool, CatBoostClassifier, CatBoostRegressor
from catboost import EFstrType

from google.protobuf import text_format

import numpy as np


from lib import catboost_options_pb2

FLAGS = flags.FLAGS

flags.DEFINE_string("mdir", None, "Model directory")
flags.DEFINE_string("svml", None, "libsvm input file")
flags.DEFINE_boolean("classification", False, "Build a classification model")
flags.DEFINE_string("option_file", None, "File containing options for learning. CatboostOptions text proto")

# Different options need different conversions.
flags.DEFINE_multi_string("str_option", None, "File containing string options for learning. CatboostOptions text proto")
flags.DEFINE_multi_string("int_option", None, "File containing integer options for learning. CatboostOptions text proto")
flags.DEFINE_multi_string("flt_option", None, "File containing float options for learning. CatboostOptions text proto")

flags.DEFINE_boolean("verbose", False, "Verbose output")

def get_catboost_options(fname: str) -> Dict[str, str]:
  """Read the CatboostOptions proto in `fname` and
  return a Dict containing the contents.
  Args:
   fname: text file containing a CatboostOptions proto
  Returns
    Dict yformed from the options map of the proto.
  """
  with open(fname, "rb") as inp:
    serialised = inp.read()
  proto = text_format.Parse(serialised, catboost_options_pb2.CatboostOptions())

  result = {}
  for k,v in proto.string_option.items():
    result[k] = v
  for k,v in proto.int_option.items():
    result[k] = v
  for k,v in proto.float_option.items():
    result[k] = v

  return result

def write_feature_importance(model: int,
                             verbose: bool,
                             fname: str):
  """Write sorted feature importance data to `fname`.
  Args:
   model: a Catboost model
   verbose: verbosity
   fname: destination
  """
  importance = model.get_feature_importance(type=EFstrType.FeatureImportance)
  if verbose:
    print(f"Feature importance returned {importance.size} values", file=sys.stderr)

  with open(fname, "w") as outp:
    nz = np.nonzero(importance)
    ordered = np.argsort(importance[nz[0]])
    print("Feature,Importance", file=outp)
    for i in ordered[::-1]:
      print(f"{nz[0][i]},{round(importance[[nz[0][i]]][0], 4)}", file=outp)

def add_into_kwargs(options: List[str],
                    converter: Callable,
                    kwargs: Dict):
  """
  """
  if options is None:
    return

  for opt in options:
    f = opt.split('=')
    kwargs[f[0]] = converter(f[1])

def catboost_build(unused_argv):
  """Builds a catboost model using a provided svml file.
  """
  del unused_argv

  mdir = FLAGS.mdir
  svml = FLAGS.svml
  classification = FLAGS.classification
  options_file = FLAGS.option_file
  verbose = FLAGS.verbose

  # Fails if this is an existing file (not dir).
  if not os.path.isdir(mdir):
    pathlib.Path(mdir).mkdir(parents=True, exist_ok=True)

  kwargs = {}
  if options_file:
    kwargs = get_catboost_options(options_file)

  add_into_kwargs(FLAGS.str_option, str, kwargs)
  add_into_kwargs(FLAGS.int_option, int, kwargs)
  add_into_kwargs(FLAGS.flt_option, float, kwargs)
  if not "train_dir" in kwargs:
    kwargs["train_dir"] = mdir

  if not "min_data_in_leaf" in kwargs:
    kwargs["min_data_in_leaf"] = 2

  if verbose:
    print(f"Options {kwargs}", file=sys.stderr)

  train_data = Pool(f"libsvm://{svml}")

  if classification:
    model = CatBoostClassifier(**kwargs)
  else:
    model = CatBoostRegressor(**kwargs)

  model.fit(train_data)

  fname = os.path.join(mdir, "Catboost.model.bin")
  model.save_model(fname)
  fname = os.path.join(mdir, "Catboost.model.cpp")
  model.save_model(fname, format="cpp")

  fname = os.path.join(mdir, "feature_importance.csv")
  write_feature_importance(model, verbose, fname)

if __name__ == '__main__':
  flags.mark_flag_as_required("mdir")
  flags.mark_flag_as_required("svml")
  app.run(catboost_build)

