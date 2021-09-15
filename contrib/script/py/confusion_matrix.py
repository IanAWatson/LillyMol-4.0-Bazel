# Compute various statistics on a classification result.
# INput is a Pandas dataframe that must contain at least two
# columns:
#  Truth label
#  Predicted label
# May optionally contain Precited_Score which enables various computations

import os
import sys
from typing import List

from absl import app
from absl import flags
from absl import logging
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import precision_score, recall_score, f1_score, roc_auc_score, confusion_matrix
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import accuracy_score, cohen_kappa_score, matthews_corrcoef
from sklearn.metrics import classification_report
from tabulate import tabulate

from lib import class_label_translation_pb2

FLAGS = flags.FLAGS
flags.DEFINE_integer("tcol", 2, "Column in which truth is found")
flags.DEFINE_integer("pcol", 3, "Column in which prediction is found")
flags.DEFINE_integer("scol", -1, "Column in which score is found")
flags.DEFINE_string("roc_curve" , None, "Input token separator")
flags.DEFINE_string("sep", ' ', "Input token separator")
flags.DEFINE_string("cltrans", None, "Name of a ClassLabelTranslation binary proto file")

def make_roc_curve(response: str,
                   y_test: str,
                   y_score: str,
                   fname: str):
  """Generate an ROC curve and plot to `fname`.
  Taken from
  https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc.html
  Args:
    response: name of response
    y_test: truth values for the test set
    y_score: score values for the test set
    fname: .png file to be created
  """

  fpr, tpr, _ = roc_curve(y_test, y_score)
  roc_auc = auc(fpr, tpr)

  plt.figure()
  lw = 2
  plt.plot(fpr, tpr, color="darkorange", lw=lw, label="ROC curve (area = %0.3f)" % roc_auc)
  plt.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
  plt.xlim([0.0, 1.0])
  plt.ylim([0.0, 1.05])
  plt.xlabel("False Positive Rate")
  plt.ylabel("True Positive Rate")
  plt.title(f"Receiver operating characteristic {response}")
  plt.legend(loc="lower right")
  plt.show()


def get_class_label_translation(fname: str) -> class_label_translation_pb2.ClassLabelTranslation:
  """Read a binary ClassLabelTranslation proto from `fname`.
  Args:
    fname: file containing proto
  Returns:
    ClassLabelTranslation
  """
  with open(fname, "rb") as inp:
    result = class_label_translation_pb2.ClassLabelTranslation()
    if not result.ParseFromString(inp.read()):
      raise ValueError(f"Cannot read ClassLabelTranslation {fname}")

    return result

def do_class_label_translation(class_label_translation: class_label_translation_pb2.ClassLabelTranslation,
                               data:int) -> List[int]:
  """Apply the string->int transformation in `class_label_translation` to
  the values in `data`.
  Args:
    class_label_translation: mapping from class names to numbers
  Returns:
    transformed `data` string->int
  """
  result = []
  for value in data:
    as_int = class_label_translation.to_numeric.get(value)
    if as_int is None:
      raise ValueError("No class number for '%s'", value)
    result.append(as_int)

  return result

def classification_statistics(argv):
  """Read argv[0] as Pandas DataFrame and compute classification metrics."""
  if len(argv) == 1:
    logging.fatal("Must specify file to process")

  tcol = FLAGS.tcol - 1
  pcol = FLAGS.pcol - 1
  scol = FLAGS.scol - 1

  if tcol == pcol or tcol == scol or tcol == scol:
    logging.fatal("Duplicated column numbers tcol %d pcol %d scol %d", tcol, pcol, scol)

  logging.info( "Columns tcol %d pcol %d", tcol+1, pcol+1)

  data = pd.read_csv(argv[1], header=0, sep=FLAGS.sep)
  nrows, ncols = data.shape
  logging.info("Shape %d rows %d cols", nrows, ncols)

  if tcol >= ncols or pcol >= ncols:
    logging.fatal("Column(s) out of range truth %d pred %d ncols %d", tcol+1, pcol+1, ncols)

  truth = data.columns[tcol]
  logging.info("Truth is '%s' pred is '%s'", truth, data.columns[pcol])
  data.reset_index(inplace=True, drop=True)
  if scol >= 0:
    logging.info("Predicted score in column %d", scol)

  y_true = data.iloc[:,tcol]
  y_pred = data.iloc[:,pcol]

  if FLAGS.cltrans:
    class_label_translation = get_class_label_translation(FLAGS.cltrans)
    y_true = do_class_label_translation(class_label_translation, y_true)
    y_pred = do_class_label_translation(class_label_translation, y_pred)

  precision = precision_score(y_true, y_pred, average='binary')
  recall = recall_score(y_true, y_pred, average='binary')
  f1 = f1_score(y_true, y_pred, average='binary')
  accuracy = accuracy_score(y_true, y_pred)
  roc = roc_auc_score(y_true, data.iloc[:,scol])
  kappa = cohen_kappa_score(y_true, y_pred)
  mcc = matthews_corrcoef(y_true, y_pred)

  table = [[accuracy, roc, recall, precision, f1, kappa, mcc]]
  print(tabulate(table, headers=["Accuracy", "AUC", "Recall", "Prec", "F1", "Kappa", "MCC"], floatfmt=".4f"))

# print(f"Acc {accuracy:.3f} precision {precision} recall {recall} f1 {f1} roc {roc}")

  cfm = confusion_matrix(y_true, y_pred).ravel()
  print(cfm)
  table = [["positive", cfm[0], cfm[2]], ["negative", cfm[1], cfm[3]]]
  print(tabulate(table, headers = ["", "Pred.Pos", "Pred.Neg"]))

  print("classification_report")
  print(classification_report(y_true, y_pred, target_names=[data.columns[tcol], data.columns[pcol]]))

  # Plot ROC curve
  if FLAGS.roc_curve:
    make_roc_curve(truth, y_true, data.iloc[:,scol], FLAGS.roc_curve)


if __name__ == '__main__':
  app.run(classification_statistics)
