"""Summarises information extracted from one or more iwstats_pb2.Stats
protos.
If given one argument, generates summary statistics on the data in there.
If given two arguments, generates t-test comparisons between the values.
"""

from dataclasses import dataclass
from dataclasses import field, fields
from itertools import count
import math
from typing import Dict, List, Set

from scipy import stats

from absl import app
from absl import flags
from absl import logging

from contrib.script.py.lib import iwstats_pb2

FLAGS = flags.FLAGS
flags.DEFINE_multi_string('pattern', [], 'pattern definining each series of results, will contain %d')
flags.DEFINE_string('osep', ' ', 'output token separator')

@dataclass
class Options:
  """Built from command line options.
  """
  # Patterns that define subset of files to be processed.
  # Not implemented.
  pattern: List[str] = field(default_factory=list)
  # output token separator.
  osep: str = ' '

# The various measures that are extracted from the protos.
@dataclass
class Results:
  """Hold results extracted from iwstats_pb2.Stats protos"""

  ave_abs_error: List[float] = field(default_factory=list)
  rms_error: List[float] = field(default_factory=list)
  bias: List[float] = field(default_factory=list)
  rsquared: List[float] = field(default_factory=list)
  qsquared: List[float] = field(default_factory=list)
  ae50: List[float] = field(default_factory=list)
  ae95: List[float] = field(default_factory=list)
  bsquared: List[float] = field(default_factory=list)

def get_responses(proto: iwstats_pb2.IWStats) -> Dict[str, Results]:
  """Return Dict of response name -> Results from `proto`.
  Args:
    proto: contains statistics on responses.
  Returns:
    dict from response name to a Results
  """
  result = {}
  for file in proto.stats:
    for res in file.stats:
      if not res.pred in result:
        result[res.pred] = Results()
        logging.info('New response %s', res.pred)

  return result

def write_summary_header(options: Options) -> None:
  """Write a header for `write_summary`.
  """
  tokens = ['pred', 'measure', 'nobs', 'min', 'max', 'ave', 'stddev']
  print(options.osep.join(tokens))


def write_summary(response: str, results: Results,
                  options: Options) -> None:
  """Write summary data for `response`.
  Args:
    response: name of the response
    results:  accumulated data for this response
    options:  options
  For each field name in `results` write summary stats.
  """
  field_names = [field.name for field in fields(Results)]
  for field_name in field_names:
    s = stats.describe(getattr(results, field_name))  # pylint: disable=C0103
    tokens = [response, field_name, str(s.nobs), f'{s.minmax[0]:.3f}', f'{s.minmax[1]:.3f}',
              f'{s.mean:.3f}', f'{math.sqrt(s.variance):.3f}']
    print(options.osep.join(tokens))


def do_append(src: iwstats_pb2.Stats, destination: Results):
  """Append the data in `src` to `destination`.
  Args:
    src:
    destination:
  """
  destination.ave_abs_error.append(src.ave_abs_error)
  destination.rms_error.append(src.rms_error)
  destination.bias.append(src.bias)
  destination.rsquared.append(src.rsquared)
  destination.qsquared.append(src.qsquared)
  destination.ae50.append(src.ae50)
  destination.ae95.append(src.ae95)
  destination.bsquared.append(src.bsquared)


def process_all_files(proto: iwstats_pb2.IWStats, responses: Dict[str, Results],
                      options: Options):
  """Process all data in `proto`.

  Args:
    proto:
    responses: mapping from prediction name to Results
    options:
  """
  for file in proto.stats:
    for result in file.stats:
      pred = result.pred
      if not pred in responses.keys():
        logging.fatal('Unknown response %s', pred)
      do_append(result, responses[pred])

  # Now that the data has been collected, write the results.
  write_summary_header(options)
  for response, data in responses.items():
    write_summary(response, data, options)

def get_patterns(patterns: List[str], files: Set[str]) -> List[List[str]]:
  """NOT IMPLEMENTED
  Return the list of files implied by each pattern.

  Args:
    patterns: Format strings describing a series of files. Must contain %d
    files: list of file names
  Returns:
    For each pattern, the list of files that match the enumerated pattern.
  """
  result: List[List[str]] = []
  # For each pattern, determine how many files are present.
  # Ignore 0 because maybe things start at 1 rather than zero.
  for pattern in patterns:
    files = []
    for i in count(start=0, step=1):
      fname = pattern%(i)
      if fname in files:
        files.append(fname)
        continue
      if i == 0:
        continue
      result.append(files)
      break

  return result

def summary(options: Options, proto: iwstats_pb2.IWStats):
  """proto is an IWStats::IWStats proto.
  Produce summary stats. If there are multiple types of predictions in
  there, do comparisons.
  """
  responses = get_responses(proto)

  if len(options.pattern) == 0:
    process_all_files(proto, responses, options)
    return

  logging.info('Functionality not implemented')
  return

  # Get the set of all file names in the proto
  files = [p.fname for p in proto.stats]
  # For each pattern work out how many files present.
  pattern_files = get_patterns(options.pattern, files)

  # Must be same number of values for each
  # lengths = set([len(p) for p in pattern_files])
  lengths = set(len(p) for p in pattern_files)
  if len(lengths) != 1:
    logging.fatal('Disparate counts %v', pattern_files)

  logging.info('Patterns contain %d files eacg', lengths.pop())

def ltgt_symbol(measure: str, means: List[float]) -> str:
  """Return < or > depending in which of `means` is better.
  Args:
  Returns:
  """
  if measure in ['ave_abs_error', 'rms_error', 'ae50', 'ae95']:
    if means[0] < means[1]:
      return '<'
    else:
      return '>'

  # Need to take abs value before comparing
  if measure == 'bias':
    if abs(means[0]) < abs(means[1]):
      return '<'
    else:
      return '>'

  # One of the cases where high values are better.
  if means[0] < means[1]:
    return '>'
  else:
    return '<'


def make_comparisons(results: List[Results],
                     fnames: List[str],
                     response: str, options: Options) -> None:
  """Write pairwise tests between the two items in `results`.
  Args:
    results: data extracted from protos
    fnames: list of file names that generates `results`
    response: name of the response
    options:
  """

  header = ['Measure', 'Response', f'{fnames[0]}_ave', f'{fnames[0]}_std',
            f'{fnames[1]}_ave', f'{fnames[1]}_std', 't', 'p', 'gtlt']
  print(options.osep.join(header))

  field_names = [field.name for field in fields(Results)]
  for field_name in field_names:
    values: List[List[float]] = []
    tokens = [field_name, response]
    means: List[float] = []
    for res in results:
      data = getattr(res, field_name)
      values.append(data)
      s = stats.describe(data)  # pylint: disable=C0103
      tokens.append(f'{s.mean:.3f}')
      means.append(s.mean)
      tokens.append(f'{math.sqrt(s.variance):.3f}')

    ttr = stats.ttest_rel(values[0], values[1])
    tokens.extend([f'{ttr.statistic:.4f}', f'{ttr.pvalue:.5f}'])
    tokens.append(ltgt_symbol(field_name, means))
    print(options.osep.join(tokens))

def comparison(protos: List[iwstats_pb2.IWStats],
               fnames: List[str],
               options: Options) -> None:
  """
  Args:
    protos:
    fnames:
    options:
  """
  results: List[Results] = []
  # All results must be for the same response
  common_pred = ''
  for proto in protos:
    res = Results()
    for file in proto.stats:
      for data in file.stats:
        if len(common_pred) == 0:
          common_pred = data.pred
        elif data.pred != common_pred:
          logging.fatal('Predicted name mismatch %s vs %s', common_pred, data.pred)

        do_append(data, res)

    results.append(res)

  make_comparisons(results, fnames, common_pred, options)

def main(argv):
  """Examine an IWStats::IWStats proto and provide summary and comparison info
    The binary proto must come in as the argument
  """
  options = Options()
  options.pattern = FLAGS.pattern
  options.osep = FLAGS.osep

  if len(argv) == 1:
    logging.fatal('Must specify IWStats::IWStats binary proto as argument')

  if len(argv) == 2:
    proto = iwstats_pb2.IWStats()
    with open(argv[1], 'rb') as reader:
      proto.ParseFromString(reader.read())

    summary(options, proto)
  elif len(argv) != 3:
    logging.fatal('Can only perform pairwise comparisons')
  else:
    protos = []
    fnames = []
    for fname in argv[1:]:
      fnames.append(fname)
      proto = iwstats_pb2.IWStats()
      with open(fname, 'rb') as reader:
        proto.ParseFromString(reader.read())
        protos.append(proto)

    comparison(protos, fnames, options)

if __name__ == '__main__':
  app.run(main)
