### Provide summary info from an IWStats::IWStats proto
###

from dataclasses import dataclass
from dataclasses import field, fields
from itertools import count
import math 
from os.path import exists
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
  pattern: List[str] = field(default_factory=list)
  osep: str = ' '

@dataclass
class Results:
  ave_abs_error: List[float] = field(default_factory=list)
  rms_error: List[float] = field(default_factory=list)
  rsquared: List[float] = field(default_factory=list)
  qsquared: List[float] = field(default_factory=list)
  ae50: List[float] = field(default_factory=list)
  ae95: List[float] = field(default_factory=list)
  bsquared: List[float] = field(default_factory=list)

def get_responses(proto: iwstats_pb2.IWStats) -> Dict[str, Results]:
  """Return Dict of response name -> Results from `proto`.
  Args:
    proto: contains statistics on responses.
  """
  result = {}
  for file in proto.stats:
    for r in file.stats:
      if not r.pred in result.keys():
        result[r.pred] = Results()
        logging.info('New response %s', r.pred)

  return result

def write_summary_header(options:Options) -> None:
  """Write a header for `write_summary`.
  """
  s = ['pred', 'measure', 'nobs', 'min', 'max', 'ave', 'stddev']
  print(options.osep.join(s))
  

def write_summary(response: str, results: Results,
                  options: Options) -> None:
  """
  Args:
    response: name of the response
    proto:    accumulated data for this response
    options:  options
  """
  n = len(results.ae50)
  field_names = [field.name for field in fields(Results)]
  for field_name in field_names:
    s = stats.describe(getattr(results, field_name))
    tokens = [response, field_name, str(s.nobs), f'{s.minmax[0]:.3f}', f'{s.minmax[1]:.3f}',
              f'{s.mean:.3f}', f'{math.sqrt(s.variance):.3f}']
    print(options.osep.join(tokens))
  

def do_append(src: iwstats_pb2.Stats, destination: Results):
  destination.ave_abs_error.append(src.ave_abs_error)
  destination.rms_error.append(src.rms_error)
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
    responses:
    options:
  """
  for file in proto.stats:
    for result in file.stats:
      pred = result.pred
      if not pred in responses.keys():
        logging.fatal('Unknown response %s', pred)
      do_append(result, responses[pred])

  write_summary_header(options)
  for k, v in responses.items():
    write_summary(k, v, options)

def get_results(files:List[str]) -> None:
  """
  Args:
    files:
  Returns:
  """
  result: list[Results] = []


def get_series(files:Set[str], options:Options) -> List[str]:
  """Discern series present in `files`.
  Args:
    files: Set of file names.

  Files is a set of file names that is assumed to consist of one or
    more sequences of files. For example 'abc1.dat' 'abc2.dat', 
    'def1.dat', 'def2.dat' would return ['abc', 'def']
    Or 'abc.1', 'abc.2' ...
  """
  if len(options.pattern) == 0:
    return get_series_all(files)
  return get_series_via_pattern(files, options.pattern)

def get_patterns(patterns:List[str], files:Set[str]) -> List[List[str]]:
  """Return the list of files implied by each pattern.

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
      if exists(fname):
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

  # Get the set of all file names in the proto
  files = [p.fname for p in proto.stats]
  # For each pattern work out how many files present.
  pattern_files = get_patterns(options.pattern, files)

  # Must be same number of values for each
  s = set([len(p) for p in pattern_files])
  if len(s) != 1:
    logging.fatal('Disparate counts %v', pattern_files)

  logging.info('Patterns contain %d files eacg', s.pop())

  results = []
  for files in pattern_files:
    results.append(get_results(proto, files))

def main(argv):
  """Examine an IWStats::IWStats proto and provide summary and comparison info
    The binary proto must come in as the argument
  """
  if len(argv) == 1:
    logging.fatal('Must specify IWStats::IWStats binary proto as argument')

  proto = iwstats_pb2.IWStats()
  with open(argv[1], 'rb') as reader:
    proto.ParseFromString(reader.read())
    
  options = Options()
  options.pattern = FLAGS.pattern
  options.osep = FLAGS.osep

  return summary(options, proto)

if __name__ == '__main__':
  app.run(main)
