# There have been multiple invocations of ring_extraction, each generating
# text_format protos. Collate those into a unified set of replacement rings.

import re, os, sys

from google.protobuf import text_format

from absl import app
from absl import flags
from absl import logging

from replacement_ring_pb2 import ReplacementRing

FLAGS = flags.FLAGS

flags.DEFINE_string('stem', 'ALL', 'file name step for results')


class Ring:
  def __init__(self, smi, usmi, smt, ex, c):
    self.smiles = smi
    self.usmi = usmi
    self.smarts = smt
    self.exemplar = ex
    self.count = c


def get_rings(fname, rings):
  print(f'Opening {fname}')
  with open(fname, 'r') as reader:
    for line in reader:
      proto = text_format.Parse(line, ReplacementRing())
      if proto.usmi in rings:
        rings[proto.usmi].count += proto.n
      else:
        rings[proto.usmi] = Ring(proto.smi, proto.usmi, proto.smt, proto.id, proto.n)

# Scan the directory for `stem` and accumulate all the rings.
# Returns a hash from ring
def get_rings_from_stem(stem, rings):
  dirname = os.path.dirname(stem)
  fname  = os.path.basename(stem)
  if dirname == '':
    dirname = '.'
  rx = re.compile(f"/{fname}_(\\S+)\\.smi")
  print(f'dirname {dirname} fname {fname} fname {fname}')
  print(f'Pattern {rx.pattern}')
  for fname in os.listdir(dirname):
    fname = os.path.join(dirname, fname)
    if os.path.isdir(fname):
      continue
    # print(f'Examining {fname}')
    m = rx.search(fname)
    if m is None:
      continue

    logging.info('%s', fname)
    rtype = m[1]
    if not rtype in rings:
      rings[rtype] = {}

    get_rings(fname, rings[rtype])

def main(argv):
  """Aggregate multiple replacement ring datasets.
  """

  output_stem = FLAGS.stem

  if len(argv) == 1:
    logging.info('Must specify file name stems for one or more existing result sets')
    sys.exit(1)

  rings = {}
  for stem in argv[1:]:
    get_rings_from_stem(stem, rings)

  for (label, data) in rings.items():
    fname = f'{output_stem}_{label}.smi'
    print(f'Writing {fname}')
    with open(fname, 'w') as writer:
      for (usmi, ringdata) in data.items():
        proto = ReplacementRing()
        proto.usmi = usmi
        proto.smi = ringdata.smiles
        proto.smt = ringdata.smarts
        proto.id = ringdata.exemplar
        proto.label = label
        proto.n = ringdata.count
        print(text_format.MessageToString(proto, as_one_line=True), file=writer)


if __name__ == '__main__':
  app.run(main)
