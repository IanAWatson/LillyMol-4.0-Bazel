# Python version of unique molecules

from unique_molecules import unique_molecules

from absl import app
from absl import flags
from absl import logging

FLAGS = flags.FLAGS

flags.DEFINE_boolean('verbose', False, 'verbose output')
flags.DEFINE_boolean('ignore_chiral', False, 'Ignore chirality')
flags.DEFINE_boolean('largest_frag', False, 'strip to largest fragment')
flags.DEFINE_boolean('ignore_cis_trans', False, 'ignore cis trans bonds')
flags.DEFINE_boolean('ignore_isotopes', False, 'ignore isotopic labels')
flags.DEFINE_boolean('standardize', False, 'apply chemical standardization (charges, tautomers, ..)')
flags.DEFINE_string('dupefile', '', 'write duplicates to <dupefile>')

def main(argv):
  if len(argv) == 1:
    logging.error("Must specify input file")
    return 1

  verbose = FLAGS.verbose

  um = unique_molecules()
  if FLAGS.ignore_chiral:
    um.set_exclude_chiral_info(1)
  if FLAGS.largest_frag:
    um.set_reduce_to_largest_fragment(1)
  if FLAGS.ignore_cis_trans:
    um.set_exclude_cis_trans_bonding_info(1)
  if FLAGS.ignore_isotopes:
    um.set_ignore_isotopes(1)
  if FLAGS.standardize:
    um.activate_standardization('all')

  if len(FLAGS.dupefile) > 0:
    dupefile = open(FLAGS.dupefile, 'w')

  with open(argv[1], 'r') as reader:
    lines = reader.read().splitlines()

  if verbose:
    logging.info("Read %d lines from %s", len(lines), argv[1])

  dupes = 0
  for line in lines:
    f = line.split()
    if not um.is_unique(f[0]):
      dupes += 1
      if dupefile is not None:
        print(line, file=dupefile)

  logging.info("Found %d duplicate structures", dupes)
    

if __name__ == '__main__':
  app.run(main)
