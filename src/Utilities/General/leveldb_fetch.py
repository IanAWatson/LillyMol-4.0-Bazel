# Fetch keys from a leveldb database

from absl import app
from absl import flags
from absl import logging

import plyvel

FLAGS = flags.FLAGS

flags.DEFINE_string('d', '', 'Database to open')
flags.DEFINE_multi_string('K', [], 'Key(s) to fetch')
flags.DEFINE_boolean('verbose', False, 'Verbose output')

def main(argv):
  """List the contents of a leveldb database
  """
  verbose = FLAGS.verbose

  if FLAGS.d == '':
    logging.error('Must specify database name via the -d option')
    return 1

  dbname = FLAGS.d

  db = plyvel.DB(dbname, create_if_missing=False)

  keys_looked_up = 0
  key_values_retrieved = 0
  for key in FLAGS.K:
    keys_looked_up += 1
    print(int(key).to_bytes(4, 'little'))
    value = db.get(int(key).to_bytes(4, 'little'))
    if value is None:
      pass
    else:
      key_values_retrieved += 1
      print(value.decode())

  if verbose:
    logging.info('Looked up %d found %d', keys_looked_up, key_values_retrieved)

if __name__ == '__main__':
  app.run(main)
