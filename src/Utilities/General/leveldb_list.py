# List the contents of a leveldb database

from absl import app
from absl import flags
from absl import logging

import plyvel

FLAGS = flags.FLAGS

flags.DEFINE_boolean('verbose', False, 'Verbose output')

def main(argv):
  """List the contents of a leveldb database
  """
  if len(argv) == 1:
    logging.error("Must specify name of database")

  verbose = FLAGS.verbose

  db = plyvel.DB(argv[1], create_if_missing=False)
  for key, value in db:
    k = int.from_bytes(key, "little", signed=False)
    if k == 1021:
      print('Got 1021')
    continue
    print(f'{key} {len(key)} {k}')
    print(f"Key {k} value {value.decode()}")


if __name__ == '__main__':
  app.run(main)
