# Generate various rdkit fingerprint types. Designed to be part
# of a gfp_make like utility

from dataclasses import dataclass, field
import re
import sys
from typing import List

from absl import app
from absl import flags
from absl import logging

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions

FLAGS = flags.FLAGS

flags.DEFINE_list('fp', ['rdkit'],  'which fingerprint to use.')
flags.DEFINE_boolean('filter', False, "Operate as a TDT filter")
flags.DEFINE_boolean("verbose", False, "verbose output")

# Information gathered from the command line, passed to functions.
@dataclass
class Options:
  # The fingerprints being generated
  fptype: List

# Several functions that generate a specific fingerprint type and
# write to stdout.
def rdkit_fingerprint(mol):
    print(f'ASRDKLIN<{Chem.RDKFingerprint(mol).ToBitString()}>')

def ec_fingerprint(mol):
    print(f'ASRDKEC<{AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=2048).ToBitString()}>')

def ap_fingerprint(mol):
    print(f'ASRDKAP<{Pairs.GetAtomPairFingerprintAsBitVect(mol).ToBitString()}>')

def tt_fingerprint(mol):
    print(f'ASRDKTT<{Torsions.GetTopologicalTorsionFingerprintAsIntVect(mol).ToBitString()}>')

def mk_fingerprint(mol):
    print(f'ASCFPRDKMK<{MACCSkeys.GenMACCSKeys(mol).ToBitString()}>')


# Counted fingerprints are different.
# The post-processing is done as text, which is not great.
def ec_counted_fingerprint(mol):
    nzero = AllChem.GetMorganFingerprint(mol, 3).GetNonzeroElements()
    s = str(nzero)
    s = re.sub(r': ([0-9][0-9]*)[^0-9]', ':\\1', s)
    print(f'ASCNCRDEC<{s[1:]}>')


def fingerprint_and_write(fp: str, mol):
  """generate fingerprint `fp` for `mol` and write
    Args:
      fp: a fingerprint type.
      mol: a molecule.
  """
  if fp == "rdkit":
    return rdkit_fingerprint(mol)
  elif fp == "ec":
    return ec_fingerprint(mol)
  elif fp == "ap":
    return ap_fingerprint(mol)
  elif fp == "tt":
    return tt_fingerprint(mol)
  elif fp == "mk":
    return mk_fingerprint(mol)
  elif fp == "ecc":
    return ec_counted_fingerprint(mol);

  logging.fatal('Unrecognised fingerprint form %s', fp)


def write_all_fingerprints(options: Options, mol):
  """Write all fingerprints for `mol`.

  Args:
    options: from the command line
    mol: a molecule
  """

  for fp in options.fptype:
    fingerprint_and_write(fp, mol)

def rdkit_fingerprint_smiles(options: Options, argv: List[str]):
  """Generate fingerprints for the smiles in `argv`
  
    Args:
      options:
      argv:
  """
  for fname in argv[1:]:
    with open(fname, 'r') as reader:
      for line in reader:
         f = line.rstrip().split(' ')
         mol = Chem.MolFromSmiles(f[0])
         if mol is None:
           logging.warning('Cannot read smiles %s', f[0])
           continue

         print(f'$SMI<{f[0]}>')
         print(f'PCN<{f[1]}>')
         write_all_fingerprints(options, mol)
         print("|")

def rdkit_fingerprint_filter(options: Options, argv:List[str]):
  """Read from a TDT file and insert fingerprints"
    Args:
      options:
      argv: ignored, we only recognise '-' as the input.
  """
  for line in sys.stdin:
    print(line.rstrip())
    regex_match = re.match('^\$SMI<(\S+)>', line)
    if regex_match is None:
      continue

    # print(f'regex_match {regex_match} {regex_match.group(1)}')
    mol = Chem.MolFromSmiles(regex_match.group(1))
    if mol is None:
      logging.warning('Cannot interpret smiles %s', regex_match.group(1))
      continue;

    write_all_fingerprints(options, mol)

def rdkit_fingerprint(argv):
  """Generates RDKit fingerprints"""
  if len(argv) == 1:
    logging.error("Must specify input smiles file")

  verbose = FLAGS.verbose
  options = Options(fptype=FLAGS.fp)
  if FLAGS.filter:
    rdkit_fingerprint_filter(options, argv)
  else:
    rdkit_fingerprint_smiles(options, argv)

if __name__ == '__main__':
  app.run(rdkit_fingerprint)
