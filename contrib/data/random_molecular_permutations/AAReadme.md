# Data files for random_molecular_permutations

This directory contains files that can be used with the -L option to
`random_molecular_permutations`. 

## single.smi
Some small fragments that are attached via the first atom. No attempt at anything
comprehensive, just some things that occurred to me.

## double.smi
Atoms that can be attached via a double bond to most atoms with 2 Hydrogen atoms
available to be displaced.

## aromR.smi
Scan through all of Chembl and identify all aromatic rings. This file contains
all aromatic rings found in Chembl, 10 or more instances. The tool will attach
them via any atom with a Hydrogen available. Better would be to preserve the
substitution points and only allow substitution there. That will hopefully make
its way to this repo.

## Invocation

Here is an invocation that uses these files, assuming that the variable FRAGS is set to
this directory - the one containing single.smi...
```
random_molecular_permutations -M nowarnvalence -y 1 -Y 5 -c 14 -C 50 -v -L single=$FRAGS/single.smi -L double=$FRAGS/double.smi -L singleR=$FRAGS/singleR.smi -L aromR=$FRAGS/aromR.smi file.smi
```
