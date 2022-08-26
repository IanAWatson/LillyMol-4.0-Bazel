# Msort

msort is a tool that sorts a file of molecules.

Sorting can be done on the basis of a variety of molecular attributes. The most
common atrribute is typically the number of atoms. A typical invocation 
might be
```
msort -l -k natoms -k amw file.smi > file.sorted.smi
```
This sorts `file.smi` by atom count, within the largest frgament since
the `-l` option was specified. A secondary sort criterion is molecular weight,
so this will place the Fluroine version of molecule before the Bromine.

When `msort` development commenced, sorting options were specified by
letter options. As more sort criteria were added, a more flexible way was
needed, and so the -k option was used. That should be used going forward.

## Attributes
The  following sorting attributes are implemented.

### natoms
The number of atoms in the connection table. If there are explicit Hydrogen atoms
those are also counted.

### nrings
The number of rings in the molecule.

### amw
The average molecular weight of the whole molecule.

###
