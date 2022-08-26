# Unique Molecules
A common task in Cheminformatics is the identification of unique or
duplicate structures. The unique smiles is commonly used for this.

`unique_molecules` reads a stream of molecules, and if any have been
seen before, that molecule is discarded.

The usage message is
```
  -l             strip to largest fragment
  -a             compare as tautomers - skeleton and Hcount
  -c             exclude chiral info - optical isomers will be duplicates
  -z             exclude cis/trans bonding information
  -I             ignore isotopic labels
  -f             function as filter (TDT input)
  -G <tag>       identifier tag when working as a filter
  -p <fname>     specify previously collected molecules
  -s <size>      specify primary hash size (default 1000)
  -S <name>      specify output file name stem
  -D <name>      write duplicate structures to <name>
  -R <rxn>       perform reaction(s) on molecules before comparing
  -T             discard molecular changes after comparison
  -r <number>    report progress every <number> molecules
  -e             report all molecules together with counts
  -j             items are the same only if both structure and name match
  -t E1=E2       element transformations, enter '-t help' for details
  -i <type>      specify input type
  -o <type>      specify output type(s)
  -A <qualifier> Aromaticity, enter "-A help" for options
  -K ...         standard smiles options, enter '-K help' for info
  -g <qualifier> chemical standardisations, enter "-g help" for usage
  -v             verbose output
```

There are a great many ways by which two molecules can be considered identical.
`unique_molecules` provides a wide variaty of structural modifiers that can
be applied to molecules before their unique smiles are generated.

* -l strip to the largest fragment
* -c discard chirality
* -z discard cis/trans bonding information
* -I discard isotopic labels.
* -R <rxn> apply a transformation reaction to the molecules.
* -t E1=E2 transform all elements E1 to E2 (suggest -T I=Cl -T Br=Cl)

Once these transformations are applied, the unique smiles is generate
and compared to what has been encountered previously. If this is
the first instance of the smiles, that molecule is written to the output
strea, otherwise it is classified as a duplicate, and if specified,
written to the -D file for duplicates.

By default, the changed molecule will be written, but the original form
can be written if the -T option is used. That should probably be the
default.

If there is a previous set of molecules that should be considered,
which will establish the unique smiles hash, without writing anything.

by default, smiles are accumulated in a C++ set, which has only keys. If
you wish to get a summary of smiles and counts, use the -e option,
which will use a map for the unique smiles, and counts can be accumulated.
