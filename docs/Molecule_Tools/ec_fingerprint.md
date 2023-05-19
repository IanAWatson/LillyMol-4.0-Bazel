# ec_fingerprint

This tool generates extended connectivity fingerprints. The radius of
the fingerprint can be controlled, as can the atom typing used.

Output can be in one of several possible forms.

## Usage
The usage message is
```
  -r <rad>      minimum shell radius (def 0)
  -R <rad>      maximum shell radius (def 3)
  -P ...        atom type specification
  -J <tag>      tag for fingerprints
  -f            function as a TDT filter
  -D ...        specify what operation is to be performed
  -D file=<fname>  For those actions that require a file.
  -D fp=<Tag>      Generate fingerprints with a given tag
  -D explain       Look for bits in <fname> and provide explanations
  -D collision     Look for bit collisions and write
  -D pgen          Gather precedent data and write
  -D puse          Read previously generated precedent data and use
  -D writeall      Write all bits generated to <fname>
  -D filter        Write smiles of molecules that do NOT contain bits in <fname>
  -D args=...      Pass arbitrary text directives to the class doing the work
  -n            suppress fingerprint output (useful with the some -D options)
  -s            gather statistics on molecules processed
  -p            precise fingerprints - includes attachment point in inner shell
  -m            bit formation is multiplicative. shells not differentiated
  -q <query>    query specification of atoms to be processed
  -w <nbits>    generate fixed width binary fingerprints, `nbits` bits
  -l            reduce to largest fragment
  -i <type>     input specification
  -g ...        chemical standardisation options
  -E ...        standard element specifications
  -A ...        standard aromaticity specifications
  -v            verbose output
```
The radius of the fingerprint defaults to 3, but can be changed.

Atom typing is specified via the -P option.

The -J option controls output and is described below.

The -f option is used when this tool is part of a fingerprint
generation pipeline.

The -D option describes a variety of more complex functioning, described
below (TODO:ianwatson).

The -n option suppresses output, but some kind of accumulation of data
must be underway.

The -s option can be useful when characterising collections.

The -p option is important. By default, fingerprints are generated with
each shell independent of the inner shell - but with the -p option,
how the outer shell is connected to the inner shell is considered.
Generally this produces better performing fingerprints.

The -q option allows restricting the atoms being fingerprinted to a
subset of the atoms in the molecule. The usual syntax for specifying
a substructure query applies here. For example to use a smarts
```
-q 'SMARTS:c1ccccc1'
```
The -w option can be used to generate a fixed width array. Obsolete, use
the -J option instead.

The other options are all standard LillyMol options on all tools.

## The -J option.

This is complex, because the tool can generate several different output types.

The following directives are recognised

### -J array
Output is a fixed width tabular file.

### -J nbits=nb
When generating tabular output, generate <nb> columns. The smaller the
number, the smaller the output, but the greater the number of bit collisions.

### -J sep=,
Specify the output token separator when generating a fixed width tabular file. Default
is space. Certain abbreviations are recognised `-J sep=tab` for example.

### -J name=nn
A prefix for the bit number column headers when writing tabular output.

### -J fixed
Generate a fingerprint as a fixed TDT form.

### -J sparse
Generate a fingerprint as s sparse, non colliding TDT form.

## Examples
Here are some example invocations

```
ec_fingerprint -i ICTE -i smi -J array -J nbits=2048 -P UST:APT file.smi > file.txt
```

The `-i ICTE` option combination tells it to ignore connection table errors. In any
external set of molecules, there is overwhelming probability that some will be
unable to be consumed. This is usually due to aromaticity differences between
different software packages.
