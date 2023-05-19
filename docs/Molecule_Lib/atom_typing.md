# Atom Typing

## Background.

When generating a molecular fingerprint, the result is a combination of
two effects
* the shape of the feature being perceived
* the numeric description of the atoms involved.
Common choices for the numeric descriptions of the atoms include
atomic number and pharmacaphore type. But there are many other possibilities.

Experience shows that atom types beyond atomic number and pharmacaphore
can perform much better in svmfp models.

Unfortunately there is a very large number of possible atom types.

## Implementation.
Many LillyMol tools support a command line specified atom typing. Most commonly
this will be via the `-P` option, followed by a single token that specifies
what kind of atom typing to apply.

When development began, it was anticipated that developers would come up
with interesting atom types, which would be named. Some atom types are
specified that way.

But as more atom type information became available, it became obvious that
naming the combinations would not scale, so the idea of a User Specified Type, UST,
came into form.

Therefore the atom type specified by the `-P` option might be something like
```
-P sb
```
which means use Sybyl atom types. But if you wanted an atom type that consisted
of atomic number and whether or not the atom was aromatic,
```
-P UST:AZ
```
would be used.

## Named atom types
The following named atom types are recognised with the `-P` option.

In all cases, the description here is cursory, describing only the
atomic properties included in how that atom typing is computed. See
the source code to see how those are actually implemented.

### sb
Sybyl (now defunct molecular modeling software) defined an atom typing
scheme, and this atom type replicates that as best we could based on
public data available at the time.

### tt
An atom typing based on the atom typing described in one of the early papers
on Topological Torsions from Merck. This implementation combines the
number of connections to the atom, the number of pi electrons and the
atomic number.

### complex
This is one of the most successful atom typing schemes, frequently found
to be an excellent choice with svm fingerprint models. It combines
the hydrogen count, pi electrons, aromaticity and number of connections
on each atom. Note that it does not explicitly include the atomic
number, although due to valence rules presumably only certain combinations
can exist.

### bs
This was an experimental type, called 'basic'. It has seldom been one of
the best choices in QSAR models. The atomic property is defined as.
* If the atom is aromatic, type is aromatic.
* Saturated or unsaturated carbon
* Saturated atom
* Hydrogen count.
How these are combined can best be explained by looking at the 
source code in the function *assign_atom_types_basic*.

### expt
an experimental type, used for testing new ideas. Currently it includes
aromaticity, atomic number and pi electrons, combined in an unusual way.
Generally use this type for testing new ideas, it is not designed to
be permanent.

### none
all atoms have the same type. This can be useful for comparing
molecular skeletons.

### sf
An atom typing used in the synthetic feasibility tool. It combines the
atomic number, aromaticity, whether or not the atom is in a ring and
whether or not the atom is saturated. This was later superceeded by the
`sfx` type.

### sfx
this is the atom type currently used in the synthetic feasibility tools.
It consists of 
* atomic number
* ring membership
* number of connections 
* number of hydrogens
* formal charge.

### ch
Carbon or heteroatom is the only differentiation.

### cc
A variant on `complex` that adds influence of formal charges.

### pp
A pharmacaphore atom typing.
If this is used, there are two ways by which the charge
assigner and donor/acceptor queries can be specified.
Either set the shell variable
*IW_PHARMACAPHORE*, or
```
-P PP:/path/to/directory
```
to specify a directory in which these files are found.

### za
Atomic number and aromaticity. THis should not exist as a named
type, the same effect is available with
```
-P UST:AZ
```

### zp
Atomic numbers encoded as prime numbers. Obsolete.

## UST
User Specified Types provide the most flexibility in constructing
atom types. More than 20 atomic properties can be combined.
Some combinations set the atomic invariant, while others alter
a value that has already been computed.

### A (alters)
Aromaticity.

### B (alters)
Unsaturation, but aromatic *not* included

### C (alters)
Number of connections

### E (sets)
Carbon or heteroatom.

### F (alters)
Only ring atoms are processed. Atoms involved in ring fusions are differentiated
from other ring atoms.

### H (alters)
Hydrogen count, implicit + explicit.

### I (alters)
The isotope on the atom.

### K (alters)
The combined atomic number of the atom and its neighbours.

### M (alters)
Starts with Y, and then all aromatic atoms then become identical.

### N (sets)
No type, all atoms identical. Note that if the fingerprint differentiates
bond types, this will not result in benzene and cyclohexane becoming
identical.

### O (alters)
Formal charge

### P (alters)
Pi electrons. The value added will depend on the number of pi electrons.

### Q (alters)
Pi electrons. The value added is the same regardless of how many pi electrons.

### R (alters)
The ring bond count for the atom.

### S (alters)
The size of the smallest ring containing the atom.

### T (sets)
A complex, but very successful atom type. Starts with the atomic number. Certain
Sulphur types are considered equivalent to Oxygen, and certain possibly
tautomeric Nitrogen atoms are combined.

### U (alters)
Unsaturation - includes aromatic.

### X (alters)
A centrality concept. Examines all atoms in the molecule and computes an
average distance. Expensive, unproven, maybe not a good idea.

### Y (sets)
Atomic number, but heavy halogens (Cl, Br and I) are equivalent. I have never
seen a model where this is worse than full atomic number (Z), and have seen
plenty of models that are improved by this.

### Z (sets)
Atomic number

## Typical Usage
Some common atom type combinations are

* -P c
* -P tt
* -P sfx
* -P UST:ACHPY
* -P UST:Y
* -P UST:ARY
* -P UST:APT
* -P UST:ACY
* -P UST:HRY

By convention atom types are combined alphabetically in order to avoid
any unfortunate combinations.

But as can be imagined, there are a very large number of combinations
possible. What might work best in any particular context is unpredictable.
