#include <iostream>

#define IWARCHIVE_IMPLEMENTATION
#define IWARCHIVE_OP_IMPLEMENTATION
#include "iwarchive.h"

template class iwarchive<double>;

template ostream & operator << (ostream &, const iwarchive<double> &);
