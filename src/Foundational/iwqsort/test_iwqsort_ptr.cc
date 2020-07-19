#include <stdlib.h>

/*
  Tests the pointer version
*/

#include "cmdline.h"
#include "iwrandom.h"
#include "ptr.h"

#include "foo_.h"

static int verbose = 0;

static void
usage (int rc)
{
  cerr << " -s <n>        size of arrays to sort\n";

  exit (rc);
}

template <typename F>
int
check_sorted (const F * const * f,
              int s)
{
//cerr << "Checking " << s << " values, between " << f[0]->zvalue () << " and " << f[s - 1]->zvalue() << endl;
  int rc = 0;

  for (int i = 1; i < s; i++)
  {
 // cerr << "i = " << i << " value " << f[i]->zvalue () << endl;
    int c = f[i]->iwqsortcompare (*(f[i - 1]));
    if (c < 0)
    {
      cerr << "Sorted items out of order, i = " << i << " " << f[i - 1]->zvalue () << " followed by " << f[i]->zvalue () << endl;
      rc = 0;
    }
  }

  return rc;
}

template <typename F, typename T>
int
test_iwqsort (F ** f, 
              int s,
              int iterations,
              T * initial_values,
              void (&new_values) (T *, int))
{
  cerr << "Testing vectors of size " << s << endl;
  for (int i = 0; i < s; i++)
  {
    f[i]->set_value (static_cast<T> (8));
  }

  if (verbose)
    cerr << "Constant values\n";

  iwqsort_ptr (f, s);

  check_sorted (f, s);

  for (int i = 0; i < s; i++)
  {
    f[i]->set_value (static_cast<T> (i));
  }

  if (verbose)
    cerr << "Increasing values\n";

  iwqsort_ptr (f, s);

  check_sorted (f, s);

  for (int i = 0; i < s; i++)
  {
    f[i]->set_value (static_cast<T> (s - i));
  }

  if (verbose)
    cerr << "Decreasing values\n";

  iwqsort_ptr (f, s);

  check_sorted (f, s);

  if (verbose)
    cerr << "Begin " << iterations << " iterations\n";

  for (int i = 0; i < iterations; i++)
  {
    new_values (initial_values, s);
    for (int j = 0; j < s; j++)
    {
      f[j]->set_value (initial_values[j]);
    }

    iwqsort_ptr (f, s);

    check_sorted (f, s);
  }

  return 1;
}

#ifdef __GNUG__

typedef Foo<int>   Foo_int;
typedef Foo<float> Foo_float;

template int test_iwqsort (Foo_int **,   int, int, int *,   void (&) (int *,   int));
template int test_iwqsort (Foo_float **, int, int, float *, void (&) (float *, int));

template int check_sorted (const Foo<int> * const *, int);
template int check_sorted (const Foo<float> * const *, int);

//template void iwqsort_ptr (Foo_int **,   int, int (Foo_int::*)   (Foo_int *));
//template void iwqsort_ptr (Foo_float **, int, int (Foo_float::*) (Foo_float *));

template void iwqsort_ptr (Foo_int **,   int);
template void iwqsort_ptr (Foo_float **, int);

#endif

static void
random_numbers_int (int * f, int s)
{
  for (int i = 0; i < s; i++)
  {
    f[i] = intbtwij (0, 10);
  }

  return;
}

static int
test_iwqsort_int (int s, int iterations)
{
  Foo<int> ** f = new Foo<int> *[s];

  for (int i = 0; i < s; i++)
  {
    f[i] = new Foo<int>;
  }

  int * initial_values = new int[s];

  int rc = test_iwqsort (f, s, iterations, initial_values, random_numbers_int);

  for (int i = 0; i < s; i++)
  {
    delete f[i];
  }

  delete [] f;
  delete initial_values;

  return rc;
}


static void
random_numbers_float (float * f, int s)
{
  for (int i = 0; i < s; i++)
  {
    f[i] = static_cast<float> (iwrandom ());

    if (i < 10)
      cerr << f[i] << endl;
  }

  return;
}

static int
test_iwqsort_float (int s, int iterations)
{
  Foo<float> ** f = new Foo<float> *[s];

  for (int i = 0; i < s; i++)
  {
    f[i] = new Foo<float>;
  }

  float * initial_values = new float[s];

  int rc = test_iwqsort (f, s, iterations, initial_values, random_numbers_float);
  
  for (int i = 0; i < s; i++)
  {
    delete f[i];
  }

  delete [] f;
  delete [] initial_values;

  return rc;
}

static int
test_iwqsort (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vs:n:if");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "unrecognised_options_encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  iw_random_seed ();

  int s;
  if (cl.option_present ('s'))
  {
    if (! cl.value ('s', s) || s < 2)
    {
      cerr << "Must sort at least 2 items\n";
      usage (3);
    }
  }
  else
  {
    s = intbtwij (2, 1000);
  }

  if (verbose)
    cerr << "Will test vectors of length " << s << " items\n";

  int iterations = 1;
  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', iterations) || iterations < 1)
    {
      cerr << "The number of iterations must be a whole number\n";
      usage (5);
    }

    if (verbose)
      cerr << "Will perform " << iterations << " iterations\n";
  }

  if (cl.option_present ('f') && cl.option_present ('i'))
  {
    cerr << "Can do either floating point or integer tests\n";
    usage (6);
  }

  if (cl.option_present ('f'))
  {
    test_iwqsort_float (s, iterations);
    if (s > 2)
      test_iwqsort_float(s-1, iterations);
  }
  else
  {
    test_iwqsort_int (s, iterations);
    if (s > 2)
     test_iwqsort_int(s-1, iterations);
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = test_iwqsort (argc, argv);

  return rc;
}
