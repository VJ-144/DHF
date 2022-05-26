#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

// Computation of an integral

double g (double *r, size_t dim, void *params) {
  (void)(dim); /* avoid unused parameter warnings */
  (void)(params);
  return r[0]*pow(r[1],2);

}

void display_results ( double result, double error) {
  // printf ("%s ==================\n", title);
  printf ("result = % .6f\n", result);
  printf ("sigma  = % .6f\n", error);
}

int main (void) {

  double res, err;
  double xl[2] = { 0, 0 };
  double xu[2] = { 2, 1 };

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G = { &g, 2, 0 };
  size_t calls = 500000;
  gsl_rng_env_setup ();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  gsl_monte_plain_state *s = gsl_monte_plain_alloc (2);
  gsl_monte_plain_integrate (&G, xl, xu, 2, calls, r, s, &res, &err);
  gsl_monte_plain_free (s);

  display_results (res, err);

  return 0;
}


  // g++ MonCarlo_IntTest.cpp -lgsl -o Int
