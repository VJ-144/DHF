#include <iostream>
#include <gsl/gsl_integration.h>

using namespace std;

// Simple RAII wrapper 
class IntegrationWorkspace {
  gsl_integration_workspace * wsp;

  public:
  IntegrationWorkspace(const size_t n=1000):
    wsp(gsl_integration_workspace_alloc(n)) {}
  ~IntegrationWorkspace() { gsl_integration_workspace_free(wsp); }

  operator gsl_integration_workspace*() { return wsp; }
};

// Build gsl_function from lambda
template <typename F>
class gsl_function_pp: public gsl_function {
  const F func;
  static double invoke(double x, void *params) {
    return static_cast<gsl_function_pp*>(params)->func(x);
  }
  public:
  gsl_function_pp(const F& f) : func(f) {
    function = &gsl_function_pp::invoke; //inherited from gsl_function
    params   = this;                     //inherited from gsl_function
  }
  operator gsl_function*(){return this;}
};

// Helper function for template construction
template <typename F>
gsl_function_pp<F> make_gsl_function(const F& func) {
  return gsl_function_pp<F>(func);
}

int main() {
    double epsabs = 1e-8;
    double epsrel = 1e-8;
    size_t limit = 100;
    double result, abserr, inner_result, inner_abserr;

    IntegrationWorkspace wsp1(limit);
    IntegrationWorkspace wsp2(limit);

    auto outer = make_gsl_function( [&](double x) {
    auto inner = make_gsl_function( [&](double y) {return exp(-x*x-y*y);} );
    gsl_integration_qagi(inner, epsabs, epsrel, limit, wsp1,
                         &inner_result, &inner_abserr);
    return inner_result;
    } );
    gsl_integration_qagi(outer, epsabs, epsrel, limit, wsp2, &result, &abserr);

    std::cout << result << std::endl;
}

// g++ doubleIntegrationTest.cpp -lgsl -o Int
