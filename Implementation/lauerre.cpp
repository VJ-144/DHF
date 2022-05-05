#include <iostream>
#include <armadillo>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf.h>
#include <cmath>

using namespace std;
using namespace arma;

const int zeta = 1;

// Wrapper to convert lambda with capture to gsl function for integration
template< typename F >  class gsl_function_pp : public gsl_function {
public:
    gsl_function_pp(const F& func) : _func(func) {
        function = &gsl_function_pp::invoke;
        params = this;
    }
private:
    const F& _func;
    static double invoke(double x, void* params) {
        return static_cast<gsl_function_pp*>(params)->_func(x);
    }
};

// Function Declerations
double laguerre_func_r(double x, double n, double l, double zeta);
double norm_r(double n1, double l1, double n2, double l2);
double coef_rinv(double n, double l);

double laguerre_func_r(double x, double n, double l, double zeta){
    // Laguerre function, see [A. E. McCoy and M. A. Caprio, J. Math. Phys. 57, (2016).] for details
    double eta = 2.0 * x / zeta;
    return sqrt(2 * tgamma(n+1) / (zeta * tgamma(n+2*l+3) ) ) * 2.0 *pow(eta, l) * exp(-0.5 * eta) * gsl_sf_laguerre_n( n, 2*l+2, eta) / zeta;
}

double norm_r(double n1, double l1, double n2, double l2){
    if (l1 != l2) {return 0.0;}
    
    double result, error;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);

    auto lambda = [&](double x) {return (laguerre_func_r(x, n1, l1, zeta) * laguerre_func_r(x, n2, l2, zeta) * pow(x,x)); };
    gsl_function_pp<decltype(lambda)> Fp(lambda);
    gsl_function* F = static_cast<gsl_function*>(&Fp);

    gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 5000, w, &result, &error);
    gsl_integration_workspace_free(w);

    return result;
}

double coef_rinv(double n, double l){

    double result, error;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);

    auto lambda = [&](double x) {return (laguerre_func_r(x, n, l, zeta)); };
    gsl_function_pp<decltype(lambda)> Fp(lambda);
    gsl_function* F = static_cast<gsl_function*>(&Fp);

    gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 5000, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result;
}

int main(){
    int nmax = 60;
    vector<double> cnlist(nmax);
    // for testing use vec
    // vec cnlist(nmax);
    for (int n = 0; n < nmax; n++){
        cnlist[n] = coef_rinv(n, 0);
    }
    // cnlist.print();
    vec x = linspace(0,2,0.1);
    vector<double> f(x.size());

    for (int i = 1; i < nmax; i++){
        f[i] += cnlist[i] * laguerre_func_r(x[i], i, 0, zeta);
        cout << f[i] << endl;
    }

    return(0);
}
