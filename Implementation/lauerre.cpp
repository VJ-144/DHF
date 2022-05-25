#include <iostream>
#include <iomanip>
#include <fstream>

#include <armadillo>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf.h>
#include <cmath>
#include <boost/range/irange.hpp>

using namespace std;
using namespace arma;
using namespace boost;

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

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Function Declerations
double laguerre_func_r(double x, double n, double l, double zeta);
double norm_r(double n1, double l1, double n2, double l2);
double coef_rinv(double n, double l);

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double laguerre_func_r(double x, double n, double l, double zeta){
    // Laguerre function, see [A. E. McCoy and M. A. Caprio, J. Math. Phys. 57, (2016).] for details
    double eta = 2.0 * x / zeta;
    return sqrt(2 * tgamma(n+1) / (zeta * tgamma(n+2*l+3) ) ) * 2.0 * pow(eta, l) * exp(-0.5 * eta) * gsl_sf_laguerre_n( n, 2*l+2, eta) / zeta;
}

double norm_r(double n1, double l1, double n2, double l2){
    if (l1 != l2) {return 0.0;}
    
    double result, error;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);

    auto lambda = [&](double x) {return (laguerre_func_r(x, n1, l1, zeta) * laguerre_func_r(x, n2, l2, zeta) * pow(x,2)); };
    gsl_function_pp<decltype(lambda)> Fp(lambda);
    gsl_function* F = static_cast<gsl_function*>(&Fp);

    gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 10000, w, &result, &error);
    gsl_integration_workspace_free(w);

    return result;
}

double coef_rinv(double n, double l){

    double result, error;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);

    auto lambda = [&](double x) {return (laguerre_func_r(x, n, l, zeta)); };
    gsl_function_pp<decltype(lambda)> Fp(lambda);
    gsl_function* F = static_cast<gsl_function*>(&Fp);

    gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 10000, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result;
}

int main(){

    int nmax = 60;
    map <int, double> cnlist;
    for (auto n : irange(nmax)) {
        cnlist[n] = coef_rinv(n, 0);
        // cout << cnlist[n] << endl;
    }

   

    vec x = arma::linspace(0,2,20);

     vec f;
    f.resize(x.size());

    ofstream output_file;
    output_file.open("lauerre_data.txt", ios_base::out);
    output_file << fixed << setw(3) << "x" << setw(17) << "1/x"  << setw(16) << "f" << endl;

    for (int i=0; i<x.size(); i++) {
        f[i] = cnlist[0] * laguerre_func_r(x[i], 0, 0, zeta);

        // for (int n : irange(1,nmax)) {

        //     f += cnlist[n] * laguerre_func_r(x[i], n, 0, zeta);
        //     // cout << f << endl;
        // }

        // output_file  << fixed << setprecision(4) << setw(3) << x[i] << setw(15) << 1/x[i] << setw(17) << f[i] << endl;

    }

    for (int i=0; i<x.size(); i++) {
        for (auto n : irange(1,nmax)) {
            f += cnlist[n] * laguerre_func_r(x[i], 0, 0, zeta);

        }
        output_file  << fixed << setprecision(4) << setw(3) << x[i] << setw(15) << 1/x[i] << setw(17) << f[i] << endl;

    }





    


    // for (int n : irange(1,nmax)) {
    //     // cout << n << endl;
    //     f += cnlist[n] * laguerre_func_r(x[n], n, 0, zeta);
    //     // cout << f[n] << endl;
    // }



    output_file.close();
    return 0;
}

// Run Command
// g++ lauerre.cpp -larmadillo -lgsl -lwignerSymbols -o lau