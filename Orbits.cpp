#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf.h>
#include <armadillo>
#include <math.h>
#include <cmath>
#include <string>
#include <regex>
#include<tuple>
#include <string_view>

#include "Orbits.h"

using namespace std;
using namespace arma;

    Orbit::Orbit(string radial_function_type1) {
        // Constructor 
        n = -1; // nodal
        l = -1; // orbital angular momentum
        j = 0;  // total angular momentum
        k = 0;  // angular
        e = 0;  // upper / lower 1 / -1
        radial_function_type = radial_function_type1;
    }

    Orbits::Orbits(int emax1, int lmax1, string radial_function_type1, bool verbose1, bool relativistic1) {
        
        norbs = -1;
        emax = -1;
        lmax = -1;
        orbits = {};
        nlje_idx = {};
        _labels_orbital_angular_momentum = { 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'o', 'q', 'r', 't', 'u', 'v', 'w', 'x', 'y', 'z' };
        
        // Default Parameters Override 
        verbose = verbose1;
        relativistic = relativistic1;
        radial_function_type = radial_function_type1;

        set_orbits(emax1, lmax1);

    }

    Orbits::Orbits(string radial_function_type1, bool verbose1, bool relativistic1) {
        
        norbs = -1;
        emax = -1;
        lmax = -1;
        orbits = {};
        nlje_idx = {};
        _labels_orbital_angular_momentum = { 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'o', 'q', 'r', 't', 'u', 'v', 'w', 'x', 'y', 'z' };
        
        // Default Parameters Override 
        verbose = verbose1;
        relativistic = relativistic1;
        radial_function_type = radial_function_type1;

        
    }


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//                      End of Constructors 


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


    vector<int> Orbit::get_nljke() {  
        vector<int> result = { n, l, j, k, e };
        return result;
    }

    int Orbit::_lj_to_k(int l, int j) {
        return floor( - (4.0 + j * (j + 2.0) - 4.0 * l * (l + 1.0) - 3.0) / 4.0);
    }


    vector<int> Orbit::_k_to_lj(int k) {
        vector<int> result;
        if (k > 0) {
            result = { k, 2 * k - 1 };
            return result;
        }
        else if (k < 0) {
            result = { -k - 1, -2 * k - 1 };
            return result;
        }
        else {return {0};}
    }


    void Orbit::set_orbit_nke(vector<int> nke) {
        n = nke[0];
        k = nke[1];
        e = nke[2];
        l = _k_to_lj(k)[0];
        j = _k_to_lj(k)[1];
    }

    void Orbit::set_orbit_nlje(vector<int> nlje) {
        n = nlje[0];
        l = nlje[1];
        j = nlje[2];
        e = nlje[3];
        k = _lj_to_k(l, j);
    }

    void Orbit::set_radial_function_type(string new_radial_function_type) {
        radial_function_type = new_radial_function_type;
    }

    
    double Orbit::eval_radial_function_rspace(double x, double zeta, double par, double PQ) {
        try {
            if (radial_function_type == "default") { return _laguerre_wave_function_rspace(x, zeta); }
            else if (radial_function_type == "LO") { return _laguerre_wave_function_rspace(x, zeta); }
            else if (radial_function_type == "HO") { return _HO_wave_function_rspace(x, zeta); }
            else if (radial_function_type == "STO") { return _STO_wave_function_rspace(x, zeta); }
            else if (radial_function_type == "LSpinor") {
                if (PQ == 1) { return _LSpinor_P_wave_function_rspace(x, zeta, par); }
                if (PQ ==-1) { return _LSpinor_Q_wave_function_rspace(x, zeta, par); }
            } else { throw radial_function_type; }
        }
        catch (string radial_function_type) { cout << "Unknown radial_function_type" << endl; }
        return(0);
    }
    
    
    double Orbit::eval_radial_function_rspace_dr(double x, double zeta, double par, double PQ) {
        try {
            if (radial_function_type == "default") { return _laguerre_wave_function_pspace(x, zeta); }
            else if (radial_function_type == "LO") { return _laguerre_wave_function_pspace(x, zeta); }
            else if (radial_function_type == "HO") {return _HO_wave_function_pspace(x, zeta); }
            if (radial_function_type == "LSpinor") {
            if (PQ == 1) { return _LSpinor_P_wave_function_rspace_dr(x, zeta, par); }
            if (PQ == -1) { return _LSpinor_Q_wave_function_rspace_dr(x, zeta, par); }
            } else { throw radial_function_type; }
        }
        catch(string radial_function_type) { cout << "Unknown radial_function_type" << endl; }
        return(0);
    }
    
    
    double Orbit::eval_radial_function_pspace(double p, double zeta) {
        try {
            if (radial_function_type == "default") { return _laguerre_wave_function_pspace(p, zeta); }
            else if (radial_function_type == "LO") { return _laguerre_wave_function_pspace(p, zeta); }
            else if (radial_function_type == "HO") { return _HO_wave_function_pspace(p, zeta); }
            else if (radial_function_type == "STO") { return _STO_wave_function_pspace(p, zeta);  
            } else { throw radial_function_type; }
        }
        catch (string radial_function_type) { cout << "Unknown radial_function_type" << endl; }
        return(0);
    }
    
    
    double Orbit::_laguerre_wave_function_rspace(double x, double zeta) {
        
        //Laguerre function, see[A.E.McCoy and M.A.Caprio, J.Math.Phys. 57, (2016).] for details
        
        double eta = 2.0 * x / zeta;
        return sqrt(2.0 * tgamma(n + 1) / (zeta * tgamma(n + 2 * l + 3))) * 2.0 * pow(eta, l) * exp(-0.5 * eta) *  gsl_sf_laguerre_n(n, 2.0 * l + 2.0, eta) / zeta * x;
    }
    
    double Orbit::_laguerre_wave_function_pspace(double p, double zeta) {
        
        //Brute force fourier transformation of Laguerre function
        
        double xi = zeta * p * 0.5;
        double x = 0.5 / sqrt(0.25 + xi * xi);
        double r = 0.0;
        double prefact = sqrt(zeta * tgamma(n + 1) * tgamma(n + 2 * l + 3) / M_PI) * (0.5 * zeta) * pow((2 * x), (2 * l + 3)) * pow((2 * xi), l) * tgamma(l + 1);
        for (int k = 0; k < n + 1; k++) {
            r += pow((-1), k) * (k + 1.0) * pow((2 * x), k) * gsl_sf_gegenpoly_n(k + 1, l + 1.0, x) / (tgamma(n + 1 - k) * tgamma(2 * l + 3 + k));
        }
        return r * prefact * p;
    }
    
    double Orbit::_HO_wave_function_rspace(double x, double zeta) {
        double eta = x / zeta;
        return sqrt((2.0 / zeta) * (tgamma(n + 1) / tgamma(n + l + 1.5))) * (1.0 / zeta) * pow(eta, l) * exp(-0.5 * eta * eta) *  gsl_sf_laguerre_n( n, l + 0.5, eta * eta) * x;
    }

    double Orbit::_HO_wave_function_pspace(double p, double zeta) {
        double scale = 1.0 / zeta;
        return pow((-1), (n)) * _HO_wave_function_rspace(p, scale) * p;
    }
    

    double Orbit::_STO_wave_function_rspace(double x, double zeta) {
        double eta = x / zeta;
        return sqrt(1 / (pow(zeta, 3) * tgamma(2 * n + 2 * l + 1))) * pow(eta, (n + l - 1)) * exp(-0.5 * eta);
    }

    double Orbit::_STO_wave_function_pspace(double p, double zeta) {
        double xi = zeta * p;
        return pow(2, l) * sqrt(2 * pow(zeta, 3) / (M_PI * tgamma(2 * n + 2 * l + 1))) * tgamma(n + 1) * tgamma(l + 1) * pow(xi, l) / pow(sqrt(pow(xi, 2) + 0.25), (n + 2 * l + 2)) * gsl_sf_gegenpoly_n(n, l + 1.0, 0.5 / sqrt(pow(xi, 2) + 0.25));

    }

    double Orbit::_LSpinor_P_wave_function_rspace(double x, double zeta, double Z) {
        double eta = 2.0 * x / zeta;
        double gam = _get_pars_Lspinor(zeta, Z)[0];
        double N = _get_pars_Lspinor(zeta, Z)[1];
        double Norm = _get_pars_Lspinor(zeta, Z)[2];
        if (Norm == 0.0) { return 0.0; } 
        double T = gsl_sf_laguerre_n(n, 2 * gam, eta) * (N - k) / (n + 2 * gam);
        if (n > 0) { T -= gsl_sf_laguerre_n(n - 1, 2 * gam, eta); } 
        return T * pow(eta, gam) * exp(-0.5 * eta) * Norm;
    }

    double Orbit::_LSpinor_Q_wave_function_rspace(double x, double zeta, double Z) {
        double eta = 2.0 * x / zeta;
        double gam = _get_pars_Lspinor(zeta, Z)[0];
        double N = _get_pars_Lspinor(zeta, Z)[1];
        double Norm = _get_pars_Lspinor(zeta, Z)[2];
        if (Norm == 0.0) { return 0.0; }
        double T = -gsl_sf_laguerre_n(n, 2 * gam, eta) * (N - k) / (n + 2 * gam);
        if (n > 0) { T -= gsl_sf_laguerre_n(n - 1, 2 * gam, eta); }
        return T * pow(eta, gam) * exp(-0.5 * eta) * Norm;
    }
    
    double Orbit::_LSpinor_P_wave_function_rspace_dr(double x, double zeta, double Z) {
        double eta = 2.0 * x / zeta;
        double gam = _get_pars_Lspinor(zeta, Z)[0];
        double N = _get_pars_Lspinor(zeta, Z)[1];
        double Norm = _get_pars_Lspinor(zeta, Z)[2];
        if (Norm == 0.0) { return 0.0; }
        double T1 = _LSpinor_P_wave_function_rspace(x, zeta, gam) * ((gam / eta) - 0.5);
        double T2;
        if (n > 0) {
            T2 = -gsl_sf_laguerre_n(n - 1, 2 * gam + 1, eta) * (N - k) / (n + 2 * gam);
            if (n > 1) {
                T2 += gsl_sf_laguerre_n(n - 2, 2 * gam + 1, eta);
            }
        } else {
            T2 = 0.0;
        }
        return 2.0 / zeta * (T1 + Norm * T2 * pow(eta, gam) * exp(-0.5 * eta));
    }

    double Orbit::_LSpinor_Q_wave_function_rspace_dr(double x, double zeta, double Z) {
        double eta = 2.0 * x / zeta;
        double gam = _get_pars_Lspinor(zeta, Z)[0];
        double N = _get_pars_Lspinor(zeta, Z)[1];
        double Norm = _get_pars_Lspinor(zeta, Z)[2];
        if (Norm == 0.0) { return 0.0; }
        double T1 = _LSpinor_Q_wave_function_rspace(x, zeta, gam) * ((gam / eta) - 0.5);
        double T2;
        if (n > 0) {
            T2 = gsl_sf_laguerre_n(n - 1, 2 * gam + 1, eta) * (N - k) / (n + 2 * gam);
            if (n > 1) {
                T2 += gsl_sf_laguerre_n(n - 2, 2 * gam + 1, eta);
            }
        } else {
            T2 = 0.0;
        }
        return 2.0 / zeta * (T1 + Norm * T2 * pow(eta, gam) * exp(-0.5 * eta));
    }

    vector<double> Orbit::_get_pars_Lspinor(double zeta, double Z) {
        double gam = sqrt(pow(k, 2) - pow(Z, 2) / pow(c, 2));
        double N = sqrt(pow(n, 2) + 2.0 * n * gam + pow(k, 2));
        double Norm;
        if (N - k == 0) {
            Norm = 0;
        }
        else {
            Norm = sqrt((tgamma(n + 1) * (2 * gam + n)) / (2 * N * (N - k) * tgamma(2 * gam + n) * zeta));
        }
        vector<double> result = { gam, N, Norm };
        return result;
    }

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//          General Non-Class Functions

    double ME_overlap(Orbit o1, Orbit o2, double zeta, double Z) {
        
        //overlap
        
        double result, error;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);

        if (o1.e == 1 and o2.e == 1) {

            // gsl_function_pp Fp([&](double x)->double {return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * o2.eval_radial_function_rspace(x, zeta, Z, o2.e)); });
            auto lambda = [&](double x) { return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * o2.eval_radial_function_rspace(x, zeta, Z, o2.e)); };
            gsl_function_pp<decltype(lambda)> Fp(lambda);
            gsl_function* F = static_cast<gsl_function*>(&Fp);
            
            gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 10000, w, &result, &error);
            gsl_integration_workspace_free(w);

            return result;

        } else if ((o1.e == -1 and o2.e == -1)){

            // gsl_function_pp Fp([&](double x)->double {return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * o2.eval_radial_function_rspace(x, zeta, Z, o2.e)); });
            auto lambda = [&](double x) {return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * o2.eval_radial_function_rspace(x, zeta, Z, o2.e)); };
            gsl_function_pp<decltype(lambda)> Fp(lambda);
            
            gsl_function* F = static_cast<gsl_function*>(&Fp);

            gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 10000, w, &result, &error);
            gsl_integration_workspace_free(w);

            return result;
                
        } else { return 0.0; }
    }

    double ME_NuclPot(Orbit o1, Orbit o2, double zeta, double Z) {
        
        // 1 / r
        
        double result, error;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);

        if (o1.e == 1 and o2.e == 1) {

            // gsl_function_pp Fp([&](double x)->double {return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * o2.eval_radial_function_rspace(x, zeta, Z, o2.e) / x); });
            auto lambda = [&](double x) { return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * o2.eval_radial_function_rspace(x, zeta, Z, o2.e) / x); };
            gsl_function_pp<decltype(lambda)> Fp(lambda);
            gsl_function* F = static_cast<gsl_function*>(&Fp);

            gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 10000, w, &result, &error);
            gsl_integration_workspace_free(w);
            return result;
        } else if (o1.e == -1 and o2.e == -1) {

            // gsl_function_pp Fp([&](double x)->double {return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * o2.eval_radial_function_rspace(x, zeta, Z, o2.e) / x); });
            auto lambda = [&](double x) { return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * o2.eval_radial_function_rspace(x, zeta, Z, o2.e) / x); };
            gsl_function_pp<decltype(lambda)> Fp(lambda);
        
            gsl_function* F = static_cast<gsl_function*>(&Fp);

            gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 10000, w, &result, &error);
            gsl_integration_workspace_free(w);
            return result;
        }
        else { return 0.0; }
    }

    double ME_Kinetic(Orbit o1, Orbit o2, double zeta, double Z) {
        
        //s dot p
        
        double result, error;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);

        if (o1.e == 1 and o2.e == -1) {

            // gsl_function_pp Fp([&](double x)->double {return (-o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * (o2.eval_radial_function_rspace_dr(x, zeta, Z, o2.e) - o2.k * o2.eval_radial_function_rspace(x, zeta, Z, o2.e) / x)); });
            auto lambda = [&](double x) { return (-o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * (o2.eval_radial_function_rspace_dr(x, zeta, Z, o2.e) - o2.k * o2.eval_radial_function_rspace(x, zeta, Z, o2.e) / x)); };
            gsl_function_pp<decltype(lambda)> Fp(lambda);
            
            gsl_function* F = static_cast<gsl_function*>(&Fp);

            gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 10000, w, &result, &error);
            gsl_integration_workspace_free(w);
            return result * c;

        } else if (o1.e == -1 and o2.e == 1) {

            // gsl_function_pp Fp([&](double x)->double {return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * (o2.eval_radial_function_rspace_dr(x, zeta, Z, o2.e) + o2.k * o2.eval_radial_function_rspace(x, zeta, Z, o2.e) / x)); });
            auto lambda = [&](double x) { return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * (o2.eval_radial_function_rspace_dr(x, zeta, Z, o2.e) + o2.k * o2.eval_radial_function_rspace(x, zeta, Z, o2.e) / x)); };
            gsl_function_pp<decltype(lambda)> Fp(lambda);
            
            gsl_function* F = static_cast<gsl_function*>(&Fp);

            gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 10000, w, &result, &error);
            gsl_integration_workspace_free(w);
            return result * c;

        } else { return 0.0; }
    }

    double ME_Kinetic_nonrel(Orbit o1, Orbit o2, double zeta) {
        
        //p^2 / 2m
        
        double result, error;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);

        // gsl_function_pp Fp([&](double x)->double {return (o1.eval_radial_function_pspace(x, zeta) * o2.eval_radial_function_pspace(x, zeta) * pow(x, 2) * 0.5); });
        auto lambda = [&](double x) { return (o1.eval_radial_function_pspace(x, zeta) * o2.eval_radial_function_pspace(x, zeta) * pow(x, 2) * 0.5); };
        gsl_function_pp<decltype(lambda)> Fp(lambda);
        
        gsl_function* F = static_cast<gsl_function*>(&Fp);

        gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 10000, w, &result, &error);
        gsl_integration_workspace_free(w);
        return result;
    }

    double norm_check_rspace(Orbit o1, Orbit o2, double zeta) {
        if (o1.e != o2.e) { return 0.0; }
        if (o1.l != o2.l) { return 0.0; }
        if (o1.j != o2.j) { return 0.0; }

        double result, error;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);

        // gsl_function_pp Fp([&](double x)->double {return (o1.eval_radial_function_rspace(x, zeta) * o2.eval_radial_function_rspace(x, zeta) * x * x); });
        auto lambda = [&](double x) { return (o1.eval_radial_function_rspace(x, zeta) * o2.eval_radial_function_rspace(x, zeta) * x * x); };
        gsl_function_pp<decltype(lambda)> Fp(lambda);
        
        gsl_function* F = static_cast<gsl_function*>(&Fp);

        gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 10000, w, &result, &error);
        gsl_integration_workspace_free(w);
        return result;

    }

    double norm_check_rspace_rel(Orbit o1, Orbit o2, double zeta, double par) {
        if (o1.e != o2.e) { return 0.0; }
        if (o1.j != o2.j) { return 0.0; }
        if (o1.e != o2.e) { return 0.0; }

        int PQ;
        double result1, error1, result2, error2;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);
        gsl_integration_workspace* m = gsl_integration_workspace_alloc(10000);

        // gsl_function_pp Fp([&](double x)->double {return (o1.eval_radial_function_pspace(x, zeta) * o2.eval_radial_function_pspace(x, zeta) * x * x); });
        auto lambda1 = [&](double x) { return (o1.eval_radial_function_rspace(x, zeta, par, o1.e) * o2.eval_radial_function_rspace(x, zeta, par, o2.e)); };
        gsl_function_pp<decltype(lambda1)> Fp(lambda1);
        gsl_function* F = static_cast<gsl_function*>(&Fp);
        gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 10000, w, &result1, &error1);
        gsl_integration_workspace_free(w);

        auto lambda2 = [&](double x) { return (o1.eval_radial_function_rspace(x,zeta,par,-o1.e) * o2.eval_radial_function_rspace(x,zeta,par,-o2.e)); };
        gsl_function_pp<decltype(lambda2)> Xp(lambda2);
        gsl_function* X = static_cast<gsl_function*>(&Xp);
        gsl_integration_qagiu(X, 0, 1e-6, 1e-6, 10000, m, &result2, &error2);
        gsl_integration_workspace_free(m);

        return result1+result2;
    }

    double norm_check_pspace(Orbit o1, Orbit o2, double zeta) {
        if (o1.l != o2.l) { return 0.0; }
        if (o1.j != o2.j) { return 0.0; }
        if (o2.e != o2.e) { return 0.0; }
        
        double result, error;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);
        auto lambda = [&](double x) { return (o1.eval_radial_function_pspace(x,zeta) * o2.eval_radial_function_pspace(x,zeta) * x*x); };
        gsl_function_pp<decltype(lambda)> Fp(lambda);
        gsl_function* F = static_cast<gsl_function*>(&Fp);
        gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 10000, w, &result, &error);
        gsl_integration_workspace_free(w);
        
        return result;
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //          Orbits Class Methods


    void Orbits::add_orbit(vector<int> nlje){
        for (int i : nlje ){
            if (verbose==true) { cout << "The orbit " << nlje[0] << nlje[1] << nlje[2] << nlje[3] << " is already there." << endl; }
        }
        norbs = orbits.size();
        int idx = norbs;
        nlje_idx[nlje] = idx;
        Orbit orb = Orbit();
        orb.set_orbit_nlje(nlje);
        orb.set_radial_function_type(radial_function_type);
        orbits.push_back(orb);
        norbs = orbits.size();
        // cout << emax << nlje[0]+nlje[1] << endl;
        emax = max(emax, nlje[0]+nlje[1]);
        lmax = max(lmax, nlje[1]);
    }


    void Orbits::add_orbit_from_label(string value) {
    
        // string format should be like l0s1 => large-compoennt 0s1/2, s0s1 => small-component 0s1/2
        
        int z;
        char pn = value.at(0);
        if ( pn == 's' ) { z = -1; }
        else if ( pn == 'l' ) { z = 1; }
        else{ 
            cout << "parse error in add_orbit_from_label: " << value << endl;
            return;
        }
        string nlj_str = value.erase(0,1);

        char n_str = value.at(0);
        char l_str = value.at(1);
        char j_str = value.at(2);

        // regex regexp1("[a-z_]+");
        // smatch l_str;
        // regex_search(nlj_str, l_str, regexp1); 

        // regex regexp2("[0-9]+");
        // smatch n_str, j_str;
        // regex_search(nlj_str, n_str, regexp2); 

        // string nlj_str_v2 = nlj_str.erase(0,1);
        // regex_search(nlj_str_v2, j_str, regexp2); 

        n = int(n_str);
        int l = 0;
        for ( char l_label : _labels_orbital_angular_momentum ) {
            if ( l_str == l_label ) { break; }
            l += 1;
        }
        j = int(j_str);
        vector <int> nlje = {n,l,j,e};
        add_orbit(nlje);
    }

    void Orbits::add_orbits_from_labels(vector <string> strings){
        for (string label : strings) { add_orbit_from_label(label); }
    }

    Orbit Orbits::get_orbit(int idx) {
        return orbits[idx];
    }

    string Orbits::get_orbit_label(int idx) {
        Orbit o = get_orbit(idx);
        char pn = 'l';
        if (o.e == 1) { pn = 'u'; }
        return pn + to_string(o.n) + _labels_orbital_angular_momentum[o.l] + to_string(o.j);
    }
    
    int Orbits::get_orbit_index(vector<int> nlj) {
        // nlje_idx.push_back(nlj);
        // nlje_idx[nlj];
        // mymap.insert(pair<int,vector<int> >(10, vector<int>()));
        return nlje_idx[nlj];
    } 
    
    int Orbits::get_orbit_index_from_orbit(Orbit o) {
        vector <int> idx = {o.n, o.l, o.j, o.e};
        return get_orbit_index(idx);
    }

    // int Orbits::get_orbit_index_from_tuple(tuple<int,int,int,int> nlje){
    //     nlje_idx.insert(nlje_idx.end(), nlje.begin(), nlje.end());
    //     return nlje_idx;
    // }
    
    int Orbits::get_num_orbits() {
        return norbs;
    }

    void Orbits::set_orbits(int nmax, int lmax) {
        vector<int> klist;
        if (relativistic==true) { klist = {1,-1}; }
        if (relativistic==false) { klist = {1}; }
        for ( int li=0; li < lmax+1; li++ ){
            for ( int n =0; n < nmax+1; n++){
                for (int j : {2*li-1, 2*li+1} ) {
                    if ( j<0 ){ continue; }
                    for ( int e : klist ) { 
                        vector<int> nlje = {n,li,j,e};
                        add_orbit(nlje);
                    }
                }
            }
        }
    }

    void Orbits::append_orbits( Orbits orbs ) {
        for (Orbit o : orbs.orbits ) {
            vector<int> nlje = {n,l,j,e};
            // o.get_nljke();
            add_orbit(nlje);
        }
    }

    void Orbits::print_orbits(){
        cout << "Orbit list" << endl;
        cout << " idx, n, l, j, k, e" << endl;
        for (Orbit o : orbits){
            int idx;
            vector<int> nljke = o.get_nljke();
            idx = get_orbit_index(nljke);
            // idx = get_orbit_index_from_tuple( {nljke} );
            idx = get_orbit_index_from_orbit(o);
            cout << fixed << setw(3) << idx << setw(4) << nljke[0] << setw(3) << nljke[1] << setw(3) << nljke[2] << setw(3) << nljke[3] << setw(3) << nljke[4] << endl;
        }
    }
