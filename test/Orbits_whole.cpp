#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
//#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf.h>
#include <armadillo>
#include <math.h>
#include <cmath>
#include <string>
#include <regex>
// #include "Orbits.h"

using namespace std;
using namespace arma;
//using namespace boost::range;
// std::string;


    

        
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

    
class Orbit {

    private:


    public:
        double c = 137.035999084;
        // Member variables
        // double c;
        int n, l, j, k, e;
        string radial_function_type;
        
        Orbit() {
            // Constructor 
            // double c = 137.035999084;
            n = -1; // nodal
            l = -1; // orbital angular momentum
            j = 0;  // total angular momentum
            k = 0;  // angular
            e = 0;  // upper / lower 1 / -1
            radial_function_type = "default";
        }

        vector<int> get_nljke() {  

            vector<int> result = { n, l, j, k, e };
            return result;
        }

        int _lj_to_k(int l, int j) {

            return floor( - (4.0 + j * (j + 2.0) - 4.0 * l * (l + 1.0) - 3.0) / 4.0);
        }

        vector<int> _k_to_lj(int k) {

        
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

        void set_orbit_nke(vector<int> nke) {


            n = nke[0];
            k = nke[1];
            e = nke[2];

            l = _k_to_lj(k)[0];
            j = _k_to_lj(k)[1];

            // cout << n << k << e << l << j << endl;

            // return 0;

        }

        void set_orbit_nlje(vector<int> nlje) {
            n = nlje[0];
            l = nlje[1];
            j = nlje[2];
            e = nlje[3];

            k = _lj_to_k(l, j);
        }

        void set_radial_function_type(string new_radial_function_type) {
            radial_function_type = new_radial_function_type;
        }

        
        double eval_radial_function_rspace(double x, double zeta, double par = 0.0, double PQ = 1.0) {
            
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
        
        
        double eval_radial_function_rspace_dr(double x, double zeta, double par = 0.0, double PQ = 1.0) {
            
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
        
        
        double eval_radial_function_pspace(double p, double zeta) {
            
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
        
        
        double _laguerre_wave_function_rspace(double x, double zeta) {
            
            //Laguerre function, see[A.E.McCoy and M.A.Caprio, J.Math.Phys. 57, (2016).] for details
            
            double eta = 2.0 * x / zeta;
            return sqrt(2.0 * tgamma(n + 1) / (zeta * tgamma(n + 2 * l + 3))) * 2.0 * pow(eta, l) * exp(-0.5 * eta) *  gsl_sf_laguerre_n(n, 2.0 * l + 2.0, eta) / zeta * x;
        }
        
        double _laguerre_wave_function_pspace(double p, double zeta) {
            
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
        
        double _HO_wave_function_rspace(double x, double zeta) {
            double eta = x / zeta;
            return sqrt((2.0 / zeta) * (tgamma(n + 1) / tgamma(n + l + 1.5))) * (1.0 / zeta) * pow(eta, l) * exp(-0.5 * eta * eta) *  gsl_sf_laguerre_n( n, l + 0.5, eta * eta) * x;
        }

        double _HO_wave_function_pspace(double p, double zeta) {
            double scale = 1.0 / zeta;
            return pow((-1), (n)) * _HO_wave_function_rspace(p, scale) * p;
        }
        

        double _STO_wave_function_rspace(double x, double zeta) {
            double eta = x / zeta;
            return sqrt(1 / (pow(zeta, 3) * tgamma(2 * n + 2 * l + 1))) * pow(eta, (n + l - 1)) * exp(-0.5 * eta);
        }

        double _STO_wave_function_pspace(double p, double zeta) {
            double xi = zeta * p;
            return pow(2, l) * sqrt(2 * pow(zeta, 3) / (M_PI * tgamma(2 * n + 2 * l + 1))) * tgamma(n + 1) * tgamma(l + 1) * pow(xi, l) / pow(sqrt(pow(xi, 2) + 0.25), (n + 2 * l + 2)) * gsl_sf_gegenpoly_n(n, l + 1.0, 0.5 / sqrt(pow(xi, 2) + 0.25));

        }

        double _LSpinor_P_wave_function_rspace(double x, double zeta, double Z) {
            double eta = 2.0 * x / zeta;
            double gam = _get_pars_Lspinor(zeta, Z)[0];
            double N = _get_pars_Lspinor(zeta, Z)[1];
            double Norm = _get_pars_Lspinor(zeta, Z)[2];
            if (Norm == 0.0) { return 0.0; } 
            double T = gsl_sf_laguerre_n(n, 2 * gam, eta) * (N - k) / (n + 2 * gam);
            if (n > 0) { T -= gsl_sf_laguerre_n(n - 1, 2 * gam, eta); } 
            return T * pow(eta, gam) * exp(-0.5 * eta) * Norm;
        }

        double _LSpinor_Q_wave_function_rspace(double x, double zeta, double Z) {
            double eta = 2.0 * x / zeta;
            double gam = _get_pars_Lspinor(zeta, Z)[0];
            double N = _get_pars_Lspinor(zeta, Z)[1];
            double Norm = _get_pars_Lspinor(zeta, Z)[2];
            if (Norm == 0.0) { return 0.0; }
            double T = -gsl_sf_laguerre_n(n, 2 * gam, eta) * (N - k) / (n + 2 * gam);
            if (n > 0) { T -= gsl_sf_laguerre_n(n - 1, 2 * gam, eta); }
            return T * pow(eta, gam) * exp(-0.5 * eta) * Norm;
        }
        
        double _LSpinor_P_wave_function_rspace_dr(double x, double zeta, double Z) {
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

        double _LSpinor_Q_wave_function_rspace_dr(double x, double zeta, double Z) {
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

        vector<double> _get_pars_Lspinor(double zeta, double Z) {
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
};

    /*
    double ME_overlap(Orbit o1, Orbit o2, double zeta, double Z);
    double ME_NuclPot(Orbit o1, Orbit o2, double zeta, double Z);
    double ME_Kinetic(Orbit o1, Orbit o2, double zeta, double Z, double c);
    double ME_Kinetic_nonrel(Orbit o1, Orbit o2, double zeta);
    double norm_check_rspace(Orbit o1, Orbit o2, double zeta);
    double norm_check_rspace_rel(Orbit o1, Orbit o2, double zeta, double par);
    */

    double ME_overlap(Orbit o1, Orbit o2, double zeta, double Z) {
        
        //overlap
        
        double result, error;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

        if (o1.e == 1 and o2.e == 1) {

            // gsl_function_pp Fp([&](double x)->double {return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * o2.eval_radial_function_rspace(x, zeta, Z, o2.e)); });
            auto lambda = [&](double x) { return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * o2.eval_radial_function_rspace(x, zeta, Z, o2.e)); };
            gsl_function_pp<decltype(lambda)> Fp(lambda);
            gsl_function* F = static_cast<gsl_function*>(&Fp);
            
            gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 1000, w, &result, &error);
            gsl_integration_workspace_free(w);

            return result;

        } else if ((o1.e == -1 and o2.e == -1)){

            // gsl_function_pp Fp([&](double x)->double {return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * o2.eval_radial_function_rspace(x, zeta, Z, o2.e)); });
            auto lambda = [&](double x) {return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * o2.eval_radial_function_rspace(x, zeta, Z, o2.e)); };
            gsl_function_pp<decltype(lambda)> Fp(lambda);
            
            gsl_function* F = static_cast<gsl_function*>(&Fp);

            gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 1000, w, &result, &error);
            gsl_integration_workspace_free(w);

            return result;
                
        } else { return 0.0; }
    }

    double ME_NuclPot(Orbit o1, Orbit o2, double zeta, double Z) {
        
        // 1 / r
        

        double result, error;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

        if (o1.e == 1 and o2.e == 1) {

            // gsl_function_pp Fp([&](double x)->double {return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * o2.eval_radial_function_rspace(x, zeta, Z, o2.e) / x); });
            auto lambda = [&](double x) { return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * o2.eval_radial_function_rspace(x, zeta, Z, o2.e) / x); };
            gsl_function_pp<decltype(lambda)> Fp(lambda);
            gsl_function* F = static_cast<gsl_function*>(&Fp);

            gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 1000, w, &result, &error);
            gsl_integration_workspace_free(w);
            return result;
        } else if (o1.e == -1 and o2.e == -1) {

            // gsl_function_pp Fp([&](double x)->double {return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * o2.eval_radial_function_rspace(x, zeta, Z, o2.e) / x); });
            auto lambda = [&](double x) { return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * o2.eval_radial_function_rspace(x, zeta, Z, o2.e) / x); };
            gsl_function_pp<decltype(lambda)> Fp(lambda);
        
            gsl_function* F = static_cast<gsl_function*>(&Fp);

            gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 1000, w, &result, &error);
            gsl_integration_workspace_free(w);
            return result;
        }
        else { return 0.0; }
    }

    double ME_Kinetic(Orbit o1, Orbit o2, double zeta, double Z, double c) {
        
        //s dot p
        
        double result, error;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

        if (o1.e == 1 and o2.e == -1) {

            // gsl_function_pp Fp([&](double x)->double {return (-o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * (o2.eval_radial_function_rspace_dr(x, zeta, Z, o2.e) - o2.k * o2.eval_radial_function_rspace(x, zeta, Z, o2.e) / x)); });
            auto lambda = [&](double x) { return (-o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * (o2.eval_radial_function_rspace_dr(x, zeta, Z, o2.e) - o2.k * o2.eval_radial_function_rspace(x, zeta, Z, o2.e) / x)); };
            gsl_function_pp<decltype(lambda)> Fp(lambda);
            
            gsl_function* F = static_cast<gsl_function*>(&Fp);

            gsl_integration_qagiu(F, 0, 1e-7, 1000, GSL_INTEG_GAUSS61, w, &result, &error);
            gsl_integration_workspace_free(w);
            return result * c;

        } else if (o1.e == -1 and o2.e == 1) {

            // gsl_function_pp Fp([&](double x)->double {return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * (o2.eval_radial_function_rspace_dr(x, zeta, Z, o2.e) + o2.k * o2.eval_radial_function_rspace(x, zeta, Z, o2.e) / x)); });
            auto lambda = [&](double x) { return (o1.eval_radial_function_rspace(x, zeta, Z, o1.e) * (o2.eval_radial_function_rspace_dr(x, zeta, Z, o2.e) + o2.k * o2.eval_radial_function_rspace(x, zeta, Z, o2.e) / x)); };
            gsl_function_pp<decltype(lambda)> Fp(lambda);
            
            gsl_function* F = static_cast<gsl_function*>(&Fp);

            gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 1000, w, &result, &error);
            gsl_integration_workspace_free(w);
            return result * c;

        } else { return 0.0; }
    }

    double ME_Kinetic_nonrel(Orbit o1, Orbit o2, double zeta) {
        
        //p^2 / 2m
        
        double result, error;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

        // gsl_function_pp Fp([&](double x)->double {return (o1.eval_radial_function_pspace(x, zeta) * o2.eval_radial_function_pspace(x, zeta) * pow(x, 2) * 0.5); });
        auto lambda = [&](double x) { return (o1.eval_radial_function_pspace(x, zeta) * o2.eval_radial_function_pspace(x, zeta) * pow(x, 2) * 0.5); };
        gsl_function_pp<decltype(lambda)> Fp(lambda);
        
        gsl_function* F = static_cast<gsl_function*>(&Fp);

        gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 1000, w, &result, &error);
        gsl_integration_workspace_free(w);
        return result;
    }

    double norm_check_rspace(Orbit o1, Orbit o2, double zeta) {
        if (o1.e != o2.e) { return 0.0; }
        if (o1.l != o2.l) { return 0.0; }
        if (o1.j != o2.j) { return 0.0; }

        double result, error;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

        // gsl_function_pp Fp([&](double x)->double {return (o1.eval_radial_function_rspace(x, zeta) * o2.eval_radial_function_rspace(x, zeta) * x * x); });
        auto lambda = [&](double x) { return (o1.eval_radial_function_rspace(x, zeta) * o2.eval_radial_function_rspace(x, zeta) * x * x); };
        gsl_function_pp<decltype(lambda)> Fp(lambda);
        
        gsl_function* F = static_cast<gsl_function*>(&Fp);

        gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 1000, w, &result, &error);
        gsl_integration_workspace_free(w);
        return result;

    }

    double norm_check_rspace_rel(Orbit o1, Orbit o2, double zeta, double par) {
        if (o1.e != o2.e) { return 0.0; }
        if (o1.j != o2.j) { return 0.0; }
        if (o1.e != o2.e) { return 0.0; }

        double result, error;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

        // gsl_function_pp Fp([&](double x)->double {return (o1.eval_radial_function_pspace(x, zeta) * o2.eval_radial_function_pspace(x, zeta) * x * x); });
        auto lambda = [&](double x) { return (o1.eval_radial_function_pspace(x, zeta) * o2.eval_radial_function_pspace(x, zeta) * x * x); };
        gsl_function_pp<decltype(lambda)> Fp(lambda);
        
        gsl_function* F = static_cast<gsl_function*>(&Fp);

        gsl_integration_qagiu(F, 0, 1e-6, 1e-6, 1000, w, &result, &error);
        gsl_integration_workspace_free(w);
        return result;
    }



/*
class Orbits: public Orbit {

    public:
        vector <int> nlje_idx = {};
        vector <int> orbits = {};
        int norbs = -1;
        int emax = -1;
        int lmax = -1;
        vector <char> _labels_orbital_angular_momentum = { 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'o', 'q', 'r', 't', 'u', 'v', 'w', 'x', 'y', 'z' };
        int relativistic = 1; // True 
        string radial_function_type = "default";



        void add_orbit(vector<double> nlje) {
        
            for (auto nlje : nlje_idx) {
            
                for (auto i = nlje_idx.begin(); i != nlje_idx.end(); ++i) {
                    cout << "The orbit" << " " << *i << " " << "is already there.";
                }
                    
            }
        }

        
        //int norbs = orbits.size();
        //Orbit orb;
        //orb.set_orbit_nlje(nlje)
        

        void add_orbit_from_label(string value) {
            
            //string format should be like l0s1 => large-compoennt 0s1/2, s0s1 => small-component 0s1/2
            
        
            char pn = value.at(0);
            int z;
            try {
                if (pn == 's') {
                    z = -1;
                }
                else if (pn == 'l') {
                    z = 1;
                }
                else { throw pn;  }
            }
            catch (char pn) { cout << "parse error in add_orbit_from_label: " << value << endl; }

            string nlj_str = value.erase(0);
            
            string l_str, matchNum;
            regex regexp("[a-z]+");
            regex_search(nlj_str, l_str, regexp);

            regex regexpNum("[0-9]+");
            regex_search(nlj_str, matchNum, regexpNum);
            string n_str, j_str;

            int n = stoi(n_str);
            int l = 0;

            for (auto l_label : _labels_orbital_angular_momentum) {
                if (l_str[0] == l_label) { break; }
                l += 1;
            }
            int j = stoi(j_str);
            //self.add_orbit_nlje(n, l, j, z)
            // Line only comes up once, not sure of meaning 
        }


        void add_orbits_from_labels(vector <string> strings) {
            // vector <string> is wrong decleration of variable 
            for (auto label : strings) {
         
                return add_orbit_from_label(label);
            
            }
        }

        Orbits get_orbit(int idx) {
            orbits[idx];
        }

        string get_orbit_label(int idx) {
            Orbit o = get_orbit(idx);
            char pn = 'l';
            if (o.e = 1) { pn = 'u'; }
            return pn + to_string(o.n) + _labels_orbital_angular_momentum[o.l] + to_string(o.j);

        }


};
*/


