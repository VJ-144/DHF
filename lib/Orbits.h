#pragma once
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
#include<tuple>


using namespace std;

const double c = 137.035999084;

class Orbit {

    public:
        // Constructor
        Orbit(string radial_function_type1="default");

        // Member variables
        int n, l, j, k, e;
        double c = 137.035999084;
        string radial_function_type;

        // Class Member Declerations
        vector<int> get_nljke();
        int _lj_to_k(int l, int j);
        vector<int> _k_to_lj(int k);
        void set_orbit_nke(vector<int> nke);
        void set_orbit_nlje(vector<int> nlje);
        void set_radial_function_type(string new_radial_function_type);
        double eval_radial_function_rspace(double x, double zeta, double par = 0.0, int PQ = 1);
        double eval_radial_function_rspace_dr(double x, double zeta, double par = 0.0, int PQ = 1);
        double eval_radial_function_pspace(double p, double zeta);
        double _laguerre_wave_function_rspace(double x, double zeta);
        double _laguerre_wave_function_pspace(double p, double zeta);
        double _HO_wave_function_rspace(double x, double zeta);
        double _HO_wave_function_pspace(double p, double zeta);
        double _STO_wave_function_rspace(double x, double zeta);
        double _STO_wave_function_pspace(double p, double zeta);
        double _LSpinor_P_wave_function_rspace(double x, double zeta, double Z);
        double _LSpinor_Q_wave_function_rspace(double x, double zeta, double Z);
        double _LSpinor_P_wave_function_rspace_dr(double x, double zeta, double Z);
        double _LSpinor_Q_wave_function_rspace_dr(double x, double zeta, double Z);
        vector<double> _get_pars_Lspinor(double zeta, double Z);
        

        
};

    // Non-Class Function Declerations
    double ME_overlap(Orbit o1, Orbit o2, double zeta, double Z);
    double ME_NuclPot(Orbit o1, Orbit o2, double zeta, double Z);
    double ME_Kinetic(Orbit o1, Orbit o2, double zeta, double Z);
    double ME_Kinetic_nonrel(Orbit o1, Orbit o2, double zeta);
    double norm_check_rspace(Orbit o1, Orbit o2, double zeta);
    double norm_check_rspace_rel(Orbit o1, Orbit o2, double zeta, double par);



class Orbits : public Orbit {

    public:

        // Member Variables
        vector <Orbit> orbits;
        int emax, lmax, norbs;
        bool relativistic, verbose;
        string radial_function_type;
        map<vector<int>, int> nlje_idx;
        vector <char> _labels_orbital_angular_momentum;

        // Constructor 
        // Orbits(int emax1=0, int lmax1=0, string radial_function_type1="default", bool verbose1=false, bool relativistic1=true);
        Orbits(int emax1, int lmax1, string radial_function_type1, bool verbose1, bool relativistic1);
        Orbits(string radial_function_type1="default", bool verbose1=false, bool relativistic1=true);

        // Class Function Declerations
        void add_orbit(vector<int> nlje);
        void add_orbit_from_label(string value);
        void add_orbits_from_labels(vector <string> strings);
        Orbit get_orbit(int idx);
        string get_orbit_label(int idx);
        int get_orbit_index( vector<int> nlj);
        int get_orbit_index_from_orbit(Orbit o);
        int get_orbit_index_from_tuple(tuple<int,int,int,int> nlje); 
        int get_num_orbits();
        void set_orbits(int nmax, int lmax);
        void append_orbits( Orbits orbs );
        void print_orbits();

};
