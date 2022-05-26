#include <iostream>
#include <ostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <regex>
#include <armadillo>

#include "../lib/Orbits.h"
// #include "Orbits_whole.cpp"


using namespace std;
using namespace arma;


int main() {

	//Orbits orbs = Orbits();
	string radial_function_type="LSpinor";
	bool relativistic=true;
	bool verbose=false;
	int emax = 0;
	int lmax = 0;
	int nmax = 1;

	Orbits orbs = Orbits( radial_function_type, verbose, relativistic);

	// string oa = "l0s1";
	// orbs.add_orbit_from_label(oa); 
	// orbs.print_orbits();


	orbs.set_orbits(1,0);
	// double zeta = 0.5
	int zeta = 1;
	for (int i=0; i<orbs.get_num_orbits(); i++){
		for (int j=0; j<orbs.get_num_orbits(); j++){
			Orbit oi = orbs.get_orbit(i);
			Orbit oj = orbs.get_orbit(j);
			//cout << i << " " << j << endl;
			if (i>j) { continue; }
			if (oi.k != oj.k) { continue; }
			if (oi.e != oj.e) { continue; }
			cout << i << " " << j << " " << norm_check_rspace_rel(oi,oj,zeta,1) << endl;
		}
	}
	
	int fac;
	mat Hmat(orbs.get_num_orbits(), orbs.get_num_orbits(), fill::zeros);
	mat Smat(orbs.get_num_orbits(), orbs.get_num_orbits(), fill::zeros);
	for (int i=0; i<orbs.get_num_orbits(); i++){
		for (int j=0; j<orbs.get_num_orbits(); j++){
			Orbit oi = orbs.get_orbit(i);
			Orbit oj = orbs.get_orbit(j);
			if (oi.e == 1 and oj.e == 1) { fac = 1; }
			if (oi.e == 1 and oj.e ==-1) { fac =-1; }
			Hmat(i,j) = ME_overlap(oi, oj, zeta, 1) * pow(orbs.c, 2) * (fac-1) - ME_NuclPot(oi, oj, zeta, 1) + ME_Kinetic(oi, oj, zeta, 1);
			Smat(i,j) = ME_overlap(oi, oj, zeta, 1);
		}
	}

	cx_dvec eigval;
	Mat<cx_double> eigvec_col;
	Mat<cx_double> eigvec, eigvec_norm, leigvec, reigvec;
	eig_pair(eigval, eigvec_col, Hmat, Smat);
	Mat<cx_double> eigvec_row = eigvec_col.t();
	eigval.print("eigenvalues");
	
	// Normalisation for non-orthonormal states
	// uvec idxs = sort_index(real(eigval), "ascend");
    // int idx = idxs(Nmax);
    // complex<double> norm = as_scalar(eigvec_row.row(idxs) * Smat * eigvec_col.col(idxs));
    


	return(0);
}

// g++ Test.cpp Orbits.cpp -o Test

//  g++ orbitsTest.cpp Orbits.cpp -lgsl -larmadillo -o Test
