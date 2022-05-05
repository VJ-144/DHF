#include <iostream>
#include <armadillo>
#include <boost/range/irange.hpp>
#include <cmath>
// #include <wignerSymbols.h>
// #include <wignerSymbols>
// #include "wignerSymbols/wignerSymbols-cpp.h"
// #include<gsl>

#include "../header/TwoBodyOperator.h"
#include "../header/TwoBodySpace.h"
#include "../header/HartreeFock.h"
#include "../header/ModelSpace.h"
#include "../header/Operator.h"
#include "../header/Orbits.h"

using namespace std;
using namespace boost;
using namespace arma;


    void Monopole::set_monopole2(){
        Orbits orbs = modelspace.orbits;
        int norbs = orbs.get_num_orbits();
        for (int i1=0; i1<norbs; i1++) {
            Orbit o1 = orbs.get_orbit(i1);
            for (int i3=0; i3<i1+1; i1++) {
                Orbit o3 = orbs.get_orbit(i3);
                for (int i2=0; i2<norbs; i2++) {
                    for (int i4=0; i4<norbs; i4++) {
                        Orbit o2 = orbs.get_orbit(i2);
                        Orbit o4 = orbs.get_orbit(i4);
                        if (o2.k != o4.k) {continue;}
                        if (o1.e + o2.e != o3.e + o4.e) {continue;}
                        idx_to_ijkl.push_back({i1,i2,i3,i4});
                        double norm = 1;
                        if (i1==i2) {norm *= sqrt(2);}
                        if (i3==i4) {norm *= sqrt(2);}
                        double v;
                        for (int J : irange( floor(abs(o1.j-o2.j)/2),  floor((o1.j+o2.j)/2)+1 ) ) {
                            if (i1==i2 && J%2==1) {continue;}
                            if (i3==i4 && J%2==1) {continue;}
                            // v += (2*J+1) * Ham.two.get_2bme_orbitsJ(o1,o2,o3,o4,J,J);
                        v *= norm / (o1.j+1);
                        v2.push_back(v);
                        }
                    }
                }
            }
        }        
    }

    void Monopole::print_monopole2() {
        for (int idx=0; idx<idx_to_ijkl.size(); idx++) {
            // cout << fixed << setw(2) << idx_to_ijkl[idx] << v2[idx] << endl;
            // idx_to_ijkl[idx] contains a vector so need to find a more efficient method
        }
    }

    void HartreeFock::solve() {
        for (int n_iter=0; n_iter<100; n_iter++) {
            // not sure if this should have double?
            double r = UpdateFock(n_iter);
            DiagonalizeFock();
            UpdateDensityMatrix();
            CalcEnergy();
            // _print_status(n_iter, detail=True);
            _print_status(n_iter, false);
            if (r < 1.8e-8);
        }
    }

    void HartreeFock::CalcEnergy() {
        // not sure if this should be its own type i.e. Orbits orbs
        Orbits orbs = modelspace.orbits;
        int norbs = orbs.get_num_orbits();
        double e1 = 0.0;
        double e2 = 0.0;
        for (int i=0; i<norbs; i++) {
            for (int j=0; j<norbs; j++) {
                Orbit oi = orbs.get_orbit(i);
                e1 += Ham.one[i,j] * rho[i,j] * (oi.j+1);
                e2 += V[i,j] * rho[i,j] * 0.5 * (oi.j+1);
            }
        }
        En = e1+e2;
    }

    double HartreeFock::UpdateFock(int n_iter) {
        Mat<double> Fock_old = F;
        // not sure if this should be its own type i.e. Orbits orbs
        Orbits orbs = modelspace.orbits;
        int norbs = orbs.get_num_orbits();
        Mat<double> V(norbs, norbs, fill::zeros);
        if (n_iter != -1) {
            for (int idx=0; idx<monopole.idx_to_ijkl.size(); idx++) {
                int i = monopole.idx_to_ijkl[idx][0];
                int j = monopole.idx_to_ijkl[idx][1];
                int k = monopole.idx_to_ijkl[idx][2];
                int l = monopole.idx_to_ijkl[idx][3];
                V(i,k) += monopole.v2[idx] * rho[j,l];
            }
            for (int i=0; i<norbs; i++) {
                for (int j=0; i<i+j; j++) {
                    V(j,i) = V(i,j);
                }
            }
        }
        F = Ham.one + V;
        double r;
        for (int i=0; i<norbs; i++) {
            for (int j=0; j<norbs; j++) {
                r += sqrt( pow(( F(i,j) - Fock_old(i,j) ), 2) );
            }
        }
        return r;
    }

    void HartreeFock::DiagonalizeFock() {
        OneBodySpace one_body_space = modelspace.one;
        for (int ich=0; ich< one_body_space.number_channels; ich++) {
            vector <int> filter_idx = one_body_space.channels[ich];
            // Figure out Hashmap stuff

        }
    }

    void HartreeFock::UpdateDensityMatrix() {
        OneBodySpace one_body_space = modelspace.one;
        Orbits orbs = modelspace.orbits;
        int norbs = orbs.get_num_orbits();
        Mat<double> tmp(norbs, norbs, fill::zeros);
        Orbits orbs_tmp = Orbits();
        map<int, int> orbit_idx_to_spe_idx;
        for (int ich=0; ich<one_body_space.number_channels; ich++) {
            vector <int> filter_idx = one_body_space.channels[ich];
            // Fix SPEsCh
            vec SPEsCh;
            for (int i=0; i<one_body_space.channels[ich].size(); i++) {
                int idx_spe = one_body_space.channels[ich][i];
                if (orbs.relativistic) {
                    if (i < floor(SPEsCh.size()/2) ) {
                        n = i;
                        e =-1;
                    } else {
                        n = i - floor(SPEsCh.size()/2);
                        e = 1;
                    }
                } else {
                    n = i;
                    e = 1;
                }
                Orbit o = orbs.get_orbit(idx_spe);
                // double check its ok to ignore tuple
                int idx = orbs.get_orbit_index({n+o.l, o.l, o.j, e});
                orbit_idx_to_spe_idx[idx] = idx_spe;
            }
        }
        for ( auto hole_idx : holes) {
            // int occ = holes.first[hole_idx];
            // int spe_idx = orbit_idx_to_spe_idx[hole_idx];

        }



    }
