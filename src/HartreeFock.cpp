#include <iostream>
#include <armadillo>
#include <boost/range/irange.hpp>
#include <cmath>

#include "../lib/TwoBodyOperator.h"
#include "../lib/TwoBodySpace.h"
#include "../lib/HartreeFock.h"
#include "../lib/ModelSpace.h"
#include "../lib/Operator.h"
#include "../lib/Orbits.h"

using namespace std;
// using namespace boost;
using namespace arma;


    Monopole::Monopole(Operator Ham1) {
        Ham=Ham1;
        modelspace = Ham.modelspace;
        idx_to_ijkl;
        v2;
    };


    HartreeFock::HartreeFock(Operator Ham1, map <int, double> holes1) {
        // En; // Not sure if this belongs here
        Ham = Ham1;
        holes = holes1;
        modelspace = Ham.modelspace;
        monopole = Monopole(Ham);
        monopole.set_monopole2();

        Orbits orbs = modelspace.orbits;
        int norbs = orbs.get_num_orbits();
        S = Ham1.S;
        C.zeros(norbs, norbs);
        rho.zeros(norbs, norbs);
        F.zeros(norbs, norbs);
        V.zeros(norbs, norbs);
        SPEs.zeros(norbs);

        double r = UpdateFock();
        DiagonalizeFock();
        UpdateDensityMatrix();
        CalcEnergy();
    };

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    void Monopole::set_monopole2(){
        Orbits orbs = modelspace.orbits;
        int norbs = orbs.get_num_orbits();
        for (int i1=0; i1<norbs; i1++) {
            Orbit o1 = orbs.get_orbit(i1);
            for (int i3=0; i3<i1+1; i3++) {
                Orbit o3 = orbs.get_orbit(i3);
                if (o1.k != o3.k) {continue;}
                for (int i2=0; i2<norbs; i2++) {
                    for (int i4=0; i4<norbs; i4++) {
                        Orbit o2 = orbs.get_orbit(i2);
                        Orbit o4 = orbs.get_orbit(i4);
                        if (o2.k != o4.k) {continue;}
                        if (o1.e + o2.e != o3.e + o4.e) {continue;}
                        idx_to_ijkl.push_back({i1,i2,i3,i4});
                        double norm = 1.0;
                        if (i1==i2) {norm *= sqrt(2);}
                        if (i3==i4) {norm *= sqrt(2);}
                        double v = 0.0;
                        for (int J : boost::irange( floor(abs(o1.j-o2.j)/2),  floor((o1.j+o2.j)/2)+1 ) ) {
                            if (i1==i2 && J%2==1) {continue;}
                            if (i3==i4 && J%2==1) {continue;}
                            v += (2*J+1) * Ham.two.get_2bme_orbitsJ(o1,o2,o3,o4,J,J);
                            // cout << Ham.two.get_2bme_orbitsJ(o1,o2,o3,o4,J,J) << endl;
                        v *= norm / (o1.j+1);
                        // cout << v << endl;
                        v2.push_back(v);
                        }
                    }
                }
            }
        }        
    }

    void Monopole::print_monopole2() {
        for (int idx=0; idx<idx_to_ijkl.size(); idx++) {
            printf("%d %d %d %d", idx_to_ijkl[idx][0], idx_to_ijkl[idx][1], idx_to_ijkl[idx][2], idx_to_ijkl[idx][3]);
            // idx_to_ijkl[idx] contains a vector so need to find a more efficient method
            // need to add v2 but not sure about length of vector
        }
    }

    void HartreeFock::solve() {
        for (int n_iter=0; n_iter<100; n_iter++) {
            double r = UpdateFock(n_iter);
            // cout << r << endl;
            DiagonalizeFock();
            UpdateDensityMatrix();
            CalcEnergy();
            // _print_status(n_iter);
            if (r < 1.8e-8) {break;};
        }
    }

    void HartreeFock::CalcEnergy() {
        Orbits orbs = modelspace.orbits;
        int norbs = orbs.get_num_orbits();
        double e1 = 0.0;
        double e2 = 0.0;
        // cout << norbs << endl;
        for (int i=0; i<norbs; i++) {
            for (int j=0; j<norbs; j++) {
                Orbit oi = orbs.get_orbit(i);
                e1 += Ham.one(i,j) * rho(i,j) * (oi.j+1);
                e2 += V(i,j) * rho(i,j) * 0.5 * (oi.j+1);
            }
        }
        En = e1+e2;
    }

    double HartreeFock::UpdateFock(int n_iter) {
        Mat<double> Fock_old = F;
        Orbits orbs = modelspace.orbits;
        int norbs = orbs.get_num_orbits();
        Mat<double> V(norbs, norbs, fill::zeros);
        if (n_iter != -1) {
            // cout << monopole.idx_to_ijkl.size() << endl;
            for (int idx=0; idx<monopole.idx_to_ijkl.size(); idx++) {
                int i = monopole.idx_to_ijkl[idx][0];
                int j = monopole.idx_to_ijkl[idx][1];
                int k = monopole.idx_to_ijkl[idx][2];
                int l = monopole.idx_to_ijkl[idx][3];
                V(i,k) += monopole.v2[idx] * rho(j,l);
                // cout << monopole.v2[idx] << endl;
                // printf("  %-3d  %-3d  %-3d  %-3d \n", i, j, k, l);
            
            }
        
            for (int i=0; i<norbs; i++) {
                for (int j=0; j<i+1; j++) {
                    V(j,i) = V(i,j);
                }
            }
        }
        F = Ham.one + V;
        // V.print();
        // cout << endl;
        double r1 = 0.0;
        for (int i=0; i<norbs; i++) {
            for (int j=0; j<norbs; j++) {
                r1 += sqrt( pow(( F(i,j) - Fock_old(i,j) ), 2) );
            }
        }
        return r1;
    }

    void HartreeFock::DiagonalizeFock() {
        OneBodySpace one_body_space = modelspace.one;
        for (int ich=0; ich< one_body_space.number_channels; ich++) {
            vector <int> filter_idx = one_body_space.channels[ich];

            Mat<double> Fch, Sch;
            Fch.zeros(size(F));
            Sch.zeros(size(S));

            for (auto i : filter_idx) {
                for (auto j : filter_idx) {
                    Fch(i,j) = F(i,j);
                    Sch(i,j) = S(i,j);
                }
            }

            // Fch.print("Fch");
            // Sch.print("Sch");

            cx_vec eigval;
            cx_mat eigvec_col;
            eig_pair( eigval, eigvec_col, Fch, Sch );
            cx_mat eigvec_row = eigvec_col.t();
            uvec idxs = arma::sort_index(real(eigval), "ascend");
            
            // int idx = idxs(ich);                
            // complex<double> norm = arma::as_scalar(eigvec_row.row(idx) * Sch * eigvec_col.col(idx));
            // eigvec_col.each_col() += -norm;
            
            // eigvec_col.col(idx) = -1*eigvec_col.col(idx)/norm;


            // Normalising eigenvectors
            // for (int i=0; i<eigvec_col.n_rows; i++) {

            //     // int idx = idxs(i);                
            //     complex<double> norm = arma::as_scalar(eigvec_row.row(i) * Sch * eigvec_col.col(i));
            //     eigvec_col.col(i) = -1*eigvec_col.col(i)/norm;
            // }


            SPEs = real(eigval);
            C = real(eigvec_col);

            // C.print("C");
            // SPEs.print("SPEs");

        }
    }

    void HartreeFock::UpdateDensityMatrix() {
        OneBodySpace one_body_space = modelspace.one;
        Orbits orbs = modelspace.orbits;
        int norbs = orbs.get_num_orbits();
        Mat<int> tmp(norbs, norbs, fill::zeros);
        Orbits orbs_tmp = Orbits();
        map<int, int> orbit_idx_to_spe_idx;
        for (int ich=0; ich<one_body_space.number_channels; ich++) {
            // vector <int> filter_idx = one_body_space.channels[ich];
            // vec filter_idx1 = conv_to<vec>::from(filter_idx);

            vec SPEsCh = SPEs;

            for (int i=0; i<one_body_space.channels[ich].size(); i++) {
                int idx_spe = one_body_space.channels[ich][i];
                if (orbs.relativistic) {
                    if (i < floor(SPEsCh.n_rows/2) ) {
                        n = i;
                        e =-1;
                    } else {
                        n = i - floor(SPEsCh.n_rows/2);
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
            int occ = holes[hole_idx.first];
            int spe_idx = orbit_idx_to_spe_idx[hole_idx.first];
            tmp(spe_idx, spe_idx) = occ;
        }
        rho = C * tmp * C.t();
    }

    void HartreeFock::_print_status(int n_iter, bool detail) {
        cout << "\n";
        printf("nth iteration: %d, HF energy: %3f ", n_iter, En);
        if (detail) {
            cout << "\n" << "\n";
            F.print("Fock Matrix");
            cout << "\n";
            SPEs.t().print("SPEs");
            cout << "\n";
            C.print("Coeffs");
            cout << "\n";
            rho.print("Density Matrix");
            cout << "\n";
        }
    }
