#pragma once
#include <iostream>
#include <armadillo>
#include <cmath>
#include <gsl/gsl_integration.h>

#include "../header/Orbits.h"
#include "../header/Operator.h"
#include "../header/ModelSpace.h"
#include "../header/TwoBodySpace.h"
#include "../header/TwoBodyOperator.h"

using namespace std;
using namespace arma;

class Monopole : public Operator {

    public:
        Operator Ham;
        ModelSpace modelspace;
        vector < vector<int> > idx_to_ijkl;
        vector <double> v2;

        // Constructor
        Monopole(Operator Ham1 = Operator()) {
            Ham=Ham1;
            modelspace = Ham.modelspace;
            idx_to_ijkl = {};
            v2 = {};
        };

        void set_monopole2();
        void print_monopole2();
};

class HartreeFock : public Monopole {

    public:
        Operator Ham;
        map <int, double> holes;
        ModelSpace modelspace;
        Monopole monopole;
        // int norbs;
        // Orbits orbs;
        Mat<double> C, rho, F, V, S;
        dvec SPEs;
        // double r;
        double En;

        // Constructor
        HartreeFock(Operator Ham1, map <int, double> holes1) {
            En; // Not sure if this belongs here
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

        void solve();
        void CalcEnergy();
        double UpdateFock(int n_itr = -1);
        void DiagonalizeFock();
        void UpdateDensityMatrix();
        void _print_status(int n_itr, bool detail);
};
