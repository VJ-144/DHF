#pragma once
#include <iostream>
#include <armadillo>
#include <cmath>

#include "../lib/Orbits.h"
#include "../lib/Operator.h"
#include "../lib/ModelSpace.h"
#include "../lib/TwoBodySpace.h"
#include "../lib/TwoBodyOperator.h"

using namespace std;
using namespace arma;

class Monopole : public Operator {

    public:
        Operator Ham;
        ModelSpace modelspace;
        vector < vector<int> > idx_to_ijkl;
        vector <double> v2;

        // Constructor
        Monopole(){};
        Monopole(Operator Ham1);

        void set_monopole2();
        void print_monopole2();
};

class HartreeFock : public Monopole {

    public:
        Operator Ham;
        map <int, double> holes;
        ModelSpace modelspace;
        Monopole monopole;
        Mat<double> C, rho, F, V, S;
        vec SPEs;
        double r;
        double En;

        // Constructor
        HartreeFock(){};
        HartreeFock(Operator Ham1, map <int, double> holes1);

        void solve();
        void CalcEnergy();
        double UpdateFock(int n_itr = -1);
        void DiagonalizeFock();
        void UpdateDensityMatrix();
        void _print_status(int n_itr, bool detail=true);
};
