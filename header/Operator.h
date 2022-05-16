#pragma once
#include <iostream>
#include <armadillo>
#include <cmath>
#include <gsl/gsl_integration.h>

#include "../header/Orbits.h"
#include "../header/ModelSpace.h"
#include "../header/TwoBodySpace.h"
#include "../header/TwoBodyOperator.h"

using namespace arma;
using namespace std;

class Operator : public TwoBodyOperator {

    public:
        int rankJ, rankP, rankZ;
        ModelSpace modelspace;
        Mat<double> S, one;
        TwoBodyOperator two;

        // Constructor
        Operator(ModelSpace modelspace1 = ModelSpace(), int rankJ1=0, int rankP1=1, int rankZ1=0);

        void set_hamiltonian(bool one_body=true, bool two_body=true);
        void set_one_body_ham();
        void set_two_body_ham();
        // double g(double *r, size_t dim, void *params);
        double _coulomb(Orbit& oa, Orbit& ob, Orbit& oc, Orbit& od, int J);
        void orthogonalize(bool scalar=true);
        void print_operator();

};
