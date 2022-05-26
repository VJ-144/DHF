#pragma once
#include <iostream>
#include <armadillo>
#include <cmath>
#include <gsl/gsl_integration.h>

#include "../lib/Orbits.h"
#include "../lib/ModelSpace.h"
#include "../lib/TwoBodySpace.h"
#include "../lib/TwoBodyOperator.h"

using namespace arma;
using namespace std;

class Operator : public TwoBodyOperator {

    public:
        int rankJ, rankP, rankZ;
        ModelSpace modelspace;
        Mat<double> S, one;
        TwoBodyOperator two;

        // Constructor
        Operator(){};
        Operator(ModelSpace modelspace1, int rankJ1, int rankP1, int rankZ1);

        void set_hamiltonian(bool one_body=true, bool two_body=true);
        void set_one_body_ham();
        void set_two_body_ham();
        double _coulomb(Orbit& oa, Orbit& ob, Orbit& oc, Orbit& od, int J);
        void orthogonalize(bool scalar=true);
        void print_operator();

};
