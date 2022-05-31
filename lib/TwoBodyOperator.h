#pragma once
#include <iostream>
#include <armadillo>

#include "../lib/Orbits.h"
#include "../lib/TwoBodySpace.h"

using namespace arma;

class TwoBodyOperatorChannel : public TwoBodySpace {

    public:
        TwoBodyChannel chbra, chket;
        Mat<double> MEs;

        // Constructor
        TwoBodyOperatorChannel(){};
        // TwoBodyOperatorChannel(TwoBodyChannel chbra1=TwoBodyChannel(), TwoBodyChannel chket1=TwoBodyChannel());
        TwoBodyOperatorChannel(TwoBodyChannel chbra1, TwoBodyChannel chket1);

        int get_2bme(int idxbra, int idxket);
        int get_2bme_orbit_indices(int a, int b, int c, int d);
        int get_2bme_orbits(Orbit oa, Orbit ob, Orbit oc, Orbit od);
        void set_2bme(int idxbra, int idxket, double v);
        void set_2bme_orbit_indices(int a, int b, int c, int d, double v);
        void set_2bme_orbits(Orbit oa, Orbit ob, Orbit oc, Orbit od, double v);
};


class TwoBodyOperator : public TwoBodyOperatorChannel {

    public:
        int rankJ, rankP, rankZ;
        TwoBodySpace two_body_space;
        map <vector<int>, TwoBodyOperatorChannel > Channels;

        // Constructor
        TwoBodyOperator(){};
        TwoBodyOperator(TwoBodySpace two_body_space1, int rankJ1=0, int rankP1=1, int rankZ1=0);

        int get_2bme(int ichbra, int ichket, int idxbra, int idxket);
        int get_2bme_orbit_indices(int ichbra, int ichket, int a, int b, int c, int d);
        int get_2bme_orbits(int ichbra, int ichket, Orbit oa, Orbit ob, Orbit oc, Orbit od);
        int get_2bme_orbitsJ(Orbit oa, Orbit ob, Orbit oc, Orbit od, int Jab, int Jcd);
        int get_2bme_orbit_indicesJ(int a, int b, int c, int d, int Jab, int Jcd);

        void set_2bme(int ichbra, int ichket, int idxbra, int idxket, double v);
        void set_2bme_orbit_indices(int ichbra, int ichket, int a, int b, int c, int d, double v);
        void set_2bme_orbits(int ichbra, int ichket, Orbit oa, Orbit ob, Orbit oc, Orbit od, double v);
        void set_2bme_orbitsJ(Orbit oa, Orbit ob, Orbit oc, Orbit od, int Jab, int Jcd, double v);
        void set_2bme_orbit_indicesJ(int a, int b, int c, int d, int Jab, int Jcd, double v);
        void print_two_body_operator();
};
