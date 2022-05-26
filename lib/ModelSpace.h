#pragma once
#include "../lib/Orbits.h"
#include "../lib/TwoBodySpace.h"

#include <set>

using namespace std;

class OneBodySpace : public TwoBodySpace{

    public:

        vector<int> kappas;
        Orbits orbits;
        vector<vector<int>> channels;
        map <int, int> orbit_index_to_channel_index;

        OneBodySpace(Orbits orbits1=Orbits());
};


class ModelSpace : public OneBodySpace{

    public:
        OneBodySpace one;
        TwoBodySpace two;

        int Z;
        int Ne;
        double zeta;
        double c = 137.035999084;

        ModelSpace(int Ne1, int Z1, double zeta1, Orbits orbs1);
        ModelSpace(int N1e=1, int Z1=1, double zeta1=1);

        // Function Declerations
        void set_model_space_from_orbits(Orbits orbs);
        void set_model_space_from_truncations(int nmax, int lmax);

};
