#include "../lib/TwoBodySpace.h"
#include "../lib/ModelSpace.h"
#include "../lib/Orbits.h"

#include <set>
#include <unordered_set>
#include <bits/stdc++.h>
#include <armadillo>
using namespace arma;
using namespace std;

    OneBodySpace::OneBodySpace(Orbits orbits1){

        orbits = orbits1;
        kappas;
        channels;
        orbit_index_to_channel_index;

        for (Orbit o : orbits.orbits){ kappas.insert(o.k); }

        vector<int> kappas1(kappas.begin(), kappas.end());
        sort(kappas1.begin(), kappas1.end(), greater<int>());
        // sorting kappas into decreasing order for printing purposes - more similar to python code



        // sort( kappas.begin(), kappas.end() );
        // kappas.erase( unique( kappas.begin(), kappas.end() ), kappas.end() );


        // unordered_set<int> s( kappas.begin(), kappas.end() );
        // kappas.assign( s.begin(), s.end() );
        // sort( kappas.begin(), kappas.end(), greater<int>() );

        // vec vect = conv_to< vec >::from(kappas1);
        // vect.t().print("kappas1");

        for (int channel_idx=0; channel_idx<kappas1.size(); channel_idx++){
            vector <int> idxs;
            for (auto o : orbits.orbits){
                if (o.k != kappas1[channel_idx]) {continue;}
                idxs.push_back(orbits.get_orbit_index_from_orbit(o));
                orbit_index_to_channel_index[orbits.get_orbit_index_from_orbit(o)] = channel_idx;
            }
            channels.push_back(idxs);  
        }
        number_channels = channels.size();
    };

    ModelSpace::ModelSpace(int Ne1, double Z1, double zeta1, Orbits orbs){
        c;
        Z=Z1;
        Ne=Ne1;
        zeta=zeta1;

        set_model_space_from_orbits(orbs);
    }

    ModelSpace::ModelSpace(int Ne1, double Z1, double zeta1){
        c;
        Z=Z1;
        Ne=Ne1;
        zeta=zeta1;
    }
    
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//          End of Constructor Overloading 

    // Model Space Class Functions
    void ModelSpace::set_model_space_from_orbits(Orbits orbs){
        orbits = orbs;
        one = OneBodySpace(orbits);
        two = TwoBodySpace(orbits);
    }

    void ModelSpace::set_model_space_from_truncations(int nmax, int lmax){
        orbits = Orbits();
        one = OneBodySpace(orbits);
        orbits.set_orbits(nmax, lmax);
        two = TwoBodySpace(orbits);
    }
