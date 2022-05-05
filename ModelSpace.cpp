
#include "Orbits.h"
#include "ModelSpace.h"
#include "TwoBodySpace.h"

    OneBodySpace::OneBodySpace(Orbits orbits){

        orbits = orbits;
        kappas;
        channels;
        orbit_index_to_channel_index;

        for (auto o : orbits.orbits){ kappas.push_back(o.k); }

        for (int i=0; i < kappas.size(); i++){
            vector <int> idxs;
            for (auto o : orbits.orbits){
                if (o.k != kappas[i]) {continue;}
                idxs.push_back(orbits.get_orbit_index_from_orbit(o));
                orbit_index_to_channel_index[orbits.get_orbit_index_from_orbit(o)] = i;
            }
            channels.push_back(idxs);
        }
        number_channels = channels.size();
    }

    ModelSpace::ModelSpace(int Ne, int Z, int zeta, Orbits orbs){
        c;
        Z=1;
        Ne=1;
        zeta=1;

        set_model_space_from_orbits(orbs);
    }

    ModelSpace::ModelSpace(int Ne, int Z, int zeta){
        c;
        Z;
        Ne;
        zeta;

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
