
#include "../header/Orbits.h"
#include "../header/ModelSpace.h"
#include "../header/TwoBodySpace.h"

#include <set>
#include <armadillo>
using namespace arma;
using namespace std;

    OneBodySpace::OneBodySpace(Orbits orbits1){

        orbits = orbits1;
        kappas;
        channels;
        orbit_index_to_channel_index;

        for (Orbit o : orbits.orbits){ kappas.push_back(o.k); }
        kappas.erase ( unique(kappas.begin(), kappas.end()), kappas.end() );

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
    };

    ModelSpace::ModelSpace(int Ne1, int Z1, double zeta1, Orbits orbs){
        c;
        Z=Z1;
        Ne=Ne1;
        zeta=zeta1;

        set_model_space_from_orbits(orbs);
    }

    ModelSpace::ModelSpace(int Ne1, int Z1, double zeta1){
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
