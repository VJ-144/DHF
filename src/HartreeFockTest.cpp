#include <iostream>
#include <algorithm> 
#include <armadillo>
#include <cmath>

#include <algorithm>
#include <vector>
#include <iterator>

#include "../lib/TwoBodyOperator.h"
#include "../lib/TwoBodySpace.h"
#include "../lib/HartreeFock.h"
#include "../lib/ModelSpace.h"
#include "../lib/Operator.h"
#include "../lib/Orbits.h"


using namespace std;
using namespace arma;


int main() {
    Orbits orbs = Orbits("LSpinor", false, true);
    // vector <double> zetainv_list = {0.5, 1, 2, 4};
    vector <double> zetainv_list = {2};
    vector <double> Nmax_list = {1};    
  
    for (auto Nmax : Nmax_list) {
        vector <double> x, y;
        x = {};
        y = {};
        for (auto zeta_inv : zetainv_list) {
            orbs.set_orbits(Nmax, 0);
            map <int, double> holes;
            /*
            He2
            */
            int i = orbs.get_orbit_index({0,0,1,1}); holes[i] = 1.0;
            /*
            Ne10
            */
            // int i = orbs.get_orbit_index({0,0,1,1}); holes[i] = 1.0;
            // int i = orbs.get_orbit_index({1,0,1,1}); holes[i] = 1.0;
            // int i = orbs.get_orbit_index({1,1,1,1}); holes[i] = 1.0;
            // int i = orbs.get_orbit_index({1,1,3,1}); holes[i] = 1.0;

            orbs.print_orbits();

            ModelSpace ms = ModelSpace(1, 2, 1/zeta_inv);
            ms.set_model_space_from_orbits(orbs);
            Operator Ham = Operator(ms, 0, 1, 0);
            // Ham.set_hamiltonian(true, false);
            Ham.set_hamiltonian(true, true);
            Ham.print_operator();
            Ham.orthogonalize(true);

            // Ham.S.print("S");

            // Ham.one.print();
            // Monopole mon = Monopole(Ham);
            // the issue is initialising my operator in the header file of HartreeFock.h

            // HartreeFock HF = HartreeFock(Ham, holes);
            // HF.solve();
            // x.push_back(zeta_inv);
            // y.push_back(HF.En);
        }
        // these are vectors, cannot print this way
        // cout << x << endl;
        // cout << y << endl;
    }
    return 0;
}


// Run Test File Command
// g++ Orbits.cpp TwoBodySpace.cpp TwoBodyOperator.cpp ModelSpace.cpp Operator.cpp HartreeFock.cpp HartreeFockTest.cpp -larmadillo -lgsl -o HartreeFock