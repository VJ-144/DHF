#include <iostream>
#include <algorithm> 
#include <armadillo>
#include <cmath>

#include <algorithm>
#include <vector>
#include <iterator>
#include <iostream>

#include "../header/TwoBodyOperator.h"
#include "../header/TwoBodySpace.h"
#include "../header/ModelSpace.h"
#include "../header/Operator.h"
#include "../header/Orbits.h"


using namespace std;
using namespace arma;

int main() {

    // TwoBodySpace Test
    // Orbits orbs = Orbits();
    // orbs.set_orbits(1,0);
    // TwoBodySpace two = TwoBodySpace(orbs);
    // two.print_channels();
    
    // Model Space Test
    // Orbits orbs = Orbits("LSpinor", false, true);
    // orbs.set_orbits(2,1);
    // ModelSpace ms = ModelSpace();
    // ms.set_model_space_from_orbits(orbs);

    // TwoBodyOperator Test
    Orbits orbs = Orbits();
    orbs.set_orbits(2,1);
    TwoBodySpace two_body_space =  TwoBodySpace(orbs);
    TwoBodyOperator Op2 = TwoBodyOperator(two_body_space);
    Op2.print_two_body_operator();

    // Operator Test
    // Orbits orbs = Orbits("LSpinor", false, true);
    // orbs.set_orbits(1,1);
    // // orbs.print_orbits();
    // ModelSpace ms = ModelSpace(1,2,1);
    // ms.set_model_space_from_orbits(orbs);
    // Operator Ham = Operator(ms, 0, 1, 0);
    // Ham.set_hamiltonian(true, true);
    // // Ham.orthogonalize();
    // Ham.print_operator();


    return(0);
}


// TwoBodyTest Command
// g++ TwoBodyTest.cpp -lgsl -larmadillo -o TwoBodyTest

// ModelSpace Test Command
// g++ Orbits.cpp TwoBodySpace.cpp ModelSpace.cpp TwoBodyTest.cpp -larmadillo -lgsl -o ModelSpace

// TwoBodySpace Test Command
// g++ Orbits.cpp TwoBodySpace.cpp TwoBodyTest.cpp -larmadillo -lgsl -o TwoBodySpace

// TwoBodyOperator Test Command
// g++ Orbits.cpp TwoBodySpace.cpp TwoBodyOperator.cpp TwoBodyTest.cpp -larmadillo -lgsl -o TwoBodyOperator

// Operator Test Command
// g++ Orbits.cpp TwoBodySpace.cpp TwoBodyOperator.cpp TwoBodyTest.cpp ModelSpace.cpp Operator.cpp -larmadillo -lgsl -lwignerSymbols -o Operator
