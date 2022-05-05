#include<armadillo>

#include "TwoBodyOperator.h"
#include "TwoBodySpace.h"
#include "ModelSpace.h"
#include "Orbits.h"


using namespace arma;

    // Constructor
    TwoBodyOperatorChannel::TwoBodyOperatorChannel(TwoBodyChannel chbra1, TwoBodyChannel chket1){
        chbra = chbra1;
        chket = chket1;

        MEs.zeros(chbra.get_number_states(), chket.get_number_states());
    }

    // Constructor
    TwoBodyOperator::TwoBodyOperator(TwoBodySpace two_body_space1, int rankJ1, int rankP1, int rankZ1){
        two_body_space = two_body_space1;
        rankJ = rankJ1;
        rankP = rankP1;
        rankZ = rankZ1;
        Channels;

        vector<int> vecSize;
        for (int i=0; i<two_body_space.get_number_channels(); i++) {vecSize.push_back(i);}

        for (int ichbra=0; ichbra<vecSize.size(); ichbra++){
            for (int ichket=0; ichket<vecSize.size(); ichket++){
                if (ichbra < ichket) {continue;}
                chbra = two_body_space.get_channel(ichbra);
                chket = two_body_space.get_channel(ichket);
                if (! (abs(chbra.J-chket.J) <= rankJ && rankJ <= chbra.J+chket.J) ) {continue;}
                if (chbra.P * chket.P * rankP == -1) {continue;}
                if (abs(chbra.Z - chket.Z) != rankZ) {continue;}
                Channels[{ichbra, ichket}] = TwoBodyOperatorChannel(chbra, chket);
            }
        }
    }

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//              End of Constructor Overloading

    int TwoBodyOperatorChannel::get_2bme(int idxbra, int idxket) {
        return MEs(idxbra, idxket);
    }

    int TwoBodyOperatorChannel::get_2bme_orbit_indices(int a, int b, int c, int d) {
        vector<int> pairAB, pairCD;
        pairAB = {a, b};
        pairCD = {c, d};

        int idxbra = chbra.index_from_indices[pairAB];
        int idxket = chket.index_from_indices[pairCD];

        int phase = chbra.phase_from_indices[pairAB] * chket.phase_from_indices[pairCD];

        return get_2bme(idxbra, idxket) * phase;
    }

    int TwoBodyOperatorChannel::get_2bme_orbits(Orbit oa, Orbit ob, Orbit oc, Orbit od) {
        Orbits orbs = chket.orbits;

        int a = orbs.get_orbit_index_from_orbit(oa);
        int b = orbs.get_orbit_index_from_orbit(ob);
        int c = orbs.get_orbit_index_from_orbit(oc);
        int d = orbs.get_orbit_index_from_orbit(od);

        return get_2bme_orbit_indices(a, b, c, d);
    }

    void TwoBodyOperatorChannel::set_2bme(int idxbra, int idxket, int v) {
        MEs(idxbra, idxket) = v;
    }

    void TwoBodyOperatorChannel::set_2bme_orbit_indices(int a, int b, int c, int d, int v) {
        vector<int> pairAB, pairCD;
        pairAB = {a,b};
        pairCD = {c,d};

        int idxbra = chbra.index_from_indices[pairAB];
        int idxket = chket.index_from_indices[pairCD];
        int phase = chbra.phase_from_indices[pairAB] * chket.phase_from_indices[pairCD];
        set_2bme(idxbra, idxket, phase*v);
    }

    void TwoBodyOperatorChannel::set_2bme_orbits(Orbit oa, Orbit ob, Orbit oc, Orbit od, int v) {
        Orbits orbs = chket.orbits;

        int a = orbs.get_orbit_index_from_orbit(oa);
        int b = orbs.get_orbit_index_from_orbit(ob);
        int c = orbs.get_orbit_index_from_orbit(oc);
        int d = orbs.get_orbit_index_from_orbit(od);

        set_2bme_orbit_indices(a, b, c, d, v);
    }


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//      TwoBodyOperator Class Methods

    int TwoBodyOperator::get_2bme(int ichbra, int ichket, int idxbra, int idxket) {
        if (ichbra < ichket) {return Channels[{ichket,ichbra}].get_2bme(idxket,idxbra);}
        else {return Channels[{ichbra,ichket}].get_2bme(idxbra,idxket);}
    } 

    int TwoBodyOperator::get_2bme_orbit_indices(int ichbra, int ichket, int a, int b, int c, int d) {
        if (ichbra < ichket) {return Channels[{ichket,ichbra}].get_2bme_orbit_indices(c,d,a,b);}
        else {return Channels[{ichbra,ichket}].get_2bme_orbit_indices(a,b,c,d);}
    }

    int TwoBodyOperator::get_2bme_orbits(int ichbra, int ichket, Orbit oa, Orbit ob, Orbit oc, Orbit od) {
        if (ichbra < ichket) {return Channels[{ichket,ichbra}].get_2bme_orbits(oc,od,oa,ob);}
        else {return Channels[{ichbra,ichket}].get_2bme_orbits(oa,ob,oc,od);}
    }

    int TwoBodyOperator::get_2bme_orbitsJ(Orbit oa, Orbit ob, Orbit oc, Orbit od, int Jab, int Jcd) {
        
        int Pab = pow( -1, (oa.l + ob.l));
        int Pcd = pow( -1, (oc.l + od.l));
        int Zab = floor( (oa.e + ob.e)/2);
        int Zcd = floor( (oc.e + od.e)/2);

        vector<int> JPZ_ab = {Jab, Pab, Zab};
        vector<int> JPZ_cd = {Jcd, Pcd, Zcd};

        int ichbra = two_body_space.get_index(JPZ_ab);
        int ichket = two_body_space.get_index(JPZ_cd);

        return get_2bme_orbits(ichbra,ichket,oa,ob,oc,od);
    }

    int TwoBodyOperator::get_2bme_orbit_indicesJ(int a, int b, int c, int d, int Jab, int Jcd) {
        Orbits orbs = two_body_space.orbits;
        Orbit oa = orbs.get_orbit(a);
        Orbit ob = orbs.get_orbit(b);
        Orbit oc = orbs.get_orbit(c);
        Orbit od = orbs.get_orbit(d);
        return get_2bme_orbitsJ(oa, ob, oc, od, Jab, Jcd);
    }

    void TwoBodyOperator::set_2bme(int ichbra, int ichket, int idxbra, int idxket, int v) {
        if (ichbra < ichket) {Channels[{ichket,ichbra}].set_2bme(idxket,idxbra,v);}
        else {Channels[{ichbra,ichket}].set_2bme(idxbra,idxket,v);}
    }

    void TwoBodyOperator::set_2bme_orbit_indices(int ichbra, int ichket, int a, int b, int c, int d, int v) {
        if (ichbra < ichket) {Channels[{ichket,ichbra}].set_2bme_orbit_indices(c,d,a,b,v);}
        else {Channels[{ichbra,ichket}].set_2bme_orbit_indices(a,b,c,d,v);}
    }

    void TwoBodyOperator::set_2bme_orbits(int ichbra, int ichket, Orbit oa, Orbit ob, Orbit oc, Orbit od, int v) {
        if (ichbra < ichket) {Channels[{ichket,ichbra}].set_2bme_orbits(oc,od,oa,ob,v);}
        else {Channels[{ichbra,ichket}].set_2bme_orbits(oa,ob,oc,od,v);}
    }

    void TwoBodyOperator::set_2bme_orbitsJ(Orbit oa, Orbit ob, Orbit oc, Orbit od, int Jab, int Jcd, int v) {
        
        int Pab = pow( -1, (oa.l + ob.l) );
        int Pcd = pow( -1, (oc.l + od.l) );
        int Zab = floor( (oa.e + ob.e)/2 );
        int Zcd = floor( (oc.e + od.e)/2 );

        vector<int> JPZ_ab = {Jab, Pab, Zab};
        vector<int> JPZ_cd = {Jcd, Pcd, Zcd};

        int ichbra = two_body_space.get_index(JPZ_ab);
        int ichket = two_body_space.get_index(JPZ_cd);

        set_2bme_orbits(ichbra, ichket, oa, ob, oc, od, v);
    }

    void TwoBodyOperator::set_2bme_orbit_indicesJ(int a, int b, int c, int d, int Jab, int Jcd, int v) {
        Orbits orbs = two_body_space.orbits;

        Orbit oa = orbs.get_orbit(a);
        Orbit ob = orbs.get_orbit(b);
        Orbit oc = orbs.get_orbit(c);
        Orbit od = orbs.get_orbit(d);

        return set_2bme_orbitsJ(oa, ob, oc, od, Jab, Jcd, v);
    }

    void TwoBodyOperator::print_two_body_operator() {
        // Note sure if this is nessesary here : two_body_space;
   
        vector<int> NumTwoBodyChannel;
        for (int i=0; i<two_body_space.get_number_channels(); i++) {NumTwoBodyChannel.push_back(i);}

        for (int ichbra=0; ichbra<NumTwoBodyChannel.size(); ichbra++){
            for (int ichket=0; ichket<NumTwoBodyChannel.size(); ichket++){
                
                if (ichbra < ichket) {continue;}
                chbra = two_body_space.get_channel(ichbra);
                chket = two_body_space.get_channel(ichket);

                if (!(abs(chbra.J-chket.J) <= rankJ && rankJ <= chbra.J+chket.J) ) {continue;}
                if (chbra.P * chket.P * rankP == -1) {continue;}
                if (abs(chbra.Z - chket.Z) != rankZ) {continue;}

                // vector<int> chbraNumStates, chketNumStates;
                // for (int i=0; i<chbra.get_number_states(); i++) {chbraNumStates.push_back(i);}
                // for (int i=0; i<chket.get_number_states(); i++) {chketNumStates.push_back(i);}

                for (int idxbra=0; idxbra<chbra.get_number_states(); idxbra++){
                    for (int idxket=0; idxket<chket.get_number_states(); idxket++){

                        int a = chbra.get_indices(idxbra)[0];
                        int b = chbra.get_indices(idxbra)[1];
                        int c = chket.get_indices(idxket)[0];
                        int d = chket.get_indices(idxket)[1];

                        cout << fixed << setw(4) << a << setw(4) << b << setw(4) << c << setw(4) << d << setw(4) << chbra.J << setw(4) << chket.J << "    " << Channels[{ichbra,ichket}].MEs(idxbra,idxket) << endl; 
                    }
                }
            }
        }
    }




