#include "../header/TwoBodySpace.h"
#include "../header/Orbits.h"

    // Constructor for full parameters
    TwoBodyChannel::TwoBodyChannel(int J1, int P1, int Z1, Orbits orbits1, int e2max1) {

        J = J1;
        P = P1;
        Z = Z1;
        // e2max = e2max1;       
        orbits = orbits1;

        // Constructor Variables 
        orbit1_index;
        orbit2_index;
        phase_from_indices;
        index_from_indices;
        number_states = 0;

        _set_two_body_channel();

    }

    // General Constructor
    TwoBodyChannel::TwoBodyChannel() {

        J;
        P;
        Z;
        // e2max;        
        orbits;

        // Constructor Variables 
        orbit1_index;
        orbit2_index;
        phase_from_indices;
        index_from_indices;
        number_states = 0;

    }

    // Constructor 2 parameters
    TwoBodySpace::TwoBodySpace(Orbits orbits1, int e2max1) {

        orbits = orbits1;
        e2max = e2max1;
        index_from_JPZ;
        channels;
        number_channels = 0;

        for ( int j=0; j<e2max+2; j++) {
            for (int p : {-1,1}) {
                for ( int z : {-1,0,1} ) {
                    TwoBodyChannel channel = TwoBodyChannel(j, p, z, orbits, e2max);
                    if ( channel.get_number_states() == 0 ) { continue; }
                    channels.push_back(channel);
                    int idx = channels.size() - 1;
                    index_from_JPZ[{J,P,Z}] = idx;
                }
            }
        number_channels = channels.size();    
        }
    }

    // Constructor 1 parameter
    TwoBodySpace::TwoBodySpace(Orbits orbits1) {

        orbits = orbits1;
        index_from_JPZ;
        channels;
        number_channels = 0;
        e2max;

        e2max = 2*orbits.emax; 


        for ( int j=0; j<e2max+2; j++) {
            for (int p : {-1,1}) {
                for ( int z : {-1,0,1} ) {
                    TwoBodyChannel channel= TwoBodyChannel(j, p, z, orbits, e2max);
                    if ( channel.get_number_states() == 0 ) { continue; }
                    channels.push_back(channel);
                    int idx = channels.size() - 1;
                    index_from_JPZ[{J,P,Z}] = idx;
                }
            }
        number_channels = channels.size();    
        }
    }


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//                  End of Constructors

    void TwoBodyChannel::_set_two_body_channel(){
        Orbits orbs = orbits;
        
        vector< vector<Orbit> > OaOb = combinCalc(orbs.orbits, 2);

        for (int i = 0; i < OaOb.size(); i++) {

            Orbit oa = OaOb[i][0];
            Orbit ob = OaOb[i][1];

            int ia = orbs.get_orbit_index_from_orbit(oa);
            int ib = orbs.get_orbit_index_from_orbit(ob);

            if (ia == ib && J%2==1) {continue;}
            if ( (oa.e + ob.e) != 2*Z) {continue;}
            if ( pow((-1), (oa.l + ob.l)) != P) {continue;}
            if ( _triag(oa.j, ob.j, 2*J) ) {continue;}

            orbit1_index.push_back(ia);
            orbit2_index.push_back(ib);

            int idx = orbit1_index.size() - 1;
            
            vector<int> Pair, InvPair;
            Pair = {ia, ib};
            InvPair = {ib, ia};

            index_from_indices[Pair] = idx;  
            index_from_indices[InvPair] = idx;
            phase_from_indices[Pair] = 1;
            phase_from_indices[InvPair] = pow( -(-1), floor( (oa.j+ob.j)/2 ) - J );
        }
        number_states = orbit1_index.size();
    }

    int TwoBodyChannel::get_number_states() {
        return number_states;
    }

    vector<int> TwoBodyChannel::get_indices(int idx) {
        vector<int> orbit_indices = {orbit1_index[idx], orbit2_index[idx]};
        return orbit_indices;
    }

    vector<Orbit> TwoBodyChannel::get_orbits(int idx) {
        vector <int> orbits_indices;  
        int ia = get_indices(idx)[0];
        int ib = get_indices(idx)[1];

        vector<Orbit> orbits_list= { orbits.get_orbit(ia), orbits.get_orbit(ib) };
        return orbits_list;
    }

    vector<int> TwoBodyChannel::get_JPZ() {
        vector<int> JPZ = {J, P, Z};
        return JPZ;
    }

    bool TwoBodyChannel::_triag(int J1, int J2, int J3) {
        bool b = true;
        if ( abs(J1-J2) <= J3 && J3 <= J1+J2 ) { b = false; }
        return b;
    }

    int TwoBodySpace::get_number_channels() {
        return number_channels;
    }

    int TwoBodySpace::get_index(vector<int> JPZ) {
        return index_from_JPZ[JPZ];
    }


    TwoBodyChannel TwoBodySpace::get_channel(int idx) {
        return channels[idx];
    }


    TwoBodyChannel TwoBodySpace::get_channel_from_JPZ(vector<int> JPZ) {
        return get_channel( get_index(JPZ) );
    }


    void TwoBodySpace::print_channels() {
        cout << " Two-body channels list " << endl;
        cout << " J, par, Z, # of states " << endl;
        for ( auto channel : channels ) { 
            
            vector <int> JPZ = channel.get_JPZ();
            cout << fixed << setw(2) << JPZ[0] << setw(4) << JPZ[1] << setw(4) << JPZ[2] << setw(4) << channel.get_number_states() << endl;
         }

    }

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //              Combination Functions

    void TwoBodyChannel::combinations_r_recursive(vector<Orbit> elems, int req_len, vector<int> pos, int depth, int margin) {
	vector<Orbit> pair;
	// Have we selected the number of required elements?
	if (depth >= req_len) {
		for (int ii = 0; ii < pos.size(); ++ii)
			// cout << elems[pos[ii]];
			pair.push_back(elems[pos[ii]]);
			comb.push_back(pair);
		// cout << endl;
		return;
	}

	// Try to select new elements to the right of the last selected one.
	for (int ii = margin; ii < elems.size(); ++ii) {
		pos[depth] = ii;
		// pair.push_back(pos[depth]);
		combinations_r_recursive(elems, req_len, pos, depth + 1, ii);
		
	}
	return;

    }

    void TwoBodyChannel::combinations_r(vector<Orbit> elems, int req_len) {
        
        vector<int> positions(req_len, 0);
        combinations_r_recursive(elems, req_len, positions, 0, 0);
    }

    vector<vector<Orbit>> TwoBodyChannel::combinCalc(vector<Orbit> comboVec, int size) {

	int num_elements = comboVec.size();
	vector<Orbit> elements(num_elements);
	copy(comboVec.begin(), comboVec.begin() + num_elements, elements.begin());
	combinations_r(elements, size);

    return comb;

    }
