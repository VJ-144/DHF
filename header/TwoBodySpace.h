#pragma once
#include "../header/Orbits.h"

class TwoBodyChannel : public Orbits {

    public:
        // Member variables
        Orbits orbits;
        vector< vector<Orbit> > comb;        
        int J, P, Z, number_states;
        vector<int> orbit1_index, orbit2_index;
        map<vector<int>, int> phase_from_indices, index_from_indices;

        // Constructors
        // TwoBodyChannel(){};
        TwoBodyChannel();
        TwoBodyChannel(int J1, int P1, int Z1, Orbits orbits1, int e2max1);
        
        // Function Declerations
        vector<int> get_JPZ();
        int get_number_states();         
        void _set_two_body_channel();
        vector<int> get_indices(int idx);
        vector<Orbit> get_orbits(int idx);
        bool _triag(int J1, int J2, int J3);

        // Combination Functions
        void combinations_r(vector<Orbit> elems, int req_len);
        vector<vector<Orbit>> combinCalc(vector<Orbit> comboVec, int size);        
        void combinations_r_recursive(vector<Orbit> elems, int req_len, vector<int> pos, int depth, int margin);
};


class TwoBodySpace : public TwoBodyChannel {

    public:
        // Member variables
        Orbits orbits;
        int e2max, number_channels;
        vector<TwoBodyChannel> channels;        
        map < vector<int>, int > index_from_JPZ;

        // Constructors
        TwoBodySpace(){};
        TwoBodySpace(Orbits orbits1);
        TwoBodySpace(Orbits orbits1, int e2max1);

        // Function Declerations
        void print_channels();        
        int get_number_channels();
        int get_index(vector<int> JPZ);
        TwoBodyChannel get_channel(int idx);
        TwoBodyChannel get_channel_from_JPZ(vector<int> JPZ);
};




