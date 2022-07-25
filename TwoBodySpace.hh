#ifndef TwoBodySpace_h
#define TwoBodySpace_h 1
#include <array>
#include <vector>
#include <map>
#include "Orbits.hh"

class TwoBodyChannel
{
  public:
    // Constructors
    ~TwoBodyChannel();
    TwoBodyChannel(int, int, Orbits&);

    // Member variables
    Orbits* orbits;
    int J, Prty, number_states;
    std::vector<int> index1, index2;
    std::map<std::array<int,2>, int> phase, index;

    // Function Declerations
    int GetNumberStates() {return number_states;};
    int GetPhaseFactor(int ip, int iq) {return phase[{ip,iq}];};
    int GetIndex(int ip, int iq) {return index[{ip,iq}];};
    int GetOrbitIndex1(int idx) {return index1[idx];};
    int GetOrbitIndex2(int idx) {return index2[idx];};
    Orbit& GetOrbit1(int idx) {return orbits->GetOrbit(GetOrbitIndex1(idx));};
    Orbit& GetOrbit2(int idx) {return orbits->GetOrbit(GetOrbitIndex2(idx));};
    void PrintChannel();
};


class TwoBodySpace
{
  public:
    ~TwoBodySpace();
    TwoBodySpace();
    TwoBodySpace(Orbits&);

    // Member variables
    Orbits* orbits;
    int number_channels;
    std::vector<TwoBodyChannel> channels;
    std::map< std::array<int,2>, int > index_from_JP;

    // Function Declerations
    int GetNumberChannels() {return number_channels;};
    void PrintSpace() {for(auto tbc: channels) tbc.PrintChannel();};
    int GetChannelIndex(int J, int Prty) {return index_from_JP[{J,Prty}];};
    TwoBodyChannel& GetChannel(int idx) {return channels[idx];};
    TwoBodyChannel& GetChannel(int J, int Prty) {return GetChannel(index_from_JP[{J,Prty}]);};
};
#endif
