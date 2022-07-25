#ifndef ModelSpace_h
#define ModelSpace_h 1
#include "Orbits.hh"
#include "TwoBodySpace.hh"

#include <set>
#include <unordered_set>

class OneBodySpace
{
  public:
    ~OneBodySpace();
    OneBodySpace();
    OneBodySpace(Orbits&);

    Orbits * orbits;
    std::vector<int> kappas;
    std::vector<int> n_large_channel;
    std::vector<int> n_small_channel;
    std::vector<std::vector<int>> channels;
    std::map <int, int> orbit_index_to_channel_index;
    int number_channels;

    int GetNumberChannels() {return number_channels;};
    int GetKappa(int idx) {return kappas[idx];};
    void PrintSpace();
};


class ModelSpace : public OneBodySpace
{
  public:
    ~ModelSpace();
    ModelSpace(int, double, double, Orbits);
    ModelSpace(std::string, double, Orbits);

    Orbits orbits;
    OneBodySpace one;
    TwoBodySpace two;
    std::map<int,double> hole_occ;
    double Z;
    int Ne;
    double zeta;
    std::set<int> holes;
    std::set<int> particles;
    std::set<int> core;
    std::set<int> valence;
    std::set<int> qspace;
    std::set<int> large_components;
    std::set<int> small_components;
    std::set<int> all_orbits;
    

    int GetElectronNumber() {return Ne;};
    std::map<int,double> AssignHoles(int);
    double GetProtonNumber() {return Z;};
    double GetZeta() {return zeta;};
    void PrintModelSpace(bool, bool);
    void PrintHoleOrbits();
    Orbits& GetOrbits() {return orbits;};
    std::string GetOrbitLabel(int idx) {return orbits.GetOrbitLabel(idx);};
    Orbit& GetOrbit(int idx) {return orbits.GetOrbit(idx);};
    int GetNmax() {return orbits.GetNmax();};
    int GetLmax() {return orbits.GetLmax();};

    int GetKappaMin() {return orbits.GetKappaMin();};
    int GetKappaMax() {return orbits.GetKappaMax();};
    int GetOrbitIndex(int n, int kappa, int ls) {return orbits.GetOrbitIndex(n,kappa,ls);};
    int GetOrbitIndex(Orbit& o) {return orbits.GetOrbitIndex(o);};
    int GetOrbitIndex(int n, int l, int j2, int ls) {return orbits.GetOrbitIndex(n,l,j2,ls);};
    int GetNumberOrbits() {return orbits.GetNumberOrbits();};
    //void UpdateHoles(std::map<int,double> holes) {hole_occ = holes;};
    //void UpdateHoles(arma::vec SPEs) {hole_occ = GetElectronOccupation(SPEs);};
    //std::map<int,double> GetElectronOccupation(arma::vec SPEs);
    void UpdateOccupation(std::map<int,double> holes);
    void UpdateOrbitals(arma::vec);

    OneBodySpace& GetOneBodySpace() {return one;};
    TwoBodySpace& GetTwoBodySpace() {return two;};
    TwoBodyChannel& GetTwoBodyChannel(int idx) {return two.GetChannel(idx);};
    TwoBodyChannel& GetTwoBodyChannel(int J, int Prty) {return two.GetChannel(J, Prty);};
};
#endif
