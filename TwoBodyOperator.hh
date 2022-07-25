#ifndef TwoBodyOperator_h
#define TwoBodyOperator_h 1
#include <armadillo>
#include "Orbits.hh"
#include "TwoBodySpace.hh"
#include "ModelSpace.hh"

class TwoBodyOperatorChannel
{
  public:
    // Constructor
    ~TwoBodyOperatorChannel();
    TwoBodyOperatorChannel();
    TwoBodyOperatorChannel(TwoBodyChannel&, TwoBodyChannel&);

    TwoBodyChannel * chbra;
    TwoBodyChannel * chket;
    arma::mat MEs;
};

class TwoBodyOperator
{
  public:
    int rankJ, rankP;
    ModelSpace * modelspace;
    std::map<std::array<int,2>,TwoBodyOperatorChannel> Channels;
    int NMesh;
    double Rmax;
    std::string FileNameCoulomb;
    std::string MeshType;

    // Constructor
    ~TwoBodyOperator();
    TwoBodyOperator();
    TwoBodyOperator(ModelSpace& ms, int rankJ=0, int rankP=1);
    TwoBodyOperator& operator*=(const double);
    TwoBodyOperator& operator+=(const TwoBodyOperator&);
    TwoBodyOperator& operator-=(const TwoBodyOperator&);


    double Get2BME(int, int, int, int, int, int);
    double Get2BME(int, int, Orbit&, Orbit&, Orbit&, Orbit&);
    double Get2BME_J(int, int, int, int, int, int);
    double Get2BME_J(int, int, Orbit&, Orbit&, Orbit&, Orbit&);
    void Set2BME(int, int, int, int, int, int, double);
    void Set2BME(int, int, Orbit&, Orbit&, Orbit&, Orbit&, double);
    void Set2BME_J(int, int, int, int, int, int, double);
    void Set2BME_J(int, int, Orbit&, Orbit&, Orbit&, Orbit&, double);
    void Print();

    void SetTwoBodyCoulombTerm();
    void SetNMesh(int num) {NMesh=num;};
    void SetRmax(double d) {Rmax=d;};
    void SetMeshType(std::string str) {MeshType=str;};
    void SetFileCoulomb(std::string str) {FileNameCoulomb=str;};
    void ReadFileInteraction(std::string);
    void WriteFileInteraction(std::string);
    long long int CountInteractionMEs();
    double MECoulomb(Orbit&, Orbit&, Orbit&, Orbit&, int);
    double TestMECoulomb(Orbit&, Orbit&, Orbit&, Orbit&, int);
    void StoreCoulombIntegrals();
    double GetCoulombIntegral(int, int, int, int, int);
};
#endif
