#ifndef Operator_h
#define Operator_h 1
#include <armadillo>
#include "Orbits.hh"
#include "ModelSpace.hh"
#include "TwoBodySpace.hh"
#include "TwoBodyOperator.hh"

class Operator
{
  public:
    int rankJ, rankP;
    ModelSpace * modelspace;
    double ZeroBody;
    arma::mat OneBody, S;
    TwoBodyOperator TwoBody;
    bool orthonormalized = false;

    // Constructor
    ~Operator();
    Operator();
    Operator(ModelSpace& ms, int rankJ=0, int rankP=1);
    Operator( const Operator& rhs); ///< Copy constructor
    Operator( Operator&&);

    Operator& operator=(Operator&& rhs);
    Operator& operator=( const Operator& rhs);

    Operator& operator+=( const Operator& rhs);
    Operator operator+( const Operator& rhs) const;
    Operator& operator+=( const double& rhs);
    Operator operator+( const double& rhs) const;
    Operator& operator-=( const Operator& rhs);
    Operator operator-( const Operator& rhs) const;
    Operator operator-( ) const;
    Operator& operator-=( const double& rhs);
    Operator operator-( const double& rhs) const;
    Operator& operator*=( const double rhs);
    Operator operator*( const double rhs) const;
    Operator& operator/=( const double rhs);
    Operator operator/( const double rhs) const;


    void SetDiracCoulombHamiltonian(bool OneBody=true, bool TwoBody=true);
    void SetCoulombHamiltonian(bool OneBody=true, bool TwoBody=true);
    void SetOneBodyDiracHamiltonian();
    void SetOneBodyHamiltonian();
    void SetTwoBodyCoulombInteraction();
    void OrthoNormalize();
    double Get2BME(int ichbra, int ichket, int i, int j, int k, int l) {return TwoBody.Get2BME(ichbra,ichket,i,j,k,l);};
    double Get2BME(int ichbra, int ichket, Orbit& o1, Orbit& o2, Orbit& o3, Orbit& o4) {return TwoBody.Get2BME(ichbra,ichket,o1,o2,o3,o4);};
    double Get2BME_J(int Jbra, int Jket, int i, int j, int k, int l) {return TwoBody.Get2BME_J(Jbra,Jket,i,j,k,l);};
    double Get2BME_J(int Jbra, int Jket, Orbit& o1, Orbit& o2, Orbit& o3, Orbit& o4) {return TwoBody.Get2BME_J(Jbra,Jket,o1,o2,o3,o4);};
    ModelSpace * GetModelSpace() const {return modelspace;};
    arma::mat EmbedBasisTrans2(arma::mat, TwoBodyChannel&);
    void Print();
    Operator DoNormalOrdering(int);
    Operator DoNormalOrdering() {return DoNormalOrdering(1);};
    Operator UndoNormalOrdering() {return DoNormalOrdering(-1);};
    arma::vec DiagonalizeOneBody();
};
#endif
