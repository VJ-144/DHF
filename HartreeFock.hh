#ifndef HartreeFock_h
#define HartreeFock_h 1
#include <iostream>
#include <armadillo>
#include <cmath>

#include "Orbits.hh"
#include "Operator.hh"
#include "ModelSpace.hh"
#include "TwoBodySpace.hh"
#include "TwoBodyOperator.hh"

class Monopole
{
  public:
    Operator * H;
    ModelSpace * modelspace;
    std::vector<std::array<int,4>> idx_to_ijkl;
    std::vector<double> vmon2;

    // Constructor
    ~Monopole();
    Monopole();
    Monopole(Operator&);
    void Print();
};

class HartreeFock
{

  public:
    Operator& H;
    ModelSpace * modelspace;
    Monopole monopole;
    arma::mat C;  // (CB|HF)
    arma::mat rho;
    arma::mat F;
    arma::mat V;
    arma::mat S;
    arma::vec SPEs;
    arma::uvec ElectronStates, PositronStates;
    double r;
    double E1, E2, EHF;

    // Constructor
    ~HartreeFock();
    HartreeFock(Operator&);

    void CalcEnergy();
    double UpdateFock(int n_itr = -1);
    void DiagonalizeFock();
    void UpdateDensityMatrix();
    void Solve();
    void PrintStatus(int n_itr);
};
#endif
