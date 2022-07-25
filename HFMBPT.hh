#ifndef HFMBPT_h
#define HFMBPT_h 1

#include "HartreeFock.hh"
class HFMBPT
{
  public:
    Operator& H;
    ModelSpace * modelspace;
    arma::mat C_CB2HF;  // computational basis to HF
    arma::mat C_CB2NAT; // computational basis to NAT
    arma::mat C_HF2NAT; // HF to NAT
    arma::mat C_trans;
    arma::vec Occ;

    ~HFMBPT();
    HFMBPT(Operator &, arma::mat);
    void SetCB2HF(arma::mat mat) {C_CB2HF = mat;};
    void SetTransformationMatrix(std::string);
    void SetTransformationMatrix(arma::mat m) {C_trans=m;};
    Operator TransformBasis(Operator&);
    double GetMP2_Energy(Operator &);
};
#endif

