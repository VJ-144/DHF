#include <armadillo>
#include <boost/range/irange.hpp>
#include <cmath>
#include <wignerSymbols.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_miser.h>
// #include <wignerSymbols>
// #include "wignerSymbols/wignerSymbols-cpp.h"
// #include<gsl>

#include "../header/TwoBodyOperator.h"
#include "../header/TwoBodySpace.h"
#include "../header/ModelSpace.h"
#include "../header/Operator.h"
#include "../header/Orbits.h"

using namespace WignerSymbols;
using namespace boost;
using namespace arma;
using namespace std;

    // Wrapper to convert lambda with capture to gsl function for integration
    template< typename F >  class gsl_function_pp : public gsl_function {
    public:
        gsl_function_pp(const F& func) : _func(func) {
            function = &gsl_function_pp::invoke;
            params = this;
        }
    private:
        const F& _func;
        static double invoke(double x, void* params) {
            return static_cast<gsl_function_pp*>(params)->_func(x);
        }
    };

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    // Constructor
    Operator::Operator(ModelSpace modelspace1, int rankJ1, int rankP1, int rankZ1) {
        modelspace = modelspace1;
        rankJ = rankJ1;
        rankP = rankP1;
        rankZ = rankZ1;
        Orbits orbs = modelspace.orbits;
        
        S.zeros(orbs.get_num_orbits(), orbs.get_num_orbits());
        one.zeros(orbs.get_num_orbits(), orbs.get_num_orbits());
        TwoBodySpace twoTest = modelspace.two;
        two = TwoBodyOperator(twoTest, rankJ1, rankP1, rankZ1);
    }

    void Operator::set_hamiltonian(bool one_body, bool two_body) {
        if (one_body) {set_one_body_ham();}
        if (two_body) {set_two_body_ham();}
    }

    void Operator::set_one_body_ham() {
        ModelSpace ms = modelspace;
        Orbits orbs = ms.orbits;
        Mat<double> H, S;
        H.zeros(orbs.get_num_orbits(), orbs.get_num_orbits());
        S.zeros(orbs.get_num_orbits(), orbs.get_num_orbits());

        for (int i=0; i<orbs.get_num_orbits(); i++){
            for (int j=0; j<orbs.get_num_orbits(); j++){
                Orbit oi = orbs.get_orbit(i);
                Orbit oj = orbs.get_orbit(j);
                int fac;
                if (oi.e== 1 && oj.e== 1) {fac = 1;}
                if (oi.e==-1 && oj.e==-1) {fac =-1;}
                if (oi.k!=oj.k) {continue;}
                if (orbs.relativistic) {
                    H(i,j) = ME_Kinetic(oi,oj,ms.zeta,ms.Z) - ms.Z * ME_NuclPot(oi,oj,ms.zeta,ms.Z) + ME_overlap(oi,oj,ms.zeta,ms.Z) * pow(ms.c, 2) * (fac-1);
                } else {
                    H(i,j) = ME_Kinetic_nonrel(oi,oj,ms.zeta) - ms.Z * ME_NuclPot(oi,oj,ms.zeta,ms.Z);
                }
                S(i,j) = ME_overlap(oi,oj,ms.zeta,ms.Z);
            }
        }
        one = H;
    }

    void Operator::set_two_body_ham() {
        ModelSpace ms = modelspace;
        two_body_space = ms.two;
        Orbits orbs = ms.orbits;

        for (auto channels : two.Channels) {
            // cout << channels.first[0] << endl;
            chbra = two_body_space.get_channel(channels.first[0]);
            chket = two_body_space.get_channel(channels.first[1]);
            for (int idxbra=0; idxbra<chbra.get_number_states(); idxbra++) {
                for (int idxket=0; idxket<chket.get_number_states(); idxket++) {
                    if (idxbra > idxket) {continue;}
                    int a = chbra.get_indices(idxbra)[0];
                    int b = chbra.get_indices(idxbra)[1];
                    int c = chket.get_indices(idxket)[0];
                    int d = chket.get_indices(idxket)[1];
                    Orbit oa = orbs.get_orbit(a);
                    Orbit ob = orbs.get_orbit(b);
                    Orbit oc = orbs.get_orbit(c);
                    Orbit od = orbs.get_orbit(d);
                    double norm = 1.0;
                    if (a==b) {norm /= sqrt(2);}
                    if (c==d) {norm /= sqrt(2);}
                    double v;
                    // Commented out for debugging _coulomb
                    // v = _coulomb(oa, ob, oc, od, chket.J);
                    // v += _coulomb(oa, ob, od, oc, chket.J) * pow( (-1), ( floor( (oc.j+od.j)/2) - chket.J + 1) );
                    v *= norm;
                    // two.set_2bme(channels.first[0], channels.first[0], idxbra, idxket, v);
                    // two.set_2bme(channels.first[0], channels.first[0], idxbra, idxket, v);
                }
            }
        }
    }

    double g(double *r, size_t dim, void *params) {
        (void)(dim); /* avoid unused parameter warnings */
        int& Z = ((int *)params)[0];
        int& zeta = ((int *)params)[1];
        int& L = ((int *)params)[2];
        Orbit& oa = ((Orbit *)params)[3];
        Orbit& ob = ((Orbit *)params)[4];
        Orbit& oc = ((Orbit *)params)[5];
        Orbit& od = ((Orbit *)params)[6];

        return oa.eval_radial_function_rspace(r[0], zeta, Z, oa.e) * \
                        ob.eval_radial_function_rspace(r[1], zeta, Z, ob.e) * \
                        oc.eval_radial_function_rspace(r[0], zeta, Z, oc.e) * \
                        od.eval_radial_function_rspace(r[1], zeta, Z, od.e) * \
                        pow( min(r[0],r[1]), L) / pow( max(r[0],r[1]), (L+1) ) * oa.e * ob.e;
    }
    
    double Operator::_coulomb(Orbit oa, Orbit ob, Orbit oc, Orbit od, int J) {
        if (oa.e != oc.e) {return 0.0;}
        if (ob.e != od.e) {return 0.0;}
        int zeta = modelspace.zeta;
        int Z = modelspace.Z;
        int Lmin = floor( max( abs(oa.j-oc.j), abs(ob.j-od.j) )/2 );
        int Lmax = floor( min( (oa.j+oc.j), (ob.j+od.j) )/2 );
        double r = 0.0;
        for (auto L : irange(Lmin, Lmax+1)) {
            if ( (oa.l+oc.l+L)%2 == 1 ) {continue;}
            if ( (ob.l+od.l+L)%2 == 1 ) {continue;}
            int angular = wigner6j(oa.j*0.5, ob.j*0.5, J, od.j*0.5, oc.j*0.5, L) * \
                            wigner3j(oa.j*0.5, L, oc.j*0.5, -0.5, 0, 0.5) * \
                            wigner3j(ob.j*0.5, L, od.j*0.5, -0.5, 0, 0.5);

            // Lambda double integral
            // auto lambda = [=](double r1, double r2) {return eval_radial_function_rspace(r1, zeta, Z, oa.e) * \
            //                                                     eval_radial_function_rspace(r2, zeta, Z, ob.e) * \
            //                                                     eval_radial_function_rspace(r1, zeta, Z, oc.e) * \
            //                                                     eval_radial_function_rspace(r2, zeta, Z, od.e) * \
            //                                                     pow( min(r1,r2), L) / pow( max(r1,r2), (L+1) ) * oa.e * ob.e;};

            double integral_res, err;


            auto lower_limit = [=](double x) {return 0;};
            auto upper_limit = [=](double x) {return 8;};
            double xl[2] = { 0, 0 };
            double xu[2] = { 8, 8 };

            const gsl_rng_type *T;
            gsl_rng *re;

            gsl_monte_function G;
            G.f = &g;
            G.dim = 2;

            // Struct for Monte Carlo Integration
            struct my_f_params { int Z; int zeta; int L; Orbit oa; Orbit ob; Orbit oc; Orbit od; };
            
            struct my_f_params params = { Z, zeta, L, oa, ob, oc, od };
            G.params = &params;

            size_t calls = 5000000;
            gsl_rng_env_setup();

            T = gsl_rng_default;
            re = gsl_rng_alloc(T);

            gsl_monte_plain_state *s = gsl_monte_plain_alloc(2);
            gsl_monte_plain_integrate(&G, xl, xu, 3, calls, re, s, &integral_res, &err);
            gsl_monte_plain_free(s);

            cout << integral_res << endl;
            r += angular * integral_res;
        }
        r *= sqrt( (oa.j+1) * (ob.j+1) * (oc.j+1) * (od.j+1) * pow( (-1), floor((oa.j+oc.j)/2)+J ) );
        return r;
    }



    void Operator::orthogonalize(bool scalar) {
        ModelSpace ms = modelspace;
        two_body_space = ms.two;
        Orbits orbs = ms.orbits;
        Mat<double> L = chol(S);
        Mat<double> T = inv(L);
        one = T * one * T.t();
        if (scalar) {
            /*
            (ab:J|U1 x U2|cd:J) = [ (a|U|c) (b|U|d)
            - (phase) (a|U|d) (b|U|c)
            - (phase) (a|U|d) (b|U|c)
            + (phase) (a|U|c) (b|U|d) ] / 2 sqrt((1+del_ab) (1+del_cd))
            */
            for (auto channels : two.Channels) {
                TwoBodyChannel tbc = two_body_space.get_channel(channels.first[0]);
                Mat<double> U;
                U.zeros(tbc.get_number_states(), tbc.get_number_states());
                for (int idxbra=0; idxbra<tbc.get_number_states(); idxbra++) {
                    for (int idxket=0; idxket<tbc.get_number_states(); idxket++) {
                        int a = tbc.get_indices(idxbra)[0];
                        int b = tbc.get_indices(idxbra)[1];
                        int c = tbc.get_indices(idxket)[0];
                        int d = tbc.get_indices(idxket)[1];

                        Orbit oa = orbs.get_orbit(a);
                        Orbit ob = orbs.get_orbit(b);
                        Orbit oc = orbs.get_orbit(c);
                        Orbit od = orbs.get_orbit(d);

                        int phase_ab = pow( (-1), floor(oa.j+ob.j)/2 - tbc.J );
                        int phase_cd = pow( (-1), floor(oc.j+od.j)/2 - tbc.J );
                        U(idxbra, idxket) = ((1+phase_ab*phase_cd) * T(a,c) * T(b,d) - (phase_ab+phase_cd) * T(a,d) * T(b,c)) * 0.5;
                        if (a==b) {U(idxbra, idxket) /= sqrt(2);}
                        if (c==d) {U(idxbra, idxket) /= sqrt(2);}
                    }
                }
                two.Channels.at(channels.first).MEs = U * two.Channels.at(channels.first).MEs * U.t();
           }

        } else {
            cout << "Not impletmented yet" << endl;
            S = eye(orbs.get_num_orbits(), orbs.get_num_orbits());
        }
    }



    void Operator::print_operator(){
        Orbits orbs = modelspace.orbits;
        OneBodySpace one_body_space = modelspace.one;
        TwoBodySpace two_body_space = modelspace.two;
        for (int ichbra=0; ichbra<one_body_space.number_channels; ichbra++) {
            for (int ichket=0; ichket<one_body_space.number_channels; ichket++) {
                vector <int> oneVec = one_body_space.channels[ichbra];
                for (auto i : one_body_space.channels[ichbra]) {
                    for (auto j : one_body_space.channels[ichket]) {  
                        Orbit oi = orbs.get_orbit(i);
                        Orbit oj = orbs.get_orbit(j);
                        if (! ( abs(oi.j-oj.j) <= 2*rankJ && 2*rankJ <= oi.j+oj.j ) ) {continue;}
                        if ( pow( (-1), (oi.j+oj.l) ) * rankP == -1 ) {continue;}
                        cout << setw(2) << i << fixed  << setw(3) << j << one(i,j) << S(i,j) << endl;
                    }
                }
                two.print_two_body_operator();
            }
        }

    }




