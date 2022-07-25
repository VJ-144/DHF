#ifndef PhysicalConstants_h
#define PhysicalConstants_h
#include <math.h>
namespace PhysConst {
const double HBARC              = 197.3269718  ;                  // reduced Planck constant * speed of light in MeV * fm
const double HARTREE            = 27.21138602;                    // 1 Hartree in eV
const double c = 137.035999084;
// const double c = 137.0373
const double Alpha = 1.0 / c;
// Math constants
const double SQRT2    = sqrt(2.0);
const double INVSQRT2 = 1.0 / SQRT2;
const double PI       = 4.0*atan(1.0) ;
} // namespace PhysConst
#endif
