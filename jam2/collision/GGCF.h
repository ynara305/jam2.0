#ifndef jam2_collision_GGCF_h
#define jam2_collision_GGCF_h
#include <Pythia8/Basics.h>
#include <jam2/hadrons/JamStdlib.h>

namespace jam2 {

// Gribov-Glauber color fluctuation (GGCF) model
// M.Alvioli, M. Strikman, Phys. Lett.B722, 347-354(2013). arXiv:1301.0728
// ATLAS Eur. Phys. J. C (2016) 76 :199
class GGCF
{
private:
  Pythia8::Rndm* rndm;
  double sigmaTot;    // total cross section (input)
  double omegaSigma;  // omega_sigma (input)
  double Norm;
  double sigma0;
  double Omega;
  double sigMax;
public:
  GGCF(Pythia8::Rndm* r,double sig=94.8, double ome=0.2) : rndm(r),sigmaTot(sig),omegaSigma(ome) {
    Norm=1.0;sigma0=72.6;Omega=1.01;
    sigMax=500.0;
  }
  double Prob(double sig, double sig0, double omega) {
    return sig/(sig+sig0)*exp(-pow2((sig/sig0-1.0)/omega));
  }
  double sample() {
    double xm = peak();
    double fmax = Prob(xm,sigma0,Omega);
    double x;
    do {
      x = rndm->flat()*5*sigmaTot;
    }while(rndm->flat() > Prob(x,sigma0,Omega)/fmax);
    return x;
  }
  bool computeParameter(double sig);
  bool computeParameter2(double sig, double omeg);
  void setParam(double sig, double omeg) {sigmaTot=sig; omegaSigma=omeg;}
  void setOmega(double ome) {Omega=ome;}
  double sum();
  double sum2();
  double sigma();
  double omega();
  double norm() {return Norm;}
  double sig0() {return sigma0;}
  double omeg() {return Omega;}
  double peak() {
    double s3=Pythia8::pow3(sigma0);
    double o2 = pow2(Omega);
    double o4= o2*o2;
    double x1 = pow(sqrt(27*o4-16.0)*s3/(4*pow(3,1.5))+s3*o2/4.0, 1.0/3.0);
    return x1 + sigma0*sigma0/(3*x1);
  }

//               4        3    3  2                        2
//     sqrt(27*o  - 16) s    s  o  1/3                   s
//x = (------------------- + -----)    + ----------------------------------]
//              3/2            4                     4        3    3  2
//           4 3                            sqrt(27 o  - 16) s    s  o  1/3
//                                       3 (------------------- + -----)
//                                                   3/2            4
//                                                4 3


};
} // end namespace jam2
#endif
