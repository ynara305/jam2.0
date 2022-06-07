#ifndef jam2_meanfield_ScalarVectorPotential_h
#define jam2_meanfield_ScalarVectorPotential_h

#include <Pythia8/Settings.h>
#include <jam2/meanfield/PotentialType.h>

namespace jam2 {

class ScalarVectorPotential : public PotentialType
{
private:
  double vf1,vf2;
  // for non-linear sigma
  double mSigma, mOmega, mSigmaFM2,mSigma2;
  double gs,gv,g2,g3,g4,gsp,gvp,GS,G2,G3;
  double G4,G24,G34,GS4;

public:
  ScalarVectorPotential(Pythia8::Settings* s);
  ~ScalarVectorPotential() { };
  void setPotentialParam(EventParticle* p);
  void Vmd(double p1, double p2,double vfac1,double vfac2,double fi1,double fi2,double fj1,double fj2);
  void dVdns(double rij, double rji);
  void dVdmd(double p1, double p2,double fs,double fv, double pmi1,double pmi2,double pmj1,double pmj2);
  void dVdp();
  void dVdpm();
  double dSigma(double spotv, double rhos);
 // function used for bisec method.
  double funcSigma(double x,double scalarDens) {
      return -scalarDens*pow2(1.0+0.5*g4*x*x)
          + mSigmaFM2*x + g2*x*x + g3*x*x*x
          + 1.0/6.0*g2*g4*pow4(x)+0.25*g3*g4*pow5(x);
  }

};

} // end namespace jam2
#endif
