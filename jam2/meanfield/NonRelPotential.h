#ifndef jam2_meanfield_NonRelPotential_h
#define jam2_meanfield_NonRelPotential_h

#include <Pythia8/Settings.h>
#include <jam2/meanfield/PotentialType.h>

namespace jam2 {

class NonRelPotential : public PotentialType
{
public:
  NonRelPotential(Pythia8::Settings* s);
  ~NonRelPotential() { };
  void setPotentialParam(EventParticle* p);
  void setLambdaPotential(double* fp);
  void setSigmaPotential(double* fp);
  void Vmd(double p1, double p2,double vfac1,double vfac2,double fi1,double fi2,double fj1,double fj2) {};
  void dVdns(double rij, double rji) {};
  void dVdmd(double p1, double p2,double fs,double fv, double pmi1,double pmj1,double pmi2,double pmj2) {};
  void dVdp() {};
  void dVdpm() {};
};

} // end namespace jam2
#endif
