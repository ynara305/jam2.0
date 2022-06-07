#ifndef jam2_meanfield_ScalarPotential_h
#define jam2_meanfield_ScalarPotential_h

#include <Pythia8/Settings.h>
#include <jam2/meanfield/PotentialType.h>

namespace jam2 {

class ScalarPotential : public PotentialType
{
public:
  ScalarPotential(Pythia8::Settings* s);
  ~ScalarPotential() {};
  void Vmd(double p1, double p2,double vfac1,double vfac2,double fi1,double fi2,double fj1,double fj2);
  void dVdns(double rij, double rji);
  void dVdmd(double p1, double p2,double fs,double fv, double pmi1,double pmj1,double pmi2,double pmj2);
  void dVdp();
  void dVdpm();
};

} // end namespace jam2
#endif
