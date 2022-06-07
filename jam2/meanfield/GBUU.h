#ifndef jam2_meanfield_GBUU_h
#define jam2_meanfield_GBUU_h

#include <jam2/meanfield/MeanField.h>
#include <jam2/meanfield/PotentialType.h>

namespace jam2 {

class GBUU : public MeanField
{
private:
  std::vector<EventParticle*> part;
  double **rhom;
  double *rho, *rhog, *vmoms;
  Pythia8::Vec4 *force, *forcer;
  double t1,t3,t3f,gam,vex1,vex2,pmu1,pmu2;
  double wG,facG;
  int overSample;
  PotentialType *potential;

public:
  GBUU(Pythia8::Settings* s);
  ~GBUU() {delete potential; }
  std::vector<EventParticle*>& getParticle() {return part;}
  void add(EventParticle* p) { part.push_back(p); }
  void evolution(std::list<EventParticle*>& plist,double t,double dt,int step);
  void singleParticlePotential();
  void init(std::list<EventParticle*>& plist) {};
  void qmdMatrix();
  void computeForce();
  Pythia8::Vec4 computeEnergy(std::list<EventParticle*>& plist,int step);
};
}
#endif
