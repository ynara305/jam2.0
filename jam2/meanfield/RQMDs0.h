#ifndef jam2_meanfield_RQMDs0_h
#define jam2_meanfield_RQMDs0_h

#include <jam2/meanfield/MeanField.h>
#include <jam2/meanfield/PotentialType.h>

namespace jam2 {

class RQMDs0 : public MeanField
{
private:
  std::vector<EventParticle*> part;
  double **rhom;
  double *rho, *rhog, *vmoms;
  Pythia8::Vec4 *force, *forcer;

  double t1,t3,t3f,gam,vex1,vex2,pmu1,pmu2;
  double wG,facG;
  double psq1, psq2;
  Pythia8::Vec4 dp2ijpi, dp2jipi, dp2ijpj, dp2jipj;
  PotentialType *potential;

public:
  RQMDs0(Pythia8::Settings* s);
  ~RQMDs0() {delete potential; }
  std::vector<EventParticle*>& getParticle() {return part;}
  void add(EventParticle* p) { part.push_back(p); }
  void evolution(std::list<EventParticle*>& plist,double t, double dt,int step);
  void singleParticlePotential();
  void init(std::list<EventParticle*>& plist) {};
  void qmdMatrix1();
  void qmdMatrix2();
  void qmdMatrix3();
  void computeForce1();
  void computeForce2();
  void computeForce3();
  Pythia8::Vec4 computeEnergy(std::list<EventParticle*>& plist,int step);
  double facSDensity(Pythia8::Vec4& r, double m);

  double devVme(double psq,double pmu,double vex) {
    double fac1 = 1.0 - psq/pmu;
    return -vex/(pmu*fac1*fac1);
  }
  double devVmd(double psq,double pmu,double vex) {
    return vex/(1.0 - psq/pmu);
  }
  double devVme(double psq) {
    double fac1 = 1.0 - psq/pmu1;
    double fac2 = 1.0 - psq/pmu2;
    return -vex1/(pmu1*fac1*fac1) - vex2/(pmu2*fac2*fac2);
  }
  double devVmd(double psq) {
    return vex1/(1.0 - psq/pmu1) + vex2/(1.0 - psq/pmu2);
  }
  void distanceP1(const Pythia8::Vec4& p1,const Pythia8::Vec4& p2);
  void distanceP2(const Pythia8::Vec4& p1,const Pythia8::Vec4& p2);
  void distanceP3(const Pythia8::Vec4& p1,const Pythia8::Vec4& p2);
};
}
#endif
