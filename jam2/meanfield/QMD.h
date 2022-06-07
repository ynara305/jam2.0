#ifndef jam2_meanfield_QMD_h
#define jam2_meanfield_QMD_h

#include <jam2/meanfield/MeanField.h>
#include <jam2/meanfield/PotentialType.h>

namespace jam2 {

class QMD : public MeanField
{
protected:
  std::vector<EventParticle*> part;
  std::vector<Vec4> force;
  std::vector<Vec4> forcer;
  std::vector<double> rho1,rho2,rho3;
  std::vector<double> rhog2,rhog3;
  std::vector<double> vmoms;
  std::vector<std::vector<double> > rhom;

  double t1,t3,t5,alpha,beta,gam,vex1,vex2,pmu1,pmu2;
  double rho0,widG,wG,facG;
  int overSample;
  //double aLambda,bLambda,cLambda,gamL1,gamL2;

  //PotentialParam* potParam;
  PotentialType *potential;

public:
  QMD(Pythia8::Settings* s);
  ~QMD();
  void evolution(std::list<EventParticle*>& plist,double t, double dt,int step);
  void initMatrix(std::list<EventParticle*>& plist,double t);
  Pythia8::Vec4 computeEnergy(list<EventParticle*>& plist, int step);
  void singleParticlePotential();
  void init(std::list<EventParticle*>& plist) {};
  void qmdMatrix();
  void computeForce();
  //void setPotentialParam(EventParticle* p);

  std::vector<EventParticle*>& getParticle() {return part;}
  int particleSize() {return part.size();}
  void add(EventParticle* p) {
      potential->setPotentialParam(p);
      part.push_back(p);}

  void addRho1(int i, double a) {rho1[i] +=a;}
  void addRho2(int i, double a) {rho2[i] +=a;}
  void addRho3(int i, double a) {rho3[i] +=a;}
  void addVmoms(int i, double a) {vmoms[i] +=a;}
  void addForce(int i, const Pythia8::Vec4& f) {force[i] +=f;}
  void addForceR(int i, const Pythia8::Vec4& f) {forcer[i] +=f;}
  double getRhog2(int i) {return rhog2[i];}
  double getRhog3(int i) {return rhog3[i];}
  double getRho(int i) {return rho1[i];}
  double getRho2(int i) {return rho2[i];}
  double getRho3(int i) {return rho3[i];}
};
}
#endif
