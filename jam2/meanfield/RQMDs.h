#ifndef jam2_meanfield_RQMDs_h
#define jam2_meanfield_RQMDs_h

#include <jam2/meanfield/MeanField.h>

namespace jam2 {

class RQMDs : public MeanField
{
protected:
  std::vector<EventParticle*> part;
  std::vector<Vec4> force;
  std::vector<Vec4> forcer;
  std::vector<double> rho;
  std::vector<double> rhog;
  std::vector<double> vmoms;
  std::vector<std::vector<double> > rhom;
  Vec4 eFree,eFree0;
  static bool firstCall;

  double t1,t3,t3f,gam,vex1,vex2,pmu1,pmu2;
  double wG,facG;
  double psq1, psq2;
  Pythia8::Vec4 dp2ijpi, dp2jipi, dp2ijpj, dp2jipj;

public:
  RQMDs(Pythia8::Settings* s);
  ~RQMDs() {clearMatrix(); }
  void evolution(std::list<EventParticle*>& plist,double t,double dt,int step);
  void initMatrix(std::list<EventParticle*>& plist,double t,int step);
  Pythia8::Vec4 computeEnergy(std::list<EventParticle*>& plist, int step);
  void singleParticlePotential(bool opt);
  void init(std::list<EventParticle*>& plist);
  void clearMatrix();

  std::vector<EventParticle*>& getParticle() {return part;}
  int particleSize()                         {return part.size();}
  void add(EventParticle* p)                 {part.push_back(p);}

  void addRho(int i, double a)                 {rho[i] +=a;}
  void addVmoms(int i, double a)               {vmoms[i] +=a;}
  void addForce(int i, const Pythia8::Vec4& f) {force[i] +=f;}
  void addForceR(int i, const Pythia8::Vec4& f){forcer[i] +=f;}
  double getRhog(int i)                        {return rhog[i];}

  //double facSDensity(Vec4& p, double m);
  void qmdMatrix(double a=1.0);
  void computeForce();

  void qmdMatrix1(double a=1.0);
  void qmdMatrix2(double a=1.0);
  void qmdMatrix3(double a=1.0);
  void computeForce1();
  void computeForce2();
  void computeForce3();

  void RecoverEnergy(std::list<EventParticle*>& plist);
  void RecoverEnergy2(std::list<EventParticle*>& plist);
  double funcEnergy(std::list<EventParticle*>& plist,double a);

  double facScalar(Vec4& p, double m) {
    if(optScalarDensity==0) return p.mCalc()/p[0];
    else return m/sqrt(m*m + p.pAbs2());
  }

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
}  // namespace jam2
#endif
