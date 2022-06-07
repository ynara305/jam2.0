#ifndef jam2_meanfield_RQMDv_h
#define jam2_meanfield_RQMDv_h

#include <jam2/meanfield/MeanField.h>
#include <jam2/meanfield/PotentialType.h>

namespace jam2 {

class RQMDv : public MeanField
{
protected:
  std::vector<EventParticle*> part; // particles for potentrial interaction
  std::vector<std::vector<double> > rhom; // rho_{ij}
  std::vector<double> rho;   // invariant baryon density
  std::vector<double> rhog,rhog2;  // rho^{gam-1}
  std::vector<Vec4> vmom4; // momentum-dependent potential
  std::vector<Vec4> JB;       // baryon current
  std::vector<Vec4> Vdot;    // time derivative of vector potential
  std::vector<Vec4> force;  // force for p
  std::vector<Vec4> forcer; // force for r
  static bool firstCall;

  double t1,t3,t3f,gam,vex1,vex2,pmu1,pmu2;
  double wG,facG;
  double psq1, psq2;
  Pythia8::Vec4 dp2ijpi, dp2jipi, dp2ijpj, dp2jipj;

  PotentialType *potential;

public:
  RQMDv(Pythia8::Settings* s);
  ~RQMDv();
  void evolution(std::list<EventParticle*>& plist,double t, double dt,int step);
  void init(std::list<EventParticle*>& plist) { };
  Pythia8::Vec4 computeEnergy(std::list<EventParticle*>& plist,int step);
  void singleParticlePotential();

  std::vector<EventParticle*>& getParticle() {return part;}
  int particleSize()                         {return part.size();}
  void add(EventParticle* p) { part.push_back(p); }

  void addRho(int i, double a)                 {rho[i] +=a;}
  void addVmom4(int i, const Pythia8::Vec4& a) {vmom4[i] +=a;}
  void addJB(int i, const Pythia8::Vec4& a)    {JB[i] +=a;}
  void addForce(int i, const Pythia8::Vec4& f) {force[i] +=f;}
  void addForceR(int i, const Pythia8::Vec4& f){forcer[i] +=f;}
  void addVdot(int i, const Pythia8::Vec4& v)  {Vdot[i] +=v;}
  double getRhog(int i)                        {return rhog[i];}
  double getRho(int i)                         {return rho[i];}
  Pythia8::Vec4   getJB(int i)                 {return JB[i];}
  Pythia8::Vec4   getForce(int i)             {return force[i];}
  Pythia8::Vec4   getForceR(int i)            {return forcer[i];}

  void qmdMatrix();
  void computeForce();

  void qmdMatrix1();
  void qmdMatrix2();
  void qmdMatrix3();
  void computeForce1();
  void computeForce2();
  void computeForce3();
  void computeForce3org();
  void computeForceP();

  Pythia8::Vec4 facV(int i, Vec4& p);
  void computeVdot();
  Vec4 getPcan(Vec4 p, Vec4 v) {
    double msq=p.m2Calc();
    p += v; p[0] = sqrt(msq + p.pAbs2());
    return p;
  }
  Vec4 getPkin(Vec4 p, Vec4 v, double m) {
    p -= v; p[0] = sqrt(m*m + p.pAbs2());
    return p;
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
}
#endif
