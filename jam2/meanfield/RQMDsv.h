#ifndef jam2_meanfield_RQMDsv_h
#define jam2_meanfield_RQMDsv_h

#include <jam2/meanfield/MeanField.h>
#include <jam2/meanfield/PotentialType.h>

namespace jam2 {

class RQMDsv : public MeanField
{
protected:
  std::vector<EventParticle*> part; // particles for potentrial interaction
  std::vector<std::vector<double> > rhom; // rho_{ij}
  std::vector<double> rhos;   // scalar density
  std::vector<double> rho;   // invariant baryon density
  std::vector<double> rhog, rhog2;  // rho^{gam-1}
  std::vector<double> rhosg;  // rhos^{gam-1}
  std::vector<double> vmoms; // momentum-dependent scalar potential
  std::vector<Vec4> vmom4; // momentum-dependent vector potential
  std::vector<Vec4> JB;       // baryon current
  std::vector<Vec4> Vdot;    // time derivative of vector potential
  std::vector<Vec4> force;  // force for p
  std::vector<Vec4> forcer; // force for r
  static bool firstCall;
  int optPotentialType,transportModel;
  int overSample;
  PotentialType *potential;

  double t0,t2,t1,t3,gam;
  double wG;

  // for non-linear sigma
  //double mSigma, mOmega, mSigmaFM2,mSigma2;
  //double gs,gv,g2,g3,g4,gsp,gvp,GS,G2,G3,t1f;
  //double G4,G24,G34,GS4;
  double scalarDens;

  double vex1,vex2,pmu1,pmu2;

public:
  RQMDsv(Pythia8::Settings* s);
  ~RQMDsv();
  void evolution(std::list<EventParticle*>& plist,double t, double dt,int step);
  void init(std::list<EventParticle*>& plist) { };
  Pythia8::Vec4 computeEnergy(std::list<EventParticle*>& plist,int step);
  void singleParticlePotential(bool opt);
  void setScalarPotential();

  std::vector<EventParticle*>& getParticle() {return part;}
  int particleSize()                         {return part.size();}
  void add(EventParticle* p) { 
      potential->setPotentialParam(p);
      part.push_back(p); }

  void addRhos(int i, double a)                {rhos[i] +=a;}
  void addRho(int i, double a)                 {rho[i] +=a;}
  void addVmoms(int i, double a)               {vmoms[i] +=a;}
  void addVmom4(int i, const Pythia8::Vec4& a) {vmom4[i] +=a;}
  void addJB(int i, const Pythia8::Vec4& a)    {JB[i] +=a;}
  void addForce(int i, const Pythia8::Vec4& f) {force[i] +=f;}
  void addForceR(int i, const Pythia8::Vec4& f){forcer[i] +=f;}
  void addVdot(int i, const Pythia8::Vec4& v)  {Vdot[i] +=v;}
  double getRhog(int i)                        {return rhog[i];}
  double getRhog2(int i)                       {return rhog2[i];}
  double getRho(int i)                         {return rho[i];}
  double getRhos(int i)                        {return rhos[i];}
  double getRhosg(int i)                       {return rhosg[i];}
  double getVmoms(int i)                       {return vmoms[i];}
  Pythia8::Vec4   getJB(int i)                 {return JB[i];}
  Pythia8::Vec4   getForce(int i)              {return force[i];}
  Pythia8::Vec4   getForceR(int i)             {return forcer[i];}

  void qmdMatrix();
  void computeForce();
  void computeBUUForce();

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

  // for non-linear sigma model
  void sigmaField(double a=1.0);
  double sigma_bisec(double sigma1, int i=0);
  Pythia8::Vec4 facV(int i, Vec4& p);
  double devVme(double psq,double pmu,double vex) {
    double fac1 = 1.0 - psq/pmu;
    return -vex/(pmu*fac1*fac1);
  }
  double devVmd(double psq,double pmu,double vex) {
    return vex/(1.0 - psq/pmu);
  }

};
}
#endif
