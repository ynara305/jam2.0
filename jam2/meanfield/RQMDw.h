#ifndef jam2_meanfield_RQMDw_h
#define jam2_meanfield_RQMDw_h

#include <jam2/meanfield/MeanField.h>
#include <jam2/meanfield/PotentialType.h>

namespace jam2 {

class RQMDw : public MeanField
{
private:
  std::vector<EventParticle*> part;
  double **rhom;
  Pythia8::Vec4 *vmom4, *JB, *Vdot;
  double *rho, *rhos,*vmoms;
  Pythia8::Vec4 *force, *forcer;

  int optVdot,optBaryonCurrent,nonlinearSigma;
  double mSigma, mOmega, mSigmaFM2,mSigma2;
  double gs,gv,g2,g3,g4,gsp,gvp,GS,G2,G3;
  double G4,G24,G34,GS4;
  double wG,facG;
  double t1,t1f,t2;
  double vex1,vex2,pmu1,pmu2;
  double scalarDens, epsEnergy;
  Vec4 *vpot;
  Vec4 pFree;
  double psq1, psq2;
  Pythia8::Vec4 dp2ijpi, dp2jipi, dp2ijpj, dp2jipj;

  PotentialType *potential;
  int optMethod;

public:
  RQMDw(Pythia8::Settings* s);
  ~RQMDw() { }
  void evolution(std::list<EventParticle*>& plist,double t,double dt,int step);
  Pythia8::Vec4 computeEnergy(std::list<EventParticle*>& plist,int step);
  void singleParticlePotential();
  void init(std::list<EventParticle*>& plist) { };

  void qmdMatrix0(double a=1.0);
  void computeForce();
  void computeForceP();

  void qmdMatrix(double a=1.0);
  void qmdMatrix1(double a=1.0);
  void qmdMatrix2(double a=1.0);
  void qmdMatrix3(double a=1.0);
  void qmdMatrix4(double a=1.0);
  void computeForce1();
  void computeForce2();
  void computeForce3();
  void computeForce4();
  void scalarDensity1(double a=1.0);
  void scalarDensity2(double a=1.0);
  void scalarDensity3(double a=1.0);
  void sigmaField(double a=1.0);
  Pythia8::Vec4 facV(int i, Vec4& p);
  void computeVdot();

  void RecoverEnergy(std::list<EventParticle*>& plist);
  double funcEnergy(std::list<EventParticle*>& plist,double a);
  void RecoverEnergy2();
  double funcEnergy2(double a);
  void RecoverEnergy3();
  double funcEnergy3(double a);
  void RecoverEnergy4();
  double funcEnergy4(double a);

  std::vector<EventParticle*>& getParticle() {return part;}
  void add(EventParticle* p) { part.push_back(p); }

  double dSigma1(int i) {
    if(!nonlinearSigma) return t1f;
    //double sigma=(part[i]->pots()-vmoms[i])/(t1*part[i]->sfactor());
    double sigma=(part[i]->pots()-vmoms[i])/t1;
    return t1*GS/(mSigma2+G2*sigma+G3*sigma*sigma);
    //return t1f*mSigma2/(mSigma2+G2*sigma+G3*sigma*sigma);
  }
  double dSigma(int i) {
    if(!nonlinearSigma) return t1f;
    //double sigma=(part[i]->pots()-vmoms[i])/(t1*part[i]->sfactor());
    double sigma=(part[i]->pots()-vmoms[i])/t1;
    return t1*GS*pow2(1.0+G4*sigma*sigma)/
	( mSigma2+G2*sigma+G3*sigma*sigma
	  +G24*pow3(sigma)+G34*pow4(sigma)
	 -2*(1.0+G4*pow2(sigma))*GS4*sigma*rhos[i] ) ;
  }


  Vec4 getPcan(Vec4 p, Vec4 v) {
    double msq=p.m2Calc();
    p += v; p[0] = sqrt(msq + p.pAbs2());
    return p;
  }
  /*
  Vec4 getPkin(Vec4 p, Vec4 v,double m) {
    //double msq=p.m2Calc();
    p -= v; p[0] = sqrt(m*m + p.pAbs2());
    return p;
  }
  */

  double sigma_bisec(double sigma1, int i=0);
  double sigma_iterate();
  double sigma_newton(double sigma1=0.0);

  // function used for iteration.
  double fSigma0(double x) {
      return (scalarDens - g2*x*x - g3*x*x*x)/mSigmaFM2;
  }
  double fSigma(double x) {
      return (scalarDens*pow2(1.0+0.5*g4*x*x)
	   - g2*x*x - g3*x*x*x
	  - 1.0/6.0*g2*g4*pow4(x)-0.25*g3*g4*pow5(x) )/mSigmaFM2;
  }

  // function used for bisec method.
  double funcSigma0(double x) {
      return -scalarDens + mSigmaFM2*x + g2*x*x + g3*x*x*x;
  }
  double funcSigma(double x) {
      return -scalarDens*pow2(1.0+0.5*g4*x*x)
	  + mSigmaFM2*x + g2*x*x + g3*x*x*x
	  + 1.0/6.0*g2*g4*pow4(x)+0.25*g3*g4*pow5(x);
  }

  // for Newton method
  double funcdSigma(double x) {
    return //-scalarDens*pow2(1.0+g4*x*x) +
	 -2*(1.0+g4*pow2(x))*g4*x*scalarDens
	+ mSigmaFM2 + 2*g2*x + 3*g3*x*x
	  + 2./3.*g2*g4*pow3(x) + 5./4.*g3*g4*pow4(x);
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
