#ifndef jam2_meanfield_PotentialType_h
#define jam2_meanfield_PotentialType_h

#include <Pythia8/Basics.h>
#include <Pythia8/Settings.h>
#include <jam2/collision/EventParticle.h>

namespace jam2 {

class PotentialType
{
protected:
  Pythia8::Settings* settings;
  bool withMomDep;
  int optTwoBodyDistance,transportModel;
  int optP0dev,optPV,optScalarDensity,optBaryonCurrent,optVectorPotential;
  int eosType,optPotentialType;
  double t0,t1,t2,t2f,t3,t3f,t5,gam,rho0,alpha,beta;
  double vex1,vex2,pmu1,pmu2;
  int optLambdaPot,optPotentialDensity,optStrangeBaryonPot;
  double facPotL[10],facPotS[10];

  double pfac1,pfac2,pfac3,pfac4,pfac5,pfac6,gama1,gamb1,gama2,gamb2;
  double pFac1,pFac2,pFac3,gam1,gam2;

  double wG;
  double rhoi,rhoj;
  double rhomij, rhomji,rhogi,rhogj,rhog2i,rhog2j;
  double rhosi, rhosj, rhosgi, rhosgj;
  double spotvi,spotvj;
  Pythia8::Vec4 JBi,JBj;
  double fi, fj, fengi,fengj;
  Pythia8::Vec4 vk1, vk2, v1,v2;
  double p01, p02;
  int bar1,bar2;
  Pythia8::Vec4 devV1,devV2;
  double s1,s2;
  double fskys1,fskys2,fskyv1,fskyv2;
  static bool firstCall;
  
public:
  PotentialType(Pythia8::Settings* s);
  virtual ~PotentialType() { };
  virtual void Vmd(double p1, double p2,double vfac1,double vfac2,double fi1,double fi2,double fj1,double fj2)=0;
  virtual void dVdns(double rij, double rji)=0;
  virtual void dVdmd(double p1, double p2,double fs,double fv,
    double pmi1,double pmi2,double pmj1,double pmj2)=0;
  virtual void dVdp()=0;
  virtual void dVdpm()=0;
  bool isMomDep() const {return withMomDep;}
  double getAlpha() const {return alpha;}
  double getBeta() const {return beta;}
  double getGam() const {return gam;}
  double getT0() const {return t0;}
  double getT1() const {return t1;}
  double getT2() const {return t2;}
  double getT3() const {return t3;}
  double getT5() const {return t5;}
  double getVex1() const {return vex1;}
  double getVex2() const {return vex2;}
  double getPmu1() const {return pmu1;}
  double getPmu2() const {return pmu2;}

  virtual void setPotentialParam(EventParticle* p);
  int isHyperon(EventParticle* p);
  //void setPotentialParam(EventParticle* p,double t0,double t3,double t5,
  //  double gam, double vex1,double vex2,double pmu1,double pmu2);

  int potentialType() const {return optPotentialType;}
  double matrix(int i1, int i2) {
    // poentials depend on baryon density
    if(optPotentialDensity>0) return 1.0;

    //Lambda potential is the same form as nucleons,
    // but strength is controlled by parameters
    // for the multiplication factors for Hyperon potential:
    //if(optLambdaPot==0) return 1.0;

    if(i1==1 && i2 ==1) return 1.0; // NN
    else if(i1==1 && i2 ==3) return 0.0; // NL
    else if(i1==3 && i2==3) return 0.0; // LL
    else return 1.0;
  }

  void set1(Pythia8::Vec4& p,Pythia8::Vec4& pkin,const Vec4& JB,
      double f,double feng,double rho,double rhog,double rhog2, double spotv,
      double rhos, double rhosg,int bar,double pf1,double pf2,double g1,double pf3,double g2) {
    vk1 = pkin/pkin[0]; v1=p/p[0]; p01=p[0];
    fi=f; fengi=feng; 
    rhosi=rhos; rhosgi=rhosg;
    rhoi= rho; rhogi=rhog; rhog2i=rhog2;
    spotvi=spotv;
    JBi=JB;bar1=bar;
    pfac1=pf1;
    pfac3=pf2;
    pfac5=pf3;
    gama1=g1;
    gamb1=g2;
  }

  void set2(Pythia8::Vec4& p,Pythia8::Vec4& pkin,const Vec4& JB,
      double f,double feng,double rho,double rhog,double rhog2, double spotv,
      double rhos, double rhosg,int bar,double pf1,double pf2,double g1,double pf3,double g2) {
    vk2 = pkin/pkin[0]; v2=p/p[0]; p02=p[0];
    fj=f; fengj=feng; 
    rhosj=rhos;
    rhosgj=rhosg;
    rhoj=rho; rhogj=rhog; rhog2j=rhog2;
    spotvj=spotv;
    JBj=JB;bar2=bar;
    pfac2=pf1;
    pfac4=pf2;
    pfac6=pf3;
    gama2=g1;
    gamb2=g2;
  }

  double devVme(double psq,double pmu,double vex) {
    double fac1 = 1.0 - psq/pmu;
    return -vex/(pmu*fac1*fac1);
  }
  double devVmd(double psq,double pmu,double vex) {
    return vex/(1.0 - psq/pmu);
  }
  double devVmd2(double psq,double vf1,double vf2, double pf1,double pf2) {
    return vf1*vex1/(1.0 - psq/(pf1*pmu1)) + vf2*vex2/(1.0 - psq/(pf2*pmu2));
  }
  double devVme2(double psq,double vf1,double vf2,double pf1,double pf2) {
    double fac1 = 1.0 - psq/(pf1*pmu1);
    double fac2 = 1.0 - psq/(pf2*pmu2);
    return -vf1*vex1/(pmu1*fac1*fac1) - vf2*vex2/(pmu2*fac2*fac2);
  }

  Vec4 facV(Vec4& JB, double rho, double rhog, double rho2,Vec4& vkin);
  void devV(const Vec4& Ai, const Vec4& Aj);

  double vmomsi, vmomsj;
  double vmom4i, vmom4j;
  double fskyi, fskyj;
  double fmomdi, fmomei, fmomdj, fmomej;
  Pythia8::Vec4 forceri, forcerj;
  double fmomds1, fmomds2, fmomdv1,fmomdv2;
  double facsk,facmom;


};



}
#endif
