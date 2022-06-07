#ifndef jam2_fluid_Fluid2Particle_H
#define jam2_fluid_Fluid2Particle_H

#include <random>
#include <vector>
#include <Pythia8/Settings.h>
#include <Pythia8/Basics.h>
#include <Pythia8/ParticleData.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/collision/EventParticle.h>
#include <jam2/fluid/Fluid.h>
#include <jam2/fluid/FreezeOut.h>

namespace jam2 {

class SMHadron
{
private:
  int id0, baryon, charge, strange, spin, istat;
  double mass, width, mult;
  Pythia8::ParticleDataEntry* pdePtr;
public:
  SMHadron(int id1, int ba, int ch, int str, int sp,int is,
	  double m, double w, Pythia8::ParticleDataEntry* p):
      id0(id1), baryon(ba), charge(ch), strange(str),
      spin(sp),istat(is), mass(m), width(w), pdePtr(p) { }
  int id() const {return id0;}
  int baryonType() const {return baryon;}  // baryon number
  int chargeType() const {return charge;}  // charge
  int strangeType() const {return strange;}
  int spinType() const {return spin;}
  int stat() const {return istat;}
  Pythia8::ParticleDataEntry* particle() {return pdePtr;}
  double m0() const {return mass;}
  double mWidth() const {return width;}
  double multiplicity() const {return mult;}
  void setMult(double p) {mult=p;}
};

using Pythia8::Vec4;

class Fluid2Particle
{
private:
  Pythia8::Settings *settings;
  JamParticleData* particleTable;
  Pythia8::Rndm* rndm;
  FreezeOut *freezeout;
  Fluid *fluid;
  int optEcon, mstc144,useMTable,isDebug;
  int nTest,nFreezeout;
  double facStat,volume, volF, tFreezeCut;
  std::vector<SMHadron> smHadrons;
  Pythia8::Vec4 fConv;
  double fConvB;
  std::default_random_engine generator;
  static int nT, nBmu, nSmu;
  static double tMin,bMuMin,sMuMin;
  static double dTmp,dbMu,dsMu;
  std::vector<double> pFreezeoutTable;
public:
  std::vector<EventParticle*> outgoing;
  Fluid2Particle(Pythia8::Settings* s,JamParticleData* pd, 
	  FreezeOut* fr, Fluid* fl,Pythia8::Rndm* r);
  ~Fluid2Particle();
  void makeHadronTable();
  int generate_poisson(double a);
  void convert(double gtime,int icheck);
  int sampling4c(int check,int& bar,int& ch,int& str,Vec4& ptot,double& nftot);
  int sampling1p(int check,int& bar,int& ch,int& str,Vec4& ptot,double& nftot);
  void p2jam(SMHadron* hc,Vec4& rc,Vec4& pc,double mh, 
	     double v1,double v2, double v3);
  SMHadron* selectParticle(double& tdns);
  SMHadron* selectParticle(double& tdns,double tch,double bmu,double smu);
  double totalMultiplicity(double tch,double bmu,double smu);
  double pdensity(double pm,double tf,double mu,int istat);
  double thermaldist3(double pm,double tf,double mu,int a);
  void recover_mom(Vec4& ptot,Vec4& pconv);
  double BWmass(double m0, double wid, double mmin, double mmax);
  void makeMultiplicityTable();
};

}
#endif
