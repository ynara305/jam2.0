#ifndef jam2_meanfield_MeanField_h
#define jam2_meanfield_MeanField_h

#include <Pythia8/Basics.h>
#include <Pythia8/Settings.h>
#include <jam2/collision/EventParticle.h>
#include <jam2/meanfield/TwoBodyDistance.h>

namespace jam2 {

class MeanField
{
protected:
  Pythia8::Settings* settings;

  int NV, selfInt;
  bool withMomDep;
  int eosType,optPotential, isDebug, optP0dev, optPV;
  int optMDarg3;
  int optPotentialArg;
  int optTwoBodyDistance,optScalarDensity, optVectorDensity,optVectorPotential;
  int optPotentialDensity;
  double facMesonPot;
  int optLambdaPot;
  Pythia8::Vec4 pTot0, pTot;
  int optRecoverEnergy;
  int optVdot, optBaryonCurrent;
  bool optDerivative;
  bool firstEng;
  static bool firstCall;
  double gamCM, gamCM2;
  double globalTime;
  TwoBodyDistance *distance, *distanceP;
  double Alpha,Beta,Gam,rho0,Vex1,Vex2,pMu1,pMu2;

public:
  MeanField(Pythia8::Settings* s);
  virtual ~MeanField() {delete distance;delete distanceP;}
  virtual void evolution(std::list<EventParticle*>& plist,double t, double dt,int step)=0;
  //virtual void singleParticlePotential()=0;
  virtual void init(std::list<EventParticle*>& plist)=0;
  virtual Pythia8::Vec4 computeEnergy(std::list<EventParticle*>& plist,int step)=0;
};

}
#endif
