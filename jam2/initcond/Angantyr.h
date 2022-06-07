#ifndef jam2_initcond_Angantyr_h
#define jam2_initcond_Angantyr_h

#include <cmath>
#include <string>
#include <jam2/collision/Collision.h>
#include <jam2/initcond/InitialCondition.h>
#include <jam2/hadrons/JamParticleData.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/HIUserHooks.h>

namespace jam2 {

class Angantyr: public InitialCondition
{
protected:
  Pythia *pythia;
  std::ofstream ofsPS;
  bool optOutPhaseSpace;

 class ImpactParamGenerator : public ImpactParameterGenerator {
   private:
    double bMin, bMax;
    Pythia8::Rndm *rndm;
   public:
   ImpactParamGenerator(double bmin, double bmax,Pythia8::Rndm* r)
    :bMin(bmin),bMax(bmax),rndm(r) {}
   bool init() {return true;}
   Pythia8::Vec4  generate(double & weight) const {
     double b=sqrt(bMin*bMin + rndm->flat()*(bMax*bMax-bMin*bMin));
     return Pythia8::Vec4(b,0.0, 0.0, 0.0);
   };
 };

 class ImpactParameterHook : public HIUserHooks {
   private:
     ImpactParamGenerator* bGen; 
   public:
   ImpactParameterHook(ImpactParamGenerator* ipg): bGen(ipg) {}
   bool hasImpactParameterGenerator() const { return true; }
   ImpactParameterGenerator * impactParameterGenerator() const {return bGen;}
 };

  ImpactParameterHook * impactParHook;
  ImpactParamGenerator * bGen;

public:
  Angantyr(Pythia8::Settings* s, JamParticleData* pd,Pythia8::Rndm* r);
  virtual~Angantyr();
  void init();
  void generate(Collision* event,int mode=0);

};

}

#endif
