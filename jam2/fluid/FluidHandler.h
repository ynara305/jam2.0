#ifndef jam2_fluid_FluidHandler_h
#define jam2_fluid_FluidHandler_h

#include <vector>
#include <Pythia8/Basics.h>
#include <Pythia8/Settings.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/collision/EventParticle.h>
#include <jam2/collision/InterList.h>
#include <jam2/fluid/Fluid.h>
#include <jam2/fluid/FluidElement.h>
#include <jam2/fluid/FreezeOut.h>
#include <jam2/fluid/Particle2Fluid.h>
#include <jam2/fluid/Fluid2Particle.h>
#include <jam2/initcond/BoostedTwoNuclei.h>

namespace jam2 {

class FluidHandler
{
private:
  Pythia8::Settings* settings;
  JamParticleData* jamParticleData;
  HadronDecay* decay;
  EoS* eos;
  Fluid* fluid;
  FreezeOut *freezeout;
  Particle2Fluid* p2fluid;
  Fluid2Particle* fluid2p;
  Pythia8::Rndm *rndm;
  double eCM,dt,dtc;
  int hydroInitCond,optConvertParticle;
  bool optHadronCascade;
  double hydroStartTime, hydroSwitchTime,hydroDynamicalFreezeOut;
  bool printFluidFraction;
  std::ofstream outputFluc;

public:
  FluidHandler(Pythia8::Settings* s,JamParticleData* jdata, HadronDecay* d,
	  Pythia8::Rndm* rnd);

  ~FluidHandler();
   //int isFluid() {return fluid->isFluid();}

  int init(Collision* event, BoostedTwoNuclei* iniTwo=0);
  void convertAll(Collision* event,double step);
  bool evolution(Collision* event,int step,double ftime,int& noCollUpdate);
  bool finalize(Collision* event,int step,double gtime);

  void checkFluidConversion(InterList* inter,
	std::vector<EventParticle*>& outgoing,double ftime) {
    p2fluid->checkFluidConversion(inter,outgoing,ftime);
  }

  int  checkFluidConversion(EventParticle* ip,
	std::vector<EventParticle*>& outgoing,double ftime) {
    return p2fluid->checkFluidConversion(ip,outgoing,ftime);
  }

  void convertAform(Collision* event,EventParticle* ep,std::vector<EventParticle*>& outgoing,
	  double ctime) {
    p2fluid->convertAform(event,ep,outgoing,ctime);
  }

  void convert(Collision* event,EventParticle* ip,double finaltime) {
    p2fluid->convert(event,ip,finaltime);
    delete ip;
   }

  Pythia8::Vec4 pFluidTotal() {return fluid->pFluidTotal();}

};
} // end namespace jam2
#endif
