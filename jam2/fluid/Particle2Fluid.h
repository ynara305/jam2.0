#ifndef jam2_fluid_Particle2Fluid_h
#define jam2_fluid_Particle2Fluid_h

#include <list>
#include <jam2/collision/EventParticle.h>
#include <jam2/collision/Collision.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/fluid/Fluid.h>
#include <Pythia8/Basics.h>
#include <Pythia8/Settings.h>

namespace jam2 {

using Pythia8::Vec4;

class Particle2Fluid
{
private:
  Pythia8::Settings *settings;
  JamParticleData* particleTable;
  std::list<EventParticle*>* particleList;
  Fluid *fluid;
  int optGauss, optConv, optDens, optPreHadron, optEdmin, optCoreCorona;
  int optVectorPotential,optVdot,optPropagate;
  int n_gaussp, optTauEta, isDebug;
  int maxX, maxY, maxZ, iP, iQ, iR;
  double minEnergyDensity, widG, widG2,widCof, gVolume;
  double dX, dY, dZ, dCut, volume;
public:
  Particle2Fluid(Pythia8::Settings *s, Fluid* f, JamParticleData* d);
  void setParticleList(list<EventParticle*>& plist) {particleList=&plist;}
  int  checkFluidConversion(EventParticle* ip,
	std::vector<EventParticle*>& outgoing,double gtime);
  void checkFluidConversion(InterList* inter,
	std::vector<EventParticle*>& outgoing,double gtime);
  bool checkEnergyDensity(EventParticle* ip,double ctime,
	  Vec4& r, int& x, int& y, int& z);
  double edminfluid(double binv);
  void convert(Collision* event,EventParticle* ip,double ctime);
  void convertAform(Collision* event,EventParticle* ip,std::vector<EventParticle*>& outgoing,
	  double ctime);
  void p2fluid(EventParticle* ep, Vec4& r, int ix, int iy, int iz,int opt=1);
  int  dens(EventParticle* i1,double ctime,
	   double& rho,double& rhob,double& einv,int iopt);
  double gaussSmear(const Vec4& pv, const Vec4& dr);
  void convertAtDenseRegion(double ctime,Collision* event);
  void convertAll(double ctime,double ftime,Collision* event);
};
} // end namespace jam2
#endif
