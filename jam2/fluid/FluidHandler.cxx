#include <jam2/fluid/FluidHandler.h>

namespace jam2 {

using namespace std;

FluidHandler::FluidHandler(Pythia8::Settings* s,JamParticleData* jdata,
 HadronDecay* d, Rndm* rnd)
{
  settings = s;
  jamParticleData = jdata;
  decay = d;
  rndm=rnd;

  hydroInitCond = settings->mode("Hydro:InitialCondition");
  hydroStartTime = settings->parm("Hydro:StartTime");
  optConvertParticle = settings->mode("Hydro:optConvertParticle");
  eCM=settings->parm("Beams:eCM");

  hydroDynamicalFreezeOut=0;
  if(settings->mode("Hydro:optFreezeOut")==11 || settings->mode("Hydro:optFreezeOut")==12) 
    hydroDynamicalFreezeOut=1;

  dt = settings->parm("Cascade:TimeStepSize");
  int collisionUpdateFreq=settings->mode("Cascade:CollisionUpdateFreq");
  dtc = dt*collisionUpdateFreq;

  int islog=0;
  eos = new EoS(settings->word("Hydro:EoSFileName"),islog);
  freezeout =  new FreezeOut(settings,eos);
  fluid = new Fluid(settings,eos,freezeout);
  p2fluid = new Particle2Fluid(settings,fluid,jamParticleData);
  fluid2p = new Fluid2Particle(settings,jamParticleData,freezeout,
	  fluid,rndm);

  printFluidFraction=settings->flag("Hydro:printFluidFraction");
  //if(printFluidFrac) outputFluc.open("fluidfrac.dat",std::ios::app);
  if(printFluidFraction) outputFluc.open("fluidfrac.dat");
  optHadronCascade=settings->flag("Hydro:optHadronCascade");

}

FluidHandler::~FluidHandler()
{
  delete eos;
  delete freezeout;
  delete fluid;
  delete p2fluid;
  delete fluid2p;
  if(printFluidFraction) outputFluc.close();
}

int FluidHandler::init(Collision* event, BoostedTwoNuclei* iniTwo)
{
  double tpass = event->nuclearPassingTime();
  fluid->reset();
  fluid->setPassTime(tpass);
  p2fluid->setParticleList(event->plist);

  int hyswitch=0;

  // hydro will start after two nuclei pass through each other.
  if(hydroInitCond==1) {

    hyswitch = tpass / dt;
    fluid->seInitialTime(hyswitch*dt);

  // hydro will start at the time user specified.
  } else if(hydroInitCond==2) {

    hyswitch = max(1, (int)(hydroStartTime / dt));
    fluid->seInitialTime(hyswitch*dt);

    // Option for smooth initial condition of colliding two nuclei.
    if(hydroStartTime < 0.0) {
	if(iniTwo==0) {
	    cout << "FluidHandler iniTwo is not initizlized"<< endl;
	    exit(1);
	}
        fluid->seInitialTime(0.0);
        double A = iniTwo->getAproj();
        double B = iniTwo->getAtarg();
	double impactPar = iniTwo->getImpactPar();
        tpass=fluid->makeTwoNuclei(eCM,A,B,impactPar,1);
        fluid->setPassTime(tpass);
        hyswitch=0;
        // delete particles.
        auto plist = &event->plist;
        for(auto ip=plist->begin();ip !=plist->end();ip++) delete *ip;
        plist->clear();
        event->makeCollisionList(0.0,0.0);
    }

  }
  //cout << " pass time = " << tpass << " switching step= "<< hyswitch << endl;


  return hyswitch;

}

void FluidHandler::convertAll(Collision* event,double step)
{
  hydroSwitchTime = dt*(step-1);
  p2fluid->convertAll(hydroSwitchTime,dt*step,event);
  fluid->resetFreezeOut();
}

bool FluidHandler::evolution(Collision* event,int step,double gtime, int& noCollUpdate)
{

  if(printFluidFraction) fluid->fluidFraction(gtime,event,outputFluc);

  // Skip if there is no fluid element.
  if(fluid->isFluid() <= 0) return true;
  if(optConvertParticle) p2fluid->convertAtDenseRegion(gtime,event);

  int check=fluid->evolution(step,gtime,event);

  // check=0: when all surface is below the critical energy density,
  // convert all fluid into particles all at once.
  if(check==0 || hydroDynamicalFreezeOut) {
    fluid2p->convert(dt*step,check);
    if(fluid2p->outgoing.size() > 0) {
      //double tnext=dtc*(step+1);
      //event->collisionUpdate(fluid2p->outgoing,hydroSwitchTime,tnext);
      event->setPList(fluid2p->outgoing);
      noCollUpdate=0;
    }
    if(check==0) {
      fluid->setNC(-1);
      if(printFluidFraction && hydroDynamicalFreezeOut==0)
	fluid->fluidFraction(gtime,event,outputFluc,1);
    }
  }

  if(check==0 && optHadronCascade==false) return false;
  return true;

}

bool FluidHandler::finalize(Collision* event,int step,double gtime)
{
  // Skip if there is no fluid element.
  if(fluid->isFluid() <= 0) return true;

  if(hydroDynamicalFreezeOut) {
    fluid2p->convert(dt*step,0);
    if(fluid2p->outgoing.size() > 0) {
      //double tnext=dtc*(step+1);
      //event->collisionUpdate(fluid2p->outgoing,hydroSwitchTime,tnext);
      event->setPList(fluid2p->outgoing);
    }
  }

  return true;

}

} // namespace jam2


