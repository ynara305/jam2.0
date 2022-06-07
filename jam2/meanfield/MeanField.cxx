#include <jam2/meanfield/MeanField.h>
#include <jam2/hadrons/JamStdlib.h>

namespace jam2 {

using namespace std;

bool MeanField::firstCall=true;

MeanField::MeanField(Pythia8::Settings* s)
{
  settings = s;

  firstEng=true;
  //dT= settings->parm("Cascade:TimeStepSize");
  //widG = settings->parm("MeanField:gaussWidth");
  //facG = 1.0/pow(4.0*M_PI*widG, 1.5);
  //wG = 1.0/(4*widG);

  eosType=settings->mode("MeanField:EoS");
  optPotential=settings->mode("MeanField:optPotential");
  isDebug = settings->mode("Check:Debug");
  int self = settings->mode("MeanField:selfInteraction");
  selfInt = 1 - self;

  optScalarDensity=settings->mode("MeanField:optScalarDensity");
  optP0dev=settings->mode("MeanField:optP0dev");
  optTwoBodyDistance=settings->mode("MeanField:twoBodyDistance");
  optVectorPotential=settings->mode("MeanField:optVectorPotential");
  optVectorDensity=settings->mode("MeanField:optVectorDensity");
  optPV=settings->mode("MeanField:optMomPotential");
  optPotentialArg=settings->mode("MeanField:optPotentialArg");
  optRecoverEnergy=settings->mode("MeanField:optRecoverEnergy");
  gamCM = settings->parm("Cascade:gamCM");
  gamCM2 = gamCM*gamCM;
  optDerivative = settings->flag("MeanField:optDerivative");
  facMesonPot = settings->parm("MeanField:factorMesonPotential");

  rho0=settings->parm("MeanField:rho0");
  optLambdaPot = settings->mode("MeanField:optLambdaPotential");
  optPotentialDensity = settings->mode("MeanField:optPotentialDensity");

  // RQMDv
  optVdot=settings->mode("MeanField:optVdot");
  // In case only time component of vector potential is included,
  // we do not need to distinguish between kinetic momenta and canonical momenta
  if(optVectorPotential !=1) optVdot=0;
    
  if(optTwoBodyDistance==1) {
    distance = new NonRelDistance(settings);
  } else if(optTwoBodyDistance==2) {
    distance = new TwoBodyCM(settings);
  } else {
    distance = new RestFrame(settings);
  }

  optMDarg3=settings->mode("MeanField:twoBodyDistanceMD");
  if(optMDarg3==1) {
    distanceP = new NonRelDistance(settings);
  } else if(optMDarg3==2) {
    distanceP = new TwoBodyCM(settings);
  } else {
    distanceP = new RestFrame(settings);
  }

}

} // namespace jam2


