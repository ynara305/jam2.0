#include <cstdlib>
#include <jam2/initcond/Angantyr.h>
#include <jam2/initcond/Nucleus.h>
#include <jam2/collision/Collision.h>

// You need to include this to get access to the HIInfo object for
// HeavyIons.
#include <Pythia8/HeavyIons.h>


using namespace std;
using namespace Pythia8;

namespace jam2 {

Angantyr::Angantyr(Settings* s, JamParticleData* p, Rndm* r)
    : InitialCondition(s,p,r)
{
  pythia = new Pythia(*settings,*jamParticleData->getParticleData(),false);
}

Angantyr::~Angantyr()
{
  delete pythia;
  //delete bGen;
  delete impactParHook;

  if(optOutPhaseSpace) ofsPS.close();
}

void Angantyr::init()
{
  optOutPhaseSpace=settings->flag("Cascade:outputInitialCondition");
  if(optOutPhaseSpace) ofsPS.open("AngantyrInitCond.dat");

  eCM=settings->parm("Beams:eCM");
  pLab=settings->parm("Beams:pLab");
  eLab=settings->parm("Beams:eLab");
  bMin=settings->parm("Beams:bmin");
  bMax=settings->parm("Beams:bmax");
  zSep=settings->parm("Beams:zseparation");
  //compFrame=settings->word("Beams:compFrame");

  string nuclA = settings->word("Beams:beamA");
  string nuclB = settings->word("Beams:beamB");
  int idProj=findAZ(nuclA);
  int idTarg=findAZ(nuclB);

  pythia->settings.mode("Beams:idA",idProj);
  pythia->settings.mode("Beams:idB",idTarg);
  pythia->settings.parm("Beams:eCM",eCM);
  pythia->settings.mode("Beams:frameType",1);
  pythia->readString("HadronLevel:Decay = off");
  pythia->readString("Fragmentation:setVertices = on");

  pythia->readString("111:mayDecay = off");   // pi0


  // Initialize the Angantyr model to fit the total and semi-inclusive
  // cross sections in Pythia within some tolerance.
  pythia->readString("HeavyIon:SigFitErr = "
                    "0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");

  // These parameters are typicall suitable for sqrt(S_NN)=5TeV
  pythia->readString("HeavyIon:SigFitDefPar = "
                    "17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0");

  // A simple genetic algorithm is run for 20 generations to fit the
  // parameters.
  pythia->readString("HeavyIon:SigFitNGen = 20");

  // add HIHooks.
  bGen = new ImpactParamGenerator(bMin,bMax,rndm);
  impactParHook = new ImpactParameterHook(bGen);
  pythia->setHIHooks(impactParHook);

  // assume both are nucleous
  mA=0.938;
  mB=0.938;

  pzAcm  = 0.5 * sqrtpos( (eCM + mA + mB) * (eCM - mA - mB)
           * (eCM - mA + mB) * (eCM + mA - mB) ) / eCM;
  eA = sqrt(mA*mA + pzAcm*pzAcm);
  eB = sqrt(mB*mB + pzAcm*pzAcm);
  yCM = 0.5*log((eA + pzAcm)/(eA - pzAcm));
  gamA= eA/mA;
  gamB= eB/mB;
  masA = atoi(nuclA.c_str());
  masB = atoi(nuclB.c_str());

  cout << " nuclA= " << nuclA << " nuclB= "<< nuclB<<endl;
  cout << "massA= "<< masA << " masB= "<< masB<<endl;

  // Initialise Pythia.
  pythia->init();

}

void Angantyr::generate(Collision* event, int mode)
{
  int ntry=0;
  do {
    if ( pythia->next() ) break;
    if(ntry==100) {
      cout << " Angantyr event generation failed 100 times."<<endl;
      exit(1);
    }
  }while(ntry++ < 100);

  // Generate impact parameter.
  //impactPar = setImpactPar();

  HeavyIons* heavyIons = pythia->heavyIonsPtr;
  HIInfo  hiinfo = heavyIons->hiinfo;
  impactPar = hiinfo.b();
  nCollG = hiinfo.nCollTot()+hiinfo.nCollND()
      +hiinfo.nCollSDP()+hiinfo.nCollSDT()
      +hiinfo.nCollDD()+hiinfo.nCollCD()+hiinfo.nCollEL();

  int nPartProj = hiinfo.nPartProj();
  int nPartTarg = hiinfo.nPartTarg();
  nPartG = nPartProj + nPartTarg;
  cout << "Angantyr ncoll = "<< nCollG << " npart= "<< nPartG <<endl;

  //HIUserHooks* hiHooks = pythia->hiHooksPtr;
  //cout << " hihooks = "<< hiHooks << endl;
  //cout << " impact par = "<< impactPar << endl;

  if(optOutPhaseSpace) ofsPS << pythia->event.size() <<endl;

  // put pythia particle into jam particle list.
  for(int i=0; i<pythia->event.size();i++) {
    Particle & py = pythia->event[i];
    if(!py.isFinal()) continue;

    int id=py.id();
    double m =  py.m();
    Pythia8::Vec4 p =  py.p();
    Pythia8::Vec4 r =  py.vProd()*MM2FM;
    ParticleDataEntry* pa= &py.particleDataEntry();
    EventParticle* cp = new EventParticle(id,m,r,p,pa);
    //cout << "id = "<< id << " name= "<< pa->name() << " m= "<< m<<endl;
    //if(id > 1000000000) continue; // exclude remnant.
    if(id > 1000000000) {
	//cout << " id= "<< id <<endl;
	continue; // exclude remnant.
    }

    cp->setPID(jamParticleData->pid(abs(id)));
    cp->setOnShell();
    double e=cp->getE0();
    double dect = jamParticleData->lifeTime(pa,m,e);
    double t =  py.tProd()*MM2FM;
    cp->setLifeTime(t+dect);
    event->setPList(cp);

    if(optOutPhaseSpace) ofsPS << setw(8) << id
        << setw(12) << fixed << cp->getMass()
        << scientific
        << setw(16) << setprecision(8) << cp->getPx()
        << setw(16) << setprecision(8) << cp->getPy()
        << setw(16) << setprecision(8) << cp->getPz()
        << setw(16) << setprecision(8) << cp->getPe()
        << setw(16) << setprecision(8) << cp->getX()
        << setw(16) << setprecision(8) << cp->getY()
        << setw(16) << setprecision(8) << cp->getZ()
        << setw(16) << setprecision(8) << cp->getT()
        << endl;

  }

}

}// end namespace jam2
