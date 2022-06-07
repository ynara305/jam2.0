// Copyright (C) 2018-2020 Yasushi Nara.
// This file contains the main class for event generation of JAM.

#ifndef jam2_JAM_H
#define jam2_JAM_H

#include <Pythia8/Pythia.h>
#include <Pythia8/SigmaTotal.h>
#include <Pythia8/Basics.h>
#include <Pythia8/FragmentationFlavZpT.h>
//#include <Pythia8/SimpleTimeShower.h>
#include <jam2/collision/EventParticle.h>
#include <jam2/collision/Collision.h>
#include <jam2/initcond/InitialCondition.h>
#include <jam2/initcond/BoostedTwoNuclei.h>
#include <jam2/initcond/Angantyr.h>
#include <jam2/xsection/CrossSection.h>
#include <jam2/interaction/Scatter.h>
#include <jam2/collision/InterList.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/interaction/JPythia.h>
#include <jam2/fluid/FluidHandler.h>
#include <jam2/meanfield/MeanField.h>
#include <jam2/interaction/NuclearCluster.h>
#include <jam2/JamAna.h>

namespace jam2 {

using namespace Pythia8;

class JAM
{
public:
    // Constants: could only be changed in the code itself.
    static const char* const VERSIONSTRING;
    static const double VERSIONNUMBERHEAD, VERSIONNUMBERCODE;

    Pythia          *pythia, *hadronize, *pydecay;
    Settings        *settings;
    Info            *info;
    ParticleData    *particleData;
    Rndm            *rndm;
    Event           pyEvent;
    JamParticleData* jamParticleData;
    CrossSection* xsection;

private:
    JPythia *jpythia=0, *jpythia2=0;
    Settings        defaultSettings;
    Angantyr* angantyr;

    //JBeamShape* jBeamShape;
    InitialCondition* initcnd = 0;
    BoostedTwoNuclei*  iniTwo = 0;
    HadronDecay* decay = 0;
    Collision* event;
    Scatter* scatt;
    FluidHandler* fluidHandler;
    CollisionHistory *anaColl=0;
    AnaTimeDepParticle   *anaParticle=0;
    AnaTimeDepFlow   *anaFlow=0;
    AnaTimeDepDensity  *anaDens=0;
    AnaMeanField *anaPot=0;
    AnaOutPutPhaseSpace *anaPhase=0;
    NuclearCluster *nuclearCluster=0;

    Vec4 pTot, pTot0;
    int nStep; // total number of time step.
    int cascadeMethod, collisionUpdateFreq;
    int isDebug,printColl,optPrintDisplay;
    int numInter; // initially expected number of collisions.
    int numPart;  // initially expected number of participants.
    int overSample;
    double impactPar=0.0,bMax; // impact parameter of an event.
    double eCM,xColl,xCollBB,xCollMB,xCollMM,xCollBBar,xInter,xPauli,xDec,xCollRR2NN;
    double xElastic,xAbsorb,xCollBarBar;
    int doPauliBlock;
    int nEvent, iEvent;
    int nColl; // number of collision
    int nPrint, nDec,ncollBB, ncollMB, ncollMM,ncollBBar,ncollBarBar,ncollRR2NN;
    int nElastic,nAbsorb,nPauli;
    int printFreq,printFreqPhaseSpace;
    bool isRandomInitialized = false;
    bool doNuclearClusterFormation = false;

    bool withHydro;
    int hydroInitCond,hydroDynamicalFreezeOut;
    double hydroStartTime, hydroSwitchTime;

    int withMeanField,optVectorPotential,optVdot,optPropagate;
    MeanField* meanField = 0;

public:
  // Constructor. (See Pythia.cc file.)
  JAM(string fname = "jam.inp", string xmlDir = "../share/Pythia8/xmldoc", bool printBanner = false);
  //: Pythia(xmlDir,printBanner) {;}

  ~JAM();

  // Initialize.
  bool init(InitialCondition* init=0);

  // Generate the next event.
  bool next();
  bool next0();

  double cascade(double initime,double ftime);

  // some weak decay of harons off.
  void mayHadronDecay(Pythia8::Pythia* py);

  void append(EventParticle* p)  {event->plist.push_back(p);}
  std::list<EventParticle*>& getEvent() {return event->plist;}
  void collisionInfor(InterList* inter);
  void finalDecay();
  void finalDecay1();
  void printCollision(double ctime,int operation,InterList* inter, std::vector<EventParticle*> outgoing);
  void printDisplay(double ctime);
  void display(double tnow);
  void collisionStatistics(double ctime,InterList* inter,std::vector<EventParticle*> outgoing);
  void collisionStatistics2(double ftime);
  void checkEnergyMomentumConservation(InterList* inter,
	  std::vector<EventParticle*>& outgoing);
  void computeTotalEnergy(int step,double dt);
  void initTimeDependentAna(const double gamcm,const double ycm);
  void eventInitTimeDependentAna();
  void anaTimeDependent(const int step, const double timeNow);
  void statTimeDependentAna(const double coltime,InterList* inter,vector<EventParticle*> outgoing);
  void printTimeDependentAna();
  void finTimeDependentAna();

 // Read in settings values: shorthand, not new functionality.
  bool   flag(std::string key) {return settings->flag(key);}
  int    mode(std::string key) {return settings->mode(key);}
  double parm(std::string key) {return settings->parm(key);}
  string word(std::string key) {return settings->word(key);}
  Pythia* getHadronize() {return hadronize;}
  int    getOverSample() const {return overSample;}
  int    getGlauberColl() const {return numInter/overSample;}
  int    getGlauberParticipants() const {return numPart/overSample;}
  int    getNColl() const {return nColl/overSample;}
  int    getNCollBB() const {return ncollBB/overSample;}
  int    getNCollRR2NN() const {return ncollRR2NN/overSample;}
  double getXColl() const {return xColl/iEvent/overSample;}
  double getXCollBB() const {return xCollBB/iEvent/overSample;}
  double getXCollMB() const {return xCollMB/iEvent/overSample;}
  double getXCollMM() const {return xCollMM/iEvent/overSample;}
  double getXCollBBar() const {return xCollBBar/iEvent/overSample;}
  double getXCollBbarBar() const {return xCollBarBar/iEvent/overSample;}
  double getXCollRR2NN() const {return xCollRR2NN/iEvent/overSample;}
  double getNInter() const {return xInter/iEvent/overSample;}
  double getNElastic()const  {return xElastic/iEvent/overSample;}
  double getNAbsorb() const {return xAbsorb/iEvent/overSample;}
  double getNPauli() const {return xPauli/iEvent/overSample;}
  double getNDecay() const {return xDec/iEvent/overSample;}
  double impactParameter() const {return impactPar;}
  InitialCondition*  initcond() {return initcnd;}

};
}
#endif // jam2_JAM_H
