// Copyright (C) 2020 Yasushi Nara
// JAM is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include "config.h" /* for PACKAGE_VERSION */
#include <cstdlib>
#include <jam2/JAM.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/hadrons/JamStdlib.h>
#include <jam2/interaction/JPythia.h>
#include <jam2/collision/Collision1.h>
#include <jam2/collision/Collision3.h>
#include <jam2/meanfield/GBUU.h>
#include <jam2/meanfield/QMD.h>
#include <jam2/meanfield/RQMDs0.h>
#include <jam2/meanfield/RQMDs.h>
#include <jam2/meanfield/RQMDv.h>
#include <jam2/meanfield/RQMDsv.h>
#include <jam2/meanfield/RQMDw.h>
#include <jam2/meanfield/QMDBox.h>
#include <jam2/meanfield/RQMDsBox.h>
#include <jam2/meanfield/RQMDvBox.h>
#include <jam2/meanfield/RQMDsvBox.h>

//#include <Pythia8/ResonanceWidths.h>


namespace jam2 {

using namespace std;

// The current JAM (sub)version number, to agree with XML version.
const char* const JAM::VERSIONSTRING = PACKAGE_VERSION;
const double JAM::VERSIONNUMBERHEAD = std::atof(PACKAGE_VERSION);
const double JAM::VERSIONNUMBERCODE = 0.000;

JAM::JAM(string fname, string xmldir, bool printBanner)
{
  pythia = new Pythia(xmldir,printBanner);
  particleData = &pythia->particleData;
  settings =     &pythia->settings;
  rndm =         &pythia->rndm;
  info =         &pythia->info;

  // Add JAM default parameters and options.
  pythia->settings.addWord("Beams:beamA","197Au");
  pythia->settings.addWord("Beams:beamB","197Au");
  pythia->settings.addWord("Beams:compFrame","nn");
  pythia->settings.addParm("Beams:bmin",0.0,true,false,0.0,15.0);
  pythia->settings.addParm("Beams:bmax",2.0,true,false,0.0,15.0);
  pythia->settings.addParm("Beams:pLab",0.0,true,false,0.0,100000.0);
  pythia->settings.addParm("Beams:eLab",0.0,true,false,0.0,100000.0);
  pythia->settings.addParm("BeamRemnants:primordialKTremnantSoft",0.6,true,false,0.0,10.0);
  pythia->settings.addParm("BeamRemnants:valenceDiqEnhanceSoft",2.0,true,false,0.5,10.0);
  pythia->settings.addParm("BeamRemnants:valenceDiqEnhance",2.0,true,false,0.5,10.0);
  pythia->settings.addParm("BeamRemnants:probDiffra",0.1,true,false,0.0,1.0);
  pythia->settings.addParm("BeamRemnants:probFlavorExchange",0.8,true,false,0.0,1.0);

  pythia->settings.addParm("StringPT:sigmaBBar",0.335,true,false,0.0,1.0);

  pythia->settings.addParm("ParticleDecays:utRatio",0.0,true,false,0.0,1.0);
  pythia->settings.addParm("ParticleDecays:utRatioS",0.0,true,false,0.0,1.0);
  pythia->settings.addParm("ParticleDecays:gWidth",0.4,true,false,0.0,10.0);
  pythia->settings.addMode("ParticleDecays:potentialHandling",0,true,false,0,1);
  pythia->settings.addMode("ParticleDecays:optDecayAngle",1,true,false,0,3);
  pythia->settings.addMode("ParticleDecays:optDecayAngleSchannel",1,true,false,0,3);
  pythia->settings.addMode("ParticleDecays:optDeltaDecayAngle",2,true,false,0,3);

  pythia->settings.addParm("Beams:eCM",10.0,true,false,0.3,100000.0);

  pythia->settings.addParm("Beams:zseparation",1.5,true,false,-2.0,15.0);
  pythia->settings.addMode("Check:Debug",0,true,false,0,10);

  pythia->settings.addMode("Model:NuclearCollision",2,true,false,0,10);
  pythia->settings.addMode("Collision:SoftModel",2,true,false,1,2);

  pythia->settings.addMode("Cascade:model",3,true,false,0,3);
  pythia->settings.addFlag("Cascade:outputInitialCondition",false);
  pythia->settings.addMode("Cascade:overSample",1,true,false,1,10000);
  pythia->settings.addMode("Cascade:initialCondition",1,true,false,1,3);
  pythia->settings.addMode("Cascade:PrintCollision",0,true,false,0,10);
  pythia->settings.addMode("Cascade:PrintDisplay",0,true,false,0,10);
  pythia->settings.addParm("Cascade:displayScale",1.0,true,false,0.1,10.0);
  pythia->settings.addParm("Cascade:timeStart",0.0,true,false,-100,100.0);
  pythia->settings.addFlag("Cascade:nuclearClusterFormation",false);
  pythia->settings.addParm("Cascade:nuclearClusterR",3.8,true,false,0.0,10.0);
  pythia->settings.addParm("Cascade:nuclearClusterP",0.3,true,false,0.0,10.0);

  pythia->settings.addFlag("Cascade:InelasticOnly",false);
  pythia->settings.addFlag("Cascade:Decay",true);
  pythia->settings.addFlag("Cascade:finalDecay",true);
  pythia->settings.addFlag("Cascade:constQCollisionOnly",false);
  pythia->settings.addFlag("Cascade:BBCollisionOnly",false);
  pythia->settings.addFlag("Cascade:noMMCollision",false);
  pythia->settings.addParm("Cascade:minCMEnergyPythia",10.0,true,false,4.0,100.0);
  pythia->settings.addMode("Cascade:optCollisionOrder",4,true,false,1,21);
  pythia->settings.addMode("Cascade:optCollisionTimeLimit",0,true,false,0,1);
  pythia->settings.addMode("Cascade:optConstQuarkScatt",1,true,false,0,10);
  pythia->settings.addMode("Cascade:optConstQuarkScattHard",3,true,false,0,10);
  pythia->settings.addMode("Cascade:optConstQuarkDiffra",5,true,false,0,10);
  pythia->settings.addMode("Cascade:optConstQSChannel",0,true,false,0,1);
  pythia->settings.addFlag("Cascade:optQuarkExchangeProcess",false);
  pythia->settings.addFlag("Cascade:optQuarkAnnihilationProcess",false);
  pythia->settings.addMode("Cascade:optElasticBackWardScattering",1,true,false,-1,2);

  pythia->settings.addParm("Cascade:GaussianWidth",1.0,true,false,0.1,3.0);
  pythia->settings.addMode("Cascade:PauliBlocking",2,true,false,0,2);
  pythia->settings.addMode("Cascade:collisionUpdateFreq",1,true,false,1,1000);
  pythia->settings.addMode("Cascade:noBaBAnnihilation",0,true,false,0,1);
  pythia->settings.addMode("Cascade:optFluctuation",0,true,false,0,1);
  pythia->settings.addFlag("Cascade:softRadiation",true);
  pythia->settings.addMode("Cascade:optDecayPosition",0,true,false,0,2);
  pythia->settings.addParm("Cascade:decayPositionSmearParam",0.5,true,false,0.0,3.0);
  pythia->settings.addMode("Cascade:optFermiMomentum",1,true,false,0,2);
  pythia->settings.addFlag("Cascade:removeSpectator",false);
  pythia->settings.addFlag("Cascade:BBabsorptionXS",true);
  pythia->settings.addParm("Cascade:ecmStringBB",3.5,true,false,2.0,5.0);
  pythia->settings.addParm("Cascade:ecmStringMB",3.0,true,false,1.0,5.0);
  pythia->settings.addParm("Cascade:ecmStringMBs",3.0,true,false,1.0,5.0);
  pythia->settings.addParm("Cascade:ecmStringMBc",4.0,true,false,1.0,5.0);
  pythia->settings.addParm("Cascade:ecmStringMM",3.0,true,false,1.0,5.0);
  pythia->settings.addParm("Cascade:ecmStringABB",3.5,true,false,1.0,5.0);
  pythia->settings.addMode("Cascade:allowRescatterSameString",0,true,false,0,2);

  pythia->settings.addMode("Cascade:optInitializeCellParameters",0,true,false,0,2);
  pythia->settings.addMode("Cascade:boxBoundary",0,true,false,0,1);
  pythia->settings.addParm("Cascade:vxCell",0.2,true,false,0.0,1.0);
  pythia->settings.addParm("Cascade:vyCell",0.2,true,false,0.0,1.0);
  pythia->settings.addParm("Cascade:vzCell",0.2,true,false,0.0,1.0);
  pythia->settings.addMode("Cascade:nxCell",10,true,false,1,50);
  pythia->settings.addMode("Cascade:nyCell",10,true,false,1,50);
  pythia->settings.addMode("Cascade:nzCell",10,true,false,1,50);

  double dx=sqrt(20.0/M_PI)+1e-6;
  double dy=dx;
  double dz=dx;
  double tstart =  std::min({dx,dy,dz});

  pythia->settings.addParm("Cascade:dxCell",dx,true,false,0.0,10.0);
  pythia->settings.addParm("Cascade:dyCell",dy,true,false,0.0,10.0);
  pythia->settings.addParm("Cascade:dzCell",dz,true,false,0.0,10.0);
  pythia->settings.addParm("Cascade:cellExpansionStartTime",tstart,true,false,0.0,10.0);

  pythia->settings.addMode("Cascade:box",0,true,false,0,1);
  pythia->settings.addParm("Cascade:boxLx",5.0,true,false,1.0,100.0);
  pythia->settings.addParm("Cascade:boxLy",5.0,true,false,1.0,100.0);
  pythia->settings.addParm("Cascade:boxLz",5.0,true,false,1.0,100.0);
  pythia->settings.addParm("Cascade:TimeStepSize",10000.0,true,false,0.01,5000000.0);
  pythia->settings.addMode("Cascade:TimeStep",1,true,false,1,5000);
  pythia->settings.addMode("Cascade:optBBJAM1",0,true,false,0,1);
  pythia->settings.addWord("Cascade:bwFileName1","BWintjam2a.dat");
  pythia->settings.addWord("Cascade:bwFileName2","BWintjam2b.dat");
  pythia->settings.addMode("Cascade:optSuppressSoftCollision",0,true,false,0,1);
  pythia->settings.addParm("Cascade:paramSoftCollision",1.0, true,false,0.0,10.0);
  pythia->settings.addFlag("Cascade:constQuarkDiffractive",false);

  pythia->settings.addFlag("Hydro:mode",false);
  pythia->settings.addMode("Hydro:InitialCondition",3,true,false,1,3);
  // =1: passing time =2: specified by Hydro:StartTime, =3: dynamical.
  pythia->settings.addMode("Hydro:optFluidConversion",3,true,false,1,3);
  pythia->settings.addMode("Hydro:optCoreCorona",1,true,false,0,3);
  pythia->settings.addParm("Cascade:gamCM",1.0,true,false,1.0,5000000.0);

  pythia->settings.addMode("Hydro:optConvertParticle",0,true,false,0,3);
  pythia->settings.addMode("Hydro:optFreezeOut",2,true,false,0,22);
  pythia->settings.addMode("Hydro:optTimeLike",1,true,false,0,2);
  pythia->settings.addParm("Hydro:StartTime",0.0,true,false,-100.0,100.0);
  pythia->settings.addWord("Hydro:EoSFileName","fluid/eosB235JAMsoft.dat");
  //pythia->settings.addWord("Hydro:EoSFileName","fluid/eosB235JAMsoft1mev.dat");
  pythia->settings.addMode("Hydro:nx",100,true,false,1,1000);
  pythia->settings.addMode("Hydro:ny",100,true,false,1,1000);
  pythia->settings.addMode("Hydro:nz",100,true,false,1,1000);
  pythia->settings.addParm("Hydro:dt",0.15,true,false,0.01,0.5);
  pythia->settings.addParm("Hydro:dx",0.3,true,false,0.01,0.5);
  pythia->settings.addParm("Hydro:dy",0.3,true,false,0.01,0.5);
  pythia->settings.addParm("Hydro:dz",0.3,true,false,0.01,0.5);
  pythia->settings.addParm("Hydro:gaussWidth",0.5,true,false,0.1,3.0);
  pythia->settings.addMode("Hydro:optGaussSmear",3,true,false,0,3);
  pythia->settings.addParm("Hydro:ParticlizationEnergyDensity",0.5,true,false,0.1,1.0);
  pythia->settings.addParm("Hydro:FluidizationEnergyDensity",0.5,true,false,0.1,1.0);
  pythia->settings.addFlag("Hydro:printFluidFraction",false);
  pythia->settings.addFlag("Hydro:optHadronCascade",true);

  pythia->settings.addMode("HeavyIon:optSample",1,true,false,1,10);
  pythia->settings.addFlag("HeavyIon:histWS",false);
  pythia->settings.addParm("HeavyIon:WSR0",0.03,true,false,0.0,1.5);

  pythia->settings.addMode("MeanField:mode",0,true,false,0,20);
  pythia->settings.addMode("MeanField:transportModel",1,true,false,1,2);
  pythia->settings.addMode("MeanField:optPotential",1,true,false,0,100);
  pythia->settings.addMode("MeanField:EoS",1,true,false,1,1000);
  pythia->settings.addMode("MeanField:optScalarDensity",0,true,false,0,2);
  pythia->settings.addMode("MeanField:optP0dev",1,true,false,0,1);
  pythia->settings.addMode("MeanField:twoBodyDistance",3,true,false,1,3);
  pythia->settings.addMode("MeanField:twoBodyDistanceMD",3,true,false,1,3);
  pythia->settings.addMode("MeanField:potentialType",2,true,false,1,3);

  pythia->settings.addParm("MeanField:gaussWidth",1.0,true,false,0.01,2.5);
  pythia->settings.addParm("MeanField:rho0",0.168,true,false,0.15,0.17);
  pythia->settings.addMode("MeanField:optVectorPotential",1,true,false,1,3);
  pythia->settings.addMode("MeanField:optVectorDensity",1,true,false,0,1);
  pythia->settings.addMode("MeanField:selfInteraction",0,true,false,0,1);
  pythia->settings.addMode("MeanField:optMomPotential",1,true,false,0,1);
  pythia->settings.addMode("MeanField:optVdot",0,true,false,0,4);
  pythia->settings.addMode("MeanField:optBaryonCurrent",0,true,false,0,1);
  pythia->settings.addMode("MeanField:optPotentialArg",2,true,false,0,1);
  pythia->settings.addFlag("MeanField:optDerivative",true);
  pythia->settings.addParm("MeanField:factorMesonPotential",0.333333,true,false,0.01,1.0);

  pythia->settings.addParm("MeanField:factorWidthL",1.0,true,false,0.0,10.0);
  pythia->settings.addParm("MeanField:factorAlphaPotentialL",1.0,true,false,0.0,10.0);
  pythia->settings.addParm("MeanField:factorBetaPotentialL",1.0,true,false,0.0,10.0);
  pythia->settings.addParm("MeanField:factorGammaPotentialL",1.0,true,false,0.0,10.0);
  pythia->settings.addParm("MeanField:factorAttMomPotentialL",1.0,true,false,0.0,10.0);
  pythia->settings.addParm("MeanField:factorRepMomPotentialL",1.0,true,false,0.0,10.0);
  pythia->settings.addParm("MeanField:factorWidthS",1.0,true,false,0.0,10.0);
  pythia->settings.addParm("MeanField:factorAlphaPotentialS",1.0,true,false,0.0,10.0);
  pythia->settings.addParm("MeanField:factorBetaPotentialS",1.0,true,false,0.0,10.0);
  pythia->settings.addParm("MeanField:factorGammaPotentialS",1.0,true,false,0.0,10.0);
  pythia->settings.addParm("MeanField:factorAttMomPotentialS",1.0,true,false,0.0,10.0);
  pythia->settings.addParm("MeanField:factorRepMomPotentialS",1.0,true,false,0.0,10.0);

  pythia->settings.addMode("MeanField:optStrangeBaryonPotential",0,true,false,0,4);
  pythia->settings.addMode("MeanField:optLambdaPotential",0,true,false,0,2);
  pythia->settings.addMode("MeanField:optPotentialDensity",0,true,false,0,1);

  pythia->settings.addMode("MeanField:optRecoverEnergy",0,true,false,0,3);
  pythia->settings.addParm("MeanField:stepVelocity",0.0,true,false,0.0,0.1);
  pythia->settings.addParm("MeanField:dtExpandStartTime",1.0,true,false,0.0,50.0);
  pythia->settings.addFlag("MeanField:gammaCorrectionGaussian",false);

  // mean-field cell
  pythia->settings.addMode("MeanField:optInitializeCellParameters",0,true,false,0,2);
  pythia->settings.addMode("MeanField:boxBoundary",0,true,false,0,1);
  pythia->settings.addParm("MeanField:vxCell",0.2,true,false,0.0,1.0);
  pythia->settings.addParm("MeanField:vyCell",0.2,true,false,0.0,1.0);
  pythia->settings.addParm("MeanField:vzCell",0.2,true,false,0.0,1.0);
  pythia->settings.addMode("MeanField:nxCell",10,true,false,1,50);
  pythia->settings.addMode("MeanField:nyCell",10,true,false,1,50);
  pythia->settings.addMode("MeanField:nzCell",10,true,false,1,50);
  pythia->settings.addParm("MeanField:dxCell",3.0,true,false,0.0,10.0);
  pythia->settings.addParm("MeanField:dyCell",3.0,true,false,0.0,10.0);
  pythia->settings.addParm("MeanField:dzCell",3.0,true,false,0.0,10.0);
  pythia->settings.addParm("MeanField:cellExpansionStartTime",3.0,true,false,0.0,10.0);

  pythia->settings.addFlag("Analysis:collision",false);
  pythia->settings.addMode("Analysis:printFreq",2,true,false,1,10);
  pythia->settings.addMode("Analysis:printFreqPhaseSpace",10,true,false,1,20);
  pythia->settings.addFlag("Analysis:timeDependenceParticle",false);
  pythia->settings.addFlag("Analysis:timeDependenceFlow",false);
  pythia->settings.addFlag("Analysis:timeDependenceDensity",false);
  pythia->settings.addFlag("Analysis:outPutPhaseSpace",false);
  pythia->settings.addFlag("Analysis:Potentials",false);
  pythia->settings.addParm("Analysis:yCut",0.5,true,false,0.0,100.0);
  pythia->settings.addParm("Analysis:yCutFoward",0.75,true,false,0.0,100.0);
  pythia->settings.addParm("Analysis:yCutMax",1.5,true,false,0.0,100.0);

  pythia->settings.flag("Random:setSeed",true);

  // JAM default settings.
  // Hadron Production Vertices.  See Hadron Scattering in the Pythia8 manual.
  pythia->readString("Fragmentation:setVertices = on");
  // 0: middle point,  1:late hadron production point,  -1:early
  pythia->readString("HadronVertex:mode = 0");
  pythia->settings.addParm("HadronVertex:kappa",1.0,true,false,0.01,100.0);
  pythia->readString("HadronVertex:kappa = 1.0");
  pythia->readString("HadronVertex:smearOn = on");
  pythia->readString("HadronVertex:xySmear = 0.7");
  pythia->readString("HadronVertex:constantTau = on");

  pythia->readString("Next:numberShowEvent = 0");// print event record n times
  pythia->readString("Next:numberShowInfo = 0");
  pythia->readString("Next:numberShowProcess = 0");
  pythia->readString("Print:quiet = on");
  pythia->readString("Main:timesAllowErrors = 100");

  pythia->settings.flag("Beams:allowVariableEnergy",true);
  pythia->settings.parm("Beams:eMinPert",10.0);
  pythia->settings.parm("Beams:eWidthPert",10.0);

  //pythia->settings.mode("Tune:pp",6);

  double mesonUDL1S0J1=pythia->settings.parm("StringFlav:mesonUDL1S0J1");
  double popcornRate=pythia->settings.parm("StringFlav:popcornRate");
  double probQQtoQ=pythia->settings.parm("StringFlav:probQQtoQ");
  double probStoUD=pythia->settings.parm("StringFlav:probStoUD");
  double ecmPow=pythia->settings.parm("MultipartonInteractions:ecmPow");

  pythia->settings.parm("StringFlav:mesonUDL1S0J1", 0.70);
  pythia->settings.parm("StringFlav:popcornRate",  0.15);
  //pythia->settings.parm("StringFlav:probQQtoQ", 0.05);
  pythia->settings.parm("StringFlav:probQQtoQ", 0.07);
  //pythia->settings.parm("StringFlav:probStoUD", 0.19);
  pythia->settings.parm("MultipartonInteractions:ecmPow", 0.177);

  cout << "# JAM version "<< PACKAGE_VERSION <<endl;

  // Read in commands from external file.
  pythia->readFile(fname);
  cout << "# input file name= " << fname << endl;

  defaultSettings = pythia->settings;
  defaultSettings.parm("StringFlav:mesonUDL1S0J1", mesonUDL1S0J1);
  defaultSettings.parm("StringFlav:popcornRate",  popcornRate);
  defaultSettings.parm("StringFlav:probQQtoQ", probQQtoQ);
  defaultSettings.parm("StringFlav:probStoUD", probStoUD);
  defaultSettings.parm("MultipartonInteractions:ecmPow", ecmPow);

  // change the minimum mass of the hadron resonances because
  // some default minimum mass is below the decay threshold.
  particleData->mMin(10333,1.39); // h_1(1380)
  //particleData->mMin(13122,1.34); // Lambda(1405)
  particleData->mMin(9010221,0.33); // f_0(980)
  particleData->mMax(9010221,1.2); // f_0(980)
  particleData->mWidth(9010221,0.05); // f_0(980)
  particleData->mayDecay(9010221,true); // f_0(980)

// Until Pythia8 fixed this....
//pythia8244/share/Pythia8/xmldoc/ParticleData.xml
//<particle id="9000211" name="a_0(980)+" antiName="a_0(980)-" spinType="1" chargeType="3" colType="0"
//          m0="0.98350" mWidth="0.06000" mMin="0.70000" mMax="1.50000">
// <channel onMode="1" bRatio="0.9000000" products="221 211"/>
// <channel onMode="1" bRatio="0.1000000" products="321 311"/>
//</particle>
//
// --->  <channel onMode="1" bRatio="0.1000000" products="321 -311"/>

  ParticleDataEntry* pd=particleData->findParticle(9000211);
  DecayChannel *dch = &pd->channel(1);
  dch->product(0,321); dch->product(1,-311);

  pythia->settings.flag("SoftQCD:inelastic", false);
  pythia->settings.flag("SoftQCD:all", false);
  pythia->settings.flag("SoftQCD:elastic", false);
  pythia->settings.flag("SoftQCD:nonDiffractive",false);
  pythia->settings.flag("SoftQCD:singleDiffractive", false);
  pythia->settings.flag("SoftQCD:doubleDiffractive", false);
  pythia->settings.flag("SoftQCD:centralDiffractive", false);

  //mayHadronDecay(pythia);
  //mayHadronDecay(hadronize);

  bMax=0.0;
  iniTwo=0;
  fluidHandler=0;
  event = 0;

  pd=particleData->findParticle(3122);
  cout << " 3122 mayDec= "<< pd->mayDecay() << " can= "<<  pd->canDecay()<<endl;

  pd=particleData->findParticle(3312);
  cout << " 3212 mayDec= "<< pd->mayDecay() << " can= "<<  pd->canDecay()<<endl;

  pd=particleData->findParticle(321);
  cout << " 321 mayDec= "<< pd->mayDecay() << " can= "<<  pd->canDecay()<<endl;

  pd=particleData->findParticle(311);
  cout << " 311 mayDec= "<< pd->mayDecay() << " can= "<<  pd->canDecay()<<endl;

  pd=particleData->findParticle(211);
  cout << " 211 mayDec= "<< pd->mayDecay() << " can= "<<  pd->canDecay()<<endl;

  pd=particleData->findParticle(111);
  cout << " 111 mayDec= "<< pd->mayDecay() << " can= "<<  pd->canDecay()<<endl;

  pd=particleData->findParticle(221);
  cout << " 221 mayDec= "<< pd->mayDecay() << " can= "<<  pd->canDecay()<<endl;

  // list all particle data.
  //particleData->listAll();
}

JAM::~JAM()
{
  delete initcnd;
  delete xsection;
  delete event;
  delete scatt;
  delete pydecay;
  delete decay;
  delete hadronize;
  delete pythia;
  delete jpythia;
  delete jpythia2;
  delete jamParticleData;
  delete fluidHandler;
  delete meanField;
  delete nuclearCluster;

  cout << "max. impact par. for two-body collision "<< bMax<<endl;

  finTimeDependentAna();
}

void JAM::mayHadronDecay(Pythia* py)
{
  py->readString("111:mayDecay = off");   // pi0
  py->readString("3122:mayDecay = off");  // Lambda decay
  py->readString("3212:mayDecay = off");  // Sigma0 decay
  py->readString("3112:mayDecay = off");  // Sigma- decay
  py->readString("3212:mayDecay = off");  // Sigma0 decay
  py->readString("3222:mayDecay = off");  // Sigma+ decay
  py->readString("3312:mayDecay = off");  // Xi- decay
  py->readString("3322:mayDecay = off");  // Xi0 decay
  py->readString("3334:mayDecay = off");  // Omega decay
  py->readString("130:mayDecay = off");   // K0_L
  py->readString("310:mayDecay = off");   // K0_S
  py->readString("311:mayDecay = off");   // K0
  py->readString("321:mayDecay = off");   // K+
  py->readString("221:mayDecay = off");   // eta
  py->readString("333:mayDecay = off");   // phi

}

bool JAM::init(InitialCondition* myinitcond)
{
  if (!isRandomInitialized) {
    isRandomInitialized = true;

    // Initialize the random number generator.
    if ( settings->flag("Random:setSeed") ) {
      rndm->init( settings->mode("Random:seed") );
      cout << " random seed = " << settings->mode("Random:seed") <<endl;
    }
  }

  // add additional hadronic resonances in the pythiaData table.
  jamParticleData = new JamParticleData(particleData,rndm);

  pydecay = new Pythia(pythia->settings,pythia->particleData,false);
  pydecay->readString("ProcessLevel:all = off");
  hadronize = new Pythia(pythia->settings,pythia->particleData,false);

  // Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
  hadronize->readString("ProcessLevel:all = off");
  hadronize->readString("HadronLevel:Decay = off");
  hadronize->readString("Fragmentation:setVertices = on");

  double Hkappa = hadronize->parm("HadronVertex:kappa");
  cout << " HadronVertex:kappa= " << Hkappa
       << " Vertex:mode= " << pythia->mode("HadronVertex:mode")
    <<endl;

  isDebug = pythia->settings.mode("Check:Debug");
  printColl=pythia->settings.mode("Cascade:PrintCollision");
  optPrintDisplay=pythia->settings.mode("Cascade:PrintDisplay");
  eCM=settings->parm("Beams:eCM");
  double pLab=settings->parm("Beams:pLab");
  double eLab=settings->parm("Beams:eLab");

  // find beam mass
  double mnuc=0.9383;
  double mA=mnuc, mB=mnuc;

  string beamA = settings->word("Beams:beamA");
  string beamB = settings->word("Beams:beamB");
  int massA= atoi(beamA.c_str());
  int massB= atoi(beamB.c_str());
  if(massA==1) {
    int id = InitialCondition::findAZ(beamA);
    ParticleDataEntry* pd=particleData->findParticle(id);
    mA = pd->m0();
  }
  if(massB==1) {
    int id = InitialCondition::findAZ(beamB);
    ParticleDataEntry* pd=particleData->findParticle(id);
    mB = pd->m0();
  }

  if(eLab > 0.0) {
    pLab=sqrt(eLab*(2*mA+eLab));
    eCM = sqrt(pow2(eLab+mA+mB)-pLab*pLab);
    pythia->settings.parm("Beams:eCM",eCM);
    pythia->settings.parm("Beams:pLab",pLab);
  } else if(pLab > 0.0) {
    eLab=sqrt(mA*mA+pLab*pLab)-mA;
    eCM = sqrt(pow2(eLab+mA+mB)-pLab*pLab);
    pythia->settings.parm("Beams:eCM",eCM);
    pythia->settings.parm("Beams:eLab",eLab);
  } else {
    eLab=(eCM*eCM -mA*mA - mB*mB)/(2*mB)-mA;
    pLab=sqrt(eLab*(2*mA+eLab));
    pythia->settings.parm("Beams:pLab",pLab);
    pythia->settings.parm("Beams:eLab",eLab);
  }
  double sCM=eCM*eCM;
  double pCM=sqrt((sCM-(mA+mB)*(mA+mB))*(sCM-(mA-mB)*(mA-mB))/(4*sCM));
  double eA=sqrt(mA*mA+pCM*pCM);
  double yCM=0.5*log((eA+pCM)/(eA-pCM));
  double gamCM = eA/mA;

  cout << " ecm= "<<eCM << " pcm= "<< pCM << " ycm= "<< yCM
        << " pLab= "<< pLab << " eLab= "<< eLab
       << " gamma= "<< gamCM
       << " v= "<< pCM/eA
       <<endl;
  if(pythia->settings.flag("MeanField:gammaCorrectionGaussian")) {
    pythia->settings.parm("Cascade:gamCM",gamCM);
  }

  double pz=PCM(eCM,mnuc,mnuc);
  double e1=sqrt(mnuc*mnuc+pz*pz);
  pythia->settings.parm("Beams:eA", e1);
  pythia->settings.parm("Beams:eB", e1);
  pythia->settings.mode("Beams:idA",2212);
  pythia->settings.mode("Beams:idB",2212);


  //jpythia = new JPythia(xmlDir,false);
  jpythia = new JPythia(pythia->settings,pythia->particleData,false);
  jpythia2 = new JPythia(defaultSettings,pythia->particleData,false);
  jpythia->settings->addFlag("Cascade:constQuarkDiffractive",false);
  jpythia2->settings->addFlag("Cascade:constQuarkDiffractive",false);
  jpythia->settings->flag("Cascade:constQuarkDiffractive",settings->flag("Cascade:constQuarkDiffractive"));
  jpythia2->settings->flag("Cascade:constQuarkDiffractive",settings->flag("Cascade:constQuarkDiffractive"));

  jpythia->settings->flag("SoftQCD:inelastic", true);
  jpythia2->settings->flag("SoftQCD:inelastic", true);

  jpythia->settings->parm("Beams:eA", e1);
  jpythia->settings->parm("Beams:eB", e1);
  jpythia->settings->mode("Beams:idA",2212);
  jpythia->settings->mode("Beams:idB",2212);
  jpythia->readString("HadronLevel:Decay = off");

  jpythia2->settings->parm("Beams:eA", e1);
  jpythia2->settings->parm("Beams:eB", e1);
  jpythia2->settings->mode("Beams:idA",2212);
  jpythia2->settings->mode("Beams:idB",2212);
  jpythia2->readString("HadronLevel:Decay = off");

  hadronize->init();
  pythia->readString("ProcessLevel:all = off");
  pythia->readString("HadronLevel:Decay = off");
  pythia->init();
  pydecay->init();

  iEvent = 0;
  nEvent = settings->mode("Main:numberOfEvents");
  collisionUpdateFreq=settings->mode("Cascade:collisionUpdateFreq");
  overSample = settings->mode("Cascade:overSample");

  int const initCond = settings->mode("Cascade:initialCondition");
  if(initCond==1) {
    initcnd = new BoostedTwoNuclei(settings,jamParticleData,rndm);
    initcnd->setNumberOfOverSample(overSample);
    cout << "initial condition:colliding two nuclei initCond= "<<initCond <<endl;
    initcnd->init();
    iniTwo = dynamic_cast<BoostedTwoNuclei*>(initcnd);
  } else if(initCond==2) {
    initcnd = new Angantyr(&defaultSettings,jamParticleData,rndm);
    cout << "Angantyr initial condition initCond= "<< initCond <<endl;
    initcnd->init();
  } else if(initCond==3) {
    initcnd = myinitcond;
    initcnd->init();
  } else {
    cout << "JAM::wrong value of initCond= "<< initCond <<endl;
    exit(1);
  }

  //pyDecay = new ParticleDecays();
  //TimeShower*    timesDec    = 0;
  //TimeShower*    timesDec    = new SimpleTimeShower();
  //DecayHandler*  decayHandle = 0;
  //vector<int> handledParticles;

  StringFlav* flavSel = new StringFlav();
  flavSel->init(*settings,  particleData, rndm, info);

  //pyDecay->init(info, *settings,particleData, rndm,
  //        &pythia->couplings, timesDec, flavSel, decayHandle,handledParticles);

  decay = new HadronDecay();
  decay->init(info,*settings,jamParticleData,pydecay,rndm);

  xsection = new CrossSection(info,settings,jamParticleData,flavSel,rndm);
  scatt = new Scatter(info,settings,jamParticleData,xsection,hadronize,flavSel,rndm);

  //double eCMPythia = settings->parm("Cascade:minCMEnergyPythia");
  double eCMPythia = settings->parm("Beams:eMinPert");
  cascadeMethod = settings->mode("Cascade:model");

  if(cascadeMethod>0) {
  if(eCM*1.3 > eCMPythia) {

    jpythia->init0(jamParticleData,decay);
    scatt->setPythia(jpythia);
  }
  if(eCM > 500.0) {
    jpythia2->init0(jamParticleData,decay);
    scatt->setPythia2(jpythia2);
  }
  }

  doPauliBlock=settings->mode("Cascade:PauliBlocking");
  cout << "# pauli= " << doPauliBlock <<endl;

  doNuclearClusterFormation = settings->flag("Cascade:nuclearClusterFormation");
  double rc0 = settings->parm("Cascade:nuclearClusterR");
  double pc0 = settings->parm("Cascade:nuclearClusterP");
  if(doNuclearClusterFormation) {
    nuclearCluster = new NuclearCluster(rc0,pc0);
  }

  // Hydro
  withHydro = settings->flag("Hydro:mode");
  hydroInitCond=0;
  if(withHydro) {
    fluidHandler = new FluidHandler(settings,jamParticleData,decay,rndm);
    hydroInitCond = settings->mode("Hydro:InitialCondition");
  }

  withMeanField = settings->mode("MeanField:mode");
  optVectorPotential=settings->mode("MeanField:optVectorPotential");
  optVdot=settings->mode("MeanField:optVdot");

  // initialize RQMD part.
  optPropagate=0;
  if(withMeanField == 1) {
    meanField = new QMD(settings);
    cout << "Meanfield mode:non-relativistic QMD " << withMeanField <<endl;

  } else if(withMeanField == 2) {
    meanField = new RQMDs(settings);
    cout << "Meanfield mode:RQMDs RQMD Skyrme scalar " << withMeanField <<endl;
    cout << "obsolete do not use this"<<endl;
    exit(1);

  } else if(withMeanField == 3) {
    meanField = new RQMDv(settings);
    cout << "Meanfield mode:RQMDv RQMD Skyrme vector " << withMeanField <<endl;
    if(optVectorPotential==1 && optVdot==0) optPropagate=1;
    cout << "obsolete do not use this"<<endl;
    exit(1);

  } else if(withMeanField == 4) {
    meanField = new RQMDsv(settings);
    cout << "Meanfield mode:RQMDsv RQMD scalar-vector " << withMeanField <<endl;
    if(optVectorPotential==1 && optVdot==0) optPropagate=1;

  } else if(withMeanField == 5) {
    meanField = new RQMDw(settings);
    cout << "Meanfield mode:RQMDw RQMD.RMF " << withMeanField <<endl;
    if(optVectorPotential==1 && optVdot==0) optPropagate=1;
    cout << "obsolete do not use this"<<endl;
    exit(1);

  } else if(withMeanField == 11) {
    meanField = new QMDBox(settings,initcnd);
    cout << "Meanfield mode:QMDBox non-relativistic QMD Skyrme " << withMeanField <<endl;

  } else if(withMeanField == 12) {
    meanField = new RQMDsBox(settings,initcnd);
    cout << "Meanfield mode:RQMDsBox RQMD Skyrme scalar " << withMeanField <<endl;
    cout << "obsolete do not use this"<<endl;
    exit(1);

  } else if(withMeanField == 13) {
    meanField = new RQMDvBox(settings,initcnd);
    cout << "Meanfield mode:RQMDvBox RQMD Skyrme vector " << withMeanField <<endl;
    if(optVectorPotential==1 && optVdot==0) optPropagate=1;
    cout << "obsolete do not use this"<<endl;
    exit(1);

  } else if(withMeanField == 14) {
    meanField = new RQMDsvBox(settings,initcnd);
    cout << "Meanfield mode:RQMDsvBox RQMD scalar-vector " << withMeanField <<endl;
    if(optVectorPotential==1 && optVdot==0) optPropagate=1;

  } else if(withMeanField == 21) {
    meanField = new GBUU(settings);
    cout << "Meanfield mode:BUU " << withMeanField <<endl;
  } else if(withMeanField == 22) {
    meanField = new RQMDs0(settings);
    cout << "Meanfield mode:RQMDs0 RQMD Skyrme scalar old version " << withMeanField <<endl;
  } else if(withMeanField==0) {
    settings->mode("ParticleDecays:potentialHandling",0);
    cout << "Cascade mode: no mean-field " << withMeanField <<endl;
  } else {
    cout << "JAM wrong Menfield:mode " << withMeanField <<endl;
    exit(1);
  }

  xColl=0;
  xDec=0;
  xCollBB=0;
  xCollMB=0;
  xCollMM=0;
  xCollBBar=0;
  xCollBarBar=0;
  xElastic=0;
  xAbsorb=0;
  xInter=0;
  xCollRR2NN=0;

  cout << "# trune = "<< jpythia->settings->mode("Tune:pp") <<endl;
  cout << "#  inel = " << jpythia->settings->flag("SoftQCD:inelastic") <<endl;
  cout << " soft nondiff = " << jpythia->settings->flag("SoftQCD:nonDiffractive") <<endl;

  cout << "# sigmaTotal:mode= " << jpythia->settings->mode("SigmaTotal:mode") <<endl;
  cout << "# StringFlav:probQQtoQ " << jpythia->settings->parm("StringFlav:probQQtoQ") <<endl;
  cout << "# StringFlav:popcornRate= " << jpythia->settings->parm("StringFlav:popcornRate") <<endl;

  initTimeDependentAna(gamCM,yCM);

  if(cascadeMethod==1) {
    event = new Collision1(settings,jamParticleData,xsection,rndm);
  } else if(cascadeMethod == 3) {
    event = new Collision3(settings,jamParticleData,xsection,rndm);
  } else {
    cout << "wrong value cascade Cascade:model= "<< cascadeMethod<<endl;
    exit(1);
  }
  event->init(initcnd);

  //vector<ResonanceWidths*> resonancePtr;
  //particleData->initWidths(resonancePtr);

  // Check resonance mass.
  //int idNow=313;
  //idNow=12112;
  //particleData->mMin(idNow,0.4);
  //particleData->mWidth(idNow,0.4);
  //ParticleDataEntry* pd=particleData->findParticle(idNow);
  //cout << " useBreitWigner= " << pd->useBreitWigner() <<endl;
  //for(int i=0;i<100;i++) {
  //  double mNow = particleData->mSel(idNow);
  //  cout << " mNow= " << mNow<<endl;
  //}

  return true;
}

bool JAM::next()
{
  iEvent++;

  //if(event) delete event;
  event->clear();

  nPrint=0;
  nColl=nDec=ncollBB=ncollMB=ncollMM=ncollBBar=nElastic=nAbsorb=nPauli=0;
  ncollRR2NN=0, ncollBarBar=0;
  eventInitTimeDependentAna();

  double dt0 = settings->parm("Cascade:TimeStepSize");
  nStep = settings->mode("Cascade:TimeStep");

  // create initial condition and count Glauber type collision number.
  initcnd->generate(event,1);
  numInter = initcnd->nColl();// initially predicted collision number.
  numPart  = initcnd->nPart();
  //if(cascadeMethod==0) nColl=numInter;
  impactPar = initcnd->getImpactPar();

  //double tFirstCol=0.0;
  //if(iniTwo) tFirstCol = event->collisionOrderTime();

  int hyswitch=0;
  if(withHydro) hyswitch = fluidHandler->init(event,iniTwo);
  if(withMeanField) meanField->init(event->plist);

  double v= settings->parm("MeanField:stepVelocity");
  double stime = settings->parm("MeanField:dtExpandStartTime");
  double timeNow = pythia->settings.parm("Cascade:timeStart");

  // Time evolution loop.
  for(int step=1; step <= nStep; step++) {

    double dt = dt0 + v* std::max(0.0,timeNow-stime);

    //timeNow = dt*(step-1);
    double nextTime = timeNow + dt;
    //if(nStep>1) tFirstCol = timeNow;

    anaTimeDependent(step,timeNow);

    if(step==hyswitch) fluidHandler->convertAll(event,step);

    int noCollUpdate= step % collisionUpdateFreq;

    // collision list needs not be updated,
    if(cascadeMethod==0) noCollUpdate=1;


    // Perform fluid evolution, and check freeze-out hyper surface.
    // If necessary, fluid elements are converted into particles.
    // For no hadronic cascade option, stop time evolution after convert all fluid to particle.
    if(withHydro) {
      if(!fluidHandler->evolution(event,step,timeNow,noCollUpdate)) break;
       //if(isDebug) {cout << " after hydro evol"<<endl; computeTotalEnergy(step,dt); }
    }

    if(withMeanField) {
      meanField->evolution(event->plist,timeNow,dt,step);
    }

    // make new collision list among all particles.
    //if(noCollUpdate==0 && tFirstCol < nextTime) {
    if(noCollUpdate==0) {
      //double nextCollUpdateTime = dt*step*collisionUpdateFreq;
      double nextCollUpdateTime = timeNow + dt*collisionUpdateFreq;
      event->makeCollisionList(timeNow,nextCollUpdateTime);
    }

    if(cascadeMethod) cascade(timeNow,nextTime);
    event->propagate(nextTime,optPropagate,nStep);

    // Compute total energy and momentum.
    if(isDebug) computeTotalEnergy(step,dt);

    timeNow = nextTime;

  } // end time evolution.

  if(withHydro) fluidHandler->finalize(event,nStep+1,timeNow);

  if(settings->flag("Cascade:finalDecay")) finalDecay();
  //finalDecay1();
  //deleteDecayedParticle();

  collisionStatistics2(timeNow);

  printTimeDependentAna();
  
  if(doNuclearClusterFormation) nuclearCluster->findCluster(event->plist);

  return true;
}

void JAM::computeTotalEnergy(int step,double dt)
{
 if(withMeanField)
   pTot = meanField->computeEnergy(event->plist,step);
 else {
   pTot=0.0;
   for(auto& i: event->plist) {
        if(optPropagate==1) {
	  pTot[1] += i->getPx();
	  pTot[2] += i->getPy();
	  pTot[3] += i->getPz();
	  Vec4 pk= i->getPkin();
	  pTot[0] += pk[0];
        } else {
	  pTot += i->getP();
        }
      pTot[0] += i->potv(0);
    }
  }

  pTot /= overSample;
  if(step==1) {
    //numInitialPart = particleList.size();
    //pTot0 = pTot/numInitialPart;
    pTot0 = pTot;
  }

  if(withHydro) pTot += fluidHandler->pFluidTotal();

//  pTot /= numInitialPart;
  double econ=abs(pTot0[0]-pTot[0])/pTot0[0]*100;
  cout << "time= " << fixed << dt*(step-1)
        << " econ= " << fixed << econ << " %"
        //<< " econ= " << scientific << setprecision(8) << econ << " %"
       << scientific << setw(15) << pTot[1]
       << scientific << setw(15) << pTot[2]
       << scientific << setw(15) << pTot[3]
       <<endl;

  if(abs(pTot[1])>1e-3) {
    cout << " momentum does not conserve " << pTot[1]<<endl;
    exit(1);
  }

  if(nStep>1 && optPrintDisplay > 0) display(dt*step);
}

void JAM::collisionStatistics2(double ftime)
{
  //if(analysis) analysis->eventAnalysis(*event);

  //nColl = event->getNColl();
  xColl   += nColl;
  xCollBB += ncollBB;
  xCollMB += ncollMB;
  xCollMM += ncollMM;
  xCollBBar += ncollBBar;
  xCollBarBar += ncollBarBar;
  xDec    += nDec;
  xInter += numInter;
  xAbsorb += nAbsorb;
  xElastic += nElastic;
  xPauli += nPauli;
  xCollRR2NN += ncollRR2NN;

  if(isDebug) {
    cout << " total particle before final decay= "
	<< event->plist.size()
	<<endl;
  }

  if(isDebug) {
    cout << "# event= " << iEvent << " nColl= " << nColl
         << " nDecay= " << event->getNDecay()
         << " total particle= " << event->plist.size()
         << " finale time = " <<  ftime
         << endl;
    cout << " ncol= " << nColl << " ndec= " << nDec <<endl;
    cout << " BB= " << ncollBB << " MB= " << ncollMB
         << " MM= " << ncollMM
         << " BBar= " << ncollBBar
         << endl;
  }

}

// do cascade process; two-body collision and decay.
double JAM::cascade(double initime,double finaltime)
{
  //double difftime=0.0;
  InterList* inter = event->findNextCollision();
  vector<EventParticle*> outgoing;
  int operation = 0;
  double coltime=0.0;
  //int icoll1=0, icoll2=0;
  int idec=0, icoll=0;
  while(inter) {
    operation++;
    //double coltime0=coltime;
    // collision time.
    coltime=inter->getCollisionOrderTime();
    //if(coltime > finaltime) { coltime = coltime0; break; }

    if(printColl) {
      if(inter->getNumberOfInComing()>0) {
      inter->print(cout);
      cout << "incoming = "<< inter->getNumberOfInComing() <<endl;
      }
    }

    int iconv=0;
     // Collision.
    if(inter->getNumberOfInComing() > 1) {

      scatt->scatter(inter,outgoing,event);

      // collision was cancelled.
      if(outgoing.size() == 0) {
        event->cancelCollision(inter);
        inter = event->findNextCollision();
	continue;
      }
      icoll++;
      if(hydroInitCond==3)
        fluidHandler->checkFluidConversion(inter,outgoing,finaltime);

    // Wall collision.
    } else if(inter->getNumberOfInComing()==0) {

      double tw=inter->getCollisionOrderTime();

      if(isDebug>2) cout << "Before wall collision tw= "<< tw << " p= "<< inter->getParticle(0)
	<< " x= " << inter->getParticle(0)->getR()
        << " id= "<< inter->getParticle(0)->getID()
	  <<endl;

      event->wallCollision(*inter);
      inter = event->findNextCollision();
      continue;

    // Decay
    } else {

      EventParticle *ep= inter->getParticle(0);
      if(ep->getStatus()==-1200) {
        fluidHandler->convertAform(event,ep,outgoing,coltime);
      } else {
        idec++;
        decay->decay(ep,outgoing);

	// Check if decayed particle can be converted into fluid.
	if(hydroInitCond==3)
	  iconv = fluidHandler->checkFluidConversion(ep,outgoing,finaltime);
      }

    }

    checkEnergyMomentumConservation(inter,outgoing);

    bool block = false;
    if(doPauliBlock && outgoing.size()>0 && iconv<=0) {
      block = event->doPauliBlocking(inter,outgoing,optPropagate);
    }

    if(!block) {
      for(int i=0; i<(int)outgoing.size();i++) {
	if(outgoing[i]->getStatus()== -1000) {
	  fluidHandler->convert(event,outgoing[i],finaltime);
	} else {
	  event->setPnewList(outgoing[i]);
	}
      }
      for(auto& p : outgoing) {
	if(p->getStatus()== -1000) delete p;
      }

      collisionStatistics(coltime,inter,outgoing);
      if(printColl) printCollision(coltime,operation,inter,outgoing);
      event->collisionUpdate(inter);
      //if(analysis) analysis->collisionAnalysis(*inter,*event);

    // This event is Pauli-blocked.
    } else {

      for(auto& p : outgoing) delete p;

      event->cancelCollision(inter);
      nPauli++;
      if(printColl) cout << " collision was Pauli blocked"<<endl;
    }

    outgoing.clear();
    inter = event->findNextCollision();


  }  // end inter

  nDec += idec;
  if(isDebug>1) {
    cout << fixed << "last coltime= "<< coltime << " nexttime= " << finaltime
        << " operation= " << operation
        << " collision= " << icoll
        << " decay= " << idec
        << endl;
  }

  return coltime;
}

void JAM::printCollision(double coltime, int operation,InterList* inter, vector<EventParticle*> outgoing)
{
  int np = inter->getNumberOfInComing();
  double ecm=0.0,sigma=0.0, sigel=0.0, sigab=0.0;
  EventParticle* p1=inter->getParticle(0);
  EventParticle* p2=0;
  ParticleDataEntry *pd1= p1->getParticleDataEntry();
  ParticleDataEntry *pd2=0;
  int id2=0;
  int id1=p1->getPID();
  if(np==2) {
    TwoBodyInterList *inter2=dynamic_cast<TwoBodyInterList*>(inter);
    CollisionPair cpair = inter2->getCpair();
    sigma=cpair.getSigma();
    sigel=cpair.getSigmaElastic();
    sigab=cpair.getSigAbs();
    ecm=cpair.getCMenergy();
    p2=inter->getParticle(1);
    pd2= p2->getParticleDataEntry();
    id2=p2->getPID();
  } else {
    ecm= p1->getMass();
  }

  if(outgoing.size()==2) {
  int idp[2]={0};
    for(int i=0;i<(int)outgoing.size();i++) {
	idp[i]=outgoing[i]->getPID();
    }
    int isdelta=0;
    //if(id1==id_delt || id2==id_delt) isdelta=1;
    if(id1 !=id_nucl || id2 !=id_nucl) isdelta=1;
    int isnn=0;
    if(idp[0]==id_nucl && idp[1]==id_nucl) isnn=1;
    if(isdelta*isnn==1) {
	cout << "<<delta absorption>> "<< pd1->name() << " + " << pd2->name()
	  << " -> " << outgoing[0]->getParticleDataEntry()->name()
	  << " + " << outgoing[1]->getParticleDataEntry()->name()
	  <<endl;
    }
  }

  cout << "coltime= "<< coltime << " operation= " << operation;
  if(np==2) cout << " ecm= "<< ecm <<  " sigma= "<< sigma << " sigel= "<< sigel<< " sigab= "<< sigab;
  cout << endl;


  cout << "p1= "<<  p1->getID()
       << " status= "<< p1->getStatus()
       << " m= " << p1->getMass()
	 << " " << pd1->name()
	 << " t= " << p1->getT()
	 << " tf= " << p1->getTf()
	 << " tcol= " << inter->getCollisionTime(0)
	 <<endl;

  // collision
  if(np == 2) {
    cout << "p2= " << p2->getID()
       << " status= "<< p2->getStatus()
	 << " m= " << p2->getMass()
	 << " " << pd2->name()
	 << " t= " << p2->getT()
	 << " tf= " << p2->getTf()
	 << " tcol= " << inter->getCollisionTime(1)
         <<endl;
  }

    for(int i=0;i<(int)outgoing.size();i++) {
      cout << "p"<< i+3 << "= "<< outgoing[i]->getID()
          << " status= "<< outgoing[i]->getStatus()
	  << " " << outgoing[i]->getParticleDataEntry()->name()
	  << " m= "<< outgoing[i]->getMass()
	  << " pz= "<< outgoing[i]->getPz()
	  << " cq1= "<< outgoing[i]->constQuark(0)
	  << " cq2= "<< outgoing[i]->constQuark(1)
	  << " t= "<< outgoing[i]->getT()
	  << " tf= "<< outgoing[i]->getTf()
          <<endl;
      double t=outgoing[i]->getT();
      double z=outgoing[i]->getZ();
      if(t < abs(z)) {
	cout << " t < z t= "<< t << " z= "<< z<<endl;
	exit(1);
      }
    }



  //event->PrintCollision(cout);
  if(printColl > 3) {
    //for(auto i:event->plist)  i->print(cout);
    cout << "enter any key to proceed" << endl;
    cin.get();
   }

}

void JAM::collisionStatistics(double coltime,InterList* inter,vector<EventParticle*> outgoing)
{
  statTimeDependentAna(coltime,inter,outgoing);

  int np = inter->getNumberOfInComing();
  if(np==1) return; // decay

  EventParticle* p1=inter->getParticle(0);
  EventParticle* p2=inter->getParticle(1);
  //ParticleDataEntry *pd1= p1->getParticleDataEntry();
  //ParticleDataEntry *pd2=p2->getParticleDataEntry();
  int pid1=p1->getPID();
  int pid2=p2->getPID();
  TwoBodyInterList *inter2=dynamic_cast<TwoBodyInterList*>(inter);
  CollisionPair cpair = inter2->getCpair();
  double b=inter2->getImpactPar();
  bMax=max(bMax,b);
  //double sigma=cpair.getSigma();
  //double sigel=cpair.getSigmaElastic();
  //double sigab=cpair.getSigAbs();
  //double ecm=cpair.getCMenergy();

  int pid3=0, pid4=0;
  if(outgoing.size()==2) {
    pid3=outgoing[0]->getPID();
    pid4=outgoing[1]->getPID();
  int isr=0, isnn=0;
  if(pid1 !=id_nucl || pid2 !=id_nucl) isr=1;
  if(pid3==id_nucl && pid4==id_nucl) isnn=1;
  if(isr*isnn==1) ncollRR2NN++;
  }

  int cltype= inter->getCollType();
  nColl++;
  switch(cltype) {
  case 1: ncollBB++;break;
  case 2: ncollMB++;break;
  case 3: ncollMM++;break;
  case 4: ncollBBar++;break;
  case 5: ncollBarBar++;break;
  default:
    cout << "JAM:wrong cltype " << cltype<<endl;
    exit(1);
  }

  int channel=scatt->getChannel();
  if(channel==Scatter::ELASTIC) nElastic++;
  else if(channel==Scatter::ABSORPTION) nAbsorb++;


  bool testprint=false;
  if(testprint) {
  double ecm=0.0,sigma=0.0, sigel=0.0, sigab=0.0;
  ParticleDataEntry *pd1= p1->getParticleDataEntry();
  ParticleDataEntry *pd2= p2->getParticleDataEntry();
  int id1=p1->getID();
  int id2=p2->getID();
    sigma=cpair.getSigma();
    sigel=cpair.getSigmaElastic();
    sigab=cpair.getSigAbs();
    ecm=cpair.getCMenergy();
    p2=inter->getParticle(1);
    pd2= p2->getParticleDataEntry();
    id2=p2->getID();
  cout << " ecm= "<< ecm <<  " sigma= "<< sigma << " sigel= "<< sigel<< " sigab= "<< sigab
      <<endl;
  cout << "incomoing: id1= "<< id1 << " " << pd1->name() << " id2= "<< id2;
  if(id2 !=0) cout << " " << pd2->name();
  cout <<endl;
  cout << "outgoing:" << outgoing.size() << " ";
  for(auto p:outgoing) {
   cout << p->getID() << " " << p->getParticleDataEntry()->name() << " ";
   }
   cout <<endl;

  }

}

void JAM::checkEnergyMomentumConservation(InterList* inter,
	vector<EventParticle*>& outgoing)
{
  if(outgoing.size()==0) return;

  Vec4 ptot=0.0;
  for(int i=0; i<(int)outgoing.size();i++) {
    ptot += outgoing[i]->getP();

    /*
    if(outgoing[i]->getMass() > 4.0) {
	cout << " hadron mass too large ? m= " <<
             outgoing[i]->getMass()
	     << " cltype= "<< inter->getCollType()
	     << " chanel= "<< scatt->getChannel()
	     << " id= "<<  outgoing[i]->getID()
	     <<endl;
    }
    */
  }

  Vec4 ptot0=inter->getTotalMom();
  double diff= (ptot0-ptot).pAbs();
  if(diff > 1e-3) {
	cout << "JAM::checkEnergyMom after operation momentum does not conserve"<<endl;
	cout << " diff= " << diff <<endl;
	cout << " ptot0=" << ptot0<<endl;
	cout << " ptot =" << ptot <<endl;
        exit(1);
  inter->print(cout);
  for(int i=0; i<(int)outgoing.size();i++) {
    outgoing[i]->print(cout);
  }

  }

}

// Force decay of particles, but we do not eliminate decayed particles
// from the list.
void JAM::finalDecay1()
{
  vector<EventParticle*> outgoing;
  std::list<EventParticle*> &particleList = event->plist;

  int ntry=0;
  list<EventParticle*>::iterator iDec = particleList.begin();
  do {
    ParticleDataEntry* pd=(*iDec)->getParticleDataEntry();
    if ( (*iDec)->getStatus() > 0 && pd->canDecay() && pd->mayDecay() ) {
      decay->decay(*iDec,outgoing,true);
      for(int i=0; i<(int)outgoing.size();i++)
        particleList.push_back(outgoing[i]);
      outgoing.clear();
      (*iDec)->setStatus(-1);
    }
    if(++ntry > 30000) {
	cout << " too many final decays? " << particleList.size()<<endl;
	break;
    }
  } while (++iDec != particleList.end());

}

// Force decay of particles.
// All decayed particles are eliminated from the list.
void JAM::finalDecay()
{
  vector<EventParticle*> outgoing;
  std::list<EventParticle*> &particleList = event->plist;

  list<EventParticle*>::iterator first = particleList.begin();
  while(first != particleList.end()) {
    list<EventParticle*>::iterator next =  first;
    ++next;
    ParticleDataEntry* pd=(*first)->getParticleDataEntry();
    if ( pd->canDecay() && pd->mayDecay() ) {

	if(printColl) {
	cout << "final decay pd= " << pd->name() << " id= " << pd->id()
	    << " m= " << (*first)->getMass()
            << " lifetime= " << (*first)->lifetime()
            <<endl;
	}
	if(pd->id() == 22) {
	    cout << " photon final decay? id= "<< pd->id()<<endl;
	    continue;
	}

      decay->decay(*first,outgoing,true);

      for(int i=0; i<(int)outgoing.size();i++) {
        particleList.push_back(outgoing[i]);
      }
      outgoing.clear();
      delete *first;
      particleList.erase(first);
    }
    first = next;
  }

}

void JAM::collisionInfor(InterList* inter)
{
    EventParticle* ip[2];
    ip[0] = inter->getParticle(0);
    ip[1] = inter->getParticle(1);
    int kf1 =  ip[0]->getID();
    int kf2 =  ip[1]->getID();
    if(abs(kf1)<10 || abs(kf2) <10) {
	cout << "quark scatt";
	cout << " kf1 = " << kf1 << " kf2=" << kf2 << endl;
    }

}

void JAM::printDisplay(double ctime)
{
     double dprint=0.5;
     int iprint=ctime/dprint;
     //cout << " ctime= " << ctime << " iprint= " << iprint << " nPrint= "<< nPrint <<endl;
     if(iprint >= nPrint) {
       for(int i=nPrint; i<=iprint;i++) display(i*dprint);
       nPrint = iprint+1;
    }

    //if((int)(ctime/1.0) % 10 ==0 ) display(ctime);

}

void JAM::display(double tnow)
{
    //ofstream ofs("aaa.dat");

    double scale=1.0, dr=pythia->settings.parm("Cascade:displayScale");
    const int maxz=30, maxx=20;
    //double xmin=-30.0;
    double xmin=-15.0;
    double ymin=-10.0;
    int iaxis=2;
    int ix=iaxis%3 + 1;
    int iy=ix%3 + 1;

    int ncount[maxz][maxx];
    for(int i=0; i<maxz;i++)
    for(int j=0; j<maxx;j++)
	ncount[i][j]=0;


    list<EventParticle*>::const_iterator jp;
    list<EventParticle*>  plist = event->plist;
    int mcount=0;
    for(jp=plist.begin(); jp != plist.end(); jp++) {
	Vec4 p = (*jp)->getP();
	Vec4 r = (*jp)->getR();
	double dt = tnow - r[0];

	if(dt < 0) continue;
	double x = r[ix] + p[ix]/p[4]*dt;
	double y = r[iy] + p[iy]/p[4]*dt;

	//ofs << setw(13) << x << setw(13) << y<<endl;

//	cout << " z= " << x << " x= " << y <<endl;
//	cin.get();

	int i=(x-xmin)*scale/dr;
	int j=(y-ymin)*scale/dr;
	if((i>=0 && i<maxz) && j>=0 && j<maxx) {
	    mcount++;
	    ncount[i][j]++;
	}
    }

    cout << " time = " << tnow << " mcount= "<< mcount << endl;
    for(int j=0;j<maxx; j++) {
    for(int i=0;i<maxz; i++) {
	cout << " ";
	string s=" ";
	stringstream ss;
	if(ncount[i][j]>0) {
	    ss << ncount[i][j]; s = ss.str();
	}
	cout << s;
    }
    cout << endl;
    }

    //ofs.close();
    //cout << " mcount= "<< mcount<<endl;
    //cin.get();

}

void JAM::initTimeDependentAna(const double gamCM,const double yCM)
{
  printFreq = settings->mode("Analysis:printFreq");
  printFreqPhaseSpace = settings->mode("Analysis:printFreqPhaseSpace");
  bool anaTPart = pythia->settings.flag("Analysis:timeDependenceParticle");
  bool anaTFlow = pythia->settings.flag("Analysis:timeDependenceFlow");
  bool anaTDens = pythia->settings.flag("Analysis:timeDependenceDensity");
  bool anapot = pythia->settings.flag("Analysis:Potentials");
  bool anacoll = pythia->settings.flag("Analysis:collision");
  bool anaph = pythia->settings.flag("Analysis:outPutPhaseSpace");

  if(anaTPart || anaTFlow || anaTDens || anaPot|| anaColl || anaph) {
    double dt0 = settings->parm("Cascade:TimeStepSize");
    int nstep = settings->mode("Cascade:TimeStep");
    double v= settings->parm("MeanField:stepVelocity");
    double stime = settings->parm("MeanField:dtExpandStartTime");
    double ftime=0.0;
    for(int step=1; step <= nstep; step++) {
      ftime += dt0 + v*std::max(0.0,ftime-stime);
    }
    double ycut = pythia->settings.parm("Analysis:yCut");
    double ycutf = pythia->settings.parm("Analysis:yCutFoward");
    double ycutmax = pythia->settings.parm("Analysis:yCutMax");
    if(nstep==1) dt0=ftime/50;
    if(anaTPart) anaParticle = new AnaTimeDepParticle(ftime,dt0*printFreq,gamCM,ycut);
    if(anaTFlow) anaFlow = new AnaTimeDepFlow(ftime,dt0*printFreq,gamCM,ycut,ycutf,ycutmax);
    if(anaTDens) anaDens = new AnaTimeDepDensity(ftime,dt0*printFreq,gamCM,yCM*0.9);
    if(anapot) anaPot = new AnaMeanField(ftime,dt0*printFreq,gamCM,ycut);
    if(anacoll) anaColl = new CollisionHistory(eCM,ftime,dt0*printFreq,gamCM);
    if(anaph) anaPhase = new AnaOutPutPhaseSpace(ftime,dt0*printFreqPhaseSpace);
  }
}

void JAM::eventInitTimeDependentAna()
{
  if(anaParticle) anaParticle->init();
  if(anaFlow) anaFlow->init();
  if(anaDens) anaDens->init();
  if(anaPot) anaPot->init();
  if(anaPhase) anaPhase->init();
}

void JAM::anaTimeDependent(const int step, const double timeNow)
{
    //if(anaFlow && nStep>1 && step%printFreq==0) anaFlow->fill(timeNow,event->plist);
    //if(anaDens && nStep>1 && step%printFreq==0) anaDens->fill(timeNow,event->plist);
    if(anaParticle && nStep>1 && (step-1)%printFreq==0) anaParticle->ana(-1,timeNow,event->plist);
    if(anaFlow && nStep>1 && (step-1)%printFreq==0)     anaFlow->ana(-1,timeNow,event->plist);
    if(anaDens && nStep>1 && (step-1)%printFreq==0)     anaDens->ana(-1,timeNow,event->plist);
    if(anaPot && nStep>1 && (step-1)%printFreq==0)      anaPot->fill(timeNow,event->plist);
    if(anaPhase && nStep>1 && (step-1)%printFreqPhaseSpace==0) anaPhase->ana(iEvent,-1,timeNow,event->plist);
}

void JAM::statTimeDependentAna(const double coltime,InterList* inter,vector<EventParticle*> outgoing)
{
  if(anaColl) anaColl->fill(inter,outgoing);
  if(nStep==1 && anaParticle) anaParticle->fill(coltime,event->plist);
  if(nStep==1 && anaFlow)     anaFlow->fill(coltime,event->plist);
  if(nStep==1 && anaDens)     anaDens->fill(coltime,event->plist);
  if(nStep==1 && optPrintDisplay > 0) printDisplay(coltime);
}

void JAM::printTimeDependentAna()
{
  if(anaPot && iEvent%10==0) anaPot->print("JAMTimeDependentPotential.dat",iEvent);
  if(anaParticle && iEvent%10==0) anaParticle->print("JAMTimeDependentParticle",iEvent);
  if(anaFlow && iEvent%10==0) anaFlow->print("JAMTimeDependentFlow");
  if(anaDens && iEvent%10==0) anaDens->print("JAMTimeDependentDensity.dat",iEvent);
}

void JAM::finTimeDependentAna()
{
  if(anaColl) {
    anaColl->print(iEvent,"JAMCollision.dat");
    delete anaColl;
  }
  if(anaParticle) {
    anaParticle->print("JAMTimeDependentParticle",iEvent);
    delete anaParticle;
  }
  if(anaFlow) {
    anaFlow->print("JAMTimeDependentFlow");
    delete anaFlow;
  }
  if(anaDens) {
    anaDens->print("JAMTimeDependentDensity.dat",iEvent);
    delete anaDens;
  }

  if(anaPot) {
    anaPot->print("JAMTimeDependentPotential.dat",iEvent);
    delete anaPot;
  }

  if(anaPhase) {
    delete anaPhase;
  }
}

} // end namespace jam2
