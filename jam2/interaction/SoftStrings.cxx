// soft string excitation
#include <Pythia8/PythiaStdlib.h>
#include <jam2/hadrons/JamStdlib.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <jam2/interaction/SoftStrings.h>
#include <jam2/hadrons/JamStdlib.h>

//const double MminBar=0.4;
//const double MminMes=0.2;


const double MminBar=0.0;
const double MminMes=0.1;

const int   NTRYKINMATCH=20;

namespace jam2 {

using namespace std;

using Pythia8::sqrtpos;
using Pythia8::ParticleData;
using Pythia8::Settings;
using Pythia8::Rndm;

SoftStrings::SoftStrings(Pythia8::Info* inf,Pythia8::Settings& set,JamParticleData* jd,
        CrossSection* xs,
        Pythia8::Pythia* py,
        Pythia8::StringFlav* flav,
        Pythia8::Rndm* r) :  Scatter2(inf,&set,jd,xs,r)
{
    jamParticleData=jd;
    particleData=jamParticleData->getParticleData();
    hadronize=py;
    flavSel=flav;

    sigtot.init(&hadronize->info,set,particleData,r);

    hadronContent = new HadronContent(particleData,r);
    beam1 = new BeamHadron(set,particleData,r);
    beam2 = new BeamHadron(set,particleData,r);

    optConstFormationTime=false;
    //optConstFormationTime=true;

    primordialKTremnant = settings->parm("BeamRemnants:primordialKTremnantSoft");

    probDiffra = settings->parm("BeamRemnants:probDiffra");
    probNonDiff= 1.0 - probDiffra;
    probFlavorExchange=settings->parm("BeamRemnants:probFlavorExchange");

    //cout << " kt = " << primordialKTremnant <<endl;
    //cout << " probDiffra= "<< probDiffra<<endl;


  primordialKTsChannelString = 1;
  sigmaAnn  = settings->parm("StringPT:sigmaBBar") / sqrt(2.);

    probBBAnn1=0.3;  // B+Bbar -> string.
    probBBAnn2=0.3;  // B+Bbar -> string + string
    probBBAnn3=0.5;  // B+Bbar -> string + string + meson

    optConstQuarkDiffra=settings->mode("Cascade:optConstQuarkDiffra");
    optConstQscatt=settings->mode("Cascade:optConstQuarkScatt");
    noBaBAnnihilation=settings->mode("Cascade:noBaBAnnihilation");

    aveConstFormationTime=0; aveFormationTime=0;
    nCFormationTime=0; nFormationTime=0;


  // Longitudinal momentum sharing of valence quarks in hadrons.
  xPowMes         = settings->parm("BeamRemnants:valencePowerMeson");

  xPowBarU = settings->parm("BeamRemnants:valencePowerUinP");
  xPowBarD = settings->parm("BeamRemnants:valencePowerDinP");
  xPowBarS         = 0.5 * ( settings->parm("BeamRemnants:valencePowerUinP")
                        + settings->parm("BeamRemnants:valencePowerDinP") );

  valenceDiqEnhance = settings->parm("BeamRemnants:valenceDiqEnhanceSoft");
  eMinPert=13.0;
  eWidthPert=1.0;
  valenceDiqE=6.0;

  // Transverse momentum spread.
  sigmaQBBar      = settings->parm("StringPT:sigmaBBar") / sqrt(2.);
  sigmaQ          = primordialKTremnant/sqrt(2.);

  // Mixing for eta and eta'.
  double theta    = settings->parm("StringFlav:thetaPS");
  double alpha    = (theta + 54.7) * M_PI / 180.;
  fracEtass       = pow2(sin(alpha));
  fracEtaPss      = 1. - fracEtass;

  // time-like final stat shower of strings.
  doShower=settings->flag("Cascade:softRadiation");

  allowRescatterSameString=settings->mode("Cascade:allowRescatterSameString");
  //doShower=false;

  //osOut.open("mass.dat");
  //osOut.open("v.dat");

}

SoftStrings::~SoftStrings()
{
  cout << " ave const. formation time= "<< constFormatinotime()
       << " formation time= " << formatinotime() 
       <<endl;

  delete hadronContent;
  delete beam1;
  delete beam2;
}

// hadron_1 + hadron_2 -> string + string
bool SoftStrings::setKinematics(int typediff,vector<EventParticle*>& outgoing)
{
  typeDiff=typediff;
  double smA=SM[0];
  double smB=SM[1];
  idBeamA=idIn[0];
  idBeamB=idIn[1];
  pCM=pInican;
  eCM=eCMcan;
  MfromCM=MfromCMcan;
  sCM = eCM * eCM;
  int ibar0 = iBaryon[0] + iBaryon[1];

  for(int i=0;i<2;i++) {
    cqBeamA[i] = constQ[0][i];
    cqBeamB[i] = constQ[1][i];
  }

  xDiqEnhance=valenceDiqEnhance
      + valenceDiqE/(1.0 + exp(-(eCM-eMinPert)/eWidthPert));

  //xDiqEnhance=valenceDiqEnhance;

  partonA.clear();
  partonB.clear();
  Vec4 pA, pB;

  /*
  if(optConstQuarkDiffra>=4) {
    //if(preHadronA && preHadronB) typediff=5;
    //else if(preHadronA  &&  !preHadronB)  typediff=4;
    //else if(!preHadronA  &&  preHadronB)  typediff=3;
    //else if(preHadronA  &&  !preHadronB)  typediff=3;
    //else if(!preHadronA  &&  preHadronB)  typediff=4;
    //if(preHadronA && preHadronB) typediff=5;
    if(preHadronA  ||  preHadronB) {
      double xrand = rndm->flat();
      if(xrand < 0.5) typediff=3;
      else typediff=4;
    }
  }
  */

  if(optConstQuarkDiffra>=4) {
    if(optConstQuarkDiffra==5) {
      if(preHadronA && preHadronB) typediff=5;
      else if(preHadronA  &&  !preHadronB)  typediff=4;
      else if(!preHadronA  &&  preHadronB)  typediff=3;
    } else if (optConstQuarkDiffra==6) {
      if(preHadronA && preHadronB) typediff=-1;
      else if(preHadronA  &&  !preHadronB) typediff=4;
      else if(!preHadronA  &&  preHadronB) typediff=3;
    } else if (optConstQuarkDiffra==7) {
      if(preHadronA && preHadronB) typediff=-1;
      else if(preHadronA  &&  !preHadronB)  typediff=3;
      else if(!preHadronA  &&  preHadronB)  typediff=4;
    } else if (optConstQuarkDiffra==8) {
        if(preHadronA  ||  preHadronB) {
        double xrand = rndm->flat();
        if(xrand < 0.33333) typediff=3;
	else if(xrand < 0.66666) typediff=4;
        else typediff=5;
	}
    } else {
        if(preHadronA  ||  preHadronB) {
        double xrand = rndm->flat();
        if(xrand < 0.5) typediff=3;
        else typediff=4;
       }
    }
  }


  if(typediff >= 3) {
    //if(typediff==-1) typediff=1;
    if(!eldiff(typediff,pA, pB)) return false;
  } else if(optConstQuarkDiffra>0 && (preHadronA || preHadronB)) {
     if(!nonDiff()) return false;

    //if(!eldiff(typediff,pA, pB)) return false;
  } else {

    int type=1;
    double xprob= rndm->flat();
    double probnon=1.0;
    double xCut = 5.53;
    if(eCM>4.9 && eCM<xCut) probnon=1.0-
    (-1.7958389716086234 + 0.5840643530390713 *eCM+ -0.05160872446609755 *eCM*eCM+ 0.0015343893713130941 *pow3(eCM));
    else if(eCM >= xCut) {
      calcSigma();
      probnon = probNonDiff;
    }

    if(xprob < probnon) {
    //if(#xprob < probNonDiff) {
      //if(!nondiff()) return false;
      if(!nonDiff()) return false;
    } else {
      type=3;
      double xrand = rndm->flat();
      if(xrand<probSDA) type=4;
      else if(xrand<probSDB) type=5;
      if(!eldiff(type,pA, pB)) return false;
    }
  }

  int idRA=0;
  int idRB=0;
  if(partonA.size()>0) {
      idRA = isResonance(mA, smA,partonA);
      if(idRA == 0) {
	  cout << "Softstrings::setKinemtaics id1= "<< idIn[0] << " idRA= "<< idRA << " mA= "<< mA <<endl;
	  exit(1);
	  return false;
      }
  }
  if(partonB.size()>0) {
      idRB = isResonance(mB, smB,partonB);
      if(idRB == 0) {
	  cout << "Softstrings::setKinemtaics id2="<< idIn[1]<< " idRB= "<< idRB << " mB= "<< mB <<endl;
	  exit(1);
	  return false;
      }
  }


  if(partonA.size()>0) {
    // jet fragmentation.
    //int idRA = isResonance(mA, smA,partonA);

    if(idRA ==  92 ) {
      passStringFrag=1;
      if(!stringFragment(100,partonA,xOut[0],smA,sPot[0],vPot[0],outgoing)) {
        cout << "SoftStrings::setKinematics fragment failed Ecm=" << eCM <<endl;
	return false;
      }
    // resonance
    } else {
     makeHadron(idRA,smA,partonA,xOut[0],tFormA,sPot[0],vPot[0],outgoing);
   }
  } else {
    outgoing.push_back(putHadron(idIn[0],mA,smA,xOut[0],tFormA,pA,constQ[0],sPot[0],vPot[0]));
  }

  /*
      vector<ResolvedPartonQ> parton;
      parton.push_back( ResolvedPartonQ(0,0,1));
      parton.push_back( ResolvedPartonQ(0,0,3201));
     cout << " mjet= "<< mJet(parton) <<endl;
     cout << " res= "<< findIdRes(parton,1.332) <<endl;
	 cin.get();
	 */

  if(partonB.size()>0) {
    //int idRB = isResonance(mB,smB,partonB);
    if(idRB == 92 ) {
	if(mB < 1.0) {
	    cout << "SoftStrings mB = " << mB << "idRB= "<< idRB <<endl;
            cout << " mjet= "<< mJet(partonB) <<endl;
	    Vec4 p=0.0;
	    for(int i=0;i<(int)partonB.size();i++) {
		p += partonB[i].p();
		cout << "id= " << partonB[i].id() <<endl;
	    }
	    cout << " m= "<< p.mCalc() <<endl;
	    exit(1);
	}
      passStringFrag=2;
      if(!stringFragment(200,partonB,xOut[1],smB,sPot[1],vPot[1],outgoing)) {
        cout << "SoftStrings::setKinematics fragment failed Ecm=" << eCM <<endl;
	return false;
      }

    } else {
      makeHadron(idRB,smB,partonB,xOut[1],tFormB,sPot[1],vPot[1],outgoing);
    }
  } else {
    outgoing.push_back(putHadron(idIn[1],mB,smB,xOut[1],tFormB,pB,constQ[1],sPot[1],vPot[1]));
  }

  Vec4 ptot=0.0;
  int ibar=0;
  for(int i=0;i<(int)outgoing.size();i++) {
      ibar += outgoing[i]->baryon();
    ptot += outgoing[i]->getMomentum();
  }
  if(ibar != ibar0) {
      cout << "Softstrings::setKinemtaics  baryon number does not converge ibar0= "<< ibar0
	  << " final = "<< ibar << endl;
      cout << " idbeamA= "<< idBeamA
	  << " idbaemB= "<< idBeamB
	  << " idRA= "<< idRA
	  << " idRB= "<< idRB
	  <<endl;

  for(int i=0;i<(int)outgoing.size();i++) {
    cout << outgoing[i]->getID() << outgoing[i]->getMomentum();
  }

      exit(1);
  }
  if(smA ==1.0 && smB ==1.0 && abs(ptot.mCalc()-eCM)>1e-5) {
      cout << "Softstrings::setKinemtaics  ptot[0]=" << ptot[0] << " ecm=" << eCM<<endl;
      cout << " partonA=" << partonA.size()
	   << " partonB=" << partonB.size()
	   <<endl;
      cout << " mA= "<< mA << " m1= "<< mInEff[0] <<endl;
      cout << " mB= "<< mB << " m2= "<< mInEff[1] <<endl;
  for(int i=0;i<(int)outgoing.size();i++) {
    cout << outgoing[i]->getID() << outgoing[i]->getMomentum();
  }
      exit(1);
  }

  return true;

}

void SoftStrings::calcSigma()
{
  int ibar1 = iBaryon[0];
  int ibar2 = iBaryon[1];

  if(ibar1*ibar2 == 9) { // B+B or antiB+antiB
    //const double minCM= 3.877;
    //if(eCM < minCM) return;
    sigtot.calc(2212,2212,eCM);
  } else if(ibar1 == 0 && ibar2==0) { // M+M
    //const double minCM= 3.551;
    //if(eCM < minCM) return;
    sigtot.calc(211,211,eCM);
  } else if(ibar1*ibar2 ==0) { // M+B
    //const double minCM= 3.714;
    //if(eCM < minCM) return;
    sigtot.calc(211,2212,eCM);
  } else if(ibar1*ibar2 ==-9) { // B+antiB
    //const double minCM= 3.877;
    //if(eCM < minCM) return;
    sigtot.calc(-2212,2212,eCM);
  } else {
    cout << "SoftStrings::calcSigma ibar1= " << ibar1
         << " ibar2= "<< ibar2 << endl;
    exit(1);
  }
  double sigdiff = sigtot.sigmaAX()+sigtot.sigmaXB() + sigtot.sigmaXX()
                 + sigtot.sigmaAXB();
  double sigin = sigtot.sigmaND()+ sigdiff;
  probNonDiff= sigtot.sigmaND() / sigin;
  probSDA= sigtot.sigmaAX() / sigdiff;
  probSDB= sigtot.sigmaXB() / sigdiff;

}

EventParticle* SoftStrings::putHadron(int id1,double m,double sm,Vec4& xout,
	double tform, Vec4& p3, int cq1[2],double spot, Vec4& vpot)
{ 
  ParticleDataEntry* pdr =jamParticleData->find(id1);
  int pid=jamParticleData->pid(id1);
  EventParticle* pa = new EventParticle(id1,pdr);
  pa->setPID(pid);
  pa->setMass(m/sm);
  pa->setCoordinate(xout);
  pa->setVertex(xout);
  p3.rotbst(MfromCM);
  pa->setBoostMatrix(MfromCM);
  pa->setMomentum(p3);
  pa->setOnShell();
  pa->setNumberOfColl(numCollision);
  pa->lastColl(nOperation);
  pa->setConstQuark(cq1);
  pa->setFormationTime(tform);
  double e=pa->getE0();
  double dect = jamParticleData->lifeTime(pdr,m/sm,e);
  //double t = max(xout[0],tform);
  pa->setLifeTime(xout[0]+dect);
  if(withMeanField) {
  if(pa->isMeanField(xout[0],optPotential)) {
    //pa->setPotS(m - m/sm);
    pa->setPotS( spot );
    pa->setPotV( vpot );
  }
  }
  return pa;
}

//--------------------------------------------------------------------------

// Do an inelastic nondiffractive scattering.

bool SoftStrings::nondiff()
{
  // Maximum number of tries to split beam particles before reconnection.
  const int MAXLOOP = 100;
  // Gradually reduce assumed quark masses from their constituent values.
  const double MASSREDUCERATE = 0.02;

  int connect=1;
  if(iBaryon[0]*iBaryon[1]  < 0) connect=2;

  // Check that not stuck in infinite loop. Allow reduced quark masses.
  int    loop = 0;
  double mAbove1, mAbove2;
  Vec4   pc1, pac1, pc2, pac2;
  int itry=0;
  do {
    int ktry=0;
    do {
      if (++loop == MAXLOOP) {
        cout << "Error in LowEnergyProcess::nondiff: "
          " failed to construct valid kinematics" << endl;
        return false;
      }
      double redStep = (loop < 10) ? 1. : exp( -MASSREDUCERATE * (loop - 9));

      // Split up hadrons A  and B into q + qbar or q + qq for meson/baryon.
      splitA( redStep);
      splitB( redStep);

      // Assign relative sharing of longitudinal momentum.
      z1     = splitZ( idc1, idac1, mTc1 / eCM, mTac1 / eCM);
      z2     = splitZ( idc2, idac2, mTc2 / eCM, mTac2 / eCM);
      mT1    = sqrt( mTsc1 / z1 + mTsac1 / (1. - z1));
      mT2    = sqrt( mTsc2 / z2 + mTsac2 / (1. - z2));

      if(ktry++>100) {
	cout << " nondif failed ktry = "<< ktry
	     << " id1= "<< idBeamA << " id2= "<< idBeamB
       << " ecm= "<< eCM
	     <<endl;
	return false;
      }

    // Ensure that hadron beam remnants are not too massive.
    } while (mT1 + mT2 > eCM);

    // Set up kinematics for outgoing beam remnants.
    double e1    = 0.5 * (sCM + mT1 * mT1 - mT2 * mT2) / eCM;
    double pz1   = sqrtpos(e1 * e1 - mT1 * mT1);
    double epz1  = z1 * (e1 + pz1);
    double pzc1  = 0.5 * (epz1 - mTsc1 / epz1 );
    double ec1   = 0.5 * (epz1 + mTsc1 / epz1 );
    pc1.p(   px1,  py1,       pzc1,      ec1 );
    pac1.p( -px1, -py1, pz1 - pzc1, e1 - ec1 );
    double epz2  = z2 * (eCM - e1 + pz1);
    double pzc2  = -0.5 * (epz2 - mTsc2 / epz2 );
    double ec2   =  0.5 * (epz2 + mTsc2 / epz2 );
    pc2.p(   px2,  py2,        pzc2,            ec2 );
    pac2.p( -px2, -py2, -pz1 - pzc2, eCM - e1 - ec2 );


    // Catch reconnected systems with too small masses.
    //mAbove1 = (pc1 + pac2).mCalc() - mThreshold( idc1, idac2);
    //mAbove2 = (pc2 + pac1).mCalc() - mThreshold( idc2, idac1);
  
    if(connect==1) {
    mAbove1 = (pc2 + pac1).mCalc() - mThreshold( idc2, idac1);
    mAbove2 = (pc1 + pac2).mCalc() - mThreshold( idc1, idac2);
    } else {
    mAbove1 = (pc1 + pac1).mCalc() - mThreshold( idc1, idac1);
    mAbove2 = (pc2 + pac2).mCalc() - mThreshold( idc2, idac2);
    }

    if(itry++>100) {
	cout << " nondif failed itry = "<< itry
	     << " id1= "<< idBeamA << " id2= "<< idBeamB
     << " ecm= "<< eCM
	     <<endl;
	return false;
    }
  } while (mAbove1 < 0. || mAbove2 < 0.);

  if(isBaryon[0] && idBeamA < 0) swap(cqBeamA[0],cqBeamA[1]);
  if(isBaryon[1] && idBeamB < 0) swap(cqBeamB[0],cqBeamB[1]);


  // Store new reconnected string systems.
  //leEvent.append(  idc1, 63, 1, 0, 0, 0, 101,   0,  pc1,  mc1);
  //leEvent.append( idac2, 63, 2, 0, 0, 0,   0, 101, pac2, mac2);
  //leEvent.append(  idc2, 63, 2, 0, 0, 0, 102,   0,  pc2,  mc2);
  //leEvent.append( idac1, 63, 1, 0, 0, 0,   0, 102, pac1, mac1);

  if(connect==1) {
    mA = (pc2 + pac1).mCalc();
    mB = (pc1 + pac2).mCalc();

  ResolvedPartonQ q1 = ResolvedPartonQ(cqBeamB[0],0,idc2);
  ResolvedPartonQ q2 = ResolvedPartonQ(cqBeamA[1],0,idac1);
  q1.m( mc2 ); q2.m( mac1 );
  q1.p(pc2); q2.p(pac1);
  q1.cols( 101, 0); // quark
  q2.cols( 0 , 101); // anti-quark or diquark
  partonA.push_back(q1);
  partonA.push_back(q2);

  ResolvedPartonQ q3 = ResolvedPartonQ(cqBeamA[0],0,idc1);
  ResolvedPartonQ q4 = ResolvedPartonQ(cqBeamB[1],0,idac2);
  q3.m( mc1 );
  q4.m( mac2 );
  q3.p(pc1);
  q4.p(pac2);
  q3.cols( 102, 0); // quark
  q4.cols( 0 , 102); // anti-quark or diquark
  partonB.push_back(q3);
  partonB.push_back(q4);

  } else {

  mA = (pc1 + pac1).mCalc();
  mB = (pc2 + pac2).mCalc();
  ResolvedPartonQ q1 = ResolvedPartonQ(cqBeamA[0],0,idc1);
  ResolvedPartonQ q2 = ResolvedPartonQ(cqBeamA[1],0,idac1);
  q1.m( mc1 ); q2.m( mac1 );
  q1.p(pc1); q2.p(pac1);
  q1.cols( 101, 0); // quark
  q2.cols( 0 , 101); // anti-quark or diquark
  partonA.push_back(q1);
  partonA.push_back(q2);

  ResolvedPartonQ q3 = ResolvedPartonQ(cqBeamB[0],0,idc2);
  ResolvedPartonQ q4 = ResolvedPartonQ(cqBeamB[1],0,idac2);
  q3.m( mc2 );
  q4.m( mac2 );
  q3.p(pc2);
  q4.p(pac2);
  q3.cols( 102, 0); // quark
  q4.cols( 0 , 102); // anti-quark or diquark
  partonB.push_back(q3);
  partonB.push_back(q4);

  }

  /*
  cout << " idBeamA= " <<idBeamA << " idc1= "<< idc1 << " idac2= "<< idac2
      << " mA= "<< mA<<endl;
  cout << " idBeamB= " <<idBeamB << " idc2= "<< idc2 << " idac1= "<< idac1
      << " mA= "<< mB<<endl;
      */


  // Done.
  return true;

}

//-------------------------------------------------------------------------

// Find relative momentum of colour and anticolour constituents in hadron.

double SoftStrings::splitZ(int iq1, int iq2, double mRat1, double mRat2) {

  // Initial values.
  int iq1Abs = abs(iq1);
  int iq2Abs = abs(iq2);
  double x1, x2;

  // Handle mesons.
  if (iq1Abs < 10 && iq2Abs < 10) {
    do x1 = pow2( mRat1 + (1. - mRat1) * rndm->flat() );
    while ( pow(1. - x1, xPowMes) < rndm->flat() );
    do x2 = pow2( mRat2 + (1. - mRat2) * rndm->flat() );
    while ( pow(1. - x2, xPowMes) < rndm->flat() );

  // Handle baryons.
  } else {
    if (iq1Abs > 10) swap( mRat1, mRat2);
    double mRat2ab = 0.5 * mRat1 / xDiqEnhance;
    int id1 = (iq2Abs/1000)%10;
    int id2 = (iq2Abs/100)%10;
    int id3 = iq1Abs;
    if(iq1Abs > 10) {
      id1 = (iq1Abs/1000)%10;
      id2 = (iq1Abs/100)%10;
      id3 = iq2Abs;
    }
    double xPow1= id1 == 1 ? xPowBarD : (id1 == 2) ? xPowBarU: xPowBarS;
    double xPow2= id2 == 1 ? xPowBarD : (id2 == 2) ? xPowBarU: xPowBarS;
    double xPow3= id3 == 1 ? xPowBarD : (id3 == 2) ? xPowBarU: xPowBarS;
    //cout << " p1= "<< xPow1 << " p2= "<< xPow2 << " p3= "<< xPow3<<endl;

    double x2a, x2b;
    do x2a = pow2( mRat2ab + (1. - mRat2ab) * rndm->flat() );
    while ( pow(1. - x2a, xPow1) < rndm->flat() );

    do x2b = pow2( mRat2ab + (1. - mRat2ab) * rndm->flat() );
    while ( pow(1. - x2b, xPow2) < rndm->flat() );

    x2 = xDiqEnhance * ( x2a + x2b);

    do x1 = pow2( mRat1 + (1. - mRat1) * rndm->flat() );
    while ( pow(1. - x1, xPow3) < rndm->flat() );

    //if (iq2Abs > 10) swap( x1, x2);
  }

  // Return z value.
  return x1 / (x1 + x2);

}

//-------------------------------------------------------------------------

// Overestimate mass of lightest 2-body state for given flavour combination.
// Only account for one c or b in hadron, and do not consider diquark spin.

double SoftStrings::mThreshold( int iq1, int iq2) {

  // Initial values.
  int iq1Abs = abs(iq1);
  int iq2Abs = abs(iq2);
  if (iq2Abs > 10) swap( iq1Abs, iq2Abs);
  double mThr = 0.14;

  // Mesonic state.
  if (iq1Abs < 10) {
    if (iq2Abs > iq1Abs) swap( iq1Abs, iq2Abs);
    if      (iq1Abs < 3)  mThr += 0.14;
    else if (iq1Abs == 3) mThr += (iq2Abs < 3) ? 0.50 : 1.00;
    else if (iq1Abs == 4) mThr += (iq2Abs < 3) ? 1.90 : 2.00;
    else if (iq1Abs == 5) mThr += (iq2Abs < 3) ? 5.30 : 5.38;

  // Baryonic state.
  } else if (iq2Abs < 10) {
    int iqo1 = (iq1Abs/1000)%10;
    int iqo2 = (iq1Abs/100)%10;
    int iqo3 = iq2Abs;
    if (iqo3 > iqo2) swap( iqo2, iqo3);
    if (iqo2 > iqo1) swap( iqo1, iqo2);
    //if      (iqo1 <  3) mThr += 0.95;
    if      (iqo1 <  3) mThr += 1.5;
    //else if (iqo1 == 3) mThr += (iqo2 < 3) ? 1.12 : ((iqo3 < 3) ? 1.33 : 1.68);
    else if (iqo1 == 3) mThr += (iqo2 < 3) ? 1.28 : ((iqo3 < 3) ? 1.35 : 1.68);
    else if (iqo1 == 4) mThr += (iqo2 < 3) ? 2.30 : ((iqo3 < 3) ? 2.48 : 2.70);
    else if (iqo1 == 5) mThr += (iqo2 < 3) ? 5.62 : ((iqo3 < 3) ? 5.80 : 6.08);

    /*
    cout << " iq1= " << iq1 << " iq2= " << iq2 
	<< " iqo1= " << iqo1 << " iqo2= " << iqo2 
	<< " mThr= " << mThr
	<< endl;
	*/

  // Baryon-antibaryon state.
  } else {
    int iqo1 = (iq1Abs/1000)%10;
    int iqo2 = (iq1Abs/100)%10;
    if      (iqo1 <  3) mThr += 0.95;
    else if (iqo1 == 3) mThr += (iqo2 < 3) ? 1.28 : 1.35;
    else if (iqo1 == 4) mThr += (iqo2 < 3) ? 2.30 : 2.48;
    else if (iqo1 == 5) mThr += (iqo2 < 3) ? 5.62 : 5.80;
    iqo1 = (iq2Abs/1000)%10;
    iqo2 = (iq2Abs/100)%10;
    if      (iqo1 <  3) mThr += 0.95;
    else if (iqo1 == 3) mThr += (iqo2 < 3) ? 1.28 : 1.35;
    else if (iqo1 == 4) mThr += (iqo2 < 3) ? 2.30 : 2.48;
    else if (iqo1 == 5) mThr += (iqo2 < 3) ? 5.62 : 5.80;
  }

  // Done.
  return mThr;

}

double SoftStrings::HadronMass(int id1, int id2)
{
  FlavContainer flav1(id1);
  FlavContainer flav2(id2);

  int idHad=0;
  do {
    idHad   = flavSel->combine( flav1, flav2 );
  } while (idHad == 0);

  return particleData->m0(idHad);
}

double SoftStrings::GSHadronMass(int id)
{
  const double mq[] ={0.0, 0.33, 0.35, 0.5, 1.5, 4.8,171.0};
  int idBeamAbs=abs(id);
  int idnew=id;
  //int id0 = (idBeamAbs)%10;
  if ((idBeamAbs/1000)%10 == 0) {
    int id0 = 1;
    int id1 = (idBeamAbs/100)%10;
    int id2 = (idBeamAbs/10)%10;
    if(id1>3 || id2 > 3 ) return max(2.0,1.001 + mq[id1]+mq[id2]);
    idnew=id1*100+id2*10+id0;

  } else {
    int id1 = (idBeamAbs/1000)%10;
    int id2 = (idBeamAbs/100)%10;
    int id3 = (idBeamAbs/10)%10;
    int id4 = id1*1000+id2*100+id3*10;
    // PDG ID of 1214 does not work for diffractive scattering in pythia.
    if(id4 == 1210) id4=2110;
    if(id4 == 2120) id4=2210;
    if(id4 == 3330) return 2.0; // Omega- does not have resonance.
    int id0 = 2;
    if(id4==1110 || id4==2220 || id4==3330 || id4==4440 || id4==5550) id0=4;
    idnew=id4+id0;
    //idnew=id1*1000+id2*100+id3*10+id0;
  }
  //idnew = id > 0 ? idnew : -idnew;
  return particleData->m0(idnew);
}

//--------------------------------------------------------------------------

// Do an elastic or diffractive scattering.
// type = 2: elastic; = 3: SD (XB); = 4: SD (AX); = 5: DD.

bool SoftStrings::eldiff( int type, Vec4& pA, Vec4& pB)
{
  // Let diffractive mass spectrum begin this much above unexcited mass.
  //const double MDIFFMIN = 0.2;
  //const double MDIFFMIN = 0.25;
  const double MDIFFMIN = 0.30;


  // Gradually reduce assumed quark masses from their constituent values.
  const double MASSREDUCERATE = 0.02;
  //const double MASSREDUCERATE = 0.03;

  // Maximum number of tries to split beam particles before reconnection.
  const int MAXLOOP = 100;

  // Classify process type.
  bool excite1 = (type == 3 || type == 5);
  bool excite2 = (type == 4 || type == 5);

  // Find excited mass ranges.
  mA           = mInEff[0];
  mB           = mInEff[1];
  //double mAmin = (excite1) ? m1 + MDIFFMIN : m1;
  //double mBmin = (excite2) ? m2 + MDIFFMIN : m2;

  int idA=abs(idBeamA);
  int idA1 = (idA/1000)%10;
  int idA2 = (idA/100)%10;
  int idA3 = (idA/10)%10;
  int idB=abs(idBeamB);
  int idB1 = (idB/1000)%10;
  int idB2 = (idB/100)%10;
  int idB3 = (idB/10)%10;
  int nsA=0;
  if(idA1 == 3) nsA += 0.1;
  if(idA2 == 3) nsA += 0.1;
  if(idA3 == 3) nsA += 0.1;
  int nsB=0;
  if(idB1 == 3) nsB += 0.1;
  if(idB2 == 3) nsB += 0.1;
  if(idB3 == 3) nsB += 0.1;

  double m01 = GSHadronMass(idBeamA);
  double m02 = GSHadronMass(idBeamB);
  double mAmin = (excite1) ? m01 + MDIFFMIN + nsA : mInEff[0];
  double mBmin = (excite2) ? m02 + MDIFFMIN + nsB : mInEff[1];

  //int kfa1,kfa2, kfb1,kfb2;
  //hadronContent->findFlavor(idBeamA,kfa1,kfa2); 
  //hadronContent->findFlavor(idBeamB,kfb1,kfb2); 
  //double mq1=particleData->m0( kfa1 )+particleData->m0( kfa2 );
  //double mq2=particleData->m0( kfb1 )+particleData->m0( kfb2 );
  //double mAmin = beam1->isBaryon() ? MminBar + mq1: MminMes+mq1;
  //double mBmin = beam2->isBaryon() ? MminBar + mq2: MminMes+mq2;

  double mAmax = eCM - mBmin;
  double mBmax = eCM - mAmin;
  if (mAmin + mBmin > eCM) {
    if(nError<5) {
    cout << "Error in SoftStrings::eldiff: "
      " too low invariant mass for diffraction eCM= " << eCM
      << " id1= "<< idBeamA << " id2= "<< idBeamB
       << endl;
    cout << "excite1= "<< excite1 << " excite2= "<< excite2<< endl;
    cout << "m1= " << mA << " m01= "<< m01 << " mA= "<< mAmin<<endl;
    cout << "m2= " << mB << " m02= "<< m02 << " mB= "<< mBmin<<endl;
    cout << " mAmin+mBmin= " << mAmin+mBmin<<endl;
     nError++;
    }
    return false;
  }

  int ntry=0;
  double theta=0.0;
  do {

  // Check that not stuck in infinite loop. Allow reduced quark masses.
  int  loop  = 0;
  bool failM = false;
  do {
    failM = false;
    if (++loop == MAXLOOP) {
      cout << "Error in LowEnergyProcess::eldiff: "
        " failed to construct valid kinematics" <<endl;
      return false;
    }
    double redStep = (loop < 10) ? 1. : exp( -MASSREDUCERATE * (loop - 9));

    // Split up hadron 1 (on side A) and assign excited A mass.
    // Should contain better low-mass description??
    if (excite1) {
      splitA(redStep);
      mA     = mAmin * pow( mAmax / mAmin, rndm->flat() );
      if (mTc1 + mTac1 > mA) failM = true;
    }

    // Split up hadron 2 (on side B) and assign excited B mass.
    // Should contain better low-mass description??
    if (excite2 && !failM) {
      splitB(redStep);
      mB     = mBmin * pow( mBmax / mBmin, rndm->flat() );
      if (mTc2 + mTac2 > mB) failM = true;
    }

  // Ensure that pair of hadron masses not too large.
  if (mA + mB > eCM) failM = true;
  } while (failM);

  // Squared masses, energies and longitudinal momenta of excited hadrons.
  double s1    = mInEff[0] * mInEff[0];
  double s2    = mInEff[1] * mInEff[1];
  double sA    = mA * mA;
  double sB    = mB * mB;


  // Select t value and rotate outgoing particles accordingly.
  double lambda12 = pow2( sCM - s1 - s2) - 4. * s1 * s2;
  double lambdaAB = pow2( sCM - sA - sB) - 4. * sA * sB;
  double tLow     = -0.5 * (sCM - (s1 + s2 + sA + sB) + (s1 - s2)
    * (sA - sB) / sCM + sqrtpos(lambda12 *  lambdaAB) / sCM);
  double tUpp     = ( (sA - s1) * (sB - s2) + (s1 + sB - s2 - sA)
    * (s1 * sB - s2 * sA) / sCM ) / tLow;
  double bNow     = bSlope( type);

  if(tLow > -50.0 && tUpp > -50.0) {
    double eBtLow   = exp( bNow * tLow);
    double eBtUpp   = exp( bNow * tUpp);
    double tNow     = log( eBtLow + rndm->flat() * (eBtUpp - eBtLow) ) / bNow;
    theta    = acos( (2. * tNow - tLow - tUpp) / (tUpp - tLow) );
    break;
  }

  double tNow   = log(rndm->flat()) / bNow;
  if( tNow > tLow && tNow < tUpp) {
    theta    = acos( (2. * tNow - tLow - tUpp) / (tUpp - tLow) );
    break;
  }

  if(ntry==100) {
     cout << "SoftStrings::eldif infinite loop?"<<endl;
     exit(1);
  }
  } while(ntry++<100);

  double phi      = 2. * M_PI * rndm->flat();

  //cout << " theta= "<<  theta <<endl;

  double sA    = mA * mA;
  double sB    = mB * mB;
  double eA    = 0.5 * (sCM + sA - sB) / eCM;
  double pzA   = sqrtpos(eA * eA - sA);
  //Vec4   pA( 0., 0.,  pzA,       eA);
  //Vec4   pB( 0., 0., -pzA, eCM - eA);
  pA.p( 0., 0.,  pzA,       eA);
  pB.p( 0., 0., -pzA, eCM - eA);

  // Internal kinematics on side A, boost to CM frame and store constituents.
  if (excite1) {
    double ec1   = 0.5 * (sA + mTsc1 - mTsac1) / mA;
    double pzc1  = sqrtpos(ec1 * ec1 - mTsc1);
    // Diquark always forward. Randomize for meson.
    if ( abs(idac1) > 10 || (abs(idc1) < 10 && abs(idac1) < 10
      && rndm->flat() > 0.5) ) pzc1 = -pzc1;
    Vec4 pc1(   px1,  py1,  pzc1,      ec1);
    Vec4 pac1( -px1, -py1, -pzc1, mA - ec1);

    pc1.bst(pA);
    pac1.bst(pA);
    pc1.rot( theta, phi);
    pac1.rot( theta, phi);

    //leEvent.append(  idc1, 63, 1, 0, 0, 0, 101,   0,  pc1,  mc1);
    //leEvent.append( idac1, 63, 1, 0, 0, 0,   0, 101, pac1, mac1);

    if(isBaryon[0] && idBeamA < 0) swap(cqBeamA[0],cqBeamA[1]);
    ResolvedPartonQ q1 = ResolvedPartonQ(cqBeamA[0],0,idc1);
    ResolvedPartonQ q2 = ResolvedPartonQ(cqBeamA[1],0,idac1);
    q1.m( mc1 );
    q2.m( mac1 );
    q1.p(pc1);
    q2.p(pac1);
    q1.cols( 101, 0); // quark
    q2.cols( 0 , 101); // anti-quark or diquark
    partonA.push_back(q1);
    partonA.push_back(q2);

    /*
    //if(isBaryon[0] && idBeamA < 0) {
  //if((abs(idc1)<1000 && cqBeamA[0] > 1) || (abs(idac1)<1000 && cqBeamA[1]>1)) 
	    cout << "eldif A id1 " << idBeamA << " id2= "<<idBeamB<<endl;
	    cout << " m = "<< (pc1+pac1).mCalc() <<endl;
	    cout << " idc1= "<< idc1 << " idac1= "<< idac1<<endl;
	    cout << " idc2= "<< idc2 << " idac2= "<< idac2<<endl;
	    cout << "eldif cqBeamB0= " << cqBeamB[0]
	         << " cqBeamB1= " << cqBeamB[1]
	         << " cqBeamA0= " << cqBeamA[0]
	         << " cqBeamA1= " << cqBeamA[1]
		 <<endl;
	    cout << " pc1 = "<< pc1;
	    cout << " pac1= "<< pac1;
	    */


  // Simple copy if not excited, and set momentum as in collision frame.
  } else {
    pA.rot( theta, phi);
    //int iNew = leEvent.copy( 1, 63);
    //leEvent[iNew].p( pA);
    //leEvent[iNew].vProd( 0., 0., 0., 0.);
  }

  // Internal kinematics on side B, boost to CM frame and store constituents.
  if (excite2) {
    double ec2   = 0.5 * (sB + mTsc2 - mTsac2) / mB;
    double pzc2  = -sqrtpos(ec2 * ec2 - mTsc2);
    // Diquark always forward (on negative side). Randomize for meson.
    if ( abs(idac2) > 10 || (abs(idc2) < 10 && abs(idac2) < 10
      && rndm->flat() > 0.5) ) pzc2 = -pzc2;
    Vec4 pc2(   px2,  py2,  pzc2,      ec2);
    Vec4 pac2( -px2, -py2, -pzc2, mB - ec2);
    pc2.bst(pB);
    pac2.bst(pB);
    pc2.rot( theta, phi);
    pac2.rot( theta, phi);


    //leEvent.append(  idc2, 63, 2, 0, 0, 0, 102,   0,  pc2,  mc2);
    //leEvent.append( idac2, 63, 2, 0, 0, 0, 0,   102, pac2, mac2);

    if(isBaryon[1] && idBeamB < 0) swap(cqBeamB[0],cqBeamB[1]);
    ResolvedPartonQ q1 = ResolvedPartonQ(cqBeamB[0],0,idc2);
    ResolvedPartonQ q2 = ResolvedPartonQ(cqBeamB[1],0,idac2);
    q1.m( mc2 );
    q2.m( mac2 );
    q1.p(pc2);
    q2.p(pac2);
    q1.cols( 101, 0); // quark
    q2.cols( 0 , 101); // anti-quark or diquark
    partonB.push_back(q1);
    partonB.push_back(q2);

    /*
    //if(isBaryon[1] && idBeamB < 0) {
  //if((abs(idc2)<1000 && cqBeamB[0] > 1) || (abs(idac2)<1000 && cqBeamB[1]>1))
	    cout << "eldif id1 " << idBeamA << " id2= "<<idBeamB<<endl;
	    cout << " m = "<< (pc2+pac2).mCalc() <<endl;
	    cout << " idc1= "<< idc1 << " idac1= "<< idac1<<endl;
	    cout << " idc2= "<< idc2 << " idac2= "<< idac2<<endl;
	    cout << "eldif B cqBeamB0= " << cqBeamB[0]
	         << " cqBeamB1= " << cqBeamB[1]
	         << " cqBeamA0= " << cqBeamA[0]
	         << " cqBeamA1= " << cqBeamA[1]
		 <<endl;
	    cout << " pc2 = "<< pc2;
	    cout << " pac2= "<< pac2;
	    */

  // Simple copy if not excited, and set momentum as in collision frame.
  } else {
    pB.rot( theta, phi);
    //int iNew = leEvent.copy( 2, 63);
    //leEvent[iNew].p( pB);
    //leEvent[iNew].vProd( 0., 0., 0., 0.);
  }

  //for (int i = 3; i < leEvent.size(); ++i) leEvent[i].rot( theta, phi);

  // Done.
  return true;

}

//-------------------------------------------------------------------------

// Split up hadron A into a colour-anticolour pair, with masses and pT values.

bool SoftStrings::splitA(double redMpT) {

  // Split up flavour of hadron into a colour and an anticolour constituent.
  pair< int, int>  paircac  = splitFlav( idBeamA );
  idc1   = paircac.first;
  idac1  = paircac.second;
  if (idc1 == 0 || idac1 == 0) return false;

  // Find constituent masses and scale down to less than full mass.
  mc1    = particleData->m0( idc1);
  mac1   = particleData->m0( idac1);
  double redNow = redMpT * min( 1., mInEff[0] / (mc1 + mac1));
  mc1   *= redNow;
  mac1  *= redNow;

  if(mac1 ==0.0 || mc1==0.0) {
      cout << "id= " << idBeamA << " idc1= "<< idc1 << " idac1= "<< idac1 <<endl;
      cout << " mac1= "<< mac1 << " mc1= "<<mc1<<endl;
      cout << " m1= "<< particleData->m0(idc1)
           << " ma1= "<< particleData->m0(idac1)
	   <<endl;
      cout << " redMpt= " << redMpT<<endl;
      cout << " redNow = " << redNow<<endl;
      exit(1);
  }

  // Select Gaussian relative transverse momenta for constituents.
  pair<double, double> gauss2 = rndm->gauss2();
  px1    = redMpT * sigmaQ * gauss2.first;
  py1    = redMpT * sigmaQ * gauss2.second;
  pTs1   = px1 * px1 + py1 * py1;

  // Construct transverse masses.
  mTsc1  = pow2(mc1)  + pTs1;
  mTsac1 = pow2(mac1) + pTs1;
  mTc1   = sqrt(mTsc1);
  mTac1  = sqrt(mTsac1);

  // Done.
  return true;

}

//-------------------------------------------------------------------------

// Split up hadron B into a colour-anticolour pair, with masses and pT values.

bool SoftStrings::splitB(double redMpT) {

  // Split up flavour of hadron into a colour and an anticolour constituent.
  pair< int, int>  paircac  = splitFlav( idBeamB );
  idc2   = paircac.first;
  idac2  = paircac.second;
  if (idc2 == 0 || idac2 == 0) return false;

  // Find constituent masses and scale down to less than full mass.
  mc2    = particleData->m0( idc2);
  mac2   = particleData->m0( idac2);
  double redNow = redMpT * min( 1., mInEff[1] / (mc2 + mac2));
  mc2   *= redNow;
  mac2  *= redNow;

  if(mac2 ==0.0 || mc2==0.0) {
      cout << "id= " << idBeamB << " idc2= "<< idc2 << " idac2= "<< idac2 <<endl;
      cout << " mac2= "<< mac2 << " mc1= "<<mc2<<endl;
      cout << " m2= "<< particleData->m0(idc2)
           << " ma2= "<< particleData->m0(idac2)
	   <<endl;
      cout << " redMpt= " << redMpT<<endl;
      cout << " redNow = " << redNow<<endl;
      exit(1);
  }

  // Select Gaussian relative transverse momenta for constituents.
  pair<double, double> gauss2 = rndm->gauss2();
  px2    = redMpT * sigmaQ * gauss2.first;
  py2    = redMpT * sigmaQ * gauss2.second;
  pTs2   = px2 * px2 + py2 * py2;

  // Construct transverse masses.
  mTsc2  = pow2(mc2)  + pTs2;
  mTsac2 = pow2(mac2) + pTs2;
  mTc2   = sqrt(mTsc2);
  mTac2  = sqrt(mTsac2);

  // Done.
  return true;

}

//-------------------------------------------------------------------------

// Split up a hadron into a colour and an anticolour part, of q or qq kinds.

pair< int, int> SoftStrings::splitFlav( int id) {

  // Hadron flavour content.
  int idAbs = abs(id);
  int iq1   = (idAbs/1000)%10;
  int iq2   = (idAbs/100)%10;
  int iq3   = (idAbs/10)%10;
  int iq4, iq5;

  // see JPythia::convertID() N(1520)0 = 1214
  int id4 = iq1*1000+iq2*100+iq3*10;
  if(id4==1210) swap(iq1,iq2);
  else if(id4==2120) swap(iq2,iq3);

  // Nondiagonal mesons.
  if (iq1 == 0 && iq2 != iq3) {
    if (id != 130 && id != 310) {
      if (iq2%2 == 1) swap( iq2, iq3);
      if (id > 0) return make_pair( iq2, -iq3);
      else        return make_pair( iq3, -iq2);

    // K0S and K0L are mixes d sbar and dbar s.
    } else {
      if (rndm->flat() < 0.5) return make_pair( 3, -1);
      else                    return make_pair( 1, -3);
    }

  // Diagonal mesons: assume complete mixing ddbar and uubar.
  } else if (iq1 == 0) {
   int id2=iq2*10 + iq3;
   //if (iq2 < 3 || id == 331) {
   if (iq2 < 3 || id2 == 33) {
     iq4 = (rndm->flat() < 0.5) ? 1 : 2;
     // eta and eta' can also be s sbar.
     //if (id == 221 && rndm->flat() < fracEtass) iq4 = 3;
     //if (id == 331 && rndm->flat() < fracEtaPss) iq4 = 3;
     if (id2 == 22 && rndm->flat() < fracEtass) iq4 = 3;
     if (id2 == 33 && rndm->flat() < fracEtaPss) iq4 = 3;
     return make_pair( iq4, -iq4);
   } else {
     return make_pair( iq2, -iq3);
   }

  // Octet baryons.
  } else if (idAbs%10 == 2 || id4==3120 ) {
    // Three identical quarks: emergency in case of higher spin 1/2 multiplet.
    if (iq1 == iq2 && iq2 == iq3) {iq4 = iq1; iq5 = 1100 * iq1 + 3;}
    // Two identical quarks, like normal p or n.
    else if (iq1 == iq2 || iq2 == iq3) {
    // Use SU(6) weight proton=1/3d(uu)1 + 1/6u(ud)1 + 1/2u(ud)0
    //                 nurtron=1/3u(dd)1 + 1/6d(ud)1 + 1/2d(ud)0
      double rr6 = 6. * rndm->flat();
      //if    (iq1 == iq2 && rr6 < 2.) { iq4 = iq3; iq5 = 1100 * iq1 + 3;}
      //else if             (rr6 < 2.) { iq4 = iq1; iq5 = 1100 * iq3 + 3;}
      //else if (rr6 < 3.) { iq4 = iq2; iq5 = 1000 * iq1 + 100 * iq3 + 3;}
      //else               { iq4 = iq2; iq5 = 1000 * iq1 + 100 * iq3 + 1;}

      //double uu1=2.0, ud1=1.0, ud0=3.0;
      //double eps1=1.3, eps2=1.0;

      //double eps1=0.0, eps2=0.0; // exact SU(6)
      double eps1=1.0, eps2=1.0;
      double uu1=2.0+eps1, ud1=1.0+eps2; //ud0=3.0;
      if    (iq1 == iq2 && rr6 < uu1) { iq4 = iq3; iq5 = 1100 * iq1 + 3;}
      else if             (rr6 < uu1) { iq4 = iq1; iq5 = 1100 * iq3 + 3;}
      else if (rr6 < uu1+ud1) { iq4 = iq2; iq5 = 1000 * iq1 + 100 * iq3 + 3;}
      else               { iq4 = iq2; iq5 = 1000 * iq1 + 100 * iq3 + 1;}

    // Three nonidentical quarks, Sigma- or Lambda-like.
    } else {
      int isp = (iq2 > iq3) ? 3 : 1;
      if (iq3 > iq2) swap( iq2, iq3);
      double rr12 = 12. * rndm->flat();
      if      (rr12 < 4.) { iq4 = iq1; iq5 = 1000 * iq2 + 100 * iq3 + isp;}
      else if (rr12 < 5.) { iq4 = iq2; iq5 = 1000 * iq1 + 100 * iq3 + isp;}
      else if (rr12 < 6.) { iq4 = iq3; iq5 = 1000 * iq1 + 100 * iq2 + isp;}
      else if (rr12 < 9.) { iq4 = iq2; iq5 = 1000 * iq1 + 100 * iq3 + 4 - isp;}
      else                { iq4 = iq3; iq5 = 1000 * iq1 + 100 * iq2 + 4 - isp;}
    }
    if (id > 0) return make_pair(  iq4,  iq5);
    else        return make_pair( -iq5, -iq4);

  // Decuplet baryons.
  } else {
    double rr3 = 3. * rndm->flat();
    if (rr3 < 1.)      { iq4 = iq1; iq5 = 1000 * iq2 + 100 * iq3 + 3;}
    else if (rr3 < 2.) { iq4 = iq2; iq5 = 1000 * iq1 + 100 * iq3 + 3;}
    else               { iq4 = iq3; iq5 = 1000 * iq1 + 100 * iq2 + 3;}
    if (id > 0) return make_pair(  iq4,  iq5);
    else        return make_pair( -iq5, -iq4);
  }

  cout << " splitFlav id= " << id <<endl;
  exit(1);

  // Done. (Fake call to avoid unwarranted compiler warning.)
  return make_pair( 0, 0);

}

//-------------------------------------------------------------------------

// Pick slope b of exp(b * t) for elastic and diffractive events.

double SoftStrings::bSlope( int type) 
{

  // Pomeron trajectory alpha(t) = 1 + epsilon + alpha' * t
  const double ALPHAPRIME = 0.25;

  // Steeper slope for baryons than mesons.
  // To do: charm and bottom should have smaller slopes.
  double bA = (isBaryon[0]) ? 2.3 : 1.4;
  double bB = (isBaryon[1]) ? 2.3 : 1.4;

  // Elastic slope.
  if (type == 2)
    return 2. * bA + 2. * bB + 2. * ALPHAPRIME * log(ALPHAPRIME * sCM);

  // Single diffractive slope for XB and AX, respectively.
  if (type == 3) return 2. * bB + 2. * ALPHAPRIME * log(sCM / (mA * mA));
  if (type == 4) return 2. * bA + 2. * ALPHAPRIME * log(sCM / (mB * mB));

  // Double diffractive slope.
  return 2. * ALPHAPRIME * log(exp(4.) + sCM / (ALPHAPRIME * pow2(mA * mB)) );

}

bool SoftStrings::nonDiff()
{

  // Note: kfa1= q, kfa2= qbar or qq
  // in case of anti-baryon kfa1= qq-bar, kfa2=qbar 2019/6/29
  //int kfa1,kfa2, kfb1,kfb2;
  //hadronContent->findFlavor(idBeamA,kfa1,kfa2); 
  //hadronContent->findFlavor(idBeamB,kfb1,kfb2); 

  pair< int, int>  paircacA  = splitFlav( idBeamA );
  int kfa1  = paircacA.first;
  int kfa2  = paircacA.second;
  pair< int, int>  paircacB  = splitFlav( idBeamB );
  int kfb1  = paircacB.first;
  int kfb2  = paircacB.second;

  // how to connect partons.  =1: no flavor exchange. =2: quark exchange.
  int connect=2;
  if(rndm->flat() > probFlavorExchange)  connect=1;

  if(optConstQuarkDiffra==1) {
    connect = (preHadronA  || preHadronB ) ? 1 : 2;
  } else if (optConstQuarkDiffra==2) {
    connect = (preHadronA  && preHadronB ) ? 1 : 2;
  } else if (optConstQuarkDiffra==3) {
    //if((preHadronA || preHadronB) && rndm->flat() < 0.1  ) connect=1;
    //connect = ((preHadronA || preHadronB) && rndm->flat() > 0.1  ) ? 1:2;
    connect = ((preHadronA || preHadronB) ) ? 2:1;
  }

  // in case of flavor exchange involving anti-baryon, anti-quark
  // will be exchanged.
  if(kfa1 < -1000 && connect==2) connect=3;
  if(kfb1 < -1000 && connect==2) connect=3;

  // No quark exchange for antiB B collision 
  if(iBaryon[0]*iBaryon[1] < 0) connect = 1; 

  // produce gluon or not.
  int pg=0;
  //if(eCM > 4.7) pg=1;
  //if(rndm->flat() < 0.7) pg=1;
  // References to beams to simplify indexing.

  Vec4 p1  = Vec4(0., 0.,  pCM, eCM1);
  Vec4 p2  = Vec4(0., 0., -pCM, eCM2);
  beam1->init(eCM,idBeamA,cqBeamA,p1,mInEff[0],kfa1,kfa2,pg);
  beam2->init(eCM,idBeamB,cqBeamB,p2,mInEff[1],kfb1,kfb2,pg);

  double mqA1=particleData->m0( kfa1 );
  double mqA2=particleData->m0( kfa2 );
  double mqB1=particleData->m0( kfb1 );
  double mqB2=particleData->m0( kfb2 );
  if(abs(kfa1)>3) mqA1 += 0.3;
  if(abs(kfa2)>3) mqA2 += 0.3;
  if(abs(kfb1)>3) mqB1 += 0.3;
  if(abs(kfb2)>3) mqB2 += 0.3;

  double mq1=(connect == 2) ? mqB1 + mqA2 :
      (connect == 1) ? mqA1 + mqA2 : mqA1 + mqB2;
  double mq2 = (connect == 2) ? mqA1 + mqB2 :
      (connect == 1) ? mqB1 + mqB2 : mqB1 + mqA2;

  //double mq1=particleData->m0( kfb1 )+particleData->m0( kfa2 );
  //double mq2=particleData->m0( kfa1 )+particleData->m0( kfb2 );

  double mdminA= isBaryon[0] ? MminBar+mq1 : MminMes + mq1;
  double mdminB= isBaryon[1] ? MminBar+mq2 : MminMes + mq2;

  /*
  double mq1=particleData->m0( kfa1 )+particleData->m0( kfa2 );
  double mq2=particleData->m0( kfb1 )+particleData->m0( kfb2 );
  double mdminA=MminBar+mq1;
  double mdminB=MminBar+mq2;
  if(!beam1->isBaryon()) mdminA=MminMes + mq1;
  if(!beam2->isBaryon()) mdminB=MminMes + mq2;
  */

  /*
  double mq1a = particleData->m0( kfa1 );
  double mq2a = particleData->m0( kfa2 );
  double mq1b = particleData->m0( kfb1 );
  double mq2b = particleData->m0( kfb2 );
  double mq1 = max(mq1b, mq2a);
  double mq2 = max(mq1a, mq2b);
  if(connect==1) {
    mq1 = max(mq1a,mq2a);
    mq2 = max(mq1b,mq2b);
  }
  */

  //double mdminA = GSHadronMass(idBeamA) + mq1;
  //double mdminB = GSHadronMass(idBeamB) + mq2;

  //double mdminA = connect == 2 ? mThreshold(kfb1,kfa2):mThreshold(kfa1,kfa2);
  //double mdminB = connect == 2 ? mThreshold(kfa1,kfb2):mThreshold(kfb1,kfb2);

  if(eCM  < mdminA+mdminB + 0.001) {
      cout << " soft string ecm too low "<< eCM
	  << " mdA= " << mdminA
	  << " mdB= " << mdminB
	  << " mq1= " << mq1
	  << " mq2= " << mq2
	  << " id1= " <<idBeamA
	  << " id2= " <<idBeamB
	  <<endl;
      return false;
  }
  BeamHadron& beamA = *beam1;
  BeamHadron& beamB = *beam2;

  // Loop over all subsystems. Default values. Find invariant mass.
  double kTcompSumA   = beamA.size();
  double kTcompSumB   = beamB.size();

  // Primordial kT and compensation power among remnants.
  double kTwidthNow = primordialKTremnant;

  mA=0.0; mB=0.0;
  Vec4 p3, p4;
  int ntry=0;
  do {
  double xSum[2], xInvM[2], w2Beam[2], wPosRem, wNegRem, w2Rem;
  for (int iTry = 0; iTry < NTRYKINMATCH; ++iTry) {

    // Loop over the two beams.
    for (int iBeam = 0; iBeam < 2; ++iBeam) {

      BeamHadron& beam = (iBeam == 0) ? beamA : beamB;
      int nPar = beam.size(); // number of partons in the beam.

      // Generate Gaussian pT for initiator partons inside hadrons.
        double pxSum = 0.;
        double pySum = 0.;
        for (int iPar = 0; iPar < nPar; ++iPar) {
            pair<double, double> gauss2 = rndm->gauss2();
	    if(iTry == 5) kTwidthNow *= 0.5;
	    if(iTry == 15) kTwidthNow =0.01;
	    //kTwidthNow = primordialKTremnant;
	    //if(beam[iPar].id()==21) kTwidthNow=1.3;
            double px = kTwidthNow * gauss2.first;
            double py = kTwidthNow * gauss2.second;
            beam[iPar].px(px);
            beam[iPar].py(py);
            pxSum += px;
            pySum += py;

        }

        // Share recoil between all initiator partons, rescatterers excluded.
        double kTcompSum = (iBeam == 0) ? kTcompSumA : kTcompSumB;

        for (int iPar = 0; iPar < nPar; ++iPar) {
            beam[iPar].px( beam[iPar].px() - pxSum * 1.0 / kTcompSum);
            beam[iPar].py( beam[iPar].py() - pySum * 1.0 / kTcompSum);
        }

      // Pick unrescaled x values for remnants. Sum up (unscaled) p+ and p-.
      xSum[iBeam]  = 0.;
      xInvM[iBeam] = 0.;
      for (int iRem = 0; iRem < nPar; ++iRem) {
        double xPrel = beam.xRemnant( iRem);
        //double xPrel = xSample1( eCM );

        beam[iRem].x(xPrel);
        xSum[iBeam]  += xPrel;
        xInvM[iBeam] += beam[iRem].mT2()/xPrel;
      }

      // Squared transverse mass for each beam, using lightcone x.
      w2Beam[iBeam] = xSum[iBeam] * xInvM[iBeam];

    // End separate treatment of the two beams.
    }


    // Recalculate kinematics of initiator systems with primordial kT.
    wPosRem = eCM;
    wNegRem = eCM;
    w2Rem = wPosRem * wNegRem;

  if (sqrtpos(w2Rem) > sqrt(w2Beam[0]) + sqrt(w2Beam[1])) break;

  }

  if (sqrtpos(w2Rem) < sqrt(w2Beam[0]) + sqrt(w2Beam[1])) {
      cout << "Error in BeamRemnants::setKinematics: idBeamA= " 
	  << idBeamA << " idBeamB= " << idBeamB
	  << " ecm= " << eCM
	  << " w2Rem= " << sqrtpos(w2Rem)
	  << " w2beam1= " << sqrt(w2Beam[0])
	  << " w2beam2= " << sqrt(w2Beam[1])
	  << " minA= " << mdminA
	  << " minB= " << mdminB
	  <<endl;
      return false;
  }

  // Construct x rescaling factors for the two remnants.
  double lambdaRoot = sqrtpos( pow2(w2Rem - w2Beam[0] - w2Beam[1])
    - 4. * w2Beam[0] * w2Beam[1] );
  double rescaleA   = (w2Rem + w2Beam[0] - w2Beam[1] + lambdaRoot)
    / (2. * w2Rem * xSum[0]) ;
  double rescaleB   = (w2Rem + w2Beam[1] - w2Beam[0] + lambdaRoot)
    / (2. * w2Rem * xSum[1]) ;

  // Construct energy and pz for remnants in first beam.
  for (int iRem = 0; iRem < beamA.size(); ++iRem) {
    double pPos = rescaleA * beamA[iRem].x() * wPosRem;
    double pNeg = beamA[iRem].mT2() / pPos;
    beamA[iRem].e( 0.5 * (pPos + pNeg) );
    beamA[iRem].pz( 0.5 * (pPos - pNeg) );
  }

  // Construct energy and pz for remnants in second beam.
  for (int iRem = 0; iRem < beamB.size(); ++iRem) {
    double pNeg = rescaleB * beamB[iRem].x() * wNegRem;
    double pPos = beamB[iRem].mT2() / pNeg;
    beamB[iRem].e( 0.5 * (pPos + pNeg) );
    beamB[iRem].pz( 0.5 * (pPos - pNeg) );
  }

  if(connect==1) {
      p3 = beamA[0].p();
      p4 = beamB[0].p();
      for(int i=1;i<beamA.size();i++) p3 += beamA[i].p();
      for(int i=1;i<beamB.size();i++) p4 += beamB[i].p();
  } else if(connect==2) {

      p3 = beamB[0].p();
      p4 = beamA[0].p();
      for(int i=1;i<beamA.size();i++) p3 += beamA[i].p();
      for(int i=1;i<beamB.size();i++) p4 += beamB[i].p();

      // back reaction?
      /*
      p3 = beamB[kTcompSumB-1].p();
      p4 = beamA[kTcompSumA-1].p();
      for(int i=0;i<beamA.size()-1;i++) p3 += beamA[i].p();
      for(int i=0;i<beamB.size()-1;i++) p4 += beamB[i].p();
      */

  } else {
     p3 = beamB[kTcompSumB-1].p();
     p4 = beamA[kTcompSumA-1].p();
     for(int i=0;i<beamA.size()-1;i++) p3 += beamA[i].p();
     for(int i=0;i<beamB.size()-1;i++) p4 += beamB[i].p();
  }

  mA = p3.mCalc();
  mB = p4.mCalc();

  if(++ntry>1000) {
    cout << "SoftStrings::nonDiff infinite loop?"<<endl;
    cout << "idA= "<< idBeamA << " mA= "<< mA << " idB= "<< idBeamB << " mB= "<< mB<<endl;
    cout << "Ecm= "<< eCM << " mdminA= "<< mdminA << " mdminB= "<< mdminB<<endl;
    return false;
  }

  } while (mA < mdminA || mB < mdminB);

  //osOut << mA << endl;
  //osOut << mB <<endl;

  for(int i=0;i<beamA.size();i++) {
      if(i==0 && connect==2) 
        partonA.push_back(beamB[i]);
      else if(i==kTcompSumA-1 && connect==3) 
        partonA.push_back(beamB[i]);
      else
        partonA.push_back(beamA[i]);
  }
  for(int i=0;i<beamB.size();i++) {
      if(i==0 && connect==2) 
        partonB.push_back(beamA[i]);
      else if(i==kTcompSumA-1 && connect==3) 
        partonB.push_back(beamA[i]);
      else
        partonB.push_back(beamB[i]);
  }

  if(isDebug==5) {
  //if(ibar1 < 0 || ibar2 < 0) {
      cout << "connect= " << connect << endl;
      cout << " id1= " << idBeamA << " id2= " << idBeamB<<endl;
      cout << " cq1= " << cqBeamA[0] << " cq11= " << cqBeamA[1]<< endl;
      cout << " cq2= " << cqBeamB[0] << " cq21= " << cqBeamB[1]<< endl;
      cout << " mA= " << mA;
      Vec4 pA=0.0;
      int nq=0;
      for(int i=0;i<(int)partonA.size();i++)  {
	  pA += partonA[i].p();
	  cout << "idA= " << partonA[i].id()
	       << " cq= " << partonA[i].getConstQ();
	       nq += partonA[i].getConstQ();
      }
      double ma=pA.mCalc();
      cout << " mA= " << mA
	  << " ma= " << ma
	  << " minA= " << mdminA
	  <<endl;

      cout << " mB= " << mB;
      Vec4 pB=0.0;
      for(int i=0;i<(int)partonB.size();i++){
	  pB += partonB[i].p();
	  cout << "idB= " << partonB[i].id()
	       << " cq= " << partonB[i].getConstQ();
	       nq += partonB[i].getConstQ();
      }
      double mb=pB.mCalc();
      cout << " mB= " << mB
	  << " mm= " << mb
	  << " minB= " << mdminB
	  <<endl;

      if(abs(mA-ma)>1e-5 || abs(mB-mb)>1e-5) {
	  cout << " strange mass " <<endl;
	  exit(1);
      }
      if(nq != cqBeamA[0]+cqBeamA[1]+cqBeamB[0]+cqBeamB[1]) {
	  cout << "SoftString::Kintematics const q? " << nq
	      << " nq0 = " <<cqBeamA[0]+cqBeamA[1]+cqBeamB[0]+cqBeamB[1]
	      <<endl;
	  exit(1);
      }

      if(abs(pA[0]+pB[0]-eCM1-eCM2) >1e-8)  {
	  cout << " energy does not conserve " << pA[0]+pB[0]
	      << " e1+e2= "<< eCM1+eCM2
	      <<endl;
	  exit(1);
      }

  }

  return true;
}


void SoftStrings::makeHadron(int idhad,double sm,vector<ResolvedPartonQ>& q,
	Vec4& xout,double tform,double spot, Vec4& vpot,
	vector<EventParticle*>& outgoing)
{
  Vec4 p=0.0;
  for(int i=0;i<(int)q.size();i++)  p += q[i].p();
  double m=p.mCalc();

  ParticleDataEntry* pdr =jamParticleData->find(idhad);

  int pid=jamParticleData->pid(idhad);
  EventParticle* pa = new EventParticle(idhad,pdr);
  pa->setPID(pid);
  pa->setMass(m/sm);
  pa->setCoordinate(xout);
  pa->setVertex(xout);
  p.rotbst(MfromCM);
  pa->setBoostMatrix(MfromCM);
  pa->setMomentum(p);
  pa->setOnShell();
  pa->setNumberOfColl(numCollision);
  pa->lastColl(nOperation);
  int q1=q[0].getConstQ();
  int q2=q.back().getConstQ();
  if(q1==2) std::swap(q1,q2);

  pa->setConstQuark(0,q1);
  pa->setConstQuark(1,q2);
  pa->setFormationTime(tform);
  double e=pa->getE0();
  double dect = jamParticleData->lifeTime(pdr,m/sm,e);
  pa->setLifeTime(max(xout[0],tform)+dect);
  //pa->setLifeTime(xout[0]+dect);

  if(pa->isMeanField(xout[0]+tform,optPotential)) {
    //pa->setPotS(m - m/sm);
    pa->setPotS( spot );
    pa->setPotV( vpot );
  }
  outgoing.push_back(pa);
}

int SoftStrings::findIdRes(vector<ResolvedPartonQ>& parton,double m)
{
  int id1=parton[0].id();
  int id2=parton.back().id();



  FlavContainer flav1(id1);
  FlavContainer flav2(id2);

  int idHad=0;
  do {
    idHad   = flavSel->combine( flav1, flav2 );
  } while (idHad == 0);

  /*
  double emin=particleData->mMin( idHad );
  double wid=particleData->mWidth( idHad );
  cout << " id1= "<< id1 << " id2= "<<id2<< endl;
  cout << " kf= "<< idHad << " m= "<< m << " emin= "<< emin
       << " wid= " << wid
       <<endl;
       */

  int idA=abs(idHad);
  ParticleTable *table=0;
  if(abs(id1)<10 && abs(id2)<10) {
    int kf2 = (idA/100)%10;
    int kf1 = (idA/10)%10;
    int kfm = kf2*10 + kf1;
    //int pid=jamParticleData->pid(id);
    switch (kfm) {
      case 22: table=jamParticleData->getLight0Meson(); break;
      case 33: table=jamParticleData->getLight0Meson(); break;
      case 11: table=jamParticleData->getLight1Meson0(); break;
      case 21: table=jamParticleData->getLight1Mesonp(); break;
      case 31: table=jamParticleData->getStrMeson0(); break;
      case 32: table=jamParticleData->getStrMesonp(); break;
      default: cout << " kfm= "<< kfm << " m= " << m << endl; exit(1);
	       return 0;
    }

  // baryons.
  } else {

    switch (idA) {
      case 2212: table=jamParticleData->getPstar();break;
      case 2112: table=jamParticleData->getNstar();break;
      case 1114: table=jamParticleData->getDmstar();break;
      case 2114: table=jamParticleData->getD0star();break;
      case 2214: table=jamParticleData->getDpstar();break;
      case 2224: table=jamParticleData->getDppstar();break;
      case 3124:
      case 3122: table=jamParticleData->getLambda();break;
      case 3114:
      case 3112: table=jamParticleData->getSmstar();break;
      case 3214:
      case 3212: table=jamParticleData->getS0star();break;
      case 3224:
      case 3222: table=jamParticleData->getSpstar();break;
      case 3314:
      case 3312: table=jamParticleData->getXmstar();break;
      case 3324:
      case 3322: table=jamParticleData->getX0star();break;
      case 3334: return 92;
      default:cout << "findIdRes idA= "<< idA << " m= "<< m << endl; exit(1);
	      return 0;
    }
    // Charm baryon??
  }

  int ntable=table->size();
  vector<double> prob(ntable);
  double total=0.0;
  for(int i=0;i<ntable;i++) {
    ParticleDataEntry* p=table->getParticle(i);
    prob[i]=0.0;
    if(m > p->mMin()) {
     prob[i]=p->spinType();
     total += prob[i];
    }
  }
  if(total == 0.0 ) {
    cout << "SoftStrings::findIdRes No hadron ? id= " << idHad 
         << " m= " << m << endl;
    for(int i=0;i<ntable;i++) {
      ParticleDataEntry* p=table->getParticle(i);
      cout << " m= "<< m << " min= " << p->mMin();
       cout << "id= " << p->id() << " spin= " << p->spinType()<<endl;
    }
    return 0;
  }

  ParticleDataEntry* pdr=0;
  double x =rndm->flat()*total;
  for(int i=0;i<ntable;i++) {
      x -= prob[i];
      if(x<0) {
	  pdr=table->getParticle(i);
	  break;
      }
  }

  int idHadR = idHad > 0 ? pdr->id() : pdr->id()*-1;

  if(idHadR==0) {
      cout << " idHadR= "<< idHadR<<endl;
      exit(1);
  }

  return idHadR;
}

int SoftStrings::isResonance(double m,double sm, vector<ResolvedPartonQ>& q)
{
  if(m >  mJet(q)) return 92;

  int id1=q[0].id();
  int id2=q.back().id();
  int id1a=abs(id1);
  int id2a=abs(id2);
  int ic1=0;
  if((id1a/1000)%10 > 3) ic1++;
  if((id1a/100)%10 > 3) ic1++;
  if(id1a > 3) ic1++;
  if((id2a/1000)%10 > 3) ic1++;
  if((id2a/100)%10 > 3) ic1++;
  if(id2a > 3) ic1++;
  if(ic1 >0) return  92;

  return findIdRes(q, m/sm);
}

// jamemjet()
double SoftStrings::mJet(vector<ResolvedPartonQ>& parton)
{
  const double mMinBar=2.0;
  const double mMinMes=2.0;

  int id1=parton[0].id();
  int id2=parton.back().id();
  double mq = parton[0].m() + parton.back().m();

  if(abs(id1)<10 && abs(id2)<10) {
    return std::max(mMinMes, 1.001 + mq);
  }

  return std::max(mMinBar, 1.0+mq);
  
}

// Deal with the reaction type of M + B -> string or M + M -> string
// B + antiB -> string.
bool SoftStrings::absorbS(int id1, int id2, int id3, const int cq1[2], const int cq2[2],
	 vector<EventParticle*>& outgoing)
{
  pCM=pInican;
  eCM=eCMcan;
  MfromCM=MfromCMcan;
  sCM = eCM * eCM;

  double sm = eCM/(eCM1/SM[0] + eCM2/SM[1]);
  //tFormA = inComing[0]->getTf();
  //tFormB = inComing[1]->getTf();

  int kf1,kf2,cqA,cqB;
  Vec4 xout=xOut[0];
  double spot = sPot[0];
  Vec4 vpot = vPot[0];
  int side=100;

  // B + antiB -> meson string.
  if(id3 ==93) {
    kf1=id1;
    kf2=id2;
    cqA=cq1[0];
    cqB=cq1[1];
    xout = 0.5*(xOut[0] + xOut[1]);
    spot = 0.0;
    vpot = 0.0;
    if(rndm->flat() < 0.5) side=200;	

  } else if(id3 != 92) {
    hadronContent->findFlavor(id3,kf1,kf2); 

    if(iBaryon[1] !=0) {
      side=200;
      xout = xOut[1];
      spot = sPot[1];
      vpot = vPot[1];
    }
    cqA=cq1[0];
    cqB=cq2[1];
    // anti-baryon.
    if(id3<0 && abs(kf1)>100) {
      cqA=cq1[1];
      cqB=cq2[0];
    }

  } else {

    // Note: kfa1= q, kfa2= qbar or qq
    int kfa1,kfa2, kfb1,kfb2;
    int iconv=0, ntry=0;
    while(iconv==0) {
      hadronContent->findFlavor(id1,kfa1,kfa2); 
      hadronContent->findFlavor(id2,kfb1,kfb2); 
      if(kfa1+kfb2==0) iconv=1;
      if(kfa2+kfb1==0) iconv += 2;

      if(ntry++ > 20) {
	cout << "SoftString::absorb no quark-antiquark pair?"<<endl;
	cout << " id1= " << id1 << " kfa1= " << kfa1 << " kfa2= " << kfa2<<endl;
	cout << " id2= " << id2 << " kfb1= " << kfb1 << " kfb2= " << kfb2<<endl;
	exit(1);
	return false;
      }
    }

    if(iBaryon[1] !=0) {
      side=200;
      xout = xOut[1];
      spot = sPot[1];
      vpot = vPot[1];
    }

    kf1=kfb1;
    kf2=kfa2;
    cqA=cq2[1];
    cqB=cq1[0];
    if(iconv==2) {
      kf1=kfa1;
      kf2=kfb2;
      cqA=cq1[0];
      cqB=cq2[1];
    } else if( iconv==3) {
      if(rndm->flat() < 0.5) {
        kf1=kfa1;
        kf2=kfb2;
        cqA=cq1[0];
        cqB=cq2[1];
      }
    }

  }

  double m1=particleData->m0( kf1 );
  double m2=particleData->m0( kf2 );
  ResolvedPartonQ q1 = ResolvedPartonQ(cqA,0,kf1);
  ResolvedPartonQ q2 = ResolvedPartonQ(cqB,0,kf2);
  q1.m( m1 );
  q2.m( m2 );
  double mt1sq=m1*m1;
  double mt2sq=m2*m2;
  double lam0 = pow2(sCM - mt1sq - mt2sq)-4.* mt1sq*mt2sq;
  if(lam0 <= 0.0001) {
    cout << "SoftString::absorbS too small energy" << endl;
    cout << "ecm= " << eCM << " kf1= " << kf1 << " kf2= " << kf2 <<endl;
    return false;
  }

  double px=0.0,py=0.0;
  double lam2=lam0;
  if(primordialKTsChannelString) {
    for (int itry=0; itry<10; itry++) {
      pair<double, double> gauss2 = rndm->gauss2();
      px = sigmaAnn * gauss2.first;
      py = sigmaAnn * gauss2.second;
      //mt1sq = q1.mT2();
      //mt2sq = q2.mT2();
      mt1sq = m1*m1 + px*px + py*py;
      mt2sq = m2*m2 + px*px + py*py;
      lam2 = pow2(sCM - mt1sq - mt2sq)-4.* mt1sq*mt2sq;
      if(lam2 > 0.001) break;
    }
    if(lam2 < 0.001) {
      px=py=0.0;
      lam2=lam0;
      mt1sq=m1*m1;
      mt2sq=m2*m2;
    }
  }

  q1.px(px);  q1.py(py);
  q2.px(-px); q2.py(-py);
  double pz = sqrt(lam2)/(2*eCM);
  if(iBaryon[0] ==0) pz *= -1;
  q1.pz(-pz);
  q2.pz(pz);
  q1.e(sqrt(mt1sq+pz*pz));
  q2.e(sqrt(mt2sq+pz*pz));
  //if(kf1>0 && abs(kf1) < 1000) {
    q1.cols( 101, 0); // quark (anti-diquark)
    q2.cols( 0 , 101); // anti-quark or diquark
  //} else {
  //  q1.cols( 0 , 101); // anti-quark or diquark
  //  q2.cols( 101, 0); // quark
  //}

  Vec4 ptot=q1.p()+q2.p();
  if(abs(eCM-ptot.mCalc())>1e-5) {
  cout << "Absorb ecm= " << eCM
       << " e1+e2= "<< sqrt(mt1sq+pz*pz) + sqrt(mt2sq+pz*pz)
       << " m= "<<ptot.mCalc()
       << " lam0= "<< lam0
       << " lam2= "<< lam2
       <<endl;
      exit(1);
  }

  vector<ResolvedPartonQ> q;
  q.push_back(q1);
  q.push_back(q2);
  passStringFrag=3;
  if(!stringFragment(side,q,xout,sm,spot,vpot,outgoing)) {
    cout << "SoftStrings::absorbS fragmentation faild id1= "<<  id1 << " id2= "<< id2 << " id3= "<< id3
     << " ecm= "<< eCM  <<endl;
    outgoing.clear();
    return false;
  }

  return true;

}



// Deal with the reaction type of B + Bar annihilations.
bool SoftStrings::BaBAnnihilation(vector<EventParticle*>& outgoing)
{
  // suppress Baryon-antiBaryon annihilation: changed to elastic scattering.
  if(noBaBAnnihilation) return false;

  pCM=pInican;
  eCM=eCMcan;
  MfromCM=MfromCMcan;
  sCM = eCM * eCM;

  Vec4 xout = 0.5*(xOut[0]+xOut[1]);
  idBeamA=idIn[0];
  idBeamB=idIn[1];
  int id1 = idBeamA;
  int id2 = idBeamB;
  int cq1[2], cq2[2];
  for(int i=0;i<2;i++) {
    cq1[i]=constQ[0][i];
    cq2[i]=constQ[1][i];
  }

//...Quark contents.
  int kfl1[3],kfl2[3];
  int kf1a=abs(idIn[0]);
  kfl1[2]=(kf1a/1000)%10;
  kfl1[1]=(kf1a/100)%10;
  kfl1[0]=(kf1a/10)%10;

  int kf2a=abs(idIn[1]);
  kfl2[2]=(kf2a/1000)%10;
  kfl2[1]=(kf2a/100)%10;
  kfl2[0]=(kf2a/10)%10;

  int kdiq1[3],kdiq2[3];
  kdiq1[0]= max(kfl1[1],kfl1[2])*10 + min(kfl1[1],kfl1[2]);
  kdiq1[1]= max(kfl1[0],kfl1[2])*10 + min(kfl1[0],kfl1[2]);
  kdiq1[2]= max(kfl1[0],kfl1[1])*10 + min(kfl1[0],kfl1[1]);

  kdiq2[0]= max(kfl2[1],kfl2[2])*10 + min(kfl2[1],kfl2[2]);
  kdiq2[1]= max(kfl2[0],kfl2[2])*10 + min(kfl2[0],kfl2[2]);
  kdiq2[2]= max(kfl2[0],kfl2[1])*10 + min(kfl2[0],kfl2[1]);
  int iann=0;
  int idq1[9],idq2[9];
  for(int i=0;i<3;i++)
  for(int j=0;j<3;j++)
  if(kdiq1[i] ==  kdiq2[j]) {
    idq1[iann]=i;
    idq2[iann]=j;
    iann++;
  }

  int ianti1 = id1 < 0 ?  -1 : 1;
  int ianti2 = id2 < 0 ?  -1 : 1;

  // Set constituent quarks.
  int cqA[3],cqB[3];
  cqA[0]=cq1[0]; cqA[1]=cqA[2]=0;
  if(cq1[1]==2) {
    cqA[1]=1; cqA[2]=1;
  } else if(cq1[1]==1){
    cqA[1]=1; cqA[2]=0;
  }
  if(cq1[0]==2) {
      cout << "SoftStrings::BaBAnnihilation wrong const q1 = " << cq1[0]
	  << " q1[1]= "<< cq1[1]
	  <<endl;
      exit(1);
  }
  if(cq2[0]==2) {
      cout << "SoftStrings::BaBAnnihilation wrong const q2 = "
	  << cq2[0] 
	  << " q2[1]= "<< cq2[1]
	  <<endl;
      exit(1);
  }

  cqB[0]=cq2[0]; cqB[1]=cqB[2]=0;
  if(cq2[1]==2) {
    cqB[1]=1; cqB[2]=1;
  } else if(cq2[1]==1){
    cqB[1]=1; cqB[2]=0;
  }


  // BBar -> one meson string.
  bool oneMeson=false;
  if(iann>0) {
    if(rndm->flat() <= probBBAnn1) oneMeson=true;
    if(eCM < 2.5) oneMeson=true;
  }

  if(oneMeson) {
    int ia=rndm->flat()*iann;
    int kfla=kfl1[idq1[ia]]*ianti1;
    int kflb=kfl2[idq2[ia]]*ianti2;
    //if(kfla < 0) swap(kfla,kflb);

    /*
    if(isDebug) {
	  cout << "SoftStrings::BaBAnnihilation iann= " << iann << " ia= " << ia<< endl;
	  cout << "idq1= "<< idq1[ia] <<endl;
	  cout << "idq2= "<< idq2[ia] <<endl;
	  cout << "kfla= "<< kfla <<endl;
	  cout << "kflb= "<< kflb <<endl;
    }
    */

    int cq[2];
    if(kfla>0) {
      cq[0]=cqA[2];cq[1]=cqB[2];
      return absorbS(kfla,kflb,93,cq, cq2,outgoing);
    } else {
      cq[0]=cqB[2];cq[1]=cqA[2];
      return absorbS(kflb,kfla,93, cq, cq1,outgoing);
    }
  }


  // Bbar + B -> 3 string 
  // Bbar + B -> 2 string + meson
  // Bbar + B -> 2 string 

  // Quark contents and color flow.
  ResolvedPartonQ qA[3], qB[3];
  double sm = eCM/(eCM1/SM[0] + eCM2/SM[1]);
  double emm;
  int kfm=0;
  Vec4 Pmeson=0.0;
  int ns;
  int ntry=0, ntry1=0,ntry2=0,ntry3=0,ntry4=0,ntry5=0,ntry6=0;
  int ntry7=0, ntry8=0;
  do {

    // find flavors.
    ns=findStringBBar(id1,id2,ianti1,ianti2,eCM,kfl1,kfl2,kfm,emm);

    for(int i=0;i<3;i++) {
      qA[i].id(kfl1[i]);
      qB[i].id(kfl2[i]);
    }

    // find momenta of the quarks inside the hadron.
    if(!BBarQuarkMom(ns,eCM1, pCM,emm,1,qA)) {ntry1++;continue;}
    if(!BBarQuarkMom(ns,eCM2,-pCM,emm,2,qB)) {ntry2++;continue;}

    // Convert meson string into hadron
    if(kfm != 0) {
      Pmeson = qA[2].p() + qB[2].p();
      double m3 = Pmeson.m2Calc();
      if(m3 <= 0.0) {ntry3++;continue;}
      m3 = sqrt(m3);
      if(m3 < emm) {ntry4++;continue;}

      double e = sqrt(emm*emm + Pmeson.pAbs2() );
      double de=Pmeson.e()-e;
      if(de <  0.0) {ntry5++;continue;}
      Pmeson.e(e);
      double e1=qA[0].e();
      double e2=qB[0].e();
      qA[0].e(e1+0.5*de);
      qB[0].e(e2+0.5*de);
    }

      // too small string mass.
      if((qA[0].p()+qB[0].p()).mCalc() < mThreshold(qA[0].id(),qB[0].id()))
      {ntry6++;continue;}
      if((qA[1].p()+qB[1].p()).mCalc() < mThreshold(qA[1].id(),qB[1].id())) 
      {ntry7++;continue;}
      if(ns==3 && kfm==0) {
        if((qA[2].p()+qB[2].p()).mCalc() < mThreshold(qA[2].id(),qB[2].id()))
        {ntry8++;continue;}
      }
      break;

  } while (++ntry < 50);

  if(ntry >= 50) {
      cout << "BBAnihilation infinite loop? eCM= " << eCM
           << " id1= "<< id1 << " id2= "<< id2 <<endl;
      cout << "ntry1= "<< ntry1 << " ntry2= "<< ntry2
           << " ntry3= "<< ntry3 << " ntry4= "<< ntry4
           << " ntry5= "<< ntry5 << " ntry6= "<< ntry6
           << " ntry7= "<< ntry7 << " ntry8= "<< ntry8
	  <<endl;
      cout << "ns= "<< ns << " emm= "<< emm << " m1= " << mInEff[0] << " m2= "<< mInEff[1] <<endl;
      return false;
  }


    for(int i=0;i<3;i++) {
      qA[i].setConstQ(cqA[i]);
      qB[i].setConstQ(cqB[i]);

      if(id1>0)  qA[i].cols( 101, 0); // quark
      else qA[i].cols( 0, 101); // anti-quark

      if(id2>0)  qB[i].cols( 101, 0); // quark
      else qB[i].cols( 0, 101); // anti-quark
    }

  //Vec4 ptot=Pmeson;
    // put meson information.
    if(kfm !=0) {
	ParticleDataEntry* pdr =jamParticleData->find(kfm);
	int pid=jamParticleData->pid(kfm);
	EventParticle* pa = new EventParticle(kfm,pdr);
	pa->setPID(pid);
	pa->setMass(emm);
        pa->setCoordinate(xout);
        pa->setVertex(xout);
        Pmeson.rotbst(MfromCM);
	pa->setBoostMatrix(MfromCM);
	pa->setMomentum(Pmeson);
	pa->setOnShell();
	pa->setNumberOfColl(numCollision);
	pa->lastColl(nOperation);

	//pa->setConstQuark(0,cqA[2]);
	//pa->setConstQuark(1,cqB[2]);

	double m=Pmeson.mCalc();
	double e=Pmeson.e();
	double dect = jamParticleData->lifeTime(pdr,m,e);
	double t = 0.5*(xOut[0][0]+xOut[1][0]);
	pa->setLifeTime(t+dect);
        outgoing.push_back(pa);
	ns=2;
    }

    if(isDebug) {
	cout << "BBar annihilation ns= " << ns << endl;
    }

/*
//-----------------------------------------------
  for(int i=0;i<ns;i++) {
    ptot += qA[i].p()+qB[i].p();
  }
  if(abs(eCM-ptot.mCalc())>1e-5) {
  cout << "BaBAnnihilation ecm= " << eCM
       << " m= "<<ptot.mCalc()
       << " kfm= "<< kfm
       << " ns= "<< ns
       <<endl;
      exit(1);
  }
  ptot = 0.0;
//-----------------------------------------------
*/

    // Loop over meson strings
    vector<ResolvedPartonQ> q;
    for(int i=0;i<ns;i++) {
      q.push_back(qA[i]);
      q.push_back(qB[i]);
      int side=100;
      if(rndm->flat() < 0.5) side=200;
      if((qA[i].p()+qB[i].p()).mCalc() < mThreshold(qA[i].id(),qB[i].id())) {
	  cout << "BBar annihilation mass = "<< (qA[i].p()+qB[i].p()).mCalc()
	      << endl;
	  cout << " qA= "<< qA[i].id() << " qB= "<< qB[i].id()<<endl;
	  exit(1);
      }
      passStringFrag=4;
      double spot = 0.0;
      Vec4 vpot = 0.0;
      if(!stringFragment(side,q,xout,sm,spot,vpot,outgoing)) {
        cout <<"SoftStrings::BaBAnnihilation string fragment faild Ecm="<< eCM
	     << " q1= "<< qA[i].id()
	     << " q1.p= "<< qA[i].p()
	     << " q2= "<< qB[i].id()
	     << " q2.p= "<< qB[i].p()
	     << " m= "<< (qA[i].p()+qB[i].p()).mCalc()
	     << endl;
        q.clear();
	return false;
      }
      q.clear();
    }

/*
//-----------------------------------------------
  ptot = 0.0;
  for(auto i : outgoing) {
    ptot += i->getP();
  }
  if(abs(eCM-ptot.mCalc())>1e-5) {
  cout << "BaBAnnihilation outgoing ecm= " << eCM
       << " m= "<<ptot.mCalc()
       << " pass= "<< passStringFrag
       <<endl;
      exit(1);
  }
//-----------------------------------------------
*/

    return true;
}

// Find flavors of meson strings for baryon-anti-baryon annihilation.
int SoftStrings::findStringBBar(int id1, int id2,int ianti1,int ianti2,
	double ecm,
	int kfl1[3],int kfl2[3],int& kfm, double& emm)
{
  int ipion=1; // option for B+antiB -> 2 string + meson.

  int ns=3;  // number of strings.
  kfm=0;    // meson id if antiB B -> strings + meson.
  emm=10.0; // meson mass.

    // Split baryon flavor for side A.
    double xr=rndm->flat();
    int l1,l2,l3;
    if(xr <= 0.33333) {
        l1=1000;
        l2=100;
        l3=10;
    } else if(xr <= 0.66666) {
        l1=100;
        l2=10;
        l3=1000;
    } else {
        l1=10;
        l2=1000;
        l3=100;
    }
    int kf1a=abs(id1);
    kfl1[2]=(kf1a/l1)%10;
    kfl1[1]=(kf1a/l2)%10; 
    kfl1[0]=(kf1a/l3)%10; 

    // Split baryon flavor for side B.
    xr=rndm->flat();
    if(xr <= 0.33333) {
        l1=1000;
        l2=100;
        l3=10;
      } else if(xr <= 0.66666) {
        l1=100;
        l2=10;
        l3=1000;
      } else {
        l1=10;
        l2=1000;
        l3=100;
      }
    int kf2a=abs(id2);
    kfl2[2]=(kf2a/l1)%10; 
    kfl2[1]=(kf2a/l2)%10; 
    kfl2[0]=(kf2a/l3)%10; 

    // Find possible q-qbar pairs for annihilation.
    int nann=0;
    int iann[3][3];
    for(int i=0;i<3;i++)
    for(int j=0;j<3;j++) {
       iann[i][j]=0;
       if(kfl1[i] == kfl2[j]) {
         nann++;
         iann[i][j]=1;
       }
    }

  // put negative sign for anti-quark.
  for(int i=0;i<3;i++) {
    kfl1[i] *= ianti1;
    kfl2[i] *= ianti2;
  }
      
  // Annihilation into two meson strings.
  if(nann >= 1) {
    if(ecm <= 2.3) ns=2;
    if(eCM1<1.5 || eCM2<1.5) ns=2;
    if(rndm->flat() <= probBBAnn2) ns=2;
  }

  if(ns == 2) {
      int kfla=0;
      int kflb=0;
      int ir=1+rndm->flat()*(nann-1);
      int l=0;
      for(int i=0;i<3;i++)
      for(int j=0;j<3;j++) {
       l += iann[i][j];
       if(ir == l) {
           kfla=kfl1[i];
           kflb=kfl2[j];
         if(i == 0) {
           kfl1[0]=kfl1[1];
           kfl1[1]=kfl1[2];
        } else if(i == 1) {
           kfl1[0]=kfl1[0];
           kfl1[1]=kfl1[2];
        }
        if(j == 0) {
             kfl2[0]=kfl2[1];
             kfl2[1]=kfl2[2];
        } else if(j == 1) {
             kfl2[0]=kfl2[0];
             kfl2[1]=kfl2[2];
	}
        kfl1[2]=0;
        kfl2[2]=0;
	goto L100;
	}
      }

L100:
     if(kfla*kflb == 0 || abs(kfla) != abs(kflb)) {

         cout << "SoftStrings::findStringBBar error kfla kflb" << endl;
	 cout << "kfla= " << kfla
	      << " kflb= " << kflb
	      << " ir= " << ir
	      << " nann= " << nann
	      << endl;
         exit(1);
     }
    return ns;
  }

  // antiB + B -> 2 strings + meson.
  //if(ipion==1 && (ecm <= 2.15 || rndm->flat() <= probBBAnn3)) {
  if(ipion==1 && (ecm <= 3.0 || rndm->flat() <= probBBAnn3)) {
        ns=3;

	// Find (smallest) mass of the hadron
       int ilf=4;
       for(int i=0;i<3;i++) {
         //call kfcnst(kfla,kflb,kfh,0.0)
         //call jamdmass(kfh,kfm1,kfd,emin,emdn,ist)

	 FlavContainer flav1(kfl1[i]);
	 FlavContainer flav2(kfl2[i]);
	 int kfh=flavSel->combine(flav1,flav2);
	 double emin=particleData->m0( kfh );
         if(emin <  emm) {
           kfm=kfh;
           emm=emin;
           ilf=i;
	 }
       }
       if(ilf == 4) {
	   cout <<  " funny ilf=0 kfm= " << kfm <<endl; 
	   cout << " emm= " << emm << endl;
	   exit(1);
       }

       if(ilf != 2) {
	  int ktmp=kfl1[2];
	  kfl1[2]=kfl1[ilf];
	  kfl1[ilf]=ktmp;
	  ktmp=kfl2[2];
	  kfl2[2]=kfl2[ilf];
	  kfl2[ilf]=ktmp;
       }

    }

    return ns;
}

//***********************************************************************
//      subroutine jamdivide(ns,pe,pz,jt,kfl,p,icon)

// see pythia6 PYREMN

bool SoftStrings::BBarQuarkMom(int ns, double pe, double pz, 
	double emm, int jt,
	ResolvedPartonQ q[3])
{
    double m1=particleData->m0( q[0].id() );
    double m2=particleData->m0( q[1].id() );
    q[0].m( m1 );
    q[1].m( m2 );
    double mq = m1 + m2;
    double m3 =0.0;
    if(ns==3) {
	m3 = emm < 10.0 ?  particleData->m0( q[2].id() ) : 0.0;
	q[2].m( m3 );
      mq += m3;
    }
    if(pe < mq + 0.001 ) {
	cout << "BBarQuarkMom energy too small pe= "<< pe
	    << " eCM= " << eCM
	    << " ns= " << ns
	    << " jt= " << jt
	    << " idA= " << idBeamA
	    << " idB= " << idBeamB;
        cout << " m= " << mq
	    << " m1= "<< m1
	    << " m2= "<< m2
	    << " q0= "<< q[0].id()
	    << " q1= "<< q[1].id();
	if(ns==3) {
	    cout << " m3= "<< m3
		<< " q2= " << q[2].id();
	}
	cout  <<endl;

	return false;
    }

    double kt =sigmaQBBar;
    double chi1=0.0,chi2=0.0;
    double mt,mt1sq,mt2sq;
    double esq=pe*pe;
    int ntry=0;
    do {
	// first sample transverse momentum for string 1.
        pair<double, double> gauss2 = rndm->gauss2();
        double px = kt * gauss2.first;
        double py = kt * gauss2.second;
        q[0].px(px);  q[0].py(py);
        mt1sq = q[0].mT2();

	// transverse momentum for string 2.
	if(ns==2) {
          q[1].px(-px); q[1].py(-py);
          mt2sq = q[1].mT2();
          chi1=1.0-sqrt(rndm->flat());
          mt=mt1sq/chi1 + mt2sq/(1.0-chi1);

	// transverse momentum for string 2 and 3.
	} else {
          pair<double, double> gauss3 = rndm->gauss2();
          double px3 = kt * gauss3.first;
          double py3 = kt * gauss3.second;
          q[1].px(px3);     q[1].py(py3);
          q[2].px(-px-px3); q[2].py(-py-py3);
          mt2sq = q[1].mT2();
          double mt3sq = q[2].mT2();
	  do {
            chi1=1.0-sqrt(rndm->flat());
            chi2=1.0-sqrt(rndm->flat());
	  } while(1.0-chi1-chi2 < 0.0);
          mt=mt1sq/chi1+mt2sq/chi2+mt3sq/(1.0-chi1-chi2);
	}
	if(ntry==5) kt *= 0.5;
	if(ntry==10) kt = 0.0;
	if(ntry++ > 20) {
	    cout << " BBarQuarkMom infinite loop? kt= " << kt
		 << " ns= " << ns
		 << " mq= " << mq
		 << " mt= " << sqrt(mt)
		 << " mt1= " << sqrt(mt1sq)
		 << " mt2= " << sqrt(mt2sq)
		 << " pe= " << pe << " eCM= " << eCM << endl;
	    return false;
	}

    } while(mt > esq);

//...Subdivide longitudinal momentum according to value selected above.
    double pw1=chi1*(pe+abs(pz));
    q[0].e(0.5*(pw1+mt1sq/pw1));
    q[0].pz(0.5*(pw1-mt1sq/pw1)*pow(-1,jt-1));
    if(ns==2) {
	q[1].e(pe  - q[0].e());
	q[1].pz(pz - q[0].pz());
    } else {
      double pw2=chi2*(pe+abs(pz));
      q[1].e(0.5*(pw2+mt2sq/pw2));
      q[1].pz(0.5*(pw2-mt2sq/pw2)*pow(-1,jt-1));
      q[2].e(pe  - q[0].e() -q[1].e());
      q[2].pz(pz - q[0].pz() -q[1].pz());

      //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
              if(q[2].e()<0.0) {
          cout << "2 quark energy<0 in BBann q1e= "<< q[0].e()
            << " q2e= "<< q[1].e()
            << " q3e= "<< q[2].e()
            << " E= "<< pe
            <<endl;
	  return false;
        }
      //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    }

  return true;

}

double SoftStrings::xSample3(double srt,double xmin, double xmax)
{
//...p(x)=1/sqrt(x^2+(c/srt)^2)
    double cut=pow2(0.1/srt);
    double cut1=xmax+sqrt(xmax*xmax+cut);
    double cut2=xmin+sqrt(xmin*xmin+cut);
    double ar=pow(cut2*(cut1/cut2), rndm->flat());
    return 0.5*(ar-cut/ar);
}

// t-channel string excitation according to the HIJING model.
bool SoftStrings::hijingSoft(int type, vector<EventParticle*>& outgoing)
{
  int id1=idIn[0];
  int id2=idIn[1];
  Vec4 xoutA = xOut[0];
  Vec4 xoutB = xOut[1];
  int cq1[2], cq2[2];
  for(int i=0;i<2;i++) {
    cq1[i]=constQ[0][i];
    cq2[i]=constQ[1][i];
  }

  int diffra=0;
  if(type >=3) {
      diffra=type-2;
  } else {
    diffra=0;
    if(rndm->flat() < probDiffra) {
      diffra=1;
      if(rndm->flat() < 0.5) diffra=2;
    }
  }

  int kfa1,kfa2,kfb1,kfb2;
  hadronContent->findFlavor(id1,kfa1,kfa2); 
  double qma1 = particleData->m0( kfa1 );
  double qma2 = particleData->m0( kfa2 );

  hadronContent->findFlavor(id2,kfb1,kfb2); 
  double qmb1 = particleData->m0( kfb1 );
  double qmb2 = particleData->m0( kfb2 );

  double mdminA=MminBar+qma1+qma2;
  double mdminB=MminBar+qmb1+qmb2;
  if((abs(id1)/1000)%10==0) mdminA=MminMes+qma1+qma2;
  if((abs(id2)/1000)%10==0) mdminB=MminMes+qmb1+qmb2;
  if(diffra==1) mdminA=mInEff[0];
  if(diffra==2) mdminB=mInEff[1];

//...Initial light-cone momenta p+,p- of proj. and targ.
    double epp=eCM1 + pCM;
    double epm=eCM1 - pCM;
    double etp=eCM2 - pCM;
    double etm=eCM2 + pCM;

//...Total W+,W- and center-of-mass energy
    double wp=epp+etp;
    double wm=epm+etm;
    double sw=wp*wm;
    double srt=sqrt(sw);
    double pkcmx2=(pow2(sw-pow2(mdminA)-pow2(mdminB))-pow2(2.0*mdminA*mdminB))
	/(4.0*sw);

    double widg=0.65;
    double widg2=widg*widg;

    double d1,d2,bb1,bb2,pkc2;
    int ntry=0;
L100:
    ntry++;
    if(ntry>20) {widg2=0.1;}
    if(ntry>100) {
	cout << "SoftStrings::hijingSoft do not converge" <<endl;
	cout << "eCM= " << srt << " id1= " << id1 << " id2= " << id2 << endl;
	return false;
    }

    int itry=0;
    do {
	itry++;
	if(itry>20) widg2=0.1;
	if(itry>100) {
	cout << "SoftStrings::hijingSoft do not converge itry=" << itry  <<endl;
	cout << "eCM= " << srt << " id1= " << id1 << " id2= " << id2 << endl;
	return false;
	}

    pkc2=widg2*(-log(1.0- rndm->flat()*(1.0-exp(-pkcmx2/widg2))));

    //double ampx=sqrt(mdminA*mdminA+pkc2);
    //double amtx=sqrt(mdminB*mdminB+pkc2);
    //double dpx=ampx*ampx/sw;
    //double dtx=amtx*ampx/sw;

    //double dpx=(mdminA*mdminA+pkc2)/sw;
    //double dtx=(mdminB*mdminB+pkc2)/sw;

    d1=(mdminA*mdminA+pkc2)/sw;
    d2=(mdminB*mdminB+pkc2)/sw;

    /*
//...Have single diffractive collision fix projectile mass.
     if(diffra== 1) {
          d1=(m1*m1+pkc2)/sw;
          d2=dtx;
          //bb2=1.0+d2-d1;
	  
     } else if(diffra== 2) {
          d1=dpx;
          d2=(m2*m2+pkc2)/sw;
          //bb2=1.0+d2-d1;
    } else {
	d1=dpx;
	d2=dtx;
    }
    */

    bb1=1.0+d1-d2;
    bb2=1.0+d2-d1;
    } while (bb1*bb1 < 4.0*d1 ||  bb2*bb2 < 4.0*d2);


    double x1,x2;
    if(diffra==1) {
      double xmin=(bb2-sqrt(bb2*bb2-4.0*d2))/2.0;
      double xmax=(bb2+sqrt(bb2*bb2-4.0*d2))/2.0;
      x2=xSample3(srt,xmin,xmax);
      x1=d1/(1.0-x2);

      if(x2*(1.0-x1) < (d2+1.e-4/sw) ) goto L100;

    } else if(diffra==2) {
      double xmin=(bb1-sqrt(bb1*bb1-4.0*d1))/2.0;
      double xmax=(bb1+sqrt(bb1*bb1-4.0*d1))/2.0;
      x1=xSample3(srt,xmin,xmax);
      x2=d2/(1.0-x1);
      if(x1*(1.0-x2) < (d1+1.e-4/sw)) goto L100;

    }else {

      double xmin1=(bb1-sqrt(bb1*bb1-4.0*d1))/2.0;
      double xmax1=(bb1+sqrt(bb1*bb1-4.0*d1))/2.0;
      double xmin2=(bb2-sqrt(bb2*bb2-4.0*d2))/2.0;
      double xmax2=(bb2+sqrt(bb2*bb2-4.0*d2))/2.0;
      x1=xSample1(srt,xmin1,xmax1);
      x2=xSample1(srt,xmin2,xmax2);
      double xxp=x1*(1.0-x2);
      double xxt=x2*(1.0-x1);
      if(xxp < (d1+1.e-4/sw) || xxt < (d2+1.e-4/sw)) goto L100;
    }

    double  epp1=(1.0-x2)*wp;
    double  epm1=x1*wm;
    double  etp1=x2*wp;
    double  etm1=(1.0-x1)*wm;

    double mAsq=epp1*epm1-pkc2;
    double mBsq=etp1*etm1-pkc2;
    if(mAsq <= 0.0) goto L100;
    if(mBsq <= 0.0) goto L100;

    double m1=sqrt(epp1*epm1-pkc2);
    double m2=sqrt(etp1*etm1-pkc2);
    if(m1 < mdminA) goto L100;
    if(m2 < mdminB) goto L100;
    if(srt < m1 + m2) goto L100;

    Vec4 p3, p4;
//...Set proj. energy, momentum and flavor.
    p3.pz((epp1-epm1)/2.0 );
    p3.e( (epp1+epm1)/2.0 );

//...Set targ. energy, momentum and flavor.
    p4.pz((etp1-etm1)/2.0);
    p4.e( (etp1+etm1)/2.0);

    double phi=2.0*M_PI*rndm->flat();
    double pkc=sqrt(pkc2);
    double pkc11=pkc*cos(phi);
    double pkc12=pkc*sin(phi);
    p3.px(pkc11);
    p3.py(pkc12);
    p4.px(-pkc11);
    p4.py(-pkc12);

      vector<ResolvedPartonQ> q3, q4;
      int sgn=1;
      //if(rndm->flat()< 0.5) sgn=-1;
      if(diffra != 1) 
	  if(!partonMomentum(sgn,id1,kfa1,kfa2,qma1,qma2 ,p3,q3,cq1)) goto L100;
      if(diffra != 2) 
	  if(!partonMomentum(-sgn,id2,kfb1,kfb2,qmb1,qmb2,p4,q4,cq2)) goto L100;


      if(isDebug) {
      cout << " diffra = " << diffra <<endl;
      cout << " id1= " << id1 << " id2= " <<id2 <<endl;
      cout << " m1= " << m1 << " m2= " << m2 <<endl;
      }


      if(diffra==1) {

        // jet fragmentation.
        int idRB = isResonance(m2,SM[1], q4);
        if(idRB == 92 ) {
          passStringFrag=5;
          if(!stringFragment(200,q4,xoutB,SM[1],sPot[1],vPot[1],outgoing)) {
            cout << "SoftStrings::hijingSoft string fragement faild"<<endl;
	    return false;
	  }
        // resonance
        } else {
          makeHadron(idRB,SM[1],q4,xoutB,tFormB,sPot[1],vPot[1],outgoing);
        }

	   // stringFragment(q4,xoutB, outgoing);

	    ParticleDataEntry* pdr =jamParticleData->find(id1);
	    int pid=jamParticleData->pid(id1);
	    EventParticle* pa = new EventParticle(id1,pdr);
	    pa->setPID(pid);
	    pa->setMass(m1);
	    pa->setCoordinate(xoutA);
	    pa->setVertex(xoutA);
            p3.rotbst(MfromCM);
	    pa->setBoostMatrix(MfromCM);
	    pa->setMomentum(p3);
	    pa->setOnShell();
	    pa->setNumberOfColl(numCollision);
	    pa->lastColl(nOperation);
	    pa->setConstQuark(cq1);
	    double m=p3.mCalc();
	    double e=p3.e();
	    double dect = jamParticleData->lifeTime(pdr,m,e);
	    pa->setLifeTime(xoutA[0]+dect);
            outgoing.push_back(pa);

      } else if(diffra==2) {

        // jet fragmentation.
        int idRA = isResonance(m1, SM[0], q3);
        if(idRA== 92 ) {
          passStringFrag=6;
          if(!stringFragment(100,q3,xoutA,SM[0],sPot[0],vPot[0],outgoing)) {
            cout << "SoftStrings::hijingSoft string fragement faild"<<endl;
	    return false;
	  }
        // resonance
        } else {
          makeHadron(idRA,SM[0],q3,xoutA,tFormA,sPot[0],vPot[0],outgoing);
        }

	    //stringFragment(q3,xoutA, outgoing);

	    ParticleDataEntry* pdr =jamParticleData->find(id2);
	    int pid=jamParticleData->pid(id2);
	    EventParticle* pa = new EventParticle(id2,pdr);
	    pa->setPID(pid);
	    pa->setMass(m2);
	    pa->setCoordinate(xoutB);
	    pa->setVertex(xoutB);
            p4.rotbst(MfromCM);
	    pa->setBoostMatrix(MfromCM);
	    pa->setMomentum(p4);
	    pa->setOnShell();
	    pa->setNumberOfColl(numCollision);
	    pa->lastColl(nOperation);
	    pa->setConstQuark(cq2);
	    double m=p4.mCalc();
	    double e=p4.e();
	    double dect = jamParticleData->lifeTime(pdr,m,e);
	    pa->setLifeTime(xoutB[0]+dect);
            outgoing.push_back(pa);
      } else {

        int idRA = isResonance(m1,SM[0], q3);
        if(idRA== 92 ) {
            passStringFrag=7;
	    if(!stringFragment(100,q3,xoutA,SM[0],sPot[0],vPot[0],outgoing)) {
              cout << "SoftStrings::hijingSoft string fragement faild"<<endl;
	      return false;
	    }
	} else {
	    makeHadron(idRA,SM[0],q3,xoutA,tFormA,sPot[0],vPot[0],outgoing);
	}

        int idRB = isResonance(m2,SM[1], q4);
        if(idRB ==  92 ) {
            passStringFrag=8;
          if(!stringFragment(200,q4,xoutB,SM[1],sPot[1],vPot[1],outgoing)) {
              cout << "SoftStrings::hijingSoft string fragement faild"<<endl;
	      return false;
	  }
        } else {
          makeHadron(idRB,SM[1],q4,xoutB,tFormB,sPot[1],vPot[1],outgoing);
        }
	    //stringFragment(q3,xoutA,outgoing);
	    //stringFragment(q4,xoutB,outgoing);
      }

  return true;

}

bool SoftStrings::partonMomentum(int sgn,int id, int kfa1,int kfa2,
	double qm1, double qm2,
	Vec4& p4, vector<ResolvedPartonQ>& q, int cq[2])
{
    // Construct pz of quarks and diquarks
    double m2=p4.mCalc();
    double pkc2 = p4.pT2();
    double emt=sqrt(m2*m2 + pkc2);
    double emt1=sqrt(qm1*qm1 + pkc2/2.0);
    double emt2=sqrt(qm2*qm2 + pkc2/2.0);
    if(emt < emt1+emt2) {
	cout << "SoftStrings::partonMomentum emt < emt1+emt2 m= " << m2 <<endl;
	cout << " pkc= " << sqrt(pkc2)
	    <<endl;
	return false;
    }
    int cq1=cq[0];
    int cq2=cq[1];
    // anti-baryon anti-qq - q.
    if(id<0 && cq2==2) {
      cq1=cq[1];
      cq2=cq[0];
    }

    ResolvedPartonQ q1 = ResolvedPartonQ(cq1,0,kfa1);
    ResolvedPartonQ q2 = ResolvedPartonQ(cq2,0,kfa2);
    q1.m( qm1 );
    q2.m( qm2 );
    q1.px(p4.px()/2.0); q1.py(p4.py()/2.0);
    q2.px(p4.px()/2.0); q2.py(p4.py()/2.0);
    q1.cols( 101, 0); // quark
    q2.cols( 0 , 101); // anti-quark or diquark

    double pzcm=PCM(emt,emt1,emt2);
    pzcm *= sgn;
    double gg=p4.e()/emt;
    double bet=p4.pz()/p4.e();
    double e1=sqrt(emt1*emt1+pzcm*pzcm);
    q1.pz( gg*(-pzcm+bet*e1) );  // pz of quark 1
    q1.e(  gg*(e1-bet*pzcm) );   // e  of quark 1

    double e2=sqrt(emt2*emt2+pzcm*pzcm);
    q2.pz(gg*(pzcm+bet*e2));  // pz of diquark 1
    q2.e(gg*(e2+bet*pzcm));  // e  of diquark 1

    q1.m(q1.p().mCalc());
    q2.m(q2.p().mCalc());

    q.push_back(q1);
    q.push_back(q2);

    return true;
}

//----------------------------------------------------------------------
bool SoftStrings::stringFragment(int side,vector<ResolvedPartonQ>& q,
	Vec4& xout, double sm, double spot, Vec4& vpot,vector<EventParticle*>& outgoing)
{
  Event& event = hadronize->event;
  event.reset();
  vector<EventParticle*> hadrons;
  vector<pair<int,int>> constq(1,make_pair(0,0));

  sm=1.0;

  // first count the number of constituent quarks.
  int nval=0;
  int isbaryon=0;
  int nconstq=0;
  Vec4 ptot0=0.0;
  int chargetot0=0;
  for(int i=0;i<(int)q.size();i++) {
    chargetot0 += particleData->chargeType(q[i].id());

    // transform partons to lab frame.
    //Vec4 p = q[i].p();
    //p.rotbst(MfromCM);
    //q[i].p(p);

    ptot0 += q[i].p();
    //q[i].p(q[i].p()/sm);

    if(abs(q[i].id()) != 21) {
      if (abs(q[i].id()) > 1000) isbaryon = q[i].id() > 0 ? 1 : -1;
      nconstq += q[i].getConstQ();
      nval++;
    }
  }

  // Check errors. -------------------------------------------------------
  int ierror=0;
  for(int i=0;i<(int)q.size();i++) {
    if(abs(q[i].id()) < 1000 && abs(q[i].getConstQ()) > 1) {
      ierror=1;
      cout << "before fragment string constq = " << q[i].getConstQ()
	        << " id= "<< q[i].id() <<endl;
    }
  }
  if(ierror) {
    for(int i=0;i<(int)q.size();i++) {
      cout << " id= "<< q[i].id() << " q= " << q[i].getConstQ()<<endl;
    }
   //exit(1);
  }
//------------------------------------------------------------------------

  makeString(nval,q,constq,event);
  Event event0=event;

  // Time-like shower
  double strmass=ptot0.mCalc();
  if(strmass < 0.28) {
    cout << "SoftString::stringFragmentaion string mass too small? " << strmass <<endl;
    event.list();
  }
  if(doShower && isbaryon == 1 && !preHadronA && !preHadronB) {
    double scale = strmass;
    for(int i=1;i<event.size();i++) event[i].scale(scale);
    hadronize->forceTimeShower(1,event.size()-1,scale);
  }
  Event event1=event;

  bool first=true;
  for (; ; ) {
  if(!hadronize->next()) {

    Vec4 ptot1=0.0;
    int charge1=0;
    for(int i=1; i<event.size();i++)
    if(event[i].isFinal()) {
        if(event[i].isParton()) {
	  cout << "SoftStrings::stringFragment hadronzation failed parton? i= " << i << " id= "<< event[i].id() <<endl;
	  event0.list();
	  event1.list();
	  event.list();
	  if(first) {
	    event = event0;
	    first=false;
	    continue;
	  } else {
	    event.list();
	    cout << "SoftStrings::stringFragment hadronzation giving up?" <<endl;
	    return false;
	  }
	}
	ptot1 += event[i].p()*sm;
	charge1 += event[i].chargeType();
    }
    if(abs(chargetot0-charge1) !=0) return false;

    Vec4 diff= ptot1-ptot0;
    double pdiff = abs(diff.mCalc());
    if(pdiff  > 1e-5) {
	cout << " after string fragmentation pdiff = " << pdiff <<endl;
        return false;
    }
  }
    break;
  }

  // boost to the computational frame.
  event.rotbst(MfromCM);

  // Find leading hadrons.
  constq.resize(event.size());
  findConstQuark(event,isbaryon,nconstq,constq);

  if(isDebug>1) {
    event.list();
    cout << " after findConstQuark nconstq="<< nconstq <<endl;
    cout << " isbaryon= "<< isbaryon<<endl;
    cout << " q0= "<< q[0].getConstQ() << " id0= " << q[0].id()
         << " q1= "<< q.back().getConstQ() << " id1= " << q.back().id()
	 <<endl;
  }

    int kaon=0;
    int nq=0;
    int hasMeanField=0;
    double etot0=0.0;
    //EventParticle* had=0;
    Vec4 ptot = 0.0;
    for(int i=1; i<event.size();i++)
	if(event[i].isFinal()) {
	    int id=event[i].id();

	    if(id==130 || id==310) {

		if(kaon==0) {
		    if(rndm->flat() < 0.5) {
			id=311; kaon=1;
		    } else {
			id=-311; kaon=-1;
		    }
		} else if(kaon==1) {
		    id=-311; kaon=0;
		} else if(kaon== -1) {
		    id=311; kaon=0;
		}
	    }
	    ParticleDataEntry* pdr =jamParticleData->find(id);
	    //ParticleDataEntry* pd =&(event[i].particleDataEntry());

	    int pid=jamParticleData->pid(id);
	    EventParticle* pa = new EventParticle(id,pdr);
	    pa->setPID(pid);
            pa->setParent(92);
	    Vec4 p=event[i].p()*sm;
	    ptot += p;
	    double m = event[i].m();
	    pa->setMass(m);
	    pa->setMomentum(p);
	    pa->setOnShell();
	    pa->setNumberOfColl(numCollision);
            if(allowRescatterSameString==2) pa->lastColl(-1);
	    else pa->lastColl(nOperation);
	    pa->setBoostMatrix(MfromCM);
	    pa->setParent(side);

	    // original position.
	    Vec4 v=xout;
	    // formation point.
	    Vec4 xprod = event[i].vProd()*MM2FM;
            double tform=event[i].tProd()*MM2FM;
	    //double tf = v[0]+tform;

	    if(optConstFormationTime) {
	      tform=1.0 * p[0]/m;
	      xprod[1] = tform*p[1]/p[0];
	      xprod[2] = tform*p[2]/p[0];
	      xprod[3] = tform*p[3]/p[0];
	      xprod[0] = tform;
	      //tf = v[0] + tform;
	    }

	    // Find leading hadron.
	    int isq=0;
	    if(constq[i].first > 0 || constq[i].second > 0) {
		isq=1;
		pa->setConstQuark(0,constq[i].first);
		pa->setConstQuark(1,constq[i].second);
	    } else {
		pa->setConstQuark(0,0);
		pa->setConstQuark(1,0);
	    }

	    int *cqq = pa->constQuark();
	    nq+=cqq[0]+cqq[1];

	    if((pa->baryon() ==0 && (cqq[0] >1 || cqq[1]>1))||
		(pa->baryon()== -3 &&(cqq[0] >1 )))
	    {
		cout << " pass string = "<< passStringFrag <<endl;
		cout << i << " fragment strange cq id= " << pa->getID()
		    << " cq0= "<< cqq[0]
		    << " cq1= "<< cqq[1]
		    <<endl;
             for(int j=0;j<(int)q.size();j++) {
	     cout << " id= "<< q[j].id() << " q= " << q[j].getConstQ()<<endl;
             }
             for(int j=0;j<(int)constq.size();j++) cout << j
	        << " q1= " << constq[j].first
	        << " q2= " << constq[j].second
	        <<endl;
	       event.list();
		exit(1);
	    }

	    if(optConstQscatt==0) {  // all hadrons can scatter after their formation time
	      v +=  xprod;
	    } else if(optConstQscatt==1) { // leading hadron can scatter before formation time
	      if(isq == 0) v +=  xprod;

            // leading hadron can scatter before formation time, others can scatter elastically.
	    } else if(optConstQscatt==2) { 

	    } else { // all hadron can scatter without formation time.
		tform=0.0;
		pa->setConstQuark(0,1);
		int qq= pa->baryon() == 0 ? 1:2;
		pa->setConstQuark(1,qq);
	    }

	    if(v[0]<xout[0]) {
	      cout << "after string frament t prod is less than t= "<< v[0]
		<< " xout0= "<< xout[0]
		<< " id= "<< id << " r= "<< v;
	      exit(1);
	    }
	    pa->setCoordinate(v);
	    pa->setVertex(v);
	    pa->setFormationTime(tform+v[0]);

	    double meff = sqrt(max(0.0,p.m2Calc()));
	    if(abs(meff-m) > 1e-3) {
	      cout << setprecision(9) << "SoftStrings::stringFragment meff= "<< meff << " m= "<<m
		<< " id= "<< id 
		<<endl;
	      event.list();
	      return false;
	      pa->setPotS(meff-m);
	    }


	    // set decay time.
	    double e=event[i].e();
	    double dect = jamParticleData->lifeTime(pdr,m,e);
	    pa->setLifeTime(v[0]+tform+dect);

            outgoing.push_back(pa);
	    hadrons.push_back(pa);

            if(withMeanField) {
	    if(pa->isMeanField(tform+v[0],optPotential)) {
	      hasMeanField++;
	      //pa->setPotS(sqrt(max(0.0,p.m2Calc())) - m);
	      pa->setPotS(spot);
	      pa->setPotV(vpot);
	      //if(isq==1) had = hadrons.back();
	    } else {
	      etot0 += p[0];
	    }
	    }
	       
	    double gam=e/m;
	    if(isq==1) {
		aveConstFormationTime += tform/gam;
		nCFormationTime++;
	    } else {
		aveFormationTime += tform/gam;
		nFormationTime++;
	    }

	    if(isDebug>1) {
	    cout << event[i].name() << " id= "<< event[i].id() << " isq= "<< isq
	         << " cq1= " << cqq[0]
	         << " cq2= " << cqq[1]
		 <<endl;
	    cout << " tform= " << event[i].tProd()*MM2FM
	        << " tau= " << event[i].tProd()*MM2FM*m/e
	    << " t= " << v[0]
	    << " tform= " << tform <<endl;
	    cout << " x= " << xout << endl;
	    cout << " v= " << v << endl;
	    }
	    
    } // end hadron loop

    if(nconstq != nq ) {
       hadronize->event.list();
	cout << " nconstq= " << nconstq
	    << " nq= " << nq <<endl;
    cout << " isbaryon= "<< isbaryon<<endl;
    cout << " q0= "<< q[0].getConstQ() << " id0= " << q[0].id()
         << " q1= "<< q.back().getConstQ() << " id1= " << q.back().id()
	 <<endl;

	exit(1);
    }


    /*
    Vec4 diff= ptot-ptot0;
    double pdiff = abs(diff.mCalc());
    if(pdiff  > 1e-5) {
	cout << "After string fragmentation pdiff = " << pdiff <<endl;
	event.list();
	cout << " p0= "<< ptot0;
	cout << " p1= "<< ptot;
	exit(1);
        return false;
    }
    */

    /*
  if(had) {
    double e=ptot[0]-etot0;
    double mef2 = e*e-had->pAbs2();
    if(mef2>0.0) {
      had->setPotS(sqrt(mef2)-had->getMass());
    } else {
      cout << "meff<0 "<< mef2<<endl;
    }
  }
  */

  //if(hasMeanField) recoverEnergy(hadrons, ptot, sm, spot, vpot);

  // Option that only the leading hadron has an effective mass after decay.
  // Other particles are put into free.
  //if(potentialHandling==1 && abs(sm-1.0) > 1e-10) {
    
  bool pot=false;
  if(pot && abs(sm-1.0) > 1e-10) {

    double ptotcm = ptot.mCalc();

    Vec4 ptot2=0.0;
    for(int i=0; i<(int)hadrons.size(); i++) {
      hadrons[i]->bstback(ptot);
      ptot2 += hadrons[i]->getP();
      double th=hadrons[i]->getT();
      if(hadrons[i]->isMeanField(th,optPotential)) { 
	//double m = hadrons[i]->getMass();
	//double meff = m*sm;
        hadrons[i]->setPotS(spot);
        hadrons[i]->setPotV(vpot);
      }
    }

    // Find scale factor by Newton method.
    double f;
    double a = 0.0;
    int ntry = 0;
    do {
      f = -ptotcm;
      double df = 0.0;
      for(int i=0; i<(int)hadrons.size(); i++) {
        double m = hadrons[i]->getEffectiveMass();
	double e = sqrt(m*m + pow2(1.0+a)*(hadrons[i]->pAbs2()));
	f += e;
	df += (1.0+a)/e;
      }
      a -= f / df;
      if(++ntry >100) break;
    } while(abs(f) > 1e-8);

    Vec4 ptot4=0.0;
    for(int i=0; i<(int)hadrons.size(); i++) {
      hadrons[i]->multP((1.0+a));
      hadrons[i]->setOnShell();
      hadrons[i]->bst(ptot);
      ptot4 += hadrons[i]->getP();
    }
  }

  hadrons.clear();

  return true;
}

  
void SoftStrings::makeString(int nval, vector<ResolvedPartonQ>& q,vector<pair<int,int>>& constq, Event& event)
{
  // see main21.cc how to do the color configuration.
  // There must be a junction.
  // Need a colour singlet mother parton to define junction origin.
  if(nval==3) {
        event.append( 1000022, -21, 0, 0, 2, 4, 0, 0,
                  0., 0., 0.0* 20, 1.01 * 20);
        event.append(q[0].id(),23,1,0,0,0,
		q[0].col(), q[0].acol(), q[0].p(), q[0].m() );
     constq.push_back(make_pair(0,0));
     constq.push_back(make_pair(q[0].getConstQ(),0));
  } else {
    // set quark
    event.append(q[0].id(),23,q[0].col(), q[0].acol(), q[0].p(), q[0].m() );
    if(q[0].id() > 0) 
      constq.push_back(make_pair(q[0].getConstQ(),0));
    else  // anti-qq
      constq.push_back(make_pair(0,q[0].getConstQ()));
  }

  int col1=0,col2=0;
  int col=q[0].col();
  int acol=q[0].acol();
  if(col !=0 && acol == 0) {
    col1=1;
  } else if(col == 0 && acol != 0) {
    col2=1;
  }

  // set gluons between q and qbar (qq)
  int n=q.size()-1;
  if(nval==3) n--;
  col=q[n].col();
  acol=q[n].acol();
  for(int i=1;i<n;i++) {
    col=101+i-1+col1;
    acol=101+i-1+col2;
    event.append(q[i].id(), 23,col, acol, q[i].p(), q[i].m() );
    constq.push_back(make_pair(0,0));
  }

  if(nval==2) {
    // set qbar (qq)
    col=event.back().col();
    acol=event.back().acol();
    if(col !=0 && acol !=0)  {
	if(col > acol) {
	event.append(q[n].id(), 23, 0,col, q[n].p(), q[n].m() );
        constq.push_back(make_pair(0,q[n].getConstQ()));
	} else {
	event.append(q[n].id(), 23, acol,0, q[n].p(), q[n].m() );
        constq.push_back(make_pair(q[n].getConstQ(),0));
	}

    } else if(col !=0 && acol==0) {
	event.append(q[n].id(), 23,0, col, q[n].p(), q[n].m() );
        constq.push_back(make_pair(0,q[n].getConstQ()));
    } else if(acol !=0 && col==0) {
	event.append(q[n].id(), 23,acol, 0, q[n].p(), q[n].m() );
        constq.push_back(make_pair(q[n].getConstQ(),0));
    } else {
	cout << " error color ?" <<endl;
	cout << " q= " << q[n].col()
	    << " acol= " << q[n].acol()
	    <<endl;
    }

  } else {
    int m=q.size()-2;
    for(int i=m;i< m+2;i++) {
      event.append(q[i].id(),23,1,0,0,0,q[i].col(),q[i].acol(),q[i].p(), q[i].m());
        if(q[i].id()>0) constq.push_back(make_pair(q[i].getConstQ(),0));
	else  constq.push_back(make_pair(0,q[i].getConstQ()));
    }
  }

}

void SoftStrings::recoverEnergy(vector<EventParticle*>& hadrons, Vec4& ptot,double sm, double spot, Vec4& vpot)
{
  double etotcm = ptot.mCalc();

  Vec4 ptot1=0.0;
  double etot0=0.0;
  for(auto& h : hadrons) {
    /*
    double tf=h->getTf();
    if(h->isMeanField(tf,optPotential)) { 
	//double m = h.getMass();
	//double meff = m*sm;
        //h.setPotS(meff - m);
        h->setPotS(spot);
        h->setPotV(vpot);
    }
    */
    ptot1 += h->getP();
    //double mef= h->getEffectiveMass();

    //h->bstback(ptot,mef);
    //h->bstback(ptot);
    //etot0 += h->getPe();
  }

  for(auto& h : hadrons) h->bstback(ptot1);
  Vec4 ptot2=0.0;
  for(auto& h : hadrons) ptot2 += h->getP();
  double etot = ptot1.mCalc();

  cout << " ptot0= "<<ptot;
  cout << " ptot1= "<<ptot1;
  cout << " ptot2= "<<ptot2;
  cout << " etot0 = "<< etot0 << " etot2= "<< ptot1.mCalc() <<endl;

  /*
    // Find scale factor by Newton method.
    double f;
    double a = 0.0;
    int ntry = 0;
    do {
      f = -etotcm;
      double df = 0.0;
      for(int i=0; i<(int)hadrons.size(); i++) {
        double m = hadrons[i]->getEffectiveMass();
	double e = sqrt(m*m + pow2(1.0+a)*(hadrons[i]->pAbs2()));
	f += e;
	df += (1.0+a)/e;
      }
      a -= f / df;
      //cout << " a= "<< a << " f= " << scientific << f <<endl;
      if(++ntry >100) break;
    } while(abs(f) > 1e-8);

  Vec4 ptot4=0.0;
  for(auto& h : hadrons) {
    h.multP((1.0+a));
    h.setOnShell();
    h.bst(ptot);
    ptot4 += h.getP();
  }
    */

  cout << setprecision(8)<< scientific << "etotcm= "<<etotcm << " etot= "<< etot << " dif= "<< etotcm-etot<<endl;

  double a = 1.0;
  int ntry=0;
  while(abs(etotcm-etot) > etotcm*1e-2) {
    a *= etotcm/etot;
    etot=0.0;
    for(auto& h : hadrons) {
      //Vec4 p = h->getP()*a;
      h->multP(a);
      h->setOnShell();
      //double mef=h->getEffectiveMass();
      //etot += sqrt(mef*mef+p.pAbs2());
      etot += h->getPe();
    }
      //cout << " etot = "<< etot  <<" diff= "<< (etot-etotcm)/etotcm << " a= "<< a <<endl;
      if(++ntry>20) break;
  }

  Vec4 ptot3=0.0;
  for(auto& h : hadrons) {
    //h->multP(a);
    //h->setOnShell();
    ptot3 += h->getP();
  }
  cout << " etot after = "<< etot<<endl;
  cout << setprecision(8)<< scientific << "etot3 "<< ptot3;

  Vec4 ptot4=0.0;
  for(auto& h : hadrons) {
    //double mef= h->getEffectiveMass();
    //h->bst(ptot,mef);
    h->bst(ptot);
    ptot4 += h->getP();
  }

  cout << " a= "<< a << " ptot = "<< scientific << (ptot-ptot4).pAbs()
	<<endl;
  cout << " ptot a = "<< ptot;
  cout << " ptot b = "<< ptot1;
  cout << " ptot c = "<< ptot4;

  //cin.get();

}

void SoftStrings::findConstQuark(Event& event,int ibar,int nconstq,
      vector<pair<int,int>>& q)
{
  // put valence quark information for all quarks.
  for(int i=0; i<event.size(); i++) {
    if(event[i].status() != -23) continue;
    int next=event[i].daughter1();
    int ntry=0;
    do {
      if(event[next].status()>0) break;
      q[next]=q[i];
      next=event[next].daughter1();
    } while(++ntry<100);
    if(ntry>99) {
	cout << "SoftString::findConstQuark infinite loop?"<<endl;
	exit(1);
    }
  }

  // add valence quark to the produced hadrons.
  int ip=0;
  do {
    if(!event[ip].isFinal()) continue;
    if(event[ip].isLepton()) continue;
    int mo1 = event[ip].mother1();  // parton1
    int mo2 = event[ip].mother2();  // parton2
    int da1 = event[mo1].daughter1(); // parton1's daughter
    int da2 = event[mo1].daughter2(); // parton2's daughter

    // from ministring into one hadron.
    if(event[da1].status()==81) {
      q[da1].first=q[mo1].first+q[mo2].first;
      q[da1].second=q[mo1].second+q[mo2].second;
      ip=da1;
      continue;
    }

    // two-body final state.
    if(da2-da1==1) {
    if(event[da1].particleDataEntry().isBaryon()) {
      if(event[mo1].isDiquark()) { 
	q[da1]=q[mo1];
	q[da2]=q[mo2];
      } else {
	q[da1]=q[mo2];
	q[da2]=q[mo1];
      }
    } else {
      if(event[mo1].isDiquark()) { 
	q[da1]=q[mo2];
	q[da2]=q[mo1];
      } else {
	q[da1]=q[mo1];
	q[da2]=q[mo2];
      }
    }
    ip=da2;
    continue;
    }

    q[da1]=q[mo1];
    q[da2]=q[mo2];
    if(event[mo1].idAbs() > 1000) { // qq or anti-qq
      if(q[mo1].second==2 && event[da1].particleDataEntry().isBaryon() ==0) {
	q[da1]=make_pair(1,0);
	q[da1+1]=make_pair(1,0);
      }
    }
    if(event[mo2].idAbs() > 1000) {
      if(q[mo2].second==2 && event[da2].particleDataEntry().isBaryon() ==0) {
	q[da2]=make_pair(1,0);
	q[da2-1]=make_pair(1,0);
      }
    }
    ip=da2;
  } while(++ip < event.size());


// debug --------------------------------------------------------------
  int nq=0;
  for(int i=1;i<event.size();i++) {
    if(event[i].isFinal()) {
      nq += q[i].first+q[i].second;
    }
  }
  if(nconstq != nq) {
    //for_each(q.begin(),q.end(),[](int x){cout << x << endl;});
    cout << "stringFrament const quark wrong number! nconstq= "<< nconstq
	<< " n= " << nq << endl;
    //for(int i : q) cout << i << " q= "<< q[i]<<endl;
    for(int i=0;i<(int)q.size();i++) cout << i 
	<< " q1= " << q[i].first
	<< " q2= " << q[i].second
	    <<endl;
     event.list();
     exit(1);
  }
//------------------------------------------------------------------

}

void SoftStrings::findConstQuark2(Event& event,int ibar,int nconstq,
      vector<ResolvedPartonQ>& q, int* iend,int* constq)
{
  const int optDiq1=1;
  const int optDiq2=1;
    
  // Mesonic string.
  if(ibar == 0) {
    if(q[0].getConstQ()==1) {
      constq[0]=1;
      for(int i=1; i<event.size();i++) 
      if(event[i].isFinal()) {iend[0]=i; break;}
    }
    if(q.back().getConstQ()==1) {
      iend[1]=event.size()-1;
      constq[1]=1;
    }

    // check
    if(constq[0]+constq[1] != nconstq) {
	cout << " findConstQ meson nconstq= "<< nconstq
	    << " constq0= "<< constq[0]
	    << " constq1= "<< constq[1]
            << " q = " << q[0].getConstQ()
            << " qbar = " << q.back().getConstQ()
	    << endl;
	event.list();
	exit(1);
    }

    return;
  }

  // baryonic string start.
  int istart=0; // find first point of produced hadron.
  int nhadron=0; // number of hadrons.
  for(int i=1; i<event.size();i++) {
    if(event[i].isFinal()) {
      if(istart==0)  istart = i;
      nhadron++;
    }
  }

  // q-qq -> h_1 + h_2
  if(nhadron==2) {
    iend[0]=istart;
    iend[1]=event.size()-1;
    int mo1 = event[istart].mother1();
    if(event[istart].particleDataEntry().isBaryon()) {
      iend[0]=event.size()-1;
      iend[1]=istart;
      if(event[mo1].isDiquark()) { 
        constq[1]=q[0].getConstQ();
        constq[0]=q.back().getConstQ();
      } else {
        constq[0]=q[0].getConstQ();
        constq[1]=q.back().getConstQ();
      }
    } else {
      if(event[mo1].isDiquark()) { 
        constq[1]=q[0].getConstQ();
        constq[0]=q.back().getConstQ();
      } else {
        constq[0]=q[0].getConstQ();
        constq[1]=q.back().getConstQ();
      }
    }

    if(constq[0]==0 && constq[1]==2) {
      constq[0]=1;
      constq[1]=1;
    }
    if(constq[0]==0) iend[0]=0;
    if(constq[1]==0) iend[1]=0;

    return;
  }


  // check left-hand side.
  if(q[0].getConstQ()==1) {
    constq[0]=1;
    for(int i=1; i<event.size();i++) {
      if(event[i].isFinal()) {iend[0]=i; break;}
    }
  } else if(q[0].getConstQ()==2) {
    for(int i=1; i<event.size();i++) {
      if(event[i].isFinal()) {iend[0]=i; break;}
    }
    if(event[iend[0]].particleDataEntry().isBaryon()) {
      constq[0]=2;
      if(optDiq1==1 && q.back().getConstQ()==0) {
        constq[0]=1;
        constq[1]=1;
        iend[1]=event.size()-1;
      }
    } else if(event[iend[0]+1].particleDataEntry().isBaryon()) {
      if(optDiq2==0) {
        iend[0]++;
	constq[0]=2;
      } else {
        constq[0]=1;
	iend[2]=iend[0]+1;
	constq[2]=1;
      }
    } else {
      constq[0]=1;
      constq[2]=1;
      iend[2]=iend[0]+1;
      cout << "1 stringFragment::findConstQuark strange constq1= "
	  << q[0].getConstQ()
	  <<endl;
      event.list();
      //exit(1);
    }
  }

  // check right-hand side.
  if(q.back().getConstQ()==2) {
    iend[1]=event.size()-1;
    if(event[iend[1]].particleDataEntry().isBaryon()) {
      constq[1]=2;
      if(optDiq1==1 && q[0].getConstQ()==0) {
        constq[1]=1;
        constq[0]=1;
        for(int i=1; i<event.size();i++) 
        if(event[i].isFinal()) {iend[0]=i; break;}
      }
    } else if(event[iend[1]-1].particleDataEntry().isBaryon()) {
      if(optDiq2==0) {
        iend[1]--;
        constq[1]=2;
      } else {
        constq[1]=1;
        iend[2]=iend[1]-1;
        constq[2]=1;
      }
    } else {
      constq[1]=1;
      constq[2]=1;
      iend[2]=iend[1]-1;
      /*
      cout << "stringFragment strange constq1 "<<endl;
      for(int i=0;i<(int)q.size();i++) {
       cout << i << " id= " << q[i].id() << " cq= " 
	   << q[i].getConstQ()<<endl;
      }
      event.list();
      exit(1);
      */
    }

  } else if(q.back().getConstQ()==1) {
	iend[1]=event.size()-1;
	constq[1]=1;
	if(nconstq-q.back().getConstQ()-q.front().getConstQ() == 1) {
	iend[2]=event.size()-2;
	constq[2]=1;
	}
  }

  if(ibar<0) {
      std::swap(constq[0],constq[1]);
      std::swap(iend[0],iend[1]);
  }

  /*
    if(constq[0]==0 && constq[1]==2) {
      constq[0]=1;
      constq[1]=1;
    }
    if(constq[0]==0) iend[0]=0;
    if(constq[1]==0) iend[1]=0;
    */

  bool ok=true;
  if(nconstq != constq[0]+constq[1]+constq[2]) ok=false;
  //if(iend[0]==iend[1]) ok=false;
  //if(iend[0]==iend[2]) ok=false;
  //if(iend[1]==iend[2]) ok=false;

  if(ok==false) {
    cout << "stringFrament const quark wrong number!" << endl;
    cout << " q0= "<< q[0].getConstQ()
         << " q1= "<< q.back().getConstQ()
	 <<endl;
    cout << " nconstq= " << nconstq
         << " q1= " << constq[0] << " iend= "<< iend[0]
	 << " q2= " << constq[1] << " iend= "<< iend[1]
	 << " q3= " << constq[2] << " iend= "<< iend[2]
	 <<endl;
    for(int i=0;i<(int)q.size();i++) {
      cout << i << " id= " << q[i].id() << " cq= " << q[i].getConstQ()<<endl;
     }
     //event.list();
     //exit(1);
  }


}

//***********************************************************************
//      function xsamp1(xmin,xmax,srt)
//double SoftStrings::xSample1(double srt,double xmin,double xmax)
//...Sampe x of valence quarks for baryon, DPM type.
double SoftStrings::xSample1(double srt,double xmin, double xmax)
{

//...Fritiof type  distribuiton: 1/x
//     if(srt.lt.10.0d0) then
//       xsamp1=xmin*(xmax/xmin)**rn(0)
//       return
//     endif

//...p(x)=(1-x)^d/(x^2+cutt^2)^b

//      b=hipr1(46)
//      d=hipr1(44)
//      cutt=hipr1(45)/srt
    double b=0.25;
    double b1=2*b;
    double d=1.5;
    double cutt=0.1/srt;
    double p,x;
    do {
      double cut1=pow(xmin+cutt,(1.0-b1));
      double cut2=pow(xmax+cutt,(1.0-b1));
      x=pow( (cut2-cut1)*rndm->flat() + cut1 ,1.0/(1.0-b1) ) - cutt;
      p = pow(1.0-x,d)*pow( pow2(x+cutt)/(2.0*(x*x+cutt*cutt)), b);

    } while( p < rndm->flat() );

    return x;

}


} // end namespace jam2

