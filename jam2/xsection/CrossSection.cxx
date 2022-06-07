#include <iostream>
#include <cmath>
#include <jam2/xsection/CrossSection.h>
#include <jam2/xsection/XsecTable.h>
#include <jam2/xsection/SigmaBB1.h>
#include <jam2/xsection/SigmaBB2.h>

using namespace std;

namespace jam2 {

const double eKinMin=0.001; // parc(41)
const double Mn=0.9396; // 0.93957  parc(24)
const double Mpi=0.1396; //  parc(27)

CrossSection::CrossSection(Pythia8::Info* inf,Pythia8::Settings* s,JamParticleData* pd, Pythia8::StringFlav *flavesel,
    Pythia8::Rndm* r)
{
  info= inf;
  settings = s;
  optJAM1=settings->mode("Cascade:optBBJAM1");
  string bwFileName1 = settings->word("Cascade:bwFileName1");
  string bwFileName2 = settings->word("Cascade:bwFileName2");
  jamParticleData = pd;
  rndm = r;
  sampleMass = new SampleMass(settings,bwFileName2,jamParticleData,flavesel,rndm);
  if(optJAM1==1) {
    sigBB = new SigmaBB1(info,settings,jamParticleData,sampleMass,rndm);
  } else {
    sigBB = new SigmaBB2(info,settings,jamParticleData,sampleMass,rndm,bwFileName1);
  }
  sigMB = new SigmaMB(settings,jamParticleData,sampleMass);
  sigKB = new SigmaKB(settings,jamParticleData,sampleMass,rndm);
  sigMM = new SigmaMM(settings,jamParticleData,sampleMass,rndm);
  sigABB = new SigmaABB(settings,jamParticleData,sampleMass);

  ecmStringMBc=settings->parm("Cascade:ecmStringMBc");
  ecmStringBB=settings->parm("Cascade:ecmStringBB");
}

CrossSection::~CrossSection()
{
  delete sigBB;
  delete sigMB;
  delete sigKB;
  delete sigMM;
  delete sigABB;
  delete sampleMass;
}

// currently S=0 BB only
void CrossSection::selectOutGoingParticle(CollisionPair& cpair)
{
  double xsig=cpair.getXsig();
  double sig=cpair.getSigma();
  double sigel=cpair.getSigmaElastic();
   if(xsig > sig-sigel) {
     cout << "Strange xsig= "<< xsig << " sigin= "<< sig-sigel;
     cout << " has anti= "<< cpair.isAnti() <<endl;
   }
   
    int id3=id1;
    int id4=id2;
    int id5=0;
    int id6=0;
    double ma=cpair.M(0);
    double mb=cpair.M(1);
    double m5=0.0, m6=0.0;
    double m3=ma,m4=mb;
    ParticleDataEntry *pa3=pd1;
    ParticleDataEntry *pa4=pd2;
    ParticleDataEntry *pa5=0, *pa6=0;

    int ch = cpair.selectChannel();
    OutGoing out = cpair.getOutGoing();
    if(ch>=0) {
	id3= out.id1;
	id4= out.id2;
	id5= out.id3;
	id6= out.id4;
	if(out.id1 != 92)  {
	    bool ok= sigBB->sampleResonanceMass(eCM,id3,id4,ma,mb);
	    if(!ok) {
                cpair.setChannel(-1);
		return;
	    }
	    if(ok) {
	    if(rndm->flat()< 0.5) {
	      pa3 = sigBB->outParticle(0);
	      pa4 = sigBB->outParticle(1);
	      id3= pa3->id();
	      id4= pa4->id();
	      m3 = ma;
	      m4 = mb;
	    } else {
	      pa4 = sigBB->outParticle(0);
	      pa3 = sigBB->outParticle(1);
	      id4= pa4->id();
	      id3= pa3->id();
	      m4 = ma;
	      m3 = mb;
	    }
	}
	}
    }

    if(cpair.isAnti()) {
	if(pa3 !=0 && pa3->hasAnti()) id3 *= -1;
	if(pa4 !=0 && pa4->hasAnti()) id4 *= -1;
	if(pa5 !=0 && pa5->hasAnti()) id5 *= -1;
	if(pa6 !=0 && pa6->hasAnti()) id6 *= -1;
    }

    cpair.setOutGoingID(id3,id4,id5,id6);
    cpair.setOutGoingMass(m3,m4,m5,m6);
    //cpair.setOutGoingParticle(pa3,pa4,pa5,pa6);

}

//...Purpose: to compute total and elastic cross sections.
double CrossSection::sigma(CollisionPair& cpair)
{
//----------------------------------------------------------------------*
//...(Outputs)
//...sig    : total cross section (mb)
//...sigel  : elastic cross section (mb)
//----------------------------------------------------------------------*
    id1 = cpair.getID(0);
    id2 = cpair.getID(1);

    //pd1 = cpair.getParticleDataEntry(0);
    //pd2 = cpair.getParticleDataEntry(1);
    pd1=jamParticleData->find(id1);
    pd2=jamParticleData->find(id2);

    eCM=cpair.getCMenergy();
    //pCM=cpair.getCMmomentum();

    int icltype=cpair.getCollType();

    anti=false;

    sigTot=0.0;
    sigEl=0.0;

    //mste(3)=0  ! final particle

//===================================
//...Baryon + baryon or antiB + antiB collision
//===================================
    if(icltype == 1 || icltype == 5) {

	/*
	if(cpair.mode()==3 && cpair.outgoingSize()>0) {
	    selectOutGoingParticle(cpair);
	    return 0.0;
	}
	*/

	if(eCM < 2*Mn+eKinMin) return 0.0;

	SigmaBaryonBaryon(cpair);

	/*
	if(cpair.mode()==3) {
	    if(cpair.outgoingSize()==0) {
		cout << " CrossSection::BB collision channel=0?" 
		     << " id1= "<< id1 << " id2= "<< id2
		     << " eCM= "<< eCM <<endl;
		exit(1);
	    }
	    selectOutGoingParticle(cpair);

      if(anti) cpair.setAnti();

	    return 0.0;
	}
	    */
	   

//===================================
//...Baryon + meson collisions part
//===================================
    } else if(icltype == 2) {

	if(eCM < Mn+Mpi+eKinMin) return 0.0;

	SigmaMesonBaryon(cpair);


//=============================
//...Meson + meson collisions
//=============================
    } else if(icltype == 3) {

	if(eCM < 2*Mpi+eKinMin) return 0.0;

	SigmaMesonMeson(cpair);

//=============================
//...Anti-baryon + baryon  collisions.
//=============================
    } else if(icltype == 4) {

	if(eCM < 2*Mn+eKinMin) return 0.0;

	SigmaAntiBaryonBaryon(cpair);

    } else {
	cout << " CrossSection::sigma wrong icltype= " << icltype
	     << " id1= "<<  cpair.getID(0)
	     << " id2= "<<  cpair.getID(1)
	    <<endl;
	exit(1);
    }

    if(anti) cpair.setAnti();
    return sigTot;

}

double CrossSection::SigmaBaryonBaryon(CollisionPair& cpair)
{
  id1 = cpair.getID(0);
  id2 = cpair.getID(1);
  int nHeavy=getHeavyQB(id1)+getHeavyQB(id2);

//...In the case of antibaryon-antibaryon collision,
  if(id1 < 0 && id2 < 0) {  // antib-antib
    anti=true;
    id1=-id1;
    id2=-id2;
    cpair.setAnti();
  }
  int str=cpair.getTotalStrangeness();

  if(str ==0  && nHeavy == 0) {
    sigBB->calc(cpair);
    if(cpair.mode()==3) {
      selectOutGoingParticle(cpair);
      if(anti) cpair.setAnti();
      return 0.0;
    }

    //....Lambda-N/Sigma-N/Xi-N/LL...
  } else if(nHeavy == 0) {

    sigBB->calcS(cpair);
    if(cpair.mode()==3) {
      double sigin=cpair.getSigma()-cpair.getSigmaElastic();
      if(cpair.getCMenergy() > ecmStringBB) {
        cpair.setOutGoing(sigin,92,92);
      } else {
	sigmaTChannel(cpair,sigin);
      }
      return 0.0;
    }

  } else {
    cpair.setXSTotal(20.0);
    cpair.setXSElastic(20.0);
  }

  sigTot=cpair.getSigma();
  sigEl=cpair.getSigmaElastic();
  return sigTot;

}

double CrossSection::SigmaMesonBaryon(CollisionPair& cpair)
{
  //pd1 = cpair.getParticleDataEntry(0);
  //pd2 = cpair.getParticleDataEntry(1);
  int ibar1=cpair.baryon(0);
  int ibar2=cpair.baryon(1);

  bool iswap=false;
  if(ibar1 !=0) {
    cpair.swap();
    iswap=true;
  }

  //...Check if this is a anti-B + M collision.
  if(ibar1 == -3 || ibar2 == -3) {
    anti=true;
    cpair.setAnti();
  }

  id1 = cpair.getID(0);
  id2 = cpair.getID(1);
  int str=cpair.getTotalStrangeness();
  int strb=cpair.getStrange(1);

//...Light collision systems.
    int nHeavy=getHeavyQM(id1)+getHeavyQB(id2);
    if(nHeavy != 0) {
	cpair.setXSTotal(20.0);
	cpair.setXSElastic(20.0);
	//return SigmaMBHeavy();

    //...piN/rhoN lambdaK+ etc.
    } else if(str == 0) {
	sigMB->calcS0(cpair);
            //call jamcmbs0(msel,srt,pr,istr,kfb,kfm,kcb,kcm,izb,izm

    //....K-p,Lambda+pi etc
    } else if(str == -1) {
	sigMB->calcS1(cpair);
            //call jamcmbs1(msel,srt,pr,istr,kfb,kfm,kcb,kcm,izb,izm

    //...Xi-pi etc.
    } else if(str == -2) {
	sigMB->calcS2(cpair);
            //call jamcmbs2(msel,srt,pr,istr,kfb,kfm,kcb,kcm,izb,izm

    //....Kaon induced collisions.
    } else if(str == 1 && strb == 0) {

	sigKB->calcKN(cpair);
            //call jamckaon(msel,kfm,kfb,srt,pr,ijet

    } else {
    XsecTable::additiveQuarkModel(id1,id2,sigTot,sigEl);
    cpair.setXSTotal(sigTot);
    cpair.setXSElastic(sigEl);
    double srt=cpair.getCMenergy();
    if(srt != eCM) {
      cout << " srt= "<< srt << " ecm= "<< eCM<<endl;
      exit(1);
    }
    if(srt < 1.8) cpair.setXSElastic(sigTot);
    if(cpair.mode()==3) {
      double sigin=sigTot-sigEl;
      // t-channel string formation.
      if(srt > ecmStringMBc) {
        cpair.setOutGoing(sigin,92,92);
      } else {
      // t-channel resonance formation.
        sigmaTChannel(cpair,sigin);
      }
    }
  }


    if(cpair.isAnti()) {
	 cpair.setAnti();
    }
    if(iswap) {
	cpair.swap();
	cpair.swapOutgoing();
    }

    sigTot=cpair.getSigma();
    sigEl=cpair.getSigmaElastic();
    return sigTot;
}


double CrossSection::SigmaAntiBaryonBaryon(CollisionPair& cpair)
{
    sigABB->calc(cpair);
    sigTot=cpair.getSigma();
    sigEl=cpair.getSigmaElastic();
    return sigTot;
}

double CrossSection::SigmaMesonMeson(CollisionPair& cpair)
{
    sigMM->calc(cpair);
    sigTot=cpair.getSigma();
    sigEl=cpair.getSigmaElastic();
    return sigTot;
}

/*
double CrossSection::SigmaMBHeavy()
{
    XsecTable::additiveQuarkModel(id1,id2,sigTot,sigEl);
    //cpair.setXSTotal(sigTot);
    //cpair.setXSElastic(sigEl);
    sigEl=sigTot;
    return sigTot;
}
*/

void CrossSection::sigmaTChannel(CollisionPair& cpair,double sigin)
{
  // test: elastic
  //cpair.setChannel(-1);
  //return;

  // t-channel resonance formation.
  double m1=cpair.M(0), m2=cpair.M(1);
  //pd1 = cpair.getParticleDataEntry(0);
  //pd2 = cpair.getParticleDataEntry(1);
  int kfm=cpair.getID(0);
  int kfb=cpair.getID(1);
  pd1=jamParticleData->find(kfm);
  pd2=jamParticleData->find(kfb);
  int idm=cpair.getPID(0);
  int idb=cpair.getPID(1);
  int izm=cpair.getZ(0)/3;
  int izb=cpair.getZ(1)/3;
  //double pr=cpair.getCMmomentum();
  bool preh1= cpair.qFactor(0) < 1.0 ? true: false;
  bool preh2 = cpair.qFactor(1) < 1.0 ? true: false;
  int icon=sampleMass->jamrmas2(cpair.getCollType(),pd1,pd2,kfm,kfb,
		//idm,idb,izm,izb, m1,m2,preh1,preh2,srt,cpair.isAnti());
		idm,idb,izm,izb, m1,m2,preh1,preh2,eCM,cpair.isAnti());

  if(icon==1) {
    cpair.setChannel(-1);
    return;
  }
  cpair.setOutGoing(sigin,sampleMass->id(0),sampleMass->id(1));
  cpair.setOutGoingMass(sampleMass->m(0),sampleMass->m(1));

}

} // namespace jam2
