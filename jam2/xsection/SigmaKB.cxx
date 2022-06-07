#include <cstdlib>
#include <jam2/xsection/SigmaKB.h>
#include <jam2/xsection/XsecTable.h>
#include <jam2/hadrons/JamStdlib.h>

using namespace std;

namespace jam2 {

using Pythia8::ParticleDataEntry;

static const double emkc=.49360, emk0=.49770,emk=.495650,widk=0.05;
static const double empion=0.139,emdelt=1.232,widdlt=0.12;
static const double emp=0.93830 ,emn=0.93960,emnuc=0.93895;
static const double ekinmi=0.001,emdmin=1.082;
static const double bhad=2.3,eps=0.0808,facel=0.0511;

// Data for one delta production.
int kfdelt[8][2][2]={
      {{321,2214},    {311,2224}}, // 1. k+ p->k+ d+   k+ p->k0 d++
      {{311,2114}, {321,1114}},    // 2. k0 n->k0 d0    k0 n->k+ d-
      {{321,2114}, {311,2214}},    // 3. k+ n->k+ d0   k+ n->k0 d+
      {{311,2214}, {321,2114}},   // 4. k0 p->k0 d+  k0 p->k+ d0
      {{321,2224}, {0,0}},        // 5. k+ d++ -> k+ d++
      {{311,1114}, {0,0}},        // 6. k0 d- -> k0 d- 
      {{321,1114}, {311,2114}},  // 7. k+ d- -> k+ d-   k+ d- -> k0 d0
      {{311,2224}, {321,2214}}   // 8. k0 d++ -> k0 d++   k0 d+ -> k+ d+
      };

SigmaKB::SigmaKB(Settings* s,JamParticleData* table, SampleMass* sm,
	    Pythia8::Rndm* r) : SigmaMB(s,table,sm)
{ 
  rndm=r;

  // make map for the reaction in which charge of delta changes.
  ParticleTable* tablem = jamtable->getDmstar();
  ParticleTable* table0 = jamtable->getD0star();
  ParticleTable* table1 = jamtable->getDpstar();
  ParticleTable* table2 = jamtable->getDppstar();
  for(int i=0;i<(int)tablem->size();i++) {
    int kfd= tablem->getParticle(i)->id();
    int kf0= table0->getParticle(i)->id();
    int kf1= table1->getParticle(i)->id();
    int kf2= table2->getParticle(i)->id();
    deltaMto0[kfd]=kf0;  // map from Delta-  to Delta0.
    delta2toP[kf2]=kf1;  // map from Delta++ to Delta+.
  }

}

//***********************************************************************
//      subroutine jamckaon(kfm,kfb,srt,pr,ijet
//     $        ,em1,em2,kf1,kf2,sig,sigel,sigin,mchanel,mxchan,ich,icon)

void SigmaKB::calcKN(CollisionPair& cpair0)
{
//...Purpose: to treat Kaon-nonstrange baryon collisions.

//...mSel=1: calculate total and elastic cross sections.
//...mSel=2: calculate inelastic cross sections.
//...mSel=3: Monte Carlo evaluation of inelastic channel.

  cpair = &cpair0;
  mChanel=0;
  sigTot=0.0;
  sigEl=0.0;
  cpair->setXS(0.0, 0.0, 0.0);

  kfm=cpair->getID(0);
  kfb=cpair->getID(1);

  //...Check particle code.
  if(kfm <= 0 ||  kfb <= 0) {
    cout << "(jamckaon:)invalid particle kfm= " << kfm 
	 << " kfb= " << kfb << endl;
    exit(1);
  }

    izm=cpair->getZ(0)/3;
    izb=cpair->getZ(1)/3;
    izt=izm+izb;
    srt=cpair->getCMenergy();
    pr=cpair->getCMmomentum();
    istr=cpair->getTotalStrangeness();
    mSel=cpair->mode();

    int kfma=abs(kfm);
    //kflr1=(kfma/10000)%10;
    kflr1=(kfma/10000)%100;
    kflb1=(kfma/100)%10;
    kflc1=(kfma/10)%10;
    kfls1=kfma%10;
    kfm1=10*kflb1+kflc1;

    plab=0.0;
    snew=srt;
    if((kfm == 311 || kfm == 321) && (kfb == 2112 || kfb == 2212)) {
	double ema,emb;
         if(kfm == 311)  ema=emk0;
         if(kfm == 321)  ema=emkc;
         if(kfb == 2112) emb=emn;
         if(kfb == 2212) emb=emp;
         snew=srt;
         if(srt > ema+emb) plab=XsecTable::plabsr(srt,ema,emb);
    } else {
        if(kfm1 == 32) {
	    snew=sqrt(emnuc*emnuc+pr*pr)+sqrt(emkc*emkc+pr*pr);
	    if(snew > emkc+emnuc) 
		plab=XsecTable::plabsr(snew,emkc,emnuc);
	    else
		return;

	} else if(kfm1 == 31) {
	    snew=sqrt(emnuc*emnuc+pr*pr)+sqrt(emk0*emk0+pr*pr);
	    if(snew > emk0+emnuc)
		plab=XsecTable::plabsr(snew,emk0,emnuc);
	    else
		return;
	} else {
          cout << "(jamckaon:3) no path kfm= " << kfm 
	      << " kfm1 = " << kfm1 
	      << endl;
	  exit(1);
	}
    }

//------------------------------------------------
//.....Get total and elastic cross sections.
//------------------------------------------------

    if(mSel == 1) {

	//......K+ p/K0n i.e. isospin =1 channel.
	if(izm == izb) {

	    // K+ p Total/elastic
            if (plab <= 3.0) {
              sigTot=XsecTable::jamsighh(15,snew);
	    //} else if(snew < 30.0) {
            //  sigTot=XsecTable::jamchc96(8,plab);
	    } else {
              //sigTot=XsecTable::jamrgg96(snew,10);
              sigTot=XsecTable::pdg2016(snew,-4);
	    }

	    // K+ p elastic
            if (plab <= 2.0) {
              sigEl=XsecTable::jamsighh(16,snew);
              //sigTot=XsecTable::jamrgg96(snew,10);
	    } else  if(plab <= 100.0)  {
              sigEl=1.06*XsecTable::jamchc96(9,plab);
	     
	    } else {
              double bel=2*(2.3+0.8)+4.0*pow(srt*srt,0.079) - 4.2;
              sigEl=1.13*facel*sigTot*sigTot/bel;
	    }
	//.......K+n/K0p
	} else {
 
	    //... K+ n Total/elastic
            //if(plab <= 4.0) {
            if(plab <= 5.0) {

              //if(plab <= 1.05)
              //  sigTot=18.6*pow(plab,0.86);
              //else
                sigTot=XsecTable::jamsighh(17,snew);

              //if(plab <= 1.05)
              //  sigEl=7.6761*pow(plab,0.566017);
              //else if(plab <= 4.5)
              //  sigEl=4.03297+9.71556/(1.10417+exp((plab-0.674619)/0.835265));
              //else
                sigEl=XsecTable::jamsighh(18,snew);

	    } else if(snew < 15.0) {
              sigTot=XsecTable::jamchc96(10,plab);
              sigEl=1.03*XsecTable::jamchc96(9,plab);
	    } else {
              sigTot=XsecTable::pdg2016(snew,-5);
              //sigTot=XsecTable::jamrgg96(snew,11);
              double bel=2*(2.3+0.8)+4.0*pow(srt*srt,0.079)-4.2;
              sigEl=1.08*facel*sigTot*sigTot/bel;
	    }

	}
	cpair->setXS(sigTot,sigEl,sigTot-sigEl);
	return;
    }

    sigRand=cpair->getXsig();

//------------------------------------------------
//.....Get charge exchange cross sections.
//------------------------------------------------
  if(izt==1 && srt < emn+emkc+eKinMin) {
    chargeExchange();
    if(sigRand <=0.0) return;
  }

//------------------------------------------------
//.....Get inelastic cross sections.
//------------------------------------------------
//...(1) K^+ + p --> K^+  + D+   1/4*sig(t=1)
//...(2) K^+ + p --> K^0  + D++  3/4*sig(t=1)
//...(3) K^+ + p --> K^+* + p (k+* p --> k+ p k* absorption)

//...1. K^+ + n --> K^0 + p                          charge exchange
//...2. K^+ + n --> K^+ + D0  0.25*sig(t=1) isig=1   delta production
//...3. K^+ + n --> K^0 + D+  0.25*sig(t=1) isig=1   delta production
//...4. K^+ + n --> K^+* + n  isig=3                 k* production
//...5. K^+ + n --> K^0* + p  isig=4                 k* production

//...1. K^0 + p --> K^+ + n    charge exchange
//...2. K^0 + p --> K^0 + D+   0.25*sig(t=1)
//...3. K^0 + p --> K^+ + D0   0.25*sig(t=1)
//...4. K^0 + p --> K^0* + p   3
//...5. K^0 + p --> K^+* + n   4

//...1. K^0 + n --> K^0 + D0   0.25
//...2. K^0 + n --> K^+ + D-   0.75
//...3. K^0 + n --> K^0* + n   t=1  equated to k+p->k+* p
//     (K^0* + n --> K^0 + n   k* absorption )

  ibra=0;
  if(izm == 1 && izb == 1)  ibra=1;  // k+ p
  if(izm == 0 && izb == 0)  ibra=2;  // k0 n
  if(izm == 1 && izb == 0)  ibra=3;  // k+ n
  if(izm == 0 && izb == 1)  ibra=4;  // k0 p

  if(izm == 1 && izb == 2)  ibra=5;  // k+ d++
  if(izm == 0 && izb == -1) ibra=6;  // k0 d-
  if(izm == 1 && izb == -1) ibra=7;  // k+ d-
  if(izm == 0 && izb == 2)  ibra=8;  // k0 d++

// Delta production. Check if there is enough energy to produce delta
  double em1 = cpair->M(0);

  int id2=cpair->getPID(1);
  if(id2 == id_delt || id2 == id_delts) {
    if(izt >= 0 && izt <= 2) {
      deltaAbsorption();
      //cout << " after deltaAbsorption sigRand= " << sigRand<<endl;
      if(sigRand <=0.0) return;
    }
  }

  if(srt > em1+emdmin+empion+ekinmi) {
    sigKNDelta();
    //cout << " after sigKNDelta sigRand= " << sigRand<<endl;
    if(sigRand <=0.0) return;
  }

//...Kstar production. Check if there is enough energy to produce delta
  double em2 = cpair->M(1);
  if(srt > (em2+empion+emk+empion+ekinmi)) {
    sigKNKstar();
    //cout << " after sigKNKstar sigRand= " << sigRand<<endl;
    if(mSel ==3 && sigRand <=0.0) return;
  }

  //  if(kfm.eq.311.or.kfm.eq.321) goto 200
  // (5) K*(892) absorption using detailed balance.
  if(kfm != 311 && kfm != 321) {
    absorbKstar();
    //cout << " after absorbKstar sigRand= " << sigRand<<endl;
    if(sigRand <=0.0) return;
  }
 
  //...(6) K(892)+D(1232) production.
  //  if(id2.eq.id_delt.or.id2.eq.id_delts) goto 300
  if(srt > emk+2*empion+emnuc+ekinmi) {
    KstarD();
    //cout << " after KstarD sigRand= " << sigRand<<endl;
    if(sigRand <=0.0) return;
  }

  bool isKD=true;
  if(kfm == 311 || kfm == 321) isKD=false;
  if(id2 != id_delt) isKD=false;
  if(izt < 0 || izt > 2)  isKD=false;

  // (7) K*(892)D(1232) absorption using detailed balance.
  if(isKD) {
    KstarDabsorb();
    //cout << " after KstarDabsorb sigRand= " << sigRand<<endl;
    if(sigRand <=0.0) return;
  }

  if(mSel !=3) return;

  // Monte-Carlo for t-channel contribution.
  double sigtr=cpair->getSigInT()*4.0/srt;
  sigRand -= sigtr;

    // t-channel string formation.
    //if(srt > 4.0 && sigRand < 0.0) {
    if(srt > 4.0) {
      cpair->setOutGoing(sigtr,92,92);
      return;
    }

    if(srt < 1.7) {
        cpair->setChannel(-1);
        return;
    }

    // t-channel resonance formation.
    double m1=cpair->M(0), m2=cpair->M(1);
    //pd1 = cpair->getParticleDataEntry(0);
    //pd2 = cpair->getParticleDataEntry(1);
    pd1=jamtable->find(kfm);
    pd2=jamtable->find(kfb);
    int idm=cpair->getPID(0);
    int idb=cpair->getPID(1);
  bool preHadronA = cpair->qFactor(0) < 1.0 ? true: false;
  bool preHadronB = cpair->qFactor(1) < 1.0 ? true: false;
    int icon=sampleMass->jamrmas2(cpair->getCollType(),pd1,pd2,kfm,kfb,
		idm,idb,izm,izb, m1,m2,preHadronA,preHadronB,srt,cpair->isAnti());

    if(icon==1) {
        cpair->setChannel(-1);
        return;
    }
    cpair->setOutGoing(sigtr,sampleMass->id(0),sampleMass->id(1));
    cpair->setOutGoingMass(sampleMass->m(0),sampleMass->m(1));

}

bool SigmaKB::chargeExchange()
{
  mChanel=1;
  double snew1=sqrt(emn*emn+pr*pr)+sqrt(emkc*emkc+pr*pr);
  double plab1=0.0;
  if(snew1 > emkc+emn) {
    plab1=sqrt( pow2((snew1*snew1-emkc*emkc-emn*emn)/(2.0*emp))-emkc*emkc);
  } else
    return false;

  // Why do you only compute K+ n ?   K0 p?
  if(plab1 <= 2.0) {
    sigin[0]=XsecTable::jamsighh(20,snew1); // K+ n charge exchange.
  } else {
    double pl=min(6.3,plab1);
    sigin[0]=XsecTable::jamchc96(26,pl);
  }

//...Monte Calro for charge exchange.
  if(mSel != 3) return true;

  sigRand -= sigin[0];
  if(sigRand > 0.0) return true;

  double em1 = cpair->M(0);
  double em2 = cpair->M(1);
  int kf1=0,kf2=0;
    if(kfm == 321) { 
      kf1=311;
      em1=emk0;
    } else if(kfm == 311) {
      kf1=321;
      em1=emkc;
    } else if(kfm1 == 31) {
      kf1=kflr1*10000+320+kfls1;
    } else if(kfm1 == 32) {
      kf1=kflr1*10000+310+kfls1;
    }

  if(kfb == 2112) {
    kf2=2212;
    em2=emp;
  } else if(kfb == 2212) {
    kf2=2112;
    em2=emn;
  } else {

    int kf=abs(kfb);
    int kflr2=(kf/10000)%100;
    int kfls2=kf%10;
    if(izb == 1) {
      kf2=10000*kflr2 + 2210 + kfls2;
    } else if(izb == 0) {
      kf2=10000*kflr2 + 2110 + kfls2;
    }
  }
    if(cpair->isAnti()) {
	kf1 *= -1;
	kf2 *= -1;
    }
    cpair->setOutGoing(sigin[0],kf1,kf2);
    cpair->setOutGoingMass(em1,em2);
    return true;
}

// (2) D(1232) production.
void SigmaKB::sigKNDelta()
{
  double sigt1=jamsigkn(1,snew,plab);
  sigin[mChanel++]=sigt1;
  //cout << " sign= " << sigin[mChanel-1] << " mch= " << mChanel<<endl;

  //if(mSel !=3) return true;
  //if(sigRand <=0.0) return true;

// Monte Carlo evaluation of delta channel.

  double sigin1, sigin2;
  if(ibra==1 || ibra==2) {
    sigin1=sigt1*0.25;  // K+ p -> K+ D+ 1/4(t=1)
    sigin2=sigt1*0.75;  // K+ p -> K0 D0 3/4(t=1)
  } else if(ibra==3 || ibra==4 || ibra==7 || ibra==8) {
    sigin1=sigt1*0.5;
    sigin2=sigt1*0.5;
  } else {               // 5: K+ D++ 6 K0 D-
    sigin1=sigt1;
    sigin2=0.0;
  }

  int iselect=-1;
  sigRand -= sigin1;
  if(sigRand<=0.0) {
    iselect=0;
  } else {
    sigRand -= sigin2;
    if(sigRand<=0.0) iselect=1;
  }

  if(iselect == -1) return;

  //cout << "ibra= " << ibra << " iselect= " << iselect << " ran= " << sigRand<<endl;

    int kf1=kfdelt[ibra-1][iselect][0];
    int kf2=kfdelt[ibra-1][iselect][1];

  //cout << " kf1= " << kf1 << " kf2= "<< kf2 <<endl;

    pd1=jamtable->find(kf1);
    pd2=jamtable->find(kf2);
    double em1=pd1->m0();

    //emmin=max(pmas(kc2,1)-pmas(kc2,3),emdmin+parc(41))
    double emmin=emdmin+ekinmi;
    double emmax=srt-em1-ekinmi;
    double emr=pd2->m0();
    double wid=pd2->mWidth();
    //double em2 = decay->BWMass(emmin,emmax,emr,wid);
    double em2 = jamtable->BWMass(emmin,emmax,emr,wid);
    if(cpair->isAnti()) { kf1 *= -1; kf2 *= -1; }
    cpair->setOutGoing(sigt1,kf1,kf2);
    cpair->setOutGoingMass(em1,em2);

}

void SigmaKB::deltaAbsorption()
{
//...(3) D(1232) absorption using detailed balance.

//  if(id2.eq.id_delt.and.(izt.ge.0.and.izt.le.2)) then

  // xsection of D* K -> N K is assumed to be the same as D K -> N K.

  double em1 = cpair->M(0);
  if(srt <= em1+emnuc) return;
  pd2=jamtable->find(kfb);
  double emd=pd2->m0();
  double widd=pd2->mWidth();
  double widcof=max(1.0,M_PI/(atan(2.0*(srt-emk-emd)/widd)-
                    atan(2.0*(emnuc+empion-emd)/widd)));
  double pfinal=XsecTable::pawt(srt,em1,emnuc);
  double facsp=2.0/double(kfb%10);
  sigin[mChanel++]=facsp*pow2(pfinal/pr)*jamsigkn(1,srt,plab)*widcof;
  if(mSel != 3) return;

  // Monte Carlo evaluation of delta absorption channel.
  double sigab= 0.25*sigin[mChanel-1];
  if(izb == -1 || izb == 2) sigab = 0.75*sigin[mChanel-1];
  sigRand -= sigab;
  if(sigRand> 0.0)return;

  int kf1=0,kf2=0;
  switch (izb) {
    case -1:  // k+ d- -> k0 + n 3/4
        kf2=2112;
        kf1=10000*kflr1+310+kfls1;
        break;

    case 0: // d0

        if(izm==1) {
          if(rndm->flat() < 0.5) {  // k+ d0 -> k+ n 1/4
            kf2=2112;
            kf1=10000*kflr1+320+kfls1;
	  } else {  // k+ d0 -> k0 + p 1/4
            kf2=2212;
            kf1=10000*kflr1+310+kfls1;
	  }
	} else if(izm == 0) {  // d0 k0 -> k0 + n0 1/4
          kf2=2112;
          kf1=kfm;
	}
	break;

    case 1: // d+

      if(izm==0) {
        if(rndm->flat() <= 0.5) {  // d+ k0 -> n k+  1/4
          kf2=2112;
          kf1=10000*kflr1+320+kfls1;
	} else {  // d+ k0 -> p k0
          kf2=2212;
          kf1=10000*kflr1+310+kfls1;
	}
      } else if(izm==1) {  // d+ k+ -> p + k+ 1/4
        kf2=2212;
        kf1=kfm;
      }
      break;

    case 2:  // d++ k0 -> k+ p  3/4
      kf2=2212;
      kf1=10000*kflr1+320+kfls1;
      break;
  }

  if(kfm ==  311 || kfm == 321) {
      pd1=jamtable->find(kf1);
      em1=pd1->m0();
  }
  pd2=jamtable->find(kf2);
  double em2=pd2->m0();
  if(cpair->isAnti()) { kf1 *= -1; kf2 *= -1; }
  cpair->setOutGoing(sigin[mChanel-1],kf1,kf2);
  cpair->setOutGoingMass(em1,em2);

}

// (4) K(892) production.
void SigmaKB::sigKNKstar()
{
  double sigin1=0.0,sigin2=0.0;
  int jch=1;
  // K+p, K0 n       K+D++      K0 D-
  if(ibra <= 2 || ibra == 5 || ibra== 6) {
    sigin1=jamsigkn(2,snew,plab); //...2:K+ p => K+*(890) p
    sigin[mChanel]=sigin1;
  } else {
    jch=2;
    sigin1=jamsigkn(3,snew,plab); //...3:K+ n => K+*(890) n
    sigin[mChanel]=sigin1;
    sigin2=jamsigkn(4,snew,plab); //...4:K+ n => K0*(890) p
    sigin[mChanel++]=sigin2;
  }
  mChanel++;

  if(mSel !=3) return;

//....Monte Carlo evaluation of K(892) production channel.

  sigRand -= sigin1;
  int iselect=0;
  double siga=sigin1;
  if(sigRand <=0.0) {
    iselect=1;
  } else if(jch == 2) {
     sigRand -= sigin2;
     if(sigRand <=0.0) iselect=2;
  }

  if(iselect == 0) return;

  int kf1=313;
  int kf2=kfb;
  // k0 + N -> k0 + N
  // k+ + N -> k+ + N
  if(iselect ==  1) {
      if(izm==0) kf1=313;
      if(izm==1) kf1=323;
      kf2=kfb;
  } else {
    siga=sigin2;
    int kf=abs(kfb);
    int kflr2=(kf/10000)%100;
    int kfla2=(kf/1000)%10;
    int kflb2=(kf/100)%10;
    int kflc2=(kf/10)%10;
    int kb=100*kfla2+10*kflb2+kflc2;
    int kfls2=kf%10;

    // k0 + p   -> k+* + n
    // k0 + d++ -> k+* + d+
    if(izm==0) {
      kf1=323;

      // Some N* has the ID such as N(1720)0 = 31214 N(1720)+ = 32124
      if(izb==1) {
        if(kb==221) kf2=10000*kflr2 + 2110 + kfls2;
	else if(kb==212) kf2=10000*kflr2 + 1210 + kfls2;
      } else if(izb==2) {
	kf2=delta2toP[kfb];
        //kf2=10000*kflr2 + 2210 + kfls2;
      }

    // k+ + n  -> k0  + p
    // k+ + d- -> k0  + d0
    } else if(izm == 1)  {
      kf1=313;
      if(izb==0) {
        if(kb==211) kf2=10000*kflr2 + 2210 + kfls2;
	else if(kb==121) kf2=10000*kflr2 + 2120 + kfls2;
      } else if(izb==-1)  {
	kf2=deltaMto0[kfb];
        //kf2=10000*kflr2 + 2110 + kfls2;
      }
    }
  }

  pd1=jamtable->find(kf1);
  double emr=pd1->m0();
  double emmin=pd1->mMin();
  double wid=pd1->mWidth();
  double em2 = cpair->M(1);
  double emmax=srt-em2-ekinmi;
  //double em1 = decay->BWMass(emmin,emmax,emr,wid);
  double em1 = jamtable->BWMass(emmin,emmax,emr,wid);
  if(em1 <=0.1) {
      cout << "sigKNKstar K* mass is wrong em1= " << em1 <<endl;
      cout << " emmin= " << emmin << " emmax= " << emmax<<endl;
      exit(1);
  }
  if(cpair->isAnti()) { kf1 *= -1; kf2 *= -1; }
  cpair->setOutGoing(siga,kf1,kf2);
  cpair->setOutGoingMass(em1,em2);


}

// (5) K*(892) absorption using detailed balance.
void SigmaKB::absorbKstar()
{
  double pfinal=0.0;
  if(srt > emk+emnuc) {
    pfinal=XsecTable::pawt(srt,emk,emnuc);
  } else {
      return;
  }

  double widcof=max(1.0,M_PI/(atan(2.0*(srt-emnuc-0.892)/widk)-
               atan(2.0*(emk+empion-0.892)/widk)));

  int jch=1;
  double sigin1=0.0,sigin2=0.0;
  if(ibra <= 2 || ibra == 5 || ibra==6) {
    sigin1=jamsigkn(2,srt,plab)*pow2(pfinal/pr)*widcof;
    sigin[mChanel]=sigin1;
  } else {
    jch=2;
    sigin1=jamsigkn(3,srt,plab)*pow2(pfinal/pr)*widcof;
    sigin[mChanel]=sigin1;
    sigin2=jamsigkn(4,srt,plab)*pow2(pfinal/pr)*widcof;
    sigin[mChanel]=sigin2;
    mChanel++;
  }

  mChanel++;
  if(mSel !=3) return;

  //....Monte Carlo evaluation of K(892) absorption channel.
  sigRand -= sigin1;
  double siga=sigin1;
  int iselect=0;
  if(sigRand <= 0.0) {
    iselect=1;
  } else if(jch == 2) {
    sigRand -= sigin2;
    if(sigRand <= 0.0) iselect=2;
  }

  if(iselect == 0) return;

  int kf1=0,kf2=0;
  if(iselect==1) {
    if(izm==0) kf1=311;
    if(izm==1) kf1=321;
    kf2=kfb;
  } else {
    siga=sigin2;
    int kf=abs(kfb);
    //int kflr2=(kf/10000)%10;
    int kflr2=(kf/10000)%100;
    int kfla2=(kf/1000)%10;
    int kflb2=(kf/100)%10;
    int kflc2=(kf/10)%10;
    int kb=100*kfla2+10*kflb2+kflc2;
    int kfls2=kf%10;
    // k0 + p -> k+* + n
    // k0 + d++ -> k+* + d+
    if(izm==0) {
      kf1=321;
      if(izb==1) {
        if(kb==221) kf2=10000*kflr2 + 2110 + kfls2;
	else if(kb==212) kf2=10000*kflr2 + 1210 + kfls2;
      } else if(izb==2) {
	kf2=delta2toP[kfb];
        //kf2=10000*kflr2 + 2210 + kfls2;
      } else {
	  cout << "absorbKstar 1 izb= " << izb <<endl;
	  exit(1);
      }
    // k+ + n -> k0+ + p
    // k+ + d- -> k0  + d0
    } else if(izm == 1)  {
      kf1=311;
      if(izb==0) {
        if(kb==211) kf2=10000*kflr2 + 2210 + kfls2;
	else if(kb==121) kf2=10000*kflr2 + 2120 + kfls2;
      } else if(izb==-1) {
	kf2=deltaMto0[kfb];
        //kf2=10000*kflr2 + 2110 + kfls2;
      } else {
	  cout << "absorbKstar 2 izb= " << izb <<endl;
	  exit(1);
      }
    }
  }

  double em1=0.0;
  if(kf1==311) em1=emk0;
  else if(kf1==321) em1=emkc;
  else {
      cout << "absorbKstar kf1=0? kf1= " << kf1<<endl;
      exit(1);
  }
  double em2 = cpair->M(1);
  if(cpair->isAnti()) { kf1 *= -1; kf2 *= -1; }
  cpair->setOutGoing(siga,kf1,kf2);
  cpair->setOutGoingMass(em1,em2);

}

void SigmaKB::KstarD()
{
  // (6) K(892)+D(1232) production.

  double sigt1=jamsigkn(5,snew,plab);
  sigin[mChanel]=sigt1;
  mChanel++;
  if(mSel !=3) return;

  double sigin1,sigin2;
  if(ibra==1 || ibra == 2) {
    sigin1=sigt1*0.25;
    sigin2=sigt1*0.75;
  } else if(ibra==3 || ibra==4 || ibra==7 || ibra==8) {
    sigin1=sigt1*0.5;
    sigin2=sigt1*0.5;
  } else {
    sigin1=sigt1;
    sigin2=0.0;
  }

  sigRand -= sigin1;
  int iselect=-1;
  if(sigRand <= 0.00) {
    iselect=0;
  } else {
    sigRand -= sigin2;
    if(sigRand <= 0.00) iselect=1;
  }

  if(iselect == -1) return;

  int kf1=kfdelt[ibra-1][iselect][0] + 2;  // K* = 313/323
  int kf2=kfdelt[ibra-1][iselect][1];

  if(kf1==0 || kf2==0) {
    cout << " kfm= " << kfm << " kfb= " << kfb<<endl;
    cout << "kf1= " << kf1 << " kf2= " << kf2 <<endl;
    cout << " ibra= " << ibra << " iselect= " << iselect<<endl;
    exit(1);
  }

  pd1=jamtable->find(kf1);
  pd2=jamtable->find(kf2);
  double emmin1=pd1->mMin();
  double emmin2=pd2->mMin();
  if(cpair->isAnti()) { kf1 *= -1; kf2 *= -1; }

  if(srt < emmin1+emmin2+ekinmi) {
    cout << " KstarD srt ? " << srt << " th= " << emmin1+emmin2+ekinmi<<endl;
    cpair->setOutGoing(sigt1,kf1,kf2);
    cpair->setOutGoingMass(emmin1,emmin2);
    return;
  }

  double emmax1=srt-emmin2-ekinmi;
  double emr1=pd1->m0();
  double wid1=pd1->mWidth();
  double emr2=pd2->m0();
  double wid2=pd2->mWidth();
  double em1=0.0, em2=0.0;

  int itry=0;
  do {
    //em1=decay->BWMass(emmin1,emmax1,emr1,wid1);
    em1=jamtable->BWMass(emmin1,emmax1,emr1,wid1);
    double emmax2=srt-em1-ekinmi;
    if(emmax2 > emmin2) {
      //em2=decay->BWMass(emmin2,emmax2,emr2,wid2);
      em2=jamtable->BWMass(emmin2,emmax2,emr2,wid2);
      break;
    } else {
      if(itry <= 50) continue;
      em1=emmin1;
      em2=emmin2;
      break;
    }
  } while(itry++ < 100);

  cpair->setOutGoing(sigin[mChanel-1],kf1,kf2);
  cpair->setOutGoingMass(em1,em2);
     
}

// (7) K*(892)D(1232) absorption using detailed balance.
void SigmaKB::KstarDabsorb()
{
  double widcof1=max(1.0,M_PI/(atan(2.0*(srt-emnuc-0.892)/widk)-
               atan(2.0*(emk+empion-0.892)/widk)));

  double widcof2=max(1.0,M_PI/(atan(2.0*(srt-emk-1.232)/widdlt)-
               atan(2.0*(emnuc+empion-1.232)/widdlt)));

  double pfinal=0.0;
  if(srt > emk+emnuc) {
    pfinal=XsecTable::pawt(srt,emk,emnuc);
  } else
      return;

  double facsp=4.0/(double(kfm%10)*double(kfb%10));
  sigin[mChanel]=facsp*pow2(pfinal/pr)*jamsigkn(5,srt,plab)*widcof1*widcof2;
  mChanel++;
  if(mSel != 3) return;

  sigRand -= sigin[mChanel-1];
  if(sigRand > 0.0) return;

  int kf1=0,kf2=0;
  switch (izb) {
     case -1: // D-
       kf1=311;
       kf2=2112;
       break;
     case 0: // D0
       if(izm==1) {        // d0 k+
         if(rndm->flat() <= 0.5) {
           kf2=2112;
           kf1=321;
	 } else {
           kf2=2212;
           kf1=311;
	 }
       } else if(izm == 0) { // d0 k0
         kf2=2112;
         kf1=311l;
       }
       break;

     case 1: // D+
       if(izm == 0) {       // d+ k0
         if(rndm->flat() <= 0.5) {
           kf2=2112;
           kf1=321;
	 } else {
           kf2=2212;
           kf1=311;
	 }
       } else if(izm == 1) {  // d+ k+
         kf2=2212;
         kf1=321;
       }
       break;

     case 2: // D++
       kf2=2212;
       kf1=321;
       break;
  }

  double em1,em2;
  if(kf1 == 311) em1=emk0;
  else if(kf1 == 321) em1=emkc;
  else {
      cout << "SigmaKB::KstarDabsorb kf1=0? kf1= " << kf1 <<endl;

      exit(1);
  }
  if(kf2 == 2112) em2=emn;
  else if(kf2 == 2212) em2=emp;
  else {
      cout << "SigmaKB::KstarDabsorb kf2=0? kf2= " << kf2 <<endl;
      exit(1);
  }
     
  if(cpair->isAnti()) { kf1 *= -1; kf2 *= -1; }
  cpair->setOutGoing(sigin[mChanel-1],kf1,kf2);
  cpair->setOutGoingMass(em1,em2);

}

//***********************************************************************

double SigmaKB::jamsigkn(int isig,double s,double pl)
{
//...Purpose: to give kaon-nucleon one/two-pion production x-sections.
//...1:K+ N => K delta  isospin=1
//...2:K+ p => K+*(890) p
//...3:K+ n => K+*(890) n
//...4:K+ n => K0*(890) p
//...5:K+ N => D(1232)+K(892) T=1

  double sth[5]={1.5719, 1.5719, 1.5719, 1.5719, 1.71231};
//double a1[5]={0.584417,0.519173,0.309147,0.465779,.262122};
  double a1[5]={0.584417,0.519173,0.309147,0.465779,.349496};
  double a2[5]={1.09807,1.50663,1.37519,0.720443,0.929574};
  double a3[5]={1.83136,1.95501,1.97189,1.99835,2.50283};
  double a4[5]={0.0365562,0.0531106,0.0972668,0.0681271,0.143572};

  double plth[5]={2.5, 2.5, 2.5, 2.5, 4.0};
// dobule b1[5]={5.44791,6.85573,2.40975,14.0587,11.1962};
  double b1[5]={5.44791,6.85573,2.40975,14.0587,14.92827};
  double b2[5]={1.66538,1.70302,1.18874,2.28464,1.79165};

  if(s < sth[isig-1]) {
    return 0.0;
  } else if(pl <= plth[isig-1]) {
    int i=isig-1;
    return a1[i]*pow(s-sth[i],a2[i])/(pow2(s-a3[i]) + a4[i]);
  } else {
    return b1[isig-1]*pow(pl,-b2[isig-1]);
  }

}


} // end of namespace jam2


