#include <jam2/xsection/SigmaABB.h>
#include <jam2/xsection/XsecTable.h>
#include <jam2/hadrons/JamStdlib.h>

using namespace std;

namespace jam2 {

//       call jamcabb(msel,srt,pr,kf1,kf2,kc1,kc2,
//    $         em1,em2,sig,sigel,sigin,mchanel,mabsrb,mxchan,ijet)
//
//     subroutine jamcabb(msel,srt,pr,kf1,kf2,kc1,kc2,
//    $            em1,em2,sig,sigel,sigin,mchanel,mabsrb,mxchan,ijet)

const double Mnucl=0.93895; //  parc(28)
const double Mp= 0.93827;
const double Mn= 0.93957;
const double sMin=2.018;

//***********************************************************************
void SigmaABB::calc(CollisionPair& cpair)
{
//...Purpose: to treat baryon-antibaryon cross sections.

    mSel=cpair.mode();

    double srt = cpair.getCMenergy();
    double pr  = cpair.getCMmomentum();
    double em1=cpair.M(0);
    double em2=cpair.M(1);
    int iz1=cpair.getZ(0)/3;
    int iz2=cpair.getZ(1)/3;
    int kf1=cpair.getID(0);
    int kf2=cpair.getID(1);
    int id1=cpair.getPID(0);
    int id2=cpair.getPID(1);
    int kf1a=abs(kf1);
    int kf2a=abs(kf2);
    //pa1=cpair.getJParticle(0);
    //pa2=cpair.getJParticle(1);
    //int ipair=cpair.getPairID();
  bool preHadronA = cpair.qFactor(0) < 1.0 ? true: false;
  bool preHadronB = cpair.qFactor(1) < 1.0 ? true: false;

    sigTot=0.0;
    sigEl=0.0;
    sigAnn=0.0;
    double snew=2*sqrt(Mnucl*Mnucl+pr*pr);
    double plab=0.0;
    if(snew > 2*Mnucl) {
        plab=XsecTable::plabsr(snew,Mnucl,Mnucl);
    } else {
        cout << "SigmaABB::calc (jamcabb:)plab<0 " << snew
	     << " srt= "<< srt
	     << " pr= "<< pr
	     << " kf1  "<<kf1
	     << " kf2  "<<kf2
	     << " m1  "<< em1
	     << " m2  "<< em2
	     << endl;
        return;
    }

//...Get total and elastic cross sections.
    getXsection(snew,plab,iz1,iz2);

//...Rescale the cross section.
    double sig1,sigel1,sig2,sigel2;
    XsecTable::additiveQuarkModel(2212,2212,sig1,sigel1);
    XsecTable::additiveQuarkModel(kf1,kf2,sig2,sigel2);
    sigTot *= sig2/sig1;
    sigEl  *= sigel2/sigel1;


  if(!preHadronA && !preHadronB) {

//...Calculate annihilation cross section.
//...Parametrization by
//...J.Gugnon and J. Vandermeulen, Ann. Phys. (Paris) 14, (1989)49.
    //sigAnn=min(sigTot-sigEl,sig2/sig1*(24/pow(plab,1.1)+38/sqrt(plab)));
    //sigAnn=sig2/sig1*(24/pow(plab,1.1)+38/sqrt(plab));
    sigAnn=sig2/sig1*(21/pow(plab,1.1)+38/sqrt(plab));

    if(kf1a==kf2a) {
      sigAnn = sig2/sig1*sigPPAnn(plab);
    } else {
      if(plab < 0.382) {
        sigAnn = sig2/sig1*( 29.0/plab + 41.4 );
      } else {
        sigAnn = sig2/sig1*sigPPAnn(plab);
      }
    }
    
  }

    /*
//...Quark contents.
    int kfl1[3],kfl2[3];
    kfl1[2]=(kf1a/1000)%10;
    kfl1[1]=(kf1a/100)%10;
    kfl1[0]=(kf1a/10)%10;

    kfl2[2]=(kf2a/1000)%10;
    kfl2[1]=(kf2a/100)%10;
    kfl2[0]=(kf2a/10)%10;

    int kfl1a=max(kfl1[0],max(kfl1[1],kfl1[2]));
    int kfl1c=min(kfl1[0],min(kfl1[1],kfl1[2]));
    int kfl1b=kfl1[0]+kfl1[1]+kfl1[2]-kfl1a-kfl1c;

    int kfl2a=max(kfl2[0],max(kfl2[1],kfl2[2]));
    int kfl2c=min(kfl2[0],min(kfl2[1],kfl2[2]));
    int kfl2b=kfl2[0]+kfl2[1]+kfl2[2]-kfl2a-kfl2c;

    int kfla,kflb;
    if(kfl1a == kfl2a && kfl1b == kfl2b) {
        kfla=kfl1c;
        kflb=kfl2c;
    } else if(kfl1a == kfl2a && kfl1c == kfl2c) {
        kfla=kfl1b;
        kflb=kfl2b;
    } else if(kfl1b == kfl2b && kfl1c == kfl2c) {
        kfla=kfl1a;
        kflb=kfl2a;
    } else {
        sigTot = sigTot-sigAnn;
        sigAnn=0.0;
    }
    */


    //if(mSel == 1) return;

    mChanel=1;
    mAbsrb=1;
    sigin[0]=max(0.0,sigTot-sigEl-sigAnn);
    sigin[1]=sigAnn;

    sigRand=cpair.getXsig();

    int ianti1 = kf1 < 0 ?  -1 : 1;
    int ianti2 = kf2 < 0 ?  -1 : 1;

    // Charge exchange.
    sigin[2]=0.0;
    if(id1==id_nucl && id2==id_nucl) {
    if(kf1a==kf2a) {
	mChanel++;
	if(plab < 0.5) {
	    sigin[2]=max(0.0,10.9*(plab-0.1)*pow(plab,-1.6));
	} else {
           sigin[2]=7.1*pow(plab,-0.9);
	}
    }
    }

    // t-channel resonance or string formation
    sigin[3]=sigin[0]-sigin[2];

    cpair.setChargeEx(sigin[2]);
    if(srt < sMin) sigTot = sigEl + sigAnn + sigin[2];

    cpair.setXSTotal(sigTot);
    //cpair.setXSElastic(sigTot); // test:elastic only
    cpair.setXSElastic(sigEl);
    cpair.setSigAbs(sigAnn);

    if(mSel == 1) return;

    sigRand -=sigin[2];
    if(sigRand < 0.0) {
	if(kf1a==2212) {
	  cpair.setOutGoing(sigin[2],2112*ianti1,2112*ianti2);
	  cpair.setOutGoingMass(Mn,Mn);
	} else {
	  cpair.setOutGoing(sigin[2],2212*ianti1,2212*ianti2);
	  cpair.setOutGoingMass(Mp,Mp);
	    
	}
	return;
    }

    // Annihilation 
    sigRand -= sigAnn;
    if(sigRand < 0.0) {
	cpair.setOutGoing(sigAnn,93);
	return;
    }

    // t-channel string formation.
    //if(srt >= 3.5) {
    if(srt >= ecmStringABB) {
	cpair.setOutGoing(sigin[3],92,92);

    // t-channel resonance formation.
    } else {
	Pythia8::ParticleDataEntry *pd1=jamtable->find(kf1);
	Pythia8::ParticleDataEntry *pd2=jamtable->find(kf2);
  //bool preHadronA = cpair.qFactor(0) < 1.0 ? true: false;
  //bool preHadronB = cpair.qFactor(1) < 1.0 ? true: false;
        int icon=sampleMass->jamrmas2(cpair.getCollType(),pd1,pd2,kf1,kf2,
                id1,id2,iz1,iz2,em1,em2,preHadronA,preHadronB,srt,cpair.isAnti());

        if(icon==1) {
            cpair.setChannel(-1);
            return;
        }
        cpair.setOutGoing(sigin[3],sampleMass->id(0),sampleMass->id(1));
        cpair.setOutGoingMass(sampleMass->m(0),sampleMass->m(1));
    }

}

//***********************************************************************

void SigmaABB::getXsection(double srt,double plab,int iz1,int iz2)
{
//...Calculate baryon-antibaryon total and elastic cross sections.
    static const double bhad=2.3,eps=0.0808,facel=0.0511;

    int az1=abs(iz1);
    int az2=abs(iz2);

    // S.A.Bass el.al. Prog.Part.Nucl.Phys.41(1998)255
    // M.Belicher, el. al. J.Phys.G Nucl.Part.Phys.25 (1999)1859

    if(srt < 2.5) {
	if(plab < 0.3) {
          //sigTot=271.6*exp(-1.1*plab*plab);
          //sigEl=78.6;
          //sigTot = 355.72 + 193.5*srt;

          sigTot=XsecTable::jamsighh(5,srt);
	  sigEl = 78.6 + 2.0/plab;


	} else {
          sigTot=75.0 + 43.1/plab + 2.6/(plab*plab) - 3.9*plab;
          //sigEl=31.6 + 18.3/plab - 1.1/(plab*plab) -3.8*plab;
          sigEl=30.6 + 18.3/plab - 1.1/(plab*plab) -3.8*plab;
	}

    // JAM table fit
    //if(srt < 4.0) {
    //    sigTot=XsecTable::jamsighh(5,srt);
    //    sigEl=XsecTable::jamsighh(6,srt);

    } else if(srt < 5.0) {
        if(az1 == az2) {
          sigTot=XsecTable::jamchc96(21,plab);
          sigEl=XsecTable::jamchc96(22,plab);
	} else {
          sigTot=XsecTable::jamchc96(23,plab);
          sigEl=XsecTable::jamchc96(22,plab);
	}

    } else if(srt < 50.0) {
        if(az1 == az2) {
          //sigTot=XsecTable::jamchc96(21,plab);
          sigTot=XsecTable::pdg2016(srt,1);
          //sigTot=XsecTable::pdg2016(srt,1);
          sigEl=XsecTable::jamchc96(22,plab);
	} else {
          //sigTot=XsecTable::jamchc96(23,plab);
          sigTot=XsecTable::pdg2016(srt,2);
          sigEl=XsecTable::jamchc96(22,plab);
	}

    } else {
        double bel=4.0*bhad+4.0*pow(srt*srt,eps)-4.2;
        if(az1 == az2) {
          //sigTot=XsecTable::jamrgg96(srt,2);
          sigTot=XsecTable::pdg2016(srt,1);
          sigEl=facel*sigTot*sigTot/bel;
	} else {
          //sigTot=XsecTable::jamrgg96(srt,4);
          sigTot=XsecTable::pdg2016(srt,2);
          sigEl=0.98*facel*sigTot*sigTot/bel;
	}
    }

}

double SigmaABB::sigPPAnn(double plab)
{
  //return 24/pow(plab,1.1)+38/sqrt(plab);
  return 21/pow(plab,1.1)+38/sqrt(plab);

  if(plab <0.51) {
    return 51.52/pow(plab,0.85) + 0.034/pow(plab,2.94);
  } else if(plab < 6.34) {
    return 88.8/pow(plab,0.4) - 24.2;
  } else {
    return 24/pow(plab,1.1)+38/sqrt(plab);
  }

} 

}
