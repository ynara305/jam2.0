#include <cstdlib>
#include <jam2/xsection/SigmaMB.h>
#include <jam2/xsection/XsecTable.h>
#include <jam2/hadrons/JamStdlib.h>

using namespace std;

namespace jam2 {

using Pythia8::ParticleDataEntry;

SigmaMB::SigmaMB(Pythia8::Settings *s,JamParticleData* table, SampleMass* sm)
{
  settings=s;
  jamtable=table;
  decay=jamtable->getDecayWidth();
  sampleMass = sm;
  pRes=0;
  optSChannel=settings->mode("Cascade:optConstQSChannel");
  optJAM1=settings->mode("Cascade:optBBJAM1");
  ecmStringMB=settings->parm("Cascade:ecmStringMB");
  ecmStringMBs=settings->parm("Cascade:ecmStringMBs");
}

//...Treat strangeness S=0 meson-baryon cross section.
void SigmaMB::calcS0(CollisionPair& cpair)
{
    int izm=cpair.getZ(0)/3;
    int izb=cpair.getZ(1)/3;
    double srt=cpair.getCMenergy();
    double pr=cpair.getCMmomentum();
    int kfm=cpair.getID(0);
    int kfb=cpair.getID(1);
    int istr=cpair.getTotalStrangeness();
    int idm=cpair.getPID(0);
    int idb=cpair.getPID(1);

    if(kfb==0 || kfm==0) {
	cout << " kfb=0? " << kfb
	    << " kfm= " << kfm
	    <<endl;
	exit(1);
    }

    int izt=izm+izb;

    mSel=cpair.mode();
    int msel=mSel;
    sigRand=cpair.getXsig();
    //const double sres=3.0, srt0=2.2;
    const double sres=5.0, srt0=2.2;

    double sigres=0.0; // s-channel resonance formation
    double sigin_t=0.0;
    bool absorb=true;
    sigEl = 0.0;
    sigElBW=0.0;
    if(izt < -1 || izt > 2) absorb=false;

    // If meson is within the formation time.
  bool preHadronA = cpair.qFactor(0) < 1.0 ? true: false;
  bool preHadronB = cpair.qFactor(1) < 1.0 ? true: false;
  if(optSChannel==1 && (preHadronA || preHadronB))
	absorb=false;

//....Get total cross section.
//==================================================================

    //ParticleDataEntry* pm;

    //.....pion-N->Delta(1232)
   if( CollisionPair::Pair(idb,idm) ==
       CollisionPair::Pair(id_nucl,id_pi) ) {

	jamxpin(izm,izb,srt,pr); // compute sigTot and sigEl

	// compute total absorption cross section.
	//if(absorb && srt < sres) {
	if(srt < sres) {
	    mSel=1;
	    sigres=jamxbw1(istr,srt,pr,izm,izb,kfm,kfb);
	    mSel=msel;
	}

	//double sigel=sigEl;
	//sigEl=sigel - sigElBW;

	//....t-channel elastic.
	if(optJAM1==0)
	  sigEl = piNtchannelElastic(srt,izt);
	else
	  sigEl = piNtchannelElasticJAM1(srt,izt);


	//sigEl=sigel-sigres;
	//sigEl=0.0;
	//sigEl = sigTot - sigres;

	double srtR=1.7;
        if(izt == 2 || izt == -1) srtR=1.85;
        //if(izt == 2 || izt == -1) srtR=1.95;
	//srt0=1.0;
	// s-channel only
	if(srt <= srtR) {
	    sigTot = sigres + sigEl;
	    if(absorb) {
	      cpair.setXS(sigres+sigEl, sigEl);
	    } else {
	      cpair.setXS(sigEl, sigEl);
	    }
	} else {

	    if(!absorb) sigTot -= sigres;
	    cpair.setXS(max(0.0,sigTot),sigEl);

	}

	/*
	cout << " q1= "<< cpair.qFactor(0) 
	     << " q2= "<< cpair.qFactor(1) 
	    << " abs= "<< absorb << " srt= "<<
	    srt << " sig= "<< sigTot << " sigel= "<< sigEl
	    << " sigab= " << sigres
	     <<endl; 
	cin.get();
	*/

    //....Higher resonance formation.
   } else {

       // Find s-channel resonance formation cross section.
	if(absorb && srt < sres)  {
	    mSel=1;
	    sigres=jamxbw1(istr,srt,pr,izm,izb,kfm,kfb);
	    mSel=msel;
	}
	sigres=min(300.0,sigres);

	double sigt=0.0;
	//if(sigres ==0 || srt > 2.0)
	XsecTable::additiveQuarkModel(kfm,kfb,sigt,sigEl);

	sigTot=sigt+sigres;
	sigin_t=max(0.0,sigt - sigEl); // t-channel inel.

	cpair.setXS(sigTot,sigEl);

    }

    cpair.setElasticBW(sigElBW);
    // test
    //cpair.setXS(sigTot,sigTot);

     //...End in case of total cross section.
    if(mSel==1) return;

//==================================================================

    double sigin_s = 0.0;
    if(srt > srt0) sigin_s = sigres*srt0/srt;  // s-channel string formation

    //if(mSel !=3 ) return;

    // Monte-Carlo for s-channel collision.
    if(absorb) {

	//....Monte Carlo for resonance formation.
	if(srt < sres)  {
	    double sig=jamxbw1(istr,srt,pr,izm,izb,kfm,kfb,cpair.isAnti());
	    if(sigRand <= 0.0) {
		cpair.setOutGoing(sig,kfm);
		cpair.setOutGoingMass(srt);
		//cpair.setOutGoingParticle(pm);

	//cout << " pi-N -> R kfm= " << kfm <<" srt= " << srt <<endl;
	 //     cout << cpair.getOutGoing().id1 <<endl;

		return;
	    }
	}

	//....Monte Carlo for s-channel string formation.
	sigRand -= sigin_s;
	if(sigRand<=0.0) {
	//if(srt > 2.0) {

	  int kf;
          if(izt == -1) kf=1114;
          if(izt ==  0) kf=12112;
          if(izt ==  1) kf=12212;
          if(izt ==  2) kf=2224;
	  if(cpair.isAnti()) kf *= -1;

	  //cpair.setOutGoing(sigin_s,92);
	  cpair.setOutGoing(sigin_s,kf);
	  cpair.setOutGoingMass(srt);
          cpair.setChannel(-2);

	  if(srt < 2.0) {
	      cout << "SigmaMB::calcS0 S=0 s-channel string formation srt="<< srt <<endl;
	      cout << " izt = " << izt << " kf= " << kf <<endl;
	      exit(1);
	  }
	  return;
	}
    }


    // test: elastic
   //cpair.setChannel(-1); return;

    double stmin=1.3;
    if(cpair.getStrange(0) !=0 && cpair.getStrange(1) !=0) stmin=1.9;
    // make it elastic
    if(srt < stmin) { cpair.setChannel(-1); return; }
	

    // Monte-Carlo for t-channel contribution.
    double sTh=3.0;
    // phi induced etc.
    if((abs(kfm)/100)%10 ==3 && (abs(kfm)/10)%10==3) sTh = 3.2;
    // K+ Lambda  etc.
    if(cpair.getStrange(0)*cpair.getStrange(1) !=0) sTh = 3.2;
    double sigtr=sigin_t*sTh/srt;
    sigRand -= sigtr;

    //if(srt > 4.0 || (srt > sTh && sigRand < 0.0)) {
    //if(srt > 3.0 || (srt > sTh && sigRand < 0.0)) {
    if(srt > ecmStringMB || (srt > sTh && sigRand < 0.0)) {
    // t-channel string formation.
      cpair.setOutGoing(sigin_t-sigtr,92,92);
      return;
    }

    // test: elastic
    //cpair.setChannel(-1);return;

    // t-channel resonance formation.
    //if(srt <= 3.0 && sigRand <= 0) {
	double m1=cpair.M(0), m2=cpair.M(1);
	//pd1 = cpair.getParticleDataEntry(0);
	//pd2 = cpair.getParticleDataEntry(1);
	pd1=jamtable->find(kfm);
	pd2=jamtable->find(kfb);
	int icon=sampleMass->jamrmas2(cpair.getCollType(),pd1,pd2,kfm,kfb,
		idm,idb,izm,izb, m1,m2,preHadronA,preHadronB,srt,cpair.isAnti());

	if(icon==1) { cpair.setChannel(-1); return; }

        cpair.setOutGoing(sigtr,sampleMass->id(0),sampleMass->id(1));
	cpair.setOutGoingMass(sampleMass->m(0),sampleMass->m(1));
	//cpair.setOutGoingParticle(sampleMass->p(0),sampleMass->p(1));

	/*
       cout << " MB t-channel collision srt= " << srt 
	   << " sigr= " << sigtr <<endl;
       cout << " idm= " << idm << " " << pd1->name() << " m1= " << m1
	   << " + " << "idb= " << idb << " " << pd2->name()
	   << " m2= " << m2 <<endl;
       cout << " id3= " << sampleMass->id(0)
	    << " " << sampleMass->p(0)->name()
	    << " m3= " <<  sampleMass->m(0)
            << " id4= " << sampleMass->id(1)
	    << " " << sampleMass->p(1)->name()
	    << " m4= " <<  sampleMass->m(1)
	    <<endl;
	    */


    //}


}

// pi-n t-channel elastic.
double SigmaMB::piNtchannelElastic(double srt, int izt)
{
  if(izt == 2 || izt == -1) {
    // JAM2
    if(srt < 1.5) return 0.0;
    else if( srt <=  4.0) 
    return 11.9969*pow(srt-1.5,1.11007 )/(pow2(srt-1.33585)+0.222022);
    else return sigEl;

  } else {
    // JAM2
    if(srt < 1.65) return 0.0;
    else if(srt <= 4.0)
    return 27.0714*pow(srt-1.65,0.871151)/(pow2(srt-0.281746)-1.67436);
    else return sigEl;
  }
}

// pi-n t-channel elastic for JAM version1.
double SigmaMB::piNtchannelElasticJAM1(double srt, int izt)
{
  if(izt == 2 || izt == -1) {
    // JAM1
    if(srt <= 1.4) return 0.0;
    else if( srt <=  4.0) 
    return 0.90*3.81*pow(srt-1.4,2.156)/(pow2(srt-1.66)+0.0352);
    else return sigEl;

  } else {
    // JAM1
    if(srt <= 1.4) return 0.0;
    else if( srt <=  4.0) 
    return 45.49*pow(srt-1.4,0.854)/(pow2(srt+1.22)-5.585);
    else return sigEl;
  }
}


//***********************************************************************

// Purpose: to give pi-N total and elastic cross sections.
void SigmaMB::jamxpin(int izpi,int izn,double srt,double pr)
{
    static const double emnuc=0.93957,empi=0.13957;
    static const double bhad1=2.3,bhad2=1.4,eps=0.0808,facel=0.0511;

    sigTot=0.0;
    sigEl=0.0;

//...Too slow. no reaction
    //if(srt < emnuc+empi+0.003) return;
    if(srt < emnuc+empi+eKinMin) return;

//....Cross section for delta formation from fitted Breit-Wigner.
    //if(srt <= 1.3) {
    //   sigTot=jamxdelt(izpi,izn,srt,pr);
    //   sigEl=0.0;
    //   return;
    //}
 
    int icha;
//...pi- n/pi+ p
    if((izpi == -1 && izn == 0) || (izpi == 1 && izn == 1))
        icha=1;
//...pi+ n/pi- p
    else if((izpi == -1 && izn == 1) || (izpi == 1 && izn == 0))
        icha=2;
//...pi0 n/pi0 p
    else if((izpi == 0 && izn == 0) || (izpi == 0 && izn == 1))
        icha=3;
    else {
        cout << "(jamxpin:)error izpi= " << izpi << " izn= " << izn << endl;
        return;
    }

//...Cross section from table fit.
    if(srt <= 3.0) {
        if(icha == 1) {
          sigTot=XsecTable::jamsighh(7,srt);
          sigEl=XsecTable::jamsighh(8,srt);
	} else if(icha == 2) {
          sigTot=XsecTable::jamsighh(9,srt);
          sigEl=XsecTable::jamsighh(10,srt);
	} else {
          double sig1 = XsecTable::jamsighh(7,srt);
          double sig2 = XsecTable::jamsighh(9,srt);
          double sigel1 = XsecTable::jamsighh(8,srt);
          double sigel2 = XsecTable::jamsighh(10,srt);
          sigTot=0.5*(sig1+sig2);
          sigEl=0.5*(sigel1+sigel2);
	}

	if(sigTot <=0.0) {
	    cout <<icha << " srt= "<< srt << " sigtot= "<< sigTot<<endl;
	    exit(1);
	}

//...Cross section from HERA fit.
    } else if(srt <= 5.6) {
	//if(srt<emnuc+empi) return;
        double plab=XsecTable::plabsr(srt,empi,emnuc);
	//double a=0.92;
	double a=0.88;
        if(icha == 1) {
          sigTot=XsecTable::jamchc96(3,plab);
          sigEl=a*XsecTable::jamchc96(4,plab);
	} else if(icha == 2) {
          sigTot=XsecTable::jamchc96(5,plab);
          sigEl=XsecTable::jamchc96(6,plab);
	} else {
          sigTot=0.5*(XsecTable::jamchc96(5,plab)+XsecTable::jamchc96(3,plab));
          sigEl=0.5*(XsecTable::jamchc96(6,plab)+a*XsecTable::jamchc96(4,plab));
	}

//...Regge fit.
    } else {
	double a=1.0;
        if(icha == 1) {
          //sigTot=XsecTable::jamrgg96(srt,7);
          sigTot=XsecTable::pdg2016(srt,-3);
	  a=1.1;
	} else if(icha == 2) {
          //sigTot=XsecTable::jamrgg96(srt,8);
          sigTot=XsecTable::pdg2016(srt,3);
	} else {
          //sigTot=0.5*(XsecTable::jamrgg96(srt,7)+XsecTable::jamrgg96(srt,8));
          sigTot=0.5*(XsecTable::pdg2016(srt,-3)+XsecTable::pdg2016(srt,3));
	}
        double bel=2.0*bhad1+2.0*bhad2+4.0*pow(srt*srt,eps)-4.2;
        sigEl=a*facel*sigTot*sigTot/bel;
    }

}

//***********************************************************************

double SigmaMB::jamxdelt(int izpi,int izn,double srt,double pr)
{
//...Purpose: to calculate pion plus nucleon to delta cross section.
//... ref. Gy.Wolf et al., Nucl.Phys. A517 (1990) 615
//...      A.Engel et al., Nucl.Phys. A572 (1994) 657

    static const double sig0=135.0, sig1=200.0, sig2=70.0;
    static const double gmr=0.11, be=0.3, be2=be*be;
    //static const double empion=0.138,emnuc=0.9383,emdelt=1.232;
    static const double emdelt=1.232;
    static const double  qr=0.227894932;
    static const double v2=be2/(be2+qr*qr);
//  qr=sqrt((emdelt**2-(emnuc+empion)**2)
//       *(emdelt**2-(emnuc-empion)**2))/(2*emdelt) )

    int izzz=(izn+izpi+2)/2;
    double  weight=sig0;
    if(izzz != 1) weight=sig1;
    if((izzz == 1) && (izpi != 0)) weight=sig2;
    double v1=be2/(be2+pr*pr);
    double gamma=pow2(pow3(pr/qr)*(emdelt/srt)*pow2(v1/v2)*gmr)/4;
    return weight*pow2(qr/pr)*gamma/(pow2(srt-emdelt)+gamma);

}



//***********************************************************************

//      subroutine jamcmbs1(msel,srt,pr,istr,kfb,kfm,kcb,kcm,izb,izm,
//     $ emmes,embar,kf1,kf2,em1,em2,sig,sigel,sigin,
//     $ mabsrb,mchanel,mxchan,ijet,jswap,icon)

void SigmaMB::calcS1(CollisionPair& cpair)
{
//...Strangeness S=-1 meson-baryon(K-p,..) collisions.

//...Local arrays.
    double sigy[4],sigyv[4];
    //static const double srt0=1.78;
    //int kfo[2][20];

    double srt=cpair.getCMenergy();
    double pr=cpair.getCMmomentum();
    int izm=cpair.getZ(0)/3;
    int izb=cpair.getZ(1)/3;
    int kfm=cpair.getID(0);
    int kfb=cpair.getID(1);
    int istr=cpair.getTotalStrangeness();
    //mSel=1;
    mSel=cpair.mode();

    /*
//...k-n
      data ((kfy(1,i,j),i=1,2),j=1,4)/
     $   -211,3122, 111,3112, -211,3212, 0,0/
//...k0p
      data ((kfy(2,i,j),i=1,2),j=1,4)/
     $    211,3122, 211,3212, 111,3222, 0,0/
//...k-p
      data ((kfy(3,i,j),i=1,2),j=1,4)/
     $     111,3122, 211,3112, 111,3212, -211,3222/
//...k0n
      data ((kfy(4,i,j),i=1,2),j=1,4)/
     $    111,3122, 211,3112, 111,3212, -211,3222/
     */

    //data (kfkp(1,i),i=1,2)/-321,2112/
    //data (kfkp(2,i),i=1,2)/-311,2212/
    //data (kfkp(3,i),i=1,2)/-321,2212/
    //data (kfkp(4,i),i=1,2)/-311,2112/

    int kfy[4][4][2]={
        {{-211,3122},{111,3112}, {-211,3212},{0,0}},  // K-N
	{{211,3122}, {211,3212}, {111,3222}, {0,0}},    // K0p
        {{111,3122}, {211,3112}, {111,3212}, {-211,3222}}, // K-P
	{{111,3122}, {211,3112}, {111,3212}, {-211,3222}}}; // K0n

    static const double emkc=.49360, emk0=.49770;
    static const double emp=0.93830 ,emn=0.93960;

    double kfkp[4][2]={{-321,2112},{-311,2212},{-321,2212},{-311,2112}};
    double mkp[4][2]={{emkc,emn},{emk0,emp},{emkc,emp},{emk0,emn}};

    //int kf1=kfb;
    //int kf2=kfm;
    sigin[0]=0.0;  // t-channel inel.
    sigin[1]=0.0;  // s-channel resonance formation
    sigin[2]=0.0;  // s-channel string formation
    mChanel=1;
    mAbsrb=2;
    int izt=izb+izm;
    bool absorb=false;
    if(izt >=-1 && izt <= 1) absorb=true;

    // If meson is in the formation time.
    //if(optSChannel==1 && cpair.qFactor(0) < 1.0) absorb=false;
    if(optSChannel==1 && (cpair.qFactor(0) < 1.0 || cpair.qFactor(1)<1.0))
	absorb=false;

    double srtR=1.8;
    double emmes,embar;
    int ik=-1;
//......K- n/antiK0p
    if(kfm == -321 && kfb == 2112) {ik=0;emmes=emkc;embar=emn;srtR=1.8;}
    if(kfm == -311 && kfb == 2212) {ik=1;emmes=emk0;embar=emp;srtR=1.8;}
    if(kfm == -321 && kfb == 2212) {ik=2;emmes=emkc;embar=emp;srtR=1.8;}
    if(kfm == -311 && kfb == 2112) {ik=3;emmes=emk0;embar=emn;srtR=1.8;}

//...Initialize cross sections.
    //int noutpa=0;
    double siga=0.0;
    double sigela=0.0;
    double sigch=0.0;
    for(int i=0;i<4;i++) sigy[i]=0.0;

//....Get t-channel elastic, charge exchange and piY cross section.

    if(ik >= 0) {
	jamxkp(kfm,kfb,srt,emmes,embar,siga,sigela,sigch,sigy);

        // Save outgoing particle types.
        /*
	for(int i=0;i<4;i++) {
          kfo[0][noutpa]=kfy[ik][i][0];
          kfo[1][noutpa]=kfy[ik][i][1];
          noutpa++;
	}
	*/
    }

//...Find resonance formation cross sections.
    double sigr=0.0, sigs=0.0;
    //ParticleDataEntry* pm;
    sigRand=cpair.getXsig();
    if(absorb && srt <= 5.0) {
        sigr=jamxbw1(istr,srt,pr,izm,izb,kfm,kfb,cpair.isAnti());
	cpair.setSigAbs(sigr);
	if(mSel == 3 && sigRand <= 0.0) {
	  cpair.setOutGoing(sigr,kfm);
	  cpair.setOutGoingMass(srt);
	  //cpair.setOutGoingParticle(pRes);
          return;
	}
  }

    sigin[1]=sigr;
    sigin[0]=sigch+sigy[0]+sigy[1]+sigy[2]+sigy[3];
    double sig=sigela+sigin[0]+sigin[1];
    double sigel=sigela;                // background t-channel elastic.

//...t-channel piY->antKN cross section from detailed balance.
    int kfi1,kfi2;
    double sigyp=0.0;
    double emi1,emi2;
    if(absorb) {

//.....Find whether this is pi Y ingoing channel.
	for(int i=0;i<4;i++) {
	    for(int j=0;j<4;j++) {
		kfi1=kfy[i][j][0];
		kfi2=kfy[i][j][1];
		if(kfi1 == kfm && kfi2 == kfb) {
		    //.....Outgoing particle codes and masses.
		    kfi1=kfkp[i][0];
		    kfi2=kfkp[i][1];
		    emi1=mkp[i][0];
		    emi2=mkp[i][1];
		    if(srt < emi1+emi2+0.0001)  goto L200;
		    double prf=XsecTable::pawt(srt,emi1,emi2);
		    double tmp1,tmp2,tmp3;
		    jamxkp(kfi1,kfi2,srt,emi1,emi2,tmp1,tmp2,tmp3,sigyv);
		    int kfam=abs(kfm);
		    int kfab=abs(kfb);
		    double spini=max(1,kfam%10)*max(1,kfab%10);
		    sigyp=sigyv[j]*pow2(prf/pr)*2.0/spini;
		    goto L200;
		}
	    }
	}
L200:
    sigin[0] += sigyp;

    } // end piY->KN

//...Additive quark cross section except KN ingoing above resonance
//...region.
//     if(ik.eq.0.and.srt.ge.srt0) then
      if(ik == -1)
	    XsecTable::additiveQuarkModel(kfm,kfb,siga,sigel);

//...Calculate total and elastic cross sections.
      if(mSel == 1) {
        if(ik >= 0) {
         if(srt <= srtR) {
           sig=sigr+sigela+sigch+sigy[0]+sigy[1]+sigy[2]+sigy[3];
	} else {
           sig=siga;
	}
	} else {
          sig=max(sig+sigyp,siga);
          sigel=sigela;
	}
        sig=min(sig,300.0);
	sigTot=sig;
	sigEl=sigel;
	cpair.setXS(sigTot,sigEl);
	cpair.setElasticBW(sigElBW);
        return;
      }

//=================================================================

  if(absorb && srt > 3.0) sigs = sigr*2.2/srt;
  sigin[3]=sigs;


  if(mSel !=3) return;

//...Monte Carlo for t-channel back ground reactions. Now KN only.
  if(ik >= 0) {

  // First charge exchange reaction.
    sigRand -= sigch;
    if(sigRand <= 0.0) {
      if(izm == -1 && izb == 1) {
            kfm=-311;
            kfb=2112;
	    //emmes=0.49761;
      } else if(izm == 0 && izb == 0) {
            kfm=-321;
            kfb=2212;
      }
      pd1=jamtable->find(kfm);
      pd2=jamtable->find(kfb);
      emmes=pd1->m0();
      embar=pd2->m0();
      if(srt <= emmes+embar+0.001) {
        sigRand += sigch;
	cpair.setChannel(-1);
	return;
      }
      if(cpair.isAnti()) {
	if(pd1->hasAnti()) kfm *= -1;
	if(pd2->hasAnti()) kfb *= -1;
      }
      cpair.setOutGoing(sigch,kfm,kfb);
      cpair.setOutGoingMass(emmes,embar);
      return;
    }

    //...Save outgoing particle types for charge exchange reaction.
    //noutpa++;
    //kfo[0][noutpa]=kfm;
    //kfo[1][noutpa]=kfb;

    //....Monte-Calro for KN->piY.
    for(int i=0;i<4;i++) {
      sigRand -= sigy[i];
      if(sigRand <= 0.0) {
        kfm=kfy[ik][i][0];
        kfb=kfy[ik][i][1];
        pd1=jamtable->find(kfm);
        pd2=jamtable->find(kfb);
        emmes=pd1->m0();
        embar=pd2->m0();
        if(cpair.isAnti()) {
  	  if(pd1->hasAnti()) kfm *= -1;
	  if(pd2->hasAnti()) kfb *= -1;
        }
        cpair.setOutGoing(sigy[i],kfm,kfb);
        cpair.setOutGoingMass(emmes,embar);

	//cout << "KN->piY kfm= " << kfm << " kfb= " << kfb <<endl;
	//cin.get();

	return;
      }
    }
  }

  //....Monte-Carlo for piY->KN.
    sigRand -= sigyp;
    if(sigRand <= 0.0) {
        if(cpair.isAnti()) {
          pd1=jamtable->find(kfi1);
          pd2=jamtable->find(kfi2);
  	  if(pd1->hasAnti()) kfi1 *= -1;
	  if(pd2->hasAnti()) kfi2 *= -1;
        }
      cpair.setOutGoing(sigyp,kfi1,kfi2);
      cpair.setOutGoingMass(emi1,emi2);
      return;

          //noutpa=noutpa+1
          //kfo(1,noutpa)=kfi1
          //kfo(2,noutpa)=kfi2
    }


  // Monte Carlo for s-channel string formation.
  if(absorb) {
      sigRand -= sigs;
      if(sigRand <= 0.0) {
	int kf=3112;
        if(izt == -1) kf=3112;
        if(izt == 0)  kf=3212;
        if(izt == 1)  kf=3222;
	if(cpair.isAnti()) kf *= -1;
	cpair.setOutGoing(sigs,kf);
	cpair.setOutGoingMass(srt);
        cpair.setChannel(-2);
	return;
    }
  }

  //...Gap of cross section from experimental total xsection.
  double sig1=max(0.0,siga-sig);

 //...t-channel resonance formation xsection.
 //double sigtc=sig1*2.5/srt;
 // string impossible due to low energy.
 //if(srt >= srt0 && srt <= 2.5) sigtc=sig1;
  //sigin[1] += sigtc;
  //sig1=max(0.0,sig1-sigtc);

    // Monte-Carlo for t-channel contribution.
  double sigtr=sig1*3.5/srt;
  sigRand -= sigtr;

  // t-channel string formation.
  //if(srt > 3.5 && sigRand < 0.0) {
  //if(srt > 3.5) {
  if(srt > ecmStringMBs) {
    cpair.setOutGoing(sigtr,92,92);
    return;
  }

  if(srt < 1.6) {
    cpair.setChannel(-1);
    return;
  }

  // t-channel resonance.
      pd1=jamtable->find(kfm);
      pd2=jamtable->find(kfb);
      int idm=cpair.getPID(0);
      int idb=cpair.getPID(1);
      double m1=cpair.M(0);
      double m2=cpair.M(1);
  bool preHadronA = cpair.qFactor(0) < 1.0 ? true: false;
  bool preHadronB = cpair.qFactor(1) < 1.0 ? true: false;
      int icon=sampleMass->jamrmas2(cpair.getCollType(),pd1,pd2,kfm,kfb,
		idm,idb,izm,izb, m1,m2,preHadronA,preHadronB,srt,cpair.isAnti());
      if(icon==1) {
        cpair.setChannel(-1);
        return;
      }
      cpair.setOutGoing(sigtr,sampleMass->id(0),sampleMass->id(1));
      cpair.setOutGoingMass(sampleMass->m(0),sampleMass->m(1));
      return;

}

//***********************************************************************

void SigmaMB::jamxkp(int kfm,int kfb,double srt,double emmes,double embar,
	double& sig,double& sigel,double& sigch,double* sigy)
{
//....Purpose: to give anti-K-nucleon total/background elastic
//...background charge exchange and background KN->Ypi cross sections.

    //static const double bhad1=2.3,bhad2=1.4,eps=0.0808,facel=0.0511;
    static const double facel=0.0511;

    sig=0.0;
    sigel=0.0;
    sigch=0.0;
    for(int i=0;i<4;i++) sigy[i]=0.0;
    if(srt < emmes+embar) return;

    double plab=XsecTable::plabsr(srt,emmes,embar);

    //......K- n/antiK0 p i.e. isospin =1 channel.
    if((kfm == -321 && kfb == 2112)
     ||(kfm == -311 && kfb == 2212)) {

	 // K- n Total
	//if (plab <= 2.5) {
	if (plab <= 1.85) {
	  sig=0.9*XsecTable::jamsighh(13,srt);
	} else if(srt < 5.0) {
          sig=XsecTable::jamchc96(14,plab);
	} else {
          sig=XsecTable::pdg2016(srt,5);
	}

	    // K- n t-channel elastic
	//if (plab <= 2.5) {
	/*
	if (plab <= 1.85) {
	    sigel=XsecTable::jamsighh(14,srt);
	} else if(srt < 5.0) {
            //sigel=XsecTable::jamchc96(13,plab);
           double bel=2*(2.3+0.8)+4.0*pow(srt*srt,0.079)-4.2;
           sigel=1.0*facel*sig*sig/bel;
	} else {
            double bel=2*(2.3+0.8)+4.0*pow(srt*srt,0.079)-4.2;
            sigel=facel*sig*sig/bel;
	}
	*/
          double bel=2*(2.3+0.8)+4.0*pow(srt*srt,0.079)-4.2;
          sigel=facel*sig*sig/bel;

        //if(srt <= 2.2) sigel=min(100.0,8.5*pow(plab,-1.0));
        //if(srt > 2.0 && srt <= 2.2) sigel=min(100.0,8.5*pow(plab,-1.0));
        //if(srt <= 2.0) sigel=6.8;

        sigy[0]=min(60.0,2*0.617537*pow(plab,-1.98735)); // k-n->Lpi
        sigy[1]=min(60.0,0.09*pow(plab,-2.42));          // k-n->Spi
        sigy[2]=sigy[1];                                // k-n->Spi

    //.......K-p/antiK0n
    } else if((kfm == -321 && kfb == 2212)
           || (kfm == -311 && kfb == 2112)) {
 
        //if (plab <= 3.0)  {
        if (plab <= 1.3)  {
	    sig=XsecTable::jamsighh(11,srt);
            sigel=XsecTable::jamsighh(12,srt);
	//} else if(srt < 30.0) {
        //      sig=XsecTable::jamchc96(12,plab);
        //      sigel=XsecTable::jamchc96(13,plab);
	} else {
              //sig=XsecTable::jamrgg96(srt,12);
              sig=XsecTable::pdg2016(srt,4);
              double bel=2*(2.3+0.8)+4.0*pow(srt*srt,0.079)-4.2;
              sigel=0.85*facel*sig*sig/bel;
	}

	//...Background elastic.
        if(srt <= 1.55)
          sigel=min(100.0,9.0*pow(plab,-1.3));
        else if(srt <= 1.69)
          sigel=6.0*pow(plab,-1.7);
        else if(srt <= 1.9)
          sigel=13.00*pow(plab,-1.21);
        //else if(srt <= 10.0)
        //  sigel=9.88*pow(plab,-0.637);

	//...Background charge exchange.
        if(srt <= 1.65)
          sigch=min(100.0,1.0*pow(plab,-1.547));
        else if(srt <= 1.87)
         sigch=0.6;
        else
         sigch=0.005;

	//...t-channel hyperon productions(K+N->Lam/Sig+pi).
	sigy[0]=min(60.0,0.617537*pow(plab,-1.98735));
	sigy[1]=min(60.0,0.18475*pow(plab,-3.06927));
	sigy[2]=min(60.0,0.172513*pow(plab,-2.83577));
	sigy[3]=min(60.0,0.23*pow(plab,-3.32));

    }

}

//***********************************************************************

//      subroutine jamcmbs2(msel,srt,pr,istr,kfb,kfm,kcb,kcm,izb,izm,
//     $ em1,em2,kf1,kf2,sig,sigel,sigin,
//     $ mabsrb,mchanel,mxchan,ijet,icon)

void SigmaMB::calcS2(CollisionPair& cpair)
{
//...Strangeness S=-2 meson-baryon cross section. xi + pi, k- + sigma

    int izm=cpair.getZ(0)/3;
    int izb=cpair.getZ(1)/3;
    int izt=izm+izb;
    double srt=cpair.getCMenergy();
    double pr=cpair.getCMmomentum();
    int kfm=cpair.getID(0);
    int kfb=cpair.getID(1);
    int istr=cpair.getTotalStrangeness();
    bool absorb=false;
    if(izt>= -1 && izt <=0 ) absorb=true;
    double srt0=2.2, sres=6.0;

    // If meson is in the formation time.
    //if(optSChannel==1 && cpair.qFactor(0) < 1.0) absorb=false;
    if(optSChannel==1 && (cpair.qFactor(0) < 1.0 || cpair.qFactor(1)<1.0))
	absorb=false;

    //mSel=1;
    mSel=cpair.mode();
    //double parc62=2.8;
//....Get total cross section.
    double sigr=0.0;
    if(mSel == 1) {
	if(absorb && srt < sres) {
          sigr=jamxbw1(istr,srt,pr,izm,izb,kfm,kfb,cpair.isAnti());
	}
	cpair.setSigAbs(sigr);
	double siga=0.0,sigel=0.0;
	//if(srt > 4.0)
	XsecTable::additiveQuarkModel(kfm,kfb,siga,sigel);
	sigTot=siga+sigr;
	sigEl=sigel;
	cpair.setXS(sigTot,sigEl,siga-sigel,sigr);
        return;
    }

    sigRand=cpair.getXsig();
    if(absorb && srt < sres) {

	// Monte-Calro.
      double sig=jamxbw1(istr,srt,pr,izm,izb,kfm,kfb,cpair.isAnti());
      if(sigRand <= 0.0) {
	cpair.setOutGoing(sig,kfm);
	cpair.setOutGoingMass(srt);
	return;
      }

      double sigr0=cpair.getSigAbs();
      double sigin_s = 0.0;
      if(srt > srt0) sigin_s = sigr0*srt0/srt;  // s-channel string formation

	//....Monte Carlo for s-channel string formation.
	sigRand -= sigin_s;
	if(sigRand<=0.0) {
	  int kf;
          if(izt == -1) kf=3312;
          if(izt ==  0) kf=3322;
	  if(cpair.isAnti()) kf *= -1;
	  cpair.setOutGoing(sigin_s,kf);
	  cpair.setOutGoingMass(srt);
          cpair.setChannel(-2);
	  return;
    }
    }

    // Monte-Carlo for t-channel contribution.
    double sigtr=cpair.getSigInT()*3.5/srt;
    sigRand -= sigtr;

    // t-channel string formation.
    //if(srt > 3.5 && sigRand < 0.0) {
    //if(srt > 3.5) {
    if(srt > ecmStringMBs) {
      cpair.setOutGoing(sigtr,92,92);
      return;
    }

  if(srt < 2.2) {
    cpair.setChannel(-1);
    return;
  }

    // t-channel resonance formation.
    double m1=cpair.M(0), m2=cpair.M(1);
    //pd1 = cpair.getParticleDataEntry(0);
    //pd2 = cpair.getParticleDataEntry(1);
    pd1=jamtable->find(kfm);
    pd2=jamtable->find(kfb);
    int idm=cpair.getPID(0);
    int idb=cpair.getPID(1);
  bool preHadronA = cpair.qFactor(0) < 1.0 ? true: false;
  bool preHadronB = cpair.qFactor(1) < 1.0 ? true: false;
    int icon=sampleMass->jamrmas2(cpair.getCollType(),pd1,pd2,kfm,kfb,
		idm,idb,izm,izb, m1,m2,preHadronA,preHadronB,srt,cpair.isAnti());

    if(icon==1) {
        cpair.setChannel(-1);
        return;
    }
    cpair.setOutGoing(sigtr,sampleMass->id(0),sampleMass->id(1));
    cpair.setOutGoingMass(sampleMass->m(0),sampleMass->m(1));
}

//***********************************************************************

double SigmaMB::jamxbw1(int istr,double srt,double pr, int iz1,int iz2,
	int& kfm,int& kfb,bool anti)
{
//...Purpose: to calculate Breit-Wigner cross section for M-B.
// -------------------------                                            *
//  meson-baryon absorptions for no total strangeness                   *
// --------------------------                                           *
//                                                                      *
//    pion   + n        ==> d(1232)/n*/d*                               *
//    pion   + d(1232)  ==> d(1232)/n*/d*                               *
//    pion   + n*       ==> d*                                          *
//    rho    + n        ==> n*/d*                                       *
//    eta    + n        ==> n*                                          *
//    eta'   + n        ==> n*                                          *
//    omega  + n        ==> n*                                          *
//    kaon   + lambda   ==> n*                                          *
//    kaon   + sigma    ==> n*/d*                                       *
//                                                                      *
// -------------------------                                            *
//  meson-baryon absorptions for S=-1 strangeness                       *
// --------------------------                                           *
//                                                                      *
//    Kbar   + n        ==> L*/S*                                       *
//    pion   + Sigma    ==> Lambda*/Sigma*                              *
//    .....etc.                                                         *
//                                                                      *
// iz1,iz2: charge
//                                                                    *
//--------------------------------------------------------------------*

    int izt=iz1+iz2;

//============================================================
//...n + pi => d(1232)
//============================================================

//...Get N + pi => D(1232) cross section.
    double sig=0.0;
    ParticleTable *table;
    sigElBW=0.0;

    if(istr == 0) {


//============================================================
//  (1)  Delta* formation cross sections,baryon + meson => d*
//============================================================
//...offset: 1:n* 2:p*  3:d*- 4:d*0 5:d*+ 6:d*++

	//ParticleTable *table=jamtable->getDelta();
	if(izt==-1) table=jamtable->getDmstar();
	else if(izt==0) table=jamtable->getD0star();
	else if(izt==1) table=jamtable->getDpstar();
	else if(izt==2) table=jamtable->getDppstar();
	else {
	    cout << "SigmaMB::jamxbw1 Delta* wrong charge izt= " << izt
		<< endl;
	    exit(1);
	}

	sig = BreitWigner(table,0,srt,pr,izt,kfm,kfb,anti);

	//cout << " iz1= "<< iz1 << " iz2= "<< iz2 << " izt= "<< izt
	//    << " sig= "<< sig <<endl;
	//cin.get();

//============================================================
// (2) Get N* formation cross sections,   n + pi => n*
//============================================================
	if(izt==0) table=jamtable->getNstar();
	else if(izt==1) table=jamtable->getPstar();
	else return sig;

	sig += BreitWigner(table,0,srt,pr,izt,kfm,kfb,anti);

	return sig;

    } else if(istr == -1) {

//============================================================
//  (3)  Sigma* formation cross sections,baryon + meson => d*
//============================================================
        if(izt < -1 || izt > 1) return 0.0;

//...Sigma*
	//ParticleTable *table=jamtable->getSigma();
	if(izt==-1) table=jamtable->getSmstar();
	else if(izt==0) table=jamtable->getS0star();
	else if(izt==1) table=jamtable->getSpstar();
	else {
	    cout << "SigmaMB::jamxbw1 Sigma* wrong charge izt= " << izt
		<< endl;
	    exit(1);
	}
	sig = BreitWigner(table,0,srt,pr,izt,kfm,kfb,anti);

        if(izt != 0) return sig;

//============================================================
//  (4) Get L* formation cross sections
//============================================================
	table=jamtable->getLambda();
	sig += BreitWigner(table,0,srt,pr,izt,kfm,kfb,anti);

    } else if(istr == -2) {
//============================================================
//  (5) Get Xi* formation cross sections
//============================================================
        //if(izt < -1 || izt > 0) return 0.0;
	//ParticleTable *table=jamtable->getXi();

	if(izt==-1) table=jamtable->getXmstar();
	else if(izt==0) table=jamtable->getX0star();
	else return 0.0;

	sig = BreitWigner(table,0,srt,pr,izt,kfm,kfb,anti);

    } else {
        cout << "(jamxbw1:)Invalid strangness istr" << istr << endl;
	exit(1);
    }

    return sig;
}

double SigmaMB::BreitWigner(ParticleTable* table,int i0,double srt, double pr,
	int iz, int& kf1, int& kf2, bool anti)
{
    double sig=0.0;
    int spin1=max(1,abs(kf1)%10);
    int spin2=max(1,abs(kf2)%10);
    static const double hbc=0.197327;
    double factbw=4*M_PI*10.0*hbc*hbc/(pr*pr)/(spin1*spin2);

    int n=table->size();
    for(int i=i0;i<n;i++) {
	ParticleDataEntry* p=table->getParticle(i);
	if(!p->canDecay()) continue;
	if(p->mWidth() < 1e-10) continue;
	if(p->chargeType() != 3*iz) continue;
	//double totwid=decay->getTotalWidth(p,srt,kf1,kf2);        

	if(srt < p->mMin()) continue;
	double pwid=decay->getPartialWidth(p,srt,kf1,kf2);        
	if(pwid <1e-15) continue;

	//cout << "kf0= "<< p->name()
	//    << " id= "<< p->id()
	//    << " m= "<< srt
	//    << " kf1= " << kf1 << " kf2= "<< kf2 << " pwid= "<< pwid <<endl;
	double totwid=jamtable->totalWidth(p,srt);        

	//int ip=decay->getIP();
	//if(ip == -1) continue;

        double  emr=p->m0();
        int     kfr=p->id();
	int ispin = p->spinType();
        //double  pwid=decay->getPWidth(ip);

	// Non-relativistic
        //double  bw=factbw*ispin*pwid*0.25
        //       /( (srt-emr)*(srt-emr) + 0.25*totwid*totwid );

	// relativistic I
	//double ss = srt*srt;
        //double  bw=factbw*ispin*pwid*srt*emr
        //       /( (ss-emr*emr)*(ss-emr*emr) + emr*emr*totwid*totwid );

	// relativistic II
	double ss = srt*srt;
        double  bw=factbw*ispin*pwid*ss /
	    ( pow2(ss-emr*emr) + ss*totwid*totwid );

	double sigbw   = bw*totwid;
	double sigbwel = bw*pwid;

	//cout << " candecay= " << p->canDecay() << " id= " << p->id() 
	//     << " sig= " << sigbw << " sigel= "<< sigbwel << endl;

	/*
	cout << " candecay= " << p->canDecay()
	     << " id= " << p->id()
	     << " w= " << p->mWidth()
	     << " iz= " << iz << " z= "<< p->chargeType()
	     << " wid= " << totwid
	     << " kf1= " << kf1 << " kf2= " << kf2
	     << " ip= " << ip
	     << " sig= " << sigbw
	     <<endl;
	     */

        if(mSel == 3) {
            sigRand -= sigbw;
            if(sigRand <= 0.0) {
              kf1=kfr;
              kf2=0;
	      pRes=p;
              if(anti && p->hasAnti()) kf1 *= -1;
              return sigbw;
	    }
	}
        sig += sigbw;
	sigElBW += sigbwel;

    }
    
    return sig;
}

} // end of namespace jam2
