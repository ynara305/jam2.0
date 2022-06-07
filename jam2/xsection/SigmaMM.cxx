#include <cmath>
#include <jam2/xsection/SigmaMM.h>
#include <jam2/xsection/XsecTable.h>

using namespace std;

namespace jam2 {

SigmaMM::SigmaMM(Pythia8::Settings *s,JamParticleData* table, SampleMass* sm,
	    Pythia8::Rndm* r)
{
  settings=s;
  jamtable=table;
  //decay=dc;
  decay=jamtable->getDecayWidth();
  particleData = jamtable->getParticleData();
  mesonTable = jamtable->getMeson();
  sampleMass = sm;
  ecmStringMM=settings->parm("Cascade:ecmStringMM");
  rndm=r;
  hadronContent = new HadronContent(particleData,rndm);
  optSChannel=settings->mode("Cascade:optConstQSChannel");
}

void SigmaMM::calc(CollisionPair& cpair)
{
  //pd1=cpair.getParticleDataEntry(0);
  //pd2=cpair.getParticleDataEntry(1);

  int iz1=cpair.getZ(0);
  int iz2=cpair.getZ(1);
  double srt=cpair.getCMenergy();
  double pr=cpair.getCMmomentum();
  int kf1=cpair.getID(0);
  int kf2=cpair.getID(1);
  int id1=cpair.getPID(0);
  int id2=cpair.getPID(1);
  bool preHadronA = cpair.qFactor(0) < 1.0 ? true: false;
  bool preHadronB = cpair.qFactor(1) < 1.0 ? true: false;

  pd1=jamtable->find(kf1);
  pd2=jamtable->find(kf2);

  bool absorb=true;
  // If meson is in the formation time.
  if(optSChannel==1 && (preHadronA || preHadronB))
    absorb=false;

  // Absorption impossible.
  int izt=iz1+iz2;
  if(izt <= -6 || izt >= 6) absorb=false;

  mChanel=1;
  mAbsrb=2;
  double sigres=0.0;
  sigRand=cpair.getXsig();
  mSel=cpair.mode();

//.....Low energy pion - pion collisions.
//       id1=kchg(kc1,5)
//       id2=kchg(kc2,5)
//       ipair=jamcpair(id1,id2)
//       if(ipair.eq.jamcpair(id_pi,id_pi).and.srt.le.0.9) then
//          call jamcpipi(msel,srt,kf1,kf2,iz1,iz2,
//    $      em1,em2,sig,sigel,sigab,mchanel,mabsrb,ijet)
//          if(msel.eq.3.and.pare(3).le.0.0d0) return
//       endif

  int ifla1, iflb1, ifla2,iflb2;
  hadronContent->findFlavor(kf1,ifla1,iflb1);
  hadronContent->findFlavor(kf2,ifla2,iflb2);
  double mth = (abs(ifla1) < 3) ? 0.14 : (abs(ifla1) == 3) ? 0.6 : 1.5;
  mth       += (abs(iflb1) < 3) ? 0.14 : (abs(iflb1) == 3) ? 0.6 : 1.5;
  mth       += (abs(ifla2) < 3) ? 0.14 : (abs(ifla2) == 3) ? 0.6 : 1.5;
  mth       += (abs(iflb2) < 3) ? 0.14 : (abs(iflb2) == 3) ? 0.6 : 1.5;

  int nc1 =(abs(ifla1) > 3 || abs(iflb1) > 3) ? 1 : 0;
  int nc2 =(abs(ifla2) > 3 || abs(iflb2) > 3) ? 1 : 0;
  if(nc1+nc2> 0) mth = 4.0;

  /*
  if((abs(kf1)==421 || abs(kf1)==413)  ||
     (abs(kf2)==421 || abs(kf2)==413) ) {
      cout << " id1= " << kf1 << " id2= "<< kf2
	   << " nc1= "<< nc1 << " nc2= "<< nc2 <<endl;
      cin.get();
  }
  */

  //double mth = 2*0.14 + 2*ns1 + 2*ns2;

  // Calculate resonance cross sections.
  //emRf=1.0 + ns1 + ns2;
  emRf=0.0;
  if(absorb && srt <= 3.0) sigres = jamxbw2(srt,pr,kf1,kf2,iz1,iz2,mSel);
  if(sigres == 0.0) absorb = false;
  double sig=sigres;
  double sigel=0.0;

  // s-channel resonance formation: m + m -> R
  if(mSel == 3 && sigRand <=0.0) {
    cpair.setOutGoing(sig,idRes);
    cpair.setOutGoingMass(srt);
    return;
  }

  // Additive quark cross section.
  double siga=0.0;
  if(srt >= emRf)  {
  //if(srt >= emRf || sigres== 0.0)  {
    XsecTable::additiveQuarkModel(kf1,kf2,siga,sigel);
    sig += siga;
    //if(sigres == 0.0 || srt < mth) sig = sigel;
    //if(srt < mth) sig = sigel;
  }

  cpair.setXSTotal(sig);
  cpair.setXSElastic(sigel);

  /*
  if(srt<1.5 && !absorb) {
  cout << " kf1= "<< kf1 << " kf2= "<< kf2 << " srt= "<< srt
      << " sig= "<< sig
      << " sigab= "<< sigres
      << " sigel= "<< sigel
      << " emRF= "<< emRf
      << " ns1= "<< ns1
      << " ns2= "<< ns2
      << endl;
  cout << " ifla1= "<< ifla1 << " iflb1= "<< iflb1 <<endl;
  cout << " ifla2= "<< ifla2 << " iflb2= "<< iflb2 <<endl;
  cin.get();
  }
  */

  if(mSel == 1) return;

  double sig3=0.0;
  if(absorb) {
  int kfla1=max(ifla1,iflb1);
  int kflb1=min(ifla1,iflb1);
  int kfla2=max(ifla2,iflb2);
  int kflb2=min(ifla2,iflb2);
  int iann=0;
  int iann1=0;
  int iann2=0;
  if(kfla1 == abs(kflb2)) iann1=1;
  if(abs(kflb1) == kfla2) iann2=1;
  if(iann1 == 1 && iann2 == 1) {
    iann=1;
    if(rndm->flat() > 0.5) iann=2;
  } else if(iann1 == 1 ) {
    iann=1;
  } else if(iann2 == 1) {
    iann=2;
  }

  double sth;
  int kfv1=0;
  if(iann == 1) {
    //kfcnst(kfla2,kflb1,kfv1,0.0);
    kfv1 = hadronContent->combine(kfla2, kflb1);
    sth=hadronContent->jamemjet(kfla2,kflb1);
  } else if(iann == 2) {
      //kfcnst(kfla1,kflb2,kfv1,0.0);
      kfv1 = hadronContent->combine(kfla1, kflb2);
      sth=hadronContent->jamemjet(kfla1,kflb2);
  }
  if(srt < sth) kfv1=0;

  // Monte Carlo for s-channel string formation.
  if(kfv1 != 0 && mSel == 3)  {
    double sigres1 = jamxbw2(srt,pr,kf1,kf2,iz1,iz2,1);
    sig3=sigres1*3.0/srt;
    sigRand -= sig3;
    if(sigRand <= 0.0) {
      cpair.setOutGoing(sig3,kfv1);
      cpair.setOutGoingMass(srt);
      cpair.setChannel(-2);
      return;
    }
  }

  } // absorb end

  // test
  //cpair.setChannel(-1);
  //return;

  // t-channel inelastic.
  if(mSel == 3) {
    //double sths=hadronContent->jamemjet(ifla1,iflb1)
    //           +hadronContent->jamemjet(ifla2,iflb2);

    double sths= max(ecmStringMM,hadronContent->jamemjet(ifla1,iflb1)
               +hadronContent->jamemjet(ifla2,iflb2));

    /*
    double s1= hadronContent->jamemjet(ifla1,iflb1);
    double s2= hadronContent->jamemjet(ifla2,iflb2);
    cout << "kf1= " <<  kf1 <<" kf2= "<< kf2
          << " s1= " << s1
          << " s2= " << s2
          << " s= " << s1+s2
	       <<endl;
    cin.get();
    */

    if(srt >= sths)  {
      // t-channel string formation.
      cpair.setOutGoing(sig,92,92);
    } else {

      if(srt < mth) {
        cpair.setChannel(-1);
        return;
      }

     // test
     //cpair.setChannel(-1);
     //return;

      double m1=cpair.M(0), m2=cpair.M(1);
      int icon=sampleMass->jamrmas2(cpair.getCollType(),pd1,pd2,kf1,kf2,id1,id2,iz1/3,iz2/3,
  	       m1,m2,preHadronA,preHadronB,srt);
      if(icon==1) {
        cpair.setChannel(-1);
	return;
      }
      cpair.setOutGoing(sig-sig3,sampleMass->id(0),sampleMass->id(1));
      cpair.setOutGoingMass(sampleMass->m(0),sampleMass->m(1));
      //cpair.setOutGoingParticle(sampleMass->p(0),sampleMass->p(1));
    }
  }


}

//***********************************************************************

double SigmaMM::jamxbw2(double srt,double pr,
	int kf1,int kf2,int iz1,int iz2,int isel)
{
//...Purpose: to calculate Berit-Wigner cross section for M-M.

    sigRand=10.0;
    int isg=1;
    //emRf=0.0;
    int izt=iz1+iz2;
    int kf01=kf1;
    int kf02=kf2;
    if(izt <= -6 || izt >= 6)  return 0.0;

    if(kf01 < 0 && pd2->hasAnti()==false) {
        kf01=-kf01;
        isg=-1;
        izt = 3;
    } else if(kf02 < 0 && pd1->hasAnti()==false) {
        kf02=-kf02;
        isg=-1;
        izt = 3;
    }
    double sig = BreitWigner(srt,pr,izt,kf01,kf02,isel);
    idRes *= isg;

    return sig;
}

double SigmaMM::BreitWigner(double srt, double pr,
	int iz, int kf1, int kf2,int isel)
{
    double xnorm= kf1 == kf2 ? 2.0 : 1.0;
    int kfa1=abs(kf1);
    int kfa2=abs(kf2);
    double spin=max(1,kfa1%10)*max(1,kfa2%10);
    static const double hbc=0.197327;
    double factbw=4*M_PI*10.0*hbc*hbc/(pr*pr)/spin;
    factbw *=xnorm;

    double sig=0.0;
    for(int i=0,n=mesonTable->size();i<n;++i) {
	Pythia8::ParticleDataEntry* p=mesonTable->getParticle(i);
	int kf=p->id();
	if(kf==130 || kf==310) continue;  // K0_L or K0_S
	if(!p->canDecay()) continue;
	if(p->mWidth() < 1e-10) continue;
	if(p->chargeType() != iz) continue;
        //if(srt.lt.pmas(ir,1)-pmas(ir,3)) goto 10

	//double totwid=decay->getTotalWidth(p,srt,kf1,kf2);        
	//int ip=decay->getIP();
	// This meson does not have decay mode to kf1 + kf2.
	//if(ip == -1) continue;
        //double  pwid=decay->getPWidth(ip);

	if(srt < p->mMin()) continue;
        double pwid=decay->getPartialWidth(p,srt,kf1,kf2);
	if(pwid <1e-15) continue;
        double totwid=jamtable->totalWidth(p,srt);


        double  emr=p->m0();
	emRf = max(emRf,emr+0.3);
        int     kfr=p->id();
	int kfra=abs(kfr);
        int  ispin=max(1,kfra%10);

	// Non-relativistic BW
        //double  bw=factbw*ispin*pwid*0.25
        //       /( (srt-emr)*(srt-emr) + 0.25*totwid*totwid );

        // relativistic I
        double ss = srt*srt;
        //double  bw=factbw*ispin*pwid*srt*emr
        //       /( (ss-emr*emr)*(ss-emr*emr) + emr*emr*totwid*totwid );

        // relativistic II
        double  bw=factbw*ispin*pwid*ss /
            ( pow2(ss-emr*emr) + ss*totwid*totwid );

	double sigbw = bw * totwid;
	// symmetry factor.
	if(kf1==kf2) sigbw *=2;

        if(isel == 3) {
	    sigRand -= sigbw;
	    if(sigRand <= 0.0) {
	      idRes=kfr;
              return sig;
	    }
	}
        sig += sigbw;

    }

    return sig;
}


} // end namespace jam2
