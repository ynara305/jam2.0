#include <iostream>
#include <cstdlib>
#include <jam2/xsection/SigmaBB.h>
#include <jam2/xsection/XsecTable.h>
#include <jam2/hadrons/ParticleTable.h>
#include <jam2/hadrons/JamStdlib.h>

using namespace std;

namespace jam2 {


//const double SigmaBB::sMin=2.0139999; // thd0 in sub. jamxnnin
const double SigmaBB::sMin=2.018;
const double SigmaBB::Mn=0.9396; // 0.93957  parc(24)
const double SigmaBB::Mp=0.9383; // 0.93827  parc(25)
const double SigmaBB::Mnucl=0.93895; //  parc(28)
const double SigmaBB::Mdelta=1.232;
const double SigmaBB::Mpion=0.14;

// maximum allowed resonance mass.
const double SigmaBB::maxRMass= 3.0;

// normalization of the BW for Delta(1232)
const double SigmaBB::normD1232= 1.163098773625961;


const double maxSigmaNN=55.0;

//const int id_nucl  = 11;
//const int id_nucls = 12;
//const int id_delt  = 13;
//const int id_delts = 14;

const int    SigmaBB::mxhlf=30;

// outgoing resonance state for pp collision
const int SigmaBB::BBpp[16][2] = {
  {2212, 2214},   {2112, 2224},     // 1:pd+ 2:nd++
  {2212, Pstar},                  // 3:pp*
  {2214, 2214}, {2114, 2224},     // 4 d+d+ 5 d0d++
  {2212, Dstarp}, {2112, Dstarpp},  // 6 pd*+ 7 nd*++
    //Pstar,2214,   Nstar,2224,     // 8 p*d+ 9 n*d++
    //Dstarp,2214,  Dstar0,2224,    // 10 d*+d+ 11 d*0d++
  {2214,Pstar},   {2224,Nstar},     // 8 p*d+ 9 n*d++
  {2214,Dstarp},  {2224,Dstar0},    // 10 d*+d+ 11 d*0d++
  {Pstar,Pstar},                  // 12 p*p*
  {Pstar,Dstarp},{Nstar,Dstarpp},   // 13 p*d*+ 14 n*d*++
  {Dstarp,Dstarp},{Dstar0,Dstarpp}  // 15 d*d*+ 16 d*d*++
    };

// outgoing resonance state for nn collision
const int SigmaBB::BBnn[16][2] = {
    2112, 2114,   2212, 1114,    // 1:nd0 2:pd-
    2112, Nstar,                 // 3:nn*
    2114, 2114,   2214, 1114,    // 4 d0d0 5 d+d-
    2112, Dstar0, 2212, Dstarm,  // 6,nd*0 7 pd*-
    //Nstar,2114,   Pstar,1114,    // 8 n*d0 9 p*d-
    //Dstar0,2114,  Dstarm,2214,   // 10 d0d*0 11 d+0d*-
    2114,Nstar,   1114,Pstar,    // 8 n*d0 9 p*d-
    2114,Dstar0,  2214,Dstarm,   // 10 d0d*0 11 d+0d*-
    Nstar,Nstar,                 // 12 n*n*
    Nstar,Dstar0,Pstar,Dstarm,   // 13 n*d*0 14 p*d*-
    Dstar0,Dstar0,Dstarp,Dstarm  // 15 d*0d*0 15 d*+d*-
    }; 

// outgoing resonance state for pn collision
const int SigmaBB::BBpn[17][2] =  {
    2112, 2214,   2212, 2114,     // 1:d+n 2 d0p
    2112, Pstar,  2212, Nstar,    // 3:np* 4 n*p
    2114, 2214,   1114, 2224,     // 5 d0d+ 6 d-d++
    2112, Dstarp, 2212, Dstar0,   // 7,nd*+ 8 pd*0
    //Nstar,2214,   Pstar,2114,     // 9 n*d+ 10 p*d0
    //Dstarp,2114,  Dstarpp,1114,   // 11 d0d*+ 12 d-0d*++
    2214,Nstar,   2114,Pstar,     // 9 n*d+ 10 p*d0
    2114,Dstarp,  1114,Dstarpp,   // 11 d0d*+ 12 d-0d*++
    Nstar,Pstar,                  // 13 n*p*
    Nstar,Dstarp,Pstar,Dstar0,    // 14 n*d*+ 15 p*d*0
    Dstar0,Dstarp,Dstarm,Dstarpp  // 16 d*0d*+ 17 d*-d*++
    };

//...1) NN -> ND  
//...2) NN -> NN*
//...3) NN -> DD
//...4) NN -> ND*
//...5) NN -> N*D
//...6) NN -> DD*
//...7) NN -> N*N*
//...8) NN -> N*D*
//...9) NN -> D*D*

// -2. D-D-
const int SigmaBB::DDm2[9][2]={
    0,0, 0,0, 1114,1114, 0,0, 0,0, 1114,Dstarm, 0,0,0,0, Dstarm,Dstarm};
// 4. D++D++
const int SigmaBB::DD4[9][2]={
    0,0, 0,0, 2224,2224, 0,0, 0,0, 2224,Dstarpp, 0,0, 0,0, Dstarpp,Dstarpp};
// -1 D0D-
const int SigmaBB::DDm1[9][2]={
    2112,1114,
    0,0,
    2114,1114,
    2112,Dstarm,
    1114,Nstar,
    2114,Dstarm,
    0,0,
    Nstar,Dstarm,
    Dstar0,Dstarm};
// 3 D+D++
const int SigmaBB::DD3[9][2]={
    2212,2224,
    0,0,
    2214,2224,
    2212,Dstarpp,
    2224,Pstar,
    2214,Dstarpp,
    0,0,
    Pstar,Dstarpp,
    Dstarp,Dstarpp};

// -1. n D-
const int SigmaBB::NDm1[9][2]={
  2112, 1114, 
  0,   0,
  2114, 1114,
  2112, Dstarm,
  1114, Nstar,
  2114, Dstarm,
  0,    0,
  Nstar, Dstarm,
  Dstar0, Dstarm};

// 3. p D++
const int SigmaBB::ND3[9][2]={
  2212, 2224,
  0,   0,
  2214, 2224,
  2212, Dstarpp,
  2224, Pstar,
  2214, Dstarpp,
  0,    0,
  Pstar, Dstarpp,
  Dstarp, Dstarpp};


SigmaBB::SigmaBB(Pythia8::Info *inf,Pythia8::Settings *s,JamParticleData* jp, SampleMass *sm, Pythia8::Rndm* r)
{
    info = inf;
    settings = s;
    jamTable=jp;
    decay=jamTable->getDecayWidth();
    sampleMass = sm;
    rndm=r;
    sigTotPy8 = new Pythia8::SigmaTotal();
    sigTotPy8->init(info,*settings,jamTable->getParticleData(),rndm);
    absorptionBB=settings->flag("Cascade:BBabsorptionXS");

    nStar[0]=jamTable->getNstar();
    nStar[1]=jamTable->getPstar();
    dStar[0]=jamTable->getDmstar();
    dStar[1]=jamTable->getD0star();
    dStar[2]=jamTable->getDpstar();
    dStar[3]=jamTable->getDppstar();
    nMaxRes = 1;
    nMaxRes=max(nStar[0]->size(),nMaxRes);
    nMaxRes=max(nStar[1]->size(),nMaxRes);
    nMaxRes=max(dStar[0]->size(),nMaxRes);
    nMaxRes=max(dStar[1]->size(),nMaxRes);
    nMaxRes=max(dStar[2]->size(),nMaxRes);
    nMaxRes=max(dStar[3]->size(),nMaxRes);

    protonT.push_back(jamTable->getProton());
    neutronT.push_back(jamTable->getNeutron());
    deltamT.push_back(jamTable->getDeltam());
    delta0T.push_back(jamTable->getDelta0());
    deltapT.push_back(jamTable->getDeltap());
    deltappT.push_back(jamTable->getDeltapp());
    N1440T.push_back(jamTable->getN1440());
    P1440T.push_back(jamTable->getP1440());
    for(int jt=0;jt<2;jt++)
    for(int i=0;i<nStar[jt]->size();i++) {
      NStar[jt].push_back(nStar[jt]->getParticle(i));
    }
    for(int jt=0;jt<4;jt++)
    for(int i=0;i<dStar[jt]->size();i++) {
      int id=dStar[jt]->getParticle(i)->id();
      if(id==1114 || id==2114 || id==2214 || id==2224) continue;
      DStar[jt].push_back(dStar[jt]->getParticle(i));
    }

    /*
    cout << " nstar= " << nStar[0]->size() <<endl;
    cout << " nstar= " << nStar[1]->size() <<endl;
    cout << " dstar= " << dStar[0]->size() <<endl;
    cout << " dstar= " << dStar[1]->size() <<endl;
    cout << " dstar= " << dStar[2]->size() <<endl;
    cout << " dstar= " << dStar[3]->size() <<endl;
    if(nMaxRes==1) {
	cout  << " nMaxRes=1? " <<nMaxRes<<endl;
	exit(1);
    }
    */

}

//***********************************************************************
//...pp/pn total and elastic cross sections.
void SigmaBB::jamxnn(double ecm,int z1, int z2)
{
//....nn/pp collision.

    if(z1 == z2) {

        //if(ecm < 4.5) {
        if(ecm < 5.0) {
          sigTot=XsecTable::jamsighh(1,ecm);
          if(ecm <= 2.4) {
            sigEl=XsecTable::nnElastic(ecm,2);
	  } else {
            sigEl=XsecTable::jamsighh(2,ecm);
	  }

            //sigEl=XsecTable::nnElastic(ecm,2);

         } else {

          double plab;
          if(z1 == 3) {
	    plab=XsecTable::plabsr(ecm,Mp,Mp);
	  } else {
	    plab=XsecTable::plabsr(ecm,Mn,Mn);
	  }
          //sigTot=XsecTable::jamchc96(16,plab);
          //sigEl=XsecTable::jamchc96(17,plab);

          sigTot=XsecTable::pdg2016(ecm,-1);

	  if(ecm < 600) {
            sigEl=XsecTable::jamchc96(17,plab);
	  } else {
	    const double bhad=2.3, eps=0.0808,facel=0.0511;
	    double bel=2.0*bhad+2.0*bhad + 4.0*pow(ecm*ecm,eps) - 4.2;
            sigEl=1.2*facel*sigTot*sigTot/bel;
	  }


	 }

//...np collision.
    } else {

        if(ecm <= 2.4) {
          sigEl=XsecTable::nnElastic(ecm,1);
	} else if(ecm < 4.5) {
          sigEl = XsecTable::jamsighh(4,ecm);
	} else if(ecm < 600) {
	  double plab=XsecTable::plabsr(ecm,Mp,Mn);
          sigEl=XsecTable::jamchc96(17,plab);
	} else {
          sigTot=XsecTable::pdg2016(ecm,-2);
	  const double bhad=2.3, eps=0.0808,facel=0.0511;
	  double bel=2.0*bhad+2.0*bhad + 4.0*pow(ecm*ecm,eps) - 4.2;
          sigEl=1.2*facel*sigTot*sigTot/bel;
	}

        if(ecm < 5.0) {
          sigTot = XsecTable::jamsighh(3,ecm);
	} else if( ecm < 6.5) {
	  sigTot = 39.2;
	} else {
          //sigTot=XsecTable::jamchc96(18,plab);
          //sigTot=XsecTable::jamchc96(18,plab);
	
          sigTot=XsecTable::pdg2016(ecm,-2);

	}
    }


    //sigTot = min(maxSigmaNN,sigTot);
    //if(ecm <= sMin+eKinMin) sigTot=min(maxSigmaNN,sigEl);
    //sigEl = min(sigTot,sigEl);
}

vector<ParticleDataEntry*> SigmaBB::pickRTable(int idn, int& idr)
{
  if(idn==2212)  {
    idr=id_nucl;
    return protonT;
  } else if (idn==2112) {
    idr=id_nucl;
    return neutronT;
  } else if(idn==1114) {
    idr=id_delt;
    return deltamT;
  } else if(idn==2114) {
    idr=id_delt;
    return delta0T;
  } else if(idn==2214) {
    idr=id_delt;
    return deltapT;
  } else if (idn==2224) {
    idr=id_delt;
    return deltappT;
  } else if (idn==12112) {
    idr=id_nucls;
    return N1440T;
  } else if (idn==12212) {
    idr=id_nucls;
    return P1440T;
  } else {
    if(idn==1) {idr=id_nucls; return NStar[0];}
    if(idn==2) {idr=id_nucls; return NStar[1];}
    if(idn==3) {idr=id_delts;return DStar[0];}
    if(idn==4) {idr=id_delts;return DStar[1];}
    if(idn==5) {idr=id_delts;return DStar[2];}
    if(idn==6) {idr=id_delts;return DStar[3];}
  }

  cout << " SigmaBB::pickRTable wrong id = " << idn <<endl;
  exit(1);
}

bool SigmaBB::sampleMassFix(double ecm,double pf0,int iex[2],double emr[2],double em0[2],
	double gam0[2],double mmin[2],double ymin[2],double ymax[2],
	double& m1, double& m2)
{
  int itry=0;
  double pf=0.0;
  do {
    if(++itry > 100) {
      cout << "SigmaBB::samplemass does not conserve ecm= " << ecm;
      cout << " mmin1+mmin2+ekin= "<< mmin[0]+mmin[1]+eKinMin;
      cout << " emr1= " << emr[0] << " emr2 " << emr[1];
      cout << " mmin1= " << mmin[0] << " mmin2= " << mmin[1] <<endl;
      return false;
    }
    for(int jt=0;jt<2;jt++) {
      if(iex[jt]==0) continue;
      double y=ymin[jt] + rndm->flat()*( ymax[jt] - ymin[jt] );
      emr[jt] = sqrt( em0[jt]*gam0[jt]*tan(y) + em0[jt]*em0[jt] );
    }
    if(emr[0]+emr[1]+eKinMin > ecm) continue;
    pf=PCM(ecm,emr[0],emr[1]);

  } while (rndm->flat()*pf0 > pf);

  m1=emr[0]; m2=emr[1];
  return true;
}

void SigmaBB::makePP(double* sig1)
{
	sigin[0]=0.5*sig1[0];
	sigin[1]=0.5*sig1[0];

	sigin[2]=sig1[1];
        //sigin[1] += sig1[9]; // include s-wave pion into pp*
	for(int i=2;i<6;i++)  {
	    sigin[2*i-1]=0.5*sig1[i];
	    sigin[2*i]=0.5*sig1[i];
	}
	sigin[11]=sig1[6];
	sigin[12]=0.5*sig1[7];
	sigin[13]=0.5*sig1[7];
	sigin[14]=0.5*sig1[8];
	sigin[15]=0.5*sig1[8];

        sigin[16] = sig1[9]; // include s-wave pion into pp(1440)*

}

void SigmaBB::makePN(double* sig1)
{
	for(int i=0;i<9;i++) {
	    if(i<=5) {
	      sigin[2*i]=0.5*sig1[i];
	      sigin[2*i+1]=0.5*sig1[i];
	    } else if(i==6) {
	      sigin[12]=sig1[i];
	    } else {
	      sigin[2*i-1]=0.5*sig1[i];
	      sigin[2*i]=0.5*sig1[i];
	    }
	}
        sigin[17] = 0.5*sig1[9]; // include s-wave pion into pp(1440)*
        sigin[18] = 0.5*sig1[9]; // include s-wave pion into pp(1440)*

}

//***********************************************************************

//      subroutine jamcbbs(msel,kf1,kf2,kc1,kc2,istr,srt,
//     $ sig,sigel,sigin,em1,em2,mchanel,mxchan,icon)

void SigmaBB::calcS(CollisionPair& cpair)
{
  int str=cpair.getTotalStrangeness();
  double srt = cpair.getCMenergy();
  //...Initialize some values.
  id1 = cpair.getPID(0);
  id2 = cpair.getPID(1);
  kf1 = cpair.getID(0);
  kf2 = cpair.getID(1);
  int ipair=cpair.getPairID();
  int iz1=cpair.getZ(0)/3;
  int iz2=cpair.getZ(1)/3;
  double thr=0.0;

//....Handle low energy Lambda-N/Sigma-N/Xi-N/LL... collisions.
  if(str == -1 && srt <= 2.50) {
    // S=-1 BB collisions.
    thr = sigmaS1(ipair,srt,iz1,iz2);

  // S=-2 BB collisions.
  } else if(str == -2 && srt <= 2.70) {
    thr = sigmaS2(ipair,srt,iz1,iz2);

  } else {
    XsecTable::additiveQuarkModel(kf1,kf2,sigTot,sigEl);
    if(str == -1) {
      if((id1==id_lambda || id1==id_lambdas) || (id2==id_lambda || id2==id_lambdas)) {
        thr = 1.11568+Mnucl+Mpion+0.003;
      } else {
        thr = 1.2 + Mnucl + Mpion + 0.003;
      }
    } else if(str == -2) {
      if((id1==id_xi || id1==id_xis) || (id2==id_xi || id2==id_xis)) {
       	thr=1.33 + Mnucl+Mpion+0.01;
      } else {
       	thr=2*1.2 + Mnucl+Mpion+0.01;
      }
    } else if(str == -3) thr = 1.67245 + 1.232 + 0.003;
    else if(str == -4) thr = 1.67245 + 1.2;
    else if(str == -5) thr = 1.67245 + 1.33;
    else thr = 2*1.67245 + 0.3;
  }

  cpair.setXSTotal(sigTot);

  /*
    double mm = cpair.M(0)+cpair.M(1);
    double thr1= str==-1 ? mm+0.14: mm+0.3;
    cout << "S= "<< str << " id1= "<< kf1 << " id2= "<< kf2 << " srt= "<< srt << " thr= "<< thr << " thr1= "<< thr1
             << " m1= " <<  cpair.M(0)
             << " m2= " <<  cpair.M(1)
             << " " << jamTable->find(kf1)->name()
             << " " << jamTable->find(kf2)->name()
	     <<endl;
    cin.get();
    */

  // elastic only.
  if(srt < thr) cpair.setXSElastic(sigTot);

}

double SigmaBB::sigmaS1(int ipair,double srt,int iz01,int iz02)
{
    int iz1=iz01;
    int iz2=iz02;
    if(id1 == id_nucl) {
        iz1=iz02;
        iz2=iz01;
    }
    int izt=iz1+iz2;
    int isigt=0;
    int isige=0;
    double thr=0.0;

//....Lambda + nucleon.
    if(ipair == CollisionPair::Pair(id_lambda,id_nucl)) {

	if(izt == 1) {           // Lambda p
            isigt=1;
            isige=2;
	} else if(izt == 0) {      // Lambda n
            isigt=3;
            isige=4;
	}
     thr = 1.11568+Mnucl+Mpion+eKinMin;

//...Sigma + nucleon.
    } else if(ipair == CollisionPair::Pair(id_sigma,id_nucl)) {
      thr =  1.2 + Mnucl+Mpion + eKinMin;

        if(izt == -1) {           // Sigma- n
            isigt=5;
            isige=0;
	    //thr =  1.19745 + 1.232 + 0.003;
	} else if(izt == 0) {
          if(iz1 == 0) { // Sigma0 n
              isigt=9;
              isige=10;
	      //thr = 1.19264 + 1.232 + 0.003;
	  } else if(iz1 == -1) {  // Sigma- p
              isigt=7;
              isige=8;
	      //thr = 1.19745 + 1.232 + 0.003;
	  }
	} else if(izt == 1) {
          if(iz1 == 1) {         // Sigma+ n
              isigt=13;
              isige=14;
	      //thr = 1.18937 + 1.232 + 0.003;
	  } else if(iz1 == 0) {      // Sigma0 p
              isigt=11;
              isige=12;
	      //thr = 1.19264 + 1.232 + 0.003;
	  }

	} else if(izt == 2) {       // Sigma+ p
            isigt=6;
            isige=0;
	    //thr = 1.18937 + 1.232 + 0.003;
	}

    } else {
      thr=1.2+Mnucl+Mpion+eKinMin;
      XsecTable::additiveQuarkModel(kf1,kf2,sigTot,sigEl);
      return thr;
    }


//...Get total and elastic cross sections.
    //if(isigt != 0) sigTot=XsecTable::sigBBS1(isigt,srt);
    sigTot=XsecTable::sigBBS1(isigt,srt);
    if(isige >= 1) sigEl=XsecTable::sigBBS1(isige,srt);
    else sigEl=sigTot;

    return thr;

}

double SigmaBB::sigmaS2(int ipair,double srt,int iz1,int iz2)
{
  //int isw=0;
  //if(id1 == id_nucl || id1 == id_lambda) isw=1;

  int isigt=0;
  int isige=0;
  int izt=iz1+iz2;
  double thr=0.0;

//...Lambda Lambda ingoing.
  if(ipair == CollisionPair::Pair(id_lambda,id_lambda)) {

	isigt=9;
        isige=10;
	//thr = 2*1.11568 + Mpion + eKinMin;
	thr = 2.3083 + eKinMin;

//.....Lambda Sigma ingoing.
  } else if(ipair == CollisionPair::Pair(id_lambda,id_sigma)) {

        if(izt == -1) {     // Lambda Sigma-
              isigt=11;
              isige=12;
	} else if(izt == 0) { // Lambda Sigma0
              isigt=13;
              isige=14;
	} else if(izt == 1) { //! Lambda Sigma+
              isigt=15;
              isige=16;
	}
    //thr = 1.11568 + 1.2 + Mpion + eKinMin;
    thr = 2.54 + eKinMin;

//....Sigma Sigma ingoing.
  } else if(ipair == CollisionPair::Pair(id_sigma,id_sigma)) {
    //thr = 2*1.2 + Mpion + 0.01;
    thr = 2.54 + eKinMin;

          if(izt == -2) {               // Sigma- Sigma-
              isigt=17;
              isige=0;
	  } else if(izt == -1) {          // Sigma- Sigma0
              isigt=18;
              isige=19;
	  } else if(izt == 0) { 
            if(iz1 == 0 && iz2 == 0) { // Sigma0 Sigma0
                isigt=20;
                isige=21;
	    } else {                        // Sigma+ Sigma-
                isigt=22;
                isige=23;
	    }
	  } else if(izt == 1) {           // Sigma+ Sigma0
              isigt=24;
              isige=25;
	  } else if(izt == 2) {           // Sigma+ Sigma+
              isigt=26;
              isige=0;
	  }

//....Xi N incoming.
  } else if(ipair == CollisionPair::Pair(id_xi,id_nucl)) {
    thr = 1.33 + Mnucl + Mpion + 0.01;
    thr = 2.4 + eKinMin;
// p1= Xi0 p2= Delta0 id1= 3322 id2= 2114 m1= 1.3149 m2= 1.0800 m1+m2= 2.3949 prob= 0.0000

          if(izt == -1) {               // Xi- n
              isigt=1;
              isige=2;
	  } else if(izt == 0) {
            if(iz1 == 0 && iz2 == 0) {  //  Xi0 n
                isigt=5;
                isige=6;
	    } else {                          //  Xi- p
                isigt=3;
                isige=4;
	    }
	  } else if(izt == 1) {            // Xi0 p
              isigt=7;
              isige=8;
	  }
  } else {
      //cout << "SigmaBB::sigmaS2 unknown particles id1= " <<  id1
	//<< " id2= "<< id2
        //<< " " << jamTable->find(kf1)->name()
        //<< " " << jamTable->find(kf2)->name()
	//<< endl;
      XsecTable::additiveQuarkModel(kf1,kf2,sigTot,sigEl);
      //return thr;
      if((id1==id_xi || id1==id_xis) || (id2==id_xi || id2==id_xis))
       	thr=1.33 + Mnucl+Mpion+0.003;
      else
       	thr=2*1.2 + Mnucl+Mpion+0.003;
      return thr;
  }

//...Get total and elastic cross sections.
  //if(isigt>0) sigTot=XsecTable::sigBBS2(isigt,srt);
  sigTot=XsecTable::sigBBS2(isigt,srt);
  if(isige >= 1) sigEl=XsecTable::sigBBS2(isige,srt);
  else sigEl=sigTot;
  return thr;

}

} // end of namespace jam2
