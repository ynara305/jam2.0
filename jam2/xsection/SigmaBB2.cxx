#include <iostream>
#include <cstdlib>
#include <jam2/xsection/SigmaBB2.h>
#include <jam2/xsection/XsecTable.h>
#include <jam2/hadrons/ParticleTable.h>
#include <jam2/hadrons/JamStdlib.h>
#include <jam2/hadrons/BaryonTable.h>

using namespace std;

namespace jam2 {

const double srtStr=3.0;   // threshold of string excitation.
//const double srtStr=3.2;   // threshold of string excitation.
//const double srtStr=2.7;   // threthold of string excitation.

SigmaBB2::SigmaBB2(Pythia8::Info *inf,Pythia8::Settings *s,JamParticleData* jp, SampleMass *sm,
	Pythia8::Rndm* r,string fname) : SigmaBB(inf,s,jp,sm,r)
{
  //optRR=1;   //=1: NN->ND, NN*, DD, DN*, DD*
  optRR=2;     //=2: NN->ND, NN*, DD, DN*, DD*, N*N*,, N*D*, D*D*

  // optDeltaWidth=1:RQMD  =2:GiBUU   =3:Randrup
  // optWid=1: RQMD =2:Giessen =3:GiBUU

  //int optDeltaWid=2, optWid=1, optBW=3;
  int optDeltaWid=2, optWid=3, optBW=3;
 
  sigmaRR = new SigmaRR(fname,optRR,optDeltaWid, optWid, optBW, jp);

}

SigmaBB2::~SigmaBB2()
{
  delete sigmaRR;
}

//***********************************************************************

void SigmaBB2::calc(CollisionPair& cpair)
{
    mMode=cpair.mode();
    srt = cpair.getCMenergy();
    pr  = cpair.getCMmomentum();
    em1=cpair.M(0);
    em2=cpair.M(1);
    iz1=cpair.getZ(0);
    iz2=cpair.getZ(1);
    kf1=cpair.getID(0);
    kf2=cpair.getID(1);
    id1=cpair.getPID(0);
    id2=cpair.getPID(1);
    //pd1=cpair.getParticleDataEntry(0);
    //pd2=cpair.getParticleDataEntry(1);
    pd1=jamTable->find(kf1);
    pd2=jamTable->find(kf2);
    int ipair=cpair.getPairID();

    if(kf1<0 || kf2<0) {
	cout << "SigmaBB2::calc kf<0? kf1= "<< kf1 << " kf2= " << kf2<<endl;
	exit(1);
    }

//.....nucleon - nucleon
    if(ipair == CollisionPair::Pair(id_nucl,id_nucl)) {
	sigmaNuclNucl(cpair);
	return;
    }

    double sigab=0.0;
    mChanel=0;
    int iabs=-1;

//...Resonance cross sections. Skip calculate resonance cross sections if energy
// is high.
    if(srt < 100.0) {


//...First get inel. cross sections at == al c.m.s. momentum.
//.......N + N* || N* + N*  inelastic collisions
        if( ( ipair == CollisionPair::Pair(id_nucl, id_nucls) ) ||
            ( ipair == CollisionPair::Pair(id_nucls,id_nucls) )  ) {

            iabs=sigmaNsNs(cpair);

//.......N(*) + Delta(*) inelastic collisions
	} else if( ( ipair == CollisionPair::Pair(id_nucl, id_delt)  ) || 
                   ( ipair == CollisionPair::Pair(id_nucl, id_delts) ) ||
                   ( ipair == CollisionPair::Pair(id_nucls,id_delt)  ) ||
                   ( ipair == CollisionPair::Pair(id_nucls,id_delts))) {

	    iabs=sigmaNsDs(cpair);

//.....Delta(*) - Delta(*)  inelastic collisions
	} else if((ipair == CollisionPair::Pair(id_delt, id_delt)  ) ||
                  (ipair == CollisionPair::Pair(id_delt, id_delts) ) || 
                  (ipair == CollisionPair::Pair(id_delts,id_delts))) {

	    iabs=sigmaDD(cpair);

	} else {

	    cout << "SigmaBB2:: BB collision something is wrong PID id1= " << id1
		 << " id2= " << id2 
		 << endl;
	    exit(1);

	}

        if(iabs>=0) {
	  if(absorptionBB) {sigab=sigin[iabs];
	  } else {
	    cpair.setOutGoing(iabs,0.0,2212, 2112);//->pn
	    sigin[iabs]=0.0;
	  }
	//cout << " iabs= "<< iabs << " abs= "<< absorptionBB<< " sig= "<< sigab << endl;
	//cin.get();
	}

	cpair.setSigAbs(sigab);

    }

    double snew=2*sqrt(Mnucl*Mnucl+pr*pr);
    jamxnn(snew,3,3);
    double sigt1 = sigTot;
    double sigel1 = sigEl;

    jamxnn(snew,3,0);
    double sigt2 = sigTot;
    double sigel2 = sigEl;
    sigTot=0.5*(sigt1+sigt2);
    sigEl=0.5*(sigel1+sigel2);

    //if(srt <= 3.6) {
    if(srt <= 3.2) {
	sigTot=sigEl;
	//for(int i=0;i<mChanel;i++) sigTot +=sigin[i];

	double sigt=0.0;
	for(int i=0;i<mChanel;i++) sigt +=sigin[i];
	sigTot += sigt;

	if(sigTot-sigEl-sigt > 0.0001) {
	     cout << " di1= " << id1 << " id2= " << id2 <<endl;
	     cout << " sigtot= " << sigTot << " sigel= " << sigEl<< endl;
	     cout << " sigt= " << sigt <<endl;
	     cout << " diff= " << sigTot-sigEl-sigt <<endl;
	     cout << " s= " << srt << " snew= " << snew<<endl;
	     exit(1);
	 }

	/*
	cout << " id1= "<< kf1 << " id2= "<< kf2
	     << " ecm= "<< srt << " snew= "<< snew
	    << " sigt= "<< sigTot << " sigel= " << sigEl
	    <<endl;
	for(int i=0;i<mChanel;i++) {
	   cout << i << " sigin= "<< sigin[i]
	       << endl;
	}
	cin.get();
	*/


    } else {
//...Add absorption cross sections.
	double siginTot=0.0;
	for( int i=0;i<mChanel;i++) siginTot += sigin[i];
	siginTot -= sigab;
	double sigsoft=sigTot-sigEl-siginTot;
	sigTot=sigEl+siginTot+sigsoft+sigab;

	// add string excitation.
	if(sigsoft>0)  {
	    mChanel++;
	    cpair.setOutGoing(sigsoft,92,92);
	}

    }

    cpair.setXSTotal(sigTot);
    cpair.setXSElastic(sigEl);


    /*
    cpair.checkChannel();

	  double sigt=0.0;
	  for( int i=0;i<mChanel;i++) sigt += sigin[i];
	  //if(sigTot-sigEl-sigt > 0.00001) {
		cout << "SigBB::RR sigsoft<0? srt= "<< srt
		    << " mMode= "<< mMode
		    << " kf1= " << kf1 << " kf2= " << kf2
		    << " siginTot= " << sigt
		    << " sig= " << sigTot
		    << " sigel= " << sigEl
	            << " diff= " << sigTot-sigEl-sigt
		    << " mChannel= " << mChanel
		    << " outgoingsize= " <<cpair.outgoingSize()
		    <<endl;
	         for( int i=0;i<mChanel;i++) 
		     cout << i << " sigin= " << sigin[i] <<endl;;
		   cout << " sigab= " << sigab<<endl;
	//	exit(1);
	   // }

	   */

}

// Compute cross sections for nucleon-nucleon scatterings.
void SigmaBB2::sigmaNuclNucl(CollisionPair& cpair)
{
  // First compute sig and sigel from cross section table.
  jamxnn(srt,iz1,iz2);
  cpair.setSigma(sigTot,sigEl);

  // Resonance cross section alone cannot fill the inelastic cross
  // section at high energy, we use parametrized sig and sigel.
  if(mMode==1 && srt >= srtStr) return;

  if(srt <= sMin+eKinMin) {
    //cpair.setSigma(sigTot,sigTot);
    cpair.setSigma(sigEl,sigEl);
    return;
  }

  if(kf1 == kf2) sigmaPP(cpair);
  else sigmaPN(cpair);

  double siginTot=0.0;
  for( int i=0;i<mChanel;i++) siginTot += sigin[i];

  if(srt < srtStr) {
    sigTot = sigEl + siginTot;
    cpair.setSigma(sigTot,sigEl);
    return;
  }

  // add string excitation cross section.
  double sigsoft=sigTot-sigEl-siginTot;
  if(sigsoft > 0) {

    mChanel++;
    cpair.setOutGoing(sigsoft,92,92);

    // test: string excitation only.
    //mChanel=1;
    //double sigsoft=sigTot-sigEl;
    //cpair.setOutGoing(sigsoft,92,92);

  } else {

      /*
    cout << "SigmaBB2::sigmaNuclNucl sigsoft<0? srt= "<< srt
         << " siginTot= " << siginTot
         << " sig= " << sigTot
         << " sigel= " << sigEl
         << " sigsoft= " << sigsoft
         << " sigt= " << sigEl + siginTot + sigsoft
         <<endl;
	 */
  }

}


//============================================================================
// Compute PN -> X
// 1) NN -> ND  2) NN -> NN*  3) NN -> DD  4) NN -> ND*  5) NN -> N*D
// 6) NN -> DD*  7) NN -> N*N*  8) NN -> N*D*  9) NN -> D*D* 
void SigmaBB2::sigmaPN(CollisionPair& cpair)
{

// 1)  np => d+(1232) + n   2) np => d0(1232) + p 
// 3)  np => np*            4) np => pn*
// 5)  np => d0d+           6) np => d-d++            
// 7)  np => nd*+           8) np => pd*0
// 9)) np=>n* d+           10) np=>p* d0
// 11) np=>d0 d*+          12) np=>d- d*++
// 13) np=>n* p*
// 14) np=>n* d*+          15) np=>p* d*0
// 16) np=>d*0 d*+         17) np=>d*- d*++
// 18,19) s-wave


  //..First compute T=1 cross section.
  double sig1[10],sig0[10];
  sigmaRR->sigma(srt,Mnucl,Mnucl,sig1,1);

  //...Compute T=0 cross section.
  sigmaRR->sigma(srt,Mnucl,Mnucl,sig0,0);

  sigin[0]=0.25*sig1[0];        // nd+
  sigin[1]=0.25*sig1[0];        // pd0
  sigin[2]=0.25*(sig1[1]+sig0[1]);        // np*
  sigin[3]=0.25*(sig1[1]+sig0[1]);        // n*p
  sigin[4]=0.05*sig1[2]+0.25*sig0[2];   // d0d+
  sigin[5]=0.45*sig1[2]+0.25*sig0[2];   // d-d++
  sigin[6]=0.25*sig1[3];                  // nd*+
  sigin[7]=0.25*sig1[3];                  // pd*0
  sigin[8]=0.25*sig1[4];                  // n*d+
  sigin[9]=0.25*sig1[4];                 // p*d0
  sigin[10]=0.05*sig1[5]+0.25*sig0[5];  // d0d*+
  sigin[11]=0.45*sig1[5]+0.25*sig0[5];  // d-d*++
  sigin[12]=0.5*(sig1[6]+sig0[6]);        // n*p*
  sigin[13]=0.25*sig1[7];                 // n*d*+
  sigin[14]=0.25*sig1[7];                 // p*d*0
  sigin[15]=0.05*sig1[8]+0.25*sig0[8];  // d*0d*+
  sigin[16]=0.45*sig1[8]+0.25*sig0[8];  // d*-d*++

  mChanel=17;
  //mChanel=2; // delta only
  cpair.setSizeOutGoing(mChanel);

  for(int i=0;i<mChanel;i++)  {
    cpair.setOutGoing(i,sigin[i],BBpn[i][0],BBpn[i][1]);

    /*
     //OutGoing out = cpair.getOutGoing();
      cout << "PN sig= " << i << " sigin= " << sigin[i]
	  << " id1= " << BBpn[i][0]
	  << " id2= " << BBpn[i][1]
	  <<endl;
	  */
  }

  /*
  double sigint=0.0;
  for(int i=0;i<mChanel;i++) sigint += sigin[i];
  //jamxnn(srt,iz1,iz2);
  double sigsoft=sigTot-sigEl-sigint;
  cout << "srt= "<< srt << " sig= " << sigTot << " sigel= " << sigEl
       << " sig-sigel=" << sigTot - sigEl
       << " sigintot= " << sigint
       << " sigsoft= " << sigsoft
       << " sigint= " << sigint + sigsoft
       <<endl;
  */

}

//============================================================================
// Compute PP/NN -> X
// 1) NN -> ND  2) NN -> NN*  3) NN -> DD  4) NN -> ND*  5) NN -> N*D
// 6) NN -> DD*  7) NN -> N*N*  8) NN -> N*D*  9) NN -> D*D* 
void SigmaBB2::sigmaPP(CollisionPair& cpair)
{
//.....Resonance production cross sections.
    double sig1[10];
    sigmaRR->sigma(srt,Mnucl,Mnucl,sig1,1);

// 1) pp=>d+p /nn=>nd0        2) pp=>d++n /nn=>pd-
// 3) pp=>pp* /nn=>nn*
// 4) pp=>d+d+/nn=>d0d0       5) pp=>d0d++/nn=>d+d-
// 6) pp=>p d*+/nn=>nd*0      7) pp=>n d*++/nn=>pd*-
// 8) pp=>p* d+/nn=>n*d0      9) pp=>n* d++/nn=>p*d-
// 10)pp=>d+ d*+/nn=>d0d*0   11) pp=>d0 d*++/nn=>d+d*-
// 12)pp=>p* p*/nn=>n*n*
// 13)pp=>p* d*+/nn=>n*d0*   14) pp=>n* d*++/nn=>p*d*-
// 15)pp=>d*+d*+/nn=>d*0d*0  16) pp=>d*0d*++/nn=>d*+d*-
// 17) s-wave

    mChanel=16;
    //mChanel=2; // delta only
    sigin[0]=0.25*sig1[0];   // pd+
    sigin[1]=0.75*sig1[0];   // nd++

//...deuteron + pi+ cross section.
    //double sigdpi=0.0;
    //if(optDeut1>= 2) sigdpi = jamxdpi1(srt);
    //if(optDeut1 == 2) sigin[1] += sigdpi;
    //else if(optDeut1 >= 3) sigin[17]=sigdpi;
    //if(optDeut2 >= 1) sigin[18]=jamxdpi2(1,srt);

    double sym = 0.5;
    sigin[2]=sig1[1];       // pp*
    sigin[3]=0.4*sig1[2]*sym;   // d+d+
    sigin[4]=0.6*sig1[2];   // d0d++
    sigin[5]=0.25*sig1[3];  // pd*+
    sigin[6]=0.75*sig1[3];  // nd*++
    sigin[7]=0.25*sig1[4];  // p*d+
    sigin[8]=0.75*sig1[4];  // n*d++
    sigin[9]=0.4*sig1[5];  // d*+d+
    sigin[10]=0.6*sig1[5];  // d*d++
    sigin[11]=sig1[6];  // p*p*  symmetry factor is already included
    sigin[12]=0.25*sig1[7]; // p*d*+
    sigin[13]=0.75*sig1[7]; // n*d*++
    sigin[14]=0.4*sig1[8];  // d*+d*+
    sigin[15]=0.6*sig1[8];  // d*0d*++

    cpair.setSizeOutGoing(mChanel);

    // nn collision
    if(iz1+iz2==0) {

	 for(int i=0;i<16;i++) 
	 cpair.setOutGoing(i,sigin[i],BBnn[i][0],BBnn[i][1]);
	 /*
	 if(optSwave==1) {
	   sigin[16]=sig1[9];
	   cpair.setOutGoing(16,sigin[16],2112,12112); // n n(1440)
	 } else {
           sigin[16]=sig1[9]/3.0;    // s-wave pp->pp pi0
           sigin[17]=sig1[9]*2.0/3.0; // s-wave pp->pn pi+
	   cpair.setOutGoing(16,sigin[16],2112,2112,111);
	   cpair.setOutGoing(17,sigin[17],2212,2112,-211);
	 }
	 */

    // pp collision
    } else {

	 for(int i=0;i<16;i++) 
	 cpair.setOutGoing(i,sigin[i],BBpp[i][0],BBpp[i][1]);
	 /*
	 if(optSwave==1) {
	   sigin[16]=sig1[9];
	   cpair.setOutGoing(16,sigin[16],2212,12212); // p p(1440)
	 } else {
           sigin[16]=sig1[9]/3.0;    // s-wave pp->pp pi0
           sigin[17]=sig1[9]*2.0/3.0; // s-wave pp->pn pi+
	   cpair.setOutGoing(16,sigin[16],2212,2212,111);
	   cpair.setOutGoing(17,sigin[17],2212,2112,211);
	 }
	 */
    }

}

//**********************************************************************
//  select N* or D* resonance.

bool SigmaBB2::sampleResonanceMass(double ecm, int idn1, int idn2, double& m1, double& m2)
{    
  //int ch = cpair.selectChannel();
  //OutGoing out = cpair.getOutGoing();
  //int idn1= out.id1;
  //int idn2= out.id2;
  //double ecm=cpair.getCMenergy();

  double emr[2],em0[2],gam0[2];
  double mmax[2],mmin[2],bwmax[2],bw[2],ymin[2],ymax[2];
  int iex[2], id[2],idv[2];

  vector<ParticleDataEntry*> table1 = pickRTable(idn1,idv[0]);
  vector<ParticleDataEntry*> table2 = pickRTable(idn2,idv[1]);
  int n1=table1.size();
  int n2=table2.size();
  if(n1==0 || n2==0) {
      cout << "SigmaBB2::sampleResonanceMass no baryon table ? n1= "<< n1
	  << " n2= "<< n2<<endl;
      exit(1);
  }

  // R R -> N N
  if(idv[0]==id_nucl && idv[1]==id_nucl) {
    m1 = table1[0]->m0();
    m2 = table2[0]->m0();
    pout[0] = table1[0];
    pout[1] = table2[0];
    return true;
  }

  //cout << "idv1= "<< idv[0] << " idv2= "<< idv[1] <<endl;
  //cin.get();

  int iso=1;
  int ip=-2;
  if(idv[0] == id_nucl && idv[1] == id_delt) { //1) NN -> ND  
    ip = 0;
  } else if(idv[0] == id_nucl && idv[1] == id_nucls) { //2) NN -> NN*
    ip = 1;
    if(idn1!=idn2) iso=0;
  }else if(idv[0] == id_delt && idv[1] == id_delt) { //3) NN -> DD
    ip = 2;
    //if(table1[0]->chargeType() + table2[0]->chargeType() == 3) {
    if( idn1 + idn2 == 9) {
      iso=0;
      //if(table2[0]->chargeType()==6) iso=2;
      if(idn2 ==6) iso=2;
    }
  } else if(idv[0] == id_nucl && idv[1] == id_delts) { //4) NN -> ND*
    ip = 3;
  } else if(idv[0] == id_delt && idv[1] == id_nucls) { //5) NN -> N*D
    ip = 4;
  } else if(idv[0] == id_delt && idv[1] == id_delts) { //6) NN -> DD*
    ip = 5;
    if(table1[0]->chargeType() + table2[0]->chargeType() == 3) {
      iso=0;
      if(table2[0]->chargeType()==6) iso=2;
    }
  } else if(idv[0] == id_nucls && idv[1] == id_nucls) { //7) NN -> N*N*
    ip = 6;
    if(idn1!=idn2) iso=0;
  } else if(idv[0] == id_nucls && idv[1] == id_delts) { //8) NN -> N*D*
    ip = 7;
  } else if(idv[0] == id_delts && idv[1] == id_delts) { //9) NN -> D*D*
    ip = 8;
    if(table1[0]->chargeType() + table2[0]->chargeType() == 3) {
      iso=0;
      if(table2[0]->chargeType()==6) iso=2;
    }
  } else {
    cout << "SigmaBB2::sampleResonanceMass wrong reaction type "<<endl;
    cout << "idv[0]= " << idv[0] << " idv[1]= "<< idv[1]<<endl;
    exit(1);
  }
  
  double totsp=0.0;
  vector<vector<double> > probr(n1, vector<double>(n2,0.0));
  for(int i=0;i<n1;i++)
  for(int j=0;j<n2;j++) { 
    probr[i][j]=0.0;
    ParticleDataEntry* pv1=table1[i];
    ParticleDataEntry* pv2=table2[j];
    double pbw = sigmaRR->prob(ip,ecm,iso,pv1->id(),pv2->id());
    if(pbw <= 0.0) continue;
    int ispin1=pv1->spinType();
    int ispin2=pv2->spinType();
    probr[i][j]= ispin1*ispin2 * pbw;
    totsp += probr[i][j];
  }

  int mtry=0;
  double pf0=0.0;
  do {
  double xrand=rndm->flat()*totsp;
  for(int i=0;i<n1;i++)
  for(int j=0;j<n2;j++) {
    xrand -= probr[i][j];
    if(xrand <= 0) {
      pout[0] = table1[i];
      pout[1] = table2[j];
      goto L100;
    }   
  }   
  cout << " SigmaBB2::sampleResonanceMass no particle? "<< srt << endl;
  exit(1);
L100:

  bool ok=false;
    // Select resonance.
    //pout[0]=selectDN(idn1,idv[0]);
    //pout[1]=selectDN(idn2,idv[1]);
    
    id[0]=pout[0]->id();
    id[1]=pout[1]->id();
    iex[0]=1; iex[1]=1;
    mmin[0]=pout[0]->mMin();
    mmin[1]=pout[1]->mMin();
    if(pout[0]->mWidth() < 1e-5 ) {mmin[0]=pout[0]->m0();iex[0]=0;}
    if(pout[1]->mWidth() < 1e-5 ) {mmin[1]=pout[1]->m0();iex[1]=0;}
    mmax[0]=mmin[0];
    mmax[1]=mmin[1];
   if(iex[0]==1) mmax[0]=min(maxRMass,ecm-mmin[1]-eKinMin);
   if(iex[1]==1) mmax[1]=min(maxRMass,ecm-mmin[0]-eKinMin);
   if(mmin[0] <= mmax[0] && mmin[1] <= mmax[1]) ok=true; 

  if(!ok) {
    cout << "SigmaBB2::ResnanceMass too small srt? srt= " << ecm <<endl;
    cout << " id1= " << id[0] << " min1= "<< mmin[0]<< " max= "<< mmax[0]
         << " iex= " << iex[0] <<endl;
    cout << " id2= " << id[1] << " min2= " << mmin[1]<< " max= "<< mmax[1]
         << " iex= " << iex[1] <<endl;
    return false;
  }

  if(ecm < mmin[0] + mmin[1] + eKinMin) {
     cout << "sampleResonanceMass too small srt?" << ecm
	   << " mmin1= " << mmin[0]
	   << " mmin2= " << mmin[1]
	   <<endl;
    return false;
  }


  // First compute maximum value.
  for(int jt=0; jt<2;jt++) {
    bwmax[jt]=1.0;
    bw[jt]=1.0;
    ymin[jt]=0.0;
    ymax[jt]=0.0;
    em0[jt]=pout[jt]->m0();
    gam0[jt]=pout[jt]->mWidth();
    emr[jt]=em0[jt];
    if(iex[jt] == 0) continue;
    ymin[jt]=atan((pow2(mmin[jt])-pow2(em0[jt]))/(em0[jt]*gam0[jt]));
    ymax[jt]=atan((pow2(mmax[jt])-pow2(em0[jt]))/(em0[jt]*gam0[jt]));

    if(mmax[jt]<em0[jt]) {
    double wid =  decay->getTotalWidth(pout[jt],mmax[jt]);
    bwmax[jt] = BW(em0[jt],mmax[jt],wid)/BW(em0[jt],mmax[jt],gam0[jt]);
    } else {
    double wid = decay->getTotalWidth(pout[jt],em0[jt]);
    bwmax[jt] = BW(em0[jt],em0[jt],wid)/BW(em0[jt],em0[jt],gam0[jt]);
    }

    if(idv[jt] == id_nucl)  bwmax[jt] *= 5.0;
    if(idv[jt] == id_delt)  bwmax[jt] *= 5.0;
    if(idv[jt] == id_nucls) bwmax[jt] *= 5.0;
    if(idv[jt] == id_delts) bwmax[jt] *= 5.0;
  }

  pf0=PCM(ecm,mmin[0],mmin[1]);
  double bwpmax=bwmax[0]*bwmax[1]*pf0;

  // Generate resonance mass according to Breit-Wigner + phase space.
  int  ntry=0;
  double pf=0.0;
  do {
    do {
      for(int jt=0;jt<2;jt++)
      if(iex[jt] == 1) {
        // First generate mass according to B-W with constant width.
        // A(m^2)=d(m^2)/((m^2 - m_R^2)^2 + (m*Gamma)^2)
        double y=ymin[jt] + rndm->flat()*( ymax[jt] - ymin[jt] );
        emr[jt] = sqrt( em0[jt]*gam0[jt]*tan(y) + em0[jt]*em0[jt] );
      }
    } while(emr[0]+emr[1]+eKinMin > ecm);

    for(int jt=0;jt<2;jt++)
    if(iex[jt] == 1) {
      double wid=decay->getTotalWidth(pout[jt],emr[jt]);
      bw[jt] = BW(em0[jt],emr[jt],wid)/BW(em0[jt],emr[jt],gam0[jt]);
    }

    //...Check final phase space.
    pf=PCM(ecm,emr[0],emr[1]);
    if(bwpmax < bw[0]*bw[1]*pf) {
      cout << "BB bwmax= " << bwpmax << " srt= " << ecm
	  << " bw= "<< bw[0]*bw[1]*pf
   	   << " id1= " << idv[0]
	   << " id2= " << idv[1]
	   << endl;
     cout << " n1= "<< n1 << " n2= "<< n2
	 << " idn1= "<< idn1 << " idn2= "<< idn2
	 <<endl;
      bwpmax *= 1.05;
      continue;
    }

    if(ntry++ > 10) break;
      
  } while (bwpmax*rndm->flat() > bw[0]*bw[1]*pf);

  if(ntry<=10) {
    m1=emr[0]; m2=emr[1];

    /*
    cout << " srt= "<< ecm << " srt= "<< srt << " n1= "<< n1
	<< " n2= "<< n2 <<endl;
  cout << "id1= "<< idn1 << " id2= "<< idn2
      << " -> id3= "<< id[0] << " m3= "<< m1 << " id4= "<< id[1]
      << " m4= "<< m2
      <<endl;
  cin.get();
  */

    return true;
  }

  } while (mtry++ < 100);

  if(nError <5) {
  cout << "warning SigBB2::sampelResonanceMass does not converge srt=" << ecm
       << " mtry= "<< mtry
       << " id1= " << id[0] << " id2= " << id[1]
       << " em1= " << emr[0] << " em2= " << emr[1]
       <<endl;
    nError++;
  }

  return sampleMassFix(ecm,pf0,iex,emr,em0,gam0,mmin,ymin,ymax,m1,m2);

}

int SigmaBB2::sigmaNsNs(CollisionPair& cpair)
{
//...Purpose : to give cross sections  || n + n* || n* + n* collisions.
//---------------------------------------------------------------------*
//...Outputs
//...sigin(): inel. cross sections.
//...mchanel: max. inel. channel.
//...at the moment dn -> x  cross sections are simply == ated to
//...t=1 cross sections at == al cms momentum
//---------------------------------------------------------------------*

//...See whether n + n* || n* + n* collision.
//...Calculate N(*)N(*) --> NN absorption cross section
//...from detailed balance.
  int iabs=-1;
  int isig=0;
  double m01, m02;
  int instar=100;
  // particle 1 is nucleon?  N_1 N*_2
  if( kf1 == 2112  ||  kf1 == 2212 ) {
    isig=1;
    instar = jamTable->nStarID(kf2)-1;
    m01=em1;
    m02=sigmaRR->getNstar(instar)->mass0();
  // particle 2 is nucleon?  N*_1 N_2
  } else if( kf2 == 2112  ||  kf2 == 2212 ) {
    isig=1;
    instar = jamTable->nStarID(kf1)-1;
    m01=sigmaRR->getNstar(instar)->mass0();
    m02=em2;
   // N*_1 N*_2
   } else {
     int ins1 = jamTable->nStarID(kf1)-1;
     int ins2 = jamTable->nStarID(kf2)-1;
     m01=sigmaRR->getNstar(ins1)->mass0();
     m02=sigmaRR->getNstar(ins2)->mass0();
     if(optRR==1) isig= 6;
  }

  double sig1[10],sig0[10];
  double sigab=0.0;
  int icc;
  //double snew=2*sqrt(Mnucl*Mnucl+pr*pr);

//...Get cross sections from nn collision.
  double faci=1.0, facf=1.0;
  if( iz1  ==  iz2 ) {
        icc=1;
	if(iz1+iz2==0) icc=2;
        faci=2.0;
	//if(isig==6) facf=2.0;
	//if(kf1==kf2) facf=2.0; // taken into account identical particle.
	if(isig>0) {
          sigab=sigmaRR->matrixElm(isig,srt,1,m01,m02,instar);
	}
        //sigmaRR->sigma(snew,em1,em2,sig1,1);
        sigmaRR->sigma(srt,em1,em2,sig1,1);
  } else {
        icc=0;
	if(isig>0) {
          //sigmaRR->sigma(srt,em1,em2,sig1,1);
          //sigmaRR->sigma(srt,em1,em2,sig0,0);
          sigab=0.25*(sigmaRR->matrixElm(isig,srt,1,m01,m02,instar)
                     +sigmaRR->matrixElm(isig,srt,0,m01,m02,instar));
	}

        sigmaRR->sigma(srt,em1,em2,sig1,1);
        sigmaRR->sigma(srt,em1,em2,sig0,0);
  }

  double prnew=sqrt(0.25*srt*srt-Mnucl*Mnucl);
  double flx=srt*srt*PCM(srt,em1,em2);
  sigab = 4.0*facf/faci*prnew/flx*sigab;

    if(icc >= 1) {
	mChanel=18;
	double sym=0.5;
        cpair.setSizeOutGoing(mChanel);
        sigin[0]=0.25*sig1[0]; // pd+
        sigin[1]=0.75*sig1[0]; // nd++
        sigin[2]=sig1[1];        // pp*+
        sigin[3]=0.4*sig1[2]*sym;  // d+d+
        sigin[4]=0.6*sig1[2];  // d0d++
        sigin[5]=0.25*sig1[3]; // pd*+
        sigin[6]=0.75*sig1[3]; // nd*++
        sigin[7]=0.25*sig1[4]; // p*d+
        sigin[8]=0.75*sig1[4]; // n*d++
        sigin[9]=0.4*sig1[5]*sym; // d*+d+
        sigin[10]=0.6*sig1[5]; // d*0d++
        sigin[11]=sig1[6]*sym;       // p*p*
        sigin[12]=0.25*sig1[7];// p*d*+
        sigin[13]=0.75*sig1[7];// n*d*++
        sigin[14]=0.4*sig1[8]; // d*+d*+
        sigin[15]=0.6*sig1[8]; // d*0d*++
        sigin[16]=sig1[9];      // s-wave
        sigin[17]=sigab;         // b*b*->nn
	iabs=17;

	// nn collision
	if(icc==2) {
	 for(int i=0;i<16;i++) {
	    cpair.setOutGoing(i,sigin[i],BBnn[i][0],BBnn[i][1]);
	 }
	    cpair.setOutGoing(16,sigin[16],2112,12112); // ->n n(1440)
	    cpair.setOutGoing(17,sigin[17],2112,2112); // ->n n
	 // pp collision
	} else {
	    for(int i=0;i<16;i++)  {
	    cpair.setOutGoing(i,sigin[i],BBpp[i][0],BBpp[i][1]);
	    }
	    cpair.setOutGoing(16,sigin[16],2212,12212); // p p(1440)

	    cpair.setOutGoing(17,sigin[17],2212,2212); // p p
	}

    } else {
        mChanel=20;
        cpair.setSizeOutGoing(mChanel);
        sigin[0]=0.25*sig1[0];                   // nd+
        sigin[1]=0.25*sig1[0];                   // pd0
        sigin[2]=0.25*(sig1[1]+sig0[1]);         // np*
        sigin[3]=0.25*(sig1[1]+sig0[1]);         // n*p
        sigin[4]=0.05*sig1[2]+0.25*sig0[2];  // d0d+
        sigin[5]=0.45*sig1[2]+0.25*sig0[2];  // d-d++
        sigin[6]=0.25*sig1[3];                   // nd*+
        sigin[7]=0.25*sig1[3];                   // pd*0
        sigin[8]=0.25*sig1[4];                   // n*d+
        sigin[9]=0.25*sig1[4];                  // p*d0
        sigin[10]=0.05*sig1[5]+0.25*sig0[5]; // d0d*+
        sigin[11]=0.45*sig1[5]+0.25*sig0[5]; // d-d*++
        sigin[12]=0.5*(sig1[6]+sig0[6]);        // n*p*
        sigin[13]=0.25*sig1[7];                  // n*d*+
        sigin[14]=0.25*sig1[7];                  // p*d*0
        sigin[15]=0.05*sig1[8]+0.25*sig0[8]; // d*0d*+
        sigin[16]=0.45*sig1[8]+0.25*sig0[8]; // d*-d*++
        sigin[17]=0.5*sig0[9];                  // s-wave
        sigin[18]=0.5*sig0[9];                  // s-wave
        sigin[19]=sigab;                           // b*b*->nn
	iabs=19;

	for(int i=0;i<17;i++) 
	cpair.setOutGoing(i,sigin[i],BBpn[i][0],BBpn[i][1]);
	cpair.setOutGoing(17,sigin[17],2212,12112);//pn->pn*(1440)
	cpair.setOutGoing(18,sigin[18],2112,12212);//pn->np*(1440)
	cpair.setOutGoing(19,sigin[19],2212, 2112);//->pn
    }

    return iabs;

}

//***********************************************************************

int SigmaBB2::sigmaNsDs(CollisionPair& cpair)
{
//...Purpose : to give cross sections  || n(*) + d(*) collisions.
//
//   srt (GeV) : c.m. energy of the colliding pair.
//   id1, id2:  particle ID
//   ires1, ires2: resonance number
//
// ... at the moment dn -> x  cross sections are simply == ated to
//     t=1 cross sections at == al cms momentum
//=======================================================================

  int itz=(iz1+iz2)/3;

  //...2 is a nucleon(*) 1 is a delta(*)
  int izn=iz2/3;
  int izd=iz1/3;
  int iresd = id1 == id_delt? 0 : 1;
  int iresn = id2 == id_nucl? 0 : 1;
  double emn=em2;
  int in = jamTable->nStarID(kf2);
  int id = jamTable->dStarID(kf1);

  //...1 is a nucleon(*) 2 is a delta(*)
  if(id1 == id_nucl  ||  id1 == id_nucls) {
    izn=iz1/3;
    izd=iz2/3;
    iresn = id1 == id_nucl? 0 : 1;
    iresd = id2 == id_delt? 0 : 1;
    emn=em1;
    in = jamTable->nStarID(kf1);
    id = jamTable->dStarID(kf2);
  }

  double elm=0.0;
  if(itz > -1 && itz < 3) {
    double em0n= in==0 ? emn   : sigmaRR->getNstar(in-1)->mass0();
    double em0d= id==0 ? Mdelta: sigmaRR->getDstar(id-1)->mass0();
    if(iresn == 0) {
      if(iresd == 0) {
	elm = sigmaRR->matrixElm(0,srt); // delta n ==> n n
      } else {
	elm = sigmaRR->matrixElm(2,srt,1,em0n,em0d); // delta* n ==> n n
      }
    } else {
      if(iresd == 0) {
	elm = sigmaRR->matrixElm(4,srt,1,em0n,em0d); // delta n* ==> n n
      } else {
	if(optRR==1)
	elm = sigmaRR->matrixElm(7,srt,1,em0n,em0d);//delta* n* ==> nn 
      }
    }
  }

//...Calculate nn cross section.
  //double snew=2*sqrt(Mnucl*Mnucl+pr*pr);
  double sig1[10];
  sigmaRR->sigma(srt,em1,em2,sig1,1);
  int iabs=-1;

  //cout << " srt= "<< srt << " em1= "<< em1 << " em2= "<< em2<<endl;
  //for(int i=0;i<9;i++)
  //cout << " sig1= "<< sig1[i]<<endl;

  // n d-
  if(itz == -1) {
    mChanel=9;
    cpair.setSizeOutGoing(mChanel);
    for(int i=0;i<9;i++)  {
      if(NDm1[i][0]==0) sig1[i]=0.0;
      sigin[i]=sig1[i];
      cpair.setOutGoing(i,sig1[i],NDm1[i][0],NDm1[i][1]);
    }

  // p d++
  } else if(itz == 3) {
    mChanel=9;
    cpair.setSizeOutGoing(mChanel);
    for(int i=0;i<9;i++)  {
      if(ND3[i][0]==0) sig1[i]=0.0;
      sigin[i]=sig1[i];
      cpair.setOutGoing(i,sig1[i],ND3[i][0],ND3[i][1]);
    }

  // p d-  n d0
  } else if(itz == 0) {
    mChanel=18;
    cpair.setSizeOutGoing(mChanel);
    makePP(sig1);
    for(int i=0;i<16;i++) 
    cpair.setOutGoing(i,sigin[i],BBnn[i][0],BBnn[i][1]);
    cpair.setOutGoing(16,sigin[16],2112,12112); // ->n n(1440)

    double faci = 0.75;
    //....pd- --> nn
    if( izn == 1 && izd == -1) faci = 0.75;
    //....nd0- --> nn  
    else if(izn == 0 && izd == 0) faci = 0.25;

    double prnew=sqrt(0.25*srt*srt-Mnucl*Mnucl);
    double flx=srt*srt*PCM(srt,em1,em2);
    double sigab=4.0*0.5*faci*prnew/flx*elm;
    sigin[17]=sigab;
    cpair.setOutGoing(17,sigab,2112,2112);
    iabs=17;

  // p d0 or n d+
  } else if(itz == 1) {
    mChanel=20;
    cpair.setSizeOutGoing(mChanel);
    makePN(sig1);
    for(int i=0;i<17;i++) 
    cpair.setOutGoing(i,sigin[i],BBpn[i][0],BBpn[i][1]);
    cpair.setOutGoing(17,sigin[17], 2212,12112);//pn->pn*(1440)
    cpair.setOutGoing(18,sigin[18], 2112,12212);//pn->np*(1440)

    //....pd0 --> pn  11:nd+ --> pn
    double prnew=sqrt(0.25*srt*srt-Mnucl*Mnucl);
    double flx=srt*srt*PCM(srt,em1,em2);
    double sigab=4.0*0.25*prnew/flx*elm;
    sigin[19]=sigab;
    cpair.setOutGoing(19,sigab,2212,2112);
    iabs=19;

  // p d+ or n d++
  } else if(itz == 2) {
    mChanel=18;
    cpair.setSizeOutGoing(mChanel);
    makePP(sig1);
    for(int i=0;i<16;i++) 
    cpair.setOutGoing(i,sigin[i],BBpp[i][0],BBpp[i][1]);
    cpair.setOutGoing(16,sigin[16],2212,12212); // p p(1440)

    //16:nd++ --> pp 0.75 //12:pd+ --> pp 0.25
    double faci = (izn == 0 && izd == 2) ? 0.75 : 0.25;
    double prnew=sqrt(0.25*srt*srt-Mnucl*Mnucl);
    double flx=srt*srt*PCM(srt,em1,em2);
    double sigab=4.0*0.5*faci*prnew/flx*elm;
    sigin[17]=sigab;
    cpair.setOutGoing(17,sigab,2212,2212);
    iabs=17;

    /*
    if(izn==0 && izd==2) {
      cout << "nd++ -> pp sig= "<< sigab << " srt= "<< srt;
    } else {
      cout << "pd+ -> pp sig= "<< sigab << " srt= "<< srt;
    }
    cout << " faci= "<< faci;
    cout << " em1= " << em1 << " em2= "<< em2
	<< " elm= "<< elm
	<< " flx= " << flx
	<< " prnew= " << pr2new
	  <<endl;
    cout << " em1= " << em1 << " em2= "<< em2
	<< " elm= "<< elm
	<< " flx= " << flx
	<< " prnew= " << pr2new
	<<endl;
    cin.get();
    */

  }

  return iabs;
}

//**********************************************************************
//...Purpose : to give cross sections  || d(*) + d(*) collisions.
//...at the moment dd -> x  cross sections are simply equated to
//...t=1 cross sections at equal cms momentum
int SigmaBB2::sigmaDD(CollisionPair& cpair)
{

  int iz01=iz1/3;
  int iz02=iz2/3;
  int iztot=iz01+iz02;
  int izmin=min(iz01,iz02);
  int izmax=max(iz01,iz02);

  int i1=0;
  int i2=0;
  int isig=0;
    // d(1232) + d(1232)
    if(id1 == id_delt && id2 == id_delt) {
        isig=3;

    // d(1232) + d*
    } else if(id1 == id_delt  && id2 == id_delts) {
        isig=5;
        i2 = jamTable->dStarID(kf2);
    // d* + d(1232)
    } else if(id1 == id_delts && id2 == id_delt ) {
        isig=5;
        i1 = jamTable->nStarID(kf1);
    // d* + d*
    } else if(id1 == id_delts && id2 == id_delts) {
        if(optRR==1) isig=8;
        i1 = jamTable->nStarID(kf1);
        i2 = jamTable->dStarID(kf2);
    } else {
	cout << "(sigmaBB2::sigmaDD) error in d+d collisions "
	<< " id1= " << id1 << " id2= " << id2 << endl;
	exit(1);
    }

  double em01= i1==0 ? Mdelta : sigmaRR->getDstar(i1-1)->mass0();
  double em02= i2==0 ? Mdelta : sigmaRR->getDstar(i2-1)->mass0();

  int nnyes=1;
  double faci=0.0;
  double elm=0.0;

  //...Identify type of collision.
  if(iztot == 4 || iztot == -2) {
        nnyes=0;
  } else if( iztot >= 3  ||  iztot <= -1 ) {
        faci=0.0;
        nnyes=0;

  //...d+d+->pp/d0d0->nn
  } else if(izmin == izmax) {

    faci=0.5;
    //if(kf1==kf2) faci=1.0;  // identical delta: symmetry factor chancel.
    if(isig!=0) elm = 0.4*sigmaRR->matrixElm(isig,srt,1,em01,em02);

  //...d0d++ -> pp /d+d- -> nn
  } else if( (izmin == 0  && izmax == 2)  || 
             (izmin == -1 && izmax == 1) ) {

    faci=0.5;
    if(isig!=0) elm = 0.6*sigmaRR->matrixElm(isig,srt,1,em01,em02);

  //...d0d+ -> np
  } else if(izmin == 0 && izmax == 1) {

    faci=1.0;
    if(isig!=0) elm = 0.05*sigmaRR->matrixElm(isig,srt,1,em01,em02)
                    + 0.25*sigmaRR->matrixElm(isig,srt,0,em01,em02);

  //...d-d++ -> np
  } else if(izmin == -1 && izmax == 2) {

    faci=1.0;
    if(isig!=0) elm = 0.45*sigmaRR->matrixElm(isig,srt,1,em01,em02)
                    + 0.25*sigmaRR->matrixElm(isig,srt,0,em01,em02);

  } else {
        cout << "(jamcbb3:) error in d+d collisions" << endl;
	exit(1);
  }

  double sigab=0.0;
  //...dd -> nn  from the detailed balance 
  if(nnyes == 1) {
    double prnew=sqrt(0.25*srt*srt-Mnucl*Mnucl);
    double flx=srt*srt*PCM(srt,em1,em2);
    sigab = 4*faci*prnew/flx*elm;
  }

  //double snew=2*sqrt(Mnucl*Mnucl+pr*pr);
  double sig1[10];
  sigmaRR->sigma(srt,em1,em2,sig1,1);
  int iabs=-1;

    // d- d-
    if(iztot == -2) {
	mChanel=9;
        cpair.setSizeOutGoing(mChanel);
	for(int i=0;i<9;i++)  {
	    if(DDm2[i][0]==0) sig1[i]=0.0;
	    sigin[i]=sig1[i];
	    cpair.setOutGoing(i,sig1[i],DDm2[i][0],DDm2[i][1]);
	}

    // d++ d++
    } else if(iztot == 4) {
	mChanel=9;
        cpair.setSizeOutGoing(mChanel);
	for(int i=0;i<9;i++)  {
	    if(DD4[i][0]==0) sig1[i]=0.0;
	    sigin[i]=sig1[i];
	    cpair.setOutGoing(i,sig1[i],DD4[i][0],DD4[i][1]);
	}
    // d0 d-
    } else if(iztot == -1) {
	mChanel=9;
        cpair.setSizeOutGoing(mChanel);
	for(int i=0;i<9;i++)  {
	    sigin[i]=sig1[i];
	    if(DDm1[i][0]==0) sigin[i]=0.0;
	    cpair.setOutGoing(i,sigin[i],DDm1[i][0],DDm1[i][1]);
	}
    // d+ d++
    } else if(iztot == 3) {
	mChanel=9;
        cpair.setSizeOutGoing(mChanel);
	for(int i=0;i<9;i++)  {
	    sigin[i]=sig1[i];
	    if(DD3[i][0]==0) sigin[i]=0.0;
	    cpair.setOutGoing(i,sigin[i],DD3[i][0],DD3[i][1]);
	}
    // d-d+   d0d0
    } else if(iztot == 0) {
	mChanel=18;
        cpair.setSizeOutGoing(mChanel);
	makePP(sig1);
	for(int i=0;i<16;i++) 
	cpair.setOutGoing(i,sigin[i],BBnn[i][0],BBnn[i][1]);
	cpair.setOutGoing(16,sigin[16],2112,12112); // ->n n(1440)
	cpair.setOutGoing(17,sigab,2112,2112);
	sigin[17]=sigab;
	iabs=17;

    // d-d++   d0d+
    } else if(iztot == 1) {
	mChanel=20;
        cpair.setSizeOutGoing(mChanel);
	makePN(sig1);
	for(int i=0;i<17;i++) 
	cpair.setOutGoing(i,sigin[i],BBpn[i][0],BBpn[i][1]);
	cpair.setOutGoing(17,sigin[17], 2212,12112);//pn->pn*(1440)
	cpair.setOutGoing(18,sigin[18], 2112,12212);//pn->np*(1440)
	cpair.setOutGoing(19,sigab,2212,2112);
	sigin[19]=sigab;
	iabs=19;

    // d0d++   d+d+
    } else if(iztot == 2) {
	mChanel=18;
        cpair.setSizeOutGoing(mChanel);
	makePP(sig1);
	for(int i=0;i<16;i++) 
	cpair.setOutGoing(i,sigin[i],BBpp[i][0],BBpp[i][1]);
	cpair.setOutGoing(16,sigin[16],2212,12212); // p p(1440)
	cpair.setOutGoing(17,sigab,2212,2212);
	sigin[17]=sigab;
	iabs=17;
    } else {
      cout << "(jamcbb3:) error in d+d collisions iztot= " << iztot << endl;
      exit(1);
    }

    return iabs;

}

} // end of namespace jam2
