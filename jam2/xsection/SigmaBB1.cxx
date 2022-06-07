#include <iostream>
#include <cstdlib>
#include <jam2/xsection/SigmaBB1.h>
#include <jam2/xsection/XsecTable.h>
#include <jam2/hadrons/ParticleTable.h>
#include <jam2/hadrons/JamStdlib.h>

using namespace std;

namespace jam2 {

const double SigmaBB1::eps = 1.0e-3;
const int    SigmaBB1::mxhlf=30;     // max. number of bin for integration

SigmaBB1::SigmaBB1(Pythia8::Info *inf,Pythia8::Settings *s,JamParticleData* jp,SampleMass *sm,
	Pythia8::Rndm* r) : SigmaBB(inf,s,jp,sm,r)
{

  bwR = new double [nMaxRes];

  optDeut1= 1; // JamCom::getMSTC(77);
  optDeut2= 0; // JamCom::getMSTC(78);
  optDetBal=22; // JamCom::getMSTC(62); // mstc(62)

  //optDetBal=20;

  optDetBal=21; // momentum independent resonance width.

  optSwave=1; // s-wave pion production is treated by N(1440) production

  optProb=0;

  //optDelta=1;   // NN -> N Delta from matrix element
  optDelta=0; // JAM1 original

  sigNDelta = new SigmaNR(2);


}
SigmaBB1::~SigmaBB1()
{
  delete [] bwR;
}

//***********************************************************************

void SigmaBB1::calc(CollisionPair& cpair)
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
	cout << "SigmaBB1::calc kf<0? kf1= "<< kf1 << " kf2= " << kf2<<endl;
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

	    cout << "SigmaBB1:: BB collision something is wrong PID id1= " << id1
		 << " id2= " << id2 
		 << endl;
	    exit(1);

	}

        if(iabs>=0 && absorptionBB) sigab=sigin[iabs];
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


    //if(snew <= 3.6) {
    if(srt <= 3.6) {
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
void SigmaBB1::sigmaNuclNucl(CollisionPair& cpair)
{
  const double srtStr=3.2;   // threthold of string excitation.
    // First compute sig and sigel from cross section table.
    jamxnn(srt,iz1,iz2);

    cpair.setSigma(sigTot,sigEl);

    //...Resonance cross section alone cannot fill the inelastic cross
    //...section at high energy, we use parametrized sig and sigel.
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

    if(sigsoft>0) {
        mChanel++;
	cpair.setOutGoing(sigsoft,92,92);
    }

}

//============================================================================
// Compute PN -> X
// 1) NN -> ND  2) NN -> NN*  3) NN -> DD  4) NN -> ND*  5) NN -> N*D
// 6) NN -> DD*  7) NN -> N*N*  8) NN -> N*D*  9) NN -> D*D* 
void SigmaBB1::sigmaPN(CollisionPair& cpair)
{
//..First compute T=1 cross section.
    double sig1[10],sig0[10];
    jamxnnin(srt,sig1,1,optDeut1);

//...deuteron + pi+ cross section.
    double sigdpi=0.0;
    if(optDeut1>= 2) sigdpi = jamxdpi1(srt);

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

    sigin[19]=0.0;
    sigin[20]=0.0;
    sigin[21]=0.0;
    mChanel=22;
    cpair.setSizeOutGoing(mChanel);
    cpair.setOutGoing(19,0.0,0,0);
    cpair.setOutGoing(20,0.0,0,0);
    cpair.setOutGoing(21,0.0,0,0);

    //...Compute T=0 cross section.
    jamxnnin(srt,sig0,0,optDeut1);

    sigin[0]=0.25*sig1[0];        // nd+
    sigin[1]=0.25*sig1[0];        // pd0
    if(optDeut1 == 2) {
        sigin[0] += 0.25*sigdpi;
        sigin[1] += 0.25*sigdpi;
    } else if(optDeut1 >= 3) {
        sigin[19]=0.5*sigdpi;       // deuteron + pi
    }
    if(optDeut2 >= 1) {
	sigin[20]=jamxdpi2(2,srt);
        sigin[21]=jamxdpi2(3,srt);
    }

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

    for(int i=0;i<17;i++)  {
      cpair.setOutGoing(i,sigin[i],BBpn[i][0],BBpn[i][1]);

      /*
      cout << "PN sig= " << i << " sigin= " << sigin[i]
	  << " id1= " << BBpn[i][0]
	  << " id2= " << BBpn[i][1]
	  << " sig= " << cpair.getOutGoing()[i].sig
	  <<endl;
	  */
    }

    if(optSwave==1) {
	sigin[17]=0.5*sig0[9];
        sigin[18]=0.5*sig0[9];
	cpair.setOutGoing(17,sigin[17], 2212,12112);//pn->pn*(1440)
	cpair.setOutGoing(18,sigin[18],12212, 2112);//pn->p*(1440)n
    } else  {
        sigin[17]=0.5*sig0[9]*1.0/3.0;                 // s-wave
        sigin[18]=0.5*sig0[9]*2.0/3.0;                 // s-wave
        sigin[19]=0.5*sig0[9]*1.0/3.0;                 // s-wave
        sigin[20]=0.5*sig0[9]*2.0/3.0;                 // s-wave
	cpair.setOutGoing(17,sigin[17],2212,2112,111);//pn->n(p pi0)
	cpair.setOutGoing(18,sigin[18],2112,2112,211);//pn->n(n pi+)
	cpair.setOutGoing(19,sigin[19],2212,2112,111);//pn->p(n pi0)
	cpair.setOutGoing(20,sigin[20],2212,2212,-211);//pn->p(p pi-)
    }

    /*
    double sigint=0.0;
    for(int i=0;i<mChanel;i++)  {
	sigint += sigin[i];
    }

    jamxnn(srt,iz1,iz2);
    cout << " sig= " << sigTot << " sigel= " << sigEl
	<< " sig-sigel=" << sigTot - sigEl
	<< " sigintot= " << sigint
	<<endl;
	*/

}

//============================================================================
// Compute PP/NN -> X
// 1) NN -> ND  2) NN -> NN*  3) NN -> DD  4) NN -> ND*  5) NN -> N*D
// 6) NN -> DD*  7) NN -> N*N*  8) NN -> N*D*  9) NN -> D*D* 
void SigmaBB1::sigmaPP(CollisionPair& cpair)
{
//.....Resonance production cross sections.
    double sig1[10];
    jamxnnin(srt,sig1,1,optDeut1);


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

    mChanel=18;
    sigin[17]=0.0;
    sigin[18]=0.0;
    sigin[0]=0.25*sig1[0];   // pd+
    sigin[1]=0.75*sig1[0];   // nd++

//...deuteron + pi+ cross section.
    double sigdpi=0.0;
    if(optDeut1>= 2) sigdpi = jamxdpi1(srt);
    if(optDeut1 == 2) sigin[1] += sigdpi;
    else if(optDeut1 >= 3) sigin[17]=sigdpi;

    if(optDeut2 >= 1) sigin[18]=jamxdpi2(1,srt);

    sigin[2]=sig1[1];       // pp*
    sigin[3]=0.4*sig1[2];   // d+d+
    sigin[4]=0.6*sig1[2];   // d0d++
    sigin[5]=0.25*sig1[3];  // pd*+
    sigin[6]=0.75*sig1[3];  // nd*++
    sigin[7]=0.25*sig1[4];  // p*d+
    sigin[8]=0.75*sig1[4];  // n*d++
    sigin[9]=0.4*sig1[5];  // d*d+
    sigin[10]=0.6*sig1[5];  // d*d++
    sigin[11]=sig1[6];        // p*p*
    sigin[12]=0.25*sig1[7]; // p*d*+
    sigin[13]=0.75*sig1[7]; // n*d*++
    sigin[14]=0.4*sig1[8];  // d*+d*+
    sigin[15]=0.6*sig1[8];  // d*0d*++

    cpair.setSizeOutGoing(mChanel);

    // nn collision
    if(iz1+iz2==0) {

	 for(int i=0;i<16;i++) 
	 cpair.setOutGoing(i,sigin[i],BBnn[i][0],BBnn[i][1]);
	 if(optSwave==1) {
	   sigin[16]=sig1[9];
	   cpair.setOutGoing(16,sigin[16],2112,12112); // n n(1440)
	 } else {
           sigin[16]=sig1[9]/3.0;    // s-wave pp->pp pi0
           sigin[17]=sig1[9]*2.0/3.0; // s-wave pp->pn pi+
	   cpair.setOutGoing(16,sigin[16],2112,2112,111);
	   cpair.setOutGoing(17,sigin[17],2212,2112,-211);
	 }

    // pp collision
    } else {

	 for(int i=0;i<16;i++) 
	 cpair.setOutGoing(i,sigin[i],BBpp[i][0],BBpp[i][1]);
	 if(optSwave==1) {
	   sigin[16]=sig1[9];
	   cpair.setOutGoing(16,sigin[16],2212,12212); // p p(1440)
	 } else {
           sigin[16]=sig1[9]/3.0;    // s-wave pp->pp pi0
           sigin[17]=sig1[9]*2.0/3.0; // s-wave pp->pn pi+
	   cpair.setOutGoing(16,sigin[16],2212,2212,111);
	   cpair.setOutGoing(17,sigin[17],2212,2112,211);
	 }
    }

}

//***********************************************************************

void SigmaBB1::jamxnnin(double x,double *siginel,int iporn,int ideut)
{
//...Output NN inelastic resonance production cross section  || given x.
//...(INPUT)
//     x: c.m.enrgy in GeV
//     i || :
//        = 1: t=1 cross sections  || pp
//        = 2: t=1 cross sections  || pn:
//               only deuteron contribution is different
//        = 0: t=0 crosss sections
//     ideut switch if deuteron pi is included || not
//   (OUTPUT)
//     siginel: inelastic cross sections
//
//...1) NN -> ND  (deuteron production is included)
//...2) NN -> NN*
//...3) NN -> DD
//...4) NN -> ND*
//...5) NN -> N*D
//...6) NN -> DD*
//...7) NN -> N*N*
//...8) NN -> N*D*
//...9) NN -> D*D*
//...10)NN-> s-wave pion

//...NN->ND(1232)
    static const double rppd0 = 0.633719, gppd0 = 2.11477, sppd0 = 0.0171405,
	   thd0 = 2.0139999;

    static const double  sgmppd0=48.64510, sppd1=2.06672,fppd1=0.480085,
		 gppd1=0.576422, thd1=2.124;
      
//...T=1
    double a1[9][5]={
     {  0.00000,  0.00000,  0.00000, 0.000000, 0.000},
     { 24.94700,  2.48150,  2.63330, 0.425358, 2.162},
     {  7.63181,  1.41140,  2.67784, 0.311722, 2.252},
     {  8.01615,  2.74161,  3.34503, 0.259703, 2.340},
     { 13.14580,  2.06775,  2.75682, 0.247810, 2.300},
     { 19.63220,  2.01946,  2.80619, 0.297073, 2.528},
     { 11.67320,  2.31682,  2.96359, 0.259223, 2.438},
     {  2.99086,  2.29380,  3.54392, 0.090438, 2.666},
     { 35.13780,  2.25498,  3.14299, 0.215611, 2.804}};
//...T=0
    double a0[9][5]={
     {  0.00000,  0.00000,  0.00000, 0.000000, 0.000},
     {166.60600,  2.10128,  2.34635, 0.284955, 2.162},
     { 39.99770,  1.83576,  2.40348, 0.288931, 2.252},
     {  0.00000,  0.00000,  0.00000, 0.000000, 0.000},
     {  0.00000,  0.00000,  0.00000, 0.000000, 0.000},
     { 56.32490,  2.00679,  2.71312, 0.362132, 2.528},
     {  2.14575,  0.21662,  3.40108, 0.252889, 2.438},
     {  0.00000,  0.00000,  0.00000, 0.000000, 0.000},
     {  4.14197,  1.67026,  3.75133, 0.476595, 2.804}};

//...s-wave pion production pp
    double as1[5]={15.644100, 1.675220,  2.07706, 0.658047, 2.014};
//...s-wave pion production pn
    double as0[5]={78.868103, 0.746742,  1.25223, 0.404072, 2.014};

    for(int i=0;i<10;i++) siginel[i]=0.0;

//...I=1 cross sections
    if(iporn == 1 || iporn== 2) {

//...Low energy.
    if(x<=thd0) return;

//...NN->ND(1232)
    double sigd;

    if(optDelta==1) {

      sigd=sigNDelta->sigma(x);

    // JAM1 parametrization.
    } else {

    double sigd0=0.0;
    if(x > thd0) sigd0 = sgmppd0*rppd0*sqrt((x-thd0)/thd0)*gppd0
              /((x-sppd0)*(x-sppd0)+gppd0*gppd0)/100;

    double sigd1=0.0;
      if(x > thd1)
       sigd1=sgmppd0*pow(x/thd1-1,fppd1)*bw(x/sppd1,gppd1);

      sigd=sigd1+sigd0; // pp->ND  (-> NNpi)

    }

//...Deuteron pi.
    double  sigdeut=0.0;
    if(ideut == 1) sigdeut = jamxdpi1(x);

//...[1]NN -> ND,dpi
    siginel[0] = sigd + sigdeut;
      
//      do i=2,9
    for(int i=1;i<9;i++) {
       if(x > a1[i][4])
       siginel[i]=a1[i][0]*pow((x/a1[i][4]-1),a1[i][1])*bw(x/a1[i][2],a1[i][3]);
    }

//...[10]s-wave pion
    if(x > as1[4])
     siginel[9]=as1[0]*pow((x/as1[4]-1),as1[1])*bw(x/as1[2],as1[3]);


//...T=0 cross sections.
    } else {

    for(int i=1;i<9;i++) {
       if(a0[i][4] > 2.0 && x > a0[i][4])
        siginel[i]=a0[i][0]*pow((x/a0[i][4]-1),a0[i][1])*bw(x/a0[i][2],a0[i][3]);
    }

//...[10] s-wave pion
    if(x > as0[4])
      siginel[9]=as0[0]*pow((x/as0[4]-1),as0[1])*bw(x/as0[2],as0[3]);

    }

}

//***********************************************************************
double SigmaBB1::jamxdpi1(double x)
{

//...Purpose: to give pp->deuteron pi+ cross section.
//     x: c.m.enrgy in GeV

//...deuteron 1
    double ad1[5]={0.14648,0.20807, 2.13072, 0.042475, 2.024};
//...deuteron 2
    double ad0[5]={0.12892,0.08448, 2.18138, 0.059207, 2.054};

//...Deuteron pi.
    double sigdeut1=0.0;
    if(x > ad1[4])
        sigdeut1 = ad1[0]*pow((x/ad1[4]-1),ad1[1])*bw(x/ad1[2],ad1[3]);

    double sigdeut2 = 0.0;
    if(x > ad0[4])
         sigdeut2 = ad0[0]*pow((x/ad0[4]-1),ad0[1])*bw(x/ad0[2],ad0[3]);

//...pp -> deuteron pi+  = 2(pn -> d + pi0)
    return sigdeut1+sigdeut2;

}

//***********************************************************************
double SigmaBB1::jamxdpi2(int iopt,double x)
{
//...Purpose: to give double pionic fusion NN -> deuteron pi + pi cross section.
//     x: c.m.enrgy in GeV
//...Ref. arXiv:1212.2881[nucl-ex]

//....iopt=1  pp -> d pi+ pi0
//....iopt=2  pn -> d pi0 pi0
//....iopt=3  pn -> d pi+ pi-

//...pp -> d pi+ pi0
    double ad1[5]={2.15249,0.0149338,1.12574,2.49005,0.091389};

//...2(pn -> d pi0 pi0)
    double ad2[5]={2.1479,0.0118512,1.74826,2.34993,0.0401094};


    double sig=0.0;
    if(iopt == 1 || iopt == 3) {
	if(x > ad1[0]) sig = 
	    ad1[1]*pow((x-ad1[0]),ad1[2])/((x-ad1[3])*(x-ad1[3])+ad1[4]*ad1[4]);
        if(iopt == 1) return sig;
    }

    double sig2=0.0;
    if(x  > ad2[0]) sig2 = 
	ad2[1]*pow((x-ad2[0]),ad2[2])/((x-ad2[3])*(x-ad2[3])+ad2[4]*ad2[4]);
    if(iopt == 2) return 0.5*sig2;
     
//...pn -> d pi+ pi-  = 2*(pn->d pi0pi0)+ 1/2(pp->dpi+pi0)
    return  sig2 + 0.5*sig;

}

//**********************************************************************
//  select N* or D* resonance.
Pythia8::ParticleDataEntry* SigmaBB1::selectDN(int idn, int& idr)
{
   if(idn==2212)  {
       idr=id_nucl;
        return jamTable->getProton();
   } else if (idn==2112) {
       idr=id_nucl;
        return jamTable->getNeutron();
   } else if(idn==1114) {
        idr=id_delt;
        return jamTable->getDeltam();
   } else if(idn==2114) {
        idr=id_delt;
        return jamTable->getDelta0();
   } else if(idn==2214) {
        idr=id_delt;
        return jamTable->getDeltap();
   } else if (idn==2224) {
        idr=id_delt;
        return jamTable->getDeltapp();
   } else if (idn==12112) {
        idr=id_nucls;
        return jamTable->getN1440();
   } else if (idn==12212) {
        idr=id_nucls;
        return jamTable->getP1440();
   } else {
    if(idn==1) {idr=id_nucls; return selectResonance(nStar[0]);}
    if(idn==2) {idr=id_nucls; return selectResonance(nStar[1]);}
    if(idn==3) {idr=id_delts;return selectResonance(dStar[0]);}
    if(idn==4) {idr=id_delts;return selectResonance(dStar[1]);}
    if(idn==5) {idr=id_delts;return selectResonance(dStar[2]);}
    if(idn==6) {idr=id_delts;return selectResonance(dStar[3]);}
   }
    cout << " SigmaBB::selectDN wrong id = " << idn <<endl;
    exit(1);
}


Pythia8::ParticleDataEntry* SigmaBB1::selectResonance(ParticleTable* table)
{
    double bwtot=0.0;
    int n=table->size();
    vector<double>  bw(n);
    for(int i=0;i<n;i++) { 
	bw[i]=0.0;
        ParticleDataEntry* p=table->getParticle(i);
	// skip delta(1232).
	int idn=p->id();
        if(idn==1114 || idn==2114 || idn==2214 || idn==2224) continue;
        bw[i] = p->spinType();
        bwtot += bw[i];
    }
    double x=rndm->flat()*bwtot;
    for(int i=0;i<n;i++) {
        x -= bw[i];
        if(x<0) {
            return table->getParticle(i);
        }   
    }   
    return table->getParticle(n-1);
}

bool SigmaBB1::sampleResonanceMass(double ecm, int idn1, int idn2, double& m1, double& m2)
{    
  double emr[2],em0[2],gam0[2];
  double mmax[2],mmin[2],bwmax[2],bw[2],ymin[2],ymax[2];
  int iex[2], id[2],idv[2];

  vector<ParticleDataEntry*> table1 = pickRTable(idn1,idv[0]);
  vector<ParticleDataEntry*> table2 = pickRTable(idn2,idv[1]);
  int n1=table1.size();
  int n2=table2.size();
  if(n1==0 || n2==0) {
      cout << " no baryon table ? n1= "<< n1
	  << " n2= "<< n2<<endl;
      exit(1);
  }

  if(idv[0]==id_nucl && idv[1]==id_nucl) {
      //cout << "id1= "<< idn1 << " id2= "<< idn2 << " m1= "<< m1
      //	  << " m2= " << m2 <<endl;
    m1 = table1[0]->m0();
    m2 = table2[0]->m0();
    pout[0] = table1[0];
    pout[1] = table2[0];
    return true;
  }

  double totsp=0.0;
  vector<vector<double> > probr(n1, vector<double>(n2,0.0));
  for(int i=0;i<n1;i++)
  for(int j=0;j<n2;j++) { 
    probr[i][j]=0.0;
    ParticleDataEntry* pv[2];
    pv[0]=table1[i];
    pv[1]=table2[j];
    double pbw=0.0;
    if(optProb==0) pbw = sampleMass->probBW2(ecm,pv[0],pv[1]);
    else pbw = sampleMass->probBW1(ecm,pv[0],pv[1]);
    if(pbw <= 0.0) continue;
    int ispin1=pv[0]->spinType();
    int ispin2=pv[1]->spinType();
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
  cout << " SigmaBB1::sampleResonanceMass no particle? "<< srt << endl;
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
    cout << "SigmaBB1::ResnanceMass too small srt? srt= " << ecm <<endl;
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
  if(idv[0]==id_nucl && idv[1]==id_nucl) {
      cout << "id1= "<< idn1 << " id2= "<< idn2 << " m1= "<< m1
	  << " m2= " << m2 <<endl;

    cout << " srt= "<< ecm << " srt= "<< srt << " n1= "<< n1
	<< " n2= "<< n2 <<endl;
  cout << "id1= "<< idn1 << " id2= "<< idn2
      << " -> id3= "<< id[0] << " m3= "<< m1 << " id4= "<< id[1]
      << " m4= "<< m2
      <<endl;
  cin.get();
  }
  */

    return true;
  }

  } while (mtry++ < 100);

  cout << "SigBB::sampelResonanceMass does not converge srt=" << ecm
       << " mtry= "<< mtry
       << " id1= " << id[0] << " id2= " << id[1]
       << " em1= " << emr[0] << " em2= " << emr[1]
       <<endl;

  return sampleMassFix(ecm,pf0,iex,emr,em0,gam0,mmin,ymin,ymax,m1,m2);

}

// compute probability to select N* or D* resonance.
double SigmaBB1::jambres2(ParticleTable* table, double mr, int id)
{
    double bwtot=0.0;
    int n=table->size();

    int ir=-1;
    //vector<double> bw(n);
    //double *bw = new double[n];
    for(int i=0;i<n;i++) {
	ParticleDataEntry* p=table->getParticle(i);
	double twid=decay->getTotalWidth(p,mr);        
        double  em0=p->m0();
	int spin=p->spinType();
	int idr=p->id();
	if(id==idr) ir=i;
	bwR[i]=spin*twid/((mr-em0)*(mr-em0) + 0.25*twid*twid);
	bwtot += bwR[i];
    }
    if(bwtot==0.0 || ir < 0) {
	cout << "SigmaBB1::jambres2 id? " << id << " mr= " << mr 
	     << " id= " << id 
	     << " n= " << table->size()
	     <<endl;

	for(int i=0;i<table->size();i++) {
	    ParticleDataEntry* p=table->getParticle(i);
	    //double  em0=p->m0();
	    int idr=p->id();
	    cout << i << " name= " << p->name() << " idr= " << idr
		 <<endl;
	}
	exit(1);
    }

    return bwR[ir]/bwtot;
}

int SigmaBB1::sigmaNsNs(CollisionPair& cpair)
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
    int iabs=-1;
    double detbal=1.0;
    double prob1=1.0, prob2=1.0;
    int isig;
    // particle 1 is nucleon?
    if( kf1 == 2112  ||  kf1 == 2212 ) {
        isig=2;
        detbal = DetBal1(1,pd2,em1,em2);
	prob2=jambres2(nStar[iz2/3],em2,kf2);
    // particle 2 is nucleon?
    } else if( kf2 == 2112  ||  kf2 == 2212 ) {
        isig=2;
        detbal = DetBal1(1,pd1,em1,em2);
	prob1=jambres2(nStar[iz1/3],em1,kf1);
    } else {
        detbal = DetBal2(1);
	 if(iz1<0 || iz2<0) {
		cout << "5 before jambres2 iz1= " << iz1/3 <<endl;
		cout << "5 before jambres2 iz2= " << iz2/3 <<endl;
		cout << " kf1= " << kf1 << " kf2= " << kf2<<endl;
		cout << " id1= " << id1 << " id2 " << id2 <<endl;
	 }
	prob1=jambres2(nStar[iz1/3],em1,kf1);
	prob2=jambres2(nStar[iz2/3],em2,kf2);
        isig=6;
    }

    double sig1[10],sig0[10];
    double sigab;
    int icc;
    double snew=2*sqrt(Mnucl*Mnucl+pr*pr);

//...Get cross sections from nn collision.
    double faci;
    if( iz1  ==  iz2 ) {
        icc=1;
	if(iz1+iz2==0) icc=2;
        faci=0.5;
	if(kf1==kf2) faci=1.0;
        jamxnnin(srt,sig1,1,optDeut1);
        sigab=sig1[isig];
        jamxnnin(snew,sig1,1,optDeut1);
    } else {
        icc=0;
        faci=1.0;
        jamxnnin(srt,sig1,2,optDeut1);
        jamxnnin(srt,sig0,0,optDeut1);
        sigab=0.25*(sig0[isig]+sig1[isig]);
        jamxnnin(snew,sig1,2,optDeut1);
        jamxnnin(snew,sig0,0,optDeut1);
    }

//...Spin factor.
    double spin= (abs(kf1) %10) *  (abs(kf2) % 10);
    double facsp=4.0/spin;

//...Calculate N(*)N(*) --> NN absorption cross section
//...from detailed balance.
    double pr2new=0.25*srt*srt-Mnucl*Mnucl;
    sigab *= facsp*faci*pr2new*detbal * prob1 * prob2;

    if(icc >= 1) {
	mChanel=18;
        cpair.setSizeOutGoing(mChanel);
        sigin[0]=0.25*sig1[0]; // pd+
        sigin[1]=0.75*sig1[0]; // nd++
        sigin[2]=sig1[1];        // pp*+
        sigin[3]=0.4*sig1[2];  // d+d+
        sigin[4]=0.6*sig1[2];  // d0d++
        sigin[5]=0.25*sig1[3]; // pd*+
        sigin[6]=0.75*sig1[3]; // nd*++
        sigin[7]=0.25*sig1[4]; // p*d+
        sigin[8]=0.75*sig1[4]; // n*d++
        sigin[9]=0.4*sig1[5]; // d*d+
        sigin[10]=0.6*sig1[5]; // d*d++
        sigin[11]=sig1[6];       // p*p*
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
	cpair.setOutGoing(17,sigin[17], 2212,12112);//pn->pn*(1440)
	cpair.setOutGoing(18,sigin[18],12212, 2112);//pn->p*(1440)n
	cpair.setOutGoing(19,sigin[19],2212, 2112);//->pn
    }

    return iabs;

}

//***********************************************************************

int SigmaBB1::sigmaNsDs(CollisionPair& cpair)
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

    double sig1[10];

    int iz01=iz1/3;
    int iz02=iz2/3;
    int itz=iz01+iz02;

//...2 is a nucleon(*) 1 is a delta(*)
    int kfn=kf2;
    int kfd=kf1;
    int izn=iz02;
    int izd=iz01;
    int iresn=1;
    int iresd=1;
    if(id1 == id_delt) iresd=0;
    if(id2 == id_nucl) iresn=0;
    Pythia8::ParticleDataEntry *kcd=pd1;
    double emn=em2;
    double emd=em1;

    //...1 is a nucleon(*) 2 is a delta(*)
    if(id1 == id_nucl  ||  id1 == id_nucls) {
	kfn=kf1;
	kfd=kf2;
        izn=iz01;
        izd=iz02;
        iresn=1;
        iresd=1;
        if(id1 == id_nucl) iresn=0;
        if(id2 == id_delt) iresd=0;
        kcd=pd2;
        emn=em1;
        emd=em2;
    }

    int idn=0;
    double faci=0.0;
    int    icc=1;
    double detbal=0.0;
    int isig;

    if(itz > -1 && itz < 3) {

	if(iresn == 0) {
	    if(iresd == 0) {
		idn = 1;   //   delta n ==> n n
		isig=0;
		detbal=DetBal1(1,kcd,emd,emn) * normD1232;

	    //.....delta* n ==> n n
	    } else {
		idn=2;
		isig=3;
		detbal=DetBal1(1,kcd,emd,emn);
		detbal *= jambres2(dStar[izd+1],emd,kfd);
	    }
	} else {

	    //.....d(d*) n* ==> n n
	    detbal=DetBal2(1);
	    detbal *= jambres2(nStar[izn],emn,kfn);
	    if(iresd == 0) {
		idn=3;   //   delta n* ==> n n
		isig=4;
	    } else {
		idn=4;   //   delta* n* ==> nn 
		isig=7;
	        detbal *= jambres2(dStar[izd+1],emd,kfd);
	    }
	}
    }

    double sigab=0.0;
    int iabs=-1;

    if(idn == 0) {
        sigab=0;
//....pd- --> nn  16:nd++ --> pp
    } else if( (izn == 1 && izd == -1)  ||  
               (izn == 0 && izd == 2)  )  {
        faci=0.5;
        icc=1;
        jamxnnin(srt,sig1,1,optDeut1);
        sigab=0.75*sig1[isig];

//....nd0- --> nn  12:pd+ --> pp
    } else if( (izn == 0 && izd == 0)  ||  
                (izn == 1 && izd == 1)  ) {

        faci=0.5;
        icc=1;
        jamxnnin(srt,sig1,1,optDeut1);
        sigab=0.25*sig1[isig];

//....pd0 --> pn  11:nd+ --> pn
    } else if( (izn == 1 && izd == 0)  ||  
                (izn == 0 && izd == 1)  ) {

        faci=1.0;
        icc=2;
        jamxnnin(srt,sig1,2,optDeut1);
        sigab=0.25*sig1[isig];
    }

//...Get absorption cross sections
      if(idn != 0) {
//       if(vfd.le.0.0d0) vfd=1.0
//       sigab=facsp*faci*vfd*pr2new*sigab/(pr*pr)
	int kfa1=abs(kf1);
	int kfa2=abs(kf2);
        double facsp=4.0/double((kfa1%10)*(kfa2%10));
        double pr2new=0.25*srt*srt-Mnucl*Mnucl;
        sigab=facsp*faci*pr2new*sigab*detbal;
      }

//...Calculate nn cross section.
    double snew=2*sqrt(Mnucl*Mnucl+pr*pr);
    jamxnnin(snew,sig1,icc,optDeut1);
    //sigin[10]=sigab;
    //for(int i=0;i<10;i++) sigin[i]=sig1[i];

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
	cpair.setOutGoing(18,sigin[18],12212, 2112);//pn->p*(1440)n
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
	sigin[17]=sigab;
	cpair.setOutGoing(17,sigab,2212,2212);
	iabs=17;

    }


//...Include s-wave pion into NN* branch
/*
      sigin[1] += sig1[9];
      sigin[9]=0.0;     // assume no s-wave pion
      if(idn == 0) {
        sigin[1]=0.0;
        sigin[6]=0.0;
      }
      mChanel=11;
      */


    return iabs;
}

//**********************************************************************
//...Purpose : to give cross sections  || d(*) + d(*) collisions.
//...at the moment dn -> x  cross sections are simply equated to
//...t=1 cross sections at equal cms momentum
int SigmaBB1::sigmaDD(CollisionPair& cpair)
{
    double sig1[10],sig0[10];

    int iz01=iz1/3;
    int iz02=iz2/3;
    int iztot=iz01+iz02;
    int izmin=min(iz01,iz02);
    int izmax=max(iz01,iz02);

    int isig;
    if(id1 == id_delt && id2 == id_delt) {
        isig=2;
    } else if( (id1 == id_delt  && id2 == id_delts)
            || (id1 == id_delts && id2 == id_delt) ) {
        isig=5;
    } else if(id1 == id_delts && id2 == id_delts) {
        isig=8;
    } else {
	cout << "(jamcbb3:) error in d+d collisions "
	<< " id1= " << id1 << " id2= " << id2 << endl;
	exit(1);
    }

    int nnyes=1;
    int ndyes=1;
    double faci=0.0;
    double sigab=0.0;
    //int    icc=1;

    //...Identify type of collision.
    if(iztot == 4 || iztot == -2) {
        nnyes=0;
        ndyes=0;
	//mChanel=0;
        //return;
    } else if( iztot >= 3  ||  iztot <= -1 ) {
        faci=0.0;
        //icc=1;
        nnyes=0;

    //...d+d+->pp/d0d0->nn
    } else if(izmin == izmax) {

        faci=0.5;
	if(kf1==kf2) faci=1.0;
        //icc=1;
        jamxnnin(srt,sig1,1,optDeut1);
        sigab=0.4*sig1[isig];

    //...d0d++ -> pp /d+d- -> nn
    } else if( (izmin == 0  && izmax == 2)  || 
               (izmin == -1 && izmax == 1) ) {

        faci=0.5;
        //icc=1;
        jamxnnin(srt,sig1,1,optDeut1);
        sigab=0.6*sig1[isig];

    //...d0d+ -> np
    } else if(izmin == 0 && izmax == 1) {

        faci=1.0;
        //icc=2;
        jamxnnin(srt,sig1,2,optDeut1);
        jamxnnin(srt,sig0,0,optDeut1);
        //sigab=0.025*sig1[isig]+0.125*sig0[isig];
        sigab=0.05*sig1[isig]+0.25*sig0[isig];

    //...d-d++ -> np
    } else if(izmin == -1 && izmax == 2) {

        faci=1.0;
        //icc=2;
        jamxnnin(srt,sig1,2,optDeut1);
        jamxnnin(srt,sig0,0,optDeut1);
        sigab=0.45*sig1[isig]+0.25*sig0[isig];

    } else {
        cout << "(jamcbb3:) error in d+d collisions" << endl;
	exit(1);
    }


    //...dd -> nn  from the detailed valance 
    if(nnyes == 1) {
        double detbal=DetBal2(1);
	if(id1==id_delts)  {
	    if(iz01+1<0 || iz01+1>3) {
		cout << "1 before jambres2 iz01= " << iz01 <<endl;
		exit(1);
	    }
	    detbal *= jambres2(dStar[iz01+1],em1,kf1);
	}
	if(id2==id_delts)  {
	    if(iz02+1<0 || iz02+1>3) {
		cout << "2 before jambres2 iz01= " << iz01 <<endl;
		exit(1);
	    }
	    detbal *= jambres2(dStar[iz02+1],em2,kf2);
	}

	int kfa1=abs(kf1);
	int kfa2=abs(kf2);
        double facsp=4.0/double((kfa1%10)*(kfa2%10));
        double pr2new=0.25*srt*srt-Mnucl*Mnucl;
        sigab=facsp*faci*pr2new*sigab*detbal;
    } else {
	sigab=0.0;
    }

    sigin[10]=sigab;
    double snew=2*sqrt(Mnucl*Mnucl+pr*pr);
    jamxnnin(snew,sig1,1,optDeut1);
    // why this was commented out?
    for(int i=0;i<10;i++) sigin[i]=sig1[i];

//...Include s-wave pion into NN* branch
    sigin[1] += sig1[9];
    sigin[9] = 0.0;  // assume no s-wave pion 
    if(nnyes == 0) {
        sigin[1]=0.0;
        sigin[6]=0.0;
    }
    if(ndyes == 0) {
        sigin[0]=0.0;
        sigin[3]=0.0;
        sigin[4]=0.0;
        sigin[7]=0.0;
    }
    mChanel=11;

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
	cpair.setOutGoing(18,sigin[18],12212, 2112);//pn->p*(1440)n
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
    }

    return iabs;

}

//***********************************************************************
//...Purpose to calculate the correction factor of BR->BN cross section.
//...Integration by using Simpson (1/3) rule.

double SigmaBB1::DetBal1(int msel,Pythia8::ParticleDataEntry* p,double m1,double m2)
{
//...m1: ingoing resonance mass.
//...m2: ingoing particle mass (not resonance assumed).
//---------------------------------------------------------------------
//...isw1=1: non-rel
//...isw1=2: rel
//...isw1=3: rel2
//...isw2=1: no momentum dependence for integrand.
//...isw2=2: momentum dependence
//...isw2=3: momentum dependence (P.Danielewicz,G.F.Bertch)
//...isw2=4: momentum-squeard dependence
//---------------------------------------------------------------------
    static const double emnuc=0.939,empion=0.138;

    int isw=optDetBal;
    int isw1=(isw/10)%10;
    int isw2=isw%10;
    if((isw1 < 1 || isw1 > 3) || (isw2 < 0 || isw2 > 4)) {
        cout << " (DetBal1:)invalid isw " << isw <<endl;
	exit(1);
    }

    double wid=p->mWidth(); // constant width.
    if(wid <= 1e-5) {
        return 1.0;
    }
    double emr0 = p->m0();
    double emin = p->mMin();
    if(emin==0.0) emin=p->m0();
    double emax=srt-m2-0.001;
    int id=p->id();
    if(id == id_nucls)
        emin=emnuc+2*empion+eKinMin;
    else if(id == id_delt)
        emin=emnuc+empion+eKinMin;
    else if(id == id_delts)
        emin=emnuc+3*empion+eKinMin;

    if(emax <= emin) {
        if(msel == 1) return 1.0/(pr*pr);
        if(msel == 2) return 0.0;
    }


//...Option for constant width and no final momentum dependence.
    if(isw2 <= 0) {
        double den1=atan((emax*emax-emr0*emr0)/(wid*emr0));
        double den2=atan((emin*emin-emr0*emr0)/(wid*emr0));
        return M_PI/max(den1-den2,1.e-4)/(pr*pr);
    }

//....Initialization for integration.
    double a=emin;
    double b=emax;
    int jmax=0;
    int nhlf=0;
    double  h=b-a;
    double s1=BW1(isw,p,a,emr0,m2)+BW1(isw,p,b,emr0,m2);
    double s2=0.0;
    double s4=0.0;
    double sum=0.5*h*s1;
    double ds=100;

//...Start the Simpson method.
    do {
      jmax = 2*jmax+1;
      nhlf++;
      h *=0.5;
      s2 += s4;
      s4=0.0;
      for(int j=1;j<=jmax;j+=2) {
      //do 110 j=1,jmax,2
       s4 += BW1(isw,p,a+j*h,emr0,m2);
      }
//110   continue
      double ss=sum;
      sum=(h/3.0)*(s1+2.0*s2+4.0*s4);
      ds=abs(sum-ss);
      if(nhlf == mxhlf) {
        cout << "simpson method does not converge. mxhlf=" << mxhlf <<endl;
        nhlf=-nhlf;
	break;
      }
    } while (ds > eps*abs(sum));


    if(msel ==  2) return sum/M_PI;

    double fac=1.0/max(sum/M_PI,1e-5);
    if(isw2 == 1) {
        fac *= 1.0/(pr*pr);
    } else if(isw2 == 2 || isw2 == 3) {
        fac *= 1.0/pr;
        if(isw2 == 3) fac *= (m1/emr0);
    }

    return fac;

}

//***********************************************************************

double SigmaBB1::BW1(int idetsw,Pythia8::ParticleDataEntry* p,double em,double mr0,double m2)
{
//...Options.
//...isw:
//...    =1  : no final-momentum dependence for integrand.
//...    =2,3: final-momentum dependence.
//...    =4  : final-momentum-squeard dependence.

    double bwf1=0.0;
    int isw= idetsw % 10;
    double gam;
    //if(idetsw >= 10)
    if(isw >= 2) {
        gam=decay->getTotalWidth(p,em);
    } else {
        gam = p->mWidth();
	if(isDelta(p->id())) gam=decay->getDeltaWidth(3,em);
    }

//...Non-rel
    if(idetsw <= 20) {
        bwf1=gam*0.5/((em-mr0)*(em-mr0)+0.25*gam*gam);
//...Rel1
    } else if(idetsw <= 30) {
	double er=em*em - mr0*mr0;
        bwf1=2*em*mr0*gam/(er*er+(mr0*gam)*(mr0*gam));
//...Rel2
    } else if(idetsw <= 40) {
	double er=em*em-mr0*mr0;
        bwf1=2*em*em*gam/(er*er+(em*gam)*(em*gam));
    }

    if(isw == 2 || isw == 3) {
        double pf;
        if(srt > em+m2+0.00001) {
	  pf=XsecTable::pawt(srt,em,m2);
        } else {
          cout << "(jambwf1:)srt<em1+m2 srt= " << srt 
	      << " em= " << em << " m2= " << m2
	      << " em+m2+eps = " << srt - (em+em2+eKinMin)
	      << " sub = " << srt-em-em2 << endl;
	  exit(1);
          pf=0.0;
        }
        bwf1 *= pf;
    } else if(isw == 4) {
        bwf1 *= XsecTable::pawt2(srt,em,m2);
    }

    return bwf1;
 
}

//***********************************************************************

double SigmaBB1::DetBal2(int msel)
{
//...Purpose to calculate the correction factor of RR->NN cross section.
//...Integration by using Simpson (1/3) rule.
//
    static double emnuc=0.939,empion=0.138;

//...Replace integral to constant width for high energy.
      if(msel ==  1 && srt >= 3.5) {
        return DetBal3()/(pr*pr);
      }

//...Check options.
    int isw=optDetBal;
    int isw1=(isw/10)%10;
    int isw2=isw%10;
    if((isw1 < 1 || isw1 > 3) || (isw2 < 0 || isw2 > 4)) {
        cout << "(DetB2al:)invalid isw " << isw <<endl;
	exit(1);
    }

    //int idbsw=isw;
//...Find min. and max. mass for integration.
    double  em01=pd1->m0();
    double  em02=pd2->m0();
    double  emin1=pd1->m0();
    double  emin2=pd2->m0();

    if(id1 == id_nucls)
        emin1=emnuc+2*empion+eKinMin;
    else if(id1 == id_delt)
        emin1=emnuc+empion+eKinMin;
    else if(id1 == id_delts)
        emin1=emnuc+3*empion+eKinMin;

    if(id2 == id_nucls)
        emin2=emnuc+2*empion+eKinMin;
    else if(id2 == id_delt)
        emin2=emnuc+empion+eKinMin;
    else if(id2 == id_delts)
        emin2=emnuc+3*empion+eKinMin;

    double emax1=srt-emin2-0.001;
    double emax2=srt-emin1-0.001;
    if(emin1 > emax1 || emin2 > emax2) {
        if(msel == 1) return 1.0/(pr*pr);
        if(msel == 2) return 0.0;
    }

//...Min. and max. for first integral.
    double  a=emin1;
    double  b=emax1;

//....Initialization
    int jmax=0;
    int nhlf=0;
    double h=b-a;
    double s1=BW2(a,emin2)+BW2(b,emin2);
    double  s2=0.0;
    double  s4=0.0;
    double  sum=0.5*h*s1;
    double ds=10000;

//...Start the Simpson method
    do {
	jmax=2*jmax+1;
	nhlf++;
	h *= 0.5;
	s2 += s4;
	s4=0.0;
	for(int j=1;j<=jmax;j+=2)
          s4 += BW2(a+j*h,emin2);
	double ss=sum;
	sum=(h/3.0)*(s1+2.0*s2+4.0*s4);
	ds=abs(sum-ss);
	if(nhlf == mxhlf) {
	    cout << "(DetBal2:)simpson does not converge. mxhlf="
		<< mxhlf <<endl;
	    nhlf=-nhlf;
	    break;
	}
    } while(ds > eps*abs(sum));

    double fac=sum/(M_PI*M_PI);
    if(msel == 2) return fac;

    fac=1.0/max(fac,1e-5);
    if(isw2 == 1) {
        fac /= (pr*pr);
    } else if(isw2 == 2 || isw2 == 3) {
        fac /= pr;
        if(isw2 == 3) fac *= (em1/em01)*(em2/em02);
    }

    return fac;

}

//***********************************************************************

double SigmaBB1::BW2(double m1,double emin2)
{
//...Provide B-W function for detbal2.

    double a=emin2;
    double b=srt-m1-0.001;
    if(a > b)  return 0.0;

//...Initialization.
    int  jmax=0;
    int  nhlf=0;
    double h=b-a;
    double s1=BW3(m1,a)+BW3(m1,b);
    double s2=0.0;
    double s4=0.0;
    double sum=0.5*h*s1;
    double ds=100;

//...Start the Simpson method.
    do {
	jmax=2*jmax+1;
	nhlf++;
	h *=0.5;
	s2 += s4;
	s4=0.0;
	for(int j=1;j<=jmax;j+=2) s4 += BW3(m1,a+j*h);
	double ss=sum;
	sum=(h/3.0)*(s1+2.0*s2+4.0*s4);
	if(nhlf == mxhlf) {
	    cout << "simpson method does not converge. mxhlf= " << mxhlf <<endl;
	    break;
	}
	ds=abs(sum-ss);
    } while(ds > eps*abs(sum));

      return sum;

}

//***********************************************************************

double SigmaBB1::BW3(double ema,double emb)
{
    int idbsw=optDetBal;
    int isw=optDetBal%10;
    double gam1,gam2;
//...Calculate momentum dependence decay width.
    //if(idbsw >= 10) {
    if(isw >= 2) {
        gam1=decay->getTotalWidth(pd1,ema);
        gam2=decay->getTotalWidth(pd2,emb);
//...Constant decay width.
    } else {
        gam1=pd1->mWidth();
        gam2=pd2->mWidth();

	if(isDelta(pd1->id())) gam1=decay->getDeltaWidth(3,ema);
	if(isDelta(pd2->id())) gam2=decay->getDeltaWidth(3,emb);

    }
    double em01=pd1->m0();
    double em02=pd2->m0();

    double jambwf3=0.0;
//...Non-rel
    if(idbsw <= 20) 
        jambwf3=0.25*gam1*gam2/(pow2(ema-em01)+0.25*gam1*gam1)
                           /(pow2(emb-em02)+0.25*gam2*gam2);
//...Rel1
    else if(idbsw <= 30) 
        jambwf3=4*ema*em01*gam1/(pow2(ema*ema-em01*em01)+pow2(em01*gam1))
               *emb*em02*gam2/(pow2(emb*emb-em02*em02)+pow2(em02*gam2));
//...Rel2
    else if(idbsw <= 40) 
        jambwf3=4*ema*ema*gam1/(pow2(ema*ema-em01*em01)+pow2(ema*gam1))
               *emb*emb*gam2/(pow2(emb*emb-em02*em02)+pow2(emb*gam2));

    double pa=XsecTable::pawt(srt,ema,emb);
    if(isw == 2 || isw == 3) {
	return jambwf3*pa;
    } else {
	return jambwf3*pa*pa;
    }

}

//***********************************************************************

double SigmaBB1::DetBal3()
{
//...Calculate the correction factor for RR->NN with constant width.
    double emnuc=0.939,empion=0.138;

//...Pole masses.
    double em01=pd1->m0();
    double em02=pd2->m0();

//...Min. masses.
    //double emin1=em01 - pa1->getBWDev();
    double emin1=pd1->m0();
    if(id1 == id_nucls) 
        emin1=emnuc+2*empion+eKinMin;
    else if(id1 ==  id_delt)
        emin1=emnuc+empion+eKinMin;
    else if(id1 == id_delts)
        emin1=emnuc+3*empion+eKinMin;

    //double emin2=em02 - pa2->getBWDev();
    double emin2=pd2->m0();
    if(id2 == id_nucls)
        emin2=emnuc+2*empion+eKinMin;
    else if(id2 == id_delt)
        emin2=emnuc+empion+eKinMin;
    else if(id2 == id_delts)
        emin2=emnuc+3*empion+eKinMin;

//...Decay width.
    double gam1=pd1->mWidth()/2.0;
    double gam2=pd2->mWidth()/2.0;

//...Max. masses.
    double  emax1=srt-emin2-0.001;
    double  emax2=srt-emin1-0.001;

    if(emin1 > emax1 || emin2 > emax2) return 1.0;

//...Min. and max. for first integral.
    double a=emin1;
    double b=emax1;

//...Initialization.
    int  jmax=0;
    int  nhlf=0;
    double  h=b-a;
    double  s1=gam1/((a-em01)*(a-em01)+gam1*gam1)
          *( atan((srt-a-em02)/gam2)-atan((emin2-em02)/gam2) )
        +gam1/(pow2(b-em01)+gam1*gam1)
          *( atan((srt-b-em02)/gam2)-atan((emin2-em02)/gam2) );
   double s2=0.0;
   double s4=0.0;
   double sum=0.5*h*s1;
   double ds=1000;

//...Start of the Simpson method.
    do {
      jmax=2*jmax+1;
      nhlf++;
      h *=0.5;
      s2 += s4;
      s4=0.0;
      for(int j=1;j<=jmax;j+=2) {
       double x=a+j*h;
       s4 += gam1/((x-em01)*(x-em01)+gam1*gam1)
          *( atan((srt-x-em02)/gam2)-atan((emin2-em02)/gam2) );
      }

      double ss=sum;
      sum=(h/3.0)*(s1+2.0*s2+4.0*s4);

      if(nhlf == mxhlf) {
        cout << "simpson method does not converge. mxhlf= " << mxhlf << endl;
      }

      ds=abs(sum-ss);
    } while(ds > eps*abs(sum));

    return max(sum/(M_PI*M_PI),1.0e-5);
}

} // end of namespace jam2
