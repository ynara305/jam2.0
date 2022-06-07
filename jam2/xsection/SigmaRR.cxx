#include <cmath>
#include <fstream>     // To use ifstream
#include <jam2/xsection/SigmaRR.h>
#include <jam2/hadrons/GaussPoints.h>
#include <jam2/hadrons/JamStdlib.h>
#include <jam2/hadrons/BaryonTable.h>

using namespace std;

namespace jam2 {

using Pythia8::pow3;
using Pythia8::pow4;

// Parameters for matrix element.
/*
const double SigmaRR::ParN1440[2][3]= {
    {19.0, 1.9, 2.0}, // T=0
    {2.4,  1.7, 0.6}  // T=1
};

const double SigmaRR::ParR[9][2][3]={
  {{0.0, 0.0, 0.0}, {68.0, 1.104, 1.951}},  // 0) NN->ND(1232)
  {{4.4, 1.9, 2.0}, {2.4, 1.7, 0.6}},       // 1) NN->NN*
  {{0.0, 0.0, 0.0}, {2.4, 1.7, 0.6}},       // 2) NN->ND*
  {{12.0, 1.5, 2.0},{10.0, 1.5, 2.0}},      // 3) NN->DD
  //{{7.0, 1.7, 1.0},{6.0, 1.7, 1.0}},      // 3) NN->DD
  {{0.0, 0.0, 0.0}, {3.3, 1.7, 2.5}},       // 4) NN->DN*
  {{6.0, 1.7, 1.7}, {3.0, 1.7, 1.7}},       // 5) NN->DD*
  {{2.0, 2.0, 2.0}, {0.7, 2.0, 2.0}},       // 6) NN->N*N*
  {{0.0, 0.0, 0.0}, {0.5, 2.0, 2.0}},       // 7) NN->N*D*
  {{1.7, 2.0, 2.0}, {2.0, 2.0, 2.0}},       // 8) NN->D*D*
  {{1.7, 2.0, 2.0}, {2.0, 2.0, 2.0}},       // 8) NN->D*D*
  //{{0.4, 1.2, 0.6}, {0.4, 1.2, 0.6}},       // 8) NN->D*D*
};
*/

/*
const double SigmaRR::ParN1440[2][3]= {
    //{15.0, 1.9, 2.0}, // T=0
    {22.0, 1.9, 2.0}, // T=0
    {2.4,  1.7, 0.6}  // T=1
};

const double SigmaRR::ParR[9][2][3]={
  {{0.0, 0.0, 0.0}, {68.0, 1.104, 1.951}},  // 0) NN->ND(1232)
  {{5.0, 1.9, 2.0}, {2.4, 1.7, 0.6}},       // 1) NN->NN*
  {{0.0, 0.0, 0.0}, {2.4, 1.7, 0.6}},       // 2) NN->ND*
  //{{10.0, 1.5, 2.0},{10.0, 1.5, 2.0}},      // 3) NN->DD
  {{12.0, 1.5, 2.0},{10.0, 1.5, 2.0}},      // 3) NN->DD
  {{0.0, 0.0, 0.0}, {3.3, 1.7, 2.5}},       // 4) NN->DN*

  {{17.0, 1.7, 1.7}, {3.0, 1.7, 1.7}},       // 5) NN->DD*

  {{1.0, 2.0, 2.0}, {0.7, 2.0, 2.0}},       // 6) NN->N*N*
  {{0.0, 0.0, 0.0}, {0.5, 2.0, 2.0}},       // 7) NN->N*D*
  {{1.7, 2.0, 2.0}, {2.0, 2.0, 2.0}},       // 8) NN->D*D*
};
*/

/*
const double SigmaRR::ParN1440[2][3]= {
    //{15.0, 1.9, 2.0}, // T=0
    {22.0, 1.9, 2.0}, // T=0
    {2.4,  1.7, 0.6}  // T=1
};
const double SigmaRR::ParR[9][2][3]={
  {{0.0, 0.0, 0.0}, {68.0, 1.104, 1.951}},  // 0) NN->ND(1232)
  {{5.2, 1.9, 2.0}, {2.4, 1.7, 0.6}},       // 1) NN->NN*
  {{0.0, 0.0, 0.0}, {2.4, 1.7, 0.6}},       // 2) NN->ND*
  {{5.0, 1.5, 2.0}, {10.0, 1.5, 2.0}},      // 3) NN->DD
  {{0.0, 0.0, 0.0}, {3.3, 1.7, 2.5}},       // 4) NN->DN*
  {{18.0, 1.7, 1.7},{3.0, 1.7, 1.7}},       // 5) NN->DD*
  {{0.7, 2.0, 2.0}, {0.7, 2.0, 2.0}},       // 6) NN->N*N*
  {{0.0, 0.0, 0.0}, {0.5, 2.0, 2.0}},       // 7) NN->N*D*
  {{1.7, 2.0, 2.0}, {2.0, 2.0, 2.0}},       // 8) NN->D*D*
};
*/

// RQMD/UrQMD width
const double SigmaRR::ParN1440_1[2][3]= {
    {8.0, 2.0, 2.3}, // T=0
    //{10.0, 1.9, 2.0}, // T=0
    {10.0, 1.1, 0.6}  // T=1
};

// with RQMD/UrQMD width
const double SigmaRR::ParR1[9][2][3]={
  {{0.0, 0.0, 0.0}, {68.0, 1.104, 1.951}},  // 0) NN->ND(1232)
  {{5.2, 1.9, 2.0}, {3.0, 1.1, 0.6}},       // 1) NN->NN*
  {{0.0, 0.0, 0.0}, {3.0, 1.1, 0.6}},       // 2) NN->ND*
  {{10.0, 1.6, 2.0},{9.0, 1.6, 2.0}},       // 3) NN->DD
  {{0.0, 0.0, 0.0}, {2.8, 1.6, 2.0}},       // 4) NN->DN*
  {{25.5,1.6, 2.0}, {3.5, 1.6, 2.0}},       // 5) NN->DD*

  {{0.7, 1.0, 2.0}, {1.0, 1.6, 2.0}},       // 6) NN->N*N*
  {{0.0, 0.0, 0.0}, {0.5, 1.6, 2.0}},       // 7) NN->N*D*
  {{1.7, 0.5, 2.0}, {1.0, 1.6, 2.0}},       // 8) NN->D*D*

  //{{0.7, 1.0, 2.0}, {1.0, 1.6, 2.0}},       // 6) NN->N*N*
  //{{0.0, 0.0, 0.0}, {0.5, 1.6, 2.0}},       // 7) NN->N*D*
  //{{1.7, 0.5, 2.0}, {2.0, 1.6, 2.0}},       // 8) NN->D*D*
};

// GiBUU/SMASH width
const double SigmaRR::ParN1440_2[2][3]= {
    //{8.0, 2.0, 2.3}, // T=0
    //{15.0, 1.90, 2.3}, // T=0
    {20.0, 1.90, 2.3}, // T=0
    //{12.0, 2.0, 2.3}, // T=0
    {3.0, 1.1, 0.6}  // T=1
};

// with GiBUU/SMASH width
const double SigmaRR::ParR2[9][2][3]={
  {{0.0, 0.0, 0.0}, {68.0, 1.104, 1.951}},  // 0) NN->ND(1232)
  {{4.3, 1.9, 2.0}, {3.0, 1.1, 0.6}},       // 1) NN->NN*
  {{0.0, 0.0, 0.0}, {2.0, 1.1, 0.6}},       // 2) NN->ND*
  {{11.0, 1.6, 2.0},{9.0, 1.6, 2.0}},       // 3) NN->DD
  {{0.0, 0.0, 0.0}, {2.8, 1.6, 2.0}},       // 4) NN->DN*
  {{24.0,1.6, 2.0}, {3.5, 1.6, 2.0}},       // 5) NN->DD*

  {{0.7, 1.0, 2.0}, {0.5, 1.6, 2.0}},       // 6) NN->N*N*
  {{0.0, 0.0, 0.0}, {0.5, 1.6, 2.0}},       // 7) NN->N*D*
  {{1.7, 0.5, 2.0}, {0.5, 1.6, 2.0}},       // 8) NN->D*D*

  //{{0.7, 1.0, 2.0}, {0.0, 1.6, 2.0}},       // 6) NN->N*N*
  //{{0.0, 0.0, 0.0}, {0.0, 1.6, 2.0}},       // 7) NN->N*D*
  //{{1.7, 0.5, 2.0}, {0.0, 1.6, 2.0}},       // 8) NN->D*D*

};

// model=1
const double SigmaRR::A1[9][2]={
      {0.0,   0.0},  // 0) N N -> N D
      {20.0, 10.0},  // 1) N N -> N N*
      {0.0,  10.0},  // 2) N N -> N D*
      {40.0, 35.0},  // 3) N N -> D(1232) D(1232) with symmetry factor
      {0.0,  25.0},  // 4) N N -> D(1232) N*
      {10.0, 25.0},  // 5) N N -> D(1232) D*
      {0.0,   0.0},  // 6) N N -> N* N*
      {0.0,   0.0},  // 7) N N -> N* D*
      {0.0,   0.0},  // 8) N N -> D* D*
};

// model=2
const double SigmaRR::A2[9][2]={
      {0.0,  0.0},  // 0) N N -> N N*
      {20.0, 10.0},  // 1) N N -> N N*
      {0.0,  13.0},  // 2) N N -> N D*
      {40.0, 35.0},  // 3) N N -> D(1232) D(1232) with symmetry factor
      {0.0,   1.8},  // 4) N N -> D(1232) N*
      {10.0,  3.0},  // 5) N N -> D(1232) D*

      //{2.0,   0.6},  // 6) N N -> N* N*
      //{0.0,   0.5},  // 7) N N -> N* D*
      //{3.5,   0.5},  // 8) N N -> D* D*

      {0.5,   0.5},  // 6) N N -> N* N*
      {0.0,   0.5},  // 7) N N -> N* D*
      {0.5,   0.5},  // 8) N N -> D* D*
};

// model=2 with GiBUU resonance width
const double SigmaRR::A3[9][2]={
      {0.0,  0.0},  // 0) N N -> N N*
      {20.0, 10.0},  // 1) N N -> N N*
      {0.0,  13.0},  // 2) N N -> N D*
      {40.0, 35.0},  // 3) N N -> D(1232) D(1232) with symmetry factor
      {0.0,   1.8},  // 4) N N -> D(1232) N*
      {10.0,  3.0},  // 5) N N -> D(1232) D*

      {0.5,   0.7},  // 6) N N -> N* N*
      {0.0,   0.7},  // 7) N N -> N* D*
      {0.5,   0.7},  // 8) N N -> D* D*

      // 12GeV
      //{0.5,   0.5},  // 6) N N -> N* N*
      //{0.0,   0.5},  // 7) N N -> N* D*
      //{0.5,   0.5},  // 8) N N -> D* D*
};

const double SigmaRR::Mnucl=0.9383,SigmaRR::Mpion=0.138;
const double SigmaRR::Mdelta=1.232,SigmaRR::Wdelta=0.117;
const double SigmaRR::eKinmi=0.0001;
//const double SigmaRR::minMass=Mnucl + Mpion + eKinmi;
const double SigmaRR::minMass=Mnucl + Mpion;

const double SigmaRR::sMin = 2*Mnucl + Mpion;
const double SigmaRR::sMax=7.0;
const int  SigmaRR::nS = 500;
const double  SigmaRR::dS = (sMax - sMin) / (nS - 1);

SHadron::SHadron(int j,int s,double m,double w,double mm,
  int ndec, double* br, double* m1, double* m2, int* s1, int* s2, int* l)
{
  idp=j; spin=s; mass=m; gam=w;
  minMass=mm;
  maxMass=6.0;
  nMass = 500;
  dMass = (maxMass - minMass) / (nMass - 1);
  width = new double [nMass];
  BWint.clear();

  nDecayMode=ndec;
  branch = br;
  mass1 = m1;
  mass2 = m2;
  spin1 = s1;
  spin2 = s2;
  angMom = l;

  facF = 0.2;  // Frankfurt
  lam2=HBARC*HBARC;

}

SHadron::~SHadron()
{ 
  delete [] width;
  for(int i=0;i<(int)BWint.size();i++) delete [] BWint[i];
}

void SHadron::computeNorm(int NG, double* xg, double* wg,int optBW,int optWidth)
{
  double ymin=2*atan(2*(minMass-mass)/gam);
  double ymax=M_PI;
  double dy = ymax - ymin;
  Norm=0.0;
  for (int i=0;i<NG;i++) {
    double y= ymin + xg[i]*dy;
    double emd=gam/2*tan(0.5*y)+mass;
    double dc= gam/( 0.25*gam*gam + pow2(emd - mass) );
    double gam1=getTotalWidth(emd,optWidth);
    Norm += SigmaRR::BW(emd,gam1,mass,optBW)*dy*wg[i]/dc;
  }
}

  inline double SHadron::getBW(double srt,int ip) {
    int i=(int)floor((srt-SigmaRR::sMin)/SigmaRR::dS);
    if(i>=0 && i<SigmaRR::nS-1) {
      double si=SigmaRR::sMin+i*SigmaRR::dS;
      double s = (srt-si)/SigmaRR::dS;
      return BWint[ip][i]*(1.0-s) + s*BWint[ip][i+1];
    } else if(i >= SigmaRR::nS ) {
      return BWint[ip][SigmaRR::nS-1] + 0.5*(srt-SigmaRR::sMax);
    } else {
      return 0.0;
    }
  }

double SHadron::deltaWidth(double emd,int optw)
{
  if(emd <= SigmaRR::Mnucl+SigmaRR::Mpion+SigmaRR::eKinmi) return 0.0;
  double pr=PCM(emd,SigmaRR::Mnucl,SigmaRR::Mpion);
  const int ldec = 1;
  // ppres= 0.227895;
  double prres=PCM(SigmaRR::Mdelta,SigmaRR::Mnucl,SigmaRR::Mpion);

  // RQMD/UrQMD
  if(optw==1) {
      return SigmaRR::Wdelta*pow(pr/prres,2*ldec+1)*(SigmaRR::Mdelta/emd)
      *(1.0 + facF)/(1.0+facF*pow(pr/prres,2*ldec));

  // GiBUU,SMASH
  } else if(optw==2)  {
      return SigmaRR::Wdelta*pow(pr/prres,2*ldec+1)*(SigmaRR::Mdelta/emd)
	  *pow((prres*prres+lam2)/(pr*pr+lam2),ldec);

  } else {

  const double pscal1= 0.238, pscal2= 0.318, p0ref=0.227, widdlt=0.12;
// Randrup NP A314 (1979) 429.Rittenberg, REV.MOD.PHYS. 43 (1971) S1.
    return widdlt*(pr*pr*pr/(1.0+pow2(pr/pscal1)+pow4(pr/pscal2)))
         /(pow3(p0ref)/(1.0+pow2(p0ref/pscal1)+pow4(p0ref/pscal2)));
  }

}

void SHadron::makeDeltaWidthTable(int optw)
{
  // make table loop over mass.
  for(int i=0;i<nMass;i++) {
    width[i]=deltaWidth(minMass + i*dMass,optw);
  }
}

double SHadron::formBlattWeisskopf(double m0, double m, double pf0, double pf, int l) 
{
  double R=1.0/HBARC;
  double bl  = BlattWeisskopf(pf*R,l);
  double bl0 = BlattWeisskopf(pf0*R,l);
  return pf/pf0 * m0/m * (bl*bl)/(bl0*bl0);
}
  
// Blatt Weisskopf barrier penetration factor
// D.M.Manley and E.M.Saleski, Phys.Rev. D45 (1992) 4002
// D.M.Manley, R.A.Arndt, Y.Goradia, and V.L.Teplitz, Phys.Rev.D30,904(1984)
// J.M.Blatt and V.F.Weisskopf, "Theoretical Nuclear Physics" 1952
// Feshbach "Theoretical Nuclear Physics Nuclear Reactions" 1992, p245
double SHadron::BlattWeisskopf(double x,int l)
{
  switch (l) {
    case 0: return 1;
          break;
    case 1: return x*x/(1+x*x);
          break;
    case 2: return x*x*x*x/(9 + 3*x*x + x*x*x*x);
          break;
    case 3: return pow(x,6)/(225 + 45*x*x + 6*pow4(x) + pow(x,6));
          break;
    case 4:
      return pow(x,8)/(11025 + 1575*x*x + 135*pow4(x) + 10*pow(x,6) + pow(x,8));
      break;
    default: cout << " angular mom > 4? " << l <<endl;
       exit(1);
  }
}

double SHadron::getTotalWidth(double emd,int opt)
{
  double totwid=0.0;
  for (int j=0;j<nDecayMode;j++) {  // loop over decay mode.

    double m1=mass1[j];
    double m2=mass2[j];
    int l=angMom[j];

    if(branch[j] < 1e-9) continue;
    if( emd < m1 + m2 + SigmaRR::eKinmi ||  mass < m1 + m2+ SigmaRR::eKinmi) continue;
    //if( emd <= m1 + m2  ||  mass <= m1 + m2 ) continue;
    double pr=PCM(emd,m1,m2);
    double pr0=PCM(mass,m1,m2);
    double form;
    if( opt==1 )
      form =  mass/emd*pow(pr/pr0,2*l+1)*(1.0+facF)/(1+facF*pow(pr/pr0,2*l));
    else if( opt==2)
      form =  mass/emd*pow(pr/pr0,2*l+1)*((pr0*pr0+lam2)/(pr*pr+lam2));
    else
      form = formBlattWeisskopf(mass, emd, pr0, pr, l);

    totwid += branch[j]*form;
  }
  return totwid*gam;
}

void SHadron::makeWidthTable(const int optW)
{
  // make table loop over mass.
  for(int i=0;i<nMass;i++) {
    double emd=minMass + i*dMass;
    width[i]=getTotalWidth(emd,optW);
  }
}

SigmaRR::SigmaRR(string bwfile,int optrr,int opt1, int opt2, int optbw, JamParticleData* jamtable)
{
  jamTable = jamtable;
  int optG=1;
  if(optG==1) {
    NG=38;
    xg = new double [NG];
    wg = new double [NG];
    for(int i=0;i<NG;i++) {
      xg[i]= xg38f[i];
      wg[i]= wg38f[i];
    }
  } else if(optG==2) {
    NG=48;
    xg = new double [NG];
    wg = new double [NG];
    for(int i=0;i<NG;i++) {
      xg[i]= xg48f[i];
      wg[i]= wg48f[i];
    }
  } else {
    NG=100;
    xg = new double [NG];
    wg = new double [NG];
    for(int i=0;i<NG;i++) {
      xg[i]= xg100f[i];
      wg[i]= wg100f[i];
    }
  }

  optDeltaWidth = opt1;
  optWidth=opt2;
  optBW=optbw;
  modelRR=optrr;
  //modelRR=1;   //=1: NN->ND, NN*, DD, DN*, DD*
  //modelRR=2;   //=2: NN->ND, NN*, DD, DN*, DD*, N*N*,, N*D*, D*D*
 
  optRR=1;    // matrix element is constant for NN->N*N*, N*D*, D*D*
  //optRR=2; // same functional form  of matrix element for all processes.

  if(modelRR==1 && optRR==2) optRR=1;

  makeNstarTable();
  makeDstarTable();

  if(!readBWTable(bwfile)) makeBWTable(bwfile);
  if(modelRR==1) {
    for(int i=0;i<9;i++) {
      A[i][0]=A1[i][0];
      A[i][1]=A1[i][1];
    }
  } else {

    if(optWidth==1) {
      for(int i=0;i<3;i++) {
        ParN1440[0][i]=ParN1440_1[0][i];
        ParN1440[1][i]=ParN1440_1[1][i];
      }
      for(int i=0;i<9;i++) {
        A[i][0]=A2[i][0];
        A[i][1]=A2[i][1];
	for(int j=0;j<3;j++) {
          ParR[i][0][j]=ParR1[i][0][j];
          ParR[i][1][j]=ParR1[i][1][j];
	}
      }
	  for(int i=0;i<9;i++) {
           for(int j=0;j<3;j++) ParR[i][1][j]=ParR2[i][1][j];
	  }

    } else {

      for(int i=0;i<3;i++) {
        ParN1440[0][i]=ParN1440_2[0][i];
        ParN1440[1][i]=ParN1440_2[1][i];
      }
      for(int i=0;i<9;i++) {
        for(int j=0;j<2;j++) {
          A[i][j]=A3[i][j];
	  for(int k=0;k<3;k++) ParR[i][j][k]=ParR2[i][j][k];
        }
      }
	  for(int k=0;k<2;k++)
	  for(int i=0;i<9;i++) {
	      for(int j=0;j<3;j++)  {
	      ParR[i][k][j]=ParR2[i][k][j];
	  }
	  }
	  for(int i=0;i<9;i++)
	      for(int j=0;j<3;j++)  {
	      //ParR[i][0][j]=ParR2[i][0][j];
	      ParR[i][1][j]=ParR2[i][1][j];
	  }

	    //ParR[i][0][0]=ParR2[i][0][0];
	    //ParR[i][0][1]=ParR2[i][0][1];
	    //ParR[i][0][2]=ParR2[i][0][2];
	    //ParR[i][1][1]=ParR2[i][1][1];
	    //ParR[i][1][2]=ParR2[i][1][2];
	    //ParR[i][1][0]=ParR2[i][1][0];

    }
  }

  /*
      for(int i=0;i<3;i++) {
        cout << i << " N1440= "<< ParN1440[0][i] << " PartN1440= "<< ParN1440_2[0][i]  <<endl;
        cout << i << " N1440= "<< ParN1440[1][i] << " PartN1440= "<< ParN1440_2[1][i]  <<endl;
      }
      for(int i=0;i<9;i++) {
        cout << "A0= "<< A[i][0] << " A3= "<< A3[i][0] <<endl;
        cout << "A1= "<< A[i][1] << " A3= "<< A3[i][1]  <<endl;
        cout << i << " ParR0= "<< ParR[i][0][0] << " " << ParR[i][0][1] << " " << ParR[i][0][2] <<endl;
        cout << i << " ParR1= "<< ParR[i][1][0] << " " << ParR[i][1][1] << " " << ParR[i][1][2]  <<endl;
        cout << i << " ParR0= "<< ParR2[i][0][0] << " " << ParR2[i][0][1] << " " << ParR2[i][0][2] <<endl;
        cout << i << " ParR1= "<< ParR2[i][1][0] << " " << ParR2[i][1][1] << " " << ParR2[i][1][2]  <<endl;
      }
      //cin.get();
      */

}

SigmaRR::~SigmaRR()
{
  delete [] xg;
  delete [] wg;
}

void SigmaRR::makeNstarTable()
{
  using namespace nucleon_resonance;

  double m1[nDecayMode],m2[nDecayMode],br[nDecayMode];
  int    sp1[nDecayMode],sp2[nDecayMode],l[nDecayMode];
  for(int j=0;j<nDecayMode;j++) {
    m1[j]=mass1[j];
    m2[j]=mass2[j];
    sp1[j]=spin1[j];
    sp2[j]=spin2[j];
    br[j]=branch[0][j];
    l[j]=  abs((spin[0] + 1 - spin1[j] - spin2[j])/2);
  }
  nucleon = new SHadron(0,2,0.9383,0.0,0.9383,0,br,m1,m2,sp1,sp2,l);

  for(int i=0;i<num;i++) {

    for(int j=0;j<nDecayMode;j++) {
      l[j]=  abs((spin[i] + 1 - spin1[j] - spin2[j])/2);
      br[j]=branch[i][j];
    }

    SHadron* p = new SHadron(pdg1[i],spin[i],mass[i],width[i],minMass,
                nDecayMode,br,m1,m2,sp1,sp2,l);

    p->makeWidthTable(optWidth);
    p->computeNorm(NG, xg, wg,optBW,optWidth);

    nStar.push_back(p);

    //cout << "# N*= "<< i+1 << " m= "<< p->mass0()
//	<< " norm= "<< p->norm()
//	<< " optW= "<< optWidth
//	<<endl;
  }

}

void SigmaRR::makeDstarTable()
{
  using namespace delta_resonance;
  double m1[nDecayMode],m2[nDecayMode],br[nDecayMode];
  int    sp1[nDecayMode],sp2[nDecayMode],l[nDecayMode];
  for(int j=0;j<nDecayMode;j++) {
    m1[j]=mass1[j];
    m2[j]=mass2[j];
    sp1[j]=spin1[j];
    sp2[j]=spin2[j];
    br[j]=branch[0][j];
    l[j]=  abs((spin[0] + 1 - spin1[j] - spin2[j])/2);
  }
  double brd[]={0.0  ,1.0 ,0.0 ,0.0 ,0.0};
  d1232 = new SHadron(1114,4,1.232,0.117,minMass,2,brd,m1,m2,sp1,sp2,l);
  d1232->makeDeltaWidthTable(optDeltaWidth);
  d1232->computeNorm(NG, xg, wg,optBW,optDeltaWidth);
  //cout << "# norm1232= "<< d1232->norm() <<endl;

  for(int i=0;i<num;i++) {

    for(int j=0;j<nDecayMode;j++) {
      br[j]=branch[i][j];
      l[j]=  abs((spin[i] + 1 - spin1[j] - spin2[j])/2);
    }
    SHadron* p = new SHadron(pdg1[i],spin[i],mass[i],width[i],minMass,
                nDecayMode,br,m1,m2,sp1,sp2,l);
    p->makeWidthTable(optWidth);
    p->computeNorm(NG, xg, wg,optBW,optWidth);
    dStar.push_back(p);
    //cout << "# D*= "<< i+1 << " m= "<< p->mass0()
//	<< " norm= "<< p->norm()
//	<< " optW= "<< optWidth
//	<<endl;
  }

}

//####### single resonance production #################################
double SigmaRR::BWint(double srt, double em2, SHadron* pa)
{
  if (srt <= 2*Mnucl + Mpion ) return 0.0;
  double emin=pa->mMin();
  double emax=srt-em2;
  double em0 = pa->mass0();
  double gam0 = pa->width0();
  double ymin=2*atan(2*(emin - em0)/gam0);
  double ymax=2*atan(2*(emax - em0)/gam0);
  double dy=ymax-ymin;
  double bwint=0.0;
  for(int i=0;i<NG;i++) {
    double y=ymin+xg[i]*dy;
    double emd=gam0/2*tan(0.5*y)+em0;
    double dc= gam0/( 0.25*gam0*gam0 + pow2(emd - em0) );
    double gam=pa->getWidth(emd);
    if(srt < emd+em2) {
      cout << "BWint srt= " << srt << " emd= "<< emd << " em2= "<< endl;
      cout << " ymin= "<< ymin << " ymax= "<< ymax<< endl;
      exit(1);
    }
    double pf=PCM(srt,emd,em2);
    bwint += BW(emd,gam,em0,optBW)*pf*dy*wg[i]/dc;
  }
  return bwint;

}

double SigmaRR::BWint2(double srt,SHadron* pa1,SHadron* pa2)
{
  if(srt < 2*Mnucl + 2*Mpion) return 0.0;
  double emin1=Mnucl+Mpion;
  double emin2=Mnucl+Mpion;
  double emax1=srt-emin2;
  double emr1=pa1->mass0();
  double gamr1=pa1->width0();
  double ymax=2*atan(2*(emax1 - emr1)/gamr1);
  double ymin=2*atan(2*(emin1 - emr1)/gamr1);
  double dy=ymax-ymin;
  if (dy <= 0.0) {
    cout << "SigmaRR:BWint2 dy=" << dy<<endl; 
    cout << "srt= "<< srt << " emin1= "<<  emin1 << " emax1= "<< emax1<<endl;
    exit(1);
  }

  double bwint=0.0;
  for(int i=0;i<NG;i++) {
    double y=ymin+xg[i]*dy;
    double emd=gamr1/2*tan(0.5*y)+emr1;
    double dc= gamr1/( 0.25*gamr1*gamr1 + pow2(emd - emr1) );
    double gam=pa1->getWidth(emd);
    bwint += BW(emd,gam,emr1,optBW) * BWint(srt,emd,pa2) * dy*wg[i]/dc;
  }
  return bwint;

}

bool SigmaRR::readBWTable(string fname)
{
  ifstream inFile(fname.c_str(),ios::in);
  if(!inFile.is_open()) {
    cout << "SigmaRR::readBWTable file does not exist " << fname
	 << " now making table" << endl;
    return false;
  }

  double *bwint = new double [nS];

    //while (count < ARRAY_SIZE && inputFile >> numbers[count]) count++;

  // N N -> N Delta
  for(int i=0;i<nS;i++) inFile >> bwint[i];
  nucleon->setBWint(nS,bwint);

  // N N -> N N*
  for(int ip=0;ip<(int)nStar.size();ip++) {
    for(int i=0;i<nS;i++) inFile >> bwint[i];
    nucleon->setBWint(nS,bwint);
  }

  // N N -> N D*
  for(int ip=0;ip<(int)dStar.size();ip++) {
    for(int i=0;i<nS;i++) inFile >> bwint[i];
    nucleon->setBWint(nS,bwint);
  }

  // N N -> D(1232) D(1232)
  for(int i=0;i<nS;i++) inFile >> bwint[i];
  d1232->setBWint(nS,bwint);

  // N N -> D(1232) N*
  for(int ip=0;ip<(int)nStar.size();ip++) {
    for(int i=0;i<nS;i++) inFile >> bwint[i];
    d1232->setBWint(nS,bwint);
  }

  // N N -> D(1232) D*
  for(int ip=0;ip<(int)dStar.size();ip++) {
    for(int i=0;i<nS;i++) inFile >> bwint[i];
    d1232->setBWint(nS,bwint);
  }

  if(modelRR==1) {
    delete [] bwint;
    return true;
  }

  // N N -> N* N*
  for(int i1=0;i1<(int)nStar.size();i1++) {
    for(int i2=i1;i2<(int)nStar.size();i2++) {
      for(int i=0;i<nS;i++) inFile >> bwint[i];
    nStar[i1]->setBWint(nS,bwint);
    }
  }

  // N N -> N* D*
  for(int i1=0;i1<(int)nStar.size();i1++) {
    for(int i2=0;i2<(int)dStar.size();i2++) {
      for(int i=0;i<nS;i++) inFile >> bwint[i];
    nStar[i1]->setBWint(nS,bwint);
    }
  }

  // N N -> D* D*
  for(int i1=0;i1<(int)dStar.size();i1++) {
    for(int i2=i1;i2<(int)dStar.size();i2++) {
      for(int i=0;i<nS;i++) inFile >> bwint[i];
    dStar[i1]->setBWint(nS,bwint);
    }
  }

  delete [] bwint;
  inFile.close();
  return true;
}

void SigmaRR::makeBWTable(string outfile)
{
  ofstream ofs(outfile.c_str());

  double *bwint = new double [nS];

  // N N -> N Delta
  for(int i=0;i<nS;i++) {
    double srt = sMin + dS * i;
    bwint[i] = BWint(srt,Mnucl,d1232)/d1232->norm();
    ofs << bwint[i]  << endl;
  }
  nucleon->setBWint(nS,bwint);

  // N N -> N N*
  for(int ip=0;ip<(int)nStar.size();ip++) {
    for(int i=0;i<nS;i++) {
      double srt = sMin + dS * i;
      bwint[i] = BWint(srt,Mnucl,nStar[ip])/nStar[ip]->norm();
      ofs << bwint[i]  << endl;
    }
    nucleon->setBWint(nS,bwint);
  }

  // N N -> N D*
  for(int ip=0;ip<(int)dStar.size();ip++) {
    for(int i=0;i<nS;i++) {
      double srt = sMin + dS * i;
      bwint[i] = BWint(srt,Mnucl,dStar[ip])/dStar[ip]->norm();
      ofs << bwint[i]  << endl;
    }
    nucleon->setBWint(nS,bwint);
  }

  // N N -> D(1232) D(1232)
  double normd=d1232->norm();
  for(int i=0;i<nS;i++) {
      double srt = sMin + dS * i;
      bwint[i] = BWint2(srt,d1232,d1232)/normd/normd;
      ofs << bwint[i]  << endl;
  }
  d1232->setBWint(nS,bwint);

  // N N -> D(1232) N*
  for(int ip=0;ip<(int)nStar.size();ip++) {
    for(int i=0;i<nS;i++) {
      double srt = sMin + dS * i;
      bwint[i] = BWint2(srt,d1232,nStar[ip])/nStar[ip]->norm()/normd;
      ofs << bwint[i]  << endl;
    }
    d1232->setBWint(nS,bwint);
  }

  // N N -> D(1232) D*
  for(int ip=0;ip<(int)dStar.size();ip++) {
    for(int i=0;i<nS;i++) {
      double srt = sMin + dS * i;
      bwint[i] = BWint2(srt,d1232,dStar[ip])/dStar[ip]->norm()/normd;
      ofs << bwint[i]  << endl;
    }
    d1232->setBWint(nS,bwint);
  }

  if(modelRR==1) {
    delete [] bwint;
    ofs.close();
    return;
  }

  // N N -> D(1232) D*
  // N N -> N* N*
  for(int i1=0;i1<(int)nStar.size();i1++) {
    for(int i2=i1;i2<(int)nStar.size();i2++) {
      for(int i=0;i<nS;i++) {
        double srt = sMin + dS * i;
        bwint[i] = BWint2(srt,nStar[i1],nStar[i2])
                   /nStar[i1]->norm()/nStar[i2]->norm();
        ofs << bwint[i]  << endl;
      }
    nStar[i1]->setBWint(nS,bwint);
    }
  }

  // N N -> N* D*
  for(int i1=0;i1<(int)nStar.size();i1++) {
    for(int i2=0;i2<(int)dStar.size();i2++) {
      for(int i=0;i<nS;i++) {
        double srt = sMin + dS * i;
        bwint[i] = BWint2(srt,nStar[i1],dStar[i2])
                   /nStar[i1]->norm()/dStar[i2]->norm();
        ofs << bwint[i]  << endl;
      }
    nStar[i1]->setBWint(nS,bwint);
    }
  }

  // N N -> D* D*
  for(int i1=0;i1<(int)dStar.size();i1++) {
    for(int i2=i1;i2<(int)dStar.size();i2++) {
      for(int i=0;i<nS;i++) {
        double srt = sMin + dS * i;
        bwint[i] = BWint2(srt,dStar[i1],dStar[i2])
                   /nStar[i1]->norm()/dStar[i2]->norm();
        ofs << bwint[i]  << endl;
      }
    dStar[i1]->setBWint(nS,bwint);
    }
  }

  delete [] bwint;
  ofs.close();

}

// N N -> N D(1232)
double SigmaRR::sigDelta2(double srt, double m1, double m2)
{
  if (srt <= 2*Mnucl+Mpion+0.0005) return 0.0;
  //double rmat = ParDelta[0]/pow(srt-ParDelta[1],ParDelta[2]);
  double rmat = ParR[0][1][0]/pow(srt-ParR[0][1][1],ParR[0][1][2]);
  double pr0=PCM(srt,m1,m2);
  int spin=2*4;
  return spin*rmat*BWint(srt,Mnucl,d1232)/(pr0*srt*srt)/d1232->norm();
}

// N N -> N D(1232)
double SigmaRR::sigDelta(double srt, double m1, double m2)
{
  if (srt <= 2*Mnucl+Mpion+0.0005) return 0.0;
  //double rmat = 71.0/pow(srt-1.104,1.951);
  //double rmat = ParDelta[0]/pow(srt-ParDelta[1],ParDelta[2]);
  double rmat = ParR[0][1][0]/pow(srt-ParR[0][1][1],ParR[0][1][2]);
  double pr0=PCM(srt,m1,m2);
  int spin=2*4;
  return spin*rmat*nucleon->getBW(srt,0)/(pr0*srt*srt);
}


// N N -> N N*
double SigmaRR::sigNNstar(double srt,double m1, double m2, int iso, int ip)
{
  if (srt <= 2*Mnucl+Mpion+0.0005) return 0.0;
  if(srt < m1 + m2 + eKinmi) return 0.0;
  double flux = PCM(srt,m1,m2)*srt*srt;
  int i1= ip==0 ? 0: ip;
  int i2= ip==0 ? nStar.size() : ip+1;
  double sig=0.0;
  for(int i=i1;i<i2;i++) {
    int spin=2*nStar[i]->getSpin();
    //double em0=nStar[i]->mass0();
    //double rmat = A[0][iso]/(2*(Mnucl*Mnucl + em0*em0 ));
    double rmat = ParR[1][iso][0]/pow(srt-ParR[1][iso][1],ParR[1][iso][2]);
    if (i==0) 
      rmat = ParN1440[iso][0]/pow(srt-ParN1440[iso][1], ParN1440[iso][2]);
    sig += spin*rmat*nucleon->getBW(srt,i+1)/flux;

    /*
    //double sig0 = spin*rmat*BWint(srt,Mnucl,nStar[i])/flux/nStar[i]->norm();
    cout << "srt= "<< srt
	<< " flux= "<< flux
	<< " rmat= "<< rmat
	<< " bw= "<< nucleon->getBW(srt,i+1)
	<< " sig= "<< rmat*nucleon->getBW(srt,i+1)/flux
	<< " i= "<< i
	<<endl;
	*/
  }
  return sig;
}


// N N -> N D*
double SigmaRR::sigNDstar(double srt,double m1, double m2, int ip)
{
  if (srt <= 2*Mnucl+Mpion+0.0005) return 0.0;
  //const double A=17.0;  // JAM2 optW=1
  //double flux = PCM(srt,Mnucl,Mnucl)*srt*srt;
  if(srt < m1 + m2 + eKinmi) return 0.0;
  double flux = PCM(srt,m1,m2)*srt*srt;
  double sig=0.0;
  int i1=ip;
  int i2=ip+1;
  if(ip==0) {
    i1=0;
    i2=dStar.size();
  }
  int shift=1+nStar.size();
  for(int i=i1; i<i2; i++) {
    int spin=2*dStar[i]->getSpin();
    //double em0=dStar[i]->mass0();
    //double rmat = A[1][1]/(2*(Mnucl*Mnucl + em0*em0 ));
    //double rmat = ParNDs[0]/pow(srt-ParNDs[1],ParNDs[2]);
    double rmat = ParR[2][1][0]/pow(srt-ParR[2][1][1],ParR[2][1][2]);
    sig += spin*rmat*nucleon->getBW(srt,i+shift)/flux;
    //sig += spin*rmat*BWint(srt,Mnucl,Dstar[i])/(pr0*srt*srt)/dStar[i]->norm();
  }
  return sig;
}


// N N -> D(1232) D(1232)
double SigmaRR::sigDeltaDelta(double srt,double m1, double m2,int iso)
{
  if (srt <= 2*Mnucl+2*Mpion+0.0005) return 0.0;
  //const double A={120.0, 50.0} // JAM2 optw=1
  double flux = PCM(srt,m1,m2)*srt*srt;
  //double rmat = A[3][iso]/( 4*Mdelta*Mdelta );
  double rmat = ParR[3][iso][0]/pow(srt-ParR[3][iso][1],ParR[3][iso][2]);
  const int spin=4*4;
  return spin*rmat*d1232->getBW(srt,0)/flux;
}


// N N -> D(1232) N*
double SigmaRR::sigDNstar(double srt,double m1, double m2, int ip)
{
  if (srt <= 2*Mnucl+2*Mpion+0.0005) return 0.0;
  //const double A=7.0;  // JAM2 optW=1
  double flux = PCM(srt,m1,m2)*srt*srt;
  int i1=ip;
  int i2=ip+1;
  if(ip==0) {
    i1=0;
    i2=nStar.size();
  }
  double sig=0.0;
  for(int i=i1;i<i2;i++) {
    int spin=4*nStar[i]->getSpin();
    //double em0=nStar[i]->mass0();
    //double rmat = A[3][1]/(2*(Mdelta*Mdelta + em0*em0 ));
    //double rmat = ParDNs[0]/pow(srt-ParDNs[1],ParDNs[2]);
    double rmat = ParR[4][1][0]/pow(srt-ParR[4][1][1],ParR[4][1][2]);
    sig += spin*rmat*d1232->getBW(srt,i+1)/flux;
  }
  return sig;
}


// N N -> D(1232) D*
double SigmaRR::sigDDstar(double srt,double m1, double m2, int iso, int ip)
{
  if (srt <= 2*Mnucl+2*Mpion+0.0005) return 0.0;
  //const double  A={7.0, 12.0}  // JAM2 optW=1
  double flux = PCM(srt,m1,m2)*srt*srt;
  int i1=ip;
  int i2=ip+1;
  if(ip==0) {
    i1=0;
    i2= dStar.size();
  }
  int shift=1+nStar.size();
  double sig=0.0;
  for(int i=i1; i<i2; i++) {
    int spin=4*dStar[i]->getSpin();
    //double em0=dStar[i]->mass0();
    //double rmat = A[4][iso]/(2*(Mdelta*Mdelta + em0*em0 ));
    //double rmat = ParDDs[iso][0]/pow(srt-ParDDs[iso][1], ParDDs[iso][2]);
    double rmat = ParR[5][iso][0]/pow(srt-ParR[5][iso][1],ParR[5][iso][2]);
    sig += spin*rmat*d1232->getBW(srt,i+shift)/flux;
  }
  return sig;
}

// N N -> N* N*
double SigmaRR::sigNstarNstar(double srt,double m1, double m2,
	int iso,int j1,int j2)
{
  if (srt <= 2*Mnucl+2*Mpion+0.0005) return 0.0;
  int i1 = std::min(j1,j2);
  int i2 = std::max(j1,j2);
  double flux = PCM(srt,m1,m2)*srt*srt;
  int spin=nStar[i1]->getSpin()*nStar[i2]->getSpin();
  double rmat;
  if(optRR==1) {
    double m01=nStar[i1]->mass0();
    double m02=nStar[i2]->mass0();
    rmat = A[6][iso]/(2*(m01*m01 + m02*m02 ));
  } else {
    rmat = ParR[6][iso][0]/pow(srt-ParR[6][iso][1],ParR[6][iso][2]);
  }
  //if(iso==1 && i1 == i2) rmat /=2;
  return spin*rmat*nStar[i1]->getBW(srt,i2-i1)/flux;
}

// N N -> N* N*
double SigmaRR::sigNstarNstar(double srt,double m1, double m2, int iso)
{
  double sig=0.0;
  for(int i=0;i<(int)nStar.size();i++) {
    for(int j=i;j<(int)nStar.size();j++) {
      sig += sigNstarNstar(srt,m1,m2,iso,i,j);
    }
  }
  return sig;
}

// N N -> N* D*
double SigmaRR::sigNstarDstar(double srt,double m1, double m2,
	int i1,int i2)
{
  if (srt <= 2*Mnucl+2*Mpion+0.0005) return 0.0;
  double flux = PCM(srt,m1,m2)*srt*srt;
  int spin=nStar[i1]->getSpin()*dStar[i2]->getSpin();
  double rmat;
  if(optRR==1) {
    double m01=nStar[i1]->mass0();
    double m02=dStar[i2]->mass0();
    rmat = A[7][1]/(2*(m01*m01 + m02*m02 ));
  } else {
    rmat = ParR[7][1][0]/pow(srt-ParR[7][1][1],ParR[7][1][2]);
  }
  int shift=nStar.size()-i1;
  return spin*rmat*nStar[i1]->getBW(srt,i2+shift)/flux;
}

// N N -> N* D*
double SigmaRR::sigNstarDstar(double srt,double m1, double m2)
{
  double sig=0.0;
  for(int i=0;i<(int)nStar.size();i++) {
    for(int j=0;j<(int)dStar.size();j++) {
      sig += sigNstarDstar(srt,m1,m2,i,j);
    }
  }
  return sig;
}

// N N -> D* D*
double SigmaRR::sigDstarDstar(double srt,double m1, double m2,
	int iso,int j1,int j2)
{
  if (srt <= 2*Mnucl+2*Mpion+0.0005) return 0.0;
  int i1 = min(j1,j2);
  int i2 = max(j1,j2);
  double flux = PCM(srt,m1,m2)*srt*srt;
  int spin=dStar[i1]->getSpin()*dStar[i2]->getSpin();
  double rmat;
  if(optRR==1) {
    double m01=dStar[i1]->mass0();
    double m02=dStar[i2]->mass0();
    rmat = A[8][iso]/(2*(m01*m01 + m02*m02 ));
  } else {
    rmat = ParR[8][iso][0]/pow(srt-ParR[8][iso][1],ParR[8][iso][2]);
  }
  //if(iso==1 && i1 == i2) rmat /=2;
  return spin*rmat*dStar[i1]->getBW(srt,i2-i1)/flux;
}

// N N -> D* D*
double SigmaRR::sigDstarDstar(double srt,double m1, double m2, int iso)
{
  double sig=0.0;
  for(int i=0;i<(int)dStar.size();i++) {
    for(int j=i;j<(int)dStar.size();j++) {
      sig += sigDstarDstar(srt,m1,m2,iso,i,j);
    }
  }
  return sig;
}



// D(1232) N -> N N
double SigmaRR::SigNDNN(double srt, double m1, double m2, int iso, int ip,
	double s12, double s34)
{
  int spin = 2*2;
  //double   rmat = ParDelta[0]/pow(srt-ParDelta[1],ParDelta[2]);
  double rmat = ParR[0][1][0]/pow(srt-ParR[0][1][1],ParR[0][1][2]);
  return spin*rmat*PCM(srt,Mnucl,Mnucl)/(PCM(srt,m1,m2)*srt*srt)*s12/s34;
}

// ip = 0) N N -> N D(1232)
//      1) N N -> N N*
//      2) N N -> N D*
//      3) N N -> D(1232) D(1232)
//      4) N N -> D(1232) N*
//      5) N N -> D(1232) D*
//      6) N N -> N* N*
//      7) N N -> N* D*
//      8) N N -> D* D*
// 1 2 -> 3 4   s12, s34 are symmetry factor: If 1 = 2, s12=2
double SigmaRR::SigDetbal(double srt, double m1, double m2, 
	double m01, double m02,
	int iso, int ip, double s12, double s34,int ir)
{
  int spin = 2*2;
  double rmat = matrixElm(ip,srt,iso,m01,m02,ir);
  return spin*rmat*PCM(srt,Mnucl,Mnucl)/(PCM(srt,m1,m2)*srt*srt)*s12/s34;

}

double SigmaRR::matrixElm(int ip, double srt, int iso,double m01, double m02,int ir)
{
  if(ip<=5) {
    if(ip==1 && ir==0) {
      return ParN1440[iso][0]/pow(srt-ParN1440[iso][1], ParN1440[iso][2]);
    } else {
      return ParR[ip][iso][0]/pow(srt-ParR[ip][iso][1],ParR[ip][iso][2]);
    }
  } else {
    if(optRR==1)
      return A[ip][iso]/(2*(m01*m01 + m02*m02 ));
    else 
      return ParR[ip][iso][0]/pow(srt-ParR[ip][iso][1],ParR[ip][iso][2]);
  }
}

double SigmaRR::prob(int ip, double srt,int iso,int id1, int id2)
{
  // NN-> ND 
  if(ip==0) {
    //double elm = ParDelta[0]/pow(srt-ParDelta[1],ParDelta[2]);
    double elm = ParR[0][1][0]/pow(srt-ParR[0][1][1],ParR[0][1][2]);
    return elm*nucleon->getBW(srt,0);

  // NN->NN*
  } else if(ip==1) {

    int i2 = jamTable->nStarID(id2) - 1;
    //double m0=nStar[i2]->mass0();
    //double elm = A[0][iso]/(2*(Mnucl*Mnucl + m0*m0 ));
    double elm = ParR[ip][1][0]/pow(srt-ParR[ip][1][1],ParR[ip][1][2]);
    if(i2==0)
      elm = ParN1440[1][0]/pow(srt-ParN1440[1][1], ParN1440[1][2]);

    if(iso==0) {
      if(i2==0)
        elm += ParN1440[0][0]/pow(srt-ParN1440[0][1], ParN1440[0][2]);
      else
        elm += ParR[ip][0][0]/pow(srt-ParR[ip][0][1],ParR[ip][0][2]);
    }
    return elm*nucleon->getBW(srt,i2+1);

  // 2) NN -> DD
  } else if(ip==2) {
    double a1=1.0, a0=0.0;
    if(iso==0) {a1=0.05; a0=0.025;}
    else if(iso==2) {a1=0.045; a0=0.025;}

    //double elm = (a1*A[3][1] + a0*A[3][0])/( 4*Mdelta*Mdelta );
    //return elm * d1232->getBW(srt,0);

    double elm1 = ParR[ip][1][0]/pow(srt-ParR[ip][1][1],ParR[ip][1][2]);
    double elm0 = ParR[ip][0][0]/pow(srt-ParR[ip][0][1],ParR[ip][0][2]);
    return (a1*elm1 + a0*elm0) * d1232->getBW(srt,0);

  //...3) NN -> ND*
  } else if(ip==3) {
    int i2=jamTable->dStarID(id2) - 1;
    //double em0=dStar[i2]->mass0();
    //double elm = A[1][1]/(2*(Mnucl*Mnucl + em0*em0 ));
    double elm = ParR[ip][1][0]/pow(srt-ParR[ip][1][1],ParR[ip][1][2]);
    int shift=1+nStar.size();
    return elm*nucleon->getBW(srt,i2+shift);

  //...4) NN -> DN*
  } else if(ip==4) {
    int i2=jamTable->nStarID(id2) - 1;
    //double em0=nStar[i2]->mass0();
    //double elm = A[3][1]/(2*(Mdelta*Mdelta + em0*em0 ));
    double elm = ParR[ip][1][0]/pow(srt-ParR[ip][1][1],ParR[ip][1][2]);
    return elm*d1232->getBW(srt,i2+1);

  //...5) NN -> DD*
  } else if(ip==5) {
    int i2=jamTable->dStarID(id2) - 1;
    //cout << " id2= "<< id2 << " i2= "<<i2 <<endl;
    double a1 = 1.0, a0 = 0.0;
    if(iso==0) {a1=0.05; a0=0.025;}
    else if(iso==2) {a1=0.045; a0=0.025;}
    //double m0=dStar[i2]->mass0();
    //double elm = (a1*A[4][1]+a0*A[4][0])/(2*(Mdelta*Mdelta + m0*m0 ));
    double elm1 = ParR[ip][1][0]/pow(srt-ParR[ip][1][1],ParR[ip][1][2]);
    double elm0 = ParR[ip][0][0]/pow(srt-ParR[ip][0][1],ParR[ip][0][2]);
    int shift=1+nStar.size();
    return (a1*elm1 + a0*elm0)*d1232->getBW(srt,i2+shift);

  //...6) NN -> N*N*
  } else if(ip==6) {
    int j1=jamTable->nStarID(id1) - 1;
    int j2=jamTable->nStarID(id2) - 1;
    int i1 = min(j1,j2);
    int i2 = max(j1,j2);
    if(optRR==1) {
      double m01=nStar[i1]->mass0();
      double m02=nStar[i2]->mass0();
      double elm = A[ip][1]/(2*(m01*m01 + m02*m02 ));
      if(iso==0) elm += A[ip][0]/(2*(m01*m01 + m02*m02 ));
      return elm*nStar[i1]->getBW(srt,i2-i1);
    } else {
      double elm = ParR[ip][1][0]/pow(srt-ParR[ip][1][1],ParR[ip][1][2]);
      if(iso==0) elm += ParR[ip][0][0]/pow(srt-ParR[ip][0][1],ParR[ip][0][2]);
      return elm*nStar[i1]->getBW(srt,i2-i1);
    }

  //...7) NN -> N*D*
  } else if(ip==7) {
    int i1=jamTable->nStarID(id1) - 1;
    int i2=jamTable->dStarID(id2) - 1;
    int shift=nStar.size()-i1;
    if(optRR==1) {
      double m01=nStar[i1]->mass0();
      double m02=dStar[i2]->mass0();
      double elm = A[ip][iso]/(2*(m01*m01 + m02*m02 ));
      return elm*nStar[i1]->getBW(srt,i2+shift);
    } else {
      double elm = ParR[ip][1][0]/pow(srt-ParR[ip][1][1],ParR[ip][1][2]);
      return elm*nStar[i1]->getBW(srt,i2+shift);
    }

  //...8) NN -> D*D*
  } else if(ip==8) {
    int j1=jamTable->dStarID(id1) - 1;
    int j2=jamTable->dStarID(id2) - 1;
    int i1 = min(j1,j2);
    int i2 = max(j1,j2);
    double a1 = 1.0, a0 = 0.0;
    if(iso==0) {a1=0.05; a0=0.025;}
    else if(iso==2) {a1=0.045; a0=0.025;}
    if(optRR==1) {
      double m01=dStar[i1]->mass0();
      double m02=dStar[i2]->mass0();
      double elm1 = A[ip][1]/(2*(m01*m01 + m02*m02 ));
      double elm0 = A[ip][0]/(2*(m01*m01 + m02*m02 ));
      return (a1*elm1+a0*elm0) * dStar[i1]->getBW(srt,i2-i1);
    } else {
      double elm1 = ParR[ip][1][0]/pow(srt-ParR[ip][1][1],ParR[ip][1][2]);
      double elm0 = ParR[ip][0][0]/pow(srt-ParR[ip][0][1],ParR[ip][0][2]);
      return (a1*elm1+a0*elm0) * dStar[i1]->getBW(srt,i2-i1);
    }
  } else {
    cout << "SigmaRR::prob ip= "<< ip <<endl;
    exit(1);
  }
}

//...0) NN -> ND
//...1) NN -> NN*
//...2) NN -> DD
//...3) NN -> ND*
//...4) NN -> DN*
//...5) NN -> DD*
//
//...6) NN -> N*N*
//...7) NN -> N*D*
//...8) NN -> D*D*
//...9) NN-> s-wave pion
void SigmaRR::sigma(double srt, double m1, double m2, double* sig, int iso)
{
  for(int i=0;i<10;i++) sig[i]=0.0;
  if(srt <= 2*Mnucl+Mpion) return;

  sig[1]=sigNNstar(srt,m1,m2,iso,0);
  sig[2]=sigDeltaDelta(srt,m1,m2,iso);
  sig[5]=sigDDstar(srt,m1,m2,iso,0);
  if(modelRR==2) {
    sig[6]=sigNstarNstar(srt,m1,m2,iso);
    sig[8]=sigDstarDstar(srt,m1,m2,iso);
  }
  if(iso==0) return;
  sig[0]=sigDelta(srt,m1,m2);
  sig[3]=sigNDstar(srt,m1,m2,0);
  sig[4]=sigDNstar(srt,m1,m2,0);
  if(modelRR==2) sig[7]=sigNstarDstar(srt,m1,m2);

}

}//end namespace jam2

//g++ SigmaRR.cxx -I.. -I/export/ynara/lib/pythia8/include

//#define MAIN 1

#ifdef MAIN
#include <iostream>
#include <iomanip>
int main() {

  using namespace jam2;
  int optRR=2, optDWid=2, optWid=3, optBW=3;
  int iso=1;
  double m1=SigmaRR::Mnucl,m2=SigmaRR::Mnucl;
  double sig[10],sig1[10];
  //SigmaRR *sigma = new SigmaRR("bwintjam2.dat",optRR,optDWid,optWid,optBW);
  SigmaRR *sigma = new SigmaRR("/export/ynara/lib/BWintjam2a.dat",optRR,optDWid,optWid,optBW);
  cout << "# norm= "<< sigma->getDelta()->norm() << endl;
  int nn=200;
  double smin=1.8;
  double smax=6.0;
  double ds=(smax-smin)/(nn-1);
  for(int i=0;i<nn;i++) {
    double  s = smin + ds*i;
    sigma->sigma(s,m1,m2,sig,iso) ;
    // pp -> d+ d+  2/5 *M_I=1 * 1/2 (final state is identical)
    // pp -> d0 d++ 3/5 *M_I=1
    // sigma(pp total) =  2/5 * 1/2 + 3/5 = 2/10 + 6/10 = 8/10
    sig[2] *= 4.0/5.0;

    if(iso==0) {
      sigma->sigma(s,m1,m2,sig1,1) ;
      for(int j=0;j<10;j++) sig[j] = 0.5*(sig[j]+sig1[j]);
    }


    double sigt=0.0;
    for(int j=0;j<10;j++) sigt += sig[j];
    cout << setw(8) << s
	<< setw(14) << sigt
	<< setw(14) << sig[0]
	<< setw(14) << sig[1]
	<< setw(14) << sig[2]
	<< setw(14) << sig[3]
	<< setw(14) << sig[4]
	<< setw(14) << sig[5]
	<< setw(14) << sig[6]
	<< setw(14) << sig[7]
	<< setw(14) << sig[8]
	<< setw(14) << sig[9]
	<<endl;
  }

  delete sigma;
  return 0;

}
#endif
