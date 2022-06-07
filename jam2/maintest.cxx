#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

#include <iomanip>
#include <sstream>
#include <ctime>

//#include <jam2/interaction/Scatter.h>

#include <jam2/JAM.h>
//#include <jam2/hadrons/JamStdlib.h>

using namespace std;
using namespace Pythia8;
using namespace jam2;

//void addhist(Analysis* ana);

//#include "book/Book1.h"
#include <Pythia8/Basics.h>
#include <Pythia8/SigmaTotal.h>



class ExclusiveReaction
{
private:
  //for pp
  static const int npdet; // number of stable particles
  static const int mpdet; // number of unstable particles
  static const int ncdet1; // number of exclusive reactions
  static const int mchan1[16][14], ktype[16];
  static const std::string rname1[16];

  // for pn
  static const int ncdet2; // number of exclusive reactions for pn
  static const int mchan2[15][14];
  static const std::string rname2[15];

  int ncdet,**mchan;
  std::vector<std::string> rname;
  std::string rtype;
  int *ncount, *xcount;
  double eCM, pLab;
public:
  ExclusiveReaction(double ecm, double plab,int mode);
  ~ExclusiveReaction();
  void analyze(JAM* jam);
  void print(ofstream& ofs1,ofstream& ofs2,ofstream& ofs3,double w,int mode);
};

const int  ExclusiveReaction::npdet=14;
const int  ExclusiveReaction::mpdet=16;
const int  ExclusiveReaction::ncdet1=16;
const int  ExclusiveReaction::ncdet2=15;

const int ExclusiveReaction::ktype[16]={
    2112,2212,3122,3112,3212,3222,
    -211,111,211,-321,-311,311,321,221,223,113};

//.... n  p  L  S- S0 S+ p- p0 p+ K- Kb K0 K+ eta
const int  ExclusiveReaction::mchan1[16][14]={
  {1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},   // np pi+
  {0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},   // pp pi0
  {1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0},   // np pi0 pi+
  {0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},   // pp pi- pi+
  {0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0},   // pp pi0 pi0
  {1, 1, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0},   // np pi- 2pi+
  {0, 2, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0},   // pp pi- pi0 pi+
  {0, 2, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0},   // pp 2pi- 2pi+
  {1, 1, 0, 0, 0, 0, 2, 0, 3, 0, 0, 0, 0, 0},   // np 2pi- 3pi+
  {0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0},   // pL  K+
  {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0},   // nS+ K+
  {0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0},   // pS0 K+
  {1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0},   // nL  pi+ K+
  {0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0},   // pL pi0 K+
  {0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},   // pp eta
  {0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0}};  // Lp pi+ K0

//.... n  p  L  S- S0 S+ p- p0 p+ K- Kb K0 K+ eta
const int ExclusiveReaction::mchan2[15][14]={
  {1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},   // np pi0
  {0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},   // pp pi-
  {1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},   // np pi- pi+
  {0, 2, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0},   // pp pi- pi0
  {1, 1, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0},   // np pi- 2pi+
  {0, 2, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0},   // pp pi- pi0 pi+
  {0, 2, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0},   // pp 2pi- 2pi+
  {1, 1, 0, 0, 0, 0, 2, 0, 3, 0, 0, 0, 0, 0},   // np 2pi- 3pi+
  {0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0},   // pL  K+
  {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0},   // nS+ K+
  {1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0},   // nL  pi+ K+
  {0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0},   // pS0 K+
  {0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0},   // pL pi0 K+
  {0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},   // pp eta
  {0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0}};  // Lp pi+ K0


const string ExclusiveReaction::rname1[16]={
        "p n pi+",  "p p pi0",
	"p n pi+ pi0", "p p pi+ pi-","p p pi0 pi0",
	"pn 2pi+ pi-",  "pp pi+ pi0 pi-",
        "pp 2pi+ 2pi-", "pn 3pi+ 2pi-",
        "LpK+",   "S+nK+", "S0pK+",  "Ln+K+",
        "Lp0K+",  "pp eta", "Lambda p pi+ K0"};

const string ExclusiveReaction::rname2[15]={
        "p n pi0",    "p p pion-",
        "p n pi^+ pi^-", "p p pi+ pi0",
        "pn 2pi+ pi-",  "pp pi+ pi0 pi-",
        "pp 2pi+ 2pi-", "pn 3pi+ 2pi-",
        "LpK+",   "S+nK+",
        "S0pK+",  "Ln+K+",
        "Lp0K+",  "pp eta",
        "Lambda p pi+ K0"};


ExclusiveReaction::ExclusiveReaction(double ecm, double plab, int mode)
{
  eCM=ecm;
  pLab=plab;

  // pp collision
  if(mode==1) {
    rtype  = "# pp exclusive";
    ncdet = ncdet1;
    for(int i=0;i<ncdet;i++) rname.push_back(rname1[i]);
    mchan = new int* [ncdet];
    for(int i=0;i<ncdet;i++) mchan[i] = new int [npdet];
    for(int i=0;i<ncdet;i++)
    for(int j=0;j<npdet;j++) { mchan[i][j] = mchan1[i][j]; }

  // pn collision
  } else {
    rtype  = "# pn exclusive";
    ncdet = ncdet2;
    mchan = new int* [ncdet];
    for(int i=0;i<ncdet;i++) rname.push_back(rname2[i]);
    for(int i=0;i<ncdet;i++) mchan[i] = new int [npdet];
    for(int i=0;i<ncdet;i++)
    for(int j=0;j<npdet;j++) { mchan[i][j] = mchan2[i][j]; }

  }

  ncount = new int[ncdet];
  for(int i=0;i<ncdet;i++) ncount[i]=0;

  xcount = new int[ncdet];
  for(int i=0;i<ncdet;i++) xcount[i]=0;

}

ExclusiveReaction::~ExclusiveReaction()
{
  for(int i=0;i<ncdet;i++) delete [] mchan[i];
  delete [] mchan;
  delete [] ncount;
  delete [] xcount;
}
void ExclusiveReaction::analyze(JAM* jam)
{
  list<EventParticle*>::const_iterator jp;
  list<EventParticle*>&  plist = jam->getEvent();

  int ntype[ncdet];
  for(int i=0;i<ncdet;i++) ntype[i]=0;

  // Loop over particles.
  for(jp=plist.begin(); jp != plist.end(); jp++) {
    int id  = (*jp)->getID();
    for(int ipdet=0;ipdet<ncdet;ipdet++) {
      if(id==ktype[ipdet]) {
	  ntype[ipdet]++;
	  xcount[ipdet]++;
      }
    }

  } // end loop over particle.

  int channel=-99;
  // Loop over reaction channel.
  for(int ic=0;ic<ncdet;ic++) {
    // Loop over particles involved this reaction channel.
    bool  next=false;
    for(int ipdet=0;ipdet<npdet;ipdet++) {
      if(ntype[ipdet] != mchan[ic][ipdet]) {next=true;break;}
    }
    if(next == true) continue;
    channel=ic;
    break;
  }

  if(channel >=0 && channel <ncdet) ncount[channel]++;

}

void ExclusiveReaction::print(ofstream& ofs1,ofstream& ofs2,ofstream& ofs3,double weight,int mode)
{

  if(mode==0) {
    ofs1 << "# (1) plab (2) eCM" <<endl;
    ofs2 << "# (1) plab (2) eCM" <<endl;
    ofs3 << rtype <<endl;
    for(int i=0;i<9;i++) ofs1 << "# (" << i+3 << ") " <<  rname[i] <<endl;
    for(int i=9;i<ncdet;i++) ofs2 << "# ("<< i-6 << ") " << rname[i] <<endl;
  }

  ofs1 << setw(5) << pLab << setw(12) << eCM;
  for(int i=0;i<9;i++) {
    ofs1 << setw(15) << ncount[i]*weight;
  }

  ofs1 << endl;

  ofs2 << setw(5) << pLab << setw(12) << eCM;
  for(int i=9;i<ncdet;i++) {
    ofs2 << setw(15) << ncount[i]*weight;
  }
  ofs2 << endl;

  ofs3 << setw(5) << pLab << setw(12) << eCM;
  for(int i=0;i<ncdet;i++) 
  ofs3 << setw(15) << xcount[i]*weight;
  ofs3 << endl;

}


class ParticleMult
{
private:
  int nevent;
  int *npar;
  const int nhist=16;
public:
  ParticleMult() {
    nevent=0;
    npar = new int [nhist];
    for(int i=0;i<nhist;i++) npar[i]=0;
  }
  void fill(JAM* jam);
  void print();
};

void ParticleMult::fill(JAM* jam)
{
  list<EventParticle*>::const_iterator jp;
  list<EventParticle*>&  plist = jam->getEvent();
  nevent++;
  for(jp=plist.begin(); jp != plist.end(); jp++) {
    int   id  = (*jp)->getID();
    int ih=-1;
    switch (id) {
     case  2212 : ih=0; break;
     case -2212 : ih=1; break;
     case -211  : ih=2; break;
     case  211  : ih=3; break;
     case  111  : ih=4; break;
     case -321  : ih=5; break;
     case  321  : ih=6; break;
     case  311  : 
     case -311  : ih=7; break;
     case  3122 : ih=8; break;
     case  3212 : ih=14; break;
     case -3122 : ih=9; break; // Lambda-bar
     case -3212 : ih=15; break;
     case  3112 : ih=10; break; // Sigma-
     case -3112 : ih=11; break; // Sigma-bar
     case  3222 : ih=12; break; // Sigma+
     case -3222 : ih=13; break; // Sigma+bar
     default: ih=-1;
    }
    if(ih>=0) npar[ih]++;
  }
}
void ParticleMult::print()
{
  double wei = 1.0/nevent;
  cout << "proton = " << npar[0]*wei << endl;;
  cout << "anti-proton = " << npar[1]*wei << endl;;
  cout << "pi- = " << npar[2]*wei << endl;
  cout << "pi+ = " << npar[3]*wei << endl;
  cout << "pi0 = " << npar[4]*wei << endl;
  cout << "K-  = " << npar[5]*wei << endl;
  cout << "K+  = " << npar[6]*wei << endl;
  cout << "K0+K0bar= " << npar[7]*wei << endl;
  cout << "Lambda= " << npar[8]*wei << endl;
  cout << "Lambdabar= " << npar[9]*wei << endl;
  cout << "Sigma-= " << npar[10]*wei << endl;
  cout << "anti-Sigma+= " << npar[11]*wei << endl;
  cout << "Sigma+= " << npar[12]*wei << endl;
  cout << "ani-Sigma+= " << npar[13]*wei << endl;
  cout << "Sigma0= " << npar[14]*wei << endl;
  cout << "ani-Sigma0= " << npar[15]*wei << endl;
}

class MyHist
{
private:
    double ymin,ymax,xmin,xmax;
    int    ny,nx;
    int    np;
    static const int    nhist1=9;
    double pmin,pmax;
    double dy, dx,dp;
    double eCM, pCM;
    int nhsit1;
    Hist  *hist1;
    Hist  *hist2, *hist3, *hist4;
    Hist  *histc, *histpt;
public:
    MyHist(double y1=-5.0, double y2=5.0);
    ~MyHist();
    void init(double ecm);
    void fill(JAM* event);
    void print(int nev);
};

MyHist::MyHist(double y1, double y2)
{
    ymin =  y1;
    ymax =  y2;
    ny = 70;
    dy = (ymax - ymin)/ny;

    xmin=0.0;
    xmax=1.0;
    nx = 50;
    dx = (xmax - xmin)/nx;

    np = 30;
    pmin = 0.0;
    pmax = 3.0;
    dp = (pmax - pmin ) /np;
}

MyHist::~MyHist()
{
    delete [] hist1;
    delete [] hist2;
    delete [] hist3;
    delete [] hist4;
    delete [] histc;
    delete [] histpt;
}

void MyHist::init(double ecm)
{
    string pa[nhist1]={"proton", "pbar", "pi-", "pi+","pi0","K-","K+","Lambda","netp"};

    eCM=ecm;
    double mA=0.938, mB=0.938;
    pCM  = 0.5 * sqrt( (eCM + mA + mB)*(eCM - mA - mB)*(eCM - mA + mB)*(eCM + mA - mB) ) / eCM;
    histc = new Hist[2];
    histc[0].book("dndeta charged",ny,ymin,ymax);
    histc[1].book("dndeta NSD charged",ny,ymin,ymax);

    hist1 = new Hist[nhist1];
    hist2 = new Hist[nhist1];
    hist3 = new Hist[nhist1];
    hist4 = new Hist[nhist1];
    histpt = new Hist[nhist1];
    for(int i=0;i<nhist1;i++) {
      hist1[i].book("dndy"+pa[i],ny,ymin,ymax);
      hist2[i].book("dndpt"+pa[i],np,pmin,pmax);
      hist3[i].book("dndmt"+pa[i],np,pmin,pmax);
      hist4[i].book("dndxF"+pa[i],nx,xmin,xmax);
      histpt[i].book("<pt> as a function of x_F"+pa[i],nx,xmin,xmax);
    }
}

void MyHist::fill(JAM* jam)
{
    Vec4 ptot=0.0;
    list<EventParticle*>::const_iterator jp;
    list<EventParticle*>&  plist = jam->getEvent();

    //cout << " JAM= " << plist.size() <<endl;

    for(jp=plist.begin(); jp != plist.end(); jp++) {
	int   id  = (*jp)->getID();
	double rap = (*jp)->getP().rap();
	double eta = (*jp)->getP().eta();
	double pt = (*jp)->getP().pT();
	double m = (*jp)->getMass();
	double mt = sqrt(m*m + pt*pt);
	double xf = (*jp)->getP().pz() /pCM;
	int ischarged = (*jp)->charge() != 0 ? 1 : 0;

	if(ischarged) histc[0].fill(eta,1.0/dy);

	ptot += (*jp)->getP();

	int ih=-1;
	double w = 1.0;
	switch (id) {
	  case  2212 : ih=0; break;
	  case -2212 : ih=1;w=-1.0; break;
	  case -211  : ih=2; break;
	  case  211  : ih=3; break;
	  case  111  : ih=4; break;
	  case -321  : ih=5; break;
	  case  321  : ih=6; break;
	  case  3122 : ih=7; break;
	  //case  3212 : ih=7; break; // include Sigma0 for Lambda
	  default: ih=-1;
	}
		       
	if(ih>=0) {
	  hist1[ih].fill(rap,1.0/dy);
	  hist2[ih].fill(pt,1.0/dp/pt/2.0);
	  //if(abs(rap)<0.1) hist3[ih].fill(mt-m,1.0/dp/mt/0.2);
	  if(rap>0.0 && rap<0.2) hist3[ih].fill(mt-m,1.0/dp/mt/0.2);
	  hist4[ih].fill(xf,1.0/dx);
	  histpt[ih].fill(xf,pt/dx);

	  if(ih==1) {
	    hist1[8].fill(rap,w/dy);
	    hist2[8].fill(pt,w/dp/pt/2.0);
	    //if(abs(rap)<0.1) hist3[8].fill(mt-m,w/dp/mt/0.2);
	    if(rap>0.0 && rap<0.2) hist3[8].fill(mt-m,w/dp/mt/0.2);
	    hist4[8].fill(xf,w/dx);
	    histpt[8].fill(xf,w*pt/dx);
	  }
	}

	// take care of net-proton.
	if(ih==0) {
	  hist1[8].fill(rap,1.0/dy);
	  hist2[8].fill(pt,1.0/dp/pt/2.0);
	  //if(abs(rap)<0.1) hist3[1].fill(mt-m,1.0/dp/mt/0.2);
	  if(rap>0.0 && rap<0.2) hist3[8].fill(mt-m,1.0/dp/mt/0.2);
	  hist4[8].fill(xf,1.0/dx);
	  histpt[8].fill(xf,pt/dx);
	}

    }

    if(ptot.pAbs()>1e-5) {
      cout << " ptot= " << ptot << " abs= " << ptot.pAbs() <<endl;
    }
    //cin.get();
}

void MyHist::print(int nev)
{
  string Dir="hist";
  string pa[nhist1]={"p", "pbar", "pi-", "pi+","pi0","K-","K+","Lambda","netp"};
  struct stat st;
  if(stat(Dir.c_str(),&st) !=0) mkdir(Dir.c_str(),0775);

  for(int i=0;i<nhist1;i++) {
    Hist *h = new Hist("<pt(x_F)> "+pa[i],nx,xmin,xmax);
    *h = histpt[i] / hist4[i];
    h->table(Dir+"/pt_"+pa[i]+".dat");
    delete h;
  }

    double wei = 1.0/(nev);
    for(int i=0;i<nhist1;i++) {
      hist1[i] *= wei;
      hist2[i] *= wei;
      hist3[i] *= wei;
      hist4[i] *= wei;
      histpt[i] *= wei;
    }
    histc[0] *= wei;
    histc[1] *= wei;
    histc[0].table(Dir+"/dndeta.dat");

    for(int i=0;i<nhist1;i++) {
      hist1[i].table(Dir+"/dndy_"+pa[i]+".dat");
      hist2[i].table(Dir+"/dndpt_"+pa[i]+".dat");
      hist3[i].table(Dir+"/dndmt_"+pa[i]+".dat");
      hist4[i].table(Dir+"/dndx_"+pa[i]+".dat");
    }

}

int npar[4];

void jamee(JAM* jam)
{
    double ecm=4.93;
    Pythia* pythia = jam->getHadronize();
    Event& event = pythia->event;
    event.reset();
    int    id1 = 2;
    int    id2 = 2101;
    double m1 = jam->particleData->m0(id1);
    double m2 = jam->particleData->m0(id2);
    //double pp = sqrtpos(ee*ee - mm*mm);
    double  pp =PCM(ecm,m1,m2);
    double ee1 = sqrt(m1*m1 + pp*pp);
    double ee2 = sqrt(m2*m2 + pp*pp);
    event.append(  id1, 23, 101,   0, 0., 0.,  pp, ee1, m1);
    event.append(  id2, 23,   0, 101, 0., 0., -pp, ee2, m2);
    pythia->next();

    for(int i=0; i<event.size();i++) {
	if(!event[i].isFinal()) continue;
	int   id  = event[i].id();
	//double rap = event[i].p().rap();
	//double pt = event[i].p().pT();
	if(id==2212) {
	    npar[0]++;
	} else if(id==-2212) {
	    npar[1]++;
	} else if(id == - 211) {
	    npar[2]++;
	} else if(id == 211) {
	    npar[3]++;
	}
    }

}

void jamevent() {

    //ofstream ofs("PCM.INFO");

    MyHist* myhist = new MyHist();

    int nev=10;
    JAM *jam = new JAM();
    jam->init();
    nev = jam->mode("Main:numberOfEvents");
    myhist->init(jam->settings->parm("Beams:eCM"));

    for(int iev=1; iev<=nev; iev++) {

	jam->next();
	myhist->fill(jam);

    } //end event loop

    myhist->print(nev);

    delete jam;
}

void jamevent_ee()
{
    npar[0]=npar[1]=npar[2]=npar[3]=0;

    int nev=10;

    JAM *jam = new JAM();
    jam->init();
    nev = jam->mode("Main:numberOfEvents");


    for(int iev=1; iev<=nev; iev++) {

	//jam->next();
	jamee(jam);

	//myhist->fill(jam);

    } //end event loop

    //myhist->print(nev);

    double wei=1.0/nev;
    cout << "proton= " << npar[0]*wei << endl;;
    cout << "anti-proton= " << npar[1]*wei << endl;;
    cout << "pi-= " << npar[2]*wei << endl;;
    cout << "pi+= " << npar[3]*wei << endl;;
}


void xsec(string file, string outfile, string outfile2,int opt)
{
  JAM *jam = new JAM(file);
  jam->settings->mode("Cascade:model",1);
  jam->init();
  CrossSection *xsec = jam->xsection;
  JamParticleData* jamParticleData = jam->jamParticleData;
  //SigmaMB *sigmb = xsec->getMB();

  ofstream ofs(outfile.c_str());
  ofstream ofs2(outfile2.c_str());

  // Specify collision type.
  //int icltype=1; // BB collision.
  int icltype=2; // MB collision.
  //int icltype=3; // MM collision.
  //int icltype=4; // B~B collision.

  //int id1=2224; // D++
  //int id1=2214; // D+
  //int id1=2212; // proton
  //int id1=111;  // pi0
  //int id1=211;  // pi+
  int id1=-211;  // pi-
  //int id1=321;  // K+
  //int id1=-321;  // K-
  //int id1=-2212; // anti-proton
  //int id1=-2112; // anti-neutron
  //int id1=3112; // sigma-
  //int id1=3122; // lambda

  int id2=2212; // proton
  if(opt==2) id2=2112; // neutron
 

  ParticleDataEntry* pd1 =jamParticleData->find(id1);
  ParticleDataEntry* pd2 =jamParticleData->find(id2);
  int pid1=jamParticleData->pid(id1);
  int pid2=jamParticleData->pid(id2);
  EventParticle* pa1 = new EventParticle(id1,pd1);
  EventParticle* pa2 = new EventParticle(id2,pd2);
  pa1->setPID(pid1);
  pa2->setPID(pid2);
  double m1 = pd1->m0();
  double m2 = pd2->m0();
  pa1->setMass(m1);
  pa2->setMass(m2);

  double smin= m1 + m2 + 0.0;
  bool islog=false;
  //islog=true;
  //double smax= 1000000;
  double smax= 10.0;
  int nn=100;
  //int nn=4000;
  double g = 1.0/nn*log10(smax/smin);
  //double ds= (smax - smin)/nn; 
  double ds=0.005;
  if(!islog) nn = int((smax-smin)/ds);

  ofs << "#  m1= " << m1 << " m2= " << m2 <<endl;
  for(int i=1;i<=nn;i++) {
    double srt = islog ? smin*pow(10,i*g) :  m1 + m2 + i*ds;
    double s = srt * srt;
    double pr=sqrt((s-(m1+m2)*(m1+m2))*(s-(m1-m2)*(m1-m2))/(4*s));
    CollisionPair cpair = CollisionPair(icltype,pa1,pa2,srt,pr);

    double sig = xsec->sigma(cpair);
    double sigel = xsec->sigmaEl();

    // inelastic cross section.
    if(srt > 0.938*2+0.138 + 0.02) {
      cpair.setMode(3);
      xsec->sigma(cpair);
    }

    //double sig = cpair.getSigma();
    //double sigel = cpair.getSigmaElastic();
    double sigelbw = cpair.getElasticBW();
    double sigabs = cpair.getSigAbs();
    double sigch = 0.0;
    if(icltype==4) {
      sigch = cpair.getChargeEx();
    }

    double plab=plabsr(srt,m1,m2);
    ofs << setw(5) << srt
	<< setw(15)<< plab
	 << setw(15) << sig
	 << setw(15) << sigel
	 << setw(15) << sigelbw
	 << setw(15) << sigabs
	 << setw(15) << sigch
	 <<endl;

    // output inelastic cross sections.
    ofs2 << setw(5) << srt << setw(15) << plab;
    double sigin[20]={0.0};
    double sigt=0.0;
    for(int j=0;i<cpair.outgoingSize();j++) {
      if(cpair.getOutGoing(j).id1 !=92) {
      sigin[j] = cpair.getOutGoing(j).sig;
      sigt += sigin[j];
      }
    }
    ofs2 << setw(15) << sigt;
    for(int j=0;j<20;j++) {
      ofs2 << setw(15) << sigin[j];
    }
    ofs2 << endl;

  }

  ofs.close();
  ofs2.close();
}

void exclusive(string inputFname,int opt)
{
  double pmin=0.9, pmax=14.0;
  int npmax=50;

  pmin=0.9;pmax=50; npmax=100;
  //pmin=0.9;pmax=15; npmax=50;

  //double dp = (pmax-pmin)/npmax;
  int nev=30000;
  ofstream ofs1, ofs2, ofs3;
  int id1=2212;
  int id2=2212;
  int icltype=1; // BB collision.
  string beam2 = "p+";
  if(opt==1) {
    ofs1.open("pppiex.dat");
    ofs2.open("ppstrangeex.dat");
    ofs3.open("ppX.dat");
  } else {
    ofs1.open("pnpiex.dat");
    ofs2.open("pnstrangeex.dat");
    ofs3.open("pnX.dat");
    id2=2112;
    beam2 = "n0";
  }

  // Loop over beam energy.
  for(int ip=0;ip<npmax;ip++) {
    JAM *jam = new JAM(inputFname,"dummy",false);

    double plab=pmin*pow(pmax/pmin, ip/(double)(npmax-1));
    jam->settings->mode("Main:numberOfEvents",nev);
    jam->settings->parm("Beams:eLab",0.0);
    jam->settings->parm("Beams:pLab",plab);
    jam->settings->word("Beams:beamA","p+");
    //jam->settings->word("Beams:beamB","p+");
    jam->settings->word("Beams:beamB",beam2);
    jam->settings->parm("Beams:bmin",0.0);
    jam->settings->parm("Beams:bmax",0.0);
    jam->settings->mode("Check:Debug",0);
    jam->settings->mode("Cascade:PrintCollision",0);
    jam->settings->parm("Beams:zseparation",1.0);
    jam->settings->flag("Cascade:InelasticOnly",true);
    jam->settings->flag("111:mayDecay", false);   // pi0
    jam->settings->flag("3122:mayDecay", true);   // Lambda

    jam->init();
    
    double ecm = jam->parm("Beams:eCM");
    ExclusiveReaction* ana = new ExclusiveReaction(ecm,plab,opt);

    CrossSection *xsec = jam->xsection;
    JamParticleData* jamParticleData = jam->jamParticleData;
    ParticleDataEntry* pd1 =jamParticleData->find(id1);
    ParticleDataEntry* pd2 =jamParticleData->find(id2);
    double m1 = pd1->m0();
    double m2 = pd2->m0();

    double s = ecm * ecm;
    double pr=sqrt((s-(m1+m2)*(m1+m2))*(s-(m1-m2)*(m1-m2))/(4*s));
    EventParticle pa1 = EventParticle(id1,pd1);
    EventParticle pa2 = EventParticle(id2,pd2);
    int pid1=jamParticleData->pid(id1);
    int pid2=jamParticleData->pid(id2);
    pa1.setPID(pid1);
    pa2.setPID(pid2);
    pa1.setMass(m1);
    pa2.setMass(m2);
    CollisionPair cpair = CollisionPair(icltype,&pa1,&pa2,ecm,pr);
    double sig = xsec->sigma(cpair);
    double sigel = xsec->sigmaEl();
    double sigin = sig - sigel;

    cout << "id1= "<< id1 << " id2= "<< id2
	<< " plab= "<< jam->parm("Beams:pLab")
         << " ecm= "<< ecm
         << " sig= "<< sig
         << " sigel= "<< sigel
         << " sigin= "<< sigin
	 <<endl;

    for(int iev=1; iev<=nev; iev++) {
      jam->next();
      ana->analyze(jam);
    } //end event loop

    ana->print(ofs1,ofs2,ofs3,sigin/nev,ip);
    delete jam;
    delete ana;

  } // end loop over beam energy.

  ofs1.close();
  ofs2.close();
  ofs3.close();

}

void sigmatot()
{
  using namespace Pythia8;
  Pythia8::Pythia *pythia = new Pythia8::Pythia();
  cout << " xmldir= "<< pythia->settings.word("xmlPath") <<endl;
  SigmaTotal sigtot;
  sigtot.init(&pythia->info,pythia->settings,&pythia->particleData,&pythia->rndm);
  double smin=4.9, smax=10.0;
  int ns=200;
  double ds=(smax-smin)/ns;
  for(int i=0;i<ns;i++) {
    double srt=smin+i*ds;
    sigtot.calc(2212,2212,srt);
    double sigin=sigtot.sigmaND()+sigtot.sigmaAX()*2+sigtot.sigmaXX();
    double sigdif=sigtot.sigmaAX()*2+sigtot.sigmaXX();
    cout << setw(6) << srt
	 << setw(12) << sigtot.sigmaTot()
	 << setw(12) << sigtot.sigmaND()
	 << setw(12) << sigtot.sigmaAX()
	 << setw(12) << sigtot.sigmaXB()
	 << setw(12) << sigtot.sigmaXX()
	 << setw(12) << sigtot.sigmaAXB()
	 << setw(12) << sigtot.sigmaND()/sigin
	 << setw(12) << (sigtot.sigmaAX()*2 + sigtot.sigmaXX())/sigin
	 << setw(12) << sigtot.sigmaAX()*2/sigdif
	 <<endl;
  }
}


void printCPU(clock_t start)
{
  const double time = static_cast<double>(clock() - start) / CLOCKS_PER_SEC;
  int isec = time;
  int imin=time/60;
  int ihrs=imin/60;
  imin -= ihrs*60;
  isec -= ihrs*3600+imin*60;
  printf("CPU time %lf[s]\n", time);
  printf("CPU time= %d h %d m %d s\n",ihrs,imin,isec);
}

int main(int argc, char* argv[]) {

  int mode = 1;
  string outputFname="phase.dat";
  string inputFname="jam.inp";
  string outfile="xsec.dat";
  string outfile2="xsec2.dat";
  //MyHist* myhist = 0;
  //ParticleMult *pmult = new ParticleMult();
  //bool hist=false;
  ofstream ofs;
  //double ymax=5.0;
  int optpp=1;

  for(int i=1; i<argc; i++) {
    if(!strcmp(argv[i],"-f")) inputFname = argv[i+1];
    if(!strcmp(argv[i],"-o")) outfile = argv[i+1];
    if(!strcmp(argv[i],"-m")) mode = atoi(argv[i+1]);
    if(!strcmp(argv[i],"-r")) optpp = atoi(argv[i+1]);
    //if(!strcmp(argv[i],"-y")) ymax = atof(argv[i+1]);
  }

  clock_t start = clock();

  switch (mode) {
   case  1 : xsec(inputFname,outfile,outfile2,optpp); break;
   case  2 : exclusive(inputFname,optpp); break;
   case  3 : sigmatot(); break;
   default: cout << " wrong mode= "<< mode<<endl;exit(1);
  }

  if( mode == 1) {
    
  }

  printCPU(start);

    return 0;
}



