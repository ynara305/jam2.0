#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

#include <iomanip>
#include <sstream>
#include <ctime>


#include <jam2/JAM.h>
//#include <jam2/initcond/InitialCondition.h>
//#include <jam2/hadrons/JamStdlib.h>

using namespace std;
using namespace Pythia8;
using namespace jam2;

//void addhist(Analysis* ana);

//#include "book/Book1.h"
#include <Pythia8/Basics.h>

//The class 'UserInitialCondition' is in the initcond/InitialCondition.h
class MyInitialCondition: public UserInitialCondition
{
public:
  MyInitialCondition(Pythia8::Settings* s, JamParticleData* pd,Pythia8::Rndm* r)
    :UserInitialCondition(s,pd,r) { 
  // open file?
    }
  ~MyInitialCondition() {}
  void init();
  void generate(Collision* event,int mode=0);
};

void MyInitialCondition::init()
{
  //some initializations.
  // If energy of two-body collision is large, we need to initialize pythia8.
  // set maximum possible collision energy.
  double eCMmax=20.0;
  settings->parm("Beams:eCM",eCMmax);
}

// put your particle list for the initial condition for JAM.
void MyInitialCondition::generate(Collision* event,int )
{
  int n=100;
  // loop over all particles.
  for(int i=0; i<n;i++) {

    // Particle PDG ID.
    int id=2212;

    // find this particle in the JAM particle list.
    ParticleDataEntry* pa= jamParticleData->find(id);

    // set particle coordinate.
    double x = -10.0 + 20*rndm->flat();
    double y = -10.0 + 20*rndm->flat();
    double z = -10.0 + 20*rndm->flat();
    double t = 0.0;
    Vec4 r(x,y,z,t);

    // set particle momentum.
    double px = -10.0 + 20*rndm->flat();
    double py = -10.0 + 20*rndm->flat();
    double pz = -10.0 + 20*rndm->flat();

    // particle mass.
    double m=0.938;

    // energy.
    double e = sqrt(m*m + px*px + py*py + pz*pz);
    Vec4 p(px,py,pz,e);

    EventParticle* cp = new EventParticle(id,m,r,p,pa);
    cp->setPID(jamParticleData->pid(abs(id)));

    // compute decay time if it is resonance.
    double dect = jamParticleData->lifeTime(pa,m,e);
    cp->setLifeTime(t+dect);
    
    // put this particle into the particle list.
    event->setPList(cp);
  }

  

}


class ParticleMult
{
private:
  int nevent;
  int *npar;
  int xpar;
  const int nhist=18;
public:
  ParticleMult() {
    nevent=0;
    xpar=0;
    npar = new int [nhist];
    for(int i=0;i<nhist;i++) npar[i]=0;
  }
  void fill(JAM* jam);
  void print(int n);
};

void ParticleMult::fill(JAM* jam)
{
  list<EventParticle*>::const_iterator jp;
  list<EventParticle*>&  plist = jam->getEvent();
  nevent++;
  xpar += plist.size();
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
     case  3122 : ih=14; break;
     case  3212 : ih=8; break;
     case -3122 : ih=15; break;
     case -3212 : ih=9; break;
     case  3112 : ih=10; break; // Sigma-
     case -3112 : ih=11; break; // Sigma-bar
     case  3222 : ih=12; break; // Sigma+
     case -3222 : ih=13; break; // Sigma+bar
     case  310  : ih=16; break; // K_S0
     case  130  : ih=17; break; // K_L0
     default: ih=-1;
    }
    if(ih>=0) npar[ih]++;
  }
}
void ParticleMult::print(int oversample)
{
  double wei = 1.0/nevent/oversample;
  cout << "total particle = " << xpar*wei << endl;;
  cout << "proton = " << npar[0]*wei << endl;;
  cout << "anti-proton = " << npar[1]*wei << endl;;
  cout << "pi- = " << npar[2]*wei << endl;
  cout << "pi+ = " << npar[3]*wei << endl;
  cout << "pi0 = " << npar[4]*wei << endl;
  cout << "K-  = " << npar[5]*wei << endl;
  cout << "K+  = " << npar[6]*wei << endl;
  cout << "K0+K0bar= " << npar[7]*wei << endl;
  cout << "K_S0= " << npar[16]*wei << endl;
  cout << "K_L0= " << npar[17]*wei << endl;
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
    static const int    nhist1=13;
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
    string pa[nhist1]={"proton", "pbar", "pi-", "pi+","pi0","K-","K+","Lambda","netp","Sigma-","Sigma0","Sigma+","K0+K0bar"};

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
          case  311  : 
          case -311  : ih=12; break;  // K0
	  case  3122 : ih=7; break; // Lambda
	  case  3112 : ih=9; break; // Sigma-
	  case  3212 : ih=10; break; // Sigma0
	  case  3222 : ih=11; break; // Sigma+
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
  string pa[nhist1]={"p", "pbar", "pi-", "pi+","pi0","K-","K+","Lambda","netp","Sigma-","Sigma0","Sigma+","K0+K0bar"};
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

void outputPhaseSpace(ostream& ofs, int iev, int mass, JAM* jam,int opt)
{
  list<EventParticle*>&  plist = jam->getEvent();
  int nv = plist.size();
  int ncoll = jam->getNColl();
  //int ncollBB = jam->getGlauberColl();
  int ncollBB = jam->getNCollBB();
  int ncollG = jam->getGlauberColl();
  int npartG = jam->getGlauberParticipants();
  int nonpart=0;
  int nmeson=0;
  for(auto jp=plist.begin(); jp != plist.end(); jp++) {
    if(abs((*jp)->getNColl())==1 && (*jp)->getParent()==0) nonpart++;
    if((*jp)->baryon() == 0) nmeson++;
  }
  //int nbaryon = nv - nmeson;
  int npart = mass - nonpart;
  double b = jam->impactParameter();
  ofs << "# " << iev << setw(12) << nv << setw(8) << ncollG
      << setw(15) << npartG
      << setw(18) << b
      << setw(8) << npart
      << setw(8) << ncoll
      << setw(8) << ncollBB
      << endl;
  Vec4 ptot=0.0;
  for(auto& p: plist) {
    Vec4 r = opt == 1 ? p->getV() : p->getR();
    ofs << setw(3) << p->getStatus()
        << setw(12) << p->getID()
        << setw(15) << scientific << p->getNColl()
        << setw(15) << fixed << p->getMass()
        << scientific
	<< setw(16) << setprecision(8) << p->getPx()
        << setw(16) << setprecision(8) << p->getPy()
        << setw(16) << setprecision(8) << p->getPz()
        << setw(16) << setprecision(8) << p->getPe()
        << setw(16) << setprecision(8) << r[1]
        << setw(16) << setprecision(8) << r[2]
        << setw(16) << setprecision(8) << r[3]
        << setw(16) << setprecision(8) << r[0]
        //<< setw(16) << setprecision(8) << p->TimeLastColl()
        << setw(16) << setprecision(8) << p->getTf()
	<< endl;
    ptot += p->getP();
  }
  //cout << " ptot = " << ptot <<endl;
}

void xsec(string file, string outfile)
{
  JAM *jam = new JAM(file);
  jam->init();
  CrossSection *xsec = jam->xsection;
  JamParticleData* jamParticleData = jam->jamParticleData;
  //SigmaMB *sigmb = xsec->getMB();

  ofstream ofs(outfile.c_str());

  //int id1=111;  // pi0
  //int id1=211;  // pi+
  //int id1=-211;  // pi-
  //int id1=321;  // K+
  //int id1=-321;  // K-
  //int id1=2212; // proton
  int id1=-2212; // anti-proton
  //int id1=-2112; // anti-neutron
  //int id1=3112; // sigma-
  //int id1=3122; // lambda

  //int id2=2212; // proton
  int id2=2112; // neutron
 
  // Specify collision type.
  //int icltype=1; // BB collision.
  //int icltype=2; // MB collision.
  //int icltype=3; // MM collision.
  int icltype=4; // B~B collision.

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

  bool islog=false;
  islog=true;
  double smin= m1 + m2 + 0.0;
  double smax= 1000000;
  //double smax= 100;
  int nn=4000;
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

    //double sig = cpair.getSigma();
    //double sigel = cpair.getSigmaElastic();
    double sigelbw = cpair.getElasticBW();
    double sigabs = cpair.getSigAbs();
    double sigch = 0.0;
    if(icltype==4) {
      sigch = cpair.getChargeEx();
    }

    double plab=plabsr(srt,m1,m2);
    ofs << setw(5) << srt << setw(15)<< plab
	 << setw(15) << sig << setw(15) << sigel << setw(15) << sigelbw
	 << setw(15) << sigabs
	 << setw(15) << sigch
	 <<endl;
  }

  ofs.close();
}

string date()
{
      /*
      time_t t=time(nullptr);
      const tm* lt = localtime(&t);
      std::stringstream s;
      s<<"20";
      s<<lt->tm_year-100; //100を引くことで20xxのxxの部分になる
      s<<"-";
      s<<lt->tm_mon+1; //月を0からカウントしているため
      s<<"-";
      s<<lt->tm_mday; //そのまま
      //result = "2015-5-19"
      std::string result = s.str();
      */

    time_t timer;        /* 時刻を取り出すための型（実際はunsigned long型） */
    struct tm *local;    /* tm構造体（時刻を扱う */

    /* 年月日と時分秒保存用 */
    int year, month, day, hour, minute, second;
    timer = time(NULL);        /* 現在時刻を取得 */
    local = localtime(&timer);    /* 地方時に変換 */
    /* 年月日と時分秒をtm構造体の各パラメタから変数に代入 */
    year = local->tm_year + 1900;        /* 1900年からの年数が取得されるため */
    month = local->tm_mon + 1;        /* 0を1月としているため */
    day = local->tm_mday;
    hour = local->tm_hour;
    minute = local->tm_min;
    second = local->tm_sec;

    std::stringstream s;
    //s<< year << month << day << hour << minute << second;
    //s<< year <<"Y" << month << "M" << day << "D" << hour <<"H" ;

    s << year <<" " << setw(2) << setfill('0') << month
      << ":" << day << ":" << hour
      << ":"<< minute
      << ":"<< second;

    //cout << resetiosflags(ios_base::floatfield);

   cout << s.str()<<endl;

    //printf("%dy%dM%dD %dH%dM%dS\n", year, month, day, hour, minute, second);
    printf("%d/%d/%d %02d:%02d:%02d\n",year,month,day,hour,minute,second);

    return s.str();
}

int main(int argc, char* argv[]) {

  string outputFname="phase.dat";
  string inputFname="jam.inp";
  string outfile="xsec.dat";
  bool outputPhaseData=true;
  MyHist* myhist = 0;
  ParticleMult *pmult = new ParticleMult();
  bool hist=false;
  ofstream ofs;
  double ymax=5.0;
  for(int i=1; i<argc; i++) {
    if(!strcmp(argv[i],"-f")) inputFname = argv[i+1];
    if(!strcmp(argv[i],"-p")) outputPhaseData = atoi(argv[i+1]);
    if(!strcmp(argv[i],"-h")) hist = atoi(argv[i+1]);
    if(!strcmp(argv[i],"-y")) ymax = atof(argv[i+1]);
    if(!strcmp(argv[i],"-o")) outfile = argv[i+1];
  }
  //cout << " phase = " << outputPhaseData <<endl;
  //cout << " hist= "<< hist << " ymax= "<< ymax <<endl;

  clock_t start = clock();

  JAM *jam = new JAM(inputFname);

  // JAM initial condtion 
  InitialCondition* initcnd=0;

  // initial condition from user?
  if(jam->settings->mode("Cascade:initialCondition")==3) {
    initcnd= new MyInitialCondition(jam->settings,jam->jamParticleData,jam->rndm);
  }
  jam->init(initcnd);
  int nev = jam->mode("Main:numberOfEvents");
  int oversample = jam->getOverSample();
  int massAB = jam->initcond()->massA() + jam->initcond()->massB();
  if(outputPhaseData) {
    //ofs.open(outputFname.c_str() + date();
    ofs.open(outputFname.c_str());
    ofs << "# " << nev << setw(8) << jam->initcond()->ycm()
        << setw(8) << jam->initcond()->gammaA()
        << setw(8) << jam->initcond()->gammaB()
        << setw(8) << 3
        <<endl;
  }

  if(hist) {
    myhist = new MyHist(-ymax,ymax);
    myhist->init(jam->settings->parm("Beams:eCM"));
  }


/*
  JamParticleData* jamParticleData = jam->jamParticleData;
  Pythia8::ParticleDataEntry *pd1=jamParticleData->find(2112);
  EventParticle* i1= new EventParticle(2112,pd1);
  EventParticle* i2= new EventParticle(2112,pd1);
  double srt=4.0;
  double pr=2.9;
  int icltyp=1;
  double ctime=1.9, tcol1=0.1, tcol2=0.2, brel=0.3;
  CollisionPair cpair = CollisionPair(icltyp,i1,i2,srt,pr);
  TwoBodyInterList* it = new TwoBodyInterList(cpair,i1,i2,
	    ctime,tcol1,tcol2,brel);

  i1->setInterList(it);
  i2->setInterList(it);
  cout << " i1= "<< i1 << " i2= "<< i2<<endl;
  i1->printCollisionTime();
  i2->printCollisionTime();

  delete i1;
  delete i2;
  delete it;
  return 0;
  */

  int isMeanField = jam->settings->mode("MeanField:mode");

  // Loop over all event.
  for(int iev=1; iev<=nev; iev++) {
    jam->next();
    if(hist) myhist->fill(jam);
    pmult->fill(jam);
    if(outputPhaseData) outputPhaseSpace(ofs,iev,massAB,jam,isMeanField);

  } //end event loop

    if(hist) myhist->print(nev*oversample);
    pmult->print(oversample);

    cout << "   number of expected Glauber type initial collision= "
	<< fixed << jam->getNInter() << endl
	<< "   number of total   collision = " << jam->getXColl() << endl
	<< "   number of BB      collision = "<< jam->getXCollBB() << endl
	<< "   number of MB      collision = "<< jam->getXCollMB()<< endl
	<< "   number of MM      collision = "<< jam->getXCollMM() << endl
	<< "   number of BBar    collision = "<< jam->getXCollBBar() << endl
	<< "   number of BbarBar collision = "<< jam->getXCollBbarBar() << endl
	<< "   number of RR->NN  collision = "<< jam->getXCollRR2NN() <<endl
	<< "   number of elastic collision = "<< jam->getNElastic() << endl
	<< "   number of absorption = "<< jam->getNAbsorb() << endl
	<< "   number of decay = "<< jam->getNDecay() << endl
	<< "   number of Pauli Blocking = "<< jam->getNPauli()
	<< endl;

    if(outputPhaseData) {
      ofs.close();
      system(("gzip -f "+outputFname).c_str());
    }

    delete jam;

    date();
    const double time = static_cast<double>(clock() - start) / CLOCKS_PER_SEC;
    int isec = time;
    int imin=time/60;
    int ihrs=imin/60;
    imin -= ihrs*60;
    isec -= ihrs*3600+imin*60;
    printf("CPU time %lf[s]\n", time);
    printf("CPU time= %d h %d m %d s\n",ihrs,imin,isec);

    return 0;
}




