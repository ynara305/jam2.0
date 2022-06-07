#ifndef jam2_xsection_SigmaRR_h
#define jam2_xsection_SigmaRR_h

#include <Pythia8/PythiaStdlib.h>
#include <jam2/hadrons/JamParticleData.h>

namespace jam2 {

class SigmaRR;


class SHadron
{
private:
  int idp;        //particle id
  int spin;      // spin= 2*J+1
  double mass;
  double gam;
  double* width; // table of width as a function of mass
  std::vector<double*> BWint;
  double Norm;
  double minMass,maxMass, dMass;
  int   nMass;
  double lam2, facF;

  int nDecayMode;
  double *mass1, *mass2, *branch;
  int    *spin1, *spin2, *angMom;
public:
  SHadron(int j=0,int s=0,double m=0,double w=0,double mm=0,
  int ndec=0, double* br=0, double* m1=0, double* m2=0,
  int* s1=0,  int* s2=0, int* l=0);
  ~SHadron();
  double mass0()    const {return mass;}
  double width0()   const {return gam;}
  double mMin() const {return minMass;}
  int    getSpin()    const {return spin;}
  int    getIP()      const {return idp;}
  double deltaWidth(double emd,int optw);
  void   makeDeltaWidthTable(int optw);
  double formBlattWeisskopf(double m0, double m, double pf0, double pf, int l); 
  double BlattWeisskopf(double x,int l);
  double getTotalWidth(double emd,int opt);
  void   makeWidthTable(const int optW);
  //void   setWidth(int i,double r) {width[i]=r;}
  double norm() {return Norm;}
  void   computeNorm(int NG, double* xg, double* wg,int optbw,int optwidth);
  void   setBWint(int n, double* bw) {
    double* bwint=new double [n];
    for(int i=0;i<n;i++) bwint[i]=bw[i];
    BWint.push_back(bwint);
  }
  int bwsize() {return BWint.size();}
  double getBW(double srt,int idp);
  double getWidth(double m) {
    //int i=(int)floor((m-minMass)/dMass);
    int i=(int)((m-minMass)/dMass);
    if(i>=0 && i<nMass-1) {
      double mi=minMass+i*dMass;
      double s = (m-mi)/dMass;
      return width[i]*(1.0-s) + s*width[i+1];
    } else if(i>=maxMass) {
      std::cout << "SHadron::getWidth mass too large "<< m <<std::endl;
      exit(1);
      //return width[nMass-1];
    } else {
      std::cout << "SHadron::getWidth mass too small "<< m
	        << " minMass= "<< minMass
	        << " i= "<< i
                << " floor= " << (int)((m-minMass)/dMass)
	        <<std::endl;
      exit(1);
      //return width[0];
    }
  }

};

class SigmaRR
{
private:
  // matrix elements.
  static const double A1[9][2],A2[9][2],A3[9][2];
  static const double ParR1[9][2][3];
  static const double ParR2[9][2][3];
  static const double ParN1440_1[2][3], ParN1440_2[2][3];
  double A[9][2];
  double ParR[9][2][3], ParN1440[2][3];

  SHadron *nucleon, *d1232;
  std::vector<SHadron*> nStar, dStar;

  int optDeltaWidth, optWidth, optBW;
  double *xg, *wg;
  int NG;
  int modelRR,optRR;
  JamParticleData* jamTable;

public:
  static const double Mnucl, Mpion,Mdelta,Wdelta, eKinmi, minMass;
  static const double sMin,sMax,dS;
  static const int    nS;

  SigmaRR(std::string bwfile,int optrr,int opt1=1, int opt2=2, int optbw=3,JamParticleData* jp=0);
  ~SigmaRR();
  void makeNstarTable();
  void makeDstarTable();
  double BWint(double srt, double em2, SHadron* pa);
  double BWint2(double srt,SHadron* pa1,SHadron* pa2);
  bool readBWTable(std::string fname);
  void makeBWTable(std::string fname);
  double sigDelta(double srt, double m1, double m2);
  double sigDelta2(double srt, double m1, double m2);
  double sigNNstar(double srt,double m1, double m2, int iso, int ip);
  double sigNDstar(double srt,double m1, double m2, int ip);
  double sigDeltaDelta(double srt,double m1, double m2,int iso);
  double sigDNstar(double srt,double m1, double m2, int ip);
  double sigDDstar(double srt,double m1, double m2, int iso, int ip);

  double sigNstarNstar(double srt,double m1, double m2,int iso,int i1,int i2);
  double sigNstarNstar(double srt,double m1, double m2, int iso);
  double sigNstarDstar(double srt,double m1, double m2,int i1,int i2);
  double sigNstarDstar(double srt,double m1, double m2);
  double sigDstarDstar(double srt,double m1, double m2,int iso,int i1,int i2);
  double sigDstarDstar(double srt,double m1, double m2, int iso);

  double SigNDNN(double srt,double m1,double m2,int iso,int ip,
	  double s12, double s34);
  double SigDetbal(double srt, double m1, double m2,double m01, double m02,
	int iso, int ip, double s12, double s34,int ir=1);
  double matrixElm(int ip,double srt,int iso=1,double m01=1.232,double m02=1.232,int ir=1);
  double prob(int ip, double srt,int iso=0,int i1=0, int i2=0);
  void sigma(double srt, double m1, double m2, double* sig, int iso);
  SHadron* getDelta() {return d1232;}
  SHadron* getNucleon() {return nucleon;}
  SHadron* getNstar(int i) {return nStar[i];}
  SHadron* getDstar(int i) {return dStar[i];}

  static double BW(double m, double gam, double emr,int optBW) {
    using Pythia8::pow2;
    if (optBW==1)
      return gam*0.5/(pow2(m-emr) + 0.25*gam*gam)/M_PI;
    else if (optBW==2)
      return 2.0/M_PI*m*emr*gam/(pow2(m*m-emr*emr)+pow2(emr*gam));
    else
      return 2.0/M_PI*m*m*gam/(pow2(m*m-emr*emr)+pow2(m*gam));
  }

};
}
#endif
