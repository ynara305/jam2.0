#include <cmath>
#include <jam2/xsection/SigmaNDelta.h>
#include <jam2/hadrons/GaussPoints.h>
#include <jam2/hadrons/JamStdlib.h>
#include <Pythia8/PythiaStdlib.h>

using namespace std;

namespace jam2 {

using Pythia8::pow3;
using Pythia8::pow4;

const double SigmaNR::emnuc=0.9383,SigmaNR::empion=0.138,SigmaNR::ekinmi=0.0001;
const double SigmaNR::emr=1.232, SigmaNR::gamr=0.117, SigmaNR::gr=emr*gamr;
const int SigmaNR::ldec=1;

SigmaNR::SigmaNR(int opt)
{

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

  optWidth=opt;
  // ppres= 0.227895;
  prres=PCM(emr,emnuc,empion);
  facF = 0.2;  // Frankfurt
  //facF = pow2(prres/HBARC);
  //cout << "# facF = "<< facF <<endl;
  lam2=HBARC*HBARC;
  Norm=computeNorm();
}

SigmaNR::~SigmaNR()
{
  delete [] xg;
  delete [] wg;
}

double SigmaNR::deltaWidth(double emd)
{
  const double pscal1= 0.238, pscal2= 0.318, p0ref=0.227, widdlt=0.12;
  if(emd <= emnuc+empion+ekinmi) return 0.0;
  double pr=PCM(emd,emnuc,empion);

  if(optWidth==1) {
      return gamr*pow(pr/prres,2*ldec+1)*(emr/emd)
      *(1.0 + facF)/(1.0+facF*pow(pr/prres,2*ldec));

  } else if(optWidth==2)  {
      return gamr*pow(pr/prres,2*ldec+1)*(emr/emd)
	  *pow((prres*prres+lam2)/(pr*pr+lam2),ldec);

  } else {
// Randrup NP A314 (1979) 429.Rittenberg, REV.MOD.PHYS. 43 (1971) S1.
    return widdlt*(pr*pr*pr/(1.0+pow2(pr/pscal1)+pow4(pr/pscal2)))
         /(pow3(p0ref)/(1.0+pow2(p0ref/pscal1)+pow4(p0ref/pscal2)));
  }

}

double SigmaNR::computeNorm()
{
  //double emin=emr;
  double  emin=emnuc+empion;
  //double ymin=atan((emin*emin-emr*emr)/gr);
  //double ymax=M_PI/2;
  //double ymin=2*atan((emin*emin-emr*emr)/gr);
  //double ymax=M_PI;
  double ymin=2*atan(2*(emin-emr)/gamr);
  double ymax=M_PI;

  double dy = ymax - ymin;
  double bwint=0.0;
  for(int i=0;i<NG;i++) {
    double y= ymin + xg[i]*dy;

    //double emd=sqrt(tan(y)*gr+ emr*emr);
    //double dc=2*emd*gr/( pow2(emd*emd-emr*emr)+pow2(emr*gamr));
    //double dc=(1+cos(2*y))*emd/gr;

    //double emd=sqrt(tan(0.5*y)*gr+ emr*emr);
    //double dc=2*emd*(1+cos(y))/gr;

    double emd=gamr/2*tan(0.5*y)+emr;
    double dc= gamr/( 0.25*gamr*gamr + pow2(emd - emr) );

    double gam=deltaWidth(emd);
    //double gam=gamr;

    //double bw= gam*0.5/(pow2(emd-emr) + 0.25*gam*gam)/M_PI;
    //double bw=2.0/M_PI*emd*emr*gam/(pow2(emd*emd-emr*emr)+pow2(emr*gam));
    double bw=2.0/M_PI*emd*emd*gam/(pow2(emd*emd-emr*emr)+pow2(emd*gam));
    bwint += bw*dy*wg[i]/dc;
  }

  return bwint;

}

double SigmaNR::sigma(double srt)
{
  if(srt <= 2*emnuc+empion) return 0.0;

  double  emin=emnuc+empion;
  double  emax=srt-emnuc;
  //double  ymax=atan((emax*emax-emr*emr)/gr);
  //double  ymin=atan((emin*emin-emr*emr)/gr);
  //double  ymax=2*atan((emax*emax-emr*emr)/gr);
  //double  ymin=2*atan((emin*emin-emr*emr)/gr);
  double  ymax=2*atan(2*(emax - emr)/gamr);
  double  ymin=2*atan(2*(emin - emr)/gamr);
  double  dy=ymax-ymin;
  if(dy <= 0.0) return 0.0;

  double bwint=0.0;
  for(int i=0;i<NG;i++) {
    double y=ymin+xg[i]*dy;

    //double emd=sqrt(tan(y)*gr+ emr*emr);
    //double dc=2*emd*gr/( pow2(emd*emd-emr*emr)+pow2(emr*gamr));
    //double dc=(1+cos(2*y))*emd/gr;

    //double emd=sqrt(tan(0.5*y)*gr+ emr*emr);
    //double dc=2*emd*(1+cos(y))/gr;

    double emd=gamr/2*tan(0.5*y)+emr;
    double dc= gamr/( 0.25*gamr*gamr + pow2(emd - emr) );

    if(srt <= emd+emnuc) continue;

    double gam=deltaWidth(emd);

    //double bw= gam*0.5/(pow2(emd-emr) + 0.25*gam*gam)/M_PI;
    //double bw=2.0/M_PI*emd*emr*gam/(pow2(emd*emd-emr*emr)+pow2(emr*gam));
    double bw=2.0/M_PI*emd*emd*gam/(pow2(emd*emd-emr*emr)+pow2(emd*gam));

    double  pf=PCM(srt,emd,emnuc);
    bwint += bw*pf*dy*wg[i]/dc;
  }

  const double spin=2*4;

  // UrQMD
  //double rmat = 4e4*pow2(0.115*emr)/(pow2(srt*srt-emr*emr)+pow2(0.115*emr));

  // AMD-JAM
  // a=2*4*7000*1.15
  //gamd2=gamr*gamr
  //rmat = 8050*srt*srt*gamd2/((s-emr*emr)**2 + srt*srt*gamd2)

  //double rmat = 68.0/pow(srt-1.104,1.951);  // SMASH
  double rmat = 71.0/pow(srt-1.104,1.951);
 
  //double rmat = 0.82*68.0/pow(srt-1.104,1.8);
  //double rmat = 1.32*68.0/pow(srt-0.95,2.1);
  //double rmat = 170.588/pow(srt-0.6276,2.48674);
  //double rmat = 455/pow(srt-0.139,2.9811);

  double pr0=PCM(srt,emnuc,emnuc);
  return spin*rmat*bwint/(pr0*srt*srt)/Norm;

}

}//end namespace jam2

//g++ SigmaNDelta.cxx -I.. -I /opt/pythia8/include
//#define MAIN 1

#ifdef MAIN
#include <iostream>
#include <iomanip>
int main() {

  using namespace jam2;
  int optWid=2;
  SigmaNR sigma = SigmaNR(optWid);
  double norm=sigma.norm();
  cout << "# norm= "<< norm << endl;
  int nn=50;
  double smin=2.0;
  double  smax=6.0;
  double ds=(smax-smin)/(nn-1);
  for(int i=0;i<nn;i++) {
    double  s = smin + ds*i;
    //cout << setw(8) << s << setw(14) << sigma.deltaWidth(s) <<endl;
    cout << setw(8) << s << setw(14) << sigma.sigma(s) <<endl;
  }
  return 0;

}
#endif
