#include <jam2/fluid/EoS.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

namespace jam2 {

using namespace std;

EoS::EoS(const std::string& fname, bool optlog)
{
  fName = fname;
  isLogScale = optlog;
  eMin=0.0;
  eMax=20.0;
  bMax=3.0;
  maxE=1000;
  maxN=150;
  dE = (eMax-eMin)/(maxE-1);
  dN = bMax/(maxN-1);
  if(isLogScale) {
   eMin=log(0.001);
   dE = (log(eMax)-eMin)/maxE;
  }

  eosTabP = new double* [maxN];   // pressure
  eosTabT = new double* [maxN];   // Temperature
  eosTabMu = new double* [maxN];  // baryon chemical potential
  eosTabSmu = new double* [maxN]; // strangness chemical potential
  eosTabLam = new double* [maxN]; // QGP fraction
  eosTabS = new double* [maxN];   // entropy
  eosTabD = new double* [maxN];   // total particle density
  for(int im=0;im<maxN;im++) {
    eosTabP[im] = new double [maxE];
    eosTabT[im] = new double [maxE];
    eosTabMu[im] = new double [maxE];
    eosTabSmu[im] = new double [maxE];
    eosTabLam[im] = new double [maxE];
    eosTabS[im] = new double [maxE];
    eosTabD[im] = new double [maxE];
  }

  readTable();
}

EoS::~EoS()
{
  for(int im=0;im<maxN;im++) {
    delete [] eosTabP[im];
    delete [] eosTabT[im];
    delete [] eosTabMu[im];
    delete [] eosTabSmu[im];
    delete [] eosTabLam[im];
    delete [] eosTabS[im];
    delete [] eosTabD[im];
  }
  delete [] eosTabP;
  delete [] eosTabT;
  delete [] eosTabMu;
  delete [] eosTabSmu;
  delete [] eosTabLam;
  delete [] eosTabS;
  delete [] eosTabD;
}


//***************************************************************************
void EoS::readTable()
{
  ifstream in;
  in.open(fName.c_str(), ios::in);
  if(!in) {
      cerr << "Error unable to open file " << fName <<endl;
      exit(1);
  }
  string templine, templine2, templine3;
  getline(in,templine);
  getline(in,templine2);
  getline(in,templine3);

  getline(in,templine);
  istringstream is(templine);
  is >> optPot >> BagConst;

  int maxe1, maxn1;
  string ctmp1, ctmp2, ctmp3;
  istringstream is2(templine2);
  is2 >> ctmp1 >>  ctmp2 >> ctmp3 >> maxe1;

  istringstream is3(templine3);
  double b, b2;
  string ctmp4,ctmp5,ctmp6;
  is3 >> ctmp1 >>  ctmp2 >> ctmp3 >> maxn1
      >> ctmp4 >> b >> ctmp5 >> b2 >> ctmp6;

  if(maxe1 != maxE || maxn1 != maxN) {
    cerr << " maxe1 = "<< maxe1
         << " maxn1 = "<< maxn1
	 <<endl;
    exit(1);
  }

  double bden,e,p,t,bmu,smu,xqgp,s,pdens;
  for(int im=0;im<maxN;im++) {
    for(int ie=1;ie<=maxE;ie++) {

      getline(in,templine);
      istringstream is4(templine);
      is4 >> bden >> e >> p >> t >> bmu >> smu >> xqgp >> s >> pdens;

      eosTabP[im][maxE-ie]=p;
      eosTabT[im][maxE-ie]=t;
      eosTabMu[im][maxE-ie]=bmu;
      eosTabSmu[im][maxE-ie]=smu;
      eosTabLam[im][maxE-ie]=xqgp;
      eosTabS[im][maxE-ie]=s;
      eosTabD[im][maxE-ie]=pdens;
    }
      getline(in,templine);
  }

  in.close();

}


double EoS::lookup(double e,double b,double** eostab) const
{
  int i=min(max(0,int(b/dN)),maxN-2);
  int j=min(max(0,int((e-eMin)/dE)),maxE-2);

  double x2=(b-(i*dN))/dN;
  double y2=(e-(eMin+j*dE))/dE;
  double x1=1.0-x2;
  double y1=1.0-y2;

  return x1*y1*eostab[i][j]
             +x2*y1*eostab[i+1][j]
             +x1*y2*eostab[i][j+1]
             +x2*y2*eostab[i+1][j+1];
}

double EoS::lookuplog(double e,double b,double** eostab) const
{
  int i=min(max(0,int(b/dN)),maxN-2);
  int j=min(max(0,int((log(e)-eMin)/dE)),maxE-2);

  double x2=(b-(i*dN))/dN;
  double y2=(log(e)-(eMin+j*dE))/dE;
  double x1=1.0-x2;
  double y1=1.0-y2;

  return x1*y1*eostab[i][j]
             +x2*y1*eostab[i+1][j]
             +x1*y2*eostab[i][j+1]
             +x2*y2*eostab[i+1][j+1];
}

//**************************************************************************
double EoS::getEoSTable(double b,double e,double** eostab,double mass_dimension) const
{
  b=abs(b);

  if(!isLogScale) return max(0.0,lookup(e,b,eostab));
  else return max(0.0,lookuplog(e,b,eostab));

  //if(isLogScale) return lookuplog(e,b,eostab);
  //else return lookup(e,b,eostab);

  if(isLogScale) e = log(e);
  int ie=int((e-eMin)/dE);


  if(ie >= maxE) ie=maxE-1;
  if(ie < 0) ie=0;
  int ib=int(b/dN);
  if(ib < 0) ib=0;
  else if(ib >= maxN) ib=maxN-1;

  if(isLogScale && (ib == 0 && ie == 0)) return eostab[0][0];

  if(ib == 0 && ie == 0) {
    double const s1 = std::pow(b/dN, mass_dimension / 3.0);
    double const t1 = std::pow(e/dE, mass_dimension / 4.0);
    double const s0 = 1.0 - s1;
    double const t0 = 1.0 - t1;
    double interpolated = 0.0;
    interpolated += s0 * t0 * eostab[0][0];
    interpolated += s0 * t1 * eostab[0][1];
    interpolated += s1 * t0 * eostab[1][0];
    interpolated += s1 * t1 * eostab[1][1];
    return interpolated;;
  }

  if(ie == maxE-1) {
    if(ib == maxN-1) {
      return eostab[ib][ie];
    } else {
      double val1=eostab[ib][ie];
      double val2=eostab[ib+1][ie];
      return max(0.0,rlinter(ib*dN,(ib+1)*dN,val1,val2,b));
    }
  }

  // first interpolation in e
  double y1=eostab[ib][ie];
  double y2=eostab[ib][ie+1];
  double x1=eMin+ie*dE;
  double x2=eMin+(ie+1)*dE;
  double val1=rlinter(x1,x2,y1,y2,e);
  if(ib == maxN-1) return max(0.0,val1);

  y1=eostab[ib+1][ie];
  y2=eostab[ib+1][ie+1];
  double val2=rlinter(x1,x2,y1,y2,e);

  /*
 double p1=max(0.0,lookup(e,b,eostab));
 double p2=max(0.0,rlinter(ib*dN,(ib+1)*dN,val1,val2,b));
 if(abs(p1-p2)>1e-9) {
 cout<< setprecision(16)<< " ie= "<< ie << " ib= "<< ib << " e= "<< e << " b= "<< b  << "p= "<<  p1 << " p2= "<< p2 <<endl;
 cin.get();
 }
 */
 
  return max(0.0,rlinter(ib*dN,(ib+1)*dN,val1,val2,b));
}

} // end namespace jam2

//#define MAIN 1
#ifdef MAIN
using namespace std;
int main() {
  //string fname = "eosB235JAMsoft.dat";
  string fname = "eosB235JAMsoft1mev.dat";
  jam2::EoS *eos = new jam2::EoS(fname,1);
  int nn=1000;
  double de=0.01;
  for(int i=1;i<nn;i++) {
    double e=i*de;
    double b = 0.0;
    double t = eos->getT(b,e/HBARC)*HBARC;
    cout << setw(13) << e
	 << setw(13) << t
	 <<endl;
  }

  return 0;
}
#endif

