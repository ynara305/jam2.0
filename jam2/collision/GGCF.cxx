#include <iostream>
#include <iomanip>
#include <cmath>
#include <jam2/hadrons/GaussPoints.h>
#include <jam2/collision/GGCF.h>

namespace jam2 {

using namespace std;

bool GGCF::computeParameter(double sig)
{
  int N2=38;  // # of integration points
  //double sigma=sigmaTot*0.1;   // sigma_in(s) [mb --> fm^2]

  //sigma0=72.6;  // starting point of iteration
  sigma0=sig;  // starting point of iteration
  double sigma00;     // holds value from previous iteration step
  int ntry=0;
  do {                                       // iterate ...
    sigma00 = sigma0;
    double f=0.0, df=0.0;
    for(int ib=0; ib<N2; ib++) {  // integral d(sigma) from 0 to infinity
      double b = xg38f[ib]*sigMax;
      double db = wg38f[ib]*sigMax;
      double P = Prob(b,sigma0,Omega);
      double dP = P*(b*2*(b/sigma0-1.0)/(Omega*Omega*sigma0*sigma0) - 1.0/(b+sigma0));
      f  += db*P*(b - sigmaTot);
      df += db*dP*(b - sigmaTot);
    }
    sigma0 -= f/df;
    //cout << "# sigma0=" << sigma0 << endl;
    if(ntry++ > 100) {
      cout << " GGCF::computeParameter infinity loop? "<< ntry << " sigmatot= "<< sigmaTot
	<< " omega_s= "<< omegaSigma<<endl;
      return false;
    }
  } while(abs(sigma0-sigma00) > 1e-6);       // ... until sigma0 has converged

  return true;
}

// compute sigma0 and Norm from the formula 
// int_0^infty dsigma P = 1, int_0^\infty dsigma sigma P = sigmaTot
bool GGCF::computeParameter2(double sig, double omeg)
{
  int N2=38;  // # of integration points
  //double sigma=sigmaTot*0.1;   // sigma_in(s) [mb --> fm^2]

  sigma0=sig;  // starting point of iteration
  Omega = omeg;
  double sigma00, Omega0;     // holds value from previous iteration step
  int ntry=0;
  //double h = 2e-8;
  do {                                       // iterate ...
    sigma00 = sigma0;
    Omega0 = Omega;
    double f1=0.0, f2=0.0;
    double df11=0.0, df12=0.0, df21=0.0, df22=0.0;
    for(int ib=0; ib<N2; ib++) {  // integral d^2b from 0 to Bmax
      double b = xg38f[ib]*sigMax;
      double db = wg38f[ib]*sigMax;
      double P = Prob(b,sigma0,Omega);
      // derivative with respect to sigma0
      double dP1 = P*(b*2*(b/sigma0-1.0)/(Omega*Omega*sigma0*sigma0) - 1.0/(b+sigma0));
      //double dP1 = (Prob(b,sigma0+h,Omega)-P)/h;
      // derivative with respect to Omega;
      double dP2 = P*2*pow2(b/sigma0-1.0)/Pythia8::pow3(Omega);
      //double dP2 = (Prob(b,sigma0,Omega+h)-P)/h;
      double sig = b - sigmaTot;
      double ome = pow2(b/sigmaTot) -(1.0+omegaSigma);
      f1   += db*P*sig;
      f2   += db*P*ome;
      df11 += db*dP1*sig;
      df12 += db*dP2*sig;
      df21 += db*dP1*ome;
      df22 += db*dP2*ome;
    }
    double det = df11*df22 - df12*df21;
    sigma0 -= ( f1*df22-df12*f2)/det;
    Omega  -= (-f1*df21+df11*f2)/det;
    //cout << "# sigma0=" << sigma0 << " Omega= " << Omega << endl;
    if(ntry++ > 100) {
      cout << " GGCF::computeParameter2 infinity loop? "<< ntry << " sigmatot= "<< sigmaTot
	<< " omega_s= "<< omegaSigma<<endl;
      return false;
    }
  } while(abs(sigma0-sigma00) > 1e-4 || abs(Omega-Omega0) > 1e-4);       // ... until sigma0 has converged

  return true;

}

double GGCF::sum2()
{
  int N=38;
  //int N=100;
  double f=0.0;
  //double xmin=0.0, a = 10.0;
  for(int i=0; i<N; i++) {
     f += wg38[i]*Prob(xg38[i],sigma0,Omega);
     //f += wg48[i]*Prob(xg48[i]);
    //double s  = (xmin + a*xg100f[i])/(1-xg100f[i]);
    //double ds = (xmin + a)/pow2(1-xg100f[i])*wg100f[i];
    //f += ds*Prob(s);
  }
  Norm = 1.0/f;
  return f;
}

double GGCF::sum()
{
  int N=38;
  double f=0.0;
  for(int ib=0; ib<N; ib++) {
     f += sigMax*wg38f[ib]*Prob(sigMax*xg38f[ib],sigma0,Omega);
    //double b = xg38[ib];
    //  f += wg38[ib]*2.0/(b*b+1.0);
  }
  Norm = 1.0/f;
  return f;
}

double GGCF::sigma()
{
  int N=38;
  double f=0.0;
  for(int ib=0; ib<N; ib++) {
    double b = xg38[ib]*sigMax;
      f += sigMax*wg38[ib]*b*Prob(b,sigma0,Omega);
  }
  return f;
}

double GGCF::omega()
{
  int N=38;
  double f=0.0;
  for(int ib=0; ib<N; ib++) {
    double b = xg38[ib]*sigMax;
      f += sigMax*wg38[ib]*b*b*Prob(b,sigma0,Omega);
  }
  return Norm*f/(sigmaTot*sigmaTot)-1.0;
}

} // end namespace jam2


//g++ SigmaRR.cxx -I.. -I /opt/pythia8/include
//g++ SigmaRR.cxx -I.. -I/export/ynara/lib/pythia8/include

//#define MAIN 1

#ifdef MAIN
int main(int argc, const char *argv[])
{
  jam2::GGCF *ggcf = new jam2::GGCF();

  double norm=ggcf->sum();
  cout << " sum= "<< norm <<endl;
  cout << " sum2= "<< ggcf->sum2()<<endl;
  cout << " sig= "<< ggcf->sigma()/norm <<endl;
  cout << " omega= "<< ggcf->omega() <<endl;
  double xm=ggcf->peak();
  cout << " peak= "<< xm <<endl;
  cout << " Pmax= "<< ggcf->Prob(xm,ggcf->sig0(),ggcf->omeg())<<endl;
  cin.get();

  /*
  ggcf->computeParameter();
  norm=ggcf->sum();
  cout << " sum= "<< norm <<endl;
  cout << " sig= "<< ggcf->sigma()/norm <<endl;
  cout << " omega= "<< ggcf->omega() <<endl;
  */

  double smin=10, smax=7000;
  int nn= 100;
  double omega=0.1;
  double Omeg=1.0;
  double sig0 = 40;
  double ds=(smax-smin)/(nn-1);
  for(int i=0;i<nn;i++) {
    double srt=smin +i*ds;
    double sigtot= 21.7*pow(srt*srt,0.0808) + 56.08*pow(srt*srt,-0.4525);
    ggcf->setParam(sigtot,omega);
    ggcf->computeParameter2(sig0,Omeg);
    double norm=ggcf->sum();
    sig0 = ggcf->sig0();
    Omeg = ggcf->omeg();
    cout << setw(7) << srt
         << setw(15) << sigtot
         << setw(15) << ggcf->sigma()/norm
         << setw(15) << omega
         << setw(15) << ggcf->omega()
         << setw(15) << sig0
         << setw(15) << Omeg
	 <<endl;
  }

  delete ggcf;

  return 0;
}
#endif
