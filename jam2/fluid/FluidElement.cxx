#include <jam2/fluid/FluidElement.h>

namespace jam2 {

using namespace std;

Vec4 FluidElement::scaleU(double& bzero, int optRescale) 
{
  const double eps=1e-12;
  double u2=ul.pAbs();
  bzero=0.0;
  if(ul[0] >= u2) return 0;

  Vec4 u0=ul;
  if(u2 < 1e-16) {
    ul=0.0;
    bzero=ql;
    return u0;
  }

  // Scale energy.
  if(optRescale == 1) {
    ul[0] = u2 + eps;

  // Scale momentum.
  } else if(optRescale == 2) {

    ul[1]=ul[0]*(1.0-eps)*ul[1]/u2;
    ul[2]=ul[0]*(1.0-eps)*ul[2]/u2;
    ul[3]=ul[0]*(1.0-eps)*ul[3]/u2;

  // Scale transverse momentum.
  } else {
    double ut2=ul[0]*ul[0]-ul[3]*ul[3];
    if(ut2 > 0.0) {
      double ut=ul.pT();
      ul[1]=sqrt(ut2)*(1.0-eps)*ul[1]/ut;
      ul[2]=sqrt(ut2)*(1.0-eps)*ul[2]/ut;
    } else {
      ul[0] = u2 + 1e-12;
    }
  }

  return u0 - ul;
}

void FluidElement::print()
{
  cout << " U= " << ul*HBARC <<endl;
  cout << " B= " << ql <<endl;
  cout << " bl= " << bl[1]
       << " el= "<< el[1]*HBARC
       << " pl= "<< pl *HBARC<< endl
       << " vx= "<< vl[1][1]
       << " vy= " <<vl[1][2]
       << " vz= "<< vl[1][3]
       << " cs= "<<csl
       << endl;
}

} // end namespace jam2
