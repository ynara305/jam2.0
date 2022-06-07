#ifndef jam2_fluid_RSolverHLLE_H
#define jam2_fluid_RSolverHLLE_H

#include <algorithm>
#include <vector>
#include <cmath>
#include <jam2/hadrons/JamStdlib.h>
#include <jam2/fluid/EoS.h>
#include <jam2/fluid/FluidElement.h>

namespace jam2 {


class RSolverHLLE
{
private:
//  FluidElement* fluid;
  double dX;
  EoS *eos;
  int Ndim;
  std::vector<double> grdQ[5], FI[5], QL[5], QR[5];
  const double Alpha = 2.0;
  int optTauEta, optTimeLike;
protected:
  double facDep;
public:
  RSolverHLLE(int n, double dx, EoS* eos0,int opttimelike=1,int opttau=0);
  ~RSolverHLLE();
   void solve(std::vector<FluidElement*>& fluid,double dt, double tau,int iz, int is);

  // Van Leer's Limiter
  double VLlimiter(double a,double b,double c) {
    //return SIGN(1.0, c)*(SIGN(0.50,a*b) + 0.50)
    return std::copysign(1.0, c)*(std::copysign(0.50,a*b) + 0.50)
         *std::min({Alpha*std::abs(a), Alpha*std::abs(b), std::abs(c)});
  }
  double SIGN(double a, double b) {return b>=0.0 ? std::abs(a) : -std::abs(a);}

  // u0 = u[i], u1=u[i-1], u2=u[i+1]
  double Limiting(double u0, double u1, double u2) {
    double dul=u0 - u1;
    double dur=u2 - u0;
    double dum=0.5*(u2 - u1 );
    return VLlimiter(dul,dur,dum);
  }
  //void HLLE(double* Qr, double* Ql, double* F, double lam);
  void HLLE(int direction);
  void QLIMIT(std::vector<FluidElement*>& f);
  void addSource(std::vector<FluidElement*>& fluid,double tau,double dt,int iz,int is);
  void limiterTimelike(std::vector<FluidElement*>& fluid, double lam,int iz);
  void urul(Vec4 uu,Vec4 uu1,double vx,double vx1,double lam, Vec4& ur,Vec4& ul);
  double scaleF(Vec4 f,Vec4 u);

};

} // end namespace jam2
#endif // jam2_fluid_RSolverHLLE_H

