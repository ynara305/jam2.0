#ifndef jam2_fluid_EoS_h
#define jam2_fluid_EoS_h

#include <iostream>
#include <string>
#include <cmath>

#ifndef HBARC
#define HBARC 0.19732698
#endif

namespace jam2 {

class EoS
{
private:
  double eMin, eMax,bMax, dE,dN;
  bool   isLogScale;
  int    maxE, maxN;
  double **eosTabP, **eosTabT, **eosTabMu, **eosTabSmu, **eosTabLam;
  double **eosTabS, **eosTabD;
  std::string fName;
  int    optPot;
  double BagConst;
public:
   EoS(const std::string& fname, bool logscale=true);
  ~EoS();

private:
  void readTable();
  double lookup(double e,double b,double** eostab) const;
  double lookuplog(double e,double b,double** eostab) const;
  static double rlinter(double x1,double x2,double y1,double y2,double x) {
      return (y2-y1)*(x-x1)/(x2-x1)+y1;
  }
  double getEoSTable(double b, double e, double** eostab, double mass_dimension) const;

public:
  double getPressure(double b, double e) const {
    if(e*HBARC >= eMax) return (e - 4*BagConst)/3.0;
    return getEoSTable(b, e*HBARC, eosTabP, 4.0)/HBARC;
    //return getPressureInGeV(b,e*HBARC)/HBARC; 
  }
  // e [GeV/fm^4] b [1/fm^3]
  double getPressureInGeV(double b,double e) const {
    if(e >= eMax) return (e - 4*BagConst)/3.0;
    return getEoSTable(b, e, eosTabP, 4.0);
  }

  double getCs(double b, double e) const {return 1.0/sqrt(3.0);}

  double getT(double b,double e) const {
      return getEoSTable(b, e*HBARC, eosTabT, 1.0)/HBARC;}

  double getMuB(double b,double e) const {
      return getEoSTable(b, e*HBARC, eosTabMu, 1.0)/HBARC;}

  double getMuS(double b,double e) const {
      return getEoSTable(b, e*HBARC, eosTabSmu, 1.0)/HBARC;}

  double getEntropy(double b,double e) const {
      return getEoSTable(b, e*HBARC, eosTabS, 3.0);}

  double getMult(double b,double e) const {
      return getEoSTable(b,e*HBARC,eosTabD, 3.0);}

};
} // end namespace jam2
#endif // jam2_fluid_EoS_h
