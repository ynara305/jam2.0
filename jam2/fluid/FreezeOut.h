#ifndef jam2_fluid_FreezeOut_H
#define jam2_fluid_FreezeOut_H

#include <vector>
#include <jam2/fluid/cornelius.h>
#include <jam2/fluid/EoS.h>
#include <jam2/fluid/FluidElement.h>
#include <Pythia8/Settings.h>

namespace jam2 {

class FreezeOut
{
private:
  Pythia8::Settings* settings;
  Cornelius cornelius;
  EoS *eos;
  int FT,FX,FY,FZ;
  const int N=2;
  double dfx[4];
  int    maxX, maxY, maxZ;
  double dT,dX,dY,dZ,dV;
  double eFreezeOut, tFreezeOutCut;
  int    nFreezeout;
  double ****HyperCube, ****dsite, ****vxsite, ****vysite, ****vzsite;
  //double ****fE;   // energy density
  //double ****fD;   // baryon density
  //double ****fVx;  // v_x
  //double ****fVy;  // v_y
  //double ****fVz;  // v_z
  std::vector<double> freezeOutTime;     // freezeout time
  std::vector<int>    freezeOutR[3];     // freezeout position
  std::vector<double> freezeOutV[4];     // flow velocity
  std::vector<double> freezeOutDsgma[4]; // surface
  std::vector<double> freezeOutT;        // Temperature
  std::vector<double> freezeOutMuB;      // mu_B
  std::vector<double> freezeOutMuS;      // mu_S
  std::vector<double> freezeOutBden;     // baryon density
  std::vector<double> freezeOutNum;      // particle density
public:
  FreezeOut(Pythia8::Settings* s,EoS* eos0);
  ~FreezeOut();
  void clear();
  void reset(std::vector<FluidElement*>& f);
  void isoenergyFreezeout(std::vector<FluidElement*>& f,int it,double htime);
  int  isochronousFreezeout(std::vector<FluidElement*>& f, double htime,int opt);
  void freezeout(std::vector<FluidElement*>& fluid,double htime);
  double interpolation(double ****qsite, double rt, double rx,
	double ry, double rz);
  double interpolation4(double ****qsite, double rt, double rx,
	double ry, double rz);

  double getFreezeOutEnergyDensity() {return eFreezeOut;}
  int site(int x, int y, int z) {return x + maxX*(y+maxY*z);}
  int    r(int i, int j) { return freezeOutR[i-1][j];}
  double t(int i)        { return freezeOutTime[i];}
  double v(int i, int j) { return freezeOutV[i][j];}
  double dsigma(int i, int j) { return freezeOutDsgma[i][j];}
  double temperature(int i)  { return freezeOutT[i];}
  double mub(int i)      { return freezeOutMuB[i];}
  double mus(int i)      { return freezeOutMuS[i];}
  double num(int i)      { return freezeOutNum[i];}
  double bd(int i)       { return freezeOutBden[i];}
  int    getN()     {return freezeOutT.size();}
  void   temperature(int i, double a ){ freezeOutT[i]=a;}
  double getTCut() {return tFreezeOutCut;}
};

}
#endif
