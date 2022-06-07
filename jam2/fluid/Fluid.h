#ifndef jam2_fluid_Fluid_h
#define jam2_fluid_Fluid_h

#include <Pythia8/Basics.h>
#include <jam2/fluid/RSolverHLLE.h>
#include <jam2/fluid/FluidElement.h>
#include <jam2/fluid/FreezeOut.h>
#include <Pythia8/Settings.h>
#include <jam2/collision/Collision.h>
#include <jam2/JamAna.h>

namespace jam2 {

class Fluid
{
private:
  Pythia8::Settings* settings;
  int maxX,maxY,maxZ, volume;
  int origX,origY,origZ;
  double dT,dX,dY,dZ, vol,volH;
  double passTime, eFreezeOut;
  double tau0;
  int    nConversion;
  FreezeOut *freezeout;
  EoS* eos;
  std::vector<FluidElement*> fluid;
  std::vector<FluidElement*> fluidx, fluidy, fluidz;
  RSolverHLLE *solverx, *solvery, *solverz;
  Pythia8::Vec4 uzero,uout,pTot;
  Pythia8::Vec4 pTotal;
  double bzero,bout;
  double baryonTot,chargeTot,strangeTot, entroTot;
  int optTauEta,jobMode,optTimeLike,optRescale,optFreezeOut,secondOrder;
  const int N=2;
  int optVectorPotential,optVdot,optPropagate;
  ParticleDensity* particleDens;
  int fluidTimeStep;

public:
  Fluid(Pythia8::Settings* s, EoS* eos0, FreezeOut* f);
  //  int nx=100,int ny=100, int nz=100,
  //  double dt=0.15, double dx=0.3,double dy=0.3, double dz=0.3,
  //  int optt=1, int optr=0, int opttau=0);

  ~Fluid();
  int  evolution(int step,double gtime, Collision* event);
  void reset();
  void evolution3d(int it,double tau);
  void evolution3d2(int it,double tau);
  int  checkout(double t);
  void update_x(double tau,double dt,int is);
  void update_y(double tau,double dt,int is);
  void update_z(double tau,double dt,int is);
  void update_local_val(std::vector<FluidElement*>& fluid,double tau,int n);
  void setLocalQuantity(int ix,int iy,int iz,double tau=1.0);
  void thermal(Vec4& uu,double u4,double& elc,double& nblc,double& plc,
	double& cslc, double& v);
  int iteratev(double u0,double u4,double m,
	double& elc,double& nblc,double& plc,double& vf);
  void hydroTotalEnergy(const int it,const double tau,
	double& etot,double& pxtot,double& pytot,double& pztot,double& btot);
  void fluidFraction(const double tau,Collision* event,std::ostream& out,int opt=0);

  void checkboundary_x();
  void checkboundary_y();
  void checkboundary_z();
  void outputhzx(int it);
  void setPassTime(double t) {passTime=t;}

  std::vector<FluidElement*>& getElement() {return fluid;}
  int    xDim()    {return maxX;}
  int    yDim()    {return maxY;}
  int    zDim()    {return maxZ;}
  double dx()      {return dX;}
  double dy()      {return dY;}
  double dz()      {return dZ;}
  double getVol()  {return vol;}
  double getVolH() {return volH;}
  double getX(int ix, double x) {return dX*(ix - origX - 0.5 + x);}
  double getY(int iy, double y) {return dY*(iy - origY - 0.5 + y);}
  double getZ(int iz, double z) {return dZ*(iz - origZ - 0.5 + z);}
  int getOptTauEta() const {return optTauEta;}
  Pythia8::Vec4 pFluidTotal() {return pTotal;}
  void updatePTotal(Pythia8::Vec4 p) {pTotal += p;}
  void updateBTotal(double b) {baryonTot += b;}
  void updateCTotal(double c) {chargeTot += c;}
  void updateSTotal(double s) {strangeTot += s;}
  void updateNC() {nConversion++;}
  void resetNC() {nConversion=0;}
  void setNC(int i)  {nConversion=i;}
  int  isFluid() {return nConversion;}
  void removeTotalU(Pythia8::Vec4 p) {pTotal -= p;}
  void removeTotalCharge(double nch) {chargeTot -= nch;}
  void removeTotalStrangeness(double s) {strangeTot -= s;}
  void removeTotalBaryon(double b) {baryonTot -= b;}

  void setU(int i, int j, int k, Vec4 u) {site(i,j,k)->setU(u);}
  void setB(int i, int j, int k, double b) {site(i,j,k)->setB(b);}
  void addU(int i, int j, int k, Vec4 u) {site(i,j,k)->addU(u);}
  void addB(int i, int j, int k, double b) {site(i,j,k)->addB(b);}

  void seInitialTime(double t) {tau0=t;}
  //int site(int x, int y, int z) {return x + maxX*(y+maxY*z);}
  //FluidElement* site(int x, int y, int z) {return fluid + x + maxX*(y+maxY*z);}
  FluidElement* site(int x, int y, int z) {return fluid[x + maxX*(y+maxY*z)];}
  Pythia8::Vec4 u(int x,int y,int z) {return site(x,y,z)->u();}
  double U(int x,int y,int z,int i) {return site(x,y,z)->U(i);}
  double energyDensity(int x, int y, int z) {return site(x,y,z)->ed();}
  double baryonDensity(int x, int y, int z) {return site(x,y,z)->bd();}
  void print(int i, int j, int k) {site(i,j,k)->print();}

  void resetFreezeOut() {freezeout->reset(fluid);}
  void reset_fluid_element(int x,int y,int z) { site(x,y,z)->reset();}
  bool inside(double x, double y, double z,int& ix, int& iy, int& iz) {
   ix = int(x/dX + 0.5 +origX);
   if(ix <0 || ix >= maxX ) return false;
   iy = int(y/dY + 0.5 +origY);
   if(iy <0 || iy >= maxY ) return false;
   iz = int(z/dZ + 0.5 +origZ);
   if(iz <0 || iz >= maxZ ) return false;
   return true;
  }
  bool inside(Vec4 r,int& x, int& y, int& z) {
   x = int(r[1]/dX + 0.5 +origX);
   if(x <0 || x >= maxX ) return false;
   y = int(r[2]/dY + 0.5 +origY);
   if(y <0 || y >= maxY ) return false;
   z = int(r[3]/dZ + 0.5 +origZ);
   if(z <0 || z >= maxZ ) return false;
   return true;
  }
  double gauss(Vec4& x, Vec4& u, double w, int i, int j, int k,int opt=1) {
    if(opt==1) {
      Vec4 y(dX*(i-origX), dY*(j-origY), dZ*(k-origZ),x[0]);
      double xtra = cmDistanceSquare(x-y,u);
      return exp( xtra/w )*vol;
    } else {
      return  gaussInter(x,u,w,i,j,k,opt);
    }
  }
  double gaussInter(Vec4& x, Vec4& u, double w, int i, int j, int k,int opt=1);

  double makeTwoNuclei(double ecm,double A, double B, double b, int opt);
  double WS(double x,double y,double z,double rad,double rho0);
  double hard_sphere(double x,double y,double z,double rad,double rho0);

};
}
#endif
