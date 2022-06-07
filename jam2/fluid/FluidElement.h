#ifndef jam2_fluid_FLuidElement_h
#define jam2_fluid_FLuidElement_h
//======================================================================
//  U(0) = g^2*(e+p)-p
//  U(1) = g^2*(e+p)*v_x
//  U(2) = g^2*(e+p)*v_y
//  U(3) = g^2*(e+p)*v_z 
//  U(4) = g*n_B
//*********************************************************************

#include <Pythia8/Basics.h>

namespace jam2 {

using Pythia8::Vec4;

class QVal
{
public:
  QVal(double q0=0.0, double q1=0.0, double q2=0.0, double q3=0.0,
	  double q4=0.0) {
      Q[0]=q0; Q[1]=q1; Q[2]=q2; Q[3]=q3; Q[4]=q4;
  }

  double Q[5];
  double& operator()(int i) {return Q[i];}
  QVal& operator-(const QVal& q) {for(int i=0;i<5;i++) Q[i] -= q.Q[i];
      return *this;}
  QVal operator-() const {QVal tmp; for(int i=0;i<5;i++) tmp.Q[i]=-Q[i];
     return tmp;}
};

class FluidElement
{
private:
  // conserved quantities. 
  Vec4   ul; // T^{0i}
  double ql; //  baryon density.
  //int optRescale;

  // Local thermodynamic variables
  double bl[2], el[2], pl, csl;
  double vl[2][4];

public:
  FluidElement(Vec4 u0=0.0, double q=0.0,
	  double b=0.0, double e=0.0, double p=0.0,
	  double vx0=0.0, double vy0=0.0,double vz0=0.0,
	  double cs0=0.0) : ul(u0), ql(q), pl(p), csl(cs0) {
	    bl[1]=b; el[1]=e;
	    bl[0]=b; el[0]=e;
	    vl[1][1]=vx0; vl[1][2]=vy0; vl[1][3]=vz0;
	    vl[0][1]=vx0; vl[0][2]=vy0; vl[0][3]=vz0;
	    vl[1][0]=sqrt(std::max(0.0,1.0-vsq())); 
	    vl[0][0]=vl[1][0];
	  }

  Vec4 u() {return ul;}
  double u(int i) {return ul[i];}
  double b() {return ql;}
  double U(int i) {if(i<4) return ul[i]; else return ql;}

  double v(int i) {return vl[1][i];}
  double vx()  {return vl[1][1];}
  double vy()  {return vl[1][2];}
  double vz()  {return vl[1][3];}
  double ed()  {return el[1];}
  double bd()  {return bl[1];}
  double vx(const int i)  {return vl[i][1];}
  double vy(const int i)  {return vl[i][2];}
  double vz(const int i)  {return vl[i][3];}
  double ed(const int i)  {return el[i];}
  double bd(const int i)  {return bl[i];}

  double pr()  {return pl;}
  double vsq() {return vl[1][1]*vl[1][1] + vl[1][2]*vl[1][2] + vl[1][3]*vl[1][3];}
  double pressure()  {return pl;}
  double cs()  {return csl;}

  void U(int i,double a) {if(i<4) ul[i]=a; else ql=a;}
  void ed(double e) {el[1]=e;}
  void bd(double e) {bl[1]=e;}
  void pr(double p) {pl=p;}
  void vx(double v) {vl[1][1]=v;}
  void vy(double v) {vl[1][2]=v;}
  void vz(double v) {vl[1][3]=v;}
  void cs(double c) {csl=c;}

  void  setU(Vec4 u)   {ul = u;}
  void  setB(double q) {ql = q;}
  void  addU(Vec4 u)   {ul += u;}
  void  addB(double q) {ql += q;}
  void  updateU(double* a) {
      for(int i=0;i<4;i++) {ul[i] += a[i];} ql+=a[4];
  }
  void  updateU(int i,double a) { if(i<4) ul[i] += a; else ql += a;}

  void  updateU(QVal a) {
      for(int i=0;i<4;i++) {ul[i] += a(i);} ql += a(4);}

  //void swapxy() {std::swap(vxl,vyl);std::swap(ul[1],ul[2]);}
  void resetU() { ul=0.0; ql=0.0;}

  void savePreviousLocalVal() {
    el[0]=el[1]; bl[0]=bl[1];
    vl[0][0]=vl[1][0];
    vl[0][1]=vl[1][1];
    vl[0][1]=vl[1][1];
    vl[0][3]=vl[1][3];
  }
    

  void reset() {
    ul=0.0; ql=0.0;
    bl[1]=0.0;el[1]=0.0;pl=0.0; vl[1][1]=vl[1][2]=vl[1][3]=0.0; vl[1][0]=1.0; csl=0.0; 
  }
  Vec4 scaleU(double& bzero, int optRescale=1); 
  void print();

};
}
#endif
