// originally taken from Nagoya fluid.
// Refs.
// D.H. Rischke et al. Nuclear Physics A 595 (1995) 346-382
// K. Okamoto and C. Nonaka, Eur. Phys. J. C (2017) 77:383.

#include <jam2/fluid/RSolverHLLE.h>

namespace jam2 {

using namespace std;

RSolverHLLE::RSolverHLLE(int n,double dx,EoS* eos0,int opttimelike,int opttau)
{
  Ndim=n;
  dX=dx;
  eos = eos0;
  facDep=0.5;
  optTauEta=opttau;
  optTimeLike=opttimelike;

  for(int i=0;i<5;i++) {
    grdQ[i].resize(n);
    FI[i].resize(n);
    QL[i].resize(n);
    QR[i].resize(n);
  }
}

RSolverHLLE::~RSolverHLLE()
{
  for(int i=0;i<5;i++) {
    grdQ[i].clear();
    FI[i].clear();
    QL[i].clear();
    QR[i].clear();
  }
}

void RSolverHLLE::solve(vector<FluidElement*>& fluid,double dt, double tau, int direction,int isource)
{
  double lam=dt/(dX*tau);

  for(int i=1;i<Ndim-1;i++) {
    grdQ[0][i]=Limiting(fluid[i]->ed(), fluid[i-1]->ed(), fluid[i+1]->ed());
    grdQ[1][i]=Limiting(fluid[i]->vx(), fluid[i-1]->vx(), fluid[i+1]->vx());
    grdQ[2][i]=Limiting(fluid[i]->vy(), fluid[i-1]->vy(), fluid[i+1]->vy());
    grdQ[3][i]=Limiting(fluid[i]->vz(), fluid[i-1]->vz(), fluid[i+1]->vz());
    grdQ[4][i]=Limiting(fluid[i]->bd(), fluid[i-1]->bd(), fluid[i+1]->bd());
  }

  QLIMIT(fluid);

  // Space average in the linear interpolation.
  for(int i=1;i<Ndim-1;i++) {
    double depend = (1.0 - fluid[i]->cs()*lam)*facDep; // Domain of dependence
    //double depend = (1.0 - abs(f->vx())*lam)*facDep; // Domain of dependence
    QL[0][i] = fluid[i]->ed() - depend * grdQ[0][i];
    QL[1][i] = fluid[i]->vx() - depend * grdQ[1][i];
    QL[2][i] = fluid[i]->vy() - depend * grdQ[2][i];
    QL[3][i] = fluid[i]->vz() - depend * grdQ[3][i];
    QL[4][i] = fluid[i]->bd() - depend * grdQ[4][i];

    QR[0][i] = fluid[i]->ed() + depend * grdQ[0][i];
    QR[1][i] = fluid[i]->vx() + depend * grdQ[1][i];
    QR[2][i] = fluid[i]->vy() + depend * grdQ[2][i];
    QR[3][i] = fluid[i]->vz() + depend * grdQ[3][i];
    QR[4][i] = fluid[i]->bd() + depend * grdQ[4][i];
  }

  HLLE(direction);

  for(int i=0;i<5;i++)
  for(int j=0;j<Ndim;j++) {
    FI[i][j] *= lam*tau;
  }

  //...Boundary condition II.
  for(int i=0;i<5;i++) {
    FI[i][0]=FI[i][1];
    FI[i][Ndim-1]=0.0;
    FI[i][Ndim-2]=FI[i][Ndim-3];
  }

  addSource(fluid,tau,dt,direction,isource);

  /*
  for(int i=1;i<Ndim-1;i++) {
    if(fluid[i]->u(0)>0.0) {
      cout << "RS " <<   i << "  before f= " <<fluid0[i]->u(0)
           << " f= " <<fluid[i]->u(0)
	   <<endl;
    }
  }
  */

}

//for(int i=1;i<Ndim-2;i++) {
//  HLLE(QL[i+1], QR[i], FI[i], lam);
//  HLLE(QR,      QL, FI[i], lam);
void RSolverHLLE::HLLE(int l )
{
// left & right
// Q = e, vx, vy, vz, d; F = numerical flux for E, Mx, My, Mz, D
// thermodynamic variables
  double UL[5], UR[5], FL[5], FR[5];
  double Ql[5], Qr[5];

  for(int i=1;i<Ndim-2;i++) {

  for(int j=0;j<5;j++) {
    Qr[j] = QL[j][i+1];
    Ql[j] = QR[j][i];
    FI[j][i]=0.0;
    UL[j]=0.0;
    FL[j]=0.0;
    UR[j]=0.0;
    FR[j]=0.0;
  }
  double c2L=0.0;
  double c2R=0.0;

  double q = pow2(Ql[1]) + pow2(Ql[2]) + pow2(Ql[3]);
  if(q < 1.0) {
    double WL = 1.0/sqrt(1.0 - q);
    double pL = eos->getPressure(Ql[4],Ql[0]);
    c2L =eos->getCs(Ql[4],Ql[0]);
    UL[0] = (pL + Ql[0])*WL*WL - pL;    // E
    UL[1] = (pL + Ql[0])*WL*WL*Ql[1];   // Mx
    UL[2] = (pL + Ql[0])*WL*WL*Ql[2];   // My
    UL[3] = (pL + Ql[0])*WL*WL*Ql[3];   // Mz
    UL[4] = WL*Ql[4];                   // D

    FL[0] = UL[l];               // flux for E
    FL[1] = UL[1]*Ql[l];         // flux for Mx       
    FL[2] = UL[2]*Ql[l];         // flux for My
    FL[3] = UL[3]*Ql[l];         // flux for Mz
    FL[4] = UL[4]*Ql[l];         // flux for D
    FL[l] += pL;
    q = UL[0]*UL[0] - UL[1]*UL[1] - UL[2]*UL[2] - UL[3]*UL[3];
    if(UL[0] <= 0.0 || q <= 0.0) {
      for(int j=0;j<5;j++) {UL[j]=0.0; FL[j]=0.0;}
      c2L=0.0;
    }
  }

  q = pow2(Qr[1]) + pow2(Qr[2]) + pow2(Qr[3]);
  if(q < 1.0) {
    double WR = 1.0/sqrt(1.0 - q);
    double pR = eos->getPressure(Qr[4],Qr[0]);
    c2R =eos->getCs(Qr[4],Qr[0]);
    UR[0] = (pR + Qr[0])*WR*WR - pR;    // E
    UR[1] = (pR + Qr[0])*WR*WR*Qr[1];   // Mx
    UR[2] = (pR + Qr[0])*WR*WR*Qr[2];   // My
    UR[3] = (pR + Qr[0])*WR*WR*Qr[3];   // Mz
    UR[4] = WR*Qr[4];                   // D

    FR[0] = UR[l];               // flux for E
    FR[1] = UR[1]*Qr[l];         // flux for Mx       
    FR[2] = UR[2]*Qr[l];         // flux for My
    FR[3] = UR[3]*Qr[l];         // flux for Mz
    FR[4] = UR[4]*Qr[l];         // flux for D
    FR[l] += pR;                 // add pressure
    q = UR[0]*UR[0] - UR[1]*UR[1] - UR[2]*UR[2] - UR[3]*UR[3];
    if(UR[0]<= 0.0 || q <= 0.0) {
      for(int j=0;j<5;j++) {UR[j]=0.0; FR[j]=0.0;}
      c2R=0.0;
    }
  }

  double eLrt = sqrt(UL[0]);
  double eRrt = sqrt(UR[0]);
  if(eLrt+eRrt ==  0.0) continue;

  // 0.5 is the suggested value
  double c2I = (eLrt*c2L*c2L + eRrt*c2R*c2R)/(eLrt + eRrt)
             + 0.5*(eLrt*eRrt)*pow2(Qr[l]-Ql[l])/pow2(eLrt + eRrt);
  c2I=sqrt(c2I);

  double vxI = (eLrt*Ql[l] + eRrt*Qr[l])/(eLrt + eRrt);

  // signal velocity
  double bL = min({0.0, (vxI - c2I)/(1.0 - vxI*c2I), 
                  (Ql[l] - c2L)/(1.0 - Ql[l]*c2L)});

  double  bR = max({0.0, (vxI + c2I)/(1.0 + vxI*c2I), 
                  (Qr[l] + c2R)/(1.0 + Qr[l]*c2R)});

  if( bR-bL == 0.0) {
      cout << i<< " bR = " << bR << " bL= "<< bL <<endl;
  }

  if(bR-bL != 0.0) {
    for(int j=0;j<5;j++)
      FI[j][i] = (bR*FL[j] - bL*FR[j] + bL*bR*(UR[j]-UL[j]))/(bR - bL);
  }

  } // end loop over fluid.

}

//======================================================================*
//           L and R state of primitive variables
//======================================================================*
void RSolverHLLE::QLIMIT(vector<FluidElement*>& fluid)
{
  for(int i=1;i<Ndim-1;i++) {
    QL[0][i] = fluid[i]->ed() - 0.5 * grdQ[0][i];
    QL[1][i] = fluid[i]->vx() - 0.5 * grdQ[1][i];
    QL[2][i] = fluid[i]->vy() - 0.5 * grdQ[2][i];
    QL[3][i] = fluid[i]->vz() - 0.5 * grdQ[3][i];
    QL[4][i] = fluid[i]->bd() - 0.5 * grdQ[4][i];

    QR[0][i] = fluid[i]->ed() + 0.5 * grdQ[0][i];
    QR[1][i] = fluid[i]->vx() + 0.5 * grdQ[1][i];
    QR[2][i] = fluid[i]->vy() + 0.5 * grdQ[2][i];
    QR[3][i] = fluid[i]->vz() + 0.5 * grdQ[3][i];
    QR[4][i] = fluid[i]->bd() + 0.5 * grdQ[4][i];
  }

  for(int i=1;i<Ndim-2;i++) {

  if(QL[4][i+1] < 0.0) {
    grdQ[4][i+1] = 0.0;
    grdQ[4][i] = 0.0;
  } else if(QR[4][i] < 0.0) {
    grdQ[4][i] = 0.0;
    grdQ[4][i+1] = 0.0;
  }
  if(QL[0][i+1] < 0.0) {
    grdQ[0][i+1] = 0.0;
    grdQ[0][i] = 0.0;
  } else if(QR[0][i] < 0.0) {
    grdQ[0][i] = 0.0;
    grdQ[0][i+1] = 0.0;
  }

  /*
  if(QL[i+1][4] < 0.0 || QR[i][4] < 0.0) {
    grdQ[i][4] = 0.0;
    grdQ[i+1][4] = 0.0;
  }
  if(QL[i+1][0] < 0.0 || QR[i][0] < 0.0) {
    grdQ[i][0] = 0.0;
    grdQ[i+1][0] = 0.0;
  }
  */
     
  double vL = sqrt(pow2(QL[1][i+1]) + pow2(QL[2][i+1]) + pow2(QL[3][i+1]));
  double vR = sqrt(pow2(QR[1][i])   + pow2(QR[2][i])   + pow2(QR[3][i]));
     
  double vmax = max({abs(fluid[i]->vx()),abs(fluid[i+1]->vx()),abs(fluid[i]->vy()),abs(fluid[i+1]->vy()),abs(fluid[i]->vz()),abs(fluid[i+1]->vz())});
  double vmin = min({abs(fluid[i]->vx()),abs(fluid[i+1]->vx()),abs(fluid[i]->vy()),abs(fluid[i+1]->vy()),abs(fluid[i]->vz()),abs(fluid[i+1]->vz())});
     
  if(vL > vmax) {
    for(int j=1;j<4;j++) {
      grdQ[j][i+1] = 0.0;
      grdQ[j][i] = 0.0;
    }
  } else if(vR > vmax) {
    for(int j=1;j<4;j++) {
      grdQ[j][i+1] = 0.0;
      grdQ[j][i] = 0.0;
    }
  }
     
  if(vL < vmin) {
    for(int j=1;j<4;j++) {
      grdQ[j][i+1] = 0.0;
      grdQ[j][i] = 0.0;
    }
  } else if(vR < vmin) {
    for(int j=1;j<4;j++) {
      grdQ[j][i+1] = 0.0;
      grdQ[j][i] = 0.0;
    }
  }

  /*
  if(vL > vmax) {
    for(int j=1;j<4;j++) {
      grdQ[i+1][j] = 0.0;
      grdQ[i][j] = 0.0;
    }
  }
  if(vR > vmax) {
    for(int j=1;j<4;j++) {
      grdQ[i+1][j] = 0.0;
      grdQ[i][j] = 0.0;
    }
  }
     
  if(vL < vmin) {
    for(int j=1;j<4;j++) {
      grdQ[i+1][j] = 0.0;
      grdQ[i][j] = 0.0;
    }
  }
  if(vR < vmin) {
    for(int j=1;j<4;j++) {
      grdQ[i+1][j] = 0.0;
      grdQ[i][j] = 0.0;
    }
  }
  */

  }

}

//***********************************************************************
//...Update conserved quantities u(:) by using flux f(:), and optionally
//...add source term in case of tau-eta coordinate.
void RSolverHLLE::addSource(vector<FluidElement*>& fluid,double tau,double dt,int iz,int is)
{
  if(optTauEta == 0 || is == 0) {

    if(optTimeLike == 1) limiterTimelike(fluid,dt/(tau*dX),iz);
    double f[5];
    for(int j=0;j<5;j++) f[j] = -FI[j][0];
    fluid[0]->updateU(f);

    for(int i=1;i<Ndim;i++) {
      for(int j=0;j<5;j++) f[j] = FI[j][i-1] - FI[j][i];
      fluid[i]->updateU(f);
    }

  // Add source term in the Milne coordinate.
  } else if(optTauEta == 1) {
      for(int i=0;i<Ndim;i++) {
	double uz=fluid[i]->u(iz);
        double vz=fluid[i]->v(iz);
	double pr=fluid[i]->pressure();
	fluid[i]->updateU(0, uz*vz/tau*dt - pr*dt);
	fluid[i]->updateU(iz,uz/tau*dt);
      }
    if(optTimeLike == 1) limiterTimelike(fluid,dt/(tau*dX),iz);
    double f[5];
    for(int j=0;j<5;j++) f[j] = -FI[j][0];
    fluid[0]->updateU(f);
    for(int i=1;i<Ndim;i++) {
      for(int j=0;j<5;j++) f[j] = FI[j][i-1] - FI[j][i];
      fluid[i]->updateU(f);
    }

  } else if(optTauEta == 2) {
//....We assume this is called only for z-direction update,and
//...assume that 3 = z-direction, 2=y-direction, 1=x-direction
    if(optTimeLike == 1) limiterTimelike(fluid,dt/(tau*dX),iz);

//...K. Murase's method. See his D-thesis for the general formulation.
    double ch=cosh(dX/2);
    double sh=sinh(dX/2);

    fluid[0]->updateU(0, -ch*FI[0][0] - sh*FI[3][0]);
    fluid[0]->updateU(3, -sh*FI[0][0] - ch*FI[3][0]);
    fluid[0]->updateU(1,-FI[1][0]);
    fluid[0]->updateU(2,-FI[2][0]);
    fluid[0]->updateU(4,-FI[4][0]);

    for(int i=1;i<Ndim;i++) {
      double f00=  ch*FI[0][i-1] - sh*FI[3][i-1];
      double f01= -sh*FI[0][i-1] + ch*FI[3][i-1];
      double f10=  ch*FI[0][i]   + sh*FI[3][i];
      double f11=  sh*FI[0][i]   + ch*FI[3][i];
      fluid[i]->updateU(0, f00-f10);
      fluid[i]->updateU(3, f01-f11);
      fluid[i]->updateU(1, FI[1][i-1]-FI[1][i]);
      fluid[i]->updateU(2, FI[2][i-1]-FI[2][i]);
      fluid[i]->updateU(4, FI[4][i-1]-FI[4][i]);
    }

  }


}

//...Flux limiter to preserve timelikeness of energy-momentum tensor
//...by K.Murase method.
//*******************************************************************
void RSolverHLLE::limiterTimelike(vector<FluidElement*>& fluid, double lam, int l)
{
  const double eps1=0.0;
  const double eps3=1e-10;

  double ch=1.0, sh=0.0;
  if(optTauEta>=1) {
    ch=cosh(dX/2);
    sh=sinh(dX/2);
  }

  for(int i=0;i<Ndim-1;i++) {
    double vx  = fluid[i]->v(l);
    double vx1 = fluid[i+1]->v(l);
    Vec4 uu  = fluid[i]->u();
    Vec4 uu1 = fluid[i+1]->u();

    //Vec4 ur, ul;
    //urul(uu,uu1,vx,vx1,lam, ur,ul);

    //...Murase new.
    Vec4 ur=min(1.0,max(0.0,0.5+vx*lam))*uu;
    Vec4 ul=min(1.0,max(0.0,0.5-vx1*lam))*uu1;

    Vec4 f4(FI[1][i],FI[2][i],FI[3][i],FI[0][i]);
    double fr0= ch*f4[0]+sh*f4[l];
    double fl0= ch*f4[0]-sh*f4[l];

    double ur0=ur[0] - fr0;
    double ul0=ul[0] + fl0;
    if(ur0 < eps1) {
      f4 = ur;
      for(int j=0;j<4;j++) FI[j][i] = f4[j];
    } else if(ul0 < eps1) {
      f4 = -ul;
      for(int j=0;j<4;j++) FI[j][i] = f4[j];
    }

    /*
    // check baryon umber
    const double eps4=1e-15;
    double b0 = fluid[i].U(4);
    double br=min(1.0,max(0.0,0.5+vx*lam))*b0;
    double bl=min(1.0,max(0.0,0.5-vx1*lam))*fluid[i+1].U(4);
    if(br - b0 < eps4 ) FI[i][4] = br;
    else if(bl + b0 < eps4 ) FI[i][4] = -bl;
    */

    //double fr3= sh*f4[0]+ch*f4[1];
    //double fl3=-sh*f4[0]+ch*f4[1];

    Vec4 fr( f4[1],f4[2],f4[3],ch*f4[0]+sh*f4[l]);
    Vec4 fl( f4[1],f4[2],f4[3],ch*f4[0]-sh*f4[l]);
    fr[l] *= ch; fr[l] += sh*f4[l];
    fl[l] *= ch; fl[l] -= sh*f4[l];

    //Vec4 fr( sh*f4[0]+ch*f4[l],f4[2],f4[3],ch*f4[0]+sh*f4[l]);
    //Vec4 fl(-sh*f4[0]+ch*f4[l],f4[2],f4[3],ch*f4[0]-sh*f4[l]);

    Vec4 ur1 = ur - fr;
    Vec4 ul1 = ul + fl;
    double ur4 = ur1.m2Calc();
    double ul4 = ul1.m2Calc();

    /*
    if(i > 0 && uu[0]> 1e-8) {
      Vec4 f3(FI[i-1][1],FI[i-1][2],FI[i-1][3],FI[i-1][0]);
      Vec4 unew =  uu + f3 - f4;
      if(unew.m2Calc() < 0.0 || unew[0]<0.0) {
      cout << i << " unew= " << unew.m2Calc() 
	  << " unew0= " << unew[0] <<endl;
	  cout << " ul4= "<< scientific <<  ul4
	       << " ul40= " << scientific << ul1[0] <<endl;
	  cout << " ur4= "<< scientific <<  ur4
	       << " ur40= " << scientific << ur1[0]
	      <<endl;
      }
    }
    */

    if(ur4 >= 0.0 && ul4 >= 0.0) continue;

    double fac=max(0.0,min(scaleF(fr,-ur),scaleF(fl,ul))-eps3);
    f4 *= fac;
    FI[0][i] = f4[0];
    FI[1][i] = f4[1];
    FI[2][i] = f4[2];
    FI[3][i] = f4[3];
    //FI[i][4] = f4[4];

    // Check.
    //Vec4 fr2( sh*f4[0]+ch*f4[1],f4[2],f4[3],ch*f4[0]+sh*f4[1]);
    //Vec4 fl2(-sh*f4[0]+ch*f4[1],f4[2],f4[3],ch*f4[0]-sh*f4[1]);
    ur1=ur-f4;
    ul1=ul+f4;
    ur4 = ur1.m2Calc();
    ul4 = ul1.m2Calc();
    if(ur4 < 0.0 || ul4 < 0.0 || ur1[0] <0.0 || ul1[0]<0.0) {
       cout << "ur<0 or ul<0?" 
	     << " ur4= "<< ur4
	     << " ul4= " << ul4
	     << " scale=  "<< fac
	     <<endl;
    }


  }

}

void RSolverHLLE::urul(Vec4 uu,Vec4 uu1,double vx,double vx1,double lam, 
	Vec4& ur,Vec4& ul)
{
 const int opt = 2;

  //....Murase original.
  if(opt == 1) {
    const double alpha=0.5;
    const double eps=0.00;
    ur=(alpha-eps)*uu;
    ul=(1.0-alpha-eps)*uu1;

  //...Murase new.
  } else if(opt == 2) {
    ur=min(1.0,max(0.0,0.5+vx*lam))*uu;
    ul=min(1.0,max(0.0,0.5-vx1*lam))*uu1;

  } else {
  //...Murase new II use the velocity of the total energy density.
    double v1=0.0;
    double v2=0.0;
    if(uu[0] != 0.0) v1=uu[1]/uu[0];
    ur = min(1.0,max(0.0,0.5+v1*lam))*uu;
    if(uu1[0] != 0.0) v2=uu1[1]/uu1[0];
    ul = min(1.0,max(0.0,0.5-v2*lam))*uu1;
  }

}

//*******************************************************************
double RSolverHLLE::scaleF(Vec4 f,Vec4 u)
{
  double a = f.m2Calc();
  double b = f * u;
  double c = u.m2Calc();

  if(a+2*b+c  >= 0.0) {
    return 1.0;
  } else if(abs(a) < 1e-10) {
    return abs(b) < 1e-10 ? 0.0 : c/(2*b);
  } else {
    b /= a;
    c /= a;
    double d4=sqrt(max(0.0,b*b-c));
    double s1=-b-d4;
    double s2=-b+d4;
    return  s2 <= 1.0 ? s2 : s1;
  }

}

} // end namespace jam2
