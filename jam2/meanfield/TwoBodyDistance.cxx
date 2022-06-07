#include <jam2/meanfield/TwoBodyDistance.h>

namespace jam2 {

using namespace std;

TwoBodyDistance::TwoBodyDistance(Pythia8::Settings* s) :settings(s)
{
  optP0dev=settings->mode("MeanField:optP0dev");
  optPV=settings->mode("MeanField:optMomPotential");
  optScalarDensity=settings->mode("MeanField:optScalarDensity");

  transportModel = settings->mode("MeanField:transportModel"); 
  widG = settings->parm("MeanField:gaussWidth");
  if(transportModel==1) {
    facG = 1.0/pow(4.0*M_PI*widG, 1.5);
    wG = 1.0/(4*widG);
  } else {
    facG = 1.0/pow(2.0*M_PI*widG, 1.5);
    wG = 1.0/(2*widG);
  }

  gamCM = settings->parm("Cascade:gamCM");
  gamCM2 = gamCM*gamCM;

  devgam1 = 0.0;
  devgam2 = 0.0;
}

void NonRelDistance::devV(const Vec4& Ai, const Vec4& Aj)
{
  double fsgam2i = optP0dev*dot3(Aj,v1)/p1[0];
  double fsgam2j = optP0dev*dot3(Ai,v2)/p2[0];
  devV1 = (fsgam2i*v1 + Aj/p1[0]);
  devV2 = (fsgam2j*v2 + Ai/p2[0]);
}

void NonRelDistance::density()
{
  double drsq = (r1 - r2).pT2() + pow2(gamCM*(r1[3]-r2[3]));
  density1 = gamCM*facG * exp(-drsq*wG);
  density2 = density1;
}

void NonRelDistance::psq()
{
  psq1 = (1-optPV) * (p1-p2)[0]*(p1-p2)[0] - (p1 - p2).pAbs2();
  psq2 = psq1;
}

void NonRelDistance::distanceR()
{
  Vec4 dR = r1 - r2;
  dR[3] *= gamCM2;
  dr2ijri =  2*dR;
  dr2jiri =  2*dR;
  dr2jirj = -2*dR;
  dr2ijrj = -2*dR;

  dr2ijpi = 0.0;
  dr2ijpi = 0.0;
  dr2jipj = 0.0;
  dr2ijpj = 0.0;
}

void NonRelDistance::distanceP()
{
  Vec4 dP = p1 - p2;
  psq1 = -dP.pAbs2() + (1-optPV) * dP[0]*dP[0];
  psq2 = psq1;
  //dp2ijpi = 2*( dP -optP0dev*(1-optPV)*dP[0]*p1/p1.e());
  //dp2ijpj = 2*(-dP -optP0dev*(1-optPV)*dP[0]*p2/p2.e());
  dp2ijpi = 2*( dP -optP0dev*(1-optPV)*dP[0]*v1);
  dp2ijpj = 2*(-dP -optP0dev*(1-optPV)*dP[0]*v2);

  dp2jipi = dp2ijpi;
  dp2jipj = dp2ijpj;
}


// relative distance and its derivatives between p1 and p2
// in the two body CM frame.
void TwoBodyCM::density()
{
  Vec4 dr = r1 - r2; dr[0] = 0.0;  // just in case;
  Vec4 pCM  = p1 + p2;
  double sInv = pCM.m2Calc();
  double drcmsq = dr.m2Calc() - pow2(dr * pCM)/sInv;
  //double g12 = optScalarDensity > 0 ? pCM[0]/sqrt(sInv) : 1.0;
  double den = facG * exp(drcmsq*wG);
  if(optScalarDensity==0) {
    density1 = pCM[0]/sqrt(sInv) * den;
    density2 = density1;
  } else if(optScalarDensity==1) {
    density1 = p2[0]/em2*den;
    density2 = p1[0]/em1*den;
  } else {
    density1 = den;
    density2 = density1;
  }
  /*
  cout << "gam/gamj= "<< pCM[0]/sqrt(sInv)/p1[0]*em1
    << " gam/gami= "<< pCM[0]/sqrt(sInv)/p2[0]*em2
    <<endl;
  cin.get();
  */
}

void TwoBodyCM::psq()
{
  Vec4 dP  = p1 - p2;
  Vec4 pCM = p1 + p2;
  double sInv = pCM.m2Calc();
  psq1 = dP.m2Calc() - optPV * pow2(dP * pCM)/sInv;
  psq2 = psq1;
}

void TwoBodyCM::distanceR()
{
  Vec4 pCM = p1 + p2;
  Vec4 bbi = pCM/pCM[0] - optP0dev*v1;
  Vec4 bbj = pCM/pCM[0] - optP0dev*v2;

  Vec4 dR  = r1 - r2;dR[0]=0.0;
  double sInv = pCM.m2Calc();
  double rbij = -dR * pCM/sInv;
  dr2ijri =  2*(dR + rbij*pCM);
  dr2jiri =  dr2ijri;
  dr2jirj = -dr2ijri;
  dr2ijrj = -dr2jiri;

  dr2ijpi = 2*(dR + pCM[0]*rbij*bbi)*rbij;
  dr2jipi = dr2ijpi;
  dr2jipj = 2*(dR + pCM[0]*rbij*bbj)*rbij;
  dr2ijpj = dr2jipj;
}

// relative distance and its derivatives between p1 and p2
// in the two body CM frame.
void TwoBodyCM::distanceP()
{
  Vec4 dP  = p1 - p2;
  Vec4 pCM = p1 + p2;
  double sInv = pCM.m2Calc();
  double pma = pow2(dP * pCM)/sInv;
  psq1 = dP.m2Calc() - optPV * pma;
  psq2 = psq1;
  Vec4 bbi = pCM/pCM[0] - optP0dev*v1;
  Vec4 bbj = pCM/pCM[0] - optP0dev*v2;
  dp2ijpi = 2*( dP - optP0dev*dP[0]*v1 + optPV * pCM[0]/sInv*pma*bbi);
  dp2ijpj = 2*(-dP + optP0dev*dP[0]*v2 + optPV * pCM[0]/sInv*pma*bbj);
  dp2jipi = dp2ijpi;
  dp2jipj = dp2ijpj;
}

void TwoBodyCM::devGamma()
{
  // Derivative of gamma_ij in front of Gaussian.
  Vec4 pCM = p1 + p2;
  double sInv = pCM.m2Calc();
  double fsgam2 = optP0dev*(1.0/pCM[0] - pCM[0]/sInv);
  devgam1 = pCM/sInv + fsgam2*v1;
  devgam2 = pCM/sInv + fsgam2*v2;
}

void TwoBodyCM::devV(const Vec4& Ai, const Vec4& Aj)
{
  devV1 = (optP0dev*dot3(Aj,v1)*v1 - Aj)/p1[0];
  devV2 = (optP0dev*dot3(Ai,v2)*v2 - Ai)/p2[0];
}

// relative distance squared and its derivatives between p1 and p2
// in the rest frame of p1 or p2.
void RestFrame::density() 
{
  Vec4 dR = r1 - r2; dR[0] = 0.0;  // just in case;
  double drsq1 = dR.m2Calc() - pow2(dR*v2);
  density1 = p2[0]/em2*facG * exp(drsq1*wG);

  double drsq2 = dR.m2Calc() - pow2(dR*v1);
  density2 = p1[0]/em1*facG * exp(drsq2*wG);
}

void RestFrame::psq()
{
  Vec4 dp = p1 - p2;
  double dpsq = dp.m2Calc();
  psq1 = dpsq - optPV * pow2(dp * v2);
  psq2 = dpsq - optPV * pow2(dp * v1);

}

void RestFrame::distanceR()
{
  Vec4 dr = r1 - r2; dr[0]=0.0;
  double rbi = dot3(dr,v1);
  double rbj = dot3(dr,v2);
  dr2ijri =  2*(dr + rbj*v2);  // R^2_{ij}/dr_i
  dr2jiri =  2*(dr + rbi*v1);  // R^2_{ji}/dr_i
  dr2ijrj =  -dr2ijri;
  dr2jirj =  -dr2jiri;

  dr2jipi =  2*dr*rbi/em1;     // R^2_{ji}/dp_i
  dr2ijpj =  2*dr*rbj/em2;     // R^2_{ij}/dp_j
  dr2ijpi = 0.0;
  dr2jipj = 0.0;

  //cout << "dr2ijri= "<< dr2ijri << " dr2jiri= "<< dr2jiri<<endl;
  
}

void RestFrame::distanceP()
{
  Vec4 dp = p1 - p2;
  double psq = dp.m2Calc();

  // distance squared in the rest frame of particle 2
  double m2 = p2.mCalc();
  double dot4j = dp * p2 / m2;
  psq2 = psq - optPV*dot4j*dot4j;

  // distance squared in the rest frame of particle 1
  double m1 = p1.mCalc();
  double dot4i = (dp * p1)/m1;
  psq1 = psq - optPV*dot4i*dot4i;

  // derivatives
  Vec4 bi = p1/p1.e();
  Vec4 bb = optP0dev * p2.e()*bi - p2;
  dp2ijpi = 2*(dp - optP0dev*dp[0]*bi + optPV*bb*dot4j/m2);
  dp2jipi = 2*(dp - optP0dev*dp[0]*bi - optPV*bb*dot4i/m1);

  Vec4 bj = p2/p2.e();
  Vec4 bb2 = optP0dev * p1.e()*bj - p1;
  dp2ijpj = 2*(-dp + optP0dev*dp[0]*bj + optPV*bb2*dot4j/m2);
  dp2jipj = 2*(-dp + optP0dev*dp[0]*bj - optPV*bb2*dot4i/m1);
}

void RestFrame::devV(const Vec4& A1, const Vec4& A2)
{
  //devV1 =  (optP0dev*A2.e()*p1/p1[0] - A2)/em1;
  //devV2 =  (optP0dev*A1.e()*p2/p2[0] - A1)/em2;
  devV1 =  (optP0dev*A2.e()*p1/p1[0] - A2)/p1[0];
  devV2 =  (optP0dev*A1.e()*p2/p2[0] - A1)/p2[0];
}


} // namespace jam2


