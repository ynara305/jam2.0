#include <jam2/meanfield/RQMDs0.h>
#include <jam2/meanfield/NonRelPotential.h>

// This is the same implementation as RQMD/S in 2005
// Relativistic quantum molecular dynamics with Skyrme potential
// potentials are implemented as scalar type.

namespace jam2 {

using namespace std;

RQMDs0::RQMDs0(Pythia8::Settings* s) : MeanField(s)
{
  //double alpha,beta,g,C1,C2,mu1,mu2;
  //MeanField::initParam(alpha,beta,g,C1,C2,mu1,mu2);

  potential = new NonRelPotential(settings);
  t1  = potential->getT1();
  t3  = potential->getT3();
  //t5  = potential->getT5();
  gam = potential->getGam();
  t3f = gam*t3;

  pmu1 = potential->getPmu1(); 
  pmu2 = potential->getPmu2();
  vex1 = potential->getVex1();
  vex2 = potential->getVex2();

  //double  rho0=settings->parm("MeanField:rho0");
  //gam = g;
  //t1 = 0.5*alpha/rho0;
  //t3 = beta/(gam+1.0)/pow(rho0,gam);
  //t3f = gam*t3;

  //pmu1= mu1*mu1;
  //pmu2= mu2*mu2;
  //vex1=C1/(2*rho0);
  //vex2=C2/(2*rho0);

  double widG = settings->parm("MeanField:gaussWidth");
  facG = 1.0/pow(4.0*M_PI*widG, 1.5);
  wG = 1.0/(4*widG);
}

void RQMDs0::evolution(list<EventParticle*>& plist,double gtime, double dt,int step)
{
  //double ctime = step * dt;
  globalTime = gtime;
  part.clear();
  for(auto& i : plist) {
    if(i->isMeanField(gtime,optPotential)) part.push_back(i);
  }

  NV = part.size();
  force = new Vec4 [NV];
  forcer = new Vec4 [NV];
  rho = new double [NV];
  rhog = new double [NV];
  vmoms = new double [NV];
  rhom = new double *[NV];
  for(int i=0;i<NV;i++) rhom[i] = new double [NV];
  for(int i=0;i<NV;i++) {
    force[i]= 0.0;
    forcer[i]= 0.0;
    rho[i]=0.0;
    vmoms[i]=0.0;
  }

  // free mass is used in the calculation of potential and force.
  if(optPotentialArg >= 1) {
    for(int i=0;i<NV;i++) part[i]->setPotS(0.0);
  }

  if(optTwoBodyDistance == 1) {
    qmdMatrix1();
    singleParticlePotential();
    computeForce1();
  } else if(optTwoBodyDistance==2) {
    qmdMatrix2();
    singleParticlePotential();
    computeForce2();
  } else {
    qmdMatrix3();
    singleParticlePotential();
    computeForce3();
  }

  for(int i=0; i<NV;i++)  {
    part[i]->updateByForce(force[i], forcer[i],dt);
  }

  delete [] rho;
  delete [] rhog;
  delete [] vmoms;
  for(int i=0;i<NV;i++) delete [] rhom[i];
  delete [] rhom;
  delete [] force;
  delete [] forcer;
}

double RQMDs0::facSDensity(Vec4& p, double m)
{
  if(optScalarDensity==0) {
    return 1.0;
  } else if(optScalarDensity==1) {
    return p.mCalc()/p[0];
  } else if (optScalarDensity==2) {
    return m/sqrt(m*m + p.pAbs2());
  } else {
    return 1.0;
  }

}

void RQMDs0::singleParticlePotential()
{
  bool optSet= optPotentialArg == 0 ? true :  false;

  for(int i=0; i< NV; i++) {
    double vsky = part[i]->baryon()/3*(t1 + t3*rhog[i])*rho[i];
    part[i]->setPotS(vsky + vmoms[i],optSet);
    part[i]->setRhoS(rho[i]);
    part[i]->setPotSm(vmoms[i]);
  }
}

// compute single particle potential energy.
Vec4 RQMDs0::computeEnergy(list<EventParticle*>& plist, int step)
{
  pTot=0.0;   
  for(auto i = plist.begin(); i != plist.end(); ++i) {
    pTot += (*i)->getP();
  }

  if(step==1) pTot0 = pTot;

  /*
  double econ=abs(pTot0[0]-pTot[0])/pTot0[0]*100;
  cout << "RQMDs0 time= " << fixed << step*dt
      << " econ= " << fixed << econ << " %"
     << scientific << setw(13) << pTot[1] 
     << scientific << setw(13) << pTot[2]
     << scientific << setw(13) << pTot[3]
     <<endl;
     */

  return pTot;
}

// non-relativistic two-body distance
// (distance at the global cm frame where the total momentum is zero)
void RQMDs0::qmdMatrix1()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;
    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;
      double drsq = (r1 - r2).pAbs2();
      double den = facG * exp(-drsq*wG)*qfac1*qfac2;
      rhom[i][j] = den;
      rhom[j][i] = den;
      rho[i] += rhom[i][j];
      rho[j] += rhom[j][i];
      if(!withMomDep) continue;
      Vec4 dp = p1 - p2;
      double ps = dp.pAbs2() - (1-optPV)*dp[0]*dp[0];
      double  pmom2ij = vex1/(1.0+ps/pmu1)+vex2/(1.0+ps/pmu2);
      vmoms[i] += pmom2ij*rhom[i][j];
      vmoms[j] += pmom2ij*rhom[j][i];
    }
  }

  for(int i=0;i<NV;i++) rhog[i] = pow(rho[i],gam-1);

}

// two-body distance at their center-of-momentum frame.
void RQMDs0::qmdMatrix2()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;
    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;
      Vec4 dr = r1 - r2; dr[0] = 0.0;  // just in case;
      Vec4 dp = p1 - p2;
      Vec4 P  = p1 + p2;
      double s = P.m2Calc();
      double drcmsq = dr.m2Calc() - pow2(dr*P)/s;
      double den = facG * exp(drcmsq*wG)*qfac1*qfac2;
      rhom[i][j] = den;
      rhom[j][i] = den;
      rho[i] += rhom[i][j];
      rho[j] += rhom[j][i];

      if(!withMomDep) continue;
      double ps = dp.m2Calc() - optPV * pow2(dp*P)/s;
      double  pmom2ij = vex1/(1.0-ps/pmu1)+vex2/(1.0-ps/pmu2);
      vmoms[i] += pmom2ij*rhom[i][j];
      vmoms[j] += pmom2ij*rhom[j][i];
    }
  }
  for(int i=0;i<NV;i++) rhog[i] = pow(rho[i],gam-1);

}

// two-body distance measured from the rest-frame of particle j.
void RQMDs0::qmdMatrix3()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double m1sq = p1.m2Calc();
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;
    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double m2sq = p2.m2Calc();
      double qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;

      Vec4 dr = r1 - r2; dr[0] = 0.0;  // just in case;
      double drsq1 = dr.m2Calc() - pow2(dr*p2)/m2sq;
      rhom[i][j] = facG * exp(drsq1*wG)*qfac1*qfac2;
      double drsq2 = dr.m2Calc() - pow2(dr*p1)/m1sq;
      rhom[j][i] = facG * exp(drsq2*wG)*qfac1*qfac2;
      rho[i] += rhom[i][j];
      rho[j] += rhom[j][i];

      if(!withMomDep) continue;

      //Vec4 dp = p1 - p2;
      //double dpsq = dp.m2Calc();
      //dot4 = dp * p2;
      //double ps1 = dpsq - optPV * dot4*dot4/m2sq;
      //dot4 = dp * p1;
      //double ps2 = dpsq - optPV * dot4*dot4/m1sq;

      Vec4 dp = p1 - p2;
      double dpsq = dp.m2Calc();
      double ps1 = dpsq;
      double ps2 = dpsq;
      if(optPV !=0 ) {
      if(optMDarg3==1)  {
        dpsq = dp.pAbs2();
	ps1 = -dpsq; ps2= -dpsq;
      } else if(optMDarg3==2) {
        double s = (p1+p2).m2Calc();
        dpsq += - optPV * pow2(dp * (p1+p2))/s;
	ps1 = dpsq; ps2= dpsq;
      } else {
        ps1 = dpsq - optPV * pow2(dp * p2)/m2sq;
        ps2 = dpsq - optPV * pow2(dp * p1)/m1sq;
      }
      }
      double  pmom2ij = vex1/(1.0-ps1/pmu1)+vex2/(1.0-ps1/pmu2);
      double  pmom2ji = vex1/(1.0-ps2/pmu1)+vex2/(1.0-ps2/pmu2);
      vmoms[i] += pmom2ij*rhom[i][j];
      vmoms[j] += pmom2ji*rhom[j][i];
    }
  }
  for(int i=0;i<NV;i++) rhog[i] = pow(rho[i],gam-1);

}

// non-relativistic distance.
void RQMDs0::computeForce1()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double emi = part[i]->getMass();
    double ei= optPotentialArg>0 ? p1[0] : sqrt(pow2(part[i]->getEffectiveMass())+p1.pAbs2());
    double fengi=emi/ei;

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double emj = part[i]->getMass();
      double ej= optPotentialArg>0 ? p2[0] : sqrt(pow2(part[j]->getEffectiveMass())+p2.pAbs2());
      double fengj=emj/ej;

      double fsky1 = fengi*(t1 + t3f*rhog[i])*rhom[i][j]; 
      double fsky2 = fengj*(t1 + t3f*rhog[j])*rhom[j][i]; 
      double fsky = -wG*(fsky1 + fsky2);

      //test
      //fsky *= 2.;

      force[i] += -2*fsky*(r1 - r2);
      force[j] += -2*fsky*(r2 - r1);


      if(!withMomDep) continue;

      Vec4 dp = p1 - p2;
      // optPV=0: distance is (p_i - p_j)^2 4-distance
      //double psq = dp.pAbs2() - (1-optPV) * dp[0]*dp[0];
      double psq = dp.pAbs2();
      double fac1 = 1.0 + psq/pmu1;
      double fac2 = 1.0 + psq/pmu2;

      // p derivative term.
      double fengij = fengi*rhom[i][j] + fengj*rhom[j][i];
      double fmome = -fengij*(vex1/(pmu1*fac1*fac1)+vex2/(pmu2*fac2*fac2));

      // r derivative term.
      double fmomd = -wG*fengij*(vex1/fac1 + vex2/fac2);

      //Vec4 dp2pi =  dp -optP0dev*(1-optPV)*dp[0]*p1/p1[0];
      //Vec4 dp2pj = -dp -optP0dev*(1-optPV)*dp[0]*p2/p2[0];

      force[i]  += -2*fmomd*(r1 - r2);
      force[j]  += -2*fmomd*(r2 - r1);
      forcer[i] +=  2*fmome*dp;
      forcer[j] += -2*fmome*dp;
    }
  }

}

// two-body distance is defined by the c.m. frame of two-particles.
void RQMDs0::computeForce2()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double emi = part[i]->getMass();
    double ei= optPotentialArg>0 ? p1[0] : sqrt(pow2(part[i]->getEffectiveMass())+p1.pAbs2());
    double fengi=emi/ei;

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double emj = part[j]->getMass();
      double ej= optPotentialArg>0 ? p2[0] : sqrt(pow2(part[j]->getEffectiveMass())+p2.pAbs2());
      double fengj=emj/ej;

      double fsky1 = fengi*(t1 + t3f*rhog[i])*rhom[i][j]; 
      double fsky2 = fengj*(t1 + t3f*rhog[j])*rhom[j][i]; 
      double fsky = -wG*(fsky1 + fsky2);

      //test
      //fsky *= 2.;

      Vec4 pcm = p1 + p2;
      Vec4 bi  = p1/p1[0];
      Vec4 bj  = p2/p2[0];
      Vec4 bbi = pcm/pcm[0] - optP0dev*bi;
      Vec4 bbj = pcm/pcm[0] - optP0dev*bj;

      Vec4 dr  = r1 - r2;
      double s = pcm.m2Calc();
      double rbij = -dr * pcm/s;
      Vec4 dr2ri = dr + rbij*pcm;
      Vec4 dr2rj = -dr2ri;
      Vec4 dr2pi = (dr + pcm[0]*rbij*bbi)*rbij;
      Vec4 dr2pj = (dr + pcm[0]*rbij*bbj)*rbij;

      force[i]  += -2*fsky*dr2ri;
      force[j]  += -2*fsky*dr2rj;
      forcer[i] +=  2*fsky*dr2pi;
      forcer[j] +=  2*fsky*dr2pj;

      if(!withMomDep) continue;

      Vec4 dp = p1 - p2;
      //double psq = dp.m2Calc() - optPV * pow2(dp*pcm)/s;
      double pma = pow2(emi*emi - emj*emj)/s;
      double psq = dp.m2Calc() - optPV * pma;
      double fac1 = 1.0 - psq/pmu1;
      double fac2 = 1.0 - psq/pmu2;

      // p derivative term.
      double fengij = fengi*rhom[i][j] + fengj*rhom[j][i];
      double fmome = -fengij*(vex1/(pmu1*fac1*fac1)+vex2/(pmu2*fac2*fac2));

      // r derivative term.
      double facmom = vex1/fac1 + vex2/fac2;
      double fmomd = -wG*fengij*facmom;
      Vec4 dp2pi =  dp -optP0dev*dp[0]*bi + optPV*pcm[0]/s*pma*bbi;
      Vec4 dp2pj = -dp +optP0dev*dp[0]*bj + optPV*pcm[0]/s*pma*bbj;

      force[i]  += -2*fmomd*dr2ri;
      force[j]  += -2*fmomd*dr2rj;
      forcer[i] +=  2*fmomd*dr2pi + 2*fmome*dp2pi;
      forcer[j] +=  2*fmomd*dr2pj + 2*fmome*dp2pj;

    }
  }

}

// two-body distance is defined by the rest frame of particle.
void RQMDs0::computeForce3()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double emi = part[i]->getMass();
    double ei= optPotentialArg>0 ? p1[0] : sqrt(pow2(part[i]->getEffectiveMass())+p1.pAbs2());
    double fengi=emi/ei;
    Vec4 vi = p1 / emi;

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double emj = part[j]->getMass();
      double ej= optPotentialArg>0 ? p2[0] : sqrt(pow2(part[j]->getEffectiveMass())+p2.pAbs2());
      double fengj=emj/ej;

      double fskyi = -wG*fengi*(t1 + t3f*rhog[i])*rhom[i][j]; 
      double fskyj = -wG*fengj*(t1 + t3f*rhog[j])*rhom[j][i]; 

      Vec4 vj = p2 / emj;
      Vec4 dr = r1 - r2; dr[0]=0.0;
      double rbi = dot3(dr,vi);
      double rbj = dot3(dr,vj);
      Vec4 dr2ri =  2*(dr + rbj*vj);
      Vec4 dr2rj =  2*(dr + rbi*vi);
      Vec4 dr2pi =  2*dr*rbi/emi;
      Vec4 dr2pj =  2*dr*rbj/emj;

      force[i]  += -fskyi*dr2ri - fskyj*dr2rj;
      force[j]  +=  fskyi*dr2ri + fskyj*dr2rj;
      forcer[i] +=  fskyj*dr2pi;
      forcer[j] +=  fskyi*dr2pj;

      if(!withMomDep) continue;

      if(optMDarg3==1) distanceP1(p1,p2);
      else if(optMDarg3==2) distanceP2(p1,p2);
      else distanceP3(p1,p2);

      //devVmom(psq2);
      double fmomei = fengi*devVme(psq2)*rhom[i][j];
      double fmomdi = -wG*fengi*devVmd(psq2)*rhom[i][j];
      double fmomej = fengj*devVme(psq1)*rhom[j][i];
      double fmomdj = -wG*fengj*devVmd(psq1)*rhom[j][i];

      force[i]  += -fmomdi*dr2ri - fmomdj*dr2rj;
      force[j]  +=  fmomdi*dr2ri + fmomdj*dr2rj;
      forcer[i] +=  fmomei*dp2ijpi + fmomej*dp2jipi + fmomdj*dr2pi;
      forcer[j] +=  fmomej*dp2jipj + fmomei*dp2ijpj + fmomdi*dr2pj;

    }
  }

}

// relative distance and its derivatives between p1 and p2
void RQMDs0::distanceP1(const Vec4& p1,const Vec4& p2)
{
  Vec4 dp = p1 - p2;
  psq1 = -dp.pAbs2() + (1-optPV) * dp[0]*dp[0];
  psq2 = psq1;
  dp2ijpi = 2*( dp -optP0dev*(1-optPV)*dp[0]*p1/p1.e());
  dp2ijpj = 2*(-dp -optP0dev*(1-optPV)*dp[0]*p2/p2.e());
  dp2jipi = dp2ijpi;
  dp2jipj = dp2ijpj;
}

// relative distance and its derivatives between p1 and p2
// in the two body CM frame.
void RQMDs0::distanceP2(const Vec4& p1,const Vec4& p2)
{
  Vec4 dp  = p1 - p2;
  Vec4 pcm = p1 + p2;
  double s = pcm.m2Calc();
  //double pma = pow2(m1*m1 - m2*m2)/s;
  double pma = pow2(dp * pcm)/s;
  psq1 = dp.m2Calc() - optPV * pma;
  psq2 = psq1;
  Vec4 b1 = p1/p1.e();
  Vec4 b2 = p2/p2.e();
  Vec4 bbi = pcm/pcm[0] - optP0dev*b1;
  Vec4 bbj = pcm/pcm[0] - optP0dev*b2;
  dp2ijpi = 2*(dp - optP0dev*dp[0]*b1 + optPV * pcm[0]/s*pma*bbi);
  dp2ijpj = 2*(-dp + optP0dev*dp[0]*b2 + optPV * pcm[0]/s*pma*bbj);
  dp2jipi = dp2ijpi;
  dp2jipj = dp2ijpj;
}

// relative distance squared and its derivatives between p1 and p2
// in the rest frame of p1 or p2.
void RQMDs0::distanceP3(const Vec4& p1,const Vec4& p2)
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



} // namespace jam2


