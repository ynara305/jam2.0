#include <jam2/meanfield/MBoxS.h>

namespace jam2 {

void NeighborS::resize(int n) 
{
  int m = box->particleSize();
  rhom.resize(n);
  for(int i=0;i<n;i++) rhom[i].resize(m);
}

MBoxS:: ~MBoxS()
{
  for(auto& i : neighbors) delete i;
  neighbors.clear();
  clearMatrix();
}

bool MBoxS::haveThisSite(MBoxS* box) const
{
  auto first = neighbors.begin();
  auto last = neighbors.end();

  while(first != last) {
    if((*first)->box == box) return true;
    ++first;
  }
  return false;
}

//std::vector<NeighborS>::iterator MBoxS::findNeighbor(MBoxS* box)
NeighborS* MBoxS::findNeighbor(MBoxS* box)
{
  auto&& first = neighbors.begin();
  auto&& last = neighbors.end();

  while(first != last) {
    if((*first)->box == box) return *first;
    ++first;
  }
  cout << "MBoxS error no such neighbor "<<endl;
  exit(1);
}

void MBoxS::clearMatrix()
{
  part.clear();
  rho.clear();
  rhog.clear();
  vmoms.clear();
  force.clear();
  forcer.clear();
  for(int i=0;i<NV;i++) rhom[i].clear();
  rhom.clear();
}

void MBoxS::initBox()
{
  NV = part.size();
  if(NV==0) return;

  rhog.resize(NV);
  rhom.resize(NV);
  for(int i=0;i<NV;i++) rhom[i].resize(NV);

  //std::fill(rho.begin(),rho.end(),0.0);
  //std::fill(vmoms.begin(),vmoms.end(),0.0);

  rho.assign(NV,0.0);
  vmoms.assign(NV,0.0);

  force.assign(NV,0.0);
  forcer.assign(NV,0.0);

  //for(int i=0;i<NV;i++) { force[i]= 0.0; forcer[i]= 0.0; }

  // free mass is used in the calculation of potential and force.
  if(optPotentialArg >= 1) {
    for(auto& i : part) i->setPotS(0.0);
  }

  for(auto& neighbor : neighbors) {
    neighbor->resize(NV);
  }

}

void MBoxS::qmdMatrix(double a)
{
  if(NV==0) return;

  RQMDs::qmdMatrix(a);
  for(auto& neighbor : neighbors) {
    if(neighbor->getAction()) qmdMatrixNeighbor(*neighbor,a);
  }

  /*
  if(optTwoBodyDistance == 1) {
    qmdMatrix1(a);
    for(auto& neighbor : neighbors) {
      if(neighbor->getAction()) qmdMatrixNeighbor1(*neighbor,a);
    }
  } else if(optTwoBodyDistance == 2) {
    qmdMatrix2(a);
    for(auto& neighbor : neighbors) {
      if(neighbor->getAction()) qmdMatrixNeighbor2(*neighbor,a);
    }
  } else if(optTwoBodyDistance == 3) {
    qmdMatrix3(a);
    for(auto& neighbor : neighbors) {
      if(neighbor->getAction()) qmdMatrixNeighbor3(*neighbor,a);
    }
  }
  */

}

void MBoxS::computeForce(double dt)
{
  if(NV==0) return;

  RQMDs::computeForce();
  for(auto& neighbor : neighbors) {
    if(neighbor->getAction()) computeForceNeighbor(*neighbor);
  }

  /*
  if(optTwoBodyDistance == 1) {
    computeForce1();
    for(auto& neighbor : neighbors) {
      if(neighbor->getAction()) computeForceNeighbor1(*neighbor);
    }
  } else if(optTwoBodyDistance==2) {
    computeForce2();
    for(auto& neighbor : neighbors) {
      if(neighbor->getAction()) computeForceNeighbor2(*neighbor);
    }
  } else {
    computeForce3();
    for(auto& neighbor : neighbors) {
      if(neighbor->getAction()) computeForceNeighbor3(*neighbor);
    }
  }
  */


  for(int i=0;i<NV;i++) {
    part[i]->updateByForce(force[i],forcer[i],dt);
  }

  if(optPotentialArg >= 1) singleParticlePotential(true);

}

void MBoxS::qmdMatrixNeighbor(NeighborS& neighbor,double a)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    p1 *= a;
    //double m1 = part[i]->getEffectiveMass();
    //p1[0]= sqrt(m1*m1 + p1.pAbs2());
    distance->setRP1(r1,p1);
    double fi=p1.mCalc()/p1[0];
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;

    for(auto j=0; j< n2; j++) {
      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      p2 *=a;
      //double m2 = part2[j]->getEffectiveMass();
      //p2[0]= sqrt(m2*m2+p2.pAbs2());
      distance->setRP2(r2,p2);
      double fj=p2.mCalc()/p2[0];
      double qfac2 = part2[j]->getTf() > globalTime ? part2[j]->qFactor() : 1.0;

      distance->density();
      double rhomij = distance->density1*qfac1*qfac2;
      double rhomji = distance->density2*qfac1*qfac2;
      neighbor.setRhom(i,j,rhomij);
      neighbor.neighbor->setRhom(j,i, rhomji);

      rho[i] += rhomij*fj;
      neighbor.box->addRho(j,rhomji*fi);

      if(!withMomDep) continue;

      distance->psq();
      double  pmom2ij = vex1/(1.0-distance->psq1/pmu1)+vex2/(1.0-distance->psq1/pmu2);
      double  pmom2ji = vex1/(1.0-distance->psq2/pmu1)+vex2/(1.0-distance->psq2/pmu2);
      vmoms[i] += pmom2ij*rhomij*fj;
      neighbor.box->addVmoms(j, pmom2ji*rhomji*fi);

    }
  }

}

void MBoxS::computeForceNeighbor(NeighborS& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double fi = p1.mCalc()/p1[0];
    double fengi=fi;
    if(optPotentialArg==3) {
      double meff1 = part[i]->getEffectiveMass();
      fengi = meff1/sqrt(meff1*meff1 + p1.pAbs2());
    }
    distance->setRP1(r1,p1);

    for(auto j=0; j< n2; j++) {
      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      double fj=p2.mCalc()/p2[0];
      double fengj=fj;
      if(optPotentialArg==3) {
        double meff2 = part2[j]->getEffectiveMass();
        fengj = meff2/sqrt(meff2*meff2 + p2.pAbs2());
      }
      distance->setRP2(r2,p2);

      double rhomij=neighbor.rhom[i][j];
      double rhomji=neighbor.neighbor->rhom[j][i];
      double fskyi = -wG*fengi*(t1 + t3f*rhog[i])*rhomij*fj; 
      double fskyj = -wG*fengj*(t1 + t3f*neighbor.box->getRhog(j))*rhomji*fi; 

      distance->distanceR();
      force[i]  += -fskyi*distance->dr2ijri - fskyj*distance->dr2jirj;
      forcer[i] +=  fskyi*distance->dr2ijpi + fskyj*distance->dr2jipi;
      neighbor.box->addForce(j,-fskyi*distance->dr2ijrj - fskyj*distance->dr2jirj);
      neighbor.box->addForceR(j,fskyi*distance->dr2ijpj + fskyj*distance->dr2jipj);

      if(optDerivative && optP0dev) {
        forcer[i] +=  fskyj*p1/(wG*p1[0]*p1[0]);
        neighbor.box->addForceR(j,fskyi*p2/(wG*p2[0]*p2[0]));

        // gamma derivative.
        distance->devGamma();
        double facsk= -(fskyi+fskyj)/wG;
        forcer[i] += distance->devgam1*facsk;
        neighbor.box->addForceR(j,distance->devgam2*facsk);
      }

      if(!withMomDep) continue;

      distance->distanceP();

      double fmomdi = -wG*fengi*devVmd(distance->psq2)*rhomij*fj;
      double fmomei =     fengi*devVme(distance->psq2)*rhomij*fj;
      double fmomdj = -wG*fengj*devVmd(distance->psq1)*rhomji*fi;
      double fmomej =     fengj*devVme(distance->psq1)*rhomji*fi;

      force[i]  += -fmomdi*distance->dr2ijri - fmomdj*distance->dr2jiri;
      neighbor.box->addForce(j,-fmomdi*distance->dr2ijrj + fmomdj*distance->dr2jirj);
      forcer[i] +=  fmomei*distance->dp2ijpi + fmomej*distance->dp2jipi
	          + fmomdi*distance->dr2ijpi + fmomdj*distance->dr2jipi;
      neighbor.box->addForceR(j,fmomei*distance->dp2ijpj + fmomej*distance->dp2jipj
	                      + fmomdi*distance->dr2ijpj + fmomdj*distance->dr2jipj);

      if(optDerivative && optP0dev) {
        forcer[i] +=  fmomdj*p1/(wG*p1[0]*p1[0]);
        neighbor.box->addForceR(j,fmomdi*p2/(wG*p2[0]*p2[0]));

        // gamma derivative.
        double facmom = -(fmomdi + fmomdj)/wG;
        forcer[i] += distance->devgam1*facmom;
        neighbor.box->addForceR(j,distance->devgam2*facmom);
      }

    }
  }

}

// non-relativistic two-body distance
// (distance at the global cm frame where the total momentum is zero)
void MBoxS::qmdMatrixNeighbor1(NeighborS& neighbor, double a)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(int i=0; i<NV; i++) {

    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double m1 = part[i]->getEffectiveMass();
    //if(optPotentialArg >=1 ) m1 = part[i]->getMass();
    p1 *= a;
    p1[0]= sqrt(m1*m1+p1.pAbs2());
    double fi=m1/p1[0];

    for(int j=0; j< n2; j++) {

      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      double m2 = part2[j]->getEffectiveMass();
      //if(optPotentialArg >=1 ) m2 = part2[j]->getMass();
      p2 *=a;
      p2[0]= sqrt(m2*m2+p2.pAbs2());
      double fj=m2/p2[0];
      //double drsq = (r1 - r2).pAbs2();
      double drsq = (r1 - r2).pT2() + pow2(gamCM*(r1[3]-r2[3]));
      double den = gamCM*facG * exp(-drsq*wG);
      neighbor.setRhom(i,j, den);
      rho[i] += den*fj;

      neighbor.neighbor->setRhom(j,i, den);
      neighbor.box->addRho(j,den*fi);

      if(!withMomDep) continue;
      Vec4 dp = p1 - p2;
      double ps = dp.pAbs2() - (1-optPV)*dp[0]*dp[0];
      double  pmom2ij = vex1/(1.0+ps/pmu1)+vex2/(1.0+ps/pmu2);
      vmoms[i] += pmom2ij*den*fj;

      neighbor.box->addVmoms(j, pmom2ij*den*fi);

    }
  }

}

// two-body distance at their center-of-momentum frame.
void MBoxS::qmdMatrixNeighbor2(NeighborS& neighbor, double a)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double m1 = part[i]->getEffectiveMass();
    //double m1 = p1.mCalc();
    //if(optPotentialArg >=1 ) m1 = part[i]->getMass();
    p1 *= a;
    p1[0]= sqrt(m1*m1+p1.pAbs2());
    double fi = m1/p1[0];
    for(auto j=0; j< n2; j++) {
      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      double m2 = part2[j]->getEffectiveMass();
      //double m2 = p2.mCalc();
      //if(optPotentialArg >=1 ) m2 = part2[j]->getMass();
      p2 *=a;
      p2[0]= sqrt(m2*m2+p2.pAbs2());
      double fj=m2/p2[0];
      Vec4 dr = r1 - r2; dr[0] = 0.0;  // just in case;
      Vec4 dp = p1 - p2;
      Vec4 P  = p1 + p2;
      double s = P.m2Calc();
      double drcmsq = dr.m2Calc() - pow2(dr*P)/s;
      double g12 = optScalarDensity > 0 ? P[0]/sqrt(s) : 1.0;
      double den = g12 * facG * exp(drcmsq*wG);
      neighbor.setRhom(i,j,den);
      neighbor.neighbor->setRhom(j,i,den);
      rho[i] += den*fj;
      neighbor.box->addRho(j, den*fi);

      if(!withMomDep) continue;
      double ps = dp.m2Calc() - optPV * pow2(dp*P)/s;
      double  pmom2ij = den*(vex1/(1.0-ps/pmu1)+vex2/(1.0-ps/pmu2));
      vmoms[i] += pmom2ij*fj;
      neighbor.box->addVmoms(j, pmom2ij*fi);
    }
  }

}

// two-body distance measured from the rest-frame of particle j.
void MBoxS::qmdMatrixNeighbor3(NeighborS& neighbor,double a)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    //double m1sq = p1.m2Calc();
    //if(optPotentialArg >=1 ) m1 = part[i]->getMass();
    double m1 = part[i]->getEffectiveMass();
    p1 *= a;
    double m1sq = m1*m1;
    p1[0]= sqrt(m1sq + p1.pAbs2());

    for(auto j=0; j< n2; j++) {
      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      //double m2sq = p2.m2Calc();
      //if(optPotentialArg >=1 ) m2 = part2[j]->getMass();
      double m2 = part2[j]->getEffectiveMass();
      double m2sq = m2*m2;
      p2 *=a;
      p2[0]= sqrt(m2*m2+p2.pAbs2());

      Vec4 dr = r1 - r2; dr[0] = 0.0;  // just in case;
      double drsq1 = dr.m2Calc() - pow2(dr * p2)/m2sq;
      double rhomij = facG * exp(drsq1*wG);
      double drsq2 = dr.m2Calc() - pow2(dr * p1)/m1sq;
      double rhomji = facG * exp(drsq2*wG);
      neighbor.setRhom(i,j,rhomij);
      neighbor.neighbor->setRhom(j,i,rhomji);

      rho[i] += rhomij;
      neighbor.box->addRho(j,rhomji);

      if(!withMomDep) continue;

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
      vmoms[i] += (vex1/(1.0-ps1/pmu1)+vex2/(1.0-ps1/pmu2))*rhomij;
      neighbor.box->addVmoms(j,(vex1/(1.0-ps2/pmu1)+vex2/(1.0-ps2/pmu2))*rhomji);
    }
  }

}

void MBoxS::computeForceNeighbor1(NeighborS& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double emi = p1.mCalc();
    double fi = emi/p1[0];
    double meff1 = part[i]->getEffectiveMass();
    double fengi = meff1/sqrt(meff1*meff1 + p1.pAbs2());

    for(int j=0; j<n2; j++) {

      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      double emj = p2.mCalc();
      double fj=emj/p2[0];
      double meff2 = part2[j]->getEffectiveMass();
      double fengj = meff2/sqrt(meff2*meff2 + p2.pAbs2());

      double rhomij=neighbor.rhom[i][j];
      double rhomji=neighbor.neighbor->rhom[j][i];

      double fsky1 = fengi*(t1 + t3f*rhog[i])*rhomij;
      double fsky2 = fengj*(t1 + t3f*neighbor.box->getRhog(j) )*rhomji; 
      double fsky = -wG*(fsky1*fj + fsky2*fi);
      Vec4 dr = r1 - r2;
      dr[3] *= gamCM2;
      force[i] += -2*fsky*dr;
      neighbor.box->addForce(j, 2*fsky*dr);

      // Derivative of m_j/p^0_j in the scalar density.
      if(optScalarDensity > 0 && optP0dev) {
        forcer[i] += -fsky2 * fi /p1[0] * p1/p1[0];
        neighbor.box->addForceR(j,-fsky1 * fj /p2[0] * p2/p2[0]);
      }

      if(!withMomDep) continue;

      // optPV=0: distance is (p_i - p_j)^2 4-distance
      Vec4 dp = p1 - p2;
      double psq = dp.pAbs2() - (1-optPV) * dp[0]*dp[0];
      double fac1 = 1.0 + psq/pmu1;
      double fac2 = 1.0 + psq/pmu2;

      // p derivative term.
      double fengij = fengi*fj*rhomij + fengj*fi*rhomji;
      double fmome = -fengij*(vex1/(pmu1*fac1*fac1)+vex2/(pmu2*fac2*fac2));

      // r derivative term.
      double facmom = vex1/fac1 + vex2/fac2;
      double fmomd = -wG*fengij*facmom;

      force[i]  += -2*fmomd*dr;
      forcer[i] +=  2*fmome*( dp -optP0dev*(1-optPV)*dp[0]*p1/p1[0]);

      neighbor.box->addForce( j,  2*fmomd*dr);
      neighbor.box->addForceR(j,  2*fmome*(-dp -optP0dev*(1-optPV)*dp[0]*p2/p2[0]));

      // Derivative of m_j/p^0_j in the scalar density.
      if(optScalarDensity > 0 && optP0dev) {
        double fs1 = -facmom*fengj * fi / p1[0] * rhomji;
        forcer[i] += fs1 * p1/p1[0];

        double fs2 = -facmom*fengi * fj / p2[0] * rhomij;
        neighbor.box->addForceR(j,fs2*p2/p2[0]);
      }
    }
  }

}

// two-body distance is defined by the c.m. frame of two-particles.
void MBoxS::computeForceNeighbor2(NeighborS& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double emi = p1.mCalc();
    double fi = emi/p1[0];
    double meff1 = part[i]->getEffectiveMass();
    double fengi = meff1/sqrt(meff1*meff1 + p1.pAbs2());

    for(auto j=0; j< n2; j++) {
      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      double emj = p2.mCalc();
      double fj=emj/p2[0];
      double meff2 = part2[j]->getEffectiveMass();
      double fengj = meff2/sqrt(meff2*meff2 + p2.pAbs2());

      double rhomij=neighbor.rhom[i][j];
      double rhomji=neighbor.neighbor->rhom[j][i];

      double fsky1 = fengi*(t1 + t3f*rhog[i])*rhomij; 
      double fsky2 = fengj*(t1 + t3f*neighbor.box->getRhog(j))*rhomji; 
      double fsky = -wG*(fsky1*fj + fsky2*fi);

      Vec4 pcm = p1 + p2;
      Vec4 bi  = p1/p1[0];
      Vec4 bj  = p2/p2[0];
      Vec4 bbi = pcm/pcm[0] - optP0dev*bi;
      Vec4 bbj = pcm/pcm[0] - optP0dev*bj;

      Vec4 dr  = r1 - r2; dr[0]=0.0;
      double s = pcm.m2Calc();
      double rbij = -dr * pcm/s;
      Vec4 dr2ri = dr + rbij*pcm;
      Vec4 dr2rj = -dr2ri;
      Vec4 dr2pi = (dr + pcm[0]*rbij*bbi)*rbij;
      Vec4 dr2pj = (dr + pcm[0]*rbij*bbj)*rbij;

      force[i]  += -2*fsky*dr2ri;
      neighbor.box->addForce(j, -2*fsky*dr2rj);
      forcer[i] +=  2*fsky*dr2pi;
      neighbor.box->addForceR(j, 2*fsky*dr2pj);

      // Derivative of gamma_ij in front of Gaussian.
      if(optScalarDensity > 0) {
        double facsk = fj * fsky1 + fi * fsky2;
	double fsgam1 = facsk / s;
	double fsgam2 = optP0dev*facsk*(1.0/pcm[0] - pcm[0]/s);
	forcer[i] += fsgam1*pcm + fsgam2*bi;
	neighbor.box->addForceR(j, fsgam1*pcm + fsgam2*bj);

        // Derivative of m_j/p^0_j in the scalar density.
        double fs1 = -fsky2 * fi / p1[0] * optP0dev;
        double fs2 = -fsky1 * fj / p2[0] * optP0dev;
        forcer[i] += fs1 * bi;
	neighbor.box->addForceR(j, fs2 * bj);
      }

      if(!withMomDep) continue;

      Vec4 dp = p1 - p2;
      double pma = pow2(emi*emi - emj*emj)/s;
      //double psq = dp.m2Calc() - optPV * pow2(dp * pcm)/s;
      double psq = dp.m2Calc() - pma;

      double fac1 = 1.0 - psq/pmu1;
      double fac2 = 1.0 - psq/pmu2;

      // p derivative term.
      double fengij = fengi*fj*rhomij + fengj*fi*rhomji;
      double fmome = -fengij*(vex1/(pmu1*fac1*fac1)+vex2/(pmu2*fac2*fac2));

      // r derivative term.
      double facmom = vex1/fac1 + vex2/fac2;
      double fmomd = -wG*fengij*facmom;
      Vec4 dp2pi =  dp -optP0dev*dp[0]*bi + optPV*pcm[0]/s*pma*bbi;
      Vec4 dp2pj = -dp +optP0dev*dp[0]*bj + optPV*pcm[0]/s*pma*bbj;

      force[i]  += -2*fmomd*dr2ri;
      neighbor.box->addForce(j, -2*fmomd*dr2rj);
      forcer[i] +=  2*fmomd*dr2pi + 2*fmome*dp2pi;
      neighbor.box->addForceR(j, 2*fmomd*dr2pj + 2*fmome*dp2pj);

      if(optScalarDensity > 0 ) {
	// derivative of gamma_ij
	double fmgam1 = fengij*facmom/s;
	double fmgam2i = optP0dev*fengij*facmom*(1.0/pcm[0] - pcm[0]/s);
	double fmgam2j = fmgam2i;

        // Derivative of m_j/p^0_j in the scalar density.
        fmgam2i += -facmom*fengj * fi / p1[0] * rhomji * optP0dev;
        fmgam2j += -facmom*fengi * fj / p2[0] * rhomij * optP0dev;

        forcer[i] += pcm*fmgam1 + bi*fmgam2i;
        neighbor.box->addForceR(j,pcm*fmgam1 + bj*fmgam2j);
      }

    }
  }

}


// two-body distance is defined by the rest frame of particle.
void MBoxS::computeForceNeighbor3(NeighborS& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double emi = p1.mCalc();
    Vec4 vi = p1 / emi;
    double meff1 = part[i]->getEffectiveMass();
    double fengi = meff1/sqrt(meff1*meff1 + p1.pAbs2());

    for(auto j=0; j< n2; j++) {
      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      double emj = p2.mCalc();
      Vec4 vj = p2 / emj;
      double meff2 = part2[j]->getEffectiveMass();
      double fengj = meff2/sqrt(meff2*meff2 + p2.pAbs2());

      double rhomij=neighbor.rhom[i][j];
      double rhomji=neighbor.neighbor->rhom[j][i];

      double fskyi = -wG*fengi*(t1 + t3f*rhog[i])*rhomij; 
      double fskyj = -wG*fengj*(t1 + t3f*neighbor.box->getRhog(j))*rhomji; 

      Vec4 dr = r1 - r2;dr[0]=0.0;
      double rbi = dot3(dr,vi);
      double rbj = dot3(dr,vj);
      Vec4 dr2ri =  2*(dr + rbj*vj);
      Vec4 dr2rj =  2*(dr + rbi*vi);
      Vec4 dr2pi =  2*dr*rbi/emi;
      Vec4 dr2pj =  2*dr*rbj/emj;

      force[i]  += -fskyi*dr2ri - fskyj*dr2rj;
      neighbor.box->addForce(j,  fskyi*dr2ri + fskyj*dr2rj);
      forcer[i] +=  fskyj*dr2pi;
      neighbor.box->addForceR(j, fskyi*dr2pj);

      if(!withMomDep) continue;

      if(optMDarg3==1) distanceP1(p1,p2);
      else if(optMDarg3==2) distanceP2(p1,p2);
      else distanceP3(p1,p2);

      double fmomdi = -wG*fengi*devVmd(psq2)*rhomij;
      double fmomei =     fengi*devVme(psq2)*rhomij;
      double fmomdj = -wG*fengj*devVmd(psq1)*rhomji;
      double fmomej =     fengj*devVme(psq1)*rhomji;

      force[i]  += -fmomdi*dr2ri   - fmomdj*dr2rj;
      neighbor.box->addForce(j, fmomdi*dr2ri   + fmomdj*dr2rj);
      forcer[i] +=  fmomei*dp2ijpi + fmomej*dp2jipi + fmomdj*dr2pi;
      neighbor.box->addForceR(j, fmomej*dp2jipj + fmomei*dp2ijpj + fmomdi*dr2pj);

    }
  }

}


} // end of namespace jam2
