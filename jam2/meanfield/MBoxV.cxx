#include <jam2/meanfield/MBoxV.h>

namespace jam2 {

void NeighborV::resize(int n) 
{
  int m = box->particleSize();
  rhom.resize(n);
  for(int i=0;i<n;i++) rhom[i].resize(m);
}

MBoxV:: ~MBoxV()
{
  for(auto& i : neighbors) delete i;
  neighbors.clear();
  clearMatrix();
}

bool MBoxV::haveThisSite(MBoxV* box) const
{
  auto first = neighbors.begin();
  auto last = neighbors.end();

  while(first != last) {
    if((*first)->box == box) return true;
    ++first;
  }
  return false;
}

NeighborV* MBoxV::findNeighbor(MBoxV* box)
{
  auto&& first = neighbors.begin();
  auto&& last = neighbors.end();

  while(first != last) {
    if((*first)->box == box) return *first;
    ++first;
  }
  cout << "MBoxV error no such neighbor "<<endl;
  exit(1);
}

void MBoxV::clearMatrix()
{
  part.clear();
  rho.clear();
  rhog.clear();
  vmom4.clear();
  force.clear();
  forcer.clear();
  Vdot.clear();
  JB.clear();
  for(int i=0;i<NV;i++) rhom[i].clear();
  rhom.clear();
}

void MBoxV::initBox()
{
  NV = part.size();
  if(NV==0) return;

  rhog.resize(NV);
  rhom.resize(NV);
  for(int i=0;i<NV;i++) rhom[i].resize(NV);

  rho.assign(NV,0.0);
  vmom4.assign(NV,0.0);
  JB.assign(NV,0.0);
  Vdot.assign(NV,0.0);

  force.assign(NV,0.0);
  forcer.assign(NV,0.0);

  // free mass is used in the calculation of potential and force.
  //if(optPotentialArg >= 1) {
  if(optVdot == 1) {
    //for(auto& i : part) i->setPotS(0.0);
    for(auto& i : part) i->setFree();
  }

  for(auto& neighbor : neighbors) {
    neighbor->resize(NV);
  }

}

void MBoxV::qmdMatrix()
{
  if(NV==0) return;

  RQMDv::qmdMatrix();
  for(auto& neighbor : neighbors) {
    if(neighbor->getAction()) qmdMatrixNeighbor(*neighbor);
  }

  /*
  if(optTwoBodyDistance == 1) {
    qmdMatrix1();
    for(auto& neighbor : neighbors) {
      if(neighbor->getAction()) qmdMatrixNeighbor1(*neighbor);
    }

  } else if(optTwoBodyDistance == 2) {
    qmdMatrix2();
    for(auto& neighbor : neighbors) {
      if(neighbor->getAction()) qmdMatrixNeighbor2(*neighbor);
    }
  } else if(optTwoBodyDistance == 3) {
    qmdMatrix3();
    for(auto& neighbor : neighbors) {
      if(neighbor->getAction()) qmdMatrixNeighbor3(*neighbor);
    }
  }
  */

}

void MBoxV::computeForce(double dt)
{
  if(NV==0) return;

  RQMDv::computeForce();
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




  // Note that optVdot=0,1: all momenta are canonical.
  // Canonical momenta update.
  if(optVdot == 0) {
    for(int i=0; i< NV; i++) {
      part[i]->updateByForce(force[i], forcer[i],dt);
    }

  // canonical momenta are used for the evaluation of forces,
  // after update canonical momenta, set kinetic momenta.
  } else if(optVdot == 1) {

    for(int i=0; i< NV; i++) {
      // canonical momenta update.
      part[i]->updateByForce(force[i], forcer[i],dt);
      // change to kinetic momenta.
      //part[i]->addP(-part[i]->potv());
      //part[i]->setOnShell();
      part[i]->setKinetic();
    }

  // compute time derivatives of vector potential V_i^\mu numerically.
  } else if(optVdot == 2) {

    // Vdot[] is already computed in singleParticlePotential().
    for(int i=0; i< NV; i++) {
      part[i]->updateByForce(force[i] - Vdot[i]/dt, forcer[i],dt);
    }

  // compute time derivatives of vector potential analytically.
  } else if(optVdot == 3) {

    RQMDv::computeVdot();
    for(int i=0; i< NV; i++) {
      part[i]->updateByForce(force[i] - Vdot[i], forcer[i],dt);
    }
  }

}

void MBoxV::qmdMatrixNeighbor(NeighborV& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(int i=0; i< NV; i++) {
    Vec4 r1  = part[i]->getR();
    Vec4 p1  = part[i]->getP();
    if(optPotentialArg>=1) p1 = part[i]->getPcan(optVdot);
    distance->setRP1(r1,p1);
    Vec4 v1 = p1/p1[0];

    Vec4 fri = optBaryonCurrent ? part[i]->forceR() : 0.0;
    int bi   = part[i]->baryon()/3;
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;

    for(int j=0; j< n2; j++) {
      Vec4 r2  = part2[j]->getR();
      Vec4 p2  = part2[j]->getP();
      if(optPotentialArg>=1) p2 = part2[j]->getPcan(optVdot);
      distance->setRP2(r2,p2);
      Vec4 v2 = p2/p2[0];

      Vec4 frj = optBaryonCurrent ? part2[j]->forceR() : 0.0;
      int bj = part2[j]->baryon()/3;
      double qfac2 = part2[j]->getTf() > globalTime ? part2[j]->qFactor() : 1.0;

      distance->density();
      double den1 = distance->density1*qfac1*qfac2;
      double den2 = distance->density2*qfac1*qfac2;
      neighbor.setRhom(i,j,den1);
      neighbor.neighbor->setRhom(j,i, den2);

      JB[i] += den1*(v2+frj)*bj;
      neighbor.box->addJB(j,den2*(v1+fri)*bi);

      if(!withMomDep) continue;

      distance->psq();
      double  pmom2ij = vex1/(1.0-distance->psq1/pmu1)+vex2/(1.0-distance->psq1/pmu2);
      double  pmom2ji = vex1/(1.0-distance->psq2/pmu1)+vex2/(1.0-distance->psq2/pmu2);
      vmom4[i] += pmom2ij*den1*v2;
      neighbor.box->addVmom4(j, pmom2ji*den2*v1);

    }
  }

}

void MBoxV::computeForceNeighbor(NeighborV& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = p1;
    if(optVdot <= 1) pk1 = part[i]->getPkin();
    if(optPotentialArg>=1) p1 = part[i]->getPcan(optVdot);
    int  bar1 = part[i]->baryon()/3;

    distance->setRP1(r1,p1);
    Vec4 v1 = p1/p1[0];

    for(auto j=0; j< n2; j++) {
      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      Vec4 pk2 = p2;
      if(optVdot <= 1) pk2 = part2[j]->getPkin();
      if(optPotentialArg>=1) p2 = part2[j]->getPcan(optVdot);
      int  bar2 = part2[j]->baryon()/3;

      distance->setRP2(r2,p2);
      Vec4 v2 = p2/p2[0];

      Vec4 Ai = RQMDv::facV(i,p2);
      Vec4 Aj = facV(j,p1,neighbor);
      double rhomij=neighbor.rhom[i][j];
      double rhomji=neighbor.neighbor->rhom[j][i];
      double fskyi = -wG*(Ai * pk1)/pk1[0]*rhomij*bar1*bar2;
      double fskyj = -wG*(Aj * pk2)/pk2[0]*rhomji*bar1*bar2;
      distance->distanceR();
      force[i]  += -fskyi*distance->dr2ijri - fskyj*distance->dr2jirj;
      forcer[i] +=  fskyi*distance->dr2ijpj + fskyj*distance->dr2jipj;
      neighbor.box->addForce(j,-fskyi*distance->dr2ijrj - fskyj*distance->dr2jirj);
      neighbor.box->addForceR(j,fskyi*distance->dr2ijpj + fskyj*distance->dr2jipj);


      // Derivative of p^\mu_j/p^0_j in the vector density.
      if(optDerivative) {
        Vec4 A1 = RQMDv::facV(i,pk1);
        Vec4 A2 = facV(j,pk2,neighbor);
        //double fsgam2i = optP0dev*rhomji*dot3(Aj,v1)/p1[0];
        //double fsgam2j = optP0dev*rhomij*dot3(Ai,v2)/p2[0];
        //double fsgam3i = -rhomji/p1[0];
        //double fsgam3j = -rhomij/p2[0];
        //forcer[i] += (fsgam2i*v1 + Aj*fsgam3i)*bar1*bar2;
        //neighbor.box->addForceR(j, (fsgam2j*v2 + Ai*fsgam3j)*bar1*bar2);

        distance->devV(A1,A2);
        forcer[i] += distance->devV1*rhomji*bar1*bar2;
        neighbor.box->addForceR(j,distance->devV2*rhomij*bar1*bar2);

        distance->devGamma();
        double facsk= -(fskyi+fskyj)/wG;
        forcer[i] += distance->devgam1*facsk;
        neighbor.box->addForceR(j,  distance->devgam2*facsk);

      }

      if(!withMomDep) continue;

      double vf1=1.0;
      double vf2=1.0;
      if(optVectorPotential==1) {
        vf1 = pk1/pk1[0] * v2;
        vf2 = pk2/pk2[0] * v1;
      }

      distance->distanceP();
      double fmomdi = -wG*devVmd(distance->psq2)*rhomij*vf1;
      double fmomdj = -wG*devVmd(distance->psq1)*rhomji*vf2;
      double fmomei =     devVme(distance->psq2)*rhomij*vf1;
      double fmomej =     devVme(distance->psq1)*rhomji*vf2;
      force[i]  += -fmomdi*distance->dr2ijri - fmomdj*distance->dr2jiri;
      neighbor.box->addForce(j,-fmomdi*distance->dr2ijrj + fmomdj*distance->dr2jirj);
      forcer[i] +=  fmomei*distance->dp2ijpi + fmomej*distance->dp2jipi
  	          + fmomdi*distance->dr2ijpi + fmomdj*distance->dr2jipi;
      neighbor.box->addForceR(j,fmomei*distance->dp2ijpj + fmomej*distance->dp2jipj
	                      + fmomdi*distance->dr2ijpj + fmomdj*distance->dr2jipj);


      // Derivative of p_j/p^0_j in the vector potential.
      if(optDerivative) {
        double facm1 = devVmd(psq1)*rhomji;
        double facm2 = devVmd(psq2)*rhomij;
        distance->devV(pk1/pk1[0],pk2/pk2[0]);
        forcer[i] += facm1*distance->devV1;
        neighbor.box->addForceR(j,facm2*distance->devV2);

        double facmom = facm1*vf1 + facm2*vf2;
        forcer[i] += distance->devgam1*facmom;
        neighbor.box->addForceR(j,distance->devgam2*facmom);


      }

    }
  }

}

// non-relativistic two-body distance
// (distance at the global cm frame where the total momentum is zero)
void MBoxV::qmdMatrixNeighbor1(NeighborV& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(int i=0; i< NV; i++) {
    Vec4 r1  = part[i]->getR();
    Vec4 p1  = part[i]->getP();
    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    Vec4 fri = optBaryonCurrent ? forcer[i] : 0.0;
    int bi   = part[i]->baryon()/3;
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;

    for(int j=0; j< n2; j++) {
      Vec4 r2  = part2[j]->getR();
      Vec4 p2  = part2[j]->getP();
      Vec4 frj = optBaryonCurrent ? forcer[j] : 0.0;
      int bj = part2[j]->baryon()/3;
      if(optPotentialArg==1) p2 = part2[j]->getPcan(optVdot);
      double qfac2 = part2[j]->getTf() > globalTime ? part2[j]->qFactor() : 1.0;

      //double den = facG * exp(-(r1 - r2).pAbs2()*wG);
      double den = gamCM*facG* exp(-wG*( (r1 - r2).pT2() + pow2(gamCM*(r1[3]-r2[3])) ) )*qfac1*qfac2;
      neighbor.setRhom(i,j, den);
      neighbor.neighbor->setRhom(j,i, den);

      JB[i] += den*(p2/p2[0]+frj)*bj;
      neighbor.box->addJB(j,den*(p1/p1[0]+fri)*bi);

      if(!withMomDep) continue;
      double ps = (p1 - p2).pAbs2() - (1-optPV) * (p1-p2)[0]*(p1-p2)[0];
      double  pmom2ij = den*(vex1/(1.0+ps/pmu1)+vex2/(1.0+ps/pmu2));
      vmom4[i] += pmom2ij*p2/p2[0];
      neighbor.box->addVmom4(j, pmom2ij*p1/p1[0]);
    }
  }

}

// two-body distance at their center-of-momentum frame.
void MBoxV::qmdMatrixNeighbor2(NeighborV& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(int i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    int bi = part[i]->baryon()/3;
    Vec4 fri=0.0;
    if(optBaryonCurrent) fri=forcer[i];
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;

    for(int j=0; j< n2; j++) {
      Vec4 frj=0.0;
      if(optBaryonCurrent) frj=forcer[j];
      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      if(optPotentialArg==1) p2 = part2[j]->getPcan(optVdot);
      int bj = part2[j]->baryon()/3;
      double qfac2 = part2[j]->getTf() > globalTime ? part2[j]->qFactor() : 1.0;
      Vec4 dr = r1 - r2; dr[0] = 0.0;  // just in case;
      Vec4 dp = p1 - p2;
      Vec4 P  = p1 + p2;
      double s = P.m2Calc();
      double dot4 = dr * P;
      double drcmsq = dr.m2Calc() - dot4*dot4/s;
      double g12 = 1.0;
      if(optScalarDensity > 0)  g12 = P[0]/sqrt(s);
      double den = g12 * facG * exp(drcmsq*wG)*qfac1*qfac2;
      //rhom[i][j] = den;
      //rhom[j][i] = den;
      neighbor.setRhom(i,j,den);
      neighbor.neighbor->setRhom(j,i,den);
      JB[i] += den*(p2/p2[0] + frj)*bj;
      neighbor.box->addJB(j,den*(p1/p1[0] + fri)*bi);

      if(!withMomDep) continue;
      dot4 = dp * P;
      double ps = dp.m2Calc() - optPV * dot4*dot4/s;
      double  pmom2ij = den*(vex1/(1.0-ps/pmu1)+vex2/(1.0-ps/pmu2));
      vmom4[i] += pmom2ij*p2/p2[0];
      neighbor.box->addVmom4(j, pmom2ij*p1/p1[0]);
    }
  }

}

// two-body distance measured from the rest-frame of particle j.
void MBoxV::qmdMatrixNeighbor3(NeighborV& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(int i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    Vec4 v1 = p1/p1.mCalc();
    int bi = part[i]->baryon()/3;
    Vec4 fri=0.0;
    if(optBaryonCurrent) fri=forcer[i];
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;
    for(int j=0; j< n2; j++) {
      Vec4 frj=0.0;
      if(optBaryonCurrent) frj=forcer[j];
      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      if(optPotentialArg == 1) p2 = part2[j]->getPcan(optVdot);
      Vec4 v2 = p2/p2.mCalc();
      int bj = part2[j]->baryon()/3;
      double qfac2 = part2[j]->getTf() > globalTime ? part2[j]->qFactor() : 1.0;

      Vec4 dr = r1 - r2; dr[0] = 0.0;  // just in case;
      double dot4 = dr * v2;
      double drsq1 = dr.m2Calc() - dot4*dot4;
      double rhomij = facG * exp(drsq1*wG)*qfac1*qfac2;
      neighbor.setRhom(i,j, rhomij);

      dot4 = dr * v1;
      double drsq2 = dr.m2Calc() - dot4*dot4;
      double rhomji= facG * exp(drsq2*wG)*qfac1*qfac2;
      neighbor.neighbor->setRhom(j,i, rhomji);

      JB[i] += rhomij*(v2 + frj)*bj;
      neighbor.box->addJB(j, rhomji*(v1 + fri)*bi);

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
        ps1 = dpsq - optPV * pow2(dp * v2);
        ps2 = dpsq - optPV * pow2(dp * v1);
      }
      }

      double  pmom2ij = vex1/(1.0-ps1/pmu1)+vex2/(1.0-ps1/pmu2);
      double  pmom2ji = vex1/(1.0-ps2/pmu1)+vex2/(1.0-ps2/pmu2);

      vmom4[i] += pmom2ij*rhomij*v2;
      neighbor.box->addVmom4(j,pmom2ji*rhomji*v1);
    }
  }

}

void MBoxV::computeForceNeighbor1(NeighborV& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = p1;
    if(optVdot <= 1) pk1 = part[i]->getPkin();
    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    Vec4 v1 = p1/p1[0];
    int  bar1 = part[i]->baryon()/3;

    for(auto j=0; j< n2; j++) {
      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      Vec4 pk2 = p2;
      if(optVdot <= 1) pk2 = part2[j]->getPkin();
      Vec4 A1 = RQMDv::facV(i,pk2);
      Vec4 A2 = facV(j,pk1,neighbor);
      if(optPotentialArg==1) p2 = part2[j]->getPcan(optVdot);
      Vec4 v2 = p2/p2[0];
      int  bar2 = part2[j]->baryon()/3;

      double rhomij=neighbor.rhom[i][j];
      double rhomji=neighbor.neighbor->rhom[j][i];
      double fengi = (A1 * v1);
      double fengj = (A2 * v2);
      double fsky = -wG*(fengi*rhomij + fengj*rhomji)*bar1*bar2;

      Vec4 dr = r1 - r2;
      dr[3] *= gamCM2;

      force[i] += -2*fsky*dr;
      neighbor.box->addForce(j, 2*fsky*dr);

      // Derivative of p^\mu_j/p^0_j in the vector density.
      if(optDerivative) {
        Vec4 Ai = RQMDv::facV(i,pk1);
        Vec4 Aj = facV(j,pk2,neighbor);
        double fsgam2i = optP0dev*rhomji*dot3(Aj,v1)/p1[0];
        double fsgam2j = optP0dev*rhomij*dot3(Ai,v2)/p2[0];
        double fsgam3i = -rhomji/p1[0];
        double fsgam3j = -rhomij/p2[0];
        forcer[i] += (fsgam2i*v1 + Aj*fsgam3i)*bar1*bar2;
        neighbor.box->addForceR(j, (fsgam2j*v2 + Ai*fsgam3j)*bar1*bar2);
      }

      if(!withMomDep) continue;

      double vf1=1.0;
      double vf2=1.0;
      if(optVectorPotential==1) {
        vf1 = pk1/pk1[0] * v2;
        vf2 = pk2/pk2[0] * v1;
      }

      Vec4 dp = p1 - p2;
      double psq = dp.pAbs2() - (1-optPV) * dp[0]*dp[0];
      double fac1 = 1.0 + psq/pmu1;
      double fac2 = 1.0 + psq/pmu2;
      double fengij = vf1*rhomij + vf2*rhomji;

      // p derivative term.
      double fmome = -fengij*(vex1/(pmu1*fac1*fac1)+vex2/(pmu2*fac2*fac2));

      // r derivative term.
      double fmomd = -wG*fengij*(vex1/fac1 + vex2/fac2);

      Vec4 dp2pi =  dp -optP0dev*(1-optPV)*dp[0]*p1/p1[0];
      Vec4 dp2pj = -dp -optP0dev*(1-optPV)*dp[0]*p2/p2[0];

      force[i]  += -2*fmomd*dr;
      neighbor.box->addForce(j,2*fmomd*dr);
      //forcer[i] +=  2*fmome*(p1 - p2);
      //forcer[j] +=  2*fmome*(p2 - p1);
      forcer[i] +=  2*fmome*dp2pi;
      neighbor.box->addForceR(j,2*fmome*dp2pj);

      // Derivative of p_j/p^0_j in the vector potential.
      if(optDerivative) {
        double fm1 = (vex1/fac1+vex2/fac2)*rhomji/p1[0];
        double fm2 = (vex1/fac1+vex2/fac2)*rhomij/p2[0];
        vf1 = dot3(pk1/pk1[0], v2);
        vf2 = dot3(pk2/pk2[0], v1);
        forcer[i] += -fm1 * pk2/pk2[0] + fm1*vf2 * v1 * optP0dev;
        neighbor.box->addForceR(j,-fm2 * pk1/pk1[0] + fm2*vf1 * v2 * optP0dev);
      }

    }
  }

}

// two-body distance is defined by the c.m. frame of two-particles.
void MBoxV::computeForceNeighbor2(NeighborV& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = p1;
    if(optVdot <= 1) pk1 = part[i]->getPkin();

    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    Vec4 bi = p1/p1[0];
    int  b1 = part[i]->baryon()/3;
    double emi = part[i]->getMass();

    for(auto j=0; j< n2; j++) {
      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      Vec4 pk2 = p2;
      if(optVdot <= 1) pk2 = part2[j]->getPkin();
      Vec4 A1 = RQMDv::facV(i,pk2);
      Vec4 A2 = facV(j,pk1,neighbor);

      if(optPotentialArg==1) p2 = part2[j]->getPcan(optVdot);
      Vec4 bj = p2/p2[0];
      int  b2 = part[i]->baryon()/3;
      double emj = part2[j]->getMass();

      double rhomij=neighbor.rhom[i][j];
      double rhomji=neighbor.neighbor->rhom[j][i];

      double fengi = (A1 * bi)*b1*b2;
      double fengj = (A2 * bj)*b2*b1;
      double fsky3 = fengi*rhomij;
      double fsky4 = fengj*rhomji;
      double fsky = -wG*( fsky3 + fsky4 );

      Vec4 pcm = p1 + p2;
      Vec4 bbi = pcm/pcm[0] - optP0dev*bi;
      Vec4 bbj = pcm/pcm[0] - optP0dev*bj;

      Vec4 dr  = r1 - r2;dr[0]=0.0;
      double s = pcm.m2Calc();
      double rbij = -dr * pcm/s;
      Vec4 dr2ri = dr + rbij*pcm;
      Vec4 dr2rj = -dr2ri;
      Vec4 dr2pi = (dr + pcm[0]*rbij*bbi)*rbij;
      Vec4 dr2pj = (dr + pcm[0]*rbij*bbj)*rbij;

      force[i]  += -2*fsky*dr2ri;
      neighbor.box->addForce(j,-2*fsky*dr2rj);
      forcer[i] +=  2*fsky*dr2pi;
      neighbor.box->addForceR(j,2*fsky*dr2pj);

      // Derivative of gamma_ij in front of Gaussian.
      //if(optScalarDensity > 0) {
      if(optDerivative) {
        double facsk = fsky3 + fsky4;
	double fsgam1 = facsk / s;
	double fsgam2 = optP0dev*facsk*(1.0/pcm[0] - pcm[0]/s);
	forcer[i] += fsgam1*pcm + fsgam2*bi;
	neighbor.box->addForceR(j,fsgam1*pcm + fsgam2*bj);
      }


        // Derivative of p^\mu_j/p^0_j in the vector density.
      if(optDerivative) {
        Vec4 Ai = RQMDv::facV(i,pk1);
        Vec4 Aj = facV(j,pk2,neighbor);
        double fsgam2i = optP0dev*rhomji*dot3(A2,bi)/p1[0];
        double fsgam2j = optP0dev*rhomij*dot3(A1,bj)/p2[0];
        double fsgam3i = -rhomji/p1[0];
        double fsgam3j = -rhomij/p2[0];
        forcer[i] += (fsgam2i*bi + Aj*fsgam3i)*b1*b2;
	neighbor.box->addForceR(j,(fsgam2j*bj + Ai*fsgam3j)*b1*b2);
      }


      if(!withMomDep) continue;

      double vf1=1.0;
      double vf2=1.0;
      if(optVectorPotential==1) {
        vf1 = pk1/pk1[0] * bj;
        vf2 = pk2/pk2[0] * bi;
      }

      Vec4 dp = p1 - p2;
      double psq = dp.m2Calc() - optPV * pow2(dp * pcm)/s;
      double fac1 = 1.0 - psq/pmu1;
      double fac2 = 1.0 - psq/pmu2;

      // p derivative term.
      double fengij = vf1*rhomij + vf2*rhomji;
      double fmome = -fengij*(vex1/(pmu1*fac1*fac1)+vex2/(pmu2*fac2*fac2));

      // r derivative term.
      double facmom = vex1/fac1 + vex2/fac2;
      double fmomd = -wG*fengij*facmom;
      double pma = pow2((emi*emi - emj*emj)/s);
      Vec4 dp2pi =  dp - optP0dev*dp[0]*bi + optPV * pcm[0]*pma*bbi;
      Vec4 dp2pj = -dp + optP0dev*dp[0]*bj + optPV * pcm[0]*pma*bbj;

      force[i]  += -2*fmomd*dr2ri;
      neighbor.box->addForce(j, -2*fmomd*dr2rj);
      forcer[i] +=  2*fmomd*dr2pi + 2*fmome*dp2pi;
      neighbor.box->addForceR(j, 2*fmomd*dr2pj + 2*fmome*dp2pj);

      //if(optScalarDensity > 0 ) {
      if(optDerivative) {
	// derivative of gamma_ij
	double fmgam1 = fengij*facmom/s;
	double fmgam2 = optP0dev*fengij*facmom*(1.0/pcm[0] - pcm[0]/s);
        forcer[i] += pcm*fmgam1 + bi*fmgam2;
        neighbor.box->addForceR(j, pcm*fmgam1 + bj*fmgam2);
      }

        // Derivative of p_i/p^0_i in the vector density.
      if(optDerivative) {
        double fm1 = -(vex1/fac1+vex2/fac2)*rhomji/p1[0];
        double fm2 = -(vex1/fac1+vex2/fac2)*rhomij/p2[0];
        vf1 = dot3(pk1/pk1[0], bj);
        vf2 = dot3(pk2/pk2[0], bi);
        forcer[i] += fm1 * pk2/pk2[0] - fm1*vf2 * bi * optP0dev;
        neighbor.box->addForceR(j,fm2 * pk1/pk1[0] - fm2*vf1 * bj * optP0dev);
      }

    }
  }


}

// two-body distance is defined by the rest frame of particle.
void MBoxV::computeForceNeighbor3(NeighborV& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;


  for(int i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = optVdot <= 1 ? part[i]->getPkin() : p1;

    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    double emi = part[i]->getMass();
    Vec4 vi = p1/emi;
    int bar1= part[i]->baryon()/3;

    for(int j=0; j< n2; j++) {
      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      Vec4 pk2 = optVdot <= 1 ? part2[j]->getPkin() : p2;

      if(optPotentialArg==1) p2 = part2[j]->getPcan(optVdot);
      double emj = part2[j]->getMass();
      Vec4 vj = p2/emj;
      int bar2 = part2[j]->baryon()/3;

      double rhomij=neighbor.rhom[i][j];
      double rhomji=neighbor.neighbor->rhom[j][i];

      Vec4 Ai = RQMDv::facV(i,pk2);
      Vec4 Aj = facV(j,pk1,neighbor);
      double fskyi = -wG*(Ai * vi)*rhomij*bar1*bar2; 
      double fskyj = -wG*(Aj * vj)*rhomji*bar1*bar2; 

      Vec4 dr = r1 - r2; dr[0]=0.0;
      double rbi = dot3(dr,vi);
      double rbj = dot3(dr,vj);
      Vec4 dr2ri =  2*(dr + rbj*vj);  // R^2_{ij}/dp_i
      Vec4 dr2rj =  2*(dr + rbi*vi);  // R^2_{ji}/dp_i
      Vec4 dr2pi =  2*dr*rbi/emi;     // R^2_{ji}/dp_i
      Vec4 dr2pj =  2*dr*rbj/emj;     // R^2_{ij}/dp_j

      force[i]  += -fskyi*dr2ri - fskyj*dr2rj;
      neighbor.box->addForce(j, fskyi*dr2ri + fskyj*dr2rj);
      forcer[i] +=  fskyj*dr2pi;
      neighbor.box->addForceR(j,fskyi*dr2pj);

      // Derivative of p^\mu/m_j term in the vector potential.
      if(optDerivative) {
        Vec4 A1 = RQMDv::facV(i,pk1)*rhomji/emi*bar1*bar2;
        Vec4 A2 = facV(j,pk2,neighbor)*rhomij/emj*bar1*bar2;
        forcer[i] += (optP0dev*A2[0]*p1/p1[0] - A2);
        neighbor.box->addForceR(j, optP0dev*A1[0]*p2/p2[0] - A1);

        //forcer[i] += (optP0dev*Aj[0]*vi/p1[0] - Aj/emi)*rhomji;
        //neighbor.box->addForceR(j, (optP0dev*Ai[0]*vj/p2[0] - Ai/emj)*rhomij);
      }


      if(!withMomDep) continue;

      double vf1=vj[0];
      double vf2=vi[0];
      if(optVectorPotential==1) {
        vf1 = pk1/pk1[0] * vj;
        vf2 = pk2/pk2[0] * vi;
      }

      if(optMDarg3==1) distanceP1(p1,p2);
      else if(optMDarg3==2) distanceP2(p1,p2);
      else distanceP3(p1,p2);

      double fmomdi = -wG*devVmd(psq2)*rhomij*vf1;
      double fmomdj = -wG*devVmd(psq1)*rhomji*vf2;
      double fmomei =     devVme(psq2)*rhomij*vf1;
      double fmomej =     devVme(psq1)*rhomji*vf2;

      force[i]  += -fmomdi*dr2ri   - fmomdj*dr2rj;
      forcer[i] +=  fmomei*dp2ijpi + fmomej*dp2jipi + fmomdj*dr2pi;
      neighbor.box->addForce(j, fmomdi*dr2ri   + fmomdj*dr2rj);
      neighbor.box->addForceR(j, fmomej*dp2jipj + fmomei*dp2ijpj + fmomdi*dr2pj);

      // Derivative of p_j / m_j part.
      if(optDerivative) {
        double facm1 = devVmd(psq1)*rhomji/emi;
        double facm2 = devVmd(psq2)*rhomij/emj;
        forcer[i] += facm1*( optP0dev*p1/p1[0] - pk2/pk2[0] );
        neighbor.box->addForceR(j,facm2*( optP0dev*p2/p2[0] - pk1/pk1[0] ));
      }

    }
  }

}

// Pre-factor from the baryon current.
Vec4 MBoxV::facV(int i, Vec4& pkin, NeighborV& neighbor)
{
  double rh = neighbor.box->getRho(i);
  if(rh<1e-15) return Vec4(0.0);
  Vec4 bi=0.0;

  // 4-components of the vector potential is fully included.
  if(optVectorPotential==1) {
    bi = pkin/pkin[0];

  // Only time-component V^0 term is included.
  } else if(optVectorPotential==2) {
    bi[0]=1.0;

  // Only time-component V^0 term is included with the form of V(rho_B).
  } else if (optVectorPotential==3) {
    return (t1 + t3f*neighbor.box->getRhog(i))/rh * neighbor.box->getJB(i);
  }

  double vj = neighbor.box->getJB(i) * bi;
  double vv = t1 + t3*neighbor.box->getRhog(i);  // V/rho_B
  double dv = (gam-1.0)*t3*pow(rh,gam-3);  // del(V/rho)/rho

  return vj*dv*neighbor.box->getJB(i) + vv*bi;

}

//...Compute time-derivatives of the vector potential.
void MBoxV::computeVdot(NeighborV& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  bool opt=true;
  //opt=false;

  Vdot.assign(NV,0.0);

  for(int i=0;i<NV; i++) {
    double vvi = t1 + t3*rhog[i];  // V_i/rho_i
    double dvi=0.0;
    // del(V_i/rho_i)/rho_i
    if(abs(rho[i]) > 1e-7) dvi = (gam-1.0)*t3*pow(rho[i],gam-3.0);

    Vec4 Bi = dvi*JB[i];
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 vk1 = p1/p1[0];
    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    double emi=p1[0];
    if(optTwoBodyDistance==3)  emi = p1.mCalc();
    Vec4 vi = p1 / emi;
    Vec4 b1 = p1 / p1[0];

    for(int j=0;j<n2; j++) {

      double vvj = t1 + t3*neighbor.box->getRhog(j);
      double dvj=0.0;
      if(abs(neighbor.box->getRho(j)) > 1e-7) dvj = (gam-1.0)*t3*pow(neighbor.box->getRho(j),gam-3.0);
      Vec4 jb = neighbor.box->getJB(j);
      Vec4 Bj = dvj*jb;
      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      Vec4 vk2 = p2/p2[0];
      if(optPotentialArg==1) p2 = part2[j]->getPcan(optVdot);
      double emj=p2[0];
      if(optTwoBodyDistance==3)  emj = p2.mCalc();
      Vec4 vj = p2 / emj;
      Vec4 b2 = p2 / p2[0];
   
      Vec4 Aj = (vj * JB[i]) * Bi;
      Vec4 Ai = (vi * jb   ) * Bj;

      // Compute derivatives

      // Non-relativistic.
      Vec4 rk = r1 - r2;
      Vec4 dgami=0.0;
      Vec4 dgamj=0.0;
      Vec4 dr2ri   = 2 * rk;
      Vec4 dr2rj   = -dr2ri;
      Vec4 drji2ri = dr2ri;
      Vec4 drji2rj = dr2rj;
      Vec4 dr2pi = 0.0;
      Vec4 dr2pj = 0.0;
      Vec4 drji2pi=0.0;
      Vec4 drji2pj=0.0;

      // two-body c.m. frame.
      if(optTwoBodyDistance==2) {

      Vec4 pcm = p1 + p2;
      Vec4 bet = pcm/pcm[0];
      Vec4 bbi = bet - optP0dev*b1;
      Vec4 bbj = bet - optP0dev*b2;
      double s = pcm.m2Calc();
      double rbij = dot3(rk,pcm)/s;
      dr2ri =  2*(rk + rbij*pcm);
      dr2rj = -dr2ri;
      dr2pi = 2*(rk + pcm[0]*rbij*bbi)*rbij;
      dr2pj = 2*(rk + pcm[0]*rbij*bbj)*rbij;
      drji2ri =  dr2ri;
      drji2rj = -drji2ri;
      drji2pi = dr2pi;
      drji2pj = dr2pj;

      // derivatives from gamma_{ij}.
      dgami = optP0dev*(1.0/pcm[0]-pcm[0]/s)*b1+pcm/s;
      dgamj = optP0dev*(1.0/pcm[0]-pcm[0]/s)*b2+pcm/s;

      // rest frame of particle i or j.
      } else if(optTwoBodyDistance == 3) {

        double rbi=dot3(rk,vi);
        double rbj=dot3(rk,vj);
        dr2ri =  2*(rk+rbj*vj);    // dR~^2_ij/dR_i
        dr2rj =  -dr2ri;           // dR~^2_ij/dR_j
        dr2pi =  0.0;              // dR~^2_ij/dP_i
        dr2pj =  2*rk*rbj/emj;     // dR~^2_ij/dP_j

        drji2rj = 2*(rk+rbi*vi);
        drji2rj =  -drji2ri;
        drji2pi =  2*rk*rbi/emi;
        drji2pj =  0.0;
      }

      Vec4 forcei=force[i];
      Vec4 forceri=forcer[i];
      Vec4 forcej=neighbor.box->getForce(j);
      Vec4 forcerj=neighbor.box->getForceR(j);
      double rhomij=neighbor.rhom[i][j];
      double rhomji=neighbor.neighbor->rhom[j][i];

      double xdoti=dot3(vk1+forceri,dr2ri)   + dot3(vk2+forcerj,dr2rj);
      double xdotj=dot3(vk2+forcerj,drji2rj) + dot3(vk1+forceri,drji2ri);

      double pdoti=dot3(forcei,dr2pi-dgami/wG)
	         + dot3(forcej,dr2pj-dgamj/wG);
      double pdotj=dot3(forcej,drji2pj-dgamj/wG)
	         + dot3(forcei,drji2pi-dgami/wG);

      double doti=-wG*(xdoti+pdoti)*rhomij;
      double dotj=-wG*(xdotj+pdotj)*rhomji;
      Vdot[i] += doti*(Aj + vvi*vj);
      neighbor.box->addVdot(j,dotj*(Ai + vvj*vi));

      // Momentum dependent potential.
      double fmomdi=0.0, fmomdj=0.0;
      if(withMomDep) {

      if(optTwoBodyDistance==1) {
        distanceP1(p1,p2);
      } else if(optTwoBodyDistance==2) {
        distanceP2(p1,p2);
      } else if(optTwoBodyDistance==3) {
        if(optMDarg3==1) distanceP1(p1,p2);
        else if(optMDarg3==2) distanceP2(p1,p2);
        else distanceP3(p1,p2);
      }
      fmomdi = devVmd(psq2);
      fmomdj = devVmd(psq1);
      Vdot[i] += doti*fmomdi*vj;
      neighbor.box->addVdot(j, dotj*fmomdj*vi);

      pdoti=dot3(forcei,dp2ijpi) + dot3(forcej,dp2ijpj);
      pdotj=dot3(forcej,dp2jipj) + dot3(forcei,dp2jipi);

      Vdot[i] += pdoti*devVme(psq2)*rhomij*vj;
      neighbor.box->addVdot(j, pdotj*devVme(psq1)*rhomji*vi);

      }


      if(opt) continue;

      //Compute derivatives of pre-factors v_j^\mu in the vector potential.
      double fai = - dot3(forcej,Bi/emj) * rhomij;
      double faj = - dot3(forcei,Bj/emi) * rhomji;
      double fai2 = dot3(forcej,vj);
      double faj2 = dot3(forcei,vi); 
      double fai3 = dot3(vj,Bi/p2[0]);
      double faj3 = dot3(vi,Bj/p1[0]);

      double  fvi2=0.0;
      double  fvj2=0.0;
      if(optTwoBodyDistance==3) {
        fai3=Bi[0]/p2[0];
        faj3=Bj[0]/p1[0];
      } else {
        fvi2=-optP0dev*fai2*rhomij*vvi/p2[0];
        fvj2=-optP0dev*faj2*rhomji*vvj/p1[0];
      }

      fai += optP0dev*fai2*fai3*rhomij;
      faj += optP0dev*faj2*faj3*rhomji;

      double fvi = vvi*rhomij/emj;
      double fvj = vvj*rhomji/emi;

      Vdot[i] += fai*JB[i] + fvi*forcej + fvi2*vj;
      neighbor.box->addVdot(j, faj*neighbor.box->getJB(j) + fvj*forcei + fvj2*vi);

      // Momentum dependent potential.
      if(withMomDep) {
        Vdot[i] += fmomdi*rhomij*forcej/emj;
        neighbor.box->addVdot(j, fmomdj*rhomji*forcei/emi);
        if(optTwoBodyDistance != 3) {
          Vdot[i] += -optP0dev*fai2*rhomij*fmomdi/p2[0]*vj;
          neighbor.box->addVdot(i, -optP0dev*faj2*rhomji*fmomdj/p1[0]*vi);
	}
      }


    } // end j loop
  }   // end i loop

  //return;

  Vec4 vdott=0.0;
  for(int i=0;i<NV;i++) {
    vdott += Vdot[i];
  }

  Vec4 dv = vdott/NV;
  for(int i=0;i<NV;i++) {
    Vdot[i] -= dv;
  }

}


} // end of namespace jam2
