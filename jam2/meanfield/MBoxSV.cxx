#include <jam2/meanfield/MBoxSV.h>

namespace jam2 {

void NeighborSV::resize(int n) 
{
  int m = box->particleSize();
  rhom.resize(n);
  for(int i=0;i<n;i++) rhom[i].resize(m);
}

MBoxSV:: ~MBoxSV()
{
  for(auto& i : neighbors) delete i;
  neighbors.clear();
  clearMatrix();
}

bool MBoxSV::haveThisSite(MBoxSV* box) const
{
  auto first = neighbors.begin();
  auto last = neighbors.end();

  while(first != last) {
    if((*first)->box == box) return true;
    ++first;
  }
  return false;
}

NeighborSV* MBoxSV::findNeighbor(MBoxSV* box)
{
  auto&& first = neighbors.begin();
  auto&& last = neighbors.end();

  while(first != last) {
    if((*first)->box == box) return *first;
    ++first;
  }
  cout << "MBoxSV error no such neighbor "<<endl;
  exit(1);
}

void MBoxSV::clearMatrix()
{
  part.clear();
  rhos.clear();
  rho.clear();
  rhog.clear();
  rhog2.clear();
  rhosg.clear();
  vmoms.clear();
  vmom4.clear();
  force.clear();
  forcer.clear();
  Vdot.clear();
  JB.clear();
  for(int i=0;i<NV;i++) rhom[i].clear();
  rhom.clear();
}

void MBoxSV::initBox()
{
  NV = part.size();
  if(NV==0) return;

  rhog.resize(NV);
  rhog2.resize(NV);
  rhosg.resize(NV);
  rhom.resize(NV);
  for(int i=0;i<NV;i++) {
    rhom[i].resize(NV);
  }
  for(auto& i : part) 
    potential->setPotentialParam(i);


  rhos.assign(NV,0.0);
  rho.assign(NV,0.0);
  vmoms.assign(NV,0.0);
  vmom4.assign(NV,0.0);
  JB.assign(NV,0.0);
  Vdot.assign(NV,0.0);

  force.assign(NV,0.0);
  forcer.assign(NV,0.0);

  // free mass is used in the calculation of potential and force.
  if(optVdot==1) {
    for(auto& i : part) i->setFree();
  }
  if(optPotentialArg >= 1) {
    for(auto& i : part) i->setPotS(0.0);
  }

  for(auto& neighbor : neighbors) {
    neighbor->resize(NV);
  }

}

void MBoxSV::qmdMatrix()
{
  if(NV==0) return;

  RQMDsv::qmdMatrix();
  for(auto& neighbor : neighbors) {
    if(neighbor->getAction()) qmdMatrixNeighbor(*neighbor);
  }

}

void MBoxSV::computeForce(double dt)
{
  if(NV==0) return;

  if(transportModel==1) {
    RQMDsv::computeForce();
    for(auto& neighbor : neighbors) {
      if(neighbor->getAction()) computeForceNeighbor(*neighbor);
    }
  } else {
    RQMDsv::computeBUUForce();
    for(auto& neighbor : neighbors) {
      if(neighbor->getAction()) computeBUUForceNeighbor(*neighbor);
    }
  }

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

    RQMDsv::computeVdot();
    for(int i=0; i< NV; i++) {
      part[i]->updateByForce(force[i] - Vdot[i], forcer[i],dt);
    }
  }

  if(optPotentialArg >= 1) setScalarPotential();

}

void MBoxSV::qmdMatrixNeighbor(NeighborSV& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(int i=0; i< NV; i++) {

    Vec4 r1  = part[i]->getR();
    Vec4 p1  = part[i]->getP();
    if(optPotentialArg>=1) p1 = part[i]->getPcan(optVdot);
    distance->setRP1(r1,p1);
    distanceP->setRP1(r1,p1);
    Vec4 v1 = p1/p1[0];
    double fi=p1.mCalc()/p1[0];

    Vec4 fri = optBaryonCurrent ? part[i]->forceR() : 0.0;
    int bi   = part[i]->baryon()/3;
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;
    if(bi==0) {
      qfac1 = facMesonPot;
      bi=1;
    }
    double* gfac1= part[i]->facPotential();
    int potid1=part[i]->potentialID();

    for(int j=0; j< n2; j++) {

      int potid2=part2[j]->potentialID();
      double pot1 = potential->matrix(potid1,potid2);
      double pot2 = potential->matrix(potid2,potid1);

      Vec4 r2  = part2[j]->getR();
      Vec4 p2  = part2[j]->getP();
      if(optPotentialArg>=1) p2 = part2[j]->getPcan(optVdot);
      distance->setRP2(r2,p2);
      distanceP->setRP2(r2,p2);
      Vec4 v2 = p2/p2[0];
      double fj=p2.mCalc()/p2[0];

      Vec4 frj = optBaryonCurrent ? part2[j]->forceR() : 0.0;
      int bj = part2[j]->baryon()/3;
      double qfac2 = part2[j]->getTf() > globalTime ? part2[j]->qFactor() : 1.0;
      if(bj==0) {
	qfac2 = facMesonPot;
	bj=1;
      }
      double* gfac2= part2[j]->facPotential();

      distance->density();
      double den1 = distance->density1*qfac1*qfac2/overSample*pot1;
      double den2 = distance->density2*qfac1*qfac2/overSample*pot2;

      neighbor.setRhom(i,j,den1);
      neighbor.neighbor->setRhom(j,i, den2);

      // scalar density
      rhos[i] += den1*fj;
      neighbor.box->addRhos(j,den2*fi);

      JB[i] += den1*(v2+frj)*bj;
      neighbor.box->addJB(j,den2*(v1+fri)*bi);

      if(!withMomDep) continue;

      distanceP->psq();
      //potential->Vmd(distanceP->psq1,distanceP->psq2);
      potential->Vmd(distanceP->psq1,distanceP->psq2,gfac1[6]*gfac2[6],gfac1[7]*gfac2[7],
	  gfac1[8],gfac1[9],gfac2[8],gfac2[9]);

      vmoms[i] += potential->vmomsi*den1*fj;
      neighbor.box->addVmoms(j, potential->vmomsj*den2*fi);

      vmom4[i] += potential->vmom4i*den1*v2;
      neighbor.box->addVmom4(j, potential->vmom4j*den2*v1);

    }
  }

}

void MBoxSV::computeForceNeighbor(NeighborSV& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = optVdot <= 1 ? part[i]->getPkin() : p1;
    if(optPotentialArg>=1) p1 = part[i]->getPcan(optVdot);
    int  bar1 = part[i]->baryon()/3;

    double fi = p1.mCalc()/p1[0];
    double fengi=fi;
    if(optPotentialArg==3) {
      double meff1 = part[i]->getEffectiveMass();
      fengi = meff1/sqrt(meff1*meff1 + pk1.pAbs2());
    }

    distance->setRP1(r1,p1);
    distanceP->setRP1(r1,p1);
    double* gf1=part[i]->facPotential();
    potential->set1(p1,pk1,JB[i],fi,fengi,rho[i],rhog[i],rhog2[i],
                    part[i]->pots()-vmoms[i],rhos[i],rhosg[i],bar1,gf1[1],gf1[2],gf1[3],gf1[4],gf1[5]);

    for(auto j=0; j< n2; j++) {

      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      Vec4 pk2 = optVdot <= 1 ? part2[j]->getPkin(): p2;
      if(optPotentialArg>=1) p2 = part2[j]->getPcan(optVdot);
      int  bar2 = part2[j]->baryon()/3;

      double fj=p2.mCalc()/p2[0];
      double fengj=fj;
      if(optPotentialArg==3) {
        double meff2 = part[j]->getEffectiveMass();
        fengj = meff2/sqrt(meff2*meff2 + pk2.pAbs2());
      }

      distance->setRP2(r2,p2);
      distanceP->setRP2(r2,p2);
      double* gf2=part2[j]->facPotential();
      potential->set2(p2,pk2,neighbor.box->getJB(j),fj,fengj,neighbor.box->getRho(j),
	  neighbor.box->getRhog(j),
	  neighbor.box->getRhog2(j),
	  part2[j]->pots()-neighbor.box->getVmoms(j),
	  neighbor.box->getRhos(j),neighbor.box->getRhosg(j),
	  bar2,gf2[1],gf2[2],gf2[3],gf2[4],gf2[5]);

      distance->distanceR();

      double rhomij=neighbor.rhom[i][j];
      double rhomji=neighbor.neighbor->rhom[j][i];
      potential->dVdns(rhomij,rhomji);

      double fsky1=potential->fskyi;
      double fsky2=potential->fskyj;

      force[i]  += -fsky1*distance->dr2ijri - fsky2*distance->dr2jiri;
      forcer[i] +=  fsky1*distance->dr2ijpi + fsky2*distance->dr2jipi;
      neighbor.box->addForce(j, -fsky1*distance->dr2ijrj - fsky2*distance->dr2jirj);
      neighbor.box->addForceR(j, fsky1*distance->dr2ijpj + fsky2*distance->dr2jipj);

      if(optDerivative) {
        potential->dVdp();
        forcer[i] +=  potential->forceri;
        neighbor.box->addForceR(j, potential->forcerj);
        // gamma derivative.
        distance->devGamma();
        forcer[i] += distance->devgam1*potential->facsk;
        neighbor.box->addForceR(j,  distance->devgam2*potential->facsk);
      }

      if(!withMomDep) continue;

      distanceP->distanceP();
      //potential->dVdmd(distanceP->psq1,distanceP->psq2,gf1[6],gf2[6],gf1[7],gf2[7]);
      potential->dVdmd(distanceP->psq1,distanceP->psq2,gf1[6]*gf2[6],gf1[7]*gf2[7],gf1[8],gf2[8],gf1[9],gf2[9]);
      double fmomdi=potential->fmomdi;
      double fmomdj=potential->fmomdj;
      double fmomei=potential->fmomei;
      double fmomej=potential->fmomej;

      force[i]  +=              -fmomdi*distance->dr2ijri   - fmomdj*distance->dr2jiri;
      neighbor.box->addForce(j, -fmomdi*distance->dr2ijrj   - fmomdj*distance->dr2jirj);

      forcer[i] +=  fmomei*distanceP->dp2ijpi + fmomej*distanceP->dp2jipi
          	  + fmomdi*distance->dr2ijpj  + fmomdj*distance->dr2jipi;
      neighbor.box->addForceR(j, fmomei*distanceP->dp2ijpj + fmomej*distanceP->dp2jipj
	                       + fmomdi*distance->dr2ijpj  + fmomdj*distance->dr2jipj);

      if(optDerivative) {
          potential->dVdpm();
          forcer[i] +=  potential->forceri;
          neighbor.box->addForceR(j, potential->forcerj);
         // gamma derivative.
          //distanceP->devGamma();
          forcer[i] += distance->devgam1*potential->facmom;
          neighbor.box->addForceR(j, distance->devgam2*potential->facmom);
      }

    }
  }

}

void MBoxSV::computeBUUForceNeighbor(NeighborSV& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = optVdot <= 1 ? part[i]->getPkin() : p1;
    if(optPotentialArg>=1) p1 = part[i]->getPcan(optVdot);
    int  bar1 = part[i]->baryon()/3;

    double fi = p1.mCalc()/p1[0];
    double fengi=fi;
    if(optPotentialArg==3) {
      double meff1 = part[i]->getEffectiveMass();
      fengi = meff1/sqrt(meff1*meff1 + pk1.pAbs2());
    }

    distance->setRP1(r1,p1);
    distanceP->setRP1(r1,p1);

    //potential->set1(p1,pk1,JB[i],fi,fengi,rho[i],rhog[i],part[i]->pots()-vmoms[i],rhos[i],rhosg[i],bar1);
    double* gf1=part[i]->facPotential();
    potential->set1(p1,pk1,JB[i],fi,fengi,rho[i],rhog[i],rhog2[i],
                    part[i]->pots()-vmoms[i],rhos[i],rhosg[i],bar1,gf1[1],gf1[2],gf1[3],gf1[4],gf1[5]);

    for(auto j=0; j< n2; j++) {

      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      Vec4 pk2 = optVdot <= 1 ? part2[j]->getPkin(): p2;
      if(optPotentialArg>=1) p2 = part2[j]->getPcan(optVdot);
      int  bar2 = part2[j]->baryon()/3;

      double fj=p2.mCalc()/p2[0];
      double fengj=fj;
      if(optPotentialArg==3) {
        double meff2 = part[j]->getEffectiveMass();
        fengj = meff2/sqrt(meff2*meff2 + pk2.pAbs2());
      }

      distance->setRP2(r2,p2);
      distanceP->setRP2(r2,p2);

      //potential->set2(p2,pk2,neighbor.box->getJB(j),fj,fengj,neighbor.box->getRho(j),
//	  neighbor.box->getRhog(j),part2[j]->pots()-neighbor.box->getVmoms(j),
//	  neighbor.box->getRhos(j),neighbor.box->getRhosg(j),bar2);

      double* gf2=part2[j]->facPotential();
      potential->set2(p2,pk2,neighbor.box->getJB(j),fj,fengj,neighbor.box->getRho(j),
	  neighbor.box->getRhog(j),
	  neighbor.box->getRhog2(j),
	  part2[j]->pots()-neighbor.box->getVmoms(j),
	  neighbor.box->getRhos(j),neighbor.box->getRhosg(j),
	  bar2,gf2[1],gf2[2],gf2[3],gf2[4],gf2[5]);

      distance->distanceR();
      double rhomij=neighbor.rhom[i][j];
      double rhomji=neighbor.neighbor->rhom[j][i];
      potential->dVdns(rhomij,rhomji);

      double fsky1=potential->fskyi;
      double fsky2=potential->fskyj;

      force[i]  += -fsky1*distance->dr2ijri;
      forcer[i] +=  fsky1*distance->dr2ijpi;
      neighbor.box->addForce(j, -fsky2*distance->dr2jirj);
      neighbor.box->addForceR(j, fsky2*distance->dr2jipj);

      if(optDerivative) {
        //potential->dVdp();
        //forcer[i] +=  potential->forceri;
        //neighbor.box->addForceR(j, potential->forcerj);

        // gamma derivative.
        distance->devGamma();
        forcer[i] += distance->devgam1*potential->facsk;
        neighbor.box->addForceR(j,  distance->devgam2*potential->facsk);
      }

      if(!withMomDep) continue;

      distanceP->distanceP();
      //potential->dVdmd(distanceP->psq1,distanceP->psq2,gf1[6],gf2[6],gf1[7],gf2[7]);
      potential->dVdmd(distanceP->psq1,distanceP->psq2,gf1[6]*gf2[6],gf1[7]*gf2[7],gf1[8],gf2[8],gf1[9],gf2[9]);
      double fmomdi=potential->fmomdi;
      double fmomdj=potential->fmomdj;
      double fmomei=potential->fmomei;
      double fmomej=potential->fmomej;

      force[i]  +=              -fmomdi*distance->dr2ijri;
      neighbor.box->addForce(j, -fmomdj*distance->dr2jirj);

      forcer[i] +=               fmomei*distanceP->dp2ijpi + fmomdi*distance->dr2ijpi;
      neighbor.box->addForceR(j, fmomej*distanceP->dp2jipj + fmomdj*distance->dr2jipj);

      if(optDerivative) {
          //potential->dVdpm();
          //forcer[i] +=  potential->forceri;
          //neighbor.box->addForceR(j, potential->forcerj);

          // gamma derivative.
          //distanceP->devGamma();
          forcer[i] += distance->devgam1*potential->facmom;
          neighbor.box->addForceR(j, distance->devgam2*potential->facmom);
      }


    }
  }

}

/*
// Pre-factor from the baryon current.
Vec4 MBoxSV::facV(int i, Vec4& pkin, NeighborSV& neighbor)
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
*/

//...Compute time-derivatives of the vector potential.
void MBoxSV::computeVdot(NeighborSV& neighbor)
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
    double* gf1=part[i]->facPotential();

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
      double* gf2=part[j]->facPotential();
   
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

      potential->dVdmd(distanceP->psq1,distanceP->psq2,gf1[6]*gf2[6],gf1[7]*gf2[7],gf1[8],gf2[8],gf1[9],gf2[9]);
      fmomdi=potential->fmomdi;
      fmomdj=potential->fmomdj;

      Vdot[i] += doti*fmomdi*vj;
      neighbor.box->addVdot(j, dotj*fmomdj*vi);

      distance->distanceP();

      pdoti=dot3(forcei,distance->dp2ijpi) + dot3(forcej,distance->dp2ijpj);
      pdotj=dot3(forcej,distance->dp2jipj) + dot3(forcei,distance->dp2jipi);

      //Vdot[i] += pdoti*devVme(psq2)*rhomij*vj;
      //neighbor.box->addVdot(j, pdotj*devVme(psq1)*rhomji*vi);

      Vdot[i] += pdoti*potential->fmomdv1*vj;
      neighbor.box->addVdot(j, pdotj*potential->fmomdv2*vi);

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
