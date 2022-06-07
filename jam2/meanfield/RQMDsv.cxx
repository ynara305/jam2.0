#include <jam2/meanfield/RQMDsv.h>
#include <jam2/meanfield/ScalarPotential.h>
#include <jam2/meanfield/VectorPotential.h>
#include <jam2/meanfield/ScalarVectorPotential.h>

// Relativistic quantum molecular dynamics with scalar and vector potentials.
// optPotentialType = 3: non-linear vector potential by
//  M I Gorenstein, D H Rischkef, H Stocker, W Greinert and K A Bugaev
//  J. Phys. G: Nucl. Part. Phys. 19 (1993) L69-L75.
// optPotentialType = 4: non-linear scalar potential
// optPotentialType = 5: non-linear sigma field

namespace jam2 {

using namespace std;

bool RQMDsv::firstCall=true;

RQMDsv::RQMDsv(Pythia8::Settings* s) : MeanField(s)
{
  optBaryonCurrent = settings->mode("MeanField:optBaryonCurrent"); 
  double widG = settings->parm("MeanField:gaussWidth");
  wG = 1.0/(4*widG);

  overSample = settings->mode("Cascade:overSample");
  transportModel = settings->mode("MeanField:transportModel"); 
  if(transportModel==2) {
    wG = 1.0/(2*widG);
    //selfInt = 0;
  }

  int mode = settings->mode("MeanField:potentialType"); 
  if(mode==1) {
    potential = new ScalarPotential(settings);
    optVdot=0;
  } else if(mode==2) {
    potential = new VectorPotential(settings);
  } else if(mode==3) {
    potential = new ScalarVectorPotential(settings);
  } else {
    cout << "RQMDsv wrong mode MeanField:potentialType="<<mode<<endl;
    exit(1);
  }

  optPotentialType = potential->potentialType();
  //=1:Skyrme scalar, =2:Skyrme vector
  //=3:scalar + nonlinear vector potential (M I Gorenstein, D H Rischke, et.al JPG 19(1993)L69)
  //=4:nonlinear scalar potential + vector
  //=5 nonlinear sigma-omega
  t0 = potential->getT0();
  t2 = potential->getT2();
  t1 = potential->getT1();
  t3 = potential->getT3();
  gam = potential->getGam();
  withMomDep = potential->isMomDep();
  if(firstCall) {
  cout << "RQMDsv potentialType= " << optPotentialType<<endl;
  cout << "t0= "<< t0
       << " t2= "<< t2
       << " t1= "<< t1
       << " t3= "<< t3
       << " gam= "<< gam
       <<endl;
  firstCall=false;
  }


}

RQMDsv::~RQMDsv()
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
  for(int i=0;i<NV;i++) rhom[i].clear();
  rhom.clear();
  JB.clear();
  Vdot.clear();
  delete potential;
}

void RQMDsv::evolution(list<EventParticle*>& plist, double t, double dt,int step)
{
  part.clear();
  globalTime = t;

  // import particles.
  for(auto& i : plist) {
    if(i->isMeanField(t,optPotential)) {
      part.push_back(i);
      potential->setPotentialParam(i);
    }
  }

  NV = part.size();
  rhom.resize(NV);
  for(int i=0;i<NV;i++) rhom[i].resize(NV);

  rhog.assign(NV,0.0);
  rhog2.assign(NV,0.0);
  rhosg.assign(NV,0.0);
  rhos.assign(NV,0.0);
  rho.assign(NV,0.0);
  vmoms.assign(NV,0.0);
  vmom4.assign(NV,0.0);
  JB.assign(NV,0.0);
  Vdot.assign(NV,0.0);
  force.assign(NV,0.0);
  forcer.assign(NV,0.0);

  // change from kinetic to canonical momenta.
  if(optVdot==1) {
    for(auto& i : part) i->setFree();
  }
  if(optPotentialType !=2 && optPotentialArg >= 1) {
    for(auto& i : part) i->setPotS(0.0);
  }

  qmdMatrix();
  bool optSet= (optPotentialType !=2 && optPotentialArg == 0) ? true :  false;
  singleParticlePotential(optSet);
  if(transportModel==1) computeForce();
  else computeBUUForce();

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

    computeVdot();
    for(int i=0; i< NV; i++) {
      part[i]->updateByForce(force[i] - Vdot[i], forcer[i],dt);
    }

  // compute time derivatives of V^\mu numerically; (V(t+dt) - V(t))/dt
  // (under construction)
  } else if(optVdot == 4 ) {

    for(int i=0;i<NV;i++) {
      part[i]->updateByForce(force[i] + part[i]->potv()/dt, forcer[i],dt);
      vmom4[i]=0.0;
      JB[i]=0.0;
    }
    vmom4.assign(NV,0.0);
    JB.assign(NV,0.0);
    qmdMatrix();
    singleParticlePotential(optSet);
    for(int i=0;i<NV;i++) {
      part[i]->addP(-part[i]->potv());
      part[i]->setOnShell();
    }

  } else {
    cout << "RQMDsv wrong option optVdot = "<< optVdot<<endl;
    exit(1);
  }

  if(!optSet) setScalarPotential();

  //qmdMatrix();
  //singleParticlePotential(true);

}

void RQMDsv::setScalarPotential()
{
  for(int i=0;i<NV;i++) {
    part[i]->setPotS(part[i]->pots());
  }
}

// Compute single particle potential.
void RQMDsv::singleParticlePotential(bool optSet)
{
  // set scalar potential.
  if(optPotentialType !=2) {  // exclude RQMDv mode
    if(optPotentialType==5) { // compute sigma-field
      sigmaField();
    } else {
      for(int i=0; i< NV; i++) {
        double* pfac=part[i]->facPotential();
        if(t2 !=0 && rhos[i]>1e-8) rhosg[i] = pow(rhos[i],pfac[3]-1);
        part[i]->setPotS((pfac[1]*t0 + t2*rhosg[i])*rhos[i] + vmoms[i],optSet);
        part[i]->setPotSm(vmoms[i]);
        part[i]->setRhoS(rhos[i]);
      }
    }
  }

  // RQMDs mode:scalar potential only.
  if(optPotentialType==1) return;

  // Compute invariant baryon density.
  for(int i=0; i< NV; i++) {
    double* pfac=part[i]->facPotential();
    rho[i] = sqrt(max(0.0,JB[i].m2Calc()));
    if(rho[i] > 1e-8) {
      rhog[i]  = pow(rho[i],pfac[3]-1);
      rhog2[i] = pow(rho[i],pfac[5]-1);
    }
    part[i]->setRhoB(rho[i]);
  }

  Vec4 vTotal=0.0;
  Vec4 vc = 0.0;
  bool optV=false;
  if(optVdot==1) optV=true;
  if(optV && optVectorPotential==1) {
    for(int i=0; i< NV; i++) {
      double* pfac=part[i]->facPotential();
      double vsky = part[i]->baryon()/3*(pfac[1]*t1 + pfac[2]*t3*rhog[i] + pfac[4]*rhog2[i]);
      vTotal += vsky * JB[i] + vmom4[i];
    }
    vc = vTotal / (NV*overSample);
    vc[0]=0.0;
  }

  Vec4 vdott = 0.0;
  for(int i=0; i< NV; i++) {

    if(rho[i] < 1e-15) continue;
    Vec4 vpot0 = part[i]->potv();
    double* pfac=part[i]->facPotential();
    double vsky = part[i]->baryon()/3*(pfac[1]*t1 + pfac[2]*t3*rhog[i] + pfac[4]*rhog2[i]);

    // four-components of vector potential are fully included.
    if(optVectorPotential==1) {
      part[i]->setPotV( vsky * JB[i] + vmom4[i] -vc);
      part[i]->setPotVm(vmom4[i] -vc);

    // only time component of vector potential is included.
    } else if(optVectorPotential==2) {
      part[i]->setPotV(0,vsky*JB[i][0] + vmom4[i][0] );
      part[i]->setPotVm(0,vmom4[i][0] );

    // only time component of vector potential is included in the form
    // of V(rho_B) where rho_B is an invariant baryon density.
    } else if(optVectorPotential==3) {
      part[i]->setPotV(0,vsky*rho[i] + vmom4[i][0] );
      part[i]->setPotVm(0,vmom4[i][0] );
    } else {
      part[i]->setPotV(0,
	  part[i]->baryon()/3*(t1*rho[i]+t3*pow(max(0.0,rho[i]),gam))
	  + vmom4[i][0] );
      part[i]->setPotVm(0,vmom4[i][0] );
    }

    if(optVdot == 2) {
      Vdot[i] = part[i]->potv() - vpot0;
      vdott += Vdot[i];
    }
  }

  if(optVdot == 2) {
    vdott /= NV;
    for(int i=0; i< NV; i++) Vdot[i] -= vdott;
  } 

}

// compute single particle potential energy.
Vec4 RQMDsv::computeEnergy(list<EventParticle*>& plist, int step)
{
  pTot=0.0;   
  for(auto& i :plist) {
    double m = i->getEffectiveMass();
    if(optVectorPotential==1 && optVdot==0) {
      Vec4 pk= i->getP() - i->potv();
      pTot[0] += sqrt( m*m + pk.pAbs2());
    } else {
      pTot[0] += sqrt( m*m + i->pAbs2());
    }
    pTot[1] += i->getP(1);
    pTot[2] += i->getP(2);
    pTot[3] += i->getP(3);
    pTot[0] += i->potv(0);
  }

  if(step==1) pTot0 = pTot/overSample;

  return pTot/overSample;
}


void RQMDsv::qmdMatrix()
{
  //double gf1[3]={1.0, 1.0, 1.0};
  //double gf2[3]={1.0, 1.0, 1.0};

  for(int i=0; i< NV; i++) {
    Vec4 r1  = part[i]->getR();
    Vec4 p1 = optPotentialArg>=1? part[i]->getPcan(optVdot) : part[i]->getP();
    distance->setRP1(r1,p1);
    distanceP->setRP1(r1,p1);
    Vec4 v1 = p1/p1[0];
    double fi=p1.mCalc()/p1[0];
    int potid1=part[i]->potentialID();

    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;
    Vec4 fri = optBaryonCurrent ? part[i]->forceR() : 0.0;
    int bi   = part[i]->baryon()/3;
    if(bi==0) {
      qfac1 = facMesonPot;
      bi=1;
    }
    double* gfac1= part[i]->facPotential();
    
    /*
    if(optPotentialDensity==1) {
      gf1[0]=gfac1[1];
      gf1[1]=gfac1[2];
      gf1[2]=gfac1[4];
    }
    */

    for(int j=i+selfInt; j< NV; j++) {

      int potid2=part[j]->potentialID();
      double pot1 = potential->matrix(potid1,potid2);
      double pot2 = potential->matrix(potid2,potid1);

      Vec4 r2  = part[j]->getR();
      Vec4 p2 = optPotentialArg>=1 ? part[j]->getPcan(optVdot) : part[j]->getP();
      distance->setRP2(r2,p2);
      distanceP->setRP2(r2,p2);
      Vec4 v2 = p2/p2[0];
      double fj=p2.mCalc()/p2[0];

      Vec4 frj = optBaryonCurrent ? part[j]->forceR() : 0.0;
      int bj = part[j]->baryon()/3;
      double qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;
      if(bj==0) {
	qfac2 = facMesonPot;
	bj=1;
      }
      double* gfac2= part[j]->facPotential();
      /*
      if(optPotentialDensity==1) {
        gf2[0]=gfac2[1];
        gf2[1]=gfac2[2];
        gf2[2]=gfac2[4];
      }
      */
      
      distance->density();
      rhom[i][j] = distance->density1*qfac1*qfac2/overSample*pot1;
      rhom[j][i] = distance->density2*qfac1*qfac2/overSample*pot2;

      // scalar density
      rhos[i] += rhom[i][j]*fj;
      rhos[j] += rhom[j][i]*fi;

      // vector current
      JB[i] += rhom[i][j]*(v2+frj)*bj;
      JB[j] += rhom[j][i]*(v1+fri)*bi;

      if(!withMomDep) continue;
      distanceP->psq();
      potential->Vmd(distanceP->psq1,distanceP->psq2,gfac1[6]*gfac2[6],gfac1[7]*gfac2[7]
	  ,gfac1[8],gfac1[9],gfac2[8],gfac2[9]);

      // momentum-dependent scalar potential
      vmoms[i] += potential->vmomsi*rhom[j][i]*fj;
      vmoms[j] += potential->vmomsj*rhom[j][i]*fi;

      // momentum-dependent vector potential
      vmom4[i] += potential->vmom4i*rhom[i][j]*v2;
      vmom4[j] += potential->vmom4j*rhom[j][i]*v1;

    }
  }

}

void RQMDsv::computeForce()
{
  for(int i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = optVdot <= 1 ? part[i]->getPkin() : p1;
    if(optPotentialArg>=1) p1 = part[i]->getPcan(optVdot);
    int bar1= part[i]->baryon()/3;

    double fi = p1.mCalc()/p1[0];
    double fengi=fi;
    if(optPotentialArg==3) {
      double meff1 = part[i]->getEffectiveMass();
      fengi = meff1/sqrt(meff1*meff1 + pk1.pAbs2());
    }

    double* gf1=part[i]->facPotential();
    distance->setRP1(r1,p1);
    distanceP->setRP1(r1,p1);
    potential->set1(p1,pk1,JB[i],fi,fengi,rho[i],rhog[i],rhog2[i],
                    part[i]->pots()-vmoms[i],rhos[i],rhosg[i],bar1,gf1[1],gf1[2],gf1[3],gf1[4],gf1[5]);

    for(int j=i+1; j< NV; j++) {

      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      Vec4 pk2 = optVdot <= 1 ? part[j]->getPkin() : p2;
      if(optPotentialArg>=1) p2 = part[j]->getPcan(optVdot);
      int bar2 = part[j]->baryon()/3;

      double fj=p2.mCalc()/p2[0];
      double fengj=fj;
      if(optPotentialArg==3) {
        double meff2 = part[j]->getEffectiveMass();
        fengj = meff2/sqrt(meff2*meff2 + pk2.pAbs2());
      }

      distance->setRP2(r2,p2);
      distanceP->setRP2(r2,p2);
      double* gf2=part[j]->facPotential();
      potential->set2(p2,pk2,JB[j],fj,fengj,rho[j],rhog[j],rhog2[j],
	  part[j]->pots()-vmoms[j],rhos[j],rhosg[j],bar2,gf2[1],gf2[2],gf2[3],gf2[4],gf2[5]);

      distance->distanceR();
      potential->dVdns(rhom[i][j],rhom[j][i]);

      double fsky1=potential->fskyi;
      double fsky2=potential->fskyj;

      force[i]  += -fsky1*distance->dr2ijri - fsky2*distance->dr2jiri;
      force[j]  += -fsky1*distance->dr2ijrj - fsky2*distance->dr2jirj;
      forcer[i] +=  fsky1*distance->dr2ijpi + fsky2*distance->dr2jipi;
      forcer[j] +=  fsky1*distance->dr2ijpj + fsky2*distance->dr2jipj;

      // Derivative of p0/m and p/p0 term.
      if(optDerivative) {

        potential->dVdp();
        forcer[i] +=  potential->forceri;
        forcer[j] +=  potential->forcerj;

        // gamma derivative.
        distance->devGamma();
        forcer[i] += distance->devgam1*potential->facsk;
        forcer[j] += distance->devgam2*potential->facsk;
      }

      if(!withMomDep) continue;

      distanceP->distanceP();
      potential->dVdmd(distanceP->psq1,distanceP->psq2,
	  gf1[6]*gf2[6],gf1[7]*gf2[7],gf1[8],gf2[8],gf1[9],gf2[9]);
      double fmomdi=potential->fmomdi;
      double fmomdj=potential->fmomdj;
      double fmomei=potential->fmomei;
      double fmomej=potential->fmomej;

      force[i]  += -fmomdi*distance->dr2ijri   - fmomdj*distance->dr2jiri;
      force[j]  += -fmomdi*distance->dr2ijrj   - fmomdj*distance->dr2jirj;

      forcer[i] +=  fmomei*distanceP->dp2ijpi + fmomej*distanceP->dp2jipi
	          + fmomdi*distance->dr2ijpi  + fmomdj*distance->dr2jipi;
      forcer[j] +=  fmomei*distanceP->dp2ijpj + fmomej*distanceP->dp2jipj
	          + fmomdi*distance->dr2ijpj  + fmomdj*distance->dr2jipj;

      if(optDerivative) {
        potential->dVdpm();
        forcer[i] +=  potential->forceri;
        forcer[j] +=  potential->forcerj;
        // gamma derivative.
        //distanceP->devGamma();
        forcer[i] += distance->devgam1*potential->facmom;
        forcer[j] += distance->devgam2*potential->facmom;
      }

    } // end loop over j
  } // end loop over i

}

void RQMDsv::computeBUUForce()
{
  for(int i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = optVdot <= 1 ? part[i]->getPkin() : p1;
    if(optPotentialArg>=1) p1 = part[i]->getPcan(optVdot);
    int bar1= part[i]->baryon()/3;

    double fi = p1.mCalc()/p1[0];
    double fengi=fi;
    if(optPotentialArg==3) {
      double meff1 = part[i]->getEffectiveMass();
      fengi = meff1/sqrt(meff1*meff1 + pk1.pAbs2());
    }

    double* gf1=part[i]->facPotential();
    distance->setRP1(r1,p1);
    distanceP->setRP1(r1,p1);
    potential->set1(p1,pk1,JB[i],fi,fengi,rho[i],rhog[i],rhog2[i],
                    part[i]->pots()-vmoms[i],rhos[i],rhosg[i],bar1,gf1[1],gf1[2],gf1[3],gf1[4],gf1[5]);

    for(int j=i+selfInt; j< NV; j++) {

      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      Vec4 pk2 = optVdot <= 1 ? part[j]->getPkin() : p2;
      if(optPotentialArg>=1) p2 = part[j]->getPcan(optVdot);
      int bar2 = part[j]->baryon()/3;

      double fj=p2.mCalc()/p2[0];
      double fengj=fj;
      if(optPotentialArg==3) {
        double meff2 = part[j]->getEffectiveMass();
        fengj = meff2/sqrt(meff2*meff2 + pk2.pAbs2());
      }

      double* gf2=part[j]->facPotential();
      distance->setRP2(r2,p2);
      distanceP->setRP2(r2,p2);
      potential->set2(p2,pk2,JB[j],fj,fengj,rho[j],rhog[j],rhog2[j],
	  part[j]->pots()-vmoms[j],rhos[j],rhosg[j],bar2,gf2[1],gf2[2],gf2[3],gf2[4],gf2[5]);

      distance->distanceR();
      potential->dVdns(rhom[i][j],rhom[j][i]);

      double fsky1=potential->fskyi;
      double fsky2=potential->fskyj;

      force[i]  += -fsky1*distance->dr2ijri;
      force[j]  += -fsky2*distance->dr2jirj;
      forcer[i] +=  fsky1*distance->dr2ijpi;
      forcer[j] +=  fsky2*distance->dr2jipj;

      // Derivative of p0/m and p/p0 term.
      if(optDerivative) {
        //potential->dVdp();
        //forcer[i] +=  potential->forceri;
        //forcer[j] +=  potential->forcerj;
        // gamma derivative.
        distance->devGamma();
        forcer[i] += distance->devgam1*potential->facsk;
        forcer[j] += distance->devgam2*potential->facsk;
      }

      if(!withMomDep) continue;

      distanceP->distanceP();
      //potential->dVdmd(distanceP->psq1,distanceP->psq2,gf1[6],gf2[6],gf1[7],gf2[7]);
      potential->dVdmd(distanceP->psq1,distanceP->psq2,gf1[6]*gf2[6],gf1[7]*gf2[7],gf1[8],gf2[8],gf1[9],gf2[9]);
      double fmomdi=potential->fmomdi;
      double fmomdj=potential->fmomdj;
      double fmomei=potential->fmomei;
      double fmomej=potential->fmomej;

      force[i]  += -fmomdi*distance->dr2ijri;
      force[j]  += -fmomdj*distance->dr2jirj;

      forcer[i] +=  fmomei*distanceP->dp2ijpi + fmomdi*distance->dr2ijpi;
      forcer[j] +=  fmomej*distanceP->dp2jipj + fmomdj*distance->dr2jipj;

      if(optDerivative) {
        //potential->dVdpm();
        //forcer[i] +=  potential->forceri;
        //forcer[j] +=  potential->forcerj;
        // gamma derivative.
        //distanceP->devGamma();
        forcer[i] += distance->devgam1*potential->facmom;
        forcer[j] += distance->devgam2*potential->facmom;
      }

    } // end loop over j
  } // end loop over i

}

//*********************************************************************
//...Compute time-derivatives of the vector potential.
void RQMDsv::computeVdot()
{
  bool opt=true;
  //opt=false;

  for(int i=0;i<NV; i++) Vdot[i]=0.0;

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

    distance->setRP1(r1,p1);
    distanceP->setRP1(r1,p1);
    double* gf1=part[i]->facPotential();

    for(int j=i+1;j<NV; j++) {

      double vvj = t1 + t3*rhog[j];
      double dvj=0.0;
      if(abs(rho[j]) > 1e-7) dvj = (gam-1.0)*t3*pow(rho[j],gam-3.0);
      Vec4 Bj = dvj*JB[j];
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      Vec4 vk2 = p2/p2[0];
      if(optPotentialArg==1) p2 = part[j]->getPcan(optVdot);
      double emj=p2[0];
      if(optTwoBodyDistance==3)  emj = p2.mCalc();
      Vec4 vj = p2 / emj;
      Vec4 b2 = p2 / p2[0];

      distance->setRP2(r2,p2);
      distanceP->setRP2(r2,p2);
      double* gf2=part[j]->facPotential();
   
      Vec4 Aj = (vj * JB[i]) * Bi;
      Vec4 Ai = (vi * JB[j]) * Bj;

      // Compute derivatives
      //distance->distanceR();

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
      //Vec4  dr2ri = distance->dr2ri; 
      //Vec4  dr2rj = distance->dr2rj;
      //Vec4  dr2pi = distance->dr2pj;
      //Vec4  dr2pj = distance->dr2pj;

      //Vec4  drji2rj = distance->drji2rj;
      //Vec4  drji2rj = distance->drji2rj;
      //Vec4  drji2pi = distance->drji2pi;
      //Vec4  drji2pj = distance->drji2pj;

      double xdoti=dot3(vk1+forcer[i],dr2ri)   + dot3(vk2+forcer[j],dr2rj);
      double xdotj=dot3(vk2+forcer[j],drji2rj) + dot3(vk1+forcer[i],drji2ri);

      double pdoti=dot3(force[i],dr2pi-dgami/wG)
	         + dot3(force[j],dr2pj-dgamj/wG);
      double pdotj=dot3(force[j],drji2pj-dgamj/wG)
	         + dot3(force[i],drji2pi-dgami/wG);

      double doti=-wG*(xdoti+pdoti)*rhom[i][j];
      double dotj=-wG*(xdotj+pdotj)*rhom[j][i];
      Vdot[i] += doti*(Aj + vvi*vj);
      Vdot[j] += dotj*(Ai + vvj*vi);

      // Momentum dependent potential.
      double fmomdi=0.0, fmomdj=0.0;
      if(withMomDep) {

      potential->dVdmd(distanceP->psq1,distanceP->psq2,gf1[6]*gf2[6],gf1[7]*gf2[7],gf1[8],gf2[8],gf1[9],gf2[9]);
      fmomdi=potential->fmomdi;
      fmomdj=potential->fmomdj;

      Vdot[i] += doti*fmomdi*vj;
      Vdot[j] += dotj*fmomdj*vi;


      //distance->distanceP();
      distanceP->distanceP();
      pdoti=dot3(force[i],distanceP->dp2ijpi) + dot3(force[j],distanceP->dp2ijpj);
      pdotj=dot3(force[j],distanceP->dp2jipj) + dot3(force[i],distanceP->dp2jipi);

      //Vdot[i] += pdoti*devVme(psq2)*rhom[i][j]*vj;
      //Vdot[j] += pdotj*devVme(psq1)*rhom[j][i]*vi;

      Vdot[i] += pdoti*potential->fmomdv1*vj;
      Vdot[j] += pdotj*potential->fmomdv2*vi;

      }


      if(opt) continue;

      //Compute derivatives of pre-factors v_j^\mu in the vector potential.
      double fai = - dot3(force[j],Bi/emj) * rhom[i][j];
      double faj = - dot3(force[i],Bj/emi) * rhom[j][i];
      double fai2 = dot3(force[j],vj);
      double faj2 = dot3(force[i],vi); 
      double fai3 = dot3(vj,Bi/p2[0]);
      double faj3 = dot3(vi,Bj/p1[0]);

      double  fvi2=0.0;
      double  fvj2=0.0;
      if(optTwoBodyDistance==3) {
        fai3=Bi[0]/p2[0];
        faj3=Bj[0]/p1[0];
      } else {
        fvi2=-optP0dev*fai2*rhom[i][j]*vvi/p2[0];
        fvj2=-optP0dev*faj2*rhom[j][i]*vvj/p1[0];
      }

      fai += optP0dev*fai2*fai3*rhom[i][j];
      faj += optP0dev*faj2*faj3*rhom[j][i];

      double fvi = vvi*rhom[i][j]/emj;
      double fvj = vvj*rhom[j][i]/emi;

      Vdot[i] += fai*JB[i] + fvi*force[j] + fvi2*vj;
      Vdot[j] += faj*JB[j] + fvj*force[i] + fvj2*vi;

      // Momentum dependent potential.
      if(withMomDep) {
        Vdot[i] += fmomdi*rhom[i][j]*force[j]/emj;
        Vdot[j] += fmomdj*rhom[j][i]*force[i]/emi;
        if(optTwoBodyDistance != 3) {
          Vdot[i] += -optP0dev*fai2*rhom[i][j]*fmomdi/p2[0]*vj;
          Vdot[i] += -optP0dev*faj2*rhom[j][i]*fmomdj/p1[0]*vi;
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

  /*
  cout << " vdott = "<< scientific
       << setw(13) << setprecision(6) << vdott[1]
       << setw(13) << setprecision(6) << vdott[2]
       << setw(13) << setprecision(6) << vdott[3]
       <<endl;
  cin.get();
  */

}

void RQMDsv::sigmaField(double a)
{
  int maxit=5;
  const double  eps = 1e-5;
  bool optSet=true;
  if(optPotentialArg >= 1) optSet=false;
  if(optPotentialArg > 0 || optVdot <= 1) maxit=1; 

  double diff=0.0;
  double gs = -2*t0;
  for(int itry=0;itry<maxit;itry++) {
    diff=0.0;
    for(int i=0; i< NV; i++) {
      double vpots0 = part[i]->pots();
      double gfac=1.0;
      scalarDens=gs*rhos[i]/gfac;// factor gr is already contained in rhos(i)
      double gss=t0*gfac*HBARC;

      double sigma1=abs((part[i]->getMass() + vmoms[i])/gss);
      //double sigma1=abs((part[i]->getMass())/gss);
      double vpots=gss*sigma_bisec(sigma1, i ) + vmoms[i];

      //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      //double mass = part[i]->getMass();
      //if(mass+vpots<=0.0) vpots=-mass*0.999;
      //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      //double vpots=gss*sigma_iterate() + vmoms[i];
      //double vpots=gss*sigma_newton() + vmoms[i];

      part[i]->setPotS(vpots,optSet);
      part[i]->setPotSm(vmoms[i]);
      part[i]->setRhoS(rhos[i]);
      diff += abs(vpots-vpots0);
    }

    //diff=abs(pot0-pot);
    if(diff/NV < eps) return;
      //cout << " diff = "<< diff << " itry= "<< itry << " mxit= "<< maxit << endl;

    if(maxit==1) return;
    qmdMatrix();

  }

    cout << "RQMDsv::sigma does not converge diff= " << diff 
	<< " itry= "<< maxit << endl;

}

// Find sigma field in 1/fm
double RQMDsv::sigma_bisec(double sigma1, int i)
{
  ScalarVectorPotential *potential3=dynamic_cast<ScalarVectorPotential*>(potential);
  double sigma0=0.0;
  //double f0 = funcSigma(sigma0);
  //double f1 = funcSigma(sigma1);
  double f0 = potential3->funcSigma(sigma0,scalarDens);
  double f1 = potential3->funcSigma(sigma1,scalarDens);

  double sig0=sigma1;
  // seeking a new initial condition.
  if(f0*f1 > 0.0) {
    do {
      sigma1 -= 0.03;
      if(sigma1 < sigma0) {
	  cout << "sigma_bisec no solution  sigma1 = "<< sigma1
	      << " rhos= "<< scalarDens
	      << " sigma0= "<< sig0
	      << " id= "<< part[i]->getID()
	      << " mass= "<< part[i]->getEffectiveMass()
	      <<endl;
	  return sigma1-1e-5;
      }
      //f1=funcSigma(sigma1);
      f1=potential3->funcSigma(sigma1,scalarDens);
      //cout << " sigma1 = " << sigma1  << " f1= "<< f1 <<endl;
      if(f0*f1 < 0.0) break;
    }while(sigma1 > sigma0);
  }

  // check if sigma field is positive.
  if(sigma1 <0.0) {
      cout << " sigma < 0 ? "<< sigma1 << " sig0= "<< sig0
	   <<endl;
      exit(1);
  }

  // check if initial condition is ok.
  if(f0*f1 > 0.0) {
    cout << "bisec no solution sigma= " << sigma1
	 << " dens= "<< scalarDens <<endl;
    cout << "func1= " << potential3->funcSigma(sigma0,scalarDens) << " f2= "<< potential3->funcSigma(sigma1,scalarDens)
	<<endl;
    return 0.0;
  }

  // bisection method start.
  double f=1.0;
  double sigma=0.0;
  int itry=0;
  do {
    sigma=0.5*(sigma0+sigma1);
    f =  potential3->funcSigma(sigma,scalarDens);
    if(f*f1 > 0.0) sigma1=sigma;
    else sigma0=sigma;
    if(++itry > 50) {
      cout << "does not converge sigma_bisec " << sigma << endl;
      exit(1);
    }
  }while (abs(f) > 1e-5);

  return sigma;
}


// Pre-factor from the baryon current.
Vec4 RQMDsv::facV(int i, Vec4& pkin)
{
  if(rho[i]<1e-15) return Vec4(0.0);
  Vec4 bj=0.0;

  // 4-components of the vector potential is fully included.
  if(optVectorPotential==1) {
    bj = pkin/pkin[0];

  // Only time-component V^0 term is included.
  } else if(optVectorPotential==2) {
    bj[0]=1.0;

  // Only time-component V^0 term is included with the form of V(rho_B).
  } else if (optVectorPotential==3) {
    return (t1 + gam*t3*rhog[i])/rho[i] * JB[i];
  }

  double vj = JB[i] * bj;
  double vv = t1 + t3*rhog[i];  // V/rho_B
  double dv = (gam-1.0)*t3*pow(rho[i],gam-3);  // del(V/rho)/rho
  return vj*dv*JB[i] + vv*bj;

}

/*
void RQMDsv::qmdMatrix0()
{
  for(int i=0; i< NV; i++) {
    Vec4 r1  = part[i]->getR();
    Vec4 p1 = optPotentialArg>=1? part[i]->getPcan(optVdot) : part[i]->getP();
    distance->setRP1(r1,p1);
    Vec4 v1 = p1/p1[0];
    double fi=p1.mCalc()/p1[0];

    Vec4 fri = optBaryonCurrent ? part[i]->forceR() : 0.0;
    int bi   = part[i]->baryon()/3;
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;

    for(int j=i+1; j< NV; j++) {
      Vec4 r2  = part[j]->getR();
      Vec4 p2 = optPotentialArg>=1 ? part[j]->getPcan(optVdot) : part[j]->getP();
      distance->setRP2(r2,p2);
      Vec4 v2 = p2/p2[0];
      double fj=p2.mCalc()/p2[0];

      Vec4 frj = optBaryonCurrent ? part[j]->forceR() : 0.0;
      int bj = part[j]->baryon()/3;
      double qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;
      
      distance->density();
      rhom[i][j] = distance->density1*qfac1*qfac2;
      rhom[j][i] = distance->density2*qfac1*qfac2;

      // scalar density
      rhos[i] += rhom[i][j]*fj;
      rhos[j] += rhom[j][i]*fi;

      // vector current
      JB[i] += rhom[i][j]*(v2+frj)*bj;
      JB[j] += rhom[j][i]*(v1+fri)*bi;

      if(!withMomDep) continue;
      distance->psq();

      // momentum-dependent scalar potential
      vmoms[i] += vex1/(1.0-distance->psq1/pmu1)*rhom[i][j]*fj;
      vmoms[j] += vex1/(1.0-distance->psq2/pmu1)*rhom[j][i]*fi;

      // momentum-dependent vector potential
      vmom4[i] += vex2/(1.0-distance->psq1/pmu2)*rhom[i][j]*v2;
      vmom4[j] += vex2/(1.0-distance->psq2/pmu2)*rhom[j][i]*v1;


    }
  }

}

void RQMDsv::computeForce0()
{
      double C1 = -0.17086291974074935;
      double mu1 = 3.1466990715061636*HBARC;
      vex1=C1/(2*0.168);
      vex2=0.0;
      pmu1=mu1*mu1;
      pmu2=1.0;

  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = optVdot <= 1 ? part[i]->getPkin() : p1;
    if(optPotentialArg>=1) p1 = part[i]->getPcan(optVdot);
    int bar1= part[i]->baryon()/3;

    distance->setRP1(r1,p1);
    Vec4 v1 = p1/p1[0];
    double fi = p1.mCalc()/p1[0];
    double fengi=fi;
    if(optPotentialArg==3) {
      double meff1 = part[i]->getEffectiveMass();
      fengi = meff1/sqrt(meff1*meff1 + p1.pAbs2());
    }

    //double s1 = optPotentialType==3 ? dSigma(i) : t0+t2*gam*rhog[i];
    double s1 = t0+t2*gam*rhog[i];

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      Vec4 pk2 = optVdot <= 1 ? part[j]->getPkin() : p2;
      if(optPotentialArg>=1) p2 = part[j]->getPcan(optVdot);
      int bar2 = part[j]->baryon()/3;

      distance->setRP2(r2,p2);
      Vec4 v2 = p2/p2[0];
      double fj=p2.mCalc()/p2[0];
      double fengj=fj;
      if(optPotentialArg==3) {
        double meff2 = part[j]->getEffectiveMass();
        fengj = meff2/sqrt(meff2*meff2 + p2.pAbs2());
      }

      distance->distanceR();
      //double s2 = optPotentialType==3 ? dSigma(j) : t0+t2*gam*rhog[j];
      double s2 = t0+t2*gam*rhog[j];

      // scalar part
      double fsky1 = -wG*fengi*s1*rhom[i][j]*fj; 
      double fsky2 = -wG*fengj*s2*rhom[j][i]*fi; 
      force[i]  += -fsky1*distance->dr2ri - fsky2*distance->dr2rj;
      force[j]  +=  fsky1*distance->dr2ri + fsky2*distance->dr2rj;
      forcer[i] +=  fsky2*distance->dr2pi;
      forcer[j] +=  fsky1*distance->dr2pj;

      // vector part
      Vec4 Ai = facV(i,p2);
      Vec4 Aj = facV(j,p1);
      double fskyi = -wG*(Ai * pk1)/pk1[0]*rhom[i][j]*bar1*bar2; 
      double fskyj = -wG*(Aj * pk2)/pk2[0]*rhom[j][i]*bar1*bar2; 
      force[i]  += -fskyi*distance->dr2ri - fskyj*distance->dr2rj;
      force[j]  +=  fskyi*distance->dr2ri + fskyj*distance->dr2rj;
      forcer[i] +=  fskyj*distance->dr2pi;
      forcer[j] +=  fskyi*distance->dr2pj;

      // Derivative of p^\mu/m_j term in the vector potential.
      if(optDerivative) {

        if(optTwoBodyDistance!=3) {
        forcer[i] +=  fsky2*p1/(wG*p1[0]*p1[0]);
        forcer[j] +=  fsky1*p2/(wG*p2[0]*p2[0]);
        // gamma derivative.
        distance->devGamma();
        double facsk= -(fsky1+fsky2)/wG;
        forcer[i] += distance->devgam1*facsk;
        forcer[j] += distance->devgam2*facsk;
	}

        Vec4 A1 = facV(i,pk1);
        Vec4 A2 = facV(j,pk2);
	distance->devV(A1,A2);
        forcer[i] += distance->devV1*rhom[j][i]*bar1*bar2;
        forcer[j] += distance->devV2*rhom[i][j]*bar1*bar2;

	double facsk= -(fskyi+fskyj)/wG;
        forcer[i] += distance->devgam1*facsk;
        forcer[j] += distance->devgam2*facsk;
      }

      if(!withMomDep) continue;

      distance->distanceP();
      double psq1=distance->psq1;
      double psq2=distance->psq2;

      // scalar part
      double fmomd1 = -wG*fengi*devVmd(psq2,pmu1,vex1)*rhom[i][j]*fj;
      double fmomd2 = -wG*fengj*devVmd(psq1,pmu1,vex1)*rhom[j][i]*fi;
      double fmome1 =     fengi*devVme(psq2,pmu1,vex1)*rhom[i][j]*fj;
      double fmome2 =     fengj*devVme(psq1,pmu1,vex1)*rhom[j][i]*fi;

      force[i]  += -fmomd1*distance->dr2ri   - fmomd2*distance->dr2rj;
      force[j]  +=  fmomd1*distance->dr2ri   + fmomd2*distance->dr2rj;
      forcer[i] +=  fmome1*distance->dp2ijpi + fmome2*distance->dp2jipi + fmomd2*distance->dr2pi;
      forcer[j] +=  fmome2*distance->dp2jipj + fmome1*distance->dp2ijpj + fmomd1*distance->dr2pj;


      // vector part
      double vf1=1.0;
      double vf2=1.0;
      if(optVectorPotential==1) {
        vf1 = pk1/pk1[0] * v2;
        vf2 = pk2/pk2[0] * v1;
      }

      double fmomdi = -wG*devVmd(psq2,pmu2,vex2)*rhom[i][j]*vf1;
      double fmomdj = -wG*devVmd(psq1,pmu2,vex2)*rhom[j][i]*vf2;
      double fmomei =     devVme(psq2,pmu2,vex2)*rhom[i][j]*vf1;
      double fmomej =     devVme(psq1,pmu2,vex2)*rhom[j][i]*vf2;

      force[i]  += -fmomdi*distance->dr2ri   - fmomdj*distance->dr2rj;
      force[j]  +=  fmomdi*distance->dr2ri   + fmomdj*distance->dr2rj;
      forcer[i] +=  fmomei*distance->dp2ijpi + fmomej*distance->dp2jipi + fmomdj*distance->dr2pi;
      forcer[j] +=  fmomej*distance->dp2jipj + fmomei*distance->dp2ijpj + fmomdi*distance->dr2pj;

      if(optDerivative) {

        if(optTwoBodyDistance!=3 && optP0dev) {
          forcer[i] +=  fmomd2*p1/(wG*p1[0]*p1[0]);
          forcer[j] +=  fmomd1*p2/(wG*p2[0]*p2[0]);
	}
        // gamma derivative.
        double facmom = -(fmomd1 + fmomd2)/wG;
        forcer[i] += distance->devgam1*facmom;
        forcer[j] += distance->devgam2*facmom;


        // Derivative of p_j / p_0j part.
        double facm1 = devVmd(psq1,pmu2,vex2)*rhom[j][i];
        double facm2 = devVmd(psq2,pmu2,vex2)*rhom[i][j];
	distance->devV(pk1/pk1[0],pk2/pk2[0]);
        forcer[i] += facm1*distance->devV1;
        forcer[j] += facm2*distance->devV2;

	// gamma derivative.
        facmom = facm1*vf1 + facm2*vf2;
        forcer[i] += distance->devgam1*facmom;
        forcer[j] += distance->devgam2*facmom;
      }

    } // end loop over j
  } // end loop over i

}
*/

} // namespace jam2


