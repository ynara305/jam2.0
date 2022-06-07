#include <jam2/meanfield/RQMDv.h>

// Relativistic quantum molecular dynamics with Skyrme potential
// potentials are implemented as vector type.

namespace jam2 {

using namespace std;

bool RQMDv::firstCall=true;

RQMDv::RQMDv(Pythia8::Settings* s) : MeanField(s)
{
  optBaryonCurrent = settings->mode("MeanField:optBaryonCurrent"); 

  double alpha,beta;
  double C1=0.0, C2=0.0, mu1=1.0, mu2=1.0;
  withMomDep=1;

  if(eosType==1) {     // Skyrme Hard Nara2019 K=380 
    alpha = -0.124838736062805;
    beta  =  0.0707373967100532;
    gam   =  2.00261915156202;
    withMomDep=0;

  } else if(eosType==2) { // Skyrme Soft Nara2019 K=210
    alpha = -0.310514700691548;
    beta  =  0.256413361338797;
    gam   =  1.20292928384297;
    withMomDep=0;

  } else if(eosType==3) {  // MH2 Skyrme Hard Nara2019 K=380 

    alpha = 0.00877512904142221;
    beta = 0.07048243772253593;
    gam = 1.8907272542475995;
    C1 = -0.3907221709140509;
    C2= 0.3341868667966711;
    mu1 = 2.02*HBARC;
    mu2= 1.0*HBARC;

  } else if(eosType==4) {  // MS2 Skyrme Hard Nara2019 K=210 

    alpha = -0.3085250472810656;
    beta = 0.3804296791388347;
    gam = 1.1142288606816626;
    C1 = -0.38616561269531274;
    C2= 0.34267680471490275;
    mu1 = 2.02*HBARC;
    mu2= 1.0*HBARC;

  } else if(eosType==5) {  // MH3 Skyrme Hard Nara2019 K=380 

    alpha = 0.031029748068250693;
    beta = 0.04755488958383087;
    gam = 2.1354953844773155;
    C1 = -0.20192630016292473;
    C2= 0.05694208245922472;
    mu1 = 2.8*HBARC;
    mu2= 1.0*HBARC;

  // 2019/1/14
  } else if(eosType==6) {  // MS3 Skyrme Soft Nara2019 K=210 
    //alpha = -0.6285006244768873;
    //beta = 0.7070852697952635;
    //gam = 1.0496537865787445;
    //C1 = -0.20192630016292473;
    //C2= 0.05694208245922472;
    //mu1 = 2.8*HBARC;
    //mu2= 1.0*HBARC;

    alpha = -0.8290165881158311;
    beta = 0.907601225767913;
    gam = 1.0386838048982703;
    C1 = -0.20192630016293525;
    C2= 0.05694208245927023;
    mu1 = 2.8*HBARC;
    mu2= 1.0*HBARC;

  } else if(eosType==7) {  // MH4 Skyrme Hard Nara2019 K=380 
    alpha = 0.03953021450755792;
    beta = 0.040627668934799493;
    gam = 2.2928923888385695;
    C1 = -0.17086291974074935;
    C2 = 0.0;
    mu1 = 3.1466990715061636*HBARC;
    mu2 = 1.0*HBARC;

  } else if(eosType==8) {  // MS4 Skyrme Soft Nara2019 K=210 
    alpha = -0.22912691249560302;
    beta = 0.30928479593796043;
    gam = 1.1087616187247877;
    C1 = -0.17086291974074935;
    mu1 = 3.1466990715061636*HBARC;
    C2 = 0.0;
    mu2 = 1.0*HBARC;

  } else {
    cout << "MeanField:EoStype not implemented " << eosType << endl;
    exit(1);
  }

  t1 = 0.5*alpha/rho0;
  t3 = beta/(gam+1.0)/pow(rho0,gam);
  t3f = gam*t3;

  pmu1= mu1*mu1;
  pmu2= mu2*mu2;
  vex1=C1/(2*rho0);
  vex2=C2/(2*rho0);

  double widG = settings->parm("MeanField:gaussWidth");
  facG = 1.0/pow(4.0*M_PI*widG, 1.5);
  wG = 1.0/(4*widG);

  if(firstCall) {
  cout << "# RQMDv mode eosType= "<< eosType  << " rho0= " << rho0
      << " alpha= " << alpha
      << " beta= "<< beta << " gamma= "<< gam
      << " t1= "<< t1 << " t3= "<< t3
      << " vex1= "<< vex1 << " pmu1= "<< pmu1
      << " vex2= "<< vex2 << " pmu2= "<< pmu2
      <<endl;
    firstCall=false;
  }

  //potential = new VectorPotential(settings);
  //t1 = potential->getT1();
  //t3 = potential->getT3();
  //gam = potential->getGam();
  //withMomDep = potential->isMomDep();
}

RQMDv::~RQMDv()
{
  part.clear();
  rho.clear();
  rhog.clear();
  rhog2.clear();
  vmom4.clear();
  force.clear();
  forcer.clear();
  for(int i=0;i<NV;i++) rhom[i].clear();
  rhom.clear();
  JB.clear();
  Vdot.clear();
}

void RQMDv::evolution(list<EventParticle*>& plist, double t, double dt,int step)
{
  part.clear();
  globalTime = t;
  for(auto& i : plist) {
    if(i->isMeanField(t,optPotential)) part.push_back(i);
  }

  NV = part.size();
  rhom.resize(NV);
  for(int i=0;i<NV;i++) rhom[i].resize(NV);

  rhog.assign(NV,0.0);
  rhog2.assign(NV,0.0);
  rho.assign(NV,0.0);
  vmom4.assign(NV,0.0);
  JB.assign(NV,0.0);
  Vdot.assign(NV,0.0);
  force.assign(NV,0.0);
  forcer.assign(NV,0.0);

  // change from kinetic to canonical momenta.
  if(optVdot==1) {
    for(auto& i : part) i->setFree();
  }

  //qmdMatrix();
  //singleParticlePotential();
  //computeForce();
  //computeForceP();

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
    singleParticlePotential();
    for(int i=0;i<NV;i++) {
      part[i]->addP(-part[i]->potv());
      part[i]->setOnShell();
    }

  } else {
    cout << "RQMDv wrong option optVdot = "<< optVdot<<endl;
    exit(1);
  }

}

// Compute single particle potential.
void RQMDv::singleParticlePotential()
{
  for(int i=0; i< NV; i++) {
    /*
    double r = JB[i].m2Calc();
    if(r> 1e-8) {
      rho[i]=sqrt(r);
      rhog[i] = pow(rho[i],gam-1);
    } else {
      rho[i]=0.0;
      rhog[i]=0.0;
    }
    */
    rho[i] = sqrt(max(0.0, JB[i].m2Calc()));
    if(rho[i]>1e-8) rhog[i] = pow(rho[i],gam-1);
  }

  Vec4 vTotal=0.0;
  Vec4 vc = 0.0;
  bool optV=false;
  if(optVdot==1) optV=true;
  if(optV && optVectorPotential==1) {
    for(int i=0; i< NV; i++) {
      double vsky = part[i]->baryon()/3*(t1 + t3*rhog[i]);
      vTotal += vsky * JB[i] + vmom4[i];
    }
    vc = vTotal / NV;
    vc[0]=0.0;
  }

  Vec4 vdott = 0.0;
  for(int i=0; i< NV; i++) {

    if(rho[i] < 1e-15) continue;
    Vec4 vpot0 = part[i]->potv();
    double vsky = part[i]->baryon()/3*(t1 + t3*rhog[i]);

    // four-components of vector potential are fully included.
    if(optVectorPotential==1) {
      part[i]->setPotV( vsky * JB[i] + vmom4[i] -vc);

    // only time component of vector potential is included.
    } else if(optVectorPotential==2) {
      part[i]->setPotV(0,vsky*JB[i][0] + vmom4[i][0] );

    // only time component of vector potential is included in the form
    // of V(rho_B) where rho_B is an invariant baryon density.
    } else {
      part[i]->setPotV(0,vsky*rho[i] + vmom4[i][0] );
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
Vec4 RQMDv::computeEnergy(list<EventParticle*>& plist, int step)
{
  pTot=0.0;   
  for(auto i = plist.begin(); i != plist.end(); ++i) {
    if(optVectorPotential==1 && optVdot==0) {
      Vec4 pk= (*i)->getP() - (*i)->potv();
      double m = (*i)->getMass();
      pTot[0] += sqrt( m*m + pk.pAbs2());
      pTot[1] += (*i)->getP(1);
      pTot[2] += (*i)->getP(2);
      pTot[3] += (*i)->getP(3);
    } else {
      pTot += (*i)->getP();
    }
    pTot[0] += (*i)->potv(0);
  }

  if(step==1) pTot0 = pTot;

  return pTot;
}

// Pre-factor from the baryon current.
Vec4 RQMDv::facV(int i, Vec4& pkin)
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
    return (t1 + t3f*rhog[i])/rho[i] * JB[i];
  }

  double vj = JB[i] * bj;
  double vv = t1 + t3*rhog[i];  // V/rho_B
  double dv = (gam-1.0)*t3*pow(rho[i],gam-3);  // del(V/rho)/rho
  return vj*dv*JB[i] + vv*bj;

}

void RQMDv::qmdMatrix()
{
  for(int i=0; i< NV; i++) {
    Vec4 r1  = part[i]->getR();
    Vec4 p1 = optPotentialArg>=1? part[i]->getPcan(optVdot) : part[i]->getP();
    distance->setRP1(r1,p1);
    Vec4 v1 = p1/p1[0];

    Vec4 fri = optBaryonCurrent ? part[i]->forceR() : 0.0;
    int bi   = part[i]->baryon()/3;
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;

    for(int j=i+1; j< NV; j++) {
      Vec4 r2  = part[j]->getR();
      Vec4 p2 = optPotentialArg>=1 ? part[j]->getPcan(optVdot) : part[j]->getP();
      distance->setRP2(r2,p2);
      Vec4 v2 = p2/p2[0];

      Vec4 frj = optBaryonCurrent ? part[j]->forceR() : 0.0;
      int bj = part[j]->baryon()/3;
      double qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;
      
      distance->density();
      rhom[i][j] = distance->density1*qfac1*qfac2;
      rhom[j][i] = distance->density2*qfac1*qfac2;

      JB[i] += rhom[i][j]*(v2+frj)*bj;
      JB[j] += rhom[j][i]*(v1+fri)*bi;

      if(!withMomDep) continue;
      distance->psq();

      //potential->Vmd(distance->psq1,distance->psq2);
      //double pmom2ij=potential->vmom4i;
      //double pmom2ji=potential->vmom4j;

      double  pmom2ij = vex1/(1.0-distance->psq1/pmu1)+vex2/(1.0-distance->psq1/pmu2);
      double  pmom2ji = vex1/(1.0-distance->psq2/pmu1)+vex2/(1.0-distance->psq2/pmu2);

      vmom4[i] += pmom2ij*rhom[i][j]*v2;
      vmom4[j] += pmom2ji*rhom[j][i]*v1;
    }
  }

}

void RQMDv::computeForceP()
{
  for(auto i=0; i< NV; i++) {
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
    potential->set1(p1,pk1,JB[i],fi,fengi,rho[i],rhog[i],rhog2[i],
                    //part[i]->pots(),rho[i],rhog[i],bar1);
                    part[i]->pots(),rho[i],rhog[i],bar1,gf1[1],gf1[2],gf1[3],gf1[4],gf1[5]);

    for(auto j=i+1; j< NV; j++) {

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
      potential->set2(p2,pk2,JB[j],fj,fengj,rho[j],rhog[j],rhog2[j],
	  part[j]->pots(),rho[j],rhog[j],bar2,gf2[1],gf2[2],gf2[3],gf2[4],gf2[5]);

      distance->distanceR();
      potential->dVdns(rhom[i][j],rhom[j][i]);

      double fsky1=potential->fskyi;
      double fsky2=potential->fskyj;

      force[i]  += -fsky1*distance->dr2ijri - fsky2*distance->dr2jirj;
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

      distance->distanceP();
      potential->dVdmd(distanceP->psq1,distanceP->psq2,gf1[6]*gf2[6],gf1[7]*gf2[7],gf1[8],gf2[8],gf1[9],gf2[9]);
      double fmomdi=potential->fmomdi;
      double fmomdj=potential->fmomdj;
      double fmomei=potential->fmomei;
      double fmomej=potential->fmomej;

      force[i]  += -fmomdi*distance->dr2ijri   - fmomdj*distance->dr2jirj;
      force[j]  += -fmomdi*distance->dr2ijrj   - fmomdj*distance->dr2jirj;
      forcer[i] +=  fmomei*distance->dp2ijpi + fmomej*distance->dp2jipi
	          + fmomdi*distance->dr2ijpi + fmomdj*distance->dr2jipi;
      forcer[j] +=  fmomej*distance->dp2jipj + fmomei*distance->dp2ijpj
	          + fmomdi*distance->dr2ijpj + fmomdj*distance->dr2jipj;

      if(optDerivative) {
        potential->dVdpm();
        forcer[i] +=  potential->forceri;
        forcer[j] +=  potential->forcerj;
        // gamma derivative.
        forcer[i] += distance->devgam1*potential->facmom;
        forcer[j] += distance->devgam2*potential->facmom;
      }

    } // end loop over j
  } // end loop over i

}

void RQMDv::computeForce()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = optVdot <= 1 ? part[i]->getPkin() : p1;
    if(optPotentialArg>=1) p1 = part[i]->getPcan(optVdot);
    int bar1= part[i]->baryon()/3;

    distance->setRP1(r1,p1);
    Vec4 v1 = p1/p1[0];

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      Vec4 pk2 = optVdot <= 1 ? part[j]->getPkin() : p2;
      if(optPotentialArg>=1) p2 = part[j]->getPcan(optVdot);
      int bar2 = part[j]->baryon()/3;

      distance->setRP2(r2,p2);
      Vec4 v2 = p2/p2[0];

      Vec4 Ai = facV(i,p2);
      Vec4 Aj = facV(j,p1);
      double fskyi = -wG*(Ai * pk1)/pk1[0]*rhom[i][j]*bar1*bar2; 
      double fskyj = -wG*(Aj * pk2)/pk2[0]*rhom[j][i]*bar1*bar2; 

      distance->distanceR();

      force[i]  += -fskyi*distance->dr2ijri - fskyj*distance->dr2jirj;
      force[j]  += -fskyi*distance->dr2ijrj - fskyj*distance->dr2jirj;
      forcer[i] +=  fskyi*distance->dr2ijpi + fskyj*distance->dr2jipi;
      forcer[j] +=  fskyi*distance->dr2ijpj + fskyj*distance->dr2jipj;

      // Derivative of p^\mu/m_j term in the vector potential.
      if(optDerivative) {
        Vec4 A1 = facV(i,pk1);
        Vec4 A2 = facV(j,pk2);
	distance->devV(A1,A2);
        forcer[i] += distance->devV1*rhom[j][i]*bar1*bar2;
        forcer[j] += distance->devV2*rhom[i][j]*bar1*bar2;

	distance->devGamma();
	double facsk= -(fskyi+fskyj)/wG;
        forcer[i] += distance->devgam1*facsk;
        forcer[j] += distance->devgam2*facsk;
      }

      if(!withMomDep) continue;

      double vf1=1.0;
      double vf2=1.0;
      if(optVectorPotential==1) {
        vf1 = pk1/pk1[0] * v2;
        vf2 = pk2/pk2[0] * v1;
      }

      distance->distanceP();

      double fmomdi = -wG*devVmd(distance->psq2)*rhom[i][j]*vf1;
      double fmomdj = -wG*devVmd(distance->psq1)*rhom[j][i]*vf2;
      double fmomei =     devVme(distance->psq2)*rhom[i][j]*vf1;
      double fmomej =     devVme(distance->psq1)*rhom[j][i]*vf2;

      force[i]  += -fmomdi*distance->dr2ijri - fmomdj*distance->dr2jiri;
      force[j]  += -fmomdi*distance->dr2ijrj - fmomdj*distance->dr2jirj;
      forcer[i] +=  fmomei*distance->dp2ijpi + fmomej*distance->dp2jipi
	          + fmomdi*distance->dr2ijpi + fmomdj*distance->dr2jipi;
      forcer[j] +=  fmomei*distance->dp2ijpj + fmomej*distance->dp2jipj
	          + fmomdi*distance->dr2ijpj + fmomdj*distance->dr2jipj;

      // Derivative of p_j / p_0j part.
      if(optDerivative) {
        double facm1 = devVmd(psq1)*rhom[j][i];
        double facm2 = devVmd(psq2)*rhom[i][j];
	distance->devV(pk1/pk1[0],pk2/pk2[0]);
        forcer[i] += facm1*distance->devV1;
        forcer[j] += facm2*distance->devV2;

	// gamma derivative.
        double facmom = facm1*vf1 + facm2*vf2;
        forcer[i] += distance->devgam1*facmom;
        forcer[j] += distance->devgam2*facmom;
      }

    } // end loop over j
  } // end loop over i

}

// non-relativistic two-body distance
// (distance at the global cm frame where the total momentum is zero)
void RQMDv::qmdMatrix1()
{
  for(int i=0; i< NV; i++) {
    Vec4 r1  = part[i]->getR();
    Vec4 p1  = part[i]->getP();
    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    Vec4 fri = optBaryonCurrent ? forcer[i] : 0.0;
    int bi   = part[i]->baryon()/3;
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;

    for(int j=i+1; j< NV; j++) {
      Vec4 r2  = part[j]->getR();
      Vec4 p2  = part[j]->getP();
      Vec4 frj = optBaryonCurrent ? forcer[j] : 0.0;
      int bj = part[j]->baryon()/3;
      if(optPotentialArg==1) p2 = part[j]->getPcan(optVdot);
      double qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;

      //double drsq = (r1 - r2).pAbs2();
      double drsq = (r1 - r2).pT2() + pow2(gamCM*(r1[3]-r2[3]));
      rhom[i][j] = gamCM*facG * exp(-drsq*wG)*qfac1*qfac2;
      rhom[j][i] = rhom[i][j];

      JB[i] += rhom[i][j]*(p2/p2[0]+frj)*bj;
      JB[j] += rhom[j][i]*(p1/p1[0]+fri)*bi;

      if(!withMomDep) continue;
      double ps = (p1 - p2).pAbs2() - (1-optPV) * (p1-p2)[0]*(p1-p2)[0];
      double  pmom2ij = vex1/(1.0+ps/pmu1)+vex2/(1.0+ps/pmu2);
      vmom4[i] += pmom2ij*rhom[i][j]*p2/p2[0];
      vmom4[j] += pmom2ij*rhom[j][i]*p1/p1[0];
    }
  }

}

// two-body distance at their center-of-momentum frame.
void RQMDv::qmdMatrix2()
{
  for(int i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    int bi = part[i]->baryon()/3;
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;
    Vec4 fri = optBaryonCurrent ? forcer[i] : 0.0;

    for(int j=i+1; j< NV; j++) {
      Vec4 frj = optBaryonCurrent ? forcer[j] : 0.0;
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      if(optPotentialArg==1) p2 = part[j]->getPcan(optVdot);
      int bj = part[j]->baryon()/3;
      double qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;

      Vec4 dr = r1 - r2; dr[0] = 0.0;  // just in case;
      Vec4 dp = p1 - p2;
      Vec4 P  = p1 + p2;
      double s = P.m2Calc();
      double drcmsq = dr.m2Calc() - pow2(dr * P)/s;
      double den = P[0]/sqrt(s) * facG * exp(drcmsq*wG)*qfac1*qfac2;
      rhom[i][j] = den;
      rhom[j][i] = den;
      JB[i] += rhom[i][j]*(p2/p2[0] + frj)*bj;
      JB[j] += rhom[j][i]*(p1/p1[0] + fri)*bi;

      if(!withMomDep) continue;
      double ps = dp.m2Calc() - optPV * pow2(dp * P)/s;
      double  pmom2ij = vex1/(1.0-ps/pmu1)+vex2/(1.0-ps/pmu2);
      vmom4[i] += pmom2ij*rhom[i][j]*p2/p2[0];
      vmom4[j] += pmom2ij*rhom[j][i]*p1/p1[0];
    }
  }

}

// two-body distance is computed from the rest-frame of particle j.
void RQMDv::qmdMatrix3()
{
  for(int i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    Vec4 v1 = p1/p1.mCalc();
    int bi = part[i]->baryon()/3;
    Vec4 fri=0.0;
    if(optBaryonCurrent) fri=forcer[i];
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;
    for(int j=i+1; j< NV; j++) {
      Vec4 frj=0.0;
      if(optBaryonCurrent) frj=forcer[j];
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      if(optPotentialArg == 1) p2 = part[j]->getPcan(optVdot);
      Vec4 v2 = p2/p2.mCalc();
      int bj = part[j]->baryon()/3;
      double qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;

      Vec4 dr = r1 - r2; dr[0] = 0.0;  // just in case;
      double drsq1 = dr.m2Calc() - pow2(dr*v2);
      rhom[i][j] = facG * exp(drsq1*wG)*qfac1*qfac2;

      double drsq2 = dr.m2Calc() - pow2(dr*v1);
      rhom[j][i] = facG * exp(drsq2*wG)*qfac1*qfac2;

      JB[i] += rhom[i][j]*(v2 + frj)*bj;
      JB[j] += rhom[j][i]*(v1 + fri)*bi;

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

      vmom4[i] += pmom2ij*rhom[i][j]*v2;
      vmom4[j] += pmom2ji*rhom[j][i]*v1;
    }
  }

}

// non-relativistic distance.
void RQMDv::computeForce1()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = optVdot <=1 ? part[i]->getPkin() : p1;
    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    Vec4 v1 = p1/p1[0];
    int  bar1 = part[i]->baryon()/3;

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      Vec4 pk2 = optVdot <=1 ? part[j]->getPkin() : p2;
      if(optPotentialArg==1) p2 = part[j]->getPcan(optVdot);
      Vec4 v2 = p2/p2[0];
      int  bar2 = part[j]->baryon()/3;

      Vec4 A1 = facV(i,pk2);
      Vec4 A2 = facV(j,pk1);
      double fsky = -wG*((A1*v1)*rhom[i][j] + (A2*v2)*rhom[j][i])*bar1*bar2;

      Vec4 dr = r1 - r2;
      dr[3] *= gamCM2;
      force[i] += -2*fsky*dr;
      force[j] +=  2*fsky*dr;

      // Derivative of p^\mu_j/p^0_j in the vector density.
      if(optDerivative) {
        Vec4 Ai = facV(i,pk1);
        Vec4 Aj = facV(j,pk2);
        double fsgam2i = optP0dev*rhom[j][i]*dot3(Aj,v1)/p1[0]*bar1*bar2;
        double fsgam2j = optP0dev*rhom[i][j]*dot3(Ai,v2)/p2[0]*bar1*bar2;
        double fsgam3i = -rhom[j][i]/p1[0];
        double fsgam3j = -rhom[i][j]/p2[0];
        forcer[i] += fsgam2i*v1 + Aj*fsgam3i;
        forcer[j] += fsgam2j*v2 + Ai*fsgam3j;
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
      double fengij = vf1*rhom[i][j] + vf2*rhom[j][i];

      // p derivative term.
      double fmome = -fengij*(vex1/(pmu1*fac1*fac1)+vex2/(pmu2*fac2*fac2));

      // r derivative term.
      double fmomd = -wG*fengij*(vex1/fac1 + vex2/fac2);

      Vec4 dp2pi =  dp -optP0dev*(1-optPV)*dp[0]*p1/p1[0];
      Vec4 dp2pj = -dp -optP0dev*(1-optPV)*dp[0]*p2/p2[0];

      forcer[i] +=  2*fmome*dp2pi;
      forcer[j] +=  2*fmome*dp2pj;

      force[i]  += -2*fmomd*dr;
      force[j]  +=  2*fmomd*dr;


      // Derivative of p_j/p^0_j in the vector potential.
      if(optDerivative) {
        double fm1 = (vex1/fac1+vex2/fac2)*rhom[j][i]/p1[0];
        double fm2 = (vex1/fac1+vex2/fac2)*rhom[i][j]/p2[0];
        vf1 = dot3(pk1/pk1[0], v2);
        vf2 = dot3(pk2/pk2[0], v1);
        forcer[i] += -fm1 * pk2/pk2[0] + fm1*vf2 * v1 * optP0dev;
        forcer[j] += -fm2 * pk1/pk1[0] + fm2*vf1 * v2 * optP0dev;
      }

    }
  }

}

// two-body distance is defined by the c.m. frame of two-particles.
void RQMDv::computeForce2()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = p1;
    if(optVdot <= 1) pk1 = part[i]->getPkin();

    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    Vec4 bi = p1/p1[0];
    int  b1 = part[i]->baryon()/3;
    double emi = part[i]->getMass();

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      Vec4 pk2 = p2;
      if(optVdot <= 1) pk2 = part[j]->getPkin();
      Vec4 A1 = facV(i,pk2);
      Vec4 A2 = facV(j,pk1);

      if(optPotentialArg==1) p2 = part[j]->getPcan(optVdot);
      Vec4 bj = p2/p2[0];
      int  b2 = part[j]->baryon()/3;
      double emj = part[j]->getMass();

      double facsk =(A1 * bi)*b1*b2*rhom[i][j] + (A2 * bj)*b2*b1*rhom[j][i];

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

      double fsky = -wG*facsk;
      force[i]  += -2*fsky*dr2ri;
      force[j]  += -2*fsky*dr2rj;
      forcer[i] +=  2*fsky*dr2pi;
      forcer[j] +=  2*fsky*dr2pj;

      // Derivative of gamma_ij in front of Gaussian.
      if(optDerivative) {
        //double facsk = fsky3 + fsky4;
	double fsgam1 = facsk / s;
	double fsgam2 = optP0dev*facsk*(1.0/pcm[0] - pcm[0]/s);
	forcer[i] += fsgam1*pcm + fsgam2*bi;
	forcer[j] += fsgam1*pcm + fsgam2*bj;

      // Derivative of p^\mu_j/p^0_j in the vector density.
        Vec4 Ai = facV(i,pk1);
        Vec4 Aj = facV(j,pk2);
        double fsgam2i = optP0dev*rhom[j][i]*dot3(Aj,bi)/p1[0];
        double fsgam2j = optP0dev*rhom[i][j]*dot3(Ai,bj)/p2[0];
        double fsgam3i = -rhom[j][i]/p1[0];
        double fsgam3j = -rhom[i][j]/p2[0];
        forcer[i] += (fsgam2i*bi + Aj*fsgam3i)*b1*b2;
        forcer[j] += (fsgam2j*bj + Ai*fsgam3j)*b1*b2;
      }


      if(!withMomDep) continue;

      double vf1=1.0;
      double vf2=1.0;
      if(optVectorPotential==1) {
        vf1 = pk1/pk1[0] * bj;
        vf2 = pk2/pk2[0] * bi;
      }

      Vec4 dp = p1 - p2;
      double psq = dp.m2Calc() - optPV * pow2(dp*pcm)/s;
      double fac1 = 1.0 - psq/pmu1;
      double fac2 = 1.0 - psq/pmu2;

      // p derivative term.
      double fengij = vf1*rhom[i][j] + vf2*rhom[j][i];
      double fmome = -fengij*(vex1/(pmu1*fac1*fac1) + vex2/(pmu2*fac2*fac2));

      // r derivative term.
      double facmom = vex1/fac1 + vex2/fac2;
      double fmomd = -wG*fengij*facmom;
      double pma = pow2((emi*emi - emj*emj)/s);
      Vec4 dp2pi =  dp - optP0dev*dp[0]*bi + optPV * pcm[0]*pma*bbi;
      Vec4 dp2pj = -dp + optP0dev*dp[0]*bj + optPV * pcm[0]*pma*bbj;

      force[i]  += -2*fmomd*dr2ri;
      force[j]  += -2*fmomd*dr2rj;
      forcer[i] +=  2*fmomd*dr2pi + 2*fmome*dp2pi;
      forcer[j] +=  2*fmomd*dr2pj + 2*fmome*dp2pj;

      // derivative of gamma_ij
      if(optDerivative) {
	double fmgam1 = fengij*facmom/s;
	double fmgam2 = optP0dev*fengij*facmom*(1.0/pcm[0] - pcm[0]/s);
        forcer[i] += pcm*fmgam1 + bi*fmgam2;
        forcer[j] += pcm*fmgam1 + bj*fmgam2;

        // Derivative of p_i/p^0_i in the vector density.
        double fm1 = -(vex1/fac1+vex2/fac2)*rhom[j][i]/p1[0];
        double fm2 = -(vex1/fac1+vex2/fac2)*rhom[i][j]/p2[0];
        vf1 = dot3(pk1/pk1[0], bj);
        vf2 = dot3(pk2/pk2[0], bi);
        forcer[i] += fm1 * pk2/pk2[0] - fm1*vf2 * bi * optP0dev;
        forcer[j] += fm2 * pk1/pk1[0] - fm2*vf1 * bj * optP0dev;
      }

    }
  }

}

// two-body distance r_ij is defined by the rest frame of particle j.
void RQMDv::computeForce3org()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = optVdot <= 1 ? part[i]->getPkin() : p1;
    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    double emi = part[i]->getMass();
    Vec4 vi = p1/emi;
    int bar1= part[i]->baryon()/3;

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      Vec4 pk2 = optVdot <= 1 ? part[j]->getPkin() : p2;
      if(optPotentialArg==1) p2 = part[j]->getPcan(optVdot);
      double emj = part[j]->getMass();
      Vec4 vj = p2/emj;
      int bar2 = part[j]->baryon()/3;

      Vec4 Ai = facV(i,pk2);
      Vec4 Aj = facV(j,pk1);
      double fskyi = -wG*(Ai * vi)*rhom[i][j]*bar1*bar2; 
      double fskyj = -wG*(Aj * vj)*rhom[j][i]*bar1*bar2; 

      Vec4 dr = r1 - r2; dr[0]=0.0;
      double rbi = dot3(dr,vi);
      double rbj = dot3(dr,vj);
      Vec4 dr2ri =  2*(dr + rbj*vj);  // R^2_{ij}/dr_i
      Vec4 dr2rj =  2*(dr + rbi*vi);  // R^2_{ji}/dr_i
      Vec4 dr2pi =  2*dr*rbi/emi;     // R^2_{ji}/dp_i
      Vec4 dr2pj =  2*dr*rbj/emj;     // R^2_{ij}/dp_j

      force[i]  += -fskyi*dr2ri - fskyj*dr2rj;
      force[j]  +=  fskyi*dr2ri + fskyj*dr2rj;
      forcer[i] +=  fskyj*dr2pi;
      forcer[j] +=  fskyi*dr2pj;

      // Derivative of p^\mu/m_j term in the vector potential.
      if(optDerivative) {
        Vec4 A1 = facV(i,pk1)*rhom[j][i]/emi*bar1*bar2;
        Vec4 A2 = facV(j,pk2)*rhom[i][j]/emj*bar1*bar2;

        forcer[i] += (optP0dev*A2[0]*p1/p1[0] - A2);
        forcer[j] += (optP0dev*A1[0]*p2/p2[0] - A1);

        //forcer[i] += (optP0dev*A2[0]*vi/p1[0] - A2/emi)*rhom[j][i];
        //forcer[j] += (optP0dev*A1[0]*vj/p2[0] - A1/emj)*rhom[i][j];
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

      double fmomdi = -wG*devVmd(psq2)*rhom[i][j]*vf1;
      double fmomdj = -wG*devVmd(psq1)*rhom[j][i]*vf2;
      double fmomei =     devVme(psq2)*rhom[i][j]*vf1;
      double fmomej =     devVme(psq1)*rhom[j][i]*vf2;

      force[i]  += -fmomdi*dr2ri   - fmomdj*dr2rj;
      force[j]  +=  fmomdi*dr2ri   + fmomdj*dr2rj;
      forcer[i] +=  fmomei*dp2ijpi + fmomej*dp2jipi + fmomdj*dr2pi;
      forcer[j] +=  fmomej*dp2jipj + fmomei*dp2ijpj + fmomdi*dr2pj;

      // Derivative of p_j / m_j part.
      if(optDerivative) {
        double facm1 = devVmd(psq1)*rhom[j][i]/emi;
        double facm2 = devVmd(psq2)*rhom[i][j]/emj;
        forcer[i] += facm1*( optP0dev*p1/p1[0] - pk2/pk2[0] );
        forcer[j] += facm2*( optP0dev*p2/p2[0] - pk1/pk1[0] );
      }

    }
  }

}

// two-body distance r_ij is defined by the rest frame of particle j.
void RQMDv::computeForce3()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = optVdot <= 1 ? part[i]->getPkin() : p1;
    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    //double emi = part[i]->getMass();
    double emi = p1.mCalc();
    Vec4 vi = p1/emi;
    int bar1= part[i]->baryon()/3;

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      Vec4 pk2 = optVdot <= 1 ? part[j]->getPkin() : p2;
      if(optPotentialArg==1) p2 = part[j]->getPcan(optVdot);
      //double emj = part[j]->getMass();
      double emj = p2.mCalc();
      Vec4 vj = p2/emj;
      int bar2 = part[j]->baryon()/3;

      Vec4 Ai = facV(i,p2);
      Vec4 Aj = facV(j,p1);
      double fskyi = -wG*(Ai * p1)/p1[0]*rhom[i][j]*bar1*bar2; 
      double fskyj = -wG*(Aj * p2)/p2[0]*rhom[j][i]*bar1*bar2; 

      Vec4 dr = r1 - r2; dr[0]=0.0;
      double rbi = dot3(dr,vi);
      double rbj = dot3(dr,vj);
      Vec4 dr2ri =  2*(dr + rbj*vj);  // R^2_{ij}/dr_i
      Vec4 dr2rj =  2*(dr + rbi*vi);  // R^2_{ji}/dr_i
      Vec4 dr2pi =  2*dr*rbi/emi;     // R^2_{ji}/dp_i
      Vec4 dr2pj =  2*dr*rbj/emj;     // R^2_{ij}/dp_j

      force[i]  += -fskyi*dr2ri - fskyj*dr2rj;
      force[j]  +=  fskyi*dr2ri + fskyj*dr2rj;
      forcer[i] +=  fskyj*dr2pi;
      forcer[j] +=  fskyi*dr2pj;

      // Derivative of p^\mu/m_j term in the vector potential.
      if(optDerivative) {
        Vec4 A1 = facV(i,pk1)*rhom[j][i]/emi*bar1*bar2;
        Vec4 A2 = facV(j,pk2)*rhom[i][j]/emj*bar1*bar2;

        forcer[i] += (optP0dev*A2[0]*p1/p1[0] - A2);
        forcer[j] += (optP0dev*A1[0]*p2/p2[0] - A1);

        //forcer[i] += (optP0dev*A2[0]*vi/p1[0] - A2/emi)*rhom[j][i];
        //forcer[j] += (optP0dev*A1[0]*vj/p2[0] - A1/emj)*rhom[i][j];
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

      double fmomdi = -wG*devVmd(psq2)*rhom[i][j]*vf1;
      double fmomdj = -wG*devVmd(psq1)*rhom[j][i]*vf2;
      double fmomei =     devVme(psq2)*rhom[i][j]*vf1;
      double fmomej =     devVme(psq1)*rhom[j][i]*vf2;

      force[i]  += -fmomdi*dr2ri   - fmomdj*dr2rj;
      force[j]  +=  fmomdi*dr2ri   + fmomdj*dr2rj;
      forcer[i] +=  fmomei*dp2ijpi + fmomej*dp2jipi + fmomdj*dr2pi;
      forcer[j] +=  fmomej*dp2jipj + fmomei*dp2ijpj + fmomdi*dr2pj;

      // Derivative of p_j / m_j part.
      if(optDerivative) {
        double facm1 = devVmd(psq1)*rhom[j][i]/emi;
        double facm2 = devVmd(psq2)*rhom[i][j]/emj;
        forcer[i] += facm1*( optP0dev*p1/p1[0] - pk2/pk2[0] );
        forcer[j] += facm2*( optP0dev*p2/p2[0] - pk1/pk1[0] );
      }

    }
  }

}

//*********************************************************************
//...Compute time-derivatives of the vector potential.
void RQMDv::computeVdot()
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
   
      Vec4 Aj = (vj * JB[i]) * Bi;
      Vec4 Ai = (vi * JB[j]) * Bj;

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
      Vdot[j] += dotj*fmomdj*vi;

      pdoti=dot3(force[i],dp2ijpi) + dot3(force[j],dp2ijpj);
      pdotj=dot3(force[j],dp2jipj) + dot3(force[i],dp2jipi);

      Vdot[i] += pdoti*devVme(psq2)*rhom[i][j]*vj;
      Vdot[j] += pdotj*devVme(psq1)*rhom[j][i]*vi;

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

// relative distance and its derivatives between p1 and p2
void RQMDv::distanceP1(const Vec4& p1,const Vec4& p2)
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
void RQMDv::distanceP2(const Vec4& p1,const Vec4& p2)
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
void RQMDv::distanceP3(const Vec4& p1,const Vec4& p2)
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


