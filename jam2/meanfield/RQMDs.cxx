#include <jam2/meanfield/RQMDs.h>

// Relativistic quantum molecular dynamics with Skyrme potential
// potentials are implemented as scalar type.

namespace jam2 {

using namespace std;

bool RQMDs::firstCall=true;

RQMDs::RQMDs(Pythia8::Settings* s) : MeanField(s)
{
  withMomDep=1;
  double alpha,beta;
  double C1=0.0, C2=0.0, mu1=1.0, mu2=1.0;
  if(eosType==1) {     // Skyrme Hard Nara2019 K=380 
    //alpha = -124.0e-3;
    //beta = 70.5e-3;
    //gam = 2.0;

    alpha = -0.12384212774318967;
    beta = 0.06875425035612921;
    gam = 2.123648826081624;
    withMomDep=0;

  } else if(eosType==2) { // Skyrme Soft Nara2019 K=210
    //alpha = -356.0e-3;
    //beta = 303.0e-3;
    //gam = 1.16667;

    alpha = -0.2621174235726939;
    beta = 0.20645616373308995;
    gam = 1.2652357572720947;
    withMomDep=0;

  } else if(eosType==3) { // MH2 Skyrme + mom.dep Nara2019 K=380
    alpha = 0.00877512904142221;
    beta = 0.07048243772253593;
    gam = 1.8907272542475995;
    C1 = -0.3907221709140509;
    C2= 0.3341868667966711;
    mu1 = 2.02*HBARC;
    mu2= 1.0*HBARC;

  } else if(eosType==4) { // MH2 Skyrme + mom.dep Nara2019 K=210
    alpha = -0.07785801989125717;
    beta = 0.1565487207865633;
    gam = 1.261622285760052;
    C1 = -0.3907221709140509;
    C2= 0.3341868667966711;
    mu1 = 2.02*HBARC;
    mu2= 1.0*HBARC;

  } else if(eosType==5) { // MH3 Skyrme + mom.dep Nara2019 K=380
    alpha = 0.03911636591740648;
    beta = 0.043064133086746524;
    gam = 2.36458860928333;
    C1 = -0.2048087947472817;
    mu1 = 2.8*HBARC;
    C2= 0.04713369448491099;
    mu2= 1.0*HBARC;

  } else if(eosType==6) { // MS3 Skyrme + mom.dep Nara2019 K=210
    alpha = -0.6285006244768873;
    beta = 0.7070852697952635;
    gam = 1.0496537865787445;
    C1 = -0.20192630016292473;
    C2= 0.05694208245922472;
    mu1 = 2.8*HBARC;
    mu2= 1.0*HBARC;

  } else if(eosType==7) { // MH4 Skyrme + mom.dep Nara2019 K=380
    alpha = 0.045166262659641;
    beta = 0.03824971687635796;
    gam = 2.509835588271163;
    C1 = -0.1787443736909779;
    mu1 = 3.0754425376340975*HBARC;

  } else if(eosType==8) { // MS4 Skyrme + mom.dep Nara2019 K=210
    alpha = -0.0896591573281629;
    beta = 0.17249044531651028;
    gam = 1.202859798739601;
    C1 = -0.1787443736909779;
    mu1 = 3.0754425376340975*HBARC;
 
  //--------------------------------------------------------------------------
  // optical potential is defined by sqrt{(m_N+S)^2+p^2} - sqrt{m_N^2 + p^2}
  // m*/m=0.88
  } else if(eosType==11) { // NH2 Skyrme + mom.dep Nara2021 K=380

    // U( Pinf= 1.7 )= 0.05  Einf= 1.0036086114353737 Pinf= 1.7
    // U( 0.685 )= 0.0 Elab= 0.22349429615474214
    // meff=  0.8823864683629444  M*/M= 0.9407105206427979 rho0= 0.168
    alpha = -0.14746618076876553;
    beta = 0.28231864726575806;
    gam = 1.3086087855271107;
    C1 = -0.8500057544566262;
    mu1 = 2.02*HBARC;
    C2= 1.0508877826863514;
    mu2= 1.0*HBARC;

  } else if(eosType==12) { // MS2 Skyrme + mom.dep Nara2021 K=210
    // U( Pinf= 1.7 )= 0.05  Einf= 1.0036086114353737 Pinf= 1.7
    // U( 0.685 )= 0.0 Elab= 0.22349429615474214
    // meff=  0.8823864683629444  M*/M= 0.9407105206427979 rho0= 0.168
    alpha = -1.7403675255660511;
    beta = 1.8746769579585103;
    gam = 1.0351617222744887;
    C1 = -0.8500057544566262;
    C2= 1.0508877826863514;
    mu1 = 2.02*HBARC;
    mu2= 1.0*HBARC;

  // optical potential is defined by sqrt{(m_N+S)^2+p^2} - sqrt{m_N^2 + p^2}
  // m*/m=0.88
  } else if(eosType==13) { // NH5 Skyrme + mom.dep Nara2019 K=380
    alpha = 0.1436791875727167;
    beta = 0.04604753299756481;
    gam = 2.3300535324726885;
    C1 = -0.272173218332104;
    mu1 = 5.077642517924502;

  } else if(eosType==14) { // Skyrme + mom.dep Ohnishi2015 K=272.59
    alpha = -0.208886959978512;
    beta  = 0.284037556284497;
    gam   = 7.0/6.0;
    mu1 = 2.02*HBARC;
    mu2 = 1.0*HBARC;
    C1  =  -0.38314;
    C2  =   0.33741;

  } else if(eosType==15) { // Skyrme  K=550
    alpha = -0.11096186950393147;
    beta = 0.05645200738039593;
    gam = 2.774173096033283;
    withMomDep=0;

  } else if(eosType==16) { // Skyrme + mom.dep K=450
    alpha = 0.050658464813240926;
    beta = 0.03300384907908334;
    gam = 3.0470954207158307;
    C1 = -0.17874437370225157;
    mu1 = 3.0754425374321817;



  } else {
    cout << "RQMDs mode MeanField:EoStype not implemented " << eosType << endl;
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

  // this mode does not have m_i/p^0_i term.
  if(optTwoBodyDistance==3) optDerivative=false;

  if(firstCall) {
  cout << "# RQMDs mode eosType= "<< eosType  << " rho0= " << rho0
      << " alpha= " << alpha
      << " beta= "<< beta << " gamma= "<< gam
      << " vex1 = "<< vex1 << " vex2 = "<< vex2
      << " pmu = "<< pmu1 << " pmu2 = "<< pmu2
      <<endl;
    firstCall=false;
    cout <<"# RQMDs t1= "<< t1<< " t3f= "<< t3f <<endl;
  }
   
}

// This is called before starting time evolution.
void RQMDs::init(std::list<EventParticle*>& plist)
{
  firstEng=true;
  return;

}

void RQMDs::initMatrix(std::list<EventParticle*>& plist,double t,int step)
{
  part.clear();
  eFree=0.0;
  for(auto& i : plist) {
    if(i->isMeanField(t,optPotential)) {
      part.push_back(i);
    }
    else eFree += i->getP();
  }
  if(step==1) eFree0 = eFree;

  NV = part.size();

  /*
  rho.resize(NV);
  rhog.resize(NV);
  vmoms.resize(NV);
  rhom.resize(NV);
  for(int i=0;i<NV;i++) rhom[i].resize(NV);
  */

  //std::fill(rho.begin(),rho.end(),0.0);
  //std::fill(vmoms.begin(),vmoms.end(),0.0);

  rho.assign(NV,0.0);
  vmoms.assign(NV,0.0);
  rhog.resize(NV);
  rhom.resize(NV);
  for(int i=0;i<NV;i++) rhom[i].resize(NV);

  //force.resize(NV);
  //forcer.resize(NV);
  force.assign(NV,0.0);
  forcer.assign(NV,0.0);

  //for(int i=0;i<NV;i++) { force[i]= 0.0; forcer[i]= 0.0; }
}

void RQMDs::clearMatrix()
{
  part.clear();
  rho.clear();
  rhog.clear();
  vmoms.clear();
  for(int i=0;i<NV;i++) rhom[i].clear();
  rhom.clear();
  force.clear();
  forcer.clear();
}



void RQMDs::evolution(list<EventParticle*>& plist, double t, double dt,int step)
{
  globalTime = t;
  initMatrix(plist,t,step);

  // free mass is used in the calculation of potential and force.
  if(optPotentialArg >= 1) {
    for(int i=0;i<NV;i++) part[i]->setPotS(0.0);
  }

  qmdMatrix();

  //if(optTwoBodyDistance == 1)      qmdMatrix1();
  //else if(optTwoBodyDistance == 2) qmdMatrix2();
  //else if(optTwoBodyDistance == 3) qmdMatrix3();

  bool optSet= optPotentialArg == 0 ? true :  false;
  singleParticlePotential(optSet);

  computeForce();

  //if(optTwoBodyDistance == 1)    computeForce1();
  //else if(optTwoBodyDistance==2) computeForce2();
  //else computeForce3();

  // update momentum 
  for(int i=0; i< NV; i++) {
    part[i]->updateByForce(force[i], forcer[i],dt);
  }

  if(!optSet) singleParticlePotential(true);

  // Save initial total energy-momentum.
  if(firstEng) computeEnergy(plist,step);

  if(optRecoverEnergy) RecoverEnergy(plist);
  //if(optRecoverEnergy) RecoverEnergy2(plist);

  if(isDebug > 1) computeEnergy(plist,step);

  //clearMatrix();

}

// compute single particle potential energy.
Vec4 RQMDs::computeEnergy(list<EventParticle*>& plist, int step)
{
  pTot=0.0;   
  for(auto i = plist.begin(); i != plist.end(); ++i) {
    pTot += (*i)->getP();
  }

  if(firstEng) {pTot0 = pTot; firstEng=false;}
  if(isDebug > 1) {
  double econ=abs(pTot0[0]-pTot[0])/pTot0[0]*100;
  cout << "RQMDs time= "
      << " econ= " << fixed << econ << " %"
     << scientific << setw(13) << pTot[1] 
     << scientific << setw(13) << pTot[2]
     << scientific << setw(13) << pTot[3]
     <<endl;
  }

  return pTot;

}

void RQMDs::singleParticlePotential(bool optSet)
{
  for(int i=0; i< NV; i++) {
    rhog[i] = pow(max(0.0,rho[i]),gam-1);
    double vsky = part[i]->baryon()/3*(t1 + t3*rhog[i])*rho[i];
    part[i]->setPotS(vsky + vmoms[i],optSet);
  }
}

void RQMDs::qmdMatrix(double a)
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    p1 *= a;
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;
    //double m1 = part[i]->getEffectiveMass();
    //p1[0]= sqrt(m1*m1 + p1.pAbs2());

    distance->setRP1(r1,p1);
    double fi=p1.mCalc()/p1[0];

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      p2 *=a;
      double qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;
      //double m2 = part[j]->getEffectiveMass();
      //p2[0]= sqrt(m2*m2 + p2.pAbs2());

      distance->setRP2(r2,p2);
      double fj=p2.mCalc()/p2[0];

      distance->density();
      rhom[i][j] = distance->density1*qfac1*qfac2;
      rhom[j][i] = distance->density2*qfac1*qfac2;

      rho[i] += rhom[i][j]*fj;
      rho[j] += rhom[j][i]*fi;

      if(!withMomDep) continue;

      distance->psq();
      double  pmom2ij = vex1/(1.0-distance->psq1/pmu1)+vex2/(1.0-distance->psq1/pmu2);
      double  pmom2ji = vex1/(1.0-distance->psq2/pmu1)+vex2/(1.0-distance->psq2/pmu2);
      vmoms[i] += pmom2ij*rhom[i][j]*fj;
      vmoms[j] += pmom2ji*rhom[j][i]*fi;

    }
  }
}

void RQMDs::computeForce()
{
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

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double fj = p2.mCalc()/p2[0];
      double fengj=fj;
      if(optPotentialArg==3) {
        double meff2 = part[j]->getEffectiveMass();
        fengj = meff2/sqrt(meff2*meff2 + p2.pAbs2());
      }
      distance->setRP2(r2,p2);

      double fskyi = -wG*fengi*(t1 + t3f*rhog[i])*rhom[i][j]*fj; 
      double fskyj = -wG*fengj*(t1 + t3f*rhog[j])*rhom[j][i]*fi; 

      distance->distanceR();
      force[i]  += -fskyi*distance->dr2ijri - fskyj*distance->dr2jiri;
      force[j]  += -fskyi*distance->dr2ijrj - fskyj*distance->dr2jirj;
      forcer[i] +=  fskyi*distance->dr2ijri + fskyj*distance->dr2jipi;
      forcer[j] +=  fskyj*distance->dr2ijpj + fskyi*distance->dr2jipj;;

      if(optDerivative && optP0dev) {
        forcer[i] +=  fskyj*p1/(wG*p1[0]*p1[0]);
        forcer[j] +=  fskyi*p2/(wG*p2[0]*p2[0]);

        // gamma derivative.
        distance->devGamma();
        double facsk= -(fskyi+fskyj)/wG;
        forcer[i] += distance->devgam1*facsk;
        forcer[j] += distance->devgam2*facsk;
      }

      if(!withMomDep) continue;

      distance->distanceP();
      double fmomdi = -wG*fengi*devVmd(distance->psq2)*rhom[i][j]*fj;
      double fmomei =     fengi*devVme(distance->psq2)*rhom[i][j]*fj;
      double fmomdj = -wG*fengj*devVmd(distance->psq1)*rhom[j][i]*fi;
      double fmomej =     fengj*devVme(distance->psq1)*rhom[j][i]*fi;

      force[i]  += -fmomdi*distance->dr2ijri - fmomdj*distance->dr2jirj;
      force[j]  += -fmomdi*distance->dr2ijri - fmomdj*distance->dr2jirj;
      forcer[i] +=  fmomei*distance->dp2ijpi + fmomej*distance->dp2jipi
	          + fmomdi*distance->dr2ijpj + fmomdj*distance->dr2jipj;
      forcer[j] +=  fmomei*distance->dp2ijpj + fmomej*distance->dp2jipj
	          + fmomdi*distance->dr2ijpj + fmomdj*distance->dr2jipj;

      if(optDerivative) {
        forcer[i] +=  fmomdj*p1/(wG*p1[0]*p1[0])*optP0dev;
        forcer[j] +=  fmomdi*p2/(wG*p2[0]*p2[0])*optP0dev;

        // gamma derivative.
        double facmom = -(fmomdi + fmomdj)/wG;
        forcer[i] += distance->devgam1*facmom;
        forcer[j] += distance->devgam2*facmom;
      }


    }
  }

}

// non-relativistic two-body distance
// (distance at the global cm frame where the total momentum is zero)
void RQMDs::qmdMatrix1(double a)
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double m1 = part[i]->getEffectiveMass();
    //double m1 = p1.mCalc();
    //if(optPotentialArg >=1 ) m1 = part[i]->getMass();
    p1 *= a;
    p1[0]= sqrt(m1*m1+p1.pAbs2());
    double fi=m1/p1[0];
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double m2 = part[j]->getEffectiveMass();
      //double m2 = p2.mCalc();
      //if(optPotentialArg >=1 ) m2 = part[j]->getMass();
      p2 *=a;
      p2[0]= sqrt(m2*m2+p2.pAbs2());
      double fj=m2/p2[0];
      double qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;
      //double drsq = (r1 - r2).pAbs2();
      double drsq = (r1 - r2).pT2() + pow2(gamCM*(r1[3]-r2[3]));
      double den = gamCM*facG * exp(-drsq*wG)*qfac1*qfac2;
      rhom[i][j] = den;
      rhom[j][i] = den;
      rho[i] += rhom[i][j]*fj;
      rho[j] += rhom[j][i]*fi;
      if(!withMomDep) continue;
      Vec4 dp = p1 - p2;
      double ps = dp.pAbs2() - (1-optPV)*dp[0]*dp[0];
      double  pmom2ij = vex1/(1.0+ps/pmu1)+vex2/(1.0+ps/pmu2);
      vmoms[i] += pmom2ij*rhom[i][j]*fj;
      vmoms[j] += pmom2ij*rhom[j][i]*fi;
    }
  }

}

// two-body distance at their center-of-momentum frame.
void RQMDs::qmdMatrix2(double a)
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    //double m1 = p1.mCalc();
    double m1 = part[i]->getEffectiveMass();
    //if(optPotentialArg >=1 ) m1 = part[i]->getMass();
    p1 *= a;
    p1[0]= sqrt(m1*m1+p1.pAbs2());
    double fi = m1/p1[0];
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      //double m2 = p2.mCalc();
      double m2 = part[j]->getEffectiveMass();
      //if(optPotentialArg >=1 ) m2 = part[j]->getMass();
      p2 *=a;
      p2[0]= sqrt(m2*m2+p2.pAbs2());
      double fj=m2/p2[0];
      double qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;

      Vec4 dr = r1 - r2; dr[0] = 0.0;  // just in case;
      Vec4 dp = p1 - p2;
      Vec4 P  = p1 + p2;
      double s = P.m2Calc();
      double drcmsq = dr.m2Calc() - pow2(dr*P)/s;
      double g12 = optScalarDensity > 0 ? P[0]/sqrt(s) : 1.0;
      rhom[i][j] = g12 * facG * exp(drcmsq*wG)*qfac1*qfac2;
      rhom[j][i] = rhom[i][j];
      rho[i] += rhom[i][j]*fj;
      rho[j] += rhom[j][i]*fi;

      if(!withMomDep) continue;
      double ps = dp.m2Calc() - optPV * pow2(dp*P)/s;
      double  pmom2ij = vex1/(1.0-ps/pmu1)+vex2/(1.0-ps/pmu2);
      vmoms[i] += pmom2ij*rhom[i][j]*fj;
      vmoms[j] += pmom2ij*rhom[j][i]*fi;
    }
  }

}

// two-body distance measured from the rest-frame of particle j.
void RQMDs::qmdMatrix3(double a)
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    //double m1sq = p1.m2Calc();
    //if(optPotentialArg >=1 ) m1 = part[i]->getMass();
    double m1 = part[i]->getEffectiveMass();
    double m1sq = m1*m1;
    p1 *= a;
    p1[0]= sqrt(m1sq + p1.pAbs2());
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      //double m2sq = p2.m2Calc();
      //if(optPotentialArg >=1 ) m2 = part[j]->getMass();
      double m2 = part[j]->getEffectiveMass();
      double m2sq = m2*m2;
      p2 *=a;
      p2[0]= sqrt(m2sq + p2.pAbs2());
      double qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;

      Vec4 dr = r1 - r2; dr[0] = 0.0;  // just in case;
      double drsq1 = dr.m2Calc() - pow2(dr * p2)/m2sq;
      rhom[i][j] = facG * exp(drsq1*wG)*qfac1*qfac2;
      double drsq2 = dr.m2Calc() - pow2(dr * p1)/m1sq;
      rhom[j][i] = facG * exp(drsq2*wG)*qfac1*qfac2;
      rho[i] += rhom[i][j];
      rho[j] += rhom[j][i];

      if(!withMomDep) continue;

      //Vec4 dp = p1 - p2;
      //double dpsq = dp.m2Calc();
      //double ps1 = dpsq - optPV * pow2(dp * p2)/m2sq;
      //double ps2 = dpsq - optPV * pow2(dp * p1)/m1sq;

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

      vmoms[i] += (vex1/(1.0-ps1/pmu1)+vex2/(1.0-ps1/pmu2))*rhom[i][j];
      vmoms[j] += (vex1/(1.0-ps2/pmu1)+vex2/(1.0-ps2/pmu2))*rhom[j][i];
    }
  }

}


// non-relativistic distance.
void RQMDs::computeForce1()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double emi = p1.mCalc();
    double fi = emi/p1[0];
    double meff1 = part[i]->getEffectiveMass();
    double fengi = meff1/sqrt(meff1*meff1 + p1.pAbs2());

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double emj = p2.mCalc();
      double fj=emj/p2[0];
      double meff2 = part[j]->getEffectiveMass();
      double fengj = meff2/sqrt(meff2*meff2 + p2.pAbs2());

      double fsky1 = fengi*(t1 + t3f*rhog[i])*rhom[i][j]; 
      double fsky2 = fengj*(t1 + t3f*rhog[j])*rhom[j][i]; 
      double fsky = -wG*(fsky1*fj + fsky2*fi);
      Vec4 dr = r1 - r2;
      dr[3] *= gamCM2;
      force[i] += -2*fsky*dr;
      force[j] +=  2*fsky*dr;

      // Derivative of m_j/p^0_j in the scalar density.
      if(optScalarDensity > 0 && optP0dev) {
        forcer[i] += -fsky2 * fi /p1[0] * p1/p1[0];
        forcer[j] += -fsky1 * fj /p2[0] * p2/p2[0];
      }

      if(!withMomDep) continue;

      // optPV=0: distance is (p_i - p_j)^2 4-distance
      Vec4 dp = p1 - p2;
      double psq = dp.pAbs2() - (1-optPV) * dp[0]*dp[0];
      double fac1 = 1.0 + psq/pmu1;
      double fac2 = 1.0 + psq/pmu2;

      // p derivative term.
      double fengij = fengi*fj*rhom[i][j] + fengj*fi*rhom[j][i];
      double fmome = -fengij*(vex1/(pmu1*fac1*fac1)+vex2/(pmu2*fac2*fac2));

      // r derivative term.
      double facmom = vex1/fac1 + vex2/fac2;
      double fmomd = -wG*fengij*facmom;

      force[i]  += -2*fmomd*dr;
      force[j]  +=  2*fmomd*dr;
      forcer[i] +=  2*fmome*( dp -optP0dev*(1-optPV)*dp[0]*p1/p1[0]);
      forcer[j] +=  2*fmome*(-dp -optP0dev*(1-optPV)*dp[0]*p2/p2[0]);

      // Derivative of m_j/p^0_j in the scalar density.
      if(optScalarDensity > 0 && optP0dev) {
        double fs1 = -facmom*fengj * fi / p1[0] * rhom[j][i];
        double fs2 = -facmom*fengi * fj / p2[0] * rhom[i][j];
        forcer[i] += fs1 * p1/p1[0];
        forcer[j] += fs2 * p2/p2[0];
      }
    }
  }

}

// two-body distance is defined by the c.m. frame of two-particles.
void RQMDs::computeForce2()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double emi = p1.mCalc();
    double fi = emi/p1[0];
    double meff1 = part[i]->getEffectiveMass();
    double fengi = meff1/sqrt(meff1*meff1 + p1.pAbs2());

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double emj = p2.mCalc();
      double fj=emj/p2[0];
      double meff2 = part[j]->getEffectiveMass();
      double fengj = meff2/sqrt(meff2*meff2 + p2.pAbs2());

      double fsky1 = fengi*(t1 + t3f*rhog[i])*rhom[i][j]; 
      double fsky2 = fengj*(t1 + t3f*rhog[j])*rhom[j][i]; 
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
      force[j]  += -2*fsky*dr2rj;
      forcer[i] +=  2*fsky*dr2pi;
      forcer[j] +=  2*fsky*dr2pj;

      // Derivative of gamma_ij in front of Gaussian.
      if(optScalarDensity > 0) {
        double facsk = fj * fsky1 + fi * fsky2;
	double fsgam1 = facsk / s;
	double fsgam2 = optP0dev*facsk*(1.0/pcm[0] - pcm[0]/s);
	forcer[i] += fsgam1*pcm + fsgam2*bi;
	forcer[j] += fsgam1*pcm + fsgam2*bj;

        // Derivative of m_j/p^0_j in the scalar density.
        double fs1 = -fsky2 * fi / p1[0] * optP0dev;
        double fs2 = -fsky1 * fj / p2[0] * optP0dev;
        forcer[i] += fs1 * bi;
        forcer[j] += fs2 * bj;
      }

      if(!withMomDep) continue;

      Vec4 dp = p1 - p2;
      double pma = pow2(emi*emi - emj*emj)/s;
      //double psq = dp.m2Calc() - optPV * pow2(dp * pcm)/s;
      double psq = dp.m2Calc() - pma;

      double fac1 = 1.0 - psq/pmu1;
      double fac2 = 1.0 - psq/pmu2;

      // p derivative term.
      double fengij = fengi*fj*rhom[i][j] + fengj*fi*rhom[j][i];
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

      if(optScalarDensity > 0 ) {
	// derivative of gamma_ij
	double fmgam1 = fengij*facmom/s;
	double fmgam2i = optP0dev*fengij*facmom*(1.0/pcm[0] - pcm[0]/s);
	double fmgam2j = fmgam2i;

        // Derivative of m_j/p^0_j in the scalar density.
        fmgam2i += -facmom*fengj * fi / p1[0] * rhom[j][i] * optP0dev;
        fmgam2j += -facmom*fengi * fj / p2[0] * rhom[i][j] * optP0dev;

        forcer[i] += pcm*fmgam1 + bi*fmgam2i;
        forcer[j] += pcm*fmgam1 + bj*fmgam2j;
      }

    }
  }

}

// two-body distance is defined by the rest frame of particle.
void RQMDs::computeForce3()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double emi = p1.mCalc();
    Vec4 vi = p1 / emi;
    double meff1 = part[i]->getEffectiveMass();
    double fengi = meff1/sqrt(meff1*meff1 + p1.pAbs2());

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double emj = p2.mCalc();
      Vec4 vj = p2 / emj;
      double meff2 = part[j]->getEffectiveMass();
      double fengj = meff2/sqrt(meff2*meff2 + p2.pAbs2());

      double fskyi = -wG*fengi*(t1 + t3f*rhog[i])*rhom[i][j]; 
      double fskyj = -wG*fengj*(t1 + t3f*rhog[j])*rhom[j][i]; 

      Vec4 dr = r1 - r2;dr[0]=0.0;
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

      double fmomdi = -wG*fengi*devVmd(psq2)*rhom[i][j];
      double fmomei =     fengi*devVme(psq2)*rhom[i][j];
      double fmomdj = -wG*fengj*devVmd(psq1)*rhom[j][i];
      double fmomej =     fengj*devVme(psq1)*rhom[j][i];

      force[i]  += -fmomdi*dr2ri   - fmomdj*dr2rj;
      force[j]  +=  fmomdi*dr2ri   + fmomdj*dr2rj;
      forcer[i] +=  fmomei*dp2ijpi + fmomej*dp2jipi + fmomdj*dr2pi;
      forcer[j] +=  fmomej*dp2jipj + fmomei*dp2ijpj + fmomdi*dr2pj;

    }
  }

}

double RQMDs::funcEnergy(list<EventParticle*>& plist,double a)
{
  std::fill(rho.begin(),rho.end(),0.0);
  std::fill(vmoms.begin(),vmoms.end(),0.0);
  qmdMatrix(a);
  bool optSet= optPotentialArg == 0 ? true :  false;
  singleParticlePotential(optSet);
  double etot=0.0;
  for(auto i = plist.begin(); i != plist.end(); ++i) {
    double m = (*i)->getMass() + (*i)->pots();
    double pa2 = (*i)->pAbs2();
    double e = sqrt(m*m + a*a*pa2);
    etot += e;
  }
  return etot;
}

void RQMDs::RecoverEnergy(list<EventParticle*>& plist)
{
  // total energy of the system 
  Vec4 ptot=pTot+eFree0;
  double etot0 = ptot[0];

  Vec4 ptotcm=0.0;
  bool boost = false;
  if(ptot.pAbs2() > 1e-9)  {
    boost = true;
    etot0 = ptot.mCalc();
    for(auto i = plist.begin(); i != plist.end(); ++i) ptotcm += (*i)->getP();
    Vec4 ptot4 = 0.0;
    for(auto i = plist.begin(); i != plist.end(); ++i) {
	(*i)->bstback(ptotcm);
	ptot4 += (*i)->getP();
    }
    cout << "RQMDs::RecoverEnergy  ptot4= "<< ptot4 ;
    cout  << " e= "<< etot0 << endl;
    exit(1);
  }

  /*
  double a0=1.0;
  double a1=1.01;
  double f0 = funcEnergy(plist,a0)-etot0;
  double f1 = funcEnergy(plist,a1)-etot0;
  if(f0*f1 > 0.0) {
    cout << " bisec in RQMDs f0= "<< f0 << " f1= "<< f1 <<endl;
    exit(1);
  }
  cout << "f0= "<< f0 << " f1= "<< f1 <<endl;
  cout << " etot0= "<< pTot0[0] << " efree0= "<< eFree0<<endl;
  int itry=0;
  double a;
  double f;
  do {
    a = 0.5*(a0 + a1);
    f = funcEnergy(plist,a) - etot0;
    if(f*f1 > 0.0) {a1 = a; }
    else {a0 = a; }
    if(++itry > 50) {
	cout << " RQMDs bisec does not converge"<<endl;
	exit(1);
    }
    cout << scientific << " a0= "<< a0 << " a1= "<< a1 << " f= "<< f <<endl;
  } while (abs(f) > 1e-5);

  for(auto i = plist.begin(); i != plist.end(); ++i) {
    (*i)->multP(a);
    (*i)->setOnShell();
  }

  cout << " etot0= "<< etot0 << " dif = "<< f <<endl;
  cin.get();
  */


  double a=1.0;
  double etot = funcEnergy(plist,a);
  int itry=0;
  while(abs(etot-etot0) > etot0*1e-3) {
    a = etot0/etot*a;
    etot = funcEnergy(plist,a);
    if(++itry>10) break;
  }
  
  double efinal=0.0;
  for(auto i = plist.begin(); i != plist.end(); ++i) {
    (*i)->multP(a);
    (*i)->setOnShell();
    efinal += (*i)->getPe();
  }

  if(boost) {
    for(auto i = plist.begin(); i != plist.end(); ++i) (*i)->bst(ptotcm);
  }

}

void RQMDs::RecoverEnergy2(list<EventParticle*>& plist)
{
  // total energy of the system 
  Vec4 ptot=pTot+eFree0;
  double etot0 = ptot[0];

  Vec4 ptotcm=0.0;
  bool boost = false;
  if(ptot.pAbs2() > 1e-9)  {
    boost = true;
    etot0 = ptot.mCalc();
    for(auto i = plist.begin(); i != plist.end(); ++i) ptotcm += (*i)->getP();
    Vec4 ptot4 = 0.0;
    for(auto i = plist.begin(); i != plist.end(); ++i) {
	(*i)->bstback(ptotcm);
	ptot4 += (*i)->getP();
    }
    cout << "RQMDs::RecoverEnergy2 Not implemented yet  ptot4= "<< ptot4 ;
    cout  << " e= "<< etot0 << endl;
    exit(1);
  }

  double a=1.0;
  double etot = funcEnergy(plist,a);
  int itry=0;
  while(abs(etot-etot0) > 1e-5) {
    a = etot0/etot*a;

    for(auto i = plist.begin(); i != plist.end(); ++i) {
      (*i)->multP(a);
      (*i)->setOnShell();
    }

    etot = funcEnergy(plist,1.0);
    if(++itry>10) break;
    cout << " a = "<< a << " etot0= " << etot0
     << " etot= "<< etot
     << scientific << " dif = "<< etot-etot0<<endl;
  }
  
  cout << " a = "<< a << " etot0= " << etot0 << " dif = "<< etot-etot0<<endl;

  double efinal=0.0;
  for(auto i = plist.begin(); i != plist.end(); ++i) {
    (*i)->multP(a);
    (*i)->setOnShell();
    efinal += (*i)->getPe();
  }

  cout << "a= "<< a << " ediff = "<< etot0-efinal <<endl;
  cin.get();

  if(boost) {
    for(auto i = plist.begin(); i != plist.end(); ++i) (*i)->bst(ptotcm);
  }

}

// relative distance and its derivatives between p1 and p2
void RQMDs::distanceP1(const Vec4& p1,const Vec4& p2)
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
void RQMDs::distanceP2(const Vec4& p1,const Vec4& p2)
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
void RQMDs::distanceP3(const Vec4& p1,const Vec4& p2)
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


