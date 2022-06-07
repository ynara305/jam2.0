#include <jam2/meanfield/QMD.h>
#include <jam2/hadrons/JamStdlib.h>
#include <jam2/meanfield/NonRelPotential.h>

// Quantum molecular dynamics with Skyrme potential

namespace jam2 {

using namespace std;

QMD::QMD(Pythia8::Settings* s) : MeanField(s)
{
  potential = new NonRelPotential(settings);
  t1  = potential->getT1();
  t3  = potential->getT3();
  t5  = potential->getT5();
  gam = potential->getGam();
  pmu1 = potential->getPmu1(); 
  pmu2 = potential->getPmu2();
  vex1 = potential->getVex1();
  vex2 = potential->getVex2();
  withMomDep = potential->isMomDep();

  rho0=settings->parm("MeanField:rho0");
  widG = settings->parm("MeanField:gaussWidth");
  facG = 1.0/pow(4.0*M_PI*widG, 1.5)/rho0;
  wG = 1.0/(4*widG);
  overSample = settings->mode("Cascade:overSample");

}

void QMD::initMatrix(std::list<EventParticle*>& plist,double gtime)
{
  part.clear();
  // Select particles which feel the potential force.
  for(auto& p : plist) {
    if(p->isMeanField(gtime,optPotential)) {
      part.push_back(p);
      potential->setPotentialParam(p);
    }
  }

  NV = part.size();
  rho1.assign(NV,0.0);
  rho2.assign(NV,0.0);
  rho3.assign(NV,0.0);
  vmoms.assign(NV,0.0);
  rhog2.resize(NV);
  rhog3.resize(NV);
  rhom.resize(NV);
  for(int i=0;i<NV;i++) rhom[i].resize(NV);

  force.assign(NV,0.0);
  forcer.assign(NV,0.0);
}

QMD::~QMD()
{
  part.clear();
  rho1.clear();
  rho2.clear();
  rho3.clear();
  rhog2.clear();
  rhog3.clear();
  vmoms.clear();
  for(int i=0;i<NV;i++) rhom[i].clear();
  rhom.clear();
  force.clear();
  forcer.clear();
  delete potential;
}

void QMD::evolution(list<EventParticle*>& plist, double gtime,double dt,int step)
{
  //std::vector<Vec4> force(NV,0.0);
  //std::vector<Vec4> forcer(NV,0.0);
  //std::vector<double> rho(NV,0.0);
  //std::vector<double> rhog(NV,0.0);
  //std::vector<double> vmoms(NV,0.0);
  //std::vector<std::vector<double> > rhom(NV, std::vector<double>(NV,0.0));

  globalTime = gtime;
  initMatrix(plist,gtime);
  qmdMatrix();
  singleParticlePotential();
  computeForce();

  for(int i=0; i<NV;i++)  {
    part[i]->updateByForce(force[i], forcer[i],dt);
  }

}

// compute single particle potential energy.
void QMD::singleParticlePotential()
{
  for(auto i=0; i< NV; i++) {
    double* pfac=part[i]->facPotential();

    int ibary = part[i]->baryon()/3;
    rhog2[i] = pfac[2]*pow(max(0.0,rho2[i]),pfac[3]-1);
    rhog3[i] = pfac[4]*pow(max(0.0,rho3[i]),pfac[5]-1);
    double  vsky = ibary*(pfac[1]*t1*rho1[i] + t3*rhog2[i]*rho2[i] + t5*rhog3[i]*rho3[i]);

    part[i]->setPotV(0,vsky + vmoms[i]);
    part[i]->setPotVm(0,vmoms[i]);
    part[i]->setRhoS(rho1[i]);
    part[i]->setRhoB(rho2[i]+rho3[i]);
  }

}

// compute total momentum and energy.
Vec4 QMD::computeEnergy(list<EventParticle*>& plist, int step)
{
  //initMatrix(plist,step);

  pTot=0.0;   
  for(auto i = plist.begin(); i != plist.end(); ++i) {
    pTot += (*i)->getP();
    pTot[0] += (*i)->potv(0);
  }
  if(step==1) pTot0 = pTot/overSample;
  return pTot/overSample;

}

void QMD::qmdMatrix()
{
  if(NV==0) return;
  double gf1[3]={1.0, 1.0, 1.0};
  double gf2[3]={1.0, 1.0, 1.0};

  for(int i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    int bi = part[i]->baryon()/3;
    double  qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;
    double* gfac1= part[i]->facPotential();
    int potid1=part[i]->potentialID();

    if(optPotentialDensity==1) {
      gf1[0]=gfac1[1];
      gf1[1]=gfac1[2];
      gf1[2]=gfac1[4];
    }

    for(int j=i+selfInt; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      int bj = part[j]->baryon()/3;
      int potid2=part[j]->potentialID();
      double  qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;
      double* gfac2= part[j]->facPotential();

      if(optPotentialDensity==1) {
        gf2[0]=gfac2[1];
        gf2[1]=gfac2[2];
        gf2[2]=gfac2[4];
      }

      // to distinguish nucleon/lambda density.
      double pot1 = potential->matrix(potid1,potid2);
      double pot2 = potential->matrix(potid2,potid1);

      double wid=2*widG*(gfac1[0]+gfac2[0]);
      double fac = pow(M_PI*wid, 1.5)*overSample;
      double den = gamCM*exp( -((r1-r2).pT2() + pow2(gamCM*(r1[3]-r2[3])) )/wid )*qfac1*qfac2/fac;
 
      rhom[i][j] = den*pot1;
      rhom[j][i] = den*pot2;

      rho1[i] += rhom[i][j]*bj*gf2[0];
      rho2[i] += rhom[i][j]*bj*gf2[1];
      rho3[i] += rhom[i][j]*bj*gf2[2];

      rho1[j] += rhom[j][i]*bi*gf1[0];
      rho2[j] += rhom[j][i]*bi*gf1[1];
      rho3[j] += rhom[j][i]*bi*gf1[2];

      if(!withMomDep) continue;
      double ps = (p1 - p2).pAbs2();
      double vx1=gfac1[6]*gfac2[6]*vex1;
      double vx2=gfac1[7]*gfac2[7]*vex2;
      double pmom2i = den*( vx1/(1.0+ps/(gfac1[8]*pmu1)) + vx2/(1.0+ps/(gfac1[9]*pmu2)) );
      double pmom2j = den*( vx1/(1.0+ps/(gfac2[8]*pmu1)) + vx2/(1.0+ps/(gfac2[9]*pmu2)) );
      vmoms[i] += pmom2i*pot1;
      vmoms[j] += pmom2j*pot2;
    }
  }

}

void QMD::computeForce()
{
  if(NV==0) return;

  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    int bar1= part[i]->baryon()/3;
    double* gf1=part[i]->facPotential();

    for(auto j=i+selfInt; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      int bar2= part[j]->baryon()/3;
      double* gf2=part[j]->facPotential();

      double vf1=gf1[1]*gf2[1]*t1;
      double vfi2=gf1[3]*t3*rhog2[i];
      double vfi3=gf1[5]*t5*rhog3[i];
      double vfj2=gf2[3]*t3*rhog2[j];
      double vfj3=gf2[5]*t5*rhog3[j];
      if(optPotentialDensity==1) {
        vfi2 *= gf2[2];
        vfi3 *= gf2[4];
        vfj2 *= gf1[2];
        vfj3 *= gf1[4];
      }

      double fsky1 = (vf1 + vfi2 + vfi3)*rhom[i][j]*bar1*bar2; 
      double fsky2 = (vf1 + vfj2 + vfj3)*rhom[j][i]*bar1*bar2; 

      //double fsky1 = (gf1[1]*t1*gf2[1] + gf1[3]*t3*rhog2[i] + t5*gf1[5]*rhog3[i])*rhom[i][j]*bar1*bar2; 
      //double fsky2 = (gf2[1]*t1*gf1[1] + gf2[3]*t3*rhog2[j] + t5*gf2[5]*rhog3[j])*rhom[j][i]*bar1*bar2; 


      double wid=1.0/(2*widG*(gf1[0]+gf2[0]));
      double fsky = -wid*(fsky1 + fsky2);
      Vec4 dr = r1 - r2;
      dr[3] *= gamCM2;
      force[i] += -2*fsky*dr;
      force[j] +=  2*fsky*dr;

      if(!withMomDep) continue;

      Vec4 p2 = part[j]->getP();
      Vec4 dp = p1 - p2;
      double psq = dp.pAbs2();

      double mui1=pmu1*gf1[8];
      double mui2=pmu2*gf1[9];
      double faci1 = 1.0 + psq/mui1;
      double faci2 = 1.0 + psq/mui2;
      double muj1=pmu1*gf2[8];
      double muj2=pmu2*gf2[9];
      double facj1 = 1.0 + psq/muj1;
      double facj2 = 1.0 + psq/muj2;

      double vx1=gf1[6]*gf2[6]*vex1;
      double vx2=gf1[7]*gf2[7]*vex2;
      double fmomei = -rhom[i][j]*(vx1/(mui1*faci1*faci1)+vx2/(mui2*faci2*faci2));
      double fmomej = -rhom[j][i]*(vx1/(muj1*facj1*facj1)+vx2/(muj2*facj2*facj2));
      double fmomdi = -wid*rhom[i][j]*(vx1/faci1 + vx2/faci2);
      double fmomdj = -wid*rhom[j][i]*(vx1/facj1 + vx2/facj2);
      force[i]  += -2*(fmomdi+fmomdj)*dr;
      forcer[i] +=  2*(fmomei+fmomej)*dp;
      force[j]  +=  2*(fmomdi+fmomdj)*dr;
      forcer[j] += -2*(fmomei+fmomej)*dp;

      /*
      //double vfac=gf1[6]*gf2[6];
      double mu1=pmu1*gf1[8];
      double mu2=pmu2*gf1[9];
      double fac1 = 1.0 + psq/mu1;
      double fac2 = 1.0 + psq/mu2;
      double rhomij = gf1[6]*gf2[6]*(rhom[i][j] + rhom[j][i]);
      double fmome = -rhomij*(vex1/(mu1*fac1*fac1)+vex2/(mu2*fac2*fac2));
      double fmomd = -wid*rhomij*(vex1/fac1 + vex2/fac2);
      force[i]  += -2*fmomd*dr;
      forcer[i] +=  2*fmome*dp;
      force[j]  +=  2*fmomd*dr;
      forcer[j] += -2*fmome*dp;
      */

    }
  }

}

} // namespace jam2


