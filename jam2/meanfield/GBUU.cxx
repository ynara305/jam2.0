#include <jam2/meanfield/GBUU.h>
#include <jam2/meanfield/NonRelPotential.h>

// BUU tranport model with Gaussian density.

namespace jam2 {

using namespace std;

GBUU::GBUU(Pythia8::Settings* s) : MeanField(s)
{ 
  potential = new NonRelPotential(settings);
  rho0=settings->parm("MeanField:rho0");
  double alpha = potential->getAlpha();
  double beta  = potential->getBeta();
  gam = potential->getGam();
  t1 = alpha/rho0;
  t3 = beta/(pow(rho0,gam));
  t3f = gam*t3;
  withMomDep = potential->isMomDep();

  pmu1 = potential->getPmu1(); 
  pmu2 = potential->getPmu2();
  vex1 = potential->getVex1()*2;
  vex2 = potential->getVex2()*2;

  overSample = settings->mode("Cascade:overSample");
  double widG = settings->parm("MeanField:gaussWidth");
  wG = 1.0/(2*widG);
  facG = 1.0/pow(2.0*M_PI*widG, 1.5);

}

void GBUU::evolution(list<EventParticle*>& plist,double t,double dt, int step)
{
  globalTime = t;
  part.clear();
  for(auto i = plist.begin(); i != plist.end(); ++i) {
    if((*i)->isMeanField(t,optPotential)) part.push_back(*i);
  }

  NV = part.size();
  force = new Vec4 [NV];
  forcer = new Vec4 [NV];
  rho = new double [NV];
  rhog = new double [NV];
  vmoms = new double [NV];
  rhom = new double *[NV];
  for(int i=0;i<NV;i++) rhom[i] = new double [NV];

  qmdMatrix();
  singleParticlePotential();
  computeForce();

  for(int i=0; i<NV;i++)  {
    part[i]->updateByForce(force[i], forcer[i],dt);
  }

  if(isDebug) {
    computeEnergy(plist,step);
    double econ=abs(pTot0[0]-pTot[0])/pTot0[0]*100;
    cout << "GBUU time= " << fixed << step*dt
      << " econ= " << fixed << econ << " %"
     << scientific << setw(13) << pTot[1] 
     << scientific << setw(13) << pTot[2]
     << scientific << setw(13) << pTot[3]
     <<endl;
  }

  delete [] rho;
  for(int i=0;i<NV;i++) delete [] rhom[i];
  delete [] rhom;
  delete [] force;
  delete [] vmoms;
  delete [] forcer;
}

void GBUU::singleParticlePotential()
{
  for(int i=0; i< NV; i++) {
    int ibary = part[i]->baryon()/3;
    double vsky = ibary*(t1 + t3*rhog[i])*rho[i];
    part[i]->setPotV(0,vsky + vmoms[i]);
  }
}

// compute single particle potential energy.
Vec4 GBUU::computeEnergy(list<EventParticle*>& plist, int step)
{
  pTot=0.0;   
  for(auto i = plist.begin(); i != plist.end(); ++i) {
    pTot += (*i)->getP();
  }
  if(step==1) pTot0 = pTot/overSample;
  return pTot/overSample;

}

void GBUU::computeForce()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();

    for(auto j=i; j< NV; j++) {

      Vec4 r2 = part[j]->getR();
      double fsky1 = -wG*(t1 + t3f*rhog[i])*rhom[i][j]; 
      double fsky2 = -wG*(t1 + t3f*rhog[j])*rhom[j][i]; 
      force[i] += -2*fsky1*(r1 - r2);
      force[j] += -2*fsky2*(r2 - r1);
      if(!withMomDep) continue;

      Vec4 p2 = part[j]->getP();
      Vec4 dp = p1 - p2;
      double psq = dp.pAbs2();
      double fac1 = 1.0 + psq/pmu1;
      double fac2 = 1.0 + psq/pmu2;
      double fmom = vex1/(pmu1*fac1*fac1) + vex2/(pmu1*fac2*fac2);
      double fmome1 = -rhom[i][j]*fmom;
      double fmome2 = -rhom[j][i]*fmom;
      double facmom=vex1/fac1 + vex2/fac2;
      double fmomd1 = -wG*rhom[i][j]*facmom;
      double fmomd2 = -wG*rhom[j][i]*facmom;
      force[i] += -2*fmomd1*(r1-r2);
      force[j] += -2*fmomd2*(r2-r1);
      forcer[i] += 2*fmome1*dp;
      forcer[j] += -2*fmome2*dp;
    }
  }

}


void GBUU::qmdMatrix()
{
  for(int i=0;i<NV;i++) {
    force[i]= 0.0;
    forcer[i]= 0.0;
    rho[i]=0.0;
    vmoms[i]=0.0;
  }

  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    for(auto j=i; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double drsq = (r1 - r2).pAbs2();
      double den = facG * exp(-drsq*wG)/overSample;
      rhom[i][j] = den;
      rhom[j][i] = den;
      rho[i] += rhom[i][j];
      rho[j] += rhom[j][i];
      if(!withMomDep) continue;
      double ps = (p1 - p2).pAbs2();
      double  pmom2ij = vex1/(1.0+ps/pmu1)+vex2/(1.0+ps/pmu2);
      vmoms[i] += pmom2ij*rhom[i][j];
      vmoms[j] += pmom2ij*rhom[j][i];
    }
  }

  for(int i=0;i<NV;i++) rhog[i] = pow(rho[i],gam-1);

}

} // namespace jam2


