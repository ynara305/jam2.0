#include <jam2/meanfield/MBoxN.h>

namespace jam2 {

void NeighborN::resize(int n) 
{
  int m = box->particleSize();
  rhom.resize(n);
  for(int i=0;i<n;i++) rhom[i].resize(m);

}

MBoxN:: ~MBoxN()
{
  for(auto& i : neighbors) delete i;
  neighbors.clear();
  clearMatrix();
}

bool MBoxN::haveThisSite(MBoxN* box) const
{
  auto first = neighbors.begin();
  auto last = neighbors.end();

  while(first != last) {
    if((*first)->box == box) return true;
    ++first;
  }
  return false;
}

//std::vector<NeighborN>::iterator MBoxN::findNeighbor(MBoxN* box)
NeighborN* MBoxN::findNeighbor(MBoxN* box)
{
  auto&& first = neighbors.begin();
  auto&& last = neighbors.end();

  while(first != last) {
    if((*first)->box == box) return *first;
    ++first;
  }
  cout << "MBoxN error no such neighbor "<<endl;
  exit(1);
}

void MBoxN::clearMatrix()
{
  part.clear();
  rho1.clear();
  rho2.clear();
  rho3.clear();
  rhog2.clear();
  rhog3.clear();
  vmoms.clear();
  force.clear();
  forcer.clear();
  for(int i=0;i<NV;i++) rhom[i].clear();
  rhom.clear();
}

void MBoxN::initBox()
{
  NV = part.size();
  if(NV==0) return;

  rho1.assign(NV,0.0);
  rho2.assign(NV,0.0);
  rho3.assign(NV,0.0);
  vmoms.assign(NV,0.0);

  rhog2.resize(NV);
  rhog3.resize(NV);
  rhom.resize(NV);
  for(int i=0;i<NV;i++) rhom[i].resize(NV);

  //std::fill(rho.begin(),rho.end(),0.0);
  //std::fill(vmoms.begin(),vmoms.end(),0.0);

  force.assign(NV,0.0);
  forcer.assign(NV,0.0);

  for(auto& neighbor : neighbors) {
    neighbor->resize(NV);
  }

}

void MBoxN::qmdMatrix()
{
  if(NV==0) return;


  QMD::qmdMatrix();
  for(auto& neighbor : neighbors) {
    if(neighbor->getAction()) qmdMatrixNeighbor(*neighbor);
  }

}

void MBoxN::computeForce(double dt)
{
  if(NV==0) return;

  QMD::computeForce();

  for(auto& neighbor : neighbors) {
    if(neighbor->getAction()) computeForceNeighbor(*neighbor);
  }

  for(int i=0;i<NV;i++) {
    part[i]->updateByForce(force[i],forcer[i],dt);
  }

}

// non-relativistic two-body distance
// (distance at the global cm frame where the total momentum is zero)
void MBoxN::qmdMatrixNeighbor(NeighborN& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;
  double gf1[3]={1.0, 1.0, 1.0};
  double gf2[3]={1.0, 1.0, 1.0};

  for(int i=0; i<NV; i++) {

    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    int bi = part[i]->baryon()/3;
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;
    double* gfac1=part[i]->facPotential();
    int potid1=part[i]->potentialID();

    if(optPotentialDensity==1) {
      gf1[0]=gfac1[1];
      gf1[1]=gfac1[2];
      gf1[2]=gfac1[4];
    }

    for(int j=0; j< n2; j++) {

      Vec4 r2 = part2[j]->getR();
      Vec4 p2 = part2[j]->getP();
      int bj = part2[j]->baryon()/3;
      double qfac2 = part2[j]->getTf() > globalTime ? part2[j]->qFactor() : 1.0;
      double* gfac2=part2[j]->facPotential();
      int potid2=part2[j]->potentialID();

      // to distinguish nucleon/lambda density.
      double pot1 = potential->matrix(potid1,potid2);
      double pot2 = potential->matrix(potid2,potid1);


      if(optPotentialDensity==1) {
        gf2[0]=gfac2[1];
        gf2[1]=gfac2[2];
        gf2[2]=gfac2[4];
      }

      double wid=2*widG*(gfac1[0]+gfac2[0]);
      //double fac = pow(M_PI*wid, 1.5)*rho0;
      double fac = pow(M_PI*wid, 1.5)*overSample;
      double den = gamCM*exp( -((r1-r2).pT2() + pow2(gamCM*(r1[3]-r2[3])) )/wid )*qfac1*qfac2/fac;

      //double den = gamCM*facG * exp( -wG * ((r1-r2).pT2() + pow2(gamCM*(r1[3]-r2[3])) ) )*qfac1*qfac2;
      //double den = facG * exp(-(r1 - r2).pAbs2()*wG);

      neighbor.setRhom(i,j, den*pot1);
      neighbor.neighbor->setRhom(j,i, den*pot2);

      rho1[i] += den*bj*pot1*gf2[0];
      rho2[i] += den*bj*pot1*gf2[1];
      rho3[i] += den*bj*pot1*gf2[2];
      neighbor.box->addRho1(j,den*bi*pot2*gf1[0]);
      neighbor.box->addRho2(j,den*bi*pot2*gf1[1]);
      neighbor.box->addRho3(j,den*bi*pot2*gf1[2]);

      /*
      rho1[i] += den*bj;
      rho2[i] += den*bj;
      rho3[i] += den*bj;
      neighbor.box->addRho1(j,den*bi);
      neighbor.box->addRho2(j,den*bi);
      neighbor.box->addRho3(j,den*bi);
      */

      if(!withMomDep) continue;
      double ps = (p1-p2).pAbs2();
      double vx1=gfac1[6]*gfac2[6]*vex1;
      double vx2=gfac1[7]*gfac2[7]*vex2;
      double pmom2i = den* ( vx1/(1.0+ps/(gfac1[8]*pmu1)) + vx2/(1.0+ps/(gfac1[9]*pmu2)) );
      double pmom2j = den* ( vx1/(1.0+ps/(gfac2[8]*pmu1)) + vx2/(1.0+ps/(gfac2[9]*pmu2)) );
      //double pmom2i = den*(vex1/(1.0+ps/(gfac1[8]*pmu1))+vex2/(1.0+ps/(gfac1[9]*pmu2)));
      //double pmom2j = den*(vex1/(1.0+ps/(gfac2[8]*pmu1))+vex2/(1.0+ps/(gfac2[9]*pmu2)));
      vmoms[i] += pmom2i*pot1;
      neighbor.box->addVmoms(j, pmom2j*pot2);

    }
  }

}

void MBoxN::computeForceNeighbor(NeighborN& neighbor)
{
  auto& part2 = neighbor.box->getParticle();
  int n2 = part2.size();
  if(n2==0) return;

  //std::vector<Vec4> force0(NV,0.0);

  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    int bar1 = part[i]->baryon()/3;
    double* gf1=part[i]->facPotential();
    for(int j=0; j<n2; j++) {

      Vec4 r2 = part2[j]->getR();
      int bar2 = part2[j]->baryon()/3;
      double* gf2=part2[j]->facPotential();
      double wid=1.0/(2*widG*(gf1[0]+gf2[0]));
      double rhomij=neighbor.rhom[i][j];
      double rhomji=neighbor.neighbor->rhom[j][i];

      double vf1=gf1[1]*gf2[1]*t1;
      double vfi2=gf1[3]*t3*rhog2[i];
      double vfi3=gf1[5]*t5*rhog3[i];
      double vfj2=gf2[3]*t3*neighbor.box->getRhog2(j);
      double vfj3=gf2[5]*t5*neighbor.box->getRhog3(j);
      if(optPotentialDensity==1) {
        vfi2 *= gf2[2];
        vfi3 *= gf2[4];
        vfj2 *= gf1[2];
        vfj3 *= gf1[4];
      }
      double fsky1 = (vf1 + vfi2 + vfi3)*rhomij*bar1*bar2;
      double fsky2 = (vf1 + vfj2 + vfj3)*rhomji*bar1*bar2; 
      double fsky = -wid*(fsky1 + fsky2);
      Vec4 dr = r1 - r2;
      dr[3] *= gamCM2;
      force[i] += -2*fsky*dr;
      neighbor.box->addForce(j,2*fsky*dr);

      if(!withMomDep) continue;

      Vec4 p2 = part2[j]->getP();
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

      double fmomei = -rhomij*(vx1/(mui1*faci1*faci1)+vx2/(mui2*faci2*faci2));
      double fmomej = -rhomji*(vx1/(muj1*facj1*facj1)+vx2/(muj2*facj2*facj2));
      double fmomdi = -wid*rhomij*(vx1/faci1 + vx2/faci2);
      double fmomdj = -wid*rhomji*(vx1/facj1 + vx2/facj2);
      double fmomd=fmomdi+fmomdj;
      double fmome=fmomei+fmomej;

      force[i]  += -2*fmomd*dr;
      forcer[i] +=  2*fmome*dp;
      neighbor.box->addForce(j,   2*fmomd*dr);
      neighbor.box->addForceR(j, -2*fmome*dp);
    }
  }


}

} // end of namespace jam2
