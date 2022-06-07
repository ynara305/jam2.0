#include <jam2/meanfield/PotentialType.h>
#include <jam2/hadrons/JamStdlib.h>

namespace jam2 {

bool PotentialType::firstCall=true;

using namespace std;

PotentialType::PotentialType(Pythia8::Settings* s) : settings(s)
{
  optP0dev=settings->mode("MeanField:optP0dev");
  optPV=settings->mode("MeanField:optMomPotential");
  optScalarDensity=settings->mode("MeanField:optScalarDensity");
  optTwoBodyDistance=settings->mode("MeanField:twoBodyDistance");
  eosType=settings->mode("MeanField:EoS");
  optBaryonCurrent = settings->mode("MeanField:optBaryonCurrent"); 
  optVectorPotential=settings->mode("MeanField:optVectorPotential");
  transportModel = settings->mode("MeanField:transportModel"); 

  //optDerivative = settings->mode("MeanField:optDerivative");
  // this mode does not have m_i/p^0_i term.
  //if(optTwoBodyDistance==3) optDerivative=false;

  rho0=settings->parm("MeanField:rho0");
  optPotentialDensity = settings->mode("MeanField:optPotentialDensity");
  optStrangeBaryonPot = settings->mode("MeanField:optStrangeBaryonPotential");
  optLambdaPot = settings->mode("MeanField:optLambdaPotential");


  double widG = settings->parm("MeanField:gaussWidth");
  //facG = 1.0/pow(4.0*M_PI*widG, 1.5);
  if(transportModel==1) 
    wG = 1.0/(4*widG);
  else {
    wG = 1.0/(2*widG);
  }

  facPotL[0]  = settings->parm("MeanField:factorWidthL");
  facPotL[1]  = settings->parm("MeanField:factorAlphaPotentialL");
  facPotL[2]  = settings->parm("MeanField:factorBetaPotentialL");
  facPotL[3] = 1.0;
  facPotL[4]  = settings->parm("MeanField:factorGammaPotentialL");
  facPotL[5] = 5.0/3.0;
  facPotL[6]  = settings->parm("MeanField:factorAttMomPotentialL");
  facPotL[7]  = settings->parm("MeanField:factorRepMomPotentialL");
  facPotL[8]  = 1.0;
  facPotL[9]  = 1.0;

  facPotS[0]  = settings->parm("MeanField:factorWidthS");
  facPotS[1]  = settings->parm("MeanField:factorAlphaPotentialS");
  facPotS[2]  = settings->parm("MeanField:factorBetaPotentialS");
  facPotS[3] = 1.0;
  facPotS[4]  = settings->parm("MeanField:factorGammaPotentialS");
  facPotS[5] = 5.0/3.0;
  facPotS[6]  = settings->parm("MeanField:factorAttMomPotentialS");
  facPotS[7]  = settings->parm("MeanField:factorRepMomPotentialS");
  facPotS[8]  = 1.0;
  facPotS[9]  = 1.0;

}

int PotentialType::isHyperon(EventParticle* p)
{
  if(p->baryon()==0) return 0;
  int pid= p->getPID();

  // non-strange baryons
  if(p->strange()==0) return 0;

  // lambda and sigma and their resonances
  else if(pid == id_lambda || pid == id_lambdas) return 1;

  // lambda and its resonances
  else if(pid == id_sigma || pid == id_sigmas) return 2;

  else if(p->strange() < 0) return 3;


  //int id= p->getID();
  // lambda and sigma0
  //else if(optStrangeBaryonPot==3 && (id==3122 || id==3212)) return true;
  // lambda only
  //else if(optStrangeBaryonPot==4 && id==3122) return true;
  //else return false;

  return 0;
}

void PotentialType::setPotentialParam(EventParticle* p)
{
  // Lambda potential U_Lambda(rho)= a*rho/rho0 + b*(rho/rho0)^4/3 + c*(rho/rho0)^5/3
 
  int islam = isHyperon(p);

  //            wid,alpha,beta,gam,  c1,  gam2,  cs, cv   mus muv
  double fp[10]={1.0, 1.0, 1.0, gam, 0.0, 5./3., 1.0,1.0, 1.0,1.0};
  int idp=1;

  if(islam==0) {
    if(t5 != 1.0) fp[4]=1.0; // three-range Skyrme potential
    p->setPotentialParam(fp);
    p->setPotentialId(idp);
    return;
  }
  idp=3;

  // Sigma potential
  if(optStrangeBaryonPot==1 && islam==2) {
    for(int i=0;i<10;i++) {
      fp[i]=facPotS[i];
    }
    if(t5 == 1.0) fp[4]=0.0; // two-range Skyrme potential
    fp[3]=gam;
    fp[5] = 5.0/3.0;
    p->setPotentialParam(fp);
    p->setPotentialId(idp);
    return;
  }


  // Hyperon potential is the same as nucleon potential multiplied by some factor.
  if(optLambdaPot==0) {
    for(int i=0;i<10;i++) {
      fp[i]=facPotL[i];
    }
    if(t5 == 1.0) fp[4]=0.0; // two-range Skyrme potential
    fp[3]=gam;
    fp[5] = 5.0/3.0;
    p->setPotentialParam(fp);
    p->setPotentialId(idp);
    return;
  }

  cout << "PotentialType::setPotentialParam ";
  cout << " something is wrong optLambaPot= "<< optLambdaPot << " islam= "<< islam<< endl;
  exit(1);

}


// Pre-factor from the baryon current.
Vec4 PotentialType::facV(Vec4& JB, double rho, double rhog, double rhog2,Vec4& vkin)
{
  if(rho < 1e-15) return Vec4(0.0);
  Vec4 bj=0.0;

  // 4-components of the vector potential is fully included.
  if(optVectorPotential==1) {
    bj = vkin;

  // Only time-component V^0 term is included.
  } else if(optVectorPotential==2) {
    bj[0]=1.0;

  // Only time-component V^0 term is included with the form of V(rho_B).
  } else if (optVectorPotential==3) {
    //return (t1 + t3f*rhog)/rho * JB;
    bj[0]= (pFac1*t1 + pFac2*gam1*t3*rhog + t5*pFac3*gam2*rhog2)/rho * JB*vkin;
    return bj;

  // fully non-relativistic
  } else {
    //bj[0]= t1 + t3f*pow(max(0.0,JB[0]),gam-1);
    bj[0]= pFac1*t1 + gam1*pFac2*t3*pow(max(0.0,JB[0]),gam1-1) + t5*gam2*pFac3*pow(max(0.0,JB[0]),gam2-1);
    return bj;
  }

  double vj = JB * bj;
  //double vv = t1 + t3*rhog;  // V/rho_B
  //double dv = (gam-1.0)*t3*pow(rho,gam-3);  // del(V/rho)/rho

  double vv = pFac1*t1 + pFac2*t3*rhog + t5*pFac3*rhog2;  // V/rho_B
  //double dv = (gam1-1.0)*pFac2*t3*pow(rho,gam1-3)+(gam2-1.0)*pFac3*pow(rho,gam2-3);  // del(V/rho)/rho
  double dv = (gam1-1.0)*pFac2*t3*rhog/(rho*rho) + t5*(gam2-1.0)*pFac3*rhog2/(rho*rho);  // del(V/rho)/rho

  return vj*dv*JB + vv*bj;

}

void PotentialType::devV(const Vec4& Ai, const Vec4& Aj)
{
  // p_i/p^0_i * p^0_i / m_i
  if(optTwoBodyDistance == 3) {
    devV1 = (optP0dev*Aj.e()*v1 - Aj)/p01;
    devV2 = (optP0dev*Ai.e()*v2 - Ai)/p02;
  // p_i/p^0_i
  } else {
    devV1 = (optP0dev*dot3(Aj,v1)*v1 - Aj)/p01;
    devV2 = (optP0dev*dot3(Ai,v2)*v2 - Ai)/p02;
  }
}





} // namespace jam2


