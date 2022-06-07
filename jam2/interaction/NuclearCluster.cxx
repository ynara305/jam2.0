#include <jam2/interaction/NuclearCluster.h>
#include <jam2/collision/EventNucleus.h>
#include <jam2/hadrons/JamStdlib.h>
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

namespace jam2 {

//#define sign(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )
template <typename T>int sign(T val){
    return (T(0) < val)-(val < T(0)) ;
}

//inline double pow2(const double& x) {return x*x;}
//const double hc=0.19732705;

double NuclearCluster::getMass(int kf)
{
  switch (abs(kf)) {
    case 2112: return 0.93957; // netron
    case 2212: return 0.93827; // proton
    case 3112: return 1.11568; // lambda
    case 3122: return 1.19744; // Sigma-
    case 3212: return 1.19255; // Sigma0
    case 3222: return 1.18937; // Sigma+
    case 3312: return 1.32130; // Xi-
    case 3322: return 1.31490; // Xi0
    case 3334: return 1.67245; // Omega-
    case 4122: return 2.28490; // Lambda_c+
    case 4132: return 2.47030; // Xi_c0
    case 4232: return 2.46560; // Xi_c+
    default:
      cout << "NuclearCluster::getMass no such baryons kf= " << kf<<endl;
      exit(1);
    }
}

// compute relative distance and momentum.
void NuclearCluster::rpCMsq(EventParticle& pa1, EventParticle& pa2,double& rr,double &pp)
{
  Vec4 p1 = pa1.getP();
  Vec4 p2 = pa2.getP();
  //Vec4 x1  = pa1.getR();
  //Vec4 x2  = pa2.getR();
  Vec4 x1  = pa1.getV();
  Vec4 x2  = pa2.getV();

  Vec4 dx  = x1 - x2;
  Vec4 dp  = p1 - p2;
  Vec4 pcm = p1 + p2;
  double s = pcm.m2Calc();
  double psq = dp.m2Calc();
  double dxpcm = dx*pcm;
  double p1sq = p1.m2Calc();
  double p2sq = p2.m2Calc();
  double dm = p1sq - p2sq;

  // relative momentum in the two-body c.m. frame.
  pp = -psq + dm*dm/s;

  // relative coordinate in the two-body c.m. frame and time of coalescence.
  double tc1,tc2;
  // When t_{cm1} >= t_{cm2}, coalescence time = t_{cm1}
  if(dxpcm >= 0.0) {
    double p2pcm=p2*pcm;
    rr = -dx.m2Calc() + dxpcm/p2pcm*(2*dx*p2 - dxpcm*p2sq/p2pcm);
    tc1 = x1.e();
    tc2 = p2.e()*dxpcm/p2pcm + x2.e();

  // When t_{cm1} < t_{cm2}, coalescence time = t_{cm2}
  } else {
    double p1pcm=p1*pcm;
    rr = -dx.m2Calc() + dxpcm/p1pcm*(2*dx*p1 - dxpcm*p1sq/p1pcm);
    tc1 = -p1.e()*dxpcm/p1pcm + x1.e();
    tc2 = x2.e();
  }

  // Set times of coalescence in the global frame.
  pa1.setFormationTime(tc1);
  pa2.setFormationTime(tc2);

  return;

  /*
  RotBstMatrix frame;
  frame.toCMframe(p1, p2);
  Vec4 rc1 = x1; rc1.rotbst(frame);
  Vec4 rc2 = x2; rc2.rotbst(frame);
  Vec4 pc1 = p1; pc1.rotbst(frame);
  Vec4 pc2 = p2; pc2.rotbst(frame);
  double v1=pc1.pz()/pc1.e();
  double v2=pc2.pz()/pc1.e();
  double t0=max(rc1.e(),rc2.e());
  rc1[3] += v1*(t0-rc1.e());
  rc2[3] += v2*(t0-rc2.e());
  rc1[0] = t0;
  rc2[0] = t0;
  double r2cm = (rc1-rc2).pAbs2();
  // go back to the original computational frame.
  frame.invert();
  rc1.rotbst(frame);
  rc2.rotbst(frame);
  double tcol1=rc1.e();
  double tcol2=rc2.e();

  cout << "rr= "<< rr << " pp= "<<pp <<endl;
  cout << "rr= "<< r2cm  <<endl;
  cout << "tc1= "<< tc1 << " tc2= "<< tc2  <<endl;
  cout << "tc1= "<< tcol1 << " tc2= "<< tcol2  <<endl;
  */

}

bool NuclearCluster::clust(EventParticle* p1, EventParticle* p2)
{
  if(p1->getStatus() >10) return false;
  if(p1->baryon() == 0) return false;

  if(p2->getStatus() >10) return false;
  if(p2->baryon() == 0) return false;

  if(p1->getID()*p2->getID() < 0) return false; // exclude B + antiB

  double rr, pp;
  rpCMsq(*p1,*p2,rr,pp);

  if(rr<= R0*R0 && pp <= P0*P0) return true;
  return false;

}

EventParticle* NuclearCluster::setClusterProperty(vector<EventParticle*>& part, vector<int>& num, int& ii,
    int nt,int isave, int& ibar,Vec4& pf)
{

  int iz=0, in=0,iy=0,ism=0, is0=0,isp=0,ixm=0,ix0=0,ig=0;
  double tmax=-100;
  // loop over particles in this cluster.
  for(int i=0;i<nt;i++) {
    int j=num[ii];
    int kfa=abs(part[j]->getID());
    tmax = max(tmax,part[j]->getTf());
    if(kfa==2212) iz++;
    else if(kfa==2112) in++;
    else if(kfa==3122) iy++;
    else if(kfa==3112) ism++;
    else if(kfa==3212) is0++;
    else if(kfa==3222) isp++;
    else if(kfa==3312) ixm++;
    else if(kfa==3322) ix0++;
    else if(kfa==3334) ig++;
    else {
	cout << "FindCluster::funny kfa? kfa= " << kfa
	     << " nt= " << nt
	     << " m= " <<part[j]->getMass()
	     <<endl;
	return 0;
    }
    ii++;
  } // end loop over particle in the cluster

  // exclude normal nucleus which is not in the mass table.
  if(nt==1) return 0;
  if(iz==0) return 0;
  if(in==0) return 0;
  //if(mastbl)  if(mastbl->inuc(iz,in,0)==0) return 0;


  Vec4 ptot=0.0;
  Vec4 rtot=0.0;
  double lang[3]={0.0, 0.0, 0.0};
  int l=isave;
  int ib=0;
  int ncoll=0;
  for(int i=0;i<nt;i++) {
    int j=num[l];
    part[j]->setStatus(30);
    ncoll = max(ncoll,part[j]->getNColl());
    part[j]->updateR(tmax);
    ptot += part[j]->getP();
    rtot += part[j]->getR();
    ib   += part[j]->baryon()/3;

    // Compute angular momentum.
    lang[0] += part[j]->getY() * part[j]->getPz() - part[j]->getZ() * part[j]->getPy();
    lang[1] += part[j]->getZ() * part[j]->getPx() - part[j]->getX() * part[j]->getPz();
    lang[2] += part[j]->getX() * part[j]->getPy() - part[j]->getY() * part[j]->getPx();
    l++;
 }

  ibar += ib;
  int istr = abs(iy)+abs(ism)+abs(is0)+abs(isp)+abs(ixm)+abs(ix0)+abs(ig);
  int idc = (1000000000 + nt*10 + iz*10000 + istr*10000000 )*sign(ib);

  pf += ptot;
  EventParticle* pc = new EventNucleus(idc,iz,in,iy,ism,is0,isp,ixm,ix0,ig);

  rtot /= nt;
  pc->setR(rtot);
  pc->setP(ptot);
  pc->setVertex(rtot);
  pc->setMass(ptot.mCalc());
  pc->setNumberOfColl(ncoll);
  int ll=int( sqrt(lang[0]*lang[0] + lang[1]*lang[1] + lang[2]*lang[2] )/HBARC + 0.5);
  dynamic_cast<EventNucleus*>(pc)->setAngMom(ll);

  return pc;

}

int NuclearCluster::findCluster(std::list<EventParticle*>& plist)
{
  // container for baryons.
  vector<EventParticle*> part;

  // Select baryons
  for(auto& p : plist) {
    if(p->getStatus() >10) continue;
    if(p->baryon() !=0) part.push_back(p);
  }
  int nv=part.size();

  // Save initial total momentum and baryon number.
  Vec4 ptot0=0.0;
  int ibar0=0;
  for(auto& p: part) {
    ptot0 += p->getP();
    ibar0 += p->baryon()/3;
  }

  // Start coalescence 
  vector<int> mscl(nv,1);
  vector<int> num(nv);
  for(int i=0;i<nv;i++) num[i]=i;

  int  nclst=0; // number of cluster
  int  icheck=0;
  // Find nuclear cluster
  for(int i=0;i<nv-1;i++) {
    int j0=icheck+1;
    int i1=num[i];
    for(int j=j0;j<nv;j++) {
      int i2=num[j];
      // Check if two-particles is close enough in phase space.
      if(clust(part[i1],part[i2])) {
	int lp=num[icheck+1];
	num[icheck+1]=i2;
	num[j]=lp;
	icheck++;
	mscl[nclst]++;
      }
    }
    if(icheck == i) {
      nclst++;
      icheck++;
    }
  }
  nclst++;

  Vec4 pf=0.0;
  int ibar=0;
  int nn=0;
  int ii=0;
  // Loop over nuclear clusters.
  for(int ic=0;ic<nclst;ic++) {
    int nt=mscl[ic]; // number of particle in the cluster
    nn += nt;
    int isave=ii;
    EventParticle* pc = setClusterProperty(part, num,ii, nt, isave, ibar,pf);
    if(pc) plist.push_back(pc);

  }  // end cluster loop


  // Check total momentum conservation.
  nv=part.size();
  for(int i=0;i<nv;i++) {
    if(part[i]->getStatus()>10) continue;
    pf += part[i]->getP();
    ibar += part[i]->baryon()/3;
  }

  if(ibar0 != ibar) {
	cout << "(NuclearCluster:) baryon number does not conserve? ibar0= "<< ibar0
	    << " ibar= " << ibar
	    << " nv " << nv
	    << " nn " << nn
	    << endl;
  }

  if(abs(pf.pAbs()-ptot0.pAbs()) >1e-2) {
    cout << "(NuclearCluster:) total momentum does not conserve" << endl;
    cout << " nclust= " << nclst << endl;
    cout << " nv= " << nv
         << " size= " << part.size()
         << " nn= " << nn
         <<endl;
    cout << " p0= " << ptot0
         << " pf= " << pf;
      return 0;
  }

  return 1;

}

} // end namespace jam2

