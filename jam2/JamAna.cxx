// Copyright (C) 2020 Yasushi Nara
// JAM is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include <jam2/JamAna.h>
#include <jam2/collision/TwoBodyInterList.h>
#include <jam2/hadrons/JamStdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>

namespace jam2 {

using namespace std;

inline bool isInside(Vec4& r,double gam) 
{
  //if((abs(r[1]) <= 6.0) && (abs(r[2]) <= 6.0) && (abs(r[3]) <= 1.0)) return true;
  //if((abs(r[1]) <= 6.0) && (abs(r[2]) <= 6.0) && (abs(r[3]) <= 6.0/gam)) return true;
  //if((abs(r[1]) <= 3.0) && (abs(r[2]) <= 3.0) && (abs(r[3]) <= 3.0/gam)) return true;
  if((abs(r[1]) <= 3.0) && (abs(r[2]) <= 3.0) && (abs(r[3]) <= 1.0)) return true;
  return false;
}

CollisionHistory::CollisionHistory(double ecm,double ftime,double dta,double g) :dt(dta),gamCM(g)
{
  sMax = ecm*2.0;
  sMin = 0.3;
  nS=50;
  dS = (sMax - sMin)/nS;
  scollBB.resize(nS,0.0);
  scollMB.resize(nS,0.0);
  scollMM.resize(nS,0.0);
  scollBBar.resize(nS,0.0);
  scollBarBar.resize(nS,0.0);
  scollBBc.resize(nS,0.0);
  scollMBc.resize(nS,0.0);
  scollMMc.resize(nS,0.0);
  scollBBarc.resize(nS,0.0);
  scollBarBarc.resize(nS,0.0);

  mprint = int(ftime/dta);
  BBcoll.assign(mprint,0);
  MBcoll.assign(mprint,0);
  MMcoll.assign(mprint,0);
  BBcolly.assign(mprint,0);
  MBcolly.assign(mprint,0);
  MMcolly.assign(mprint,0);

  BBy.assign(mprint,0);
  BBf.assign(mprint,0);
  MBy.assign(mprint,0);
  MBf.assign(mprint,0);

  BBStry.assign(mprint,0);
  BBStrf.assign(mprint,0);
  MBStry.assign(mprint,0);
  MBStrf.assign(mprint,0);
  DecStry.assign(mprint,0);
  DecStrf.assign(mprint,0);
}

void CollisionHistory::fill(InterList* inter,const vector<EventParticle*>& outgoing)
{
  int np = inter->getNumberOfInComing();
  //if(np==1) return; // decay

  int cltype= 0;

  EventParticle* p1=inter->getParticle(0);
  EventParticle* p2=0;
  if(np>1) {
  p2=inter->getParticle(1);
  //ParticleDataEntry *pd1= p1->getParticleDataEntry();
  //ParticleDataEntry *pd2=p2->getParticleDataEntry();
  //int pid1=p1->getPID();
  //int pid2=p2->getPID();
  //TwoBodyInterList *inter2=dynamic_cast<TwoBodyInterList*>(inter);
  CollisionPair cpair = inter->getCpair();
  double eCM = cpair.getCMenergy();
  //double sigma=cpair.getSigma();
  //double sigel=cpair.getSigmaElastic();
  //double sigab=cpair.getSigAbs();
  bool preHadronA = cpair.qFactor(0) < 1.0 ? true: false;
  bool preHadronB = cpair.qFactor(1) < 1.0 ? true: false;
  bool preH = preHadronA || preHadronB;

  //cout << " preA = " << preHadronA << " preB= "<< preHadronB << " pre= "<< preH<<endl;

  /*
  int pid3=0, pid4=0;
  if(outgoing.size()==2) {
    pid3=outgoing[0]->getPID();
    pid4=outgoing[1]->getPID();
  }
  */

  //int isr=0, isnn=0;
  //if(pid1 !=id_nucl || pid2 !=id_nucl) isr=1;
  //if(pid3==id_nucl && pid4==id_nucl) isnn=1;
  //if(isr*isnn==1) ncollRR2NN++;

  cltype= inter->getCollType();
  int ix =(eCM - sMin) / dS;
  if(ix <0 || ix >= nS) {
    cout << "CollisionHistory::fill srt too large srt= "<< eCM<< " ix= "<< ix << endl;
    return;
  }
  switch(cltype) {
  case 1: scollBB[ix]+=1.0/dS;
          if(preH) scollBBc[ix]+=1.0/dS;
	  break;
  case 2: scollMB[ix]+=1.0/dS;
          if(preH) scollMBc[ix]+=1.0/dS;
	  break;
  case 3: scollMM[ix]+=1.0/dS;
          if(preH) scollMMc[ix]+=1.0/dS;
	  break;
  case 4: scollBBar[ix]+=1.0/dS;
          if(preH) scollBBarc[ix]+=1.0/dS;
	  break;
  case 5: scollBarBar[ix]+=1.0/dS;
          if(preH) scollBarBarc[ix]+=1.0/dS;
	  break;
  default:
    cout << "CollisionHistory:wrong cltype " << cltype<<endl;
    exit(1);
  }
  }

  //int channel=scatt->getChannel();
  //if(channel==Scatter::ELASTIC) nElastic++;
  //else if(channel==Scatter::ABSORPTION) nAbsorb++;

  double ctime=inter->getCollisionOrderTime();
  double itime = int(ctime/dt+0.5);
  if(itime <0 || itime >= mprint) return;

  bool sBaryon=false;
  double ay= 0.0;
  for(auto& p: outgoing) {
    if(p->baryon() >0 && p->strange()!=0) {
      sBaryon=true;
    }
    ay = max(ay,abs(p->getP().rap()));
  }
  if(np==1 && sBaryon) {
    if(ay<0.5) DecStry[itime]++;
    else if(ay>0.5 && ay <1.5) DecStrf[itime]++;
  }

  if(cltype>0) {
  Vec4 r1 =p1->getR();
  Vec4 r2 =p2->getR();
  bool inside1=isInside(r1,gamCM);
  bool inside2=isInside(r2,gamCM);
  //if((abs(r1[1]) <= 2.0) && (abs(r1[2]) <= 2.0) && (abs(r1[3])/gamCM <= 2.0)) inside1=true;
  //if((abs(r2[1]) <= 2.0) && (abs(r2[2]) <= 2.0) && (abs(r2[3])/gamCM <= 2.0)) inside2=true;
  bool inside = (inside1 || inside2) ? true: false;
  switch(cltype) {
  case 1: BBcoll[itime]++;
	  if(inside) BBcolly[itime]++;
	  if(sBaryon) {
	    if(ay <0.5) BBStry[itime]++;
	    else if(ay > 0.5 && ay < 1.5) BBStrf[itime]++;
	  } else {
	    if(ay <0.5) BBy[itime]++;
	    else if(ay > 0.5 && ay < 1.5) BBf[itime]++;
	  }
	  break;
  case 2: MBcoll[itime]++;
	  if(inside) MBcolly[itime]++;
	  if(sBaryon) {
	    if(ay <0.5) MBStry[itime]++;
	    else if(ay > 0.5 && ay < 1.5) MBStrf[itime]++;
	  } else {
	    if(ay <0.5) MBy[itime]++;
	    else if(ay > 0.5 && ay < 1.5) MBf[itime]++;
	  }
	  break;
  case 3: MMcoll[itime]++;
	  if(inside) MMcolly[itime]++;
	  break;
  }
  }


}

void CollisionHistory::print(int nev,string outFile)
{
  ofstream ofs(outFile.c_str());

  ofs << "# eCM BB  MB  MM BBar BBc MBc MMc BBarc"<< endl;
  ofs << "# smin= "<< sMin << " sMax= "<< sMax << " dS= "<< dS << " nS= "<< nS<<endl;

  double w = 1.0/nev;
  for(int ix=0;ix<(int)scollBB.size();ix++) {
    double xx=sMin + dS*ix;
    ofs << scientific << setprecision(4)
      << setw(6) << xx
      << setw(12) << scollBB[ix]*w
      << setw(12) << scollMB[ix]*w
      << setw(12) << scollMM[ix]*w
      << setw(12) << scollBBar[ix]*w
      << setw(12) << scollBBc[ix]*w
      << setw(12) << scollMBc[ix]*w
      << setw(12) << scollMMc[ix]*w
      << setw(12) << scollBBarc[ix]*w
      <<endl;
  }

  ofs.close();

  ofs.open("JAMTimeEvolCollision.dat");

  ofs << "# (1) time (fm/c) (2) BB collision  (3) MB collision  (4) MM collision"<<endl;
  ofs << "# (5) BB collision(inside)  (6) MB collision (iside)  (7) MM collision (iside)"<<endl;
  ofs <<  "# (8) BB(mid) (9) BB(fb) (10) MB(mid) (11) MB(fb)"<< endl;
  ofs <<  "# (12) BBstr(mid) (13) BBstr(fb) (14) MBstr(mid) (15) MBstr(fb) (16) Decstr(mid) (17) Decstr(fb)"<< endl;

  for(int i=0;i<mprint;i++) {
    ofs << setprecision(4)
      << std::fixed
      << setw(6) << i*dt
      << scientific << setprecision(4)
      << setw(12) << BBcoll[i]*w
      << setw(12) << MBcoll[i]*w
      << setw(12) << MMcoll[i]*w
      << setw(12) << BBcolly[i]*w
      << setw(12) << MBcolly[i]*w
      << setw(12) << MMcolly[i]*w
      << setw(12) << BBy[i]*w
      << setw(12) << BBf[i]*w
      << setw(12) << MBy[i]*w
      << setw(12) << MBf[i]*w
      << setw(12) << BBStry[i]*w
      << setw(12) << BBStrf[i]*w
      << setw(12) << MBStry[i]*w
      << setw(12) << MBStrf[i]*w
      << setw(12) << DecStry[i]*w
      << setw(12) << DecStrf[i]*w
      <<endl;
  }

  ofs.close();
}

// analyze the time evolution of particles.
//--------------------------------------------------------------------------------------
AnaTimeDepParticle::AnaTimeDepParticle(double ftime, double dta,double g, double yc)
  : dt(dta), gamCM(g),yCut(yc)
{
  mprint = int(ftime/dta);
  nPrint=0;
  nP=6;

//  npart[0]:nucleon
//  npart[1]:delta
//  npart[2]:non-strange B*
//  npart[3]:Hyperon
//  npart[4]:non-strange meson
//  npart[5]:strange meson

  npart.resize(nP);

  for(int i=0;i<nP;i++) {
    npart[i].assign(mprint,0.0);
  }

}

void AnaTimeDepParticle::ana(int itime, double coltime,std::list<EventParticle*>& plist)
{
  if(itime==-1) {
    itime = int((coltime+1e-10)/dt);
    if(itime <0 || itime >= mprint) return;
  }

  for(auto& i : plist) {
    if(i->getT() > coltime) continue;
    if(i->getMass()<=0.0) continue;
    // baryon or meson?
    int ip=-1;
    if(i->baryon() !=0) {
      if(i->getPID()==id_nucl) ip=0;
      else if(i->getPID()==id_delt) ip=1;
      else if(i->strange()==0) ip=2;
      else ip=3;
    } else {
      if(i->strange()==0) ip=4;
      else ip=5;
    }

    npart[ip][itime] += 1;

  } // end particle loop
}

void AnaTimeDepParticle::fill(double coltime,std::list<EventParticle*>& plist)
{
  int iprint = int((coltime+1e-10)/dt);
  if(iprint<0 || iprint >= mprint) return;

  if(iprint >= nPrint) {
    for(int i=nPrint;i<=iprint;i++) {
      ana(i,coltime,plist);
    }
    nPrint++;
  } else {
  }
}

void AnaTimeDepParticle::print(std::string outfile, const int nevent)
{
  string f1=outfile+".dat";
  ofstream ofs(f1.c_str());

  ofs << "# (1) time (2) nucleon (3) delta (4) B* (5) hyperon (6) non-strange meson (7) strange meson"<<endl;

  for(int i=0;i<mprint;i++) {
    ofs << setprecision(4)
      << std::fixed
      << setw(6) << i*dt
      << std::scientific
      << setw(12) << npart[0][i]/nevent
      << setw(12) << npart[1][i]/nevent
      << setw(12) << npart[2][i]/nevent
      << setw(12) << npart[3][i]/nevent
      << setw(12) << npart[4][i]/nevent
      << setw(12) << npart[5][i]/nevent
      <<endl;
  }
  ofs.close();
}

//--------------------------------------------------------------------------------------
AnaTimeDepFlow::AnaTimeDepFlow(double ftime, double dta,double g, double yc, double yf,double yfmax)
  : dt(dta), gamCM(g),yCut(yc),yCutF(yf),yCutMax(yfmax)
{
  mprint = int(ftime/dta);
  nPrint=0;
  nP=3;

  v1.resize(nP);
  v2.resize(nP);
  mult.resize(nP);

  v1y.resize(nP);
  v2y.resize(nP);
  multy.resize(nP);

  v1pos.resize(nP);
  v1neg.resize(nP);

  v1r.resize(nP);
  v2r.resize(nP);
  multr.resize(nP);

  v1f.resize(nP);
  v2f.resize(nP);
  multf.resize(nP);

  for(int i=0;i<nP;i++) {
    v1[i].assign(mprint,0.0);
    v2[i].assign(mprint,0.0);
    mult[i].assign(mprint,0);

    v1y[i].assign(mprint,0.0);
    v2y[i].assign(mprint,0.0);
    multy[i].assign(mprint,0);

    v1r[i].assign(mprint,0.0);
    v2r[i].assign(mprint,0.0);
    multr[i].assign(mprint,0);

    v1pos[i].assign(mprint,0.0);
    v1neg[i].assign(mprint,0.0);

    v1f[i].assign(mprint,0.0);
    v2f[i].assign(mprint,0.0);
    multf[i].assign(mprint,0);
  }

  xave.assign(mprint,0.0);
  xavey.assign(mprint,0.0);
  xmult.assign(mprint,0);
  xmulty.assign(mprint,0);

  v1Lambda.assign(mprint,0.0);
  multLambda.assign(mprint,0);
  v1LambdaY.assign(mprint,0.0);
  multLambdaY.assign(mprint,0);
  v1LambdaF.assign(mprint,0.0);
  multLambdaF.assign(mprint,0);

}

//double sign(double A){ return (A>0)-(A<0); }

void AnaTimeDepFlow::ana(int itime, double coltime,std::list<EventParticle*>& plist)
{
  if(itime==-1) {
    itime = int((coltime+1e-10)/dt);
    if(itime <0 || itime >= mprint) return;
  }

  for(auto& i : plist) {
    if(i->getT() > coltime) continue;
    if(i->getMass()<=0.0) continue;
    if(i->baryon()< 0) continue;
    if(abs(i->getNColl()) == 1) continue;
    // baryon or meson?
    int ib = i->baryon()> 0 ? 0 : 1;
    // hyperon
    if(ib==0 && i->strange() !=0) ib=2;


    Vec4 p=i->getP();
    double pt2 = p.pAbs2();
    int sgn= (p[3] > 0) ? 1 : ((p[3] < 0) ? -1 : 0);
    double y = p.rap();
    if(abs(y)>yCutMax) continue;

    //double v1p=p[1]*sgn;  // <px>
    double v1p=p[1]*sgn/sqrt(pt2);  // v1
    double v2p=(p[1]*p[1] - p[2]*p[2])/pt2;

    int id= i->getID();
    if(id==3122 || id==3212) {
      v1Lambda[itime] += v1p;
      multLambda[itime] += 1;
      if(abs(y) <= yCut) {
        v1LambdaY[itime] += v1p;
        multLambdaY[itime] += 1;
      } else if(abs(y) >= yCutF) {
        v1LambdaF[itime] += v1p;
        multLambdaF[itime] += 1;
      }
    }

    v1[ib][itime] += v1p;
    v2[ib][itime] += v2p;
    mult[ib][itime] += 1;

    Vec4 r=i->getR();
    if(isInside(r,gamCM)){
      v1r[ib][itime] += v1p;
      v2r[ib][itime] += v2p;
      multr[ib][itime] += 1;
    }

    // mid-rapidity
    if(abs(y) <= yCut) {
      v1y[ib][itime] += v1p;
      v2y[ib][itime] += v2p;
      multy[ib][itime] += 1;
      if(y*p[1]>0) {
        v1pos[ib][itime] += v1p;
      } else {
        v1neg[ib][itime] += v1p;
      }
      if(ib >=0) {
	xavey[itime] += r[1]*sgn;
	xmulty[itime] += 1;
      }
    }

    // forward-backward
    if(abs(y) >= yCutF) {
      v1f[ib][itime] += v1p;
      v2f[ib][itime] += v2p;
      multf[ib][itime] += 1;
    }

    if(abs(i->getNColl()) != 1) {
      if(ib >=0) {
	xave[itime] += r[1]*sgn;
	xmult[itime] += 1;
      }
    }


  } // end particle loop
}

void AnaTimeDepFlow::fill(double coltime,std::list<EventParticle*>& plist)
{
  int iprint = int((coltime+1e-10)/dt);
  if(iprint<0 || iprint >= mprint) return;

  //cout << "nPrint= "<< nPrint << " iprint= "<< iprint<<endl;
  if(iprint >= nPrint) {
    for(int i=nPrint;i<=iprint;i++) {
      ana(i,coltime,plist);
    }
    nPrint++;
  } else {
  }
}

void AnaTimeDepFlow::print(std::string outfile)
{
  string f1=outfile+".dat";
  ofstream ofs(f1.c_str());

  ofs << "# (1) time (2) v1(baryon) (3) v1(meson) (4) v1(hyperon) (5) v2(baryon) (6) v2(meson) (7) v2(hyperon) (8) v1(Lambda)"<<endl;

  for(int i=0;i<mprint;i++) {
    double w1 =  mult[0][i] > 0  ? 1.0/mult[0][i]  : 0.0;
    double w2 =  mult[1][i] > 0  ? 1.0/mult[1][i]  : 0.0;
    double w3 =  mult[2][i] > 0  ? 1.0/mult[2][i]  : 0.0;
    double w4 =  xmult[i]   > 0  ? 1.0/xmult[i]    : 0.0;
    double wl =  multLambda[i] > 0  ? 1.0/multLambda[i]    : 0.0;

    ofs << setprecision(4)
      << std::fixed
      << setw(6) << i*dt
      << std::scientific
      << setw(12) << v1[0][i]*w1
      << setw(12) << v1[1][i]*w2
      << setw(12) << v1[2][i]*w3
      << setw(12) << v2[0][i]*w1
      << setw(12) << v2[1][i]*w2
      << setw(12) << v2[2][i]*w3
      << setw(12) << v1Lambda[i]*wl
      << setw(12) << xave[i]*w4
      <<endl;
  }
  ofs.close();

  string f2=outfile+"y.dat";
  ofstream ofs2(f2.c_str());
  ofs2 << "# (2) v1(baryon) at mid-rap (3) v1(meson) at mid-rap (4)v1(hyperon) at mid-rap";
  ofs2 << "(5) v2(baryon) at mid-rap (6) v2(meson) at mid-rap (7) v2(hypero) at mid-rap (8) v1(Lambda)"<<endl;
  ofs2 << "# ycut "<< yCut << " <|y|< " << yCutMax << endl;

  for(int i=0;i<mprint;i++) {
    double w1y = multy[0][i] > 0 ? 1.0/multy[0][i] : 0.0;
    double w2y = multy[1][i] > 0 ? 1.0/multy[1][i] : 0.0;
    double w3y = multy[2][i] > 0 ? 1.0/multy[2][i] : 0.0;
    double w4y = xmulty[i]   > 0  ? 1.0/xmulty[i]  : 0.0;
    double wly = multLambdaY[i]   > 0  ? 1.0/multLambdaY[i]  : 0.0;
    ofs2 << setprecision(4)
      << std::fixed
      << setw(6) << i*dt
      << std::scientific
      << setw(12) << v1y[0][i]*w1y
      << setw(12) << v1y[1][i]*w2y
      << setw(12) << v1y[2][i]*w3y
      << setw(12) << v2y[0][i]*w1y
      << setw(12) << v2y[1][i]*w2y
      << setw(12) << v2y[2][i]*w3y
      << setw(12) << v1LambdaY[i]*wly
      << setw(12) << v1pos[0][i]*w1y
      << setw(12) << v1neg[0][i]*w1y
      << setw(12) << v1pos[1][i]*w2y
      << setw(12) << v1neg[1][i]*w2y
      << setw(12) << v1pos[2][i]*w3y
      << setw(12) << v1neg[2][i]*w3y
      << setw(12) << xavey[i]*w4y
      <<endl;
  }
  ofs2.close();

  string f3=outfile+"r.dat";
  ofstream ofs3(f3.c_str());
  ofs3 << "# (6) v1(baryon) (7) v1(meson) (8) v2(baryon) (9) v1(meson) "<<endl;

  for(int i=0;i<mprint;i++) {
    double w1r = multr[0][i] > 0 ? 1.0/multr[0][i] : 0.0;
    double w2r = multr[1][i] > 0 ? 1.0/multr[1][i] : 0.0;
    double w3r = multr[2][i] > 0 ? 1.0/multr[2][i] : 0.0;
    ofs3 << setprecision(4)
      << std::fixed
      << setw(6) << i*dt
      << std::scientific
      << setw(12) << v1r[0][i]*w1r
      << setw(12) << v1r[1][i]*w2r
      << setw(12) << v1r[2][i]*w3r
      << setw(12) << v2r[0][i]*w1r
      << setw(12) << v2r[1][i]*w2r
      << setw(12) << v2r[2][i]*w3r
      <<endl;
  }
  ofs3.close();

  // output v1 and v2 at forward-backward rapidity.
  string f4=outfile+"f.dat";
  ofstream ofs4(f4.c_str());
  ofs4 << "# (2) v1(baryon) at fb-rap (3) v1(meson) at fb-rap (4)v1(hyperon) at fd-rap";
  ofs4 << "(5) v2(baryon) at fb-rap (6) v2(meson) at fb-rap (7) v2(hypero) at fb-rap (8) v1(Lambda)"<<endl;
  ofs4 << "# ycut >= "<< yCutF << " ycutmax= "<< yCutMax <<endl;

  for(int i=0;i<mprint;i++) {
    double w1f = multf[0][i] > 0 ? 1.0/multf[0][i] : 0.0;
    double w2f = multf[1][i] > 0 ? 1.0/multf[1][i] : 0.0;
    double w3f = multf[2][i] > 0 ? 1.0/multf[2][i] : 0.0;
    double wlf = multLambdaF[i] > 0 ? 1.0/multLambdaF[i] : 0.0;
    ofs4 << setprecision(4)
      << std::fixed
      << setw(6) << i*dt
      << std::scientific
      << setw(12) << v1f[0][i]*w1f
      << setw(12) << v1f[1][i]*w2f
      << setw(12) << v1f[2][i]*w3f
      << setw(12) << v2f[0][i]*w1f
      << setw(12) << v2f[1][i]*w2f
      << setw(12) << v2f[2][i]*w3f
      << setw(12) << v1LambdaF[i]*wlf
      <<endl;
  }
  ofs4.close();
}

//--------------------------------------------------------------------------------------

AnaTimeDepDensity::AnaTimeDepDensity(double ftime, double dta,double g, double yc)
  : dt(dta), gamCM(g), yCut(yc)
{
  //double widG=0.25;
  double widG=1.0;
  facG = 1.0/pow(2.0*M_PI*widG, 1.5);
  wG = 1.0/(2*widG);

  //gamCM=1.0;

  mprint = int(ftime/dta);
  nPrint=0;
  nP=2;
  rho1.assign(mprint,0.0);
  rho2.assign(mprint,0.0);
  rhos1.assign(mprint,0.0);
  rhos2.assign(mprint,0.0);
  JB.resize(nP);
  JB2.resize(nP);
  JB3.resize(nP);
  for(int i=0;i<nP;i++) {
    JB[i].assign(mprint,0.0);
    JB2[i].assign(mprint,0.0);
    JB3[i].assign(mprint,0.0);
  }
}

// compute time-dependence of the baryon current.
void AnaTimeDepDensity::ana(int itime, double coltime,std::list<EventParticle*>& plist)
{
  //int optPropagate=0;

  if(itime==-1) {
    itime = int((coltime+1e-10)/dt);
    if(itime <0 || itime >= mprint) return;
  }

  //int ncount = 0;

  Vec4 r0=0.0;
  r0[0]=coltime;
  for(auto& i : plist) {

    if(abs(i->getNColl()) == 1) continue;
    if(i->getT() > coltime) continue;
    if(i->getMass()<=0.0) continue;
    int ib = i->baryon()/3;

    if(ib==0) continue;
    if(ib==0) {
      // exclude formed meson.
      if(i->getTf() <= coltime) continue;
      // cont.quark in the meson.
      ib=0.333333;
    }

    //double t0=i->getT();

    Vec4 p=i->getP();
    //if(abs(p.rap()) > yCut) continue;
    //if(abs(p.rap()) > 2.0) continue;

    //Vec4 r = i->propagate(coltime,optPropagate);
    Vec4 r = i->getR();
    Vec4 dr = r0-r;

    double m = p.mCalc(); 
    double drsq = dr.m2Calc() - pow2(dr*p/m);
    double qfac = i->getTf() > coltime ? i->qFactor() : 1.0;
    double den = facG * exp(drsq*wG);
    double den2 = facG * exp(-wG*dr.pAbs2());
    double den3 = gamCM*facG * exp(-wG*(dr.pT2() + pow2(gamCM*dr[3]))  );

    rho1[itime] += p[0]/m*den*ib;
    rhos1[itime] += den*ib*qfac;
    JB[0][itime] += p/m*den*ib*qfac;
    JB2[0][itime] += p/p[0]*den2*ib*qfac;
    JB3[0][itime] += p/p[0]*den3*ib*qfac;


    // exclude pre-hadrons.
    if(i->getTf() > coltime) continue;

    
    if(qfac < 1.0 || i->baryon()==0)  {
      cout << " qfac < 1 ?"<< qfac << " t= "<< r[0]
	<< " tform= "<< i->getTf()
	<< " tcol= "<< coltime
	<< " id= "<< i->getID()
	<<endl;
      exit(1);
    }

    rho2[itime] += p[0]/m*den*ib;
    rhos2[itime] += den*ib;
    JB[1][itime] += p/m*den*ib;
    JB2[1][itime] += p/p[0]*den2*ib;
    JB3[1][itime] += p/p[0]*den3*ib;

    /*
    bool ism = i->isMeanField(coltime,1);
    cout << "id= "<< i->getID() << " y= "<< p.rap()
      << " ncol= "<< i->getNColl()
      << " lastcl= "<< i->lastColl()
      << " ism= "<< ism<<endl;
    ncount++;
    */

  }
  //cout <<"time= "<< coltime <<  "n= "<< ncount << endl;
}

void AnaTimeDepDensity::fill(double coltime,std::list<EventParticle*>& plist)
{
  int iprint = int((coltime+1e-10)/dt);

  if(iprint<0 || iprint >= mprint) return;

  if(iprint >= nPrint) {
  //cout << "coltime= "<< coltime << " iprint= "<<iprint << " nPrint= "<< nPrint<<endl;
    for(int i=nPrint;i<=iprint;i++) {
      ana(i,coltime,plist);
    }
    nPrint++;
  } else {
  }
}

void AnaTimeDepDensity::print(std::string outfile, int nevent)
{
  ofstream ofs(outfile.c_str());

  ofs << "# (1) time (2) baryon density (3) baryon density (formed hadron) "<<endl;
  ofs << "# ycut = "<< yCut <<endl;

  for(int i=0;i<mprint;i++) {
    ofs << setprecision(4)
      << std::fixed
      << setw(6) << i*dt
      << std::scientific
      << setw(12) << JB[0][i].mCalc()/nevent
      << setw(12) << JB[1][i].mCalc()/nevent
      << setw(12) << rho1[i]/nevent
      << setw(12) << rho2[i]/nevent
      << setw(12) << rhos1[i]/nevent
      << setw(12) << rhos2[i]/nevent
      << setw(12) << JB2[0][i].mCalc()/nevent
      << setw(12) << JB2[1][i].mCalc()/nevent
      << setw(12) << JB3[0][i].mCalc()/nevent
      << setw(12) << JB3[1][i].mCalc()/nevent
      <<endl;
  }
  ofs.close();
}

//--------------------------------------------------------------------------------------

AnaOutPutPhaseSpace::AnaOutPutPhaseSpace(double ftime, double dta) : dt(dta)
{
  //double widG=0.25;
  double widG=1.0;
  facG = 1.0/pow(2.0*M_PI*widG, 1.5);
  wG = 1.0/(2*widG);

  mprint = int(ftime/dta);
  nPrint=0;
  dir = "Phs/";
  struct stat st;
  if(stat(dir.c_str(),&st) !=0) mkdir(dir.c_str(),0775);

}

// compute time-dependence of the baryon current.
void AnaOutPutPhaseSpace::ana(int ievent,int itime, double coltime,std::list<EventParticle*>& plist)
{
  if(itime==-1) {
    itime = int((coltime+1e-10)/dt);
    if(itime <0 || itime >= mprint) return;
  }

  //Vec4 r0=0.0;
  //r0[0]=coltime;

  std::stringstream s;
  s << dir << "ph" << setfill('0') << ievent << "t" << setfill('0') << itime << ".dat";
  ofstream ofs(s.str().c_str());
  //ofstream ofs((dir+"ph"+s.str()+".dat").c_str());

  ofs <<"# time = "<< coltime<<endl;

  for(auto& i : plist) {

    //if(abs(i->getNColl()) == 1) continue;
    if(i->getT() > coltime) continue;
    if(i->getMass()<=0.0) continue;
    //int ib = i->baryon()/3;
    int ib = i->baryon();
    double qfac = i->getTf() > coltime ? i->qFactor() : 1.0;


    // exclude meson.
    //if(ib==0) continue;

    //double qfac = i->getTf() > coltime ? i->qFactor() : 1.0;
    //Vec4 r = i->propagate(coltime,optPropagate);

    ofs << setw(3) << i->getStatus()
        << setw(10) << i->getID()
        << setw(4) << int(ib*qfac)
        << setw(15) << scientific << i->getNColl()
        << setw(12) << fixed << i->getMass()
        << scientific
	<< setw(16) << setprecision(8) << i->getPx()
        << setw(16) << setprecision(8) << i->getPy()
        << setw(16) << setprecision(8) << i->getPz()
        << setw(16) << setprecision(8) << i->getPe()
        << setw(16) << setprecision(8) << i->getX()
        << setw(16) << setprecision(8) << i->getY()
        << setw(16) << setprecision(8) << i->getZ()
        << setw(16) << setprecision(8) << i->getT()
        << setw(16) << setprecision(8) << i->getTf()
	<< endl;

  }

  ofs.close();

}

void AnaOutPutPhaseSpace::fill(int iev,double coltime,std::list<EventParticle*>& plist)
{
  int iprint = int((coltime+1e-10)/dt);

  if(iprint<0 || iprint >= mprint) return;

  if(iprint >= nPrint) {
  //cout << "coltime= "<< coltime << " iprint= "<<iprint << " nPrint= "<< nPrint<<endl;
    for(int i=nPrint;i<=iprint;i++) {
      ana(iev,i,coltime,plist);
    }
    nPrint++;
  } else {
  }
}
//--------------------------------------------------------------------------------------

AnaOutPutDensity::AnaOutPutDensity(double ftime, double dta) : dt(dta)
{
  //double widG=0.25;
  double widG=1.0;
  facG = 1.0/pow(2.0*M_PI*widG, 1.5);
  wG = 1.0/(2*widG);

  mprint = int(ftime/dta);
  nPrint=0;
  dir = "Phs/";
  struct stat st;
  if(stat(dir.c_str(),&st) !=0) mkdir(dir.c_str(),0775);

  int nx=50;
  j0.resize(nx);
  jx.resize(nx);
  jy.resize(nx);
  jz.resize(nx);

  for(int i=0;i<nP;i++) {
    j0[i].assign(mprint,0.0);
    jx[i].assign(mprint,0.0);
  }

}

// compute time-dependence of the baryon current.
void AnaOutPutDensity::ana(int ievent,int itime, double coltime,std::list<EventParticle*>& plist)
{
  if(itime==-1) {
    itime = int((coltime+1e-10)/dt);
    if(itime <0 || itime >= mprint) return;
  }

  //Vec4 r0=0.0;
  //r0[0]=coltime;

  std::stringstream s;
  s << dir << "ph" << setfill('0') << ievent << "t" << setfill('0') << itime << ".dat";
  ofstream ofs(s.str().c_str());
  //ofstream ofs((dir+"ph"+s.str()+".dat").c_str());

  ofs <<"# time = "<< coltime<<endl;

  for(auto& i : plist) {

    //if(abs(i->getNColl()) == 1) continue;
    if(i->getT() > coltime) continue;
    if(i->getMass()<=0.0) continue;
    int ib = i->baryon()/3;

    // exclude meson.
    if(ib==0) continue;

    //double qfac = i->getTf() > coltime ? i->qFactor() : 1.0;
    //Vec4 r = i->propagate(coltime,optPropagate);

    ofs << setw(3) << i->getStatus()
        << setw(8) << i->getID()
        << setw(4) << ib
        << setw(12) << fixed << i->getMass()
        << scientific
	<< setw(16) << setprecision(8) << i->getPx()
        << setw(16) << setprecision(8) << i->getPy()
        << setw(16) << setprecision(8) << i->getPz()
        << setw(16) << setprecision(8) << i->getPe()
        << setw(16) << setprecision(8) << i->getX()
        << setw(16) << setprecision(8) << i->getY()
        << setw(16) << setprecision(8) << i->getZ()
        << setw(16) << setprecision(8) << i->getT()
        << setw(16) << setprecision(8) << i->getTf()
	<< endl;

  }

  ofs.close();

}

//--------------------------------------------------------------------------------------
AnaMeanField::AnaMeanField(double ftime, double dta,double g, double yc)
  : dt(dta), gamCM(g), yCut(yc)
{
  mprint = int(ftime/dta);
  nPrint=0;
  potSk.assign(mprint,0.0);
  potSm.assign(mprint,0.0);
  potVk.assign(mprint,0.0);
  potVm.assign(mprint,0.0);
  rhoS.assign(mprint,0.0);
  rhoB.assign(mprint,0.0);
}

// compute time-dependence of the baryon current.
void AnaMeanField::fill(double coltime,list<EventParticle*>& plist)
{
  int itime = int((coltime+1e-10)/dt);
  //int itime = (int)floor(coltime/dt+0.5);
  //int itime = (int)ceil(coltime/dt);
  if(itime <0 || itime >= mprint) return;

  double rhos=0.0, rhob=0.0;
  double vsk=0.0, vsm=0.0;
  Vec4 vk=0.0, vm=0.0;
  int n=0;
  for(auto& i : plist) {
    //if(abs(i->getNColl()) == 1) continue;
    //Vec4 p=i->getP();
    //if(abs(p.rap()) > yCut) continue;
 
    Vec4 r=i->getR();
    if(!isInside(r,gamCM)) continue;

    //cout << i->getID()<<  " mean= "<< i->meanFieldOn()<<endl;
    if(!(i->meanFieldOn())) continue;
    vsk += i->pots() - i->potsm();
    vsm += i->potsm();
    vk  += i->potv() - i->potvm();
    vm  += i->potvm();
    rhos += i->rhos();
    rhob += i->rhob();
    n++;
  }
  if(n==0) return;
  potSk[itime] += vsk/n;
  potSm[itime] += vsm/n;
  potVk[itime] += vk/n;
  potVm[itime] += vm/n;
  rhoS[itime] += rhos/n;
  rhoB[itime] += rhob/n;
}

void AnaMeanField::print(std::string outfile, int nevent)
{
  ofstream ofs(outfile.c_str());

  ofs << "# (1) time (2) scalar density (3) baryon density "<<endl;
  ofs << "# ycut = "<< yCut<< endl;

  for(int i=0;i<mprint;i++) {
    ofs << setprecision(4)
      << std::fixed
      << setw(6) << i*dt
      << std::scientific
      << setw(12) << rhoS[i]/nevent
      << setw(12) << rhoB[i]/nevent
      << setw(12) << potSk[i]/nevent
      << setw(12) << potSm[i]/nevent
      << setw(12) << potVk[i][0]/nevent
      << setw(12) << potVk[i][1]/nevent
      << setw(12) << potVk[i][2]/nevent
      << setw(12) << potVk[i][3]/nevent
      << setw(12) << potVm[i][0]/nevent
      << setw(12) << potVm[i][1]/nevent
      << setw(12) << potVm[i][2]/nevent
      << setw(12) << potVm[i][3]/nevent
      <<endl;
  }
  ofs.close();
}



//...Purpose: compute energy-momentum tensor and hydrodynamics velocity.
//***********************************************************************

ParticleDensity::ParticleDensity(Pythia8::Settings* s)
{
  settings = s;
  optGauss = settings->mode("Hydro:optGaussSmear");  // mstc(146)
  //optGauss = 3;  // mstc(146)
  widG = settings->parm("Hydro:gaussWidth");
  widG2=2*widG*widG; // Gaussian width
  widCof=pow(1.0/(M_PI*widG2),1.5);
  gVolume=4.0/3.0*M_PI*widG*widG*widG;

  //fluid=f;
  //dX=fluid->dx();
  //dY=fluid->dy();
  //dZ=fluid->dz();
  dX=0.3;dY=0.3;dZ=0.3;

  optDens = 0; //mstc137=0;
  optPreHadron=1;// mstc89
}


int ParticleDensity::computeEnergyMomentum(list<EventParticle*>& plist,EventParticle* i1,
	double ctime,int iopt)
{
  rho=0.0;
  rhob=0.0;
  einv=0.0;
  Vec4 u=0.0;
  Vec4 r0=0.0;
  if(i1) r0 = i1->propagate(ctime,optPropagate);
  double cur[4]={0.0}, curb[4]={0.0}, tens[4][4]={0.0};
  int ncount=0;

  // copy particle list.
  //list<EventParticle*> partlist;
  //partlist.insert(partlist.begin(),plist.begin(),plist.end());

    /*
  //CascBox* box= i1 !=0 ? i1->box(): 0;
  if(iopt==1) {
    for(auto& b : box->getNeighbors2()) {
      if(b->getParticles().size()>0) partlist.insert(partlist.begin(),b->getParticles().begin(),b->getParticles().end());
    }
  } else {
    partlist.insert(partlist.begin(),plist.begin(),plist.end());
  }

  if(iopt==1) {
  cout << " plist.size= "<< plist.size()<<endl;
  cin.get();
  }
  */

//....Loop over all particles
  for(auto ip : plist) {
    if(ip ==  i1) continue;
    if(ip->getMass()  < 1e-5) continue;

    // not yet collide
    if(abs(ip->getNColl()) <=  1) continue;  

    //pre-formed hadrons.
    if(optPreHadron == 0 && ip->getTf() > ctime) continue;

    double dt=ctime - (ip)->getT();
    //if(iopt==1) cout << "dt= "<< dt<<endl;
    if(iopt == 1 && dt <  0.0) continue;

    Vec4 r1 = ip->propagate(ctime,optPropagate);

    //int ix,iy,iz;
    //if(!fluid->inside(r1,ix,iy,iz)) continue;

    double bar=ip->baryon()/3.0;
    double facq=1.0;
    double facb=bar;
    if(iopt == 1 && ip->getTf() > ctime) {
	facq = ip->qFactor();
	facb = facq*bar;
    }

    Vec4 pv = ip->getP();
    double den = gaussSmear(pv,r1-r0);
    double bden = den*facb;

    //...Compute current and energy-momentum tensor.
    for(int im=0; im<4; im++) {
      cur[im]  += pv[im]/pv[0]*den;
      curb[im] += pv[im]/pv[0]*bden;
      for(int ik=0; ik<4;ik++) {
        tens[im][ik] += pv[im]*pv[ik]/pv[0]*den;
      }
    }

    ncount++;

  } // End loop over particles.


  if(ncount ==  0)   return 1;
  if(cur[0] <  1e-7) return 1;

  double  cc=cur[0]*cur[0] - cur[1]*cur[1]-cur[2]*cur[2]-cur[3]*cur[3];
  if(cc < 1e-7) return 1;

  // Compute local energy density and baryon density by using EoS.
  /*
  if(optDens == 1) {
    Vec4 uu(tens[1][0], tens[2][0], tens[3][0], tens[0][0]);
    double pre, sva, vel;
    fluid->thermal(uu,curb[0],einv,rhob,pre,sva,vel);
    double tau=1.0;
    double fg = uu[0] + tau*pre;
    u[1] = uu[1]/fg;
    u[2] = uu[2]/fg;
    u[3] = uu[3]/fg;
    u[0]=sqrt(1.0-u[1]*u[1]-u[2]*u[2]-u[3]*u[3]);
    return 0;
  }
  */

  cc=sqrt(cc);
  u[0]=cur[0]/cc;
  u[1]=cur[1]/cc;
  u[2]=cur[2]/cc;
  u[3]=cur[3]/cc;

//     if(mstc(47).eq.2) call landauframe(u,tens)

//...Lorentz invariant Scalar Number Density
  rho=cur[0]*u[0]-cur[1]*u[1]-cur[2]*u[2]-cur[3]*u[3];

//...Lorentz invariant Baryon density
  rhob=curb[0]*u[0]-curb[1]*u[1]-curb[2]*u[2]-curb[3]*u[3];

//...Lorentz invariant pressure and energy density
  int g[4][4]={{-1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,1}};
  int gg[4][4]={{1,1,1,-1},  {1,1,1,-1}, {1,1,1,-1}, {-1,-1,-1,1}};
  int dl[4][4]={{1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}};

  double gam=u[0];
  double vz=u[3]/u[0];
  double gamz=1.0/sqrt(1.0-vz*vz);
  double z[4];
  z[0]=gamz*vz;
  z[3]=gamz;
  z[2]=0.0;
  z[1]=0.0;
  double pfree=0.0;
  double pxx=0.0;
  double pyy=0.0;
  double pzz=0.0;
  double pperp=0.0;
  double plong=0.0;
  for(int i=0;i<4;i++)
  for(int j=0;j<4;j++) {
    double tmp=gg[i][j]*u[i]*u[j];
    einv += tens[i][j]*tmp;                     // energy density
    pfree += - 1.0/3.0*tens[i][j]*(g[i][j]-tmp); // pressure
    double pl=gg[i][j]*z[i]*z[j];
    plong += tens[i][j]*pl;
    pperp -= 0.5*tens[i][j]*(g[i][j]-tmp+pl);

//...Txx,Tyy,Tzz at the local rest frame.
    double ci=-1.0;
    double cj=-1.0;
    if(i > 0) ci=u[i]/(1+gam);
    if(j > 0) cj=u[j]/(1+gam);
    pxx += tens[i][j]*(dl[1][i]+u[1]*ci)*(dl[1][j]+u[1]*cj);
    pyy += tens[i][j]*(dl[2][i]+u[2]*ci)*(dl[2][j]+u[2]*cj);
    pzz += tens[i][j]*(dl[3][i]+u[3]*ci)*(dl[3][j]+u[3]*cj);
  }

  //cout << " einv= "<< einv << " nv= "<< partlist.size()<<endl;
  //partlist.clear();

  return 0;

}

//------------------------------------------------------------
double ParticleDensity::gaussSmear(const Vec4& pv, const Vec4& dr)
{
  static const double xg3[]={0.0, -.77459666924148337703,.77459666924148337703};
  static const double wg3[]={.8888888888888888,.55555555555555,.555555555555555};
  static const int optg=1;

  Vec4 u=0.0;
  double gam=1.0;
  double pe=pv.e();
  if(optGauss == 0 || optGauss == 3)  {
    double  emd= pv.mCalc();
    u=pv/emd;
    gam=pe/emd;
    if(optGauss == 0) {
	return cmDistanceSquare(dr,u) < widG ? 1.0/gVolume: 0.0;
    }
  } else if(optGauss == 2) {
    double pz=pv.pz();
    double emt=sqrt(pe*pe -pz*pz);
    u[3]=pz/emt;
    gam=pe/emt;
  }

  if(optg == 1) {
    double xtra = cmDistanceSquare(dr,u);
    return widCof*gam*exp(xtra/widG2);

  } else {

    Vec4 dr2 = Vec4(dX/2, dY/2, dZ/2, 0.0);
    double den=0.0;
    for(int ig=0;ig<3;ig++)
    for(int jg=0;jg<3;jg++)
    for(int kg=0;kg<3;kg++) {
      Vec4 r3 = dr - xg3[ig]*dr2;
      double xtra=cmDistanceSquare(r3,u);
      den += exp(xtra/widG2)*wg3[ig]*wg3[jg]*wg3[kg]/8;
    }
    return widCof*gam*den;

  }


}

//--------------------------------------------------------------------------------------

} // end namespace jam2
