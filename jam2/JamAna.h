// Copyright (C) 2020 Yasushi Nara.

#ifndef jam2_JamAna_H
#define jam2_JamAna_H

#include <jam2/collision/EventParticle.h>
#include <jam2/collision/InterList.h>
//#include <jam2/fluid/Fluid.h>
#include <string>
#include <vector>

namespace jam2 {

//using namespace Pythia8;

class JamAna
{
public:
    //ParticleData    *particleData;
    //JamParticleData* jamParticleData;

private:
    //Collision* event;
    //Scatter* scatt;

    //Vec4 pTot, pTot0;
    //std::list<EventParticle*>* particleList;

public:
  // Initialize.
  bool init();

};

class CollisionHistory
{
private:
  int    nS;
  double dS, sMax,sMin;
  std::vector<double> scollBB, scollMB, scollMM, scollBBar, scollBarBar;
  std::vector<double> scollBBc, scollMBc, scollMMc, scollBBarc, scollBarBarc;
  double dt,gamCM; // time slice for time dependent analysis.
  int mprint; // number of time step
  std::vector<int> BBcoll, MBcoll, MMcoll;
  std::vector<int> BBcolly, MBcolly, MMcolly;
  std::vector<int> BBy,BBf,MBy,MBf;
  std::vector<int> BBStry,BBStrf,MBStry,MBStrf,DecStry,DecStrf;
public:
  CollisionHistory(double ecm,double ftime,double dta,double gam);
  void fill(InterList* inter,const std::vector<EventParticle*>& outgoing);
  void print(int nev,std::string outFile);
};

class AnaTimeDepParticle
{
private:
  double dt; // time slice for time dependent analysis.
  int mprint; // number of time step
  int nPrint;
  int nP; // number of particles
  double gamCM,yCut;
  std::vector<std::vector<double> > npart;
  std::vector<std::vector<int> > mult;
public:
  AnaTimeDepParticle(double ftime, double dta,double gam, double ycut);
  void fill(double ctime,std::list<EventParticle*>& plist);
  void ana(int itime, double ct,std::list<EventParticle*>& plist);
  void print(std::string outfile,const int nevent);
  void init() {nPrint=0;}
};

class AnaTimeDepFlow
{
private:
  double dt; // time slice for time dependent analysis.
  int mprint; // number of time step
  int nPrint;
  int nP; // number of particles
  double gamCM,yCut,yCutF,yCutMax;
  std::vector<std::vector<double> > v1,v2, v1y,v2y, v1r,v2r,v1neg,v1pos;
  std::vector<std::vector<int> > mult,multy,multr;
  std::vector<double> xave,xavey;
  std::vector<int> xmult,xmulty;

  // flow at foward and backward direction.
  std::vector<std::vector<double> > v1f,v2f;
  std::vector<std::vector<int> > multf;

  std::vector<double > v1Lambda,v1LambdaY,v1LambdaF;
  std::vector<int> multLambda,multLambdaY,multLambdaF;

public:
  AnaTimeDepFlow(double ftime, double dta,double gam, double ycut,double yf=1.0,double ymax=3.0);
  void fill(double ctime,std::list<EventParticle*>& plist);
  void ana(int itime, double ct,std::list<EventParticle*>& plist);
  void print(std::string outfile);
  void init() {nPrint=0;}
};

class AnaTimeDepDensity
{
private:
  double dt; // time slice for time dependent analysis.
  int mprint; // number of time step
  int nPrint;
  int nP; // number of particles
  double wG,facG,gamCM,yCut;
  std::vector<double> rho1,rho2,rhos1,rhos2;
  std::vector<std::vector<Pythia8::Vec4> > JB, JB2, JB3;
public:
  AnaTimeDepDensity(double ftime, double dta,double g, double ycut);
  void fill(double ctime,std::list<EventParticle*>& plist);
  void ana(int itime, double ct,std::list<EventParticle*>& plist);
  void print(std::string outfile,int nevent);
  void init() {nPrint=0;}
};

class AnaOutPutPhaseSpace
{
private:
  std::string dir;
  double dt; // time slice for time dependent analysis.
  int mprint; // number of time step
  int nPrint;
  int nP; // number of particles
  double wG,facG;
public:
  AnaOutPutPhaseSpace(double ftime, double dta);
  void fill(int iev,double ctime,std::list<EventParticle*>& plist);
  void ana(int iev,int itime, double ct,std::list<EventParticle*>& plist);
  void init() {nPrint=0;}
};

class AnaOutPutDensity
{
private:
  std::string dir;
  double dt; // time slice for time dependent analysis.
  int mprint; // number of time step
  int nPrint;
  int nP; // number of particles
  double wG,facG;
  std::vector<std::vector<double> > j0,jx,jy,jz;
public:
  AnaOutPutDensity(double ftime, double dta);
  void ana(int iev,int itime, double ct,std::list<EventParticle*>& plist);
  void init() {nPrint=0;}
};

class AnaMeanField
{
private:
  double dt; // time slice for time dependent analysis.
  int mprint; // number of time step
  int nPrint;
  int nP; // number of particles
  double gamCM,yCut;
  std::vector<double>  rhoS,rhoB,potSk,potSm;
  std::vector<Pythia8::Vec4> potVk,potVm;
public:
  AnaMeanField(double ftime, double dta,double gam, double ycut);
  void fill(double coltime,std::list<EventParticle*>& plist);
  void print(std::string outfile,int nevent);
  void init() {nPrint=0;}
};

class ParticleDensity
{
private:
  Pythia8::Settings *settings;
  //Fluid* fluid;
  double rho,rhob,einv;
  int optGauss;
  double widG,widG2,widCof,gVolume;
  double dX,dY,dZ;
  int optPropagate,optDens,optPreHadron;
public:
  //ParticleDensity(Pythia8::Settings* s, Fluid* f);
  ParticleDensity(Pythia8::Settings* s);
  double getRho() {return rho;}
  double getRhob() {return rhob;}
  double getEdens() {return einv;}
  int  computeEnergyMomentum(std::list<EventParticle*>& plist,EventParticle* i1,
	double ctime,int iopt);
  double gaussSmear(const Vec4& pv, const Vec4& dr);

};

}
#endif // jam2_JamAna_H


