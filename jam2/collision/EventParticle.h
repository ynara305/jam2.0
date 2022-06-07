#ifndef jam2_collision_EventParticle_H
#define jam2_collision_EventParticle_H

#include <fstream>
#include <set>
#include <cmath>
#include <Pythia8/Basics.h>
#include <Pythia8/ParticleData.h>
#include <jam2/collision/CascBox.h>
#include <jam2/xsection/CollisionPair.h>

namespace jam2 {

using namespace Pythia8;

class EventParticle
{
private:
  int     status;     // status
  int     id;         // PDG flavour code
  int     pid;        // particle ID defined in JamStdlib.h
  int     parent;
  double  mass;       // rest mass in GeV
  Pythia8::Vec4 r;    // position in fm
  Pythia8::Vec4 p;    // momentum in GeV
  Pythia8::Vec4 v;    // vertex in fm
  Pythia8::Vec4 force;
  Pythia8::ParticleDataEntry* pdePtr;
  CascBox *cBox;
  double  tform;     // formation time in fm
  double  lifeTime;  // life time in fm
  double  tlastcl;   // time of last collision in fm
  double  twall;     // wall collision time in fm.
  int     numberOfCollision;
  int     lastcol;
  int     baryonType;   // 3 times the baryon number
  int     chargeType;
  int     strangeType;
  int     constQ[2];
  bool    meanField;
  double  facPot[10];   // factor to be multiplied to the Skyrme potential
  int     potentialId;
  double  rhoS, rhoB;   // scalar and  baryon density
  double  potS,potSm;        // scalar potential
  Pythia8::Vec4 potV,potVm;  // vector potential
  Pythia8::RotBstMatrix MfromCM;

public:
  EventParticle(const int id0=0, ParticleDataEntry* pd=0);
  EventParticle(int id0, double m,Vec4 r0, Vec4 p0,ParticleDataEntry* pd=0, int ip=0);
  virtual ~EventParticle() { cBox=nullptr;}

  void setProperty() {
    baryonType=pdePtr->baryonNumberType(id);
    chargeType=pdePtr->chargeType(id);
    strangeType=pdePtr->nQuarksInCode(3);
    if(baryonType==0) {
      strangeType = id > 0 ? strangeType : -strangeType;
      // eta' phi f0(1710)...
      if((abs(id)/100)%10 ==3 && (abs(id)/10)%10==3) strangeType=0;
      // Bs
      else if((abs(id)/100)%10 ==5 && (abs(id)/10)%10==3) {
        strangeType = id > 0 ? -1 : 1;
      }
    } else {
      strangeType = id < 0 ? strangeType : -strangeType;
    }
    constQ[0] = 1;
    constQ[1] = abs(baryonType) == 3 ? 2: 1; 
    meanField=false;
  }

  Pythia8::ParticleDataEntry* getParticleDataEntry() {return pdePtr;}
  void setParticleData(Pythia8::ParticleDataEntry* pd) {pdePtr=pd;}

  void setBoostMatrix(Pythia8::RotBstMatrix m) {MfromCM=m;}
  Pythia8::RotBstMatrix getMfromCM() {return MfromCM;}
  Pythia8::RotBstMatrix getMtoCM() {
    Pythia8::RotBstMatrix m=MfromCM; m.invert();
    return m;
  }

  int   getID()          const {return id;}
  int   getPID()         const {return pid;}
  int   getStatus()          const {return status;}
  int   getParent()      const {return parent;}
  int   getNColl()  const {return numberOfCollision;}
  int   lastColl() const {return lastcol;}
  void  lastColl(int i) {lastcol=i;}
  int   baryon() const {return baryonType;}
  int   charge() const {return chargeType;}
  int   strange() const {return strangeType;}
  int*  constQuark()  {return constQ;}
  int   constQuark(int i)  {return constQ[i];}
  void  setConstQuark() {
    constQ[0] = 1;
    constQ[1] = abs(baryonType) == 3 ? 2: 1; 
  }
  void  setConstQuark(int q[2]) {constQ[0]=q[0];constQ[1]=q[1];}
  void  setConstQuark(int i,int q) {constQ[i]=q;}
  double qFactor() const {
    int nq = baryonType == 0 ? 2 : 3;
    return (double)(constQ[0]+constQ[1])/nq;
  }

  CascBox* box() const {return cBox;}
  void clearBox() {cBox=nullptr;}
  void removeBox();
  void addBox(CascBox* b) {cBox=b; cBox->add(this);}
  void changeBox(CascBox* box);

  Vec4 velocity(int opt) {
    if(opt !=0) return p/p[0];
      Vec4 pk = getPkin();
      return pk/pk[0];
  }

  Vec4 propagate(double t,int opt) {
        Vec4 rp;
	Vec4 vel = p/p[0];
	// make v* = p*/e*
	if(opt==1) {
	  vel = p - potV;
	  double ee=sqrt(pow2(mass+potS)+vel.pAbs2());
          vel /= ee; 
	}
	rp = r + vel*(t - r[0]);
	rp[0] = t;
	return rp;
  }
  void updateR(double t, int opt=0) {
	Vec4 vel = p/p[0];
	// make v* = p*/e*
	if(opt==1) {
	  vel = p - potV;
	  double ee=sqrt(pow2(mass+potS)+vel.pAbs2());
          vel /= ee; 
	}
	r += vel*(t - r[0]);
	r[0] = t;
	checkConstQuark(t);
  }
  void updateByForce(const Vec4& forcep, const Vec4& fr,const double dt) {
    force=fr;
    double meff = mass + potS;
    p += forcep*dt;
    p[0] = sqrt(meff*meff + p.pAbs2());
    // displacement only with interaction term.
    r[1] += force[1]*dt;
    r[2] += force[2]*dt;
    r[3] += force[3]*dt;
  }
  void checkConstQuark(double t) {
    if(t >= tform)  {
      constQ[0]=1;
      constQ[1] = abs(baryonType) == 3 ? 2: 1; 
    }
  }

  double tWall() const {return twall;}
  void setTWall(double t) {twall = t;}

  //double collisionOrderTime() const;

  bool hasConst() {
    if( r[0] < tform ) return true;
    else return false;
  }

  void    baryon(int i)    {baryonType=i;}
  void    setID(int i)     {id = i;}
  void    setPID(int i)    {pid = i;}
  void    setStatus(int i)     {status = i;}
  void    setParent(int i) {parent = i;}

  double  mCalc()  {return p.mCalc();}
  double  getMass() const { return mass;}
  double  getEffectiveMass() const { return mass+potS;}
  double  lifetime() const {return lifeTime;}
  void    setLifeTime(double t) {lifeTime=t;}

  Vec4 forceR() const {return force;}
  double rhob() const {return rhoB;}
  double rhos() const {return rhoS;}
  void setRhoS(double rin) {rhoS=rin;}
  void setRhoB(double rin) {rhoB=rin;}
  double pots() {return potS;}
  double potsm() {return potSm;}
  Vec4 potv() {return potV;}
  Vec4 potvm() {return potVm;}
  double potv(int i) {return potV[i];}
  void setPotV(Vec4 vin) {potV=vin;}
  void setPotVm(Vec4 vin) {potVm=vin;}
  void setPotV(int i,double vin) {potV[i]=vin;}
  void setPotVm(int i,double vin) {potVm[i]=vin;}
  void setPotSm(double s) {potSm = s;}
  void setPotS(double s, bool opt=true) {
    potS=s; 
    if(!opt) return;
    double meff= mass + potS;
    if(meff<0.0) {
      cout << "EventParticle::setPotS meff < 0? "<< meff << " m= "<< mass
   	   << " potS= " << potS <<endl;
      exit(1);
    }
    p[0]=sqrt(meff*meff + p.pAbs2());
  }
  void setOnShell() {
    p[0]=sqrt((mass+potS)*(mass+potS)+p.pAbs2());
  }
  void setFree() {
    p += potV;
    p[0]=sqrt(mass*mass+p.pAbs2());
  }
  void setFreeMass() {
    p[0]=sqrt(mass*mass+p.pAbs2());
  }
  void setKinetic() {
    p -= potV;
    p[0]=sqrt((mass+potS)*(mass+potS)+p.pAbs2());
  }
  Vec4 getPkin() {
    Vec4 pk = p - potV;
    pk[0]=sqrt((mass+potS)*(mass+potS)+pk.pAbs2());
    return pk;
  }
  Vec4 getPcan(int optvdot=2) {
    Vec4 pk = optvdot > 1 ? p + potV : p ;
    pk[0]=sqrt(mass*mass+pk.pAbs2());
    return pk;
  }

  double  getPx() const { return p.px();}
  double  getPy() const { return p.py();}
  double  getPz() const { return p.pz();}
  double  getPe() const { return p.e();}
  double  getE0() const {return sqrt(mass*mass + p.pAbs2());}
  double  pAbs2() {return p.pAbs2();}
  double  getX()  const { return r.px();}
  double  getY()  const { return r.py();}
  double  getZ()  const { return r.pz();}
  double  getT()  const { return r.e();}
  double  getTf() const { return tform;}
  double  TimeLastColl() const {return tlastcl;}
  Vec4 getR() {return r;}
  Vec4 getP()  {return p;}
  Vec4 getP() const {return p;}
  Vec4 getV() {return v;}
  Vec4 getMomentum() const {return p;}
  double getR(int i)  { return r[i];}
  double getP(int i)   { return p[i];}
  double getV(int i)  { return v[i];}

  void    setNumberOfColl(int n) {numberOfCollision=n;}
  void    setFormationTime(double t)  {tform=t;}
  void    setT(double t)  {r.e(t);}
  void    setX(double x)  {r.px(x);}
  void    setY(double y)  {r.py(y);}
  void    setZ(double z)  {r.pz(z);}
  void    setCoordinate(Vec4  rr) {r=rr;}
  void    setVertex(Vec4&  rr) {v=rr;tlastcl=v[0];}
  void    setTimeLastColl() {tlastcl=r[0];}
  void    setVertex() {v=r;}
  //void    setMomentum(double p0[4]) { for (int i=0; i<4; i++) p[i] = p0[i]; }
  void    setMomentum(Vec4  pp) {p=pp;}
  void    setR(Vec4  rr) {r=rr;}
  void    setP(Vec4  pp) {p=pp;}
  void    addP(Vec4  pp) {p +=pp;}
  void    multP(double  a) {p[1] *=a; p[2] *=a; p[3] *=a;}
  void    bst(Vec4& pin) {p.bst(pin);}
  void    bst(Vec4& pin, double m) {p.bst(pin,m);}
  void    bstback(Vec4& pin) {p.bstback(pin);}
  void    bstback(Vec4& pin,double m) {p.bstback(pin,m);}

  void    addT(double t)  {r[0] += t;}
  void    addX(double x)  {r[1] += x;}
  void    addY(double y)  {r[2] += y;}
  void    addZ(double z)  {r[3] += z;}

  void    setE(double e)    {p.e(e);}
  void    setPx(double px)  {p.px(px);}
  void    setPy(double py)  {p.py(py);}
  void    setPz(double pz)  {p.pz(pz);}
  void    setMass(double m)  {mass=m;}
  void print(ostream& os = std::cout) const;
  bool isconv(const int opt);
  bool isMeanField(const double t, const int opt, const double* opt2);
  bool isMeanField(const double t, const int opt);
  bool meanFieldOn() {return meanField;}
  double facPotential(const int i) {return facPot[i];}
  double* facPotential() {return facPot;}
  int     potentialID() {return potentialId;}
  void setPotentialParam(double* p0) {for(int i=0;i<10;i++) facPot[i]=p0[i];}
  void setPotentialId(int i) {potentialId=i;}

};

}

#endif
