#ifndef jam2_collision_Collision_h
#define jam2_collision_Collision_h

#include <list>
#include <vector>
#include <set>
#include <fstream>
#include <jam2/collision/EventParticle.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/collision/TimeOrder.h>
#include <jam2/collision/TwoBodyInterList.h>
#include <jam2/xsection/CrossSection.h>
#include <jam2/initcond/InitialCondition.h>
#include <Pythia8/Settings.h>
#include <jam2/collision/GGCF.h>

namespace jam2 {

class Collision
{
protected:
  Pythia8::Settings* settings;
  JamParticleData* jamParticleData;
  CrossSection* xsection;
  Pythia8::Rndm* rndm;
  int cascadeMethod;
  double  passingTime,lastCollisionTime;
  double  initialTime,finalTime;
  int     overSample;
  int     withBox;
  double  xBox,yBox,zBox;
  int     numberOfCollision, numberOfParticipants;
  int     predictedColl;
  int     numberOfDecay, numberOfWallCollision;
  double  minKinEnergy;
  int     optCollisionOrdering, optCollisionTimeLimit;
  int     optVectorPotential,optVdot,optPropagate;
  int     optFluctuation;
  double  omegaSigma; // parameter for color fluctuation.
  bool decayOn, inelOnly, bbCollisionOnly,constQCollisionOnly,noMMCollision;
  double rMaxCutSq,rMaxCutSqBB,rMaxCutSqMB,rMaxCutSqMM,rMaxCutSqBBar; 
  bool removeSpectator;
  double rMaxCutSqNN;
  double bMax;
  bool optTau;
  GGCF *ggcf;

  std::list<EventParticle*> pnew;
  std::list<EventParticle*> pnewA, pnewB;
  double gWidth, pauliR, pauliP, pauliC;
  int hitPath;

public:

  std::list<EventParticle*> plist;
    virtual void init(const InitialCondition* inic) {}
    virtual void clear()=0;
    void clearPlist();
    virtual InterList* findNextCollision()=0;
    virtual void collisionUpdate(InterList* inter)=0;
    virtual void collisionUpdate(vector<EventParticle*>& outgoing,double it,double t)=0;
    virtual void removeInteraction(EventParticle* i1)=0;
    //virtual int  getInterListSize() const=0;
    virtual void makeCollisionList(double itime,double gtime)=0;
    virtual void cancelCollision(InterList* inter)=0;
    virtual bool checkNextCollisionTime(TwoBodyInterList* inter,double dtau1,double dtau2,
                 bool preHadronA,bool preHadronB)=0;
    virtual void wallCollision(InterList& inter) { }
    virtual void propagate(double ctime,int opt, int n=2);
    virtual void check(std::string name) {  }

  double nuclearPassingTime() const {return passingTime;}
  void makeCollisionListTwoNuclei(double itime,double gtime);
  int computeNumberOfParticipants(std::multiset<InterList*,TimeOrder>& interlist);
  int getNumberOfParticipants() {return numberOfParticipants;}

    typedef std::list<EventParticle*> listEV;
    typedef listEV::iterator          listEVIt;
    typedef listEV::const_iterator    listEVcIt;

    //std::list<EventParticle*> getPlist()   const    {return plist;}
    listEVIt  plistEnd()   { return plist.end();}
    int  getPlistSize() const {return plist.size();}

    std::list<EventParticle*> getPnewList() const   {return pnew;}
    void setPList(EventParticle* ep) {plist.push_front(ep);} // ep->setPlistIt(plist.begin());}
    void setPListB(EventParticle* ep) {plist.push_back(ep);} // ep->setPlistIt(plist.begin());}
    void setPList(std::vector<EventParticle*>& outgoing) {
     for(const auto& p : outgoing) plist.push_back(p);}
    void setPnewList(EventParticle* ep) {pnew.push_back(ep);}
    void setPnewAList(EventParticle* ep) {pnewA.push_back(ep);}
    void setPnewBList(EventParticle* ep) {pnewB.push_back(ep);}
    listEVIt  getPnewListBegin() { return pnew.begin();}
    listEVIt  getPnewListEnd()   { return pnew.end();}
    listEVcIt getPnewListBegin() const { return pnew.begin();}
    listEVcIt getPnewListEnd()   const { return pnew.end();}
    void erasePnewList(listEVIt ev) { pnew.erase(ev);}
    void erasePList(listEVIt ev) { plist.erase(ev);}
    void eraseParticle(EventParticle* ip) {
      delete ip;
      plist.erase(std::find(plist.begin(),plist.end(),ip));
      //plist.erase(i1->plistIt());
    }

    int  getPnewlistSize() const {return pnew.size();}

    static bool dectim(EventParticle* a, EventParticle* b) {
	return ( a->lifetime() < b->lifetime() );
    }

    Collision(Pythia8::Settings* s, JamParticleData* jp, CrossSection* in,Pythia8::Rndm* r);
    virtual ~Collision();

    double getFinalTime()      const { return finalTime;}
    double collisionOrderTime() const { return lastCollisionTime;}
    int getNColl()             const { return numberOfCollision;}
    int getNWallColl()         const { return numberOfWallCollision;}
    int getNDecay()            const { return numberOfDecay;}
    int getNumberOfOperation() const { return numberOfCollision+numberOfDecay;}
  int getNumberOfGlauberCollision() const {return predictedColl;}

    void setFinalTime(double t)       {finalTime = t;}
    void setNumberOfCollision(int n)  {numberOfCollision = n;}
    virtual TwoBodyInterList* hit(EventParticle* i1, EventParticle* i2);
    virtual bool doPauliBlocking(InterList* inter,std::vector<EventParticle*>& out,int opt);
    int collisionType(EventParticle* i1, EventParticle* i2);
    void deleteParticle();
    void deleteParticle(int i);
    void PrintCollision(ostream& os) const;


};
}
#endif
