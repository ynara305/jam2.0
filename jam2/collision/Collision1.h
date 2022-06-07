#ifndef jam2_collision_Collision1_h
#define jam2_collision_Collision1_h

#include <list>
#include <vector>
#include <set>
#include <fstream>
#include <jam2/collision/EventParticle.h>
#include <jam2/collision/TimeOrder.h>
#include <jam2/collision/InterList.h>
#include <jam2/xsection/CrossSection.h>
#include <Pythia8/Settings.h>
#include <jam2/collision/Collision.h>

namespace jam2 {

typedef std::list<EventParticle*>::iterator EventParIt;

class Collision1 : public Collision
{
private:
    typedef std::multiset<InterList*,TimeOrder> InterListSet;
    typedef InterListSet::const_iterator InterCIt;
    typedef InterListSet::iterator InterIt;
    InterListSet interList;

public:
    Collision1(Pythia8::Settings* s, JamParticleData* jp, CrossSection* in,Pythia8::Rndm* r)
	: Collision(s,jp,in,r) {}
    ~Collision1();
    int  getInterListSize() const {return interList.size();}
    void makeCollisionList(double itime,double gtime);
    int computeNumberOfParticipants();
    InterList* findNextCollision();
    void collisionUpdate(InterList* inter);
    void collisionUpdate(vector<EventParticle*>& outgoing,double itime,double gtime);
    void cancelCollision(InterList* inter);
    double finalCollisionTime() {
      if(interList.size() ==0) return 0;
      //InterList* lastColl = *interList.rbegin();
      //cout << " lastcoll= " << lastColl<<endl;
      //double t = lastColl->getCollisionOrderTime();
      //cout << " t=  " << t<<endl;
      //return t;
      return (*interList.rbegin())->getCollisionOrderTime();
    }
  bool checkNextCollisionTime(TwoBodyInterList* inter,double dtau1,double dtau2,
                 bool preHadronA,bool preHadronB);
  void clear();
  void removeInteraction(EventParticle* i1);

private:
  void removeInteraction2(InterList* inter);

};
}
#endif
