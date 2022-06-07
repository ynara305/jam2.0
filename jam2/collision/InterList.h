#ifndef jam2_collision_InterList_H
#define jam2_collision_InterList_H
#include <list>
#include <jam2/collision/EventParticle.h>
#include <Pythia8/Basics.h>
#include <jam2/xsection/CollisionPair.h>

namespace jam2 {

class EventParticle;

class InterList
{
protected:
    int      collType;
    double   collisionOrderTime;
    EventParticle* scatt[2];
    double       collTime[2];
public:
    InterList() { }
    virtual ~InterList() { }
    InterList(double ctime) : collisionOrderTime(ctime){
      scatt[0]=nullptr;
      scatt[1]=nullptr;
      collTime[0]=0.0;
      collTime[1]=0.0;
      collType=0;
    }


    //virtual std::list<EventParticle*>::iterator getParticle(int i) const =0;
    virtual int getNumberOfInComing() const=0;
    virtual void print(std::ostream& os=std::cout) const=0;
    virtual Pythia8::Vec4 getTotalMom() =0;
    virtual CollisionPair& getCpair() {std::cout << "InterList::getCpair wrong usage" << std::endl;exit(1);}
    virtual double getImpactPar() const {std::cout << "InterList::getImpactPar wrong usage" << std::endl;exit(1);}

  void   setCollisionTime(int i, double t) {collTime[i] = t;}
  double getCollisionTime(int i) const {return collTime[i];}
  EventParticle* getParticle(int i=0) const {return scatt[i];}
  double getCollisionOrderTime() const {return collisionOrderTime;}
  void   setCollisionOrderTime(double t) {collisionOrderTime = t;}
  int    getCollType() const {return collType;}

};
}
#endif
