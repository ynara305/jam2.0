#ifndef jam2_collision_TwoBodyInterList_H
#define jam2_collision_TwoBodyInterList_H

#include <ostream>
#include <list>
#include <jam2/collision/InterList.h>
#include <jam2/xsection/CollisionPair.h>

namespace jam2 {

class EventParticle;

class TwoBodyInterList : public InterList
{
private:
    CollisionPair cpair;
    //std::list<EventParticle*>::iterator scatt[2];
    //EventParticle* scatt[2];
    //double       collTime[2];
    double       impactPar;
public:
  TwoBodyInterList() { }
  TwoBodyInterList(CollisionPair &cp, 
        EventParticle* i1, EventParticle* i2,
	double ctime, double tcl1=0.0, double tcl2=0.0,double b=0.0) {
	cpair = cp;
	scatt[0]= i1;
	scatt[1]= i2;
	collTime[0] = tcl1;
	collTime[1] = tcl2;
	collisionOrderTime = ctime;
	impactPar = b;
	collType=cpair.getCollType();
  }
  /*
  TwoBodyInterList(const TwoBodyInterList& it) {
	cpair = it.cpair;
	scatt[0]= it.scatt[0];
	scatt[1]= it.scatt[1];
	collTime[0] = it.collTime[0]; 
	collTime[1] = it.collTime[1];
	collisionOrderTime = it.collisionOrderTime;
	impactPar = it.impactPar;
	collType= it.collType;
  }
  */

  int getNumberOfInComing() const {return 2;}
  CollisionPair& getCpair() {return cpair;}
  Pythia8::Vec4 getTotalMom() {return scatt[0]->getP()+scatt[1]->getP();}
  double getImpactPar() const {return impactPar;}
  void   print(std::ostream& os=std::cout) const;

};
}
#endif
