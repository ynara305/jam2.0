#ifndef jam2_collision_WallList_H
#define jam2_collision_WallList_H

#include <ostream>
#include <list>
#include <jam2/collision/InterList.h>

namespace jam2 {

class WallList : public InterList
{
public:
  WallList() { }
  WallList(EventParticle* i1, double t) { 
    collisionOrderTime=t;
    collType=-1;
    scatt[0] = i1;
    scatt[1] = nullptr;
    collTime[0]=collisionOrderTime;
    collTime[1]=collisionOrderTime;
  }

  int getNumberOfInComing() const {return 0;}
  Pythia8::Vec4 getTotalMom() {return scatt[0]->getP();}
  void print(std::ostream& os=std::cout) const;

};

}
#endif
