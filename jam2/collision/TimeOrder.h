#ifndef jam2_collision_TimeOrder_h
#define jam2_collision_TimeOrder_h

#include <jam2/collision/InterList.h>
//#include <jam2/collision/InterType.h>

class InterList;

namespace jam2 {
class TimeOrder {
public:
   bool operator() (const InterList* i1, const InterList* i2) const {
     return i1->getCollisionOrderTime() < i2->getCollisionOrderTime();
   }
};

class TimeOrderDescend {
public:
  bool operator() (const InterList* i1, const InterList* i2) const {
    return i1->getCollisionOrderTime() > i2->getCollisionOrderTime();
  }
};

/*
class CollisionTimeOrder {
    public:
	bool operator() (const InterType& i1, const InterType& i2) const {
	return i1.collisionTime() > i2.collisionTime();
	}
};

class PCollisionTimeOrder {
    public:
	bool operator() (const EventParticle* i1, const EventParticle* i2) const {
	return i1->collisionOrderTime() > i2->collisionOrderTime();
	}
};
*/

}
#endif
