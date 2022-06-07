#ifndef jam2_collision_DecayList_H
#define jam2_collision_DecayList_H

#include <ostream>
#include <list>
#include <jam2/collision/InterList.h>

namespace jam2 {

//class EventParticle;

class DecayList : public InterList
{
private:
    //std::list<EventParticle*>::iterator decayer;
    //EventParticle* decayer;
public:
    DecayList() { }
    DecayList(EventParticle* i1) {
	    scatt[0] = i1;
	    scatt[1] = nullptr;
	    collisionOrderTime = i1->lifetime();
	    collTime[0]=collisionOrderTime;
	    collTime[1]=collisionOrderTime;
	    collType=0;
    }
    int getNumberOfInComing() const {return 1;}
    Pythia8::Vec4 getTotalMom() {return scatt[0]->getP();}
    void   print(std::ostream& os=std::cout) const;

};

}
#endif
