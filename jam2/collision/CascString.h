#ifndef jam2_collision_cascade_CascPart_h
#define jam2_collision_cascade_CascPart_h

#include <jam2/collision/EventParticle.h>
#include <jam2/collision/CascParton.h>
#include <list>
#include <iostream>

namespace jam2 {

class CascString : public EventParticle
{
private:
    typedef std::list<CascParton*>  strList;
    strList  sparts;  // List of the particles in the string

public:
    CascString(int kf): EventParticle(kf) { }
    virtual ~CascString() {sparts.clear();}
    void printString(std::ostream& os) const;

    strList& getPartons() {return sparts;}
    void setParton(CascParton* cp) {sparts.push_back(cp);}
    int  getStringSize() {return sparts.size();}

    strList::const_iterator getStrBegin() const { return sparts.begin();}
    strList::const_iterator getStrEnd()   const { return sparts.end();}
    void replaceParton(CascParton* cold, CascParton* cnew) {
			replace(sparts.begin(),sparts.end(),cold,cnew); }

};

}
#endif /* jam2_collision_cascade_CascPart_h */
