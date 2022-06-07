#ifndef jam2_collision_cascade_CascParton_h
#define jam2_collision_cascade_CascParton_h

#include <jam2/collision/EventParticle.h>

namespace jam2 {

class CascString;

class CascParton : public EventParticle
{
private:
    // This is to maintain net color flow. 
    CascParton* colorFrom;
    CascParton* colorTo;
    CascString* stringTo;    // pointer to stringList.
    int         constQuark;  // Flag to tell constituent quark or not.

public:
    CascParton(int kf): EventParticle(kf) { 
	colorFrom  = 0;
	colorTo    = 0;
	stringTo   = 0;
	constQuark = 0;
    }
    CascParton(const CascParton& arg): EventParticle(arg) {
	colorFrom  = arg.colorFrom;
	colorTo    = arg.colorTo;
	stringTo   = arg.stringTo;
	constQuark = arg.constQuark;
    }
    CascParton*  getColorFrom() const {return colorFrom;}
    CascParton*  getColorTo()   const {return colorTo;}
    CascString*  getString()    const {return stringTo;}
    int          getConstQuark() const {return constQuark;}
    void setColorFrom(CascParton* cs) {colorFrom = cs;}
    void setColorTo(CascParton* cs)   {colorTo = cs;}
    void setString(CascString* cs)    {stringTo = cs;}
    void setConstQuark(int i)         {constQuark = i;}

};
}
#endif /* jam2_collision_cascade_CascParton_h */
