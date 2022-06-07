#include <jam2/collision/CascString.h>

namespace jam2 {

using namespace std;

void CascString::printString(ostream& os) const
{
    os << "mass of string = " << getMass()
       << " # coll. = " <<  this->getNColl()
       << endl;
    list<CascParton*>::const_iterator ci;
    int i=0;
    for (ci = sparts.begin(); ci != sparts.end(); ci++) {
	os << i << " kf=" << (*ci)->getID() ;

        if((*ci)->getColorFrom() !=0)
	    os <<  " from " << (*ci)->getColorFrom()
	       	<< " kf=" << (*ci)->getColorFrom()->getID();

        if((*ci)->getColorTo() !=0)
	    os << " to " << (*ci)->getColorTo()
	       << " kf= " << (*ci)->getColorTo()->getID();

	os << " e= " << (*ci)->getPe();
	os << " ci= " << *ci;
        os << endl;

	/*
	if((*ci)->stringTo != this) {
	    os << " partons does not indicate this string ";
	    os << " stringTo = " << (*ci)->stringTo << endl;
	}
	*/
	i++;
    }

}

}
