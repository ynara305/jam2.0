#include <jam2/collision/DecayList.h>

namespace jam2 {

using namespace std;
void DecayList::print(ostream& os) const
{
    cout << "*** Before Decay ***"<<endl;
    cout << "print: inter= " << this << endl;
	double srt = scatt[0]->getMass();

    os << "# decay time= " <<  collisionOrderTime;
    os << " m= " << srt;
    os << endl;

    (*scatt[0]).print(os);
    os << endl;

}

} // end of namespace jam2
