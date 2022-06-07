#include <jam2/collision/WallList.h>

namespace jam2 {

using namespace std;
void WallList::print(ostream& os) const
{
  cout << "*** Before Wall collision ***"<<endl;
  cout << "print: inter= " << this << " p= "<< scatt[0] << endl;

  os << "# wall colllision time= " <<  collisionOrderTime;
  os << endl;

  scatt[0]->print(os);
  os << endl;

}

} // end of namespace jam2
