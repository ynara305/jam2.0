#include <jam2/collision/TwoBodyInterList.h>
#include <jam2/collision/EventParticle.h>

namespace jam2 {

using namespace std;

void TwoBodyInterList::print(ostream& os) const
{
    cout << endl;

    double px1 = scatt[0]->getPx();
    double py1 = scatt[0]->getPy();
    double pz1 = scatt[0]->getPz();
    double pe1 = scatt[0]->getPe();

    cout << "*** Twobody inter: Before Collision ***"<<endl;
    cout << "print: inter= " << this << endl;
    cout << "scat0= " << scatt[0] << " scat1=  " << scatt[1] <<endl;

    double px2 = scatt[1]->getPx();
    double py2 = scatt[1]->getPy();
    double pz2 = scatt[1]->getPz();
    double pe2 = scatt[1]->getPe();

    double s = (pe1+pe2)*(pe1+pe2)
	       -(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2);

    double srt = sqrt(max(0.0,s));

    os << "# coll time= " <<  collisionOrderTime;
    os << " ct1=" <<  collTime[0];
    os << " ct2=" <<  collTime[1];
    os << " ecm= " << srt;
    os << " b= " << impactPar;
    os << " sigma= " << cpair.getSigma();
    os << endl;
    (*scatt[0]).print(os);
    (*scatt[1]).print(os);
    os << " box1= "<< scatt[0]->box() << " box2= "<< scatt[1]->box() <<endl;
    os << " pttox= "<< px1+px2 << " ptotoy= "<< py1+py2 << " ptotz= "<< pz1+pz2<< " etot= "<< pe1+pe2 << endl;
    os << endl;
    /*
    if(collTime[0]< scatt[0]->getT()  || collTime[1]< scatt[1]->getT() ) {
      os << " collision time strange !"<< endl;
      exit(1);
    }
    */


}

} // end of namespace jam2


