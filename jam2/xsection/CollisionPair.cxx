#include <jam2/xsection/CollisionPair.h>
#include <jam2/collision/EventParticle.h>

namespace jam2 {
using namespace std;

CollisionPair::CollisionPair(int icoltype,
	    EventParticle* p1, EventParticle* p2, double srt,double pr) 
{
  collType=icoltype;
  pid[0]  = p1->getPID();
  pid[1]  = p2->getPID();
  id[0]   = p1->getID();
  id[1]   = p2->getID();
  iz[0]   = p1->charge();
  iz[1]   = p2->charge();
  ibar[0] = p1->baryon();
  ibar[1] = p2->baryon();
  str[0]  = p1->strange();
  str[1]  = p2->strange();
  m[0]    = p1->getMass();
  m[1]    = p2->getMass();

  //eCM=srt;

  eCM = sqrt(m[0]*m[0]+pr*pr) + sqrt(m[1]*m[1]+pr*pr);

  //eCM = srt - p1->getEffectiveMass() + m[0] - p2->getEffectiveMass() + m[1];

  //double m1=p1->getParticleDataEntry()->m0();
  //double m2=p2->getParticleDataEntry()->m0();
  //eCM0=sqrt(m[0]*m[0]+pr*pr)+sqrt(m[1]*m[1]+pr*pr);

  pCM=pr;
  sig=0.0;
  sigel=0.0;
  sigabs=0.0;
  sigt=0.0;
  sigelbw=0.0;
  sigelbw=0.0;
  sigch=0.0;
  mMode=1;
  iChannel=0;
  outgoing.clear();
  anti=false;
  hasAnti[0] = p1->getParticleDataEntry()->hasAnti();
  hasAnti[1] = p2->getParticleDataEntry()->hasAnti();

  //qFac[0]=p1->qFactor();
  //qFac[1]=p2->qFactor();

  //pdata[0] = p1->getParticleDataEntry();
  //pdata[1] = p2->getParticleDataEntry();
}

int CollisionPair::selectChannel() 
{
  //checkChannel();

  if(xSig > sig-sigel) {
    cout << "(CollisionPari::selectChannel) xsig ? "<< xSig << " sigin= "<< sig-sigel<<endl;
    //exit(1);
  }

  double sigint=0.0;
  double xSigSave = xSig;
  for(int i=0, N=outgoing.size(); i<N; i++) {
    sigint += outgoing[i].sig;
    xSig -= outgoing[i].sig;
    if(xSig <= 0.0) {
      iChannel=i;
      return iChannel;
    }
  }
    std::cout << "CollisionPair::selectChannel: no channel? sum sigin= "<< sigint
	      << " xSigSave= "<< xSigSave
	      << " xSig= "<< xSig
                 << " ecm= " << eCM
		 << " id1= " << id[0]
		 << " id2= " << id[1]
		 << endl;
    std::cout << " outgoing size= " << outgoing.size();
    std::cout << " has anti= "<< isAnti() <<endl;
    std::cout<< " sig= " << sig << " sigel= " << sigel
	<< " sigintot= " << sigint<< endl;
    std::cout<< " sig-sigel-sigint= " << sig-sigel-sigint<< endl;

    for(int i=0, N=outgoing.size(); i<N; i++) {
	std::cout << " out1= " << outgoing[i].id1
	     << " out2= " << outgoing[i].id2
	     << " sig= " << outgoing[i].sig
	     <<endl;
    }
    exit(1);

     return -1;
}

void CollisionPair::checkChannel() 
{
    std::cout << "checkChannel:  "
                 << " ecm= " << eCM
		 << " id1= " << id[0]
		 << " id2= " << id[1]
		 << endl;
    double sigint=0.0;
    for(int i=0, N=outgoing.size(); i<N; i++) {
	sigint += outgoing[i].sig;
	std::cout << i << " out1= " << outgoing[i].id1
	     << " out2= " << outgoing[i].id2
	     << " sig= " << outgoing[i].sig
	     <<endl;
    }

    std::cout << " outgoing size= " << outgoing.size()<<endl;
    std::cout<< " sig= " << sig << " sigel= " << sigel
	<< " sigintot= " << sigint<< endl;
    std::cout<< " sig-sigel-sigint= " << sig-sigel-sigint<< endl;

}

} // end namespace jam2
