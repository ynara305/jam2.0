#include <jam2/collision/EventParticle.h>
#include <jam2/collision/TimeOrder.h>
#include <jam2/collision/InterList.h>

namespace jam2 {

using namespace std;

EventParticle::EventParticle(const int id0, ParticleDataEntry* pd)
{
  id = id0;
  pdePtr=pd;
  pid=0;
  status = 1;
  parent = 0;
  r=0.0; p=0.0; v=0.0;
  mass = 0;
  tform = 0;
  lifeTime = 1e+35;
  tlastcl=0.0;
  numberOfCollision=0;
  lastcol=-1;
  baryonType=0;
  chargeType=0;
  strangeType=0;
  constQ[0]=0;
  constQ[1]=0;
  meanField=false;
  potentialId=0;
  for(int i=0;i<10;i++) facPot[i]=1.0;
  MfromCM.reset();
  if(id !=0 && pdePtr !=0) setProperty();
  rhoS=0.0;
  rhoB=0.0;
  potS=0.0;
  potSm=0.0;
  potV=0.0;
  potVm=0.0;
  force=0.0;
  cBox=nullptr;
  twall=1e+35;
}

EventParticle::EventParticle(int id0, double m,Vec4 r0, Vec4 p0,
	    ParticleDataEntry* pd, int ip)
	: id(id0), mass(m), r(r0), p(p0) ,pdePtr(pd),numberOfCollision(ip)
{
  pid=0;
  status=1; parent=0;
  v=r;
  tform=r[0];
  lifeTime=1e+35;
  tlastcl=v[0];
  lastcol=-1;
  MfromCM.reset();
  meanField=false;
  potentialId=0;
  for(int i=0;i<10;i++) facPot[i]=1.0;
  MfromCM.reset();
  if(id !=0 && pdePtr !=0) setProperty();
  rhoS=0.0;
  rhoB=0.0;
  potS=0.0;
  potSm=0.0;
  potV=0.0;
  potVm=0.0;
  force=0.0;
  cBox=nullptr;
  twall=1e+35;
}

void EventParticle::removeBox() 
{
  if(cBox==nullptr) {
    cout << "EventParticle::removeBox cbox? "<< cBox << endl;
    exit(1);
  }
  cBox->removeParticle(this);
  cBox=nullptr;
}

void EventParticle::changeBox(CascBox* box) 
{
  if(box==cBox) return;
  removeBox();
  cBox=box;
  cBox->add(this);
}

// check if this particle is converted into fluid
bool EventParticle::isconv(const int opt)
{
  if( mass  < 1e-5) return false;     // photons

  //...Only mesons are converted into fluid
  if(opt ==  1) {
    if(baryonType != 0) return false;;

  //....Exclude leading baryons.
  } else if(opt == 2) {
    if(tform > r[0] && baryonType == 3) return false; ;

  //...Exclude leading hadrons.
  } else if(opt ==  3) {
    if( tform > r[0] ) return false;
  }

  return true;
}

bool EventParticle::isMeanField(const double t, const int opt, const double* opt2)
{
  if(!isMeanField(t,opt)) return false;
  if(strangeType !=0) {
    for(int i=0;i<10;i++) facPot[i]=opt2[i];
  }
  return true;
}

bool EventParticle::isMeanField(const double t, const int opt)
{
  meanField=false;
  switch (opt) {
    case 0: return false;

    // only formed baryons feel potentials.
    case 1: 
      if(baryonType == 0) return false;
      //if(r[0] < tform) return false;
      if(t < tform) return false;
      if(r[0] > t) return false;
      meanField=true;
      return true;
      break;

    // pre-formed baryons can feel potentials.
    case 2: 
      if(baryonType == 0) return false;
      if(r[0] > t) return false; // does not have const.quark
      meanField=true;
      return true;
      break;

    // constituent quarks can interact by potentials.
    case 3:
      if(r[0] > t) return false;   // not formed
      // exclude formed mesons.
      if(t > tform && baryonType == 0)  return false;
      meanField=true;
      return true;
      break;
    default:
      if(baryonType==0) return false;
      if(t < tform ) return false;
      meanField=true;
      return true;
  }
      
}

void EventParticle::print(ostream& os) const
{
  os << setw(8) << pdePtr->name(id)
     << setw(9) <<  id
     << setw(12) <<  " col= " << numberOfCollision
     << setw(6) <<  " q1= " << constQ[0]
     << setw(6) <<  " q2= " << constQ[1]
     << setw(6) <<  " t= " << r.e()
     << setw(6) <<  " z= " << r.pz()
     << setw(6) <<  " m= " << mass
     << setw(6) <<  " y= " << p.rap()
     << setw(6) <<  " r= " << r
     << setw(6) <<  " p= " << p
     <<endl;
}

} // end of namespace jam2
