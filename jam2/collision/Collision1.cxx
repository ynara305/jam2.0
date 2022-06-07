#include <vector>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <iostream>

#include <jam2/collision/Collision1.h>
#include <jam2/collision/DecayList.h>

namespace jam2 {

using namespace std;

void Collision1::clear()
{
  if(interList.size()>0) {
  for(auto& it:interList) delete it;
  interList.clear();
  }
  clearPlist();
} 

Collision1::~Collision1()
{
    InterIt it;
    for(it = interList.begin(); it != interList.end(); ++it) {
	delete *it;
    }
    interList.clear();
}



// =========================================================================//
//...Purpose: to make collision list into collision predictor arrays.
void Collision1::makeCollisionList(double itime, double gtime)
{
  initialTime=itime;
  finalTime=gtime;
  EventParIt i0 = ++plist.begin();

  // Loop over all particle pairs.
  for(EventParIt i1 = i0; i1 != plist.end(); ++i1)
  for(EventParIt i2 = plist.begin(); i2 != i1; ++i2) {
    TwoBodyInterList* it = hit(*i1,*i2);
    if(it !=nullptr) interList.insert(it);
  }

  predictedColl = interList.size();
  numberOfCollision=0;

  // check decay.
  if(decayOn) {
    for(EventParIt ip = plist.begin(); ip != plist.end(); ++ip) {
      double ct=(*ip)->lifetime();
      if(optTau) {
        Vec4 v = (*ip)->velocity(optPropagate);
        double zc= (*ip)->getZ() +v[3]*(ct-(*ip)->getT());
        if(abs(ct) < abs(zc)) {
          cout << "Collision1:makeCollisionList tc= "<< ct << " zc= "<< zc<<endl;
          exit(1);
        }
        ct = sqrt(ct*ct - zc*zc);
      }
      if(ct < finalTime) interList.insert(new DecayList(*ip));
    }
  }

}


//***********************************************************************
//...Purpose: to search decay and collision arrays to find next event
InterList* Collision1::findNextCollision()
{
    if(interList.size() <= 0 ) return 0;
    InterList* nextColl = *interList.begin();
    lastCollisionTime = nextColl->getCollisionOrderTime();
    return nextColl;
}


//...Remove collisions including i1.
void Collision1::removeInteraction(EventParticle* i1)
{

//=== delete interactions ================================================ 
  InterIt first = interList.begin();
  InterIt last  = interList.end();
  while(first != last) {
    InterIt next = first;
    ++next;
    if((*first)->getParticle(0) == i1 || (*first)->getParticle(1) == i1) {
      delete *first;
      interList.erase(first);
    }
    first = next;
  }

/*
  InterListSet list2;
  for(auto& it : interList) {
    if(it->getParticle(0) == i1 || it->getParticle(1) == i1) {
      delete it;
    } else {
      list2.insert(it);
    }
  }
  interList = list2;
  list2.clear();
  */


  /*
  // vector
  auto it = interList.begin();
  while(it != interList.end()) {
    if((*it)->getParticle(0) == i1 || (*it)->getParticle(1) == i1) {
      delete *it;
      it = interList.erase(it);
    } else ++it;
  }
  */

    /*
    for (InterIt it=interList.begin(); it != interList.end(); ) {
	if((*it)->getParticle(0) == i1 || (*it)->getParticle(1) == i1) {
	    delete *it;
	    it=interList.erase(it);
	    continue;
	}
	it++;
    }
    */

//=== end delete interactions ================================================ 

  //i1->deleteInterList();
  //delete i1;
  //plist.erase(i1->plistIt());
  eraseParticle(i1);
}

//...Remove collisions including i1,i2.
void Collision1::removeInteraction2(InterList* inter)
{
  EventParticle* i1 = inter->getParticle(0);
  EventParticle* i2 = inter->getParticle(1);

//=== delete interactions ================================================ 
  InterIt first = interList.begin();
  while(first != interList.end()) {
	InterIt next = first;
	++next;
	if((*first)->getParticle(0) == i1 || (*first)->getParticle(1) == i1
        || (*first)->getParticle(0) == i2 || (*first)->getParticle(1) == i2) {
	    delete *first;
	    interList.erase(first);
	}
	first = next;
  }

/*
  InterListSet list2;
  for(auto& it : interList) {
    if(it->getParticle(0) == i1 || it->getParticle(1) == i1
    || it->getParticle(0) == i2 || it->getParticle(1) == i2) {
      delete it;
    } else {
      list2.insert(it);
    }
  }
  interList = list2;
  list2.clear();
  */

    /*
    for (InterIt it=interList.begin(); it != interList.end(); ) {
	if((*it)->getParticle(0) == i1 || (*it)->getParticle(1) == i1
   	 ||(*it)->getParticle(0) == i2 || (*it)->getParticle(1) == i2) {
	    delete *it;
	    it=interList.erase(it);
	    continue;
	}
	it++;
    }
    */

//=== end delete interactions ================================================ 

  //delete i1;
  //plist.erase(i1->plistIt());
  //delete i2;
  //plist.erase(i2->plistIt());

  eraseParticle(i1);
  eraseParticle(i2);

}

// This is called when particles are created from fluid.
//----------------------------------------------------------------------------
void Collision1::collisionUpdate(vector<EventParticle*>& outgoing,
	double itime, double gtime)
{
  initialTime=itime;
  finalTime=gtime;

  //ofstream out("particle.dat");
  for(int i=0; i<(int)outgoing.size();i++) {
    pnew.push_back(outgoing[i]);

    /*
    out  << setw(14) << outgoing[i]->getX()
	 << setw(14) << outgoing[i]->getY()
	 << setw(14) << outgoing[i]->getZ()
	 << setw(14) << outgoing[i]->getT()
	 <<endl;
	 */
  }
  //out.close();

  EventParIt i0 = ++pnew.begin();
  // Loop over newly produced particles.
  for(EventParIt i1 = i0; i1 != pnew.end(); ++i1)
  for(EventParIt i2 = pnew.begin(); i2 != i1; ++i2) {
    TwoBodyInterList* it = hit(*i1,*i2);
    if(it !=nullptr) interList.insert(it);
  }

  // Find new collisions between newly produced particle and old particles.
  for(EventParIt ip = pnew.begin(); ip != pnew.end(); ++ip) {
    if((*ip)->getStatus()== -1200) continue;
    for(EventParIt i3 = plist.begin(); i3 != plist.end(); ++i3) {
      TwoBodyInterList* it = hit(*ip,*i3);
      if(it !=nullptr) interList.insert(it);
    }
  }

  // check decay.
  if(decayOn) {
    for(EventParIt ip = pnew.begin(); ip != pnew.end(); ++ip) {
      double ct=(*ip)->lifetime();
      if(optTau) {
        Vec4 v = (*ip)->velocity(optPropagate);
        double zc= (*ip)->getZ() +v[3]*(ct-(*ip)->getT());
        if(abs(ct) < abs(zc)) {
          cout << "Collision1::collisionUpdate fluid tc= "<< ct << " zc= "<< zc<<endl;
          exit(1);
        }
        ct = sqrt(ct*ct - zc*zc);
      }
      if(ct < finalTime) interList.insert(new DecayList(*ip));
    }
  }

  cout << "Collision1::collisionUpdate inter size (+decay) = " << interList.size()<<endl;

  //plist.merge(pnew);
  for(auto p=pnew.begin();p != pnew.end();p++) {
    plist.push_front(*p);
    //(*p)->setPlistIt(plist.begin());
  }
  pnew.clear();
  
  //cout << " plist= "<< plist.size()<<endl;

}

//***********************************************************************
//...Purpose: to update collision array.
void Collision1::collisionUpdate(InterList* inter)
{
  int np = inter->getNumberOfInComing();
  if(np == 1 ) {
    numberOfDecay++;
    removeInteraction(inter->getParticle(0));
  } else {
    numberOfCollision++;
    removeInteraction2(inter);
  }

  if(pnew.size()==0) return;

  // Find next new collisions.
  for(EventParIt ip = pnew.begin(); ip != pnew.end(); ++ip) {
    if((*ip)->getStatus()== -1200) continue;
     for(EventParIt i3 = plist.begin(); i3 != plist.end(); ++i3) {
       TwoBodyInterList* it = hit(*ip,*i3);
       if(it !=nullptr) interList.insert(it);
     }
  }

  // check decay.
  if(decayOn) {
    for(auto& p : pnew) {
      double ct=p->lifetime();
      if(optTau) {
        Vec4 v = p->velocity(optPropagate);
        double zc= p->getZ() +v[3]*(ct-p->getT());
        if(abs(ct) < abs(zc)) {
          cout << "Collision1:collisionUpdate inter tc= "<< ct << " zc= "<< zc<<endl;
          cout << "t0= "<< p->getT() << " z0= "<< p->getZ() << " vz= "<< v[3] <<endl;
          exit(1);
        }
        ct = sqrt(ct*ct - zc*zc);
      }
      if(ct < finalTime) interList.insert(new DecayList(p));
    }
  }

  //plist.merge(pnew);
  for(auto p=pnew.begin();p != pnew.end();p++) {
    plist.push_front(*p);
    //(*p)->setPlistIt(plist.begin());
  }
  pnew.clear();

}

void Collision1::cancelCollision(InterList* inter)
{
  //double ctime = inter->getCollisionOrderTime();
  //int np = inter->getNumberOfInComing();
  //for(int i=0;i<np;i++) {
   // EventParticle* i1 = inter->getParticle(i);
    //i1->deleteInterList();
    //(*i1)->updateR(ctime);
  //}

  /*
  // Remove current collision list.
  InterList* current = *interList.begin();
  if(current != inter) {
      cout << " cancelcollision " << current << " inter= "<<inter<<endl;
      exit(1);
  }
  */

  // reset decay time.
  if(inter->getNumberOfInComing() == 1) {
    EventParticle *p= inter->getParticle(0);
    double m=p->getMass();
    double e=p->getE0();
    double tdec = p->lifetime() + jamParticleData->lifeTime(p->getParticleDataEntry(),m,e);
    p->setLifeTime(tdec);
    if(tdec < finalTime) interList.insert(new DecayList(p));
  }

  delete inter;
  interList.erase(interList.begin());

  // cancel newly produced particles.
  //list<EventParticle*>::iterator cp;
  //for(cp = pnew.begin(); cp != pnew.end(); ++cp) {
  //  delete *cp;
  //}
  //pnew.clear();

}

bool Collision1::checkNextCollisionTime(TwoBodyInterList* inter,double dtau1,double dtau2,
    bool preHadronA,bool preHadronB)
{
  //double tc1 = inter->getCollisionTime(0);
  //double tc2 = inter->getCollisionTime(1);
  EventParticle* i1=inter->getParticle(0);
  EventParticle* i2=inter->getParticle(1);

  for (auto it=interList.begin(); it != interList.end();it++ ) {
    TwoBodyInterList *iw=dynamic_cast<TwoBodyInterList*>(*it);
    if(iw==0) continue;
    if(iw==inter) continue;

    if(!preHadronA) {
      if(iw->getParticle(0) == i1 && iw->getCollisionTime(0) < dtau1) return true;
      if(iw->getParticle(1) == i1 && iw->getCollisionTime(1) < dtau1) return true;
    }
    if(!preHadronB) {
      if(iw->getParticle(0) == i2 && iw->getCollisionTime(0) < dtau2) return true;
      if(iw->getParticle(1) == i2 && iw->getCollisionTime(1) < dtau2) return true;
    }

  }
  return false;
}

} // end namespace jam2


