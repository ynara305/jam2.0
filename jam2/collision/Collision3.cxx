#include <vector>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <algorithm>

#include <jam2/collision/Collision3.h>
#include <jam2/collision/InterType.h>
#include <jam2/collision/DecayList.h>
#include <jam2/collision/WallList.h>
#include <jam2/xsection/XsecTable.h>

// Box is used. interList is stored in each box.

namespace jam2 {

using namespace std;

void Collision3::clear()
{
  cascCell->clear();
  clearPlist();
}

void Collision3::init(const InitialCondition* initcnd)
{
  cascCell = new CascCell(settings);
  cascCell->init(initcnd);
  cascBox = cascCell->boxes();
}

Collision3::Collision3(Pythia8::Settings* s, JamParticleData* jp,CrossSection* in,Pythia8::Rndm* r)
	: Collision(s,jp,in,r)
{
  optBoxBoundary=settings->mode("Cascade:boxBoundary");  // =0: boxes at the edge are infinitely large  =1: finite box

}

Collision3::~Collision3()
{
  delete cascCell;
}

// =========================================================================//
//...Purpose: to make collision list into collision predictor arrays.
void Collision3::makeCollisionList(double itime, double gtime)
{
  //for(auto& box : cascBox) box->cleanInterList();

  initialTime=itime;
  finalTime=gtime;

  // assign a box to a particle and check wall collision and decay.
  for(const auto& i1 : plist ) {

    if(optBoxBoundary>0) {
    if(cascCell->boxPosition(i1->getR()) == -1) {
      cout << " makeCollisionList outside ? " << i1->getR()<<endl;
      exit(1);
    }
    }

    CascBox *box = cascCell->box(i1->getR());
    i1->addBox(box);
    checkWallCollision(*i1,*box);

  }

  plist.clear();

  // Loop over all boxes.
  for(const auto& box : cascBox) searchCollisionPair(box);

  predictedColl = cascCell->numberOfCollision();
  numberOfCollision=0;

}

void Collision3::searchCollisionPair(CascBox* box)
{
  auto& pl = box->getParticles();
  EventParIt i0 = ++pl.begin();
  for(EventParIt i1 = i0; i1 != pl.end(); ++i1)
  for(EventParIt i2 = pl.begin(); i2 != i1; ++i2) {
    TwoBodyInterList* it = hit(*i1,*i2);
    if(it !=nullptr) box->setInterList(it);
 }

  // search for collision in the neighbor grid.
  for(auto& neighbor : box->getNeighbors1()) {
    auto const& pl2=neighbor->getParticles();
    for(auto& p1:pl)
    for(auto& p2:pl2) {
      TwoBodyInterList* it = hit(p1,p2);
      if(it !=nullptr) box->setInterList(it);
    }
  }
}

//***********************************************************************
//...Purpose: to search decay and collision arrays to find next event
InterList* Collision3::findNextCollision()
{
  return cascCell->findNextCollision();


  if(cascBox.size() < plist.size()) {

    return cascCell->findNextCollision();

  } else {

  InterList* inter=0;
  double lastCollTime=1e+25;
  for(const auto& p : plist) {
    InterList* i = p->box()->sortInterList();
    if(i==nullptr) continue;
    if(i->getCollisionOrderTime() < lastCollTime) {
      lastCollTime = i->getCollisionOrderTime();
      inter = i;
    }
  }
  return inter;

  }
}

void Collision3::checkWallCollision(EventParticle& ip, CascBox& box)
{
  Vec4 r = ip.getR();
  Vec4 v = ip.velocity(optPropagate);
  double tw = cascCell->wallCollisionTime(box,r,v);

  if(tw<=r[0] || tw < initialTime) {
    Vec4 r1 = r;
    r1 += v*(tw - r[0]);
    r1[0]= tw;
    cout << "Collision3::checkWallCollision tw strange t-t0 "<< setprecision(20) << tw-r[0]<<endl;
    cout << " t-t0 "<< scientific << setprecision(20) << tw-r[0]<<endl;
    cout << scientific << setprecision(9) << " tw = "<< tw
          << " t0 = "<< r[0]
          << " ini t = "<< initialTime <<endl;
    cout << "p = "<<  &ip << " id= "<< ip.getID()<<endl;
    cout << " v= "<< v ;
    cout << "r0= "<< r <<endl;
    cout << "r1= "<< r1 <<endl;
    cout << "p= "<< ip.getP() <<endl;
    ip.setTWall(1e+6);
    exit(1);
    return;
  }

  if(optBoxBoundary>0) {
    Vec4 r1 = r;
    r1 += v*(tw - r[0]);
    r1[0]= tw;
    if(cascCell->boxPosition(r1)==-1) {
    cout << " going outside of wall? tw = "<< tw
      << " r0= "<< r
      << " r1= "<< r1 << endl;
    cout << " v= "<< v;
    cout << " p= "<< ip.getP();
    cout << " ipos orig= "<< cascCell->boxPosition(r) <<endl;
    exit(1);
    }
  }

  double tauw=tw;
  // in case of the option of tau ordering
  if(optTau) {
    double zc= r[3]+v[3]*(tw-r[0]);
    if(tw*tw < zc*zc) {
      cout << "Collision3:checkWallCollision tw= "<< tw << " zc= "<< zc<<endl;
      cout << "t0= "<< r[0] << " z0= "<< r[3] << " v= "<< v[3] <<endl;
      cout << " dt "<< tw - r[0] << " dt*v= " << (tw-r[0])*v[3] <<endl;
      cout << " t0 "<< r[0]  << " zv= " << r[3]*v[3] <<endl;
      //exit(1);
    }
    if(abs(tw) >= abs(zc)) {
      tauw = sqrt(tw*tw - zc*zc);
    }
  }

  // set time of wall collision.
  ip.setTWall(tw);

  if(decayOn) {
    double td = ip.lifetime();
    if(optTau && td < 1e+10) {
      double zd= r[3]+v[3]*(td-r[0]);
      if(abs(td) < abs(zd)) {
        cout << "Collision3:checkWallCollision td= "<< td << " zd= "<< zd <<endl;
      } else {
        td = sqrt(td*td - zd*zd);
      }
    }
    if(tauw < td) {
      if(tauw < finalTime) box.setInterList(new WallList(&ip,tauw));
    } else if(td < finalTime) {
      box.setInterList(new DecayList(&ip));
    }
  } else {
    if(tauw < finalTime) box.setInterList(new WallList(&ip,tauw));
  }

}

void Collision3::wallCollision(InterList& inter)
{
  numberOfWallCollision++;
  EventParticle* i1=inter.getParticle(0);
  // ordering time.
  //double tw=inter.getCollisionOrderTime();
  // The time a particle hits the wall.
  double tw=i1->tWall();

  CascBox *box0= i1->box();
  i1->updateR(tw,optPropagate);

  //i1->setVertex();
  i1->setTimeLastColl();

  CascBox* box1= cascCell->box(i1->getR());
  i1->changeBox(box1);
  box0->cleanInterList(i1);
  checkWallCollision(*i1,*box1);
  searchCollisionPair(i1,box1);

}

void Collision3::searchCollisionPair(EventParticle* p1, CascBox* box)
{
  //setParticle1(p1);
  // find collision with particles including the neighbor grid.
  for(const auto& neighbor : box->getNeighbors2()) {
    // loop over all particle in this box.
    for(const auto& p2:neighbor->getParticles()) {
      //TwoBodyInterList* it = hit(p2);
      TwoBodyInterList* it = hit(p1,p2);
      if(it !=nullptr) box->setInterList(it);
    }
  }

}

//...Remove collisions including i1,i2.
void Collision3::removeInteraction2(InterList* inter)
{
  EventParticle* i1 = inter->getParticle(0);
  EventParticle* i2 = inter->getParticle(1);
  CascBox* box1= i1->box();
  CascBox* box2= i2->box();
  if(box1==0 || box2==0) { 
    cout << "Collision3::removeInteraction2 box1= "<< box1 << " box2= "<< box2<<endl;
    exit(1);
  }

  for(auto& neighbor : box1->getNeighbors2()) {
    neighbor->cleanInterList(i1,i2);
  }

  if(box1 != box2) {
    for(auto& neighbor : box2->getNeighbors2()) neighbor->cleanInterList(i1,i2);
  }

  //i1->removeBox();
  box1->eraseParticle(i1);

  //i2->removeBox();
  box2->eraseParticle(i2);

}

void Collision3::removeInteraction(EventParticle* i1)
{
  for(auto& neighbor : i1->box()->getNeighbors2()) {
    neighbor->cleanInterList(i1);
  }
  //i1->removeBox();
  i1->box()->eraseParticle(i1);
}

//***********************************************************************
//...Purpose: to update collision array.
void Collision3::collisionUpdate(InterList* inter)
{
  int np = inter->getNumberOfInComing();
  if(np == 1 ) {
    numberOfDecay++;
    auto&& i1 = inter->getParticle(0);
    //removeInteraction(i1);

    for(auto& neighbor : i1->box()->getNeighbors2()) {
      neighbor->cleanInterList(i1);
    }
    //i1->removeBox();
    i1->box()->eraseParticle(i1);

    if(inter->getParticle(1) != nullptr) {
      cout << "collisonUpdate strange "<<endl;
      exit(1);
    }

  } else if(np ==2) {
    numberOfCollision++;
    removeInteraction2(inter);
  } else {
    cout << "Collision3::collisionUpdate wrong np np= "<< np<< endl;
    exit(1);
  }

  if(pnew.size()==0) return;

  collisionUpdate();

}

void Collision3::collisionUpdate()
{
  // Find next new collisions for newly produced particles.
  for(auto&& ip : pnew) {

    // check if this particle is outside the cell.
    if(optBoxBoundary >0) {
      if(cascCell->boxPosition(ip->getR())==-1) {
      Vec4 r=ip->getR();
      cout << "Collision3::collisionUpdate particle is produced outside box r = "<< r;
      exit(1);
      }
    }

    // find a box
    CascBox *box = cascCell->box(ip->getR());

    // search wall collision and particle decay time.
    checkWallCollision(*ip,*box);

    // This particle will be a fluid after formation time.
    if(ip->getStatus()== -1200) {
      ip->addBox(box);
      continue;
    }

    // search collisions
    searchCollisionPair(ip,box);

    // put this new particle into a box after search collision
    // to avoid collision between newly produced particles.
    ip->addBox(box);

  }

  pnew.clear();

  /*
  // merge pnew into plist.
  for(auto& p:pnew) {
    plist.push_front(p);
    p->setPlistIt(plist.begin());
  }
  pnew.clear();
  */

}

// This is called when particles are created from fluid.
//----------------------------------------------------------------------------
void Collision3::collisionUpdate(std::vector<EventParticle*>& outgoing,
	double itime, double gtime)
{
  initialTime=itime;
  finalTime=gtime;

  for(int i=0; i<(int)outgoing.size();i++) {
    pnew.push_back(outgoing[i]);
  }

  EventParIt i0 = ++pnew.begin();
  // Loop over newly produced particles.
  for(EventParIt i1 = i0; i1 != pnew.end(); ++i1) {
    CascBox *box = cascBox[cascCell->inside((*i1)->getR())];
    for(EventParIt i2 = pnew.begin(); i2 != i1; ++i2) {
      // if i2 is far away from i1 we should not search collision between i1 and i2.
      if(box->isNeighbor(**i2)) {
	TwoBodyInterList* it = hit(*i1,*i2);
        if(it !=nullptr) box->setInterList(it);
      }
    }
  }

  // collision search between newly created particle and old particles.
  collisionUpdate();

}


void Collision3::cancelCollision(InterList* inter)
{
  EventParticle* p = inter->getParticle(0);

  // reset decay time.
  if(inter->getNumberOfInComing() == 1) {
    double m=p->getMass();
    double e=p->getE0();
    double tdec = p->lifetime() + jamParticleData->lifeTime(p->getParticleDataEntry(),m,e);
    p->setLifeTime(tdec);
    if(tdec < finalTime) p->box()->setInterList(new DecayList(p));
  }


  p->box()->removeInter(inter);
  return;

}

bool Collision3::checkNextCollisionTime(TwoBodyInterList* inter,double dtau1,double dtau2,
    bool preHadronA,bool preHadronB)
{
  EventParticle* i1=inter->getParticle(0);
  EventParticle* i2=inter->getParticle(1);
  CascBox *box1 = i1->box();
  CascBox *box2 = i2->box();

  if (!preHadronA && !preHadronB && box1 == box2) {
    for (const auto& neighbor: box1->getNeighbors2())
      if (neighbor->checkNextCollisionTime(static_cast<InterList*>(inter), i1, dtau1, i2, dtau2))
	return true;
    return false;
  }

  if (!preHadronA) {
    for(const auto& neighbor: box1->getNeighbors2())
      if (neighbor->checkNextCollisionTime(static_cast<InterList*>(inter), i1, dtau1))
	return true;
  }

  if (!preHadronB) {
    for(const auto& neighbor: box2->getNeighbors2())
      if (neighbor->checkNextCollisionTime(static_cast<InterList*>(inter), i2, dtau2))
	return true;
  }

  return false;
}

void Collision3::propagate(double ctime, int opt, int step)
{
  for(auto& b : cascBox) {
    if(b->getParticles().size()>0) {
    plist.splice(plist.end(), b->getParticles());
    }
  }
  if(step==1) return;

  for(auto& i : plist) {
    if(ctime > i->getT()) i->updateR(ctime,opt);
    i->clearBox();
  }

}

bool Collision3::doPauliBlocking(InterList* inter,
	vector<EventParticle*>& outgoing,int opt)
{
  int np = inter->getNumberOfInComing();
  EventParticle* i1 = inter->getParticle(0);
  EventParticle* i2 = np == 2 ? inter->getParticle(1):i1;

  for(const auto& ip : outgoing) {
    int idp = ip->getID();
    if(idp != 2212 && idp != 2112) continue;
    double ctime = ip->getT();
    if(ctime < ip->getTf()) continue; // not formed
    CascBox *box = cascCell->box(ip->getR());
    double phase = 0.0;
    for(const auto& neighbor : box->getNeighbors2())
      phase += neighbor->phaseSpace(i1,i2,ip,ctime,opt);

    // Loop over all boxes.
    //for(const auto& box : cascBox) phase += box->phaseSpace(i1,i2,ip,ctime,opt);


    //cout << " phase = "<< pauliC*phase <<endl;
    //cin.get();
   if(pauliC*phase/overSample > rndm->flat()) return true; // Pauli-blocked.
  }

  // No Pauli blocking.
  return false;
}

} // end namespace jam2
