#ifndef jam2_collision_Collision3_h
#define jam2_collision_Collision3_h

#include <list>
#include <vector>
#include <set>
#include <fstream>
#include <jam2/collision/EventParticle.h>
#include <jam2/collision/TimeOrder.h>
#include <jam2/collision/InterList.h>
#include <jam2/xsection/CrossSection.h>
#include <Pythia8/Settings.h>
#include <jam2/collision/Collision.h>
#include <jam2/collision/CascBox.h>
#include <jam2/collision/CascCell.h>

namespace jam2 {

typedef std::list<EventParticle*>::iterator EventParIt;

class Collision3 : public Collision
{
private:
  //typedef std::multiset<InterList*,TimeOrder> InterListSet;
  //typedef InterListSet::iterator InterIt;
  //InterListSet interList;
 
  //typedef std::list<InterList*>::iterator InterIt;
  //std::list<InterList*> interList;

  CascCell* cascCell;
  std::vector<CascBox*> cascBox;
  int optBoxBoundary;
  /*
  int optInitCell;
  int maxX, maxY, maxZ, maxN;
  double dX, dY, dZ;
  double xMin, yMin, zMin, xMax, yMax, zMax;
  static const double infLength;
  double expansionStartTime;
  double vX, vY, vZ;
  */

public:
    Collision3(Pythia8::Settings* s,JamParticleData* jp, CrossSection* in,Pythia8::Rndm* r);
    ~Collision3();
  void clear();
  void init(const InitialCondition* initc);
    void makeCollisionList(double itime,double gtime);
    InterList* findNextCollision();
    void collisionUpdate();
    void collisionUpdate(InterList* inter);
    void collisionUpdate(vector<EventParticle*>& outgoing,double itime,double gtime);
    void cancelCollision(InterList* inter);
    /*
    double finalCollisionTime() {
      if(interList.size() ==0) return 0;
      return (*interList.begin())->getCollisionOrderTime();
    }
    */
  bool checkNextCollisionTime(TwoBodyInterList* inter,double dtau1,double dtau2,
                 bool preHadronA,bool preHadronB);
  void searchCollisionPair(CascBox* box);
  void searchCollisionPair(EventParticle* p1, CascBox* box);
  void checkWallCollision(EventParticle& ip, CascBox& box);
  void wallCollision(InterList& inter);

  //int site(int x, int y, int z) {return x + maxX*(y+maxY*z);}

  /*
  int inside(const Vec4& r) {
    double tw = std::max(0.0, r.e() - expansionStartTime);
    //int ix = std::min(std::max(0,int(floor((r.px()-xMin)/(dX+vX*tw)))),maxX-1);
    //int iy = std::min(std::max(0,int(floor((r.py()-yMin)/(dY+vY*tw)))),maxY-1);
    //int iz = std::min(std::max(0,int(floor((r.pz()-zMin)/(dZ+vZ*tw)))),maxZ-1);
    int ix = std::min(std::max(0,int(floor(r.px()/(dX+vX*tw)) + maxX/2)),maxX-1);
    int iy = std::min(std::max(0,int(floor(r.py()/(dY+vY*tw)) + maxY/2)),maxY-1);
    int iz = std::min(std::max(0,int(floor(r.pz()/(dZ+vZ*tw)) + maxZ/2)),maxZ-1);
    return  ix + maxX*(iy+maxY*iz);
    //return site(ix,iy,iz);
  }

  int boxPosition(const Vec4& r) {
    double tw = std::max(0.0, r.e() - expansionStartTime);
    int ix = int(floor(r.px()/(dX+vX*tw)) + maxX/2);
    int iy = int(floor(r.py()/(dY+vY*tw)) + maxY/2);
    int iz = int(floor(r.pz()/(dZ+vZ*tw)) + maxZ/2);
    if(ix<0 || ix>=maxX) return -1;
    if(iy<0 || iy>=maxY) return -1;
    if(iz<0 || iz>=maxZ) return -1;
    return ix + maxX*(iy+maxY*iz);
  }
  */

  void propagate(double ctime,int opt,int step=2);
  //TwoBodyInterList* hit2(EventParticle* i1, EventParticle* i2);
  bool doPauliBlocking(InterList* inter,vector<EventParticle*>& outgoing,int opt);
  void removeInteraction(EventParticle* i1);

private:
  void removeInteraction2(InterList* inter);

};
}
#endif
