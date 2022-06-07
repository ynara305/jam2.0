#ifndef jam2_collision_CascBox_h
#define jam2_collision_CascBox_h

#include <vector>
#include <list>
#include <array>
#include <memory>
#include <algorithm>
#include <Pythia8/Settings.h>

namespace jam2 {

class EventParticle;
class InterList;
class InterListContainer;

class CascBox
{
private:
  Pythia8::Settings* settings;
  std::array<double,3> xMin,xMax;
  std::array<int,3> position;
  bool edgeBox;
  std::list<EventParticle*> particles;
  std::vector<CascBox*> neighbors1;  // neighbor excluding duplicate action
  std::vector<CascBox*> neighbors2;  // all neighbor cell including myself
  std::unique_ptr<InterListContainer> m_interList; // PImpl
  int     withBox;
  double  xBox,yBox,zBox;
  double gWidth,pauliR, pauliP;

public:
  CascBox(Pythia8::Settings* s,std::array<double,3>& xmin, std::array<double,3>& xmax, std::array<int,3>& pos, bool ed);
  ~CascBox();

  std::list<EventParticle*>& getParticles()  {return particles;}
  const std::vector<CascBox*>&     getNeighbors1() {return neighbors1;}
  const std::vector<CascBox*>&     getNeighbors2() {return neighbors2;}

  bool edge() const {return edgeBox;}
  int x(int i) const {return position[i];}
  double xmin(int i) const {return xMin[i];}
  double xmax(int i) const {return xMax[i];}
  std::array<double,3> xmin() const {return xMin;}
  std::array<double,3> xmax() const {return xMax;}
  std::array<int,3> pos() const {return position;}
  void addNeighbor1(CascBox* b) {neighbors1.push_back(b);}
  void addNeighbor2(CascBox* b) {neighbors2.push_back(b);}
  bool haveThisSite(CascBox* box) const {
    if(find(neighbors1.begin(),neighbors1.end(),box) != neighbors1.end()) return true;
    else return false;
  }
  std::list<EventParticle*>::iterator add(EventParticle* p) {
    particles.push_front(p);
    return particles.begin();
  }
  void removeParticle(std::list<EventParticle*>::iterator p) {
    particles.erase(p);
  }
  void removeParticle(EventParticle* ip) {
    particles.erase(std::find(particles.begin(),particles.end(),ip));
  }
  void eraseParticle(EventParticle* ip);

  void clearParticle();
  void propagate(double ctime, int opt, int step=2);

  void setInterList(InterList* inter);
  void cleanInterList(EventParticle* i1, EventParticle* i2);
  void cleanInterList(EventParticle* i1);
  void removeInterList(EventParticle* i1);
  int cleanInterList();
  void removeInter(InterList* inter);
  InterList* sortInterList();
  void printInterList();
  int interListSize();
  bool checkNextCollisionTime(InterList const* inter,
    EventParticle const* p1, double dtau1,
    EventParticle const* p2 = nullptr, double dtau2 = 0.0) const;

  bool isNeighbor(const EventParticle& p);
  bool isNeighbor(const CascBox& box);
  double phaseSpace(EventParticle* i1, EventParticle* i2, EventParticle* ip,double ctime,int opt);

  //void searchCollisionPair();
  //void searchCollisionPair(EventParticle* p1);
};


} // end namespace jam2
#endif
