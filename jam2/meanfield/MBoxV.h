#ifndef jam2_meanfield_MBoxV_h
#define jam2_meanfield_MBoxV_h

#include <vector>
#include <list>
#include <array>
#include <algorithm>
#include <jam2/meanfield/MBox.h>
#include <jam2/meanfield/RQMDv.h>

namespace jam2 {

class MBoxV;

class NeighborV
{
public:
  MBoxV* box;    // pointer to the neighbor box
  NeighborV*  neighbor; // pointer to the neighbor of the neighbor box
  std::vector<std::vector<double> > rhom;
private:
  bool action;
public:
  NeighborV(MBoxV* b) : box(b) {neighbor=nullptr;action=false;}
  void setAction() {action=true;}
  bool getAction() {return action;}
  void resize(int n);
  void setNeighbor(NeighborV* n) {neighbor=n;}
  void setRhom(int i, int j, double a) {rhom[i][j]=a;}
  //std::vector<EventParticle*>& getParticle() {return box->getParticle();}
 
};

class MBoxV : public MBox, public RQMDv
{
private:
  int MV;
  //std::vector<MBoxV*> neighbors;
  std::vector<NeighborV*> neighbors;

public:
  MBoxV(std::array<double,3>& xmin, std::array<double,3>& xmax, std::array<int,3>& pos, bool ed,Pythia8::Settings* s) :
    MBox(xmin,xmax,pos,ed),RQMDv(s) { MV=0;}

  ~MBoxV();
  void initBox();
  void qmdMatrix();
  void computeForce(double dt);
  void clearMatrix();

  void qmdMatrixNeighbor(NeighborV& b);
  void computeForceNeighbor(NeighborV& b);

  void qmdMatrixNeighbor1(NeighborV& b);
  void qmdMatrixNeighbor2(NeighborV& b);
  void qmdMatrixNeighbor3(NeighborV& b);
  void computeForceNeighbor1(NeighborV& b);
  void computeForceNeighbor2(NeighborV& b);
  void computeForceNeighbor3(NeighborV& b);
  Pythia8::Vec4 facV(int i, Vec4& p,NeighborV& b);
  void computeVdot(NeighborV& neighbor);

  void clear() {part.clear();}
  const std::vector<NeighborV*>& getNeighbors() const {return neighbors;}
  void addNeighbor(MBoxV* b) {neighbors.push_back(new NeighborV(b));}
  void setAction() {neighbors.back()->setAction();}
  bool haveThisSite(MBoxV* box) const;
  //{if(find(neighbors.begin(),neighbors.end(),box) != neighbors.end()) return true;
  //  else return false;}

  //std::vector<NeighborV>::iterator findNeighbor(MBoxV* box);
  NeighborV* findNeighbor(MBoxV* box);

};


} // end namespace jam2
#endif
