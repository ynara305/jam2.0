#ifndef jam2_meanfield_MBoxSV_h
#define jam2_meanfield_MBoxSV_h

#include <vector>
#include <list>
#include <array>
#include <algorithm>
#include <jam2/meanfield/MBox.h>
#include <jam2/meanfield/RQMDsv.h>

namespace jam2 {

class MBoxSV;

class NeighborSV
{
public:
  MBoxSV* box;    // pointer to the neighbor box
  NeighborSV*  neighbor; // pointer to the neighbor of the neighbor box
  std::vector<std::vector<double> > rhom;
private:
  bool action;
public:
  NeighborSV(MBoxSV* b) : box(b) {neighbor=nullptr;action=false;}
  void setAction() {action=true;}
  bool getAction() {return action;}
  void resize(int n);
  void setNeighbor(NeighborSV* n) {neighbor=n;}
  void setRhom(int i, int j, double a) {rhom[i][j]=a;}
  //std::vector<EventParticle*>& getParticle() {return box->getParticle();}
 
};

class MBoxSV : public MBox, public RQMDsv
{
private:
  //int MV;
  std::vector<NeighborSV*> neighbors;

public:
  MBoxSV(std::array<double,3>& xmin, std::array<double,3>& xmax, std::array<int,3>& pos, bool ed,Pythia8::Settings* s) :
    MBox(xmin,xmax,pos,ed),RQMDsv(s) { } //MV=0;}

  ~MBoxSV();
  void initBox();
  void qmdMatrix();
  void computeForce(double dt);
  void clearMatrix();

  void qmdMatrixNeighbor(NeighborSV& b);
  void computeForceNeighbor(NeighborSV& b);
  void computeBUUForceNeighbor(NeighborSV& b);

  //Pythia8::Vec4 facV(int i, Vec4& p,NeighborSV& b);
  void computeVdot(NeighborSV& neighbor);

  void clear() {part.clear();}
  const std::vector<NeighborSV*>& getNeighbors() const {return neighbors;}
  void addNeighbor(MBoxSV* b) {neighbors.push_back(new NeighborSV(b));}
  void setAction() {neighbors.back()->setAction();}
  bool haveThisSite(MBoxSV* box) const;
  //{if(find(neighbors.begin(),neighbors.end(),box) != neighbors.end()) return true;
  //  else return false;}

  //std::vector<NeighborV>::iterator findNeighbor(MBoxSV* box);
  NeighborSV* findNeighbor(MBoxSV* box);

};


} // end namespace jam2
#endif
