#ifndef jam2_meanfield_MBoxS_h
#define jam2_meanfield_MBoxS_h

#include <vector>
#include <list>
#include <array>
#include <algorithm>
#include <jam2/meanfield/MBox.h>
#include <jam2/meanfield/RQMDs.h>

namespace jam2 {

class MBoxS;

class NeighborS
{
public:
  MBoxS* box;    // pointer to the neighbor box
  NeighborS*  neighbor; // pointer to the neighbor of the neighbor box
  std::vector<std::vector<double> > rhom;
private:
  bool action;
public:
  NeighborS(MBoxS* b) : box(b) {neighbor=nullptr;action=false;}
  void setAction() {action=true;}
  bool getAction() {return action;}
  void resize(int n);
  void setNeighbor(NeighborS* n) {neighbor=n;}
  void setRhom(int i, int j, double a) {rhom[i][j]=a;}
  //std::vector<EventParticle*>& getParticle() {return box->getParticle();}
 
};

class MBoxS : public MBox, public RQMDs
{
private:
  //int MV;
  //std::vector<MBoxS*> neighbors;
  std::vector<NeighborS*> neighbors;

public:
  MBoxS(std::array<double,3>& xmin, std::array<double,3>& xmax, std::array<int,3>& pos, bool ed,Pythia8::Settings* s) :
    MBox(xmin,xmax,pos,ed),RQMDs(s) { }//MV=0;}

  ~MBoxS();
  void initBox();
  void qmdMatrix(double a=1.0);
  void computeForce(double dt);
  void clearMatrix();

  void qmdMatrixNeighbor(NeighborS& b, double a=1.0);
  void computeForceNeighbor(NeighborS& b);

  void qmdMatrixNeighbor1(NeighborS& b, double a=1.0);
  void qmdMatrixNeighbor2(NeighborS& b, double a=1.0);
  void qmdMatrixNeighbor3(NeighborS& b, double a=1.0);
  void computeForceNeighbor1(NeighborS& b);
  void computeForceNeighbor2(NeighborS& b);
  void computeForceNeighbor3(NeighborS& b);

  void clear() {part.clear();}
  const std::vector<NeighborS*>& getNeighbors() const {return neighbors;}
  void addNeighbor(MBoxS* b) {neighbors.push_back(new NeighborS(b));}
  void setAction() {neighbors.back()->setAction();}
  bool haveThisSite(MBoxS* box) const;
  //{if(find(neighbors.begin(),neighbors.end(),box) != neighbors.end()) return true;
  //  else return false;}

  //std::vector<NeighborS>::iterator findNeighbor(MBoxS* box);
  NeighborS* findNeighbor(MBoxS* box);

};


} // end namespace jam2
#endif
