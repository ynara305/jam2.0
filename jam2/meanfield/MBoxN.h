#ifndef jam2_meanfield_MBoxN_h
#define jam2_meanfield_MBoxN_h

#include <vector>
#include <list>
#include <array>
#include <algorithm>
#include <jam2/meanfield/MBox.h>
#include <jam2/meanfield/QMD.h>

namespace jam2 {

class MBoxN;

class NeighborN
{
public:
  MBoxN* box;    // pointer to the neighbor box
  NeighborN*  neighbor; // pointer to the neighbor of the neighbor box
  std::vector<std::vector<double> > rhom;
private:
  bool action;
public:
  NeighborN(MBoxN* b) : box(b) {neighbor=nullptr;action=false;}
  void setAction() {action=true;}
  bool getAction() {return action;}
  void resize(int n);
  void setNeighbor(NeighborN* n) {neighbor=n;}
  void setRhom(int i, int j, double a) {rhom[i][j]=a;}
  //std::vector<EventParticle*>& getParticle() {return box->getParticle();}
 
};

class MBoxN : public MBox, public QMD
{
private:
  int MV;
  std::vector<NeighborN*> neighbors;

public:
  MBoxN(std::array<double,3>& xmin, std::array<double,3>& xmax, std::array<int,3>& pos, bool ed,Pythia8::Settings* s) :
    MBox(xmin,xmax,pos,ed),QMD(s) { MV=0;}

  ~MBoxN();
  void initBox();
  void qmdMatrix();
  void computeForce(double dt);
  void clearMatrix();

  void qmdMatrixNeighbor(NeighborN& b);
  void computeForceNeighbor(NeighborN& b);

  void clear() {part.clear();}
  const std::vector<NeighborN*>& getNeighbors() const {return neighbors;}
  void addNeighbor(MBoxN* b) {neighbors.push_back(new NeighborN(b));}
  void setAction() {neighbors.back()->setAction();}
  bool haveThisSite(MBoxN* box) const;
  NeighborN* findNeighbor(MBoxN* box);

};


} // end namespace jam2
#endif
