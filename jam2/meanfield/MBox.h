#ifndef jam2_meanfield_MBox_h
#define jam2_meanfield_MBox_h

#include <vector>
#include <array>
#include <algorithm>
#include <jam2/collision/EventParticle.h>

namespace jam2 {

class MBox
{
protected:
  std::array<double,3> xMin,xMax;
  std::array<int,3> position;
  //std::vector<MBox*> neighbors;
  bool edgeBox;

public:
  MBox(std::array<double,3>& xmin, std::array<double,3>& xmax, std::array<int,3>& pos, bool ed):
    xMin(xmin), xMax(xmax), position(pos), edgeBox(ed) {}

  bool edge() const {return edgeBox;}
  int x(int i) const {return position[i];}
  double xmin(int i) const {return xMin[i];}
  double xmax(int i) const {return xMax[i];}
  std::array<double,3> xmin() const {return xMin;}
  std::array<double,3> xmax() const {return xMax;}
  std::array<int,3> pos() const {return position;}

  //virtual void init()=0;
  //virtual void qmdMatrix(double a=1.0)=0;
  //virtual void computeForce()=0;
  //virtual void propagate(double dt)=0;
  //virtual void clearMatrix()=0;

  //virtual const std::vector<MBox*>& getNeighbors() =0;
  //virtual void add(EventParticle* p)=0;
  //virtual void addNeighbor(MBox* b)=0;
  //virtual bool haveThisSite(MBox* box) const=0;

};


} // end namespace jam2
#endif
