#ifndef jam2_meanfield_MCell_h
#define jam2_meanfield_MCell_h

#include <vector>
#include <Pythia8/Basics.h>
#include <Pythia8/Settings.h>
#include <jam2/initcond/InitialCondition.h>
//#include <jam2/meanfield/MBox.h>

namespace jam2 {

class MCell
{
protected:
  Pythia8::Settings* settings;
  int optInitCell;
  int optBoxBoundary;
  int maxX, maxY, maxZ, maxN;
  double dX, dY, dZ;
  double xMin, yMin, zMin, xMax, yMax, zMax;
  static const double infLength;
  double expansionStartTime;
  double vX, vY, vZ;

public:
  MCell(Pythia8::Settings* s);
  void initCell(const InitialCondition* initc);
  //void initNeighbor(std::vector<MBox*>& mBox);

  int inside(const Vec4& r) {
    double tw = std::max(0.0, r.e() - expansionStartTime);
    //int ix = std::min(std::max(0,int(floor((r.px()-xMin)/(dX+vX*tw)))),maxX-1);
    //int iy = std::min(std::max(0,int(floor((r.py()-yMin)/(dY+vY*tw)))),maxY-1);
    //int iz = std::min(std::max(0,int(floor((r.pz()-zMin)/(dZ+vZ*tw)))),maxZ-1);
    int ix = std::min(std::max(0,int(floor(r.px()/(dX+vX*tw)) + maxX/2)),maxX-1);
    int iy = std::min(std::max(0,int(floor(r.py()/(dY+vY*tw)) + maxY/2)),maxY-1);
    int iz = std::min(std::max(0,int(floor(r.pz()/(dZ+vZ*tw)) + maxZ/2)),maxZ-1);
    return  ix + maxX*(iy+maxY*iz);
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

};
}
#endif
