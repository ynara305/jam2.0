#ifndef jam2_collision_CascCell_h
#define jam2_collision_CascCell_h

#include <vector>
#include <Pythia8/Basics.h>
#include <Pythia8/Settings.h>
#include <jam2/collision/CascBox.h>
#include <jam2/initcond/InitialCondition.h>
#include <jam2/collision/InterList.h>

namespace jam2 {

class CascCell
{
private:
  Pythia8::Settings* settings;
  int optInitCell;
  int optBoxBoundary;
  std::vector<CascBox*> cascBox;
  int maxX, maxY, maxZ, maxN;
  double dX, dY, dZ;
  double xMin, yMin, zMin, xMax, yMax, zMax;
  static const double infLength;
  double expansionStartTime;
  double vX, vY, vZ;

public:
  CascCell(Pythia8::Settings* s);
  ~CascCell();
  void clear();
  void init(const InitialCondition* initc);

  std::vector<CascBox*>& boxes() {return cascBox;}
  double wallCollisionTime(const CascBox& box, const Vec4& r, const Vec4& v);
  double staticWallCollisionTime(CascBox& box, Vec4& r, Vec4& v);
  int numberOfCollision();
  InterList* findNextCollision();

  //int site(int x, int y, int z) {return x + maxX*(y+maxY*z);}

  CascBox* box(const Vec4& r) const {
    double tw = std::max(0.0, r.e() - expansionStartTime);
    int ix = std::min(std::max(0,int(floor(r.px()/(dX+vX*tw)) + maxX/2)),maxX-1);
    int iy = std::min(std::max(0,int(floor(r.py()/(dY+vY*tw)) + maxY/2)),maxY-1);
    int iz = std::min(std::max(0,int(floor(r.pz()/(dZ+vZ*tw)) + maxZ/2)),maxZ-1);
    return cascBox[ix + maxX*(iy+maxY*iz)];
  }

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

};
}
#endif
