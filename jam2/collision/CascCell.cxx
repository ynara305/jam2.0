#include <array>
#include <jam2/collision/CascCell.h>

namespace jam2 {

const double CascCell::infLength=1e+6;

CascCell::CascCell(Pythia8::Settings* s):settings(s)
{
  optInitCell = settings->mode("Cascade:optInitializeCellParameters");
  maxX = settings->mode("Cascade:nxCell");
  maxY = settings->mode("Cascade:nyCell");
  maxZ = settings->mode("Cascade:nzCell");

  vX = settings->parm("Cascade:vxCell");
  vY = settings->parm("Cascade:vyCell");
  vZ = settings->parm("Cascade:vzCell");

  dX = settings->parm("Cascade:dxCell");
  dY = settings->parm("Cascade:dyCell");
  dZ = settings->parm("Cascade:dzCell");
  expansionStartTime = settings->parm("Cascade:cellExpansionStartTime");
  optBoxBoundary=settings->mode("Cascade:boxBoundary");  // =0: boxes at the edge are infinitely large  =1: finite box

}

CascCell::~CascCell()
{
  for (auto& b : cascBox) delete b;
  cascBox.clear();
}

void CascCell::clear()
{
  int n=0;
  for(auto& box : cascBox) n += box->cleanInterList();

  if(n>0) {
    cout << "CascCell::clear() there are " << n << " interaction remain "<<endl;
  }
}

int CascCell::numberOfCollision()
{
  int ncoll = 0;
  for(auto& box : cascBox) ncoll += box->interListSize();
  return ncoll;
}

void CascCell::init(const InitialCondition* initcnd)
{
  // under construction.
  if(optInitCell==1) {
    double dx=sqrt(20.0/M_PI) + 1e-5;
    dX = dx, dY= dx, dZ = dx;

    // assume that simulation is done in the equal velocity system.
    //double eCM=initcnd->ecm();
    //double bMin=initcnd->bmin();
    double bMax=initcnd->bmax();
    //double zSep=initcnd->zseparation();

    double mA = initcnd->massA();
    double mB = initcnd->massB();
    double pzAcm  = initcnd->pzA();
    //double pzBcm  = initcnd->pzB();
    double eA = initcnd->energyA();
    double eB = initcnd->energyB();
    double gamA= eA/mA;
    double gamB= eB/mB;
    double vA = pzAcm/eA;
    //double vB = pzBcm/eB;
    int masA = initcnd->massNumberA();
    int masB = initcnd->massNumberB();
    double rA = initcnd->getRadius(masA);
    double rB = initcnd->getRadius(masB);

    double tpass = 2*rA/gamA/vA;
    expansionStartTime= std::min({dX,dY,dZ}) + tpass;
    double zmax = 2*rA/gamA + 2*rB/gamB;
    double xmax = rA + rB + bMax;
    double ymax = rA + rB;
    maxX = xmax / dX;
    maxY = ymax / dY;
    maxZ = std::max(6,(int)(zmax / dZ));
    vX = 2.0/ maxX;
    vY = 2.0/ maxY;
    vZ = 2.0/ maxZ;

  }

  int offx = maxX ==1 ? 0: maxX%2;
  int offy = maxY ==1 ? 0: maxY%2;
  int offz = maxZ ==1 ? 0: maxZ%2;
  xMin = -(maxX - offx)*dX/2;
  yMin = -(maxY - offy)*dY/2;
  zMin = -(maxZ - offz)*dZ/2;

  xMax = maxX*dX + xMin;
  yMax = maxY*dY + yMin;
  zMax = maxZ*dZ + zMin;
  maxN = maxX*maxY*maxZ -1;

  //xMin = -(maxX - maxX%2)*dX/2;
  //yMin = -(maxY - maxY%2)*dY/2;
  //zMin = -(maxZ - maxZ%2)*dZ/2;

  cout << "optBoxBoundary = "<<optBoxBoundary << endl;
  cout << " grid x= "<< maxX << " grid y= "<< maxY << " grid z= "<< maxZ<<endl;
  cout << " dx= "<< dX << " dy= "<< dY << " dz= "<< dZ<<endl;
  cout << " xMin= "<< xMin << " yMin= "<< yMin << " zMin= "<< zMin <<endl;
  cout << " xMax= "<< xMax << " yMax= "<< yMax << " zMax= "<< zMax <<endl;

  cout << "Box expansion start at "<< expansionStartTime << " with the velocity Vx= "<< vX << " Vy= "<< vY << " vZ= "<< vZ << endl;

  std::array<double,3> xmin,xmax;
  std::array<int,3> pos;
  double infLx =  infLength, infLy=  infLength, infLz=  infLength;
  double imfLx = -infLength, imfLy= -infLength, imfLz= -infLength;
  if(optBoxBoundary>0) {
    imfLx=xMin; imfLy=yMin; imfLz=zMin;
    infLx=xMax; infLy=yMax; infLz=zMax;
  }

  // make boxes for collision search.
  for(int z=0;z<maxZ;z++)
  for(int y=0;y<maxY;y++)
  for(int x=0;x<maxX;x++) {
    xmin[0]= x==0 ? imfLx : dX*x+xMin;
    xmin[1]= y==0 ? imfLy : dY*y+yMin;
    xmin[2]= z==0 ? imfLz : dZ*z+zMin;
    xmax[0]= x==maxX-1 ? infLx : dX*(x+1)+xMin;
    xmax[1]= y==maxY-1 ? infLy : dY*(y+1)+yMin;
    xmax[2]= z==maxZ-1 ? infLz : dZ*(z+1)+zMin;
    pos[0]=x; pos[1]=y; pos[2]=z;
    bool edgex = (x == 0) ? true : (x== maxX-1) ? true :false;
    bool edgey = (y == 0) ? true : (y== maxY-1) ? true :false;
    bool edgez = (z == 0) ? true : (z== maxZ-1) ? true :false;
    bool edge = false;
    if(edgex|| edgey || edgez ) edge = true;

    cascBox.push_back(new CascBox(settings,xmin,xmax,pos,edge));

  }

  // set neighbor box for collision search.
  for(auto& box : cascBox ) {
    int x = box->x(0);
    int y = box->x(1);
    int z = box->x(2);
    //cout << x + maxX*(y + maxY*z);
    for(int k=max(0,-1+z);k<min(2+z,maxZ);k++)
    for(int j=max(0,-1+y);j<min(2+y,maxY);j++)
    for(int i=max(0,-1+x);i<min(2+x,maxX);i++) {
      int ib = i + maxX*(j + maxY*k);
      box->addNeighbor2(cascBox[ib]);
      if(i==x && j==y && k==z) continue;
      if(cascBox[ib]->haveThisSite(box)) continue;
      box->addNeighbor1(cascBox[ib]);
    }
  }
}

double CascCell::staticWallCollisionTime(CascBox& box, Vec4& r, Vec4& v)
{
  double xbox = v[1] >0.0 ? box.xmax(0) : box.xmin(0);
  double ybox = v[2] >0.0 ? box.xmax(1) : box.xmin(1);
  double zbox = v[3] >0.0 ? box.xmax(2) : box.xmin(2);
  double tx = v[1] !=0.0 ? (xbox - r[1])/v[1] : infLength;
  double ty = v[2] !=0.0 ? (ybox - r[2])/v[2] : infLength;
  double tz = v[3] !=0.0 ? (zbox - r[3])/v[3] : infLength;
  double tw = std::min({tx,ty,tz}) +1e-6;
  return tw + r[0];
}


double CascCell::wallCollisionTime(const CascBox& box, const Vec4& r, const Vec4& v)
{
  double tp = r.e();
  double hX= tp > expansionStartTime ? vX/dX : 0.0;
  double hY= tp > expansionStartTime ? vY/dY : 0.0;
  double hZ= tp > expansionStartTime ? vZ/dZ : 0.0;
  double xB = 1.0 - hX*expansionStartTime;
  double yB = 1.0 - hY*expansionStartTime;
  double zB = 1.0 - hZ*expansionStartTime;

  double vxl=v.px() - hX*box.xmin(0);
  double vxr=v.px() - hX*box.xmax(0);
  double vyl=v.py() - hY*box.xmin(1);
  double vyr=v.py() - hY*box.xmax(1);
  double vzl=v.pz() - hZ*box.xmin(2);
  double vzr=v.pz() - hZ*box.xmax(2);

  double xp = -r.px() + v.px()*tp;
  double yp = -r.py() + v.py()*tp;
  double zp = -r.pz() + v.pz()*tp;

  double tx = vxr > 0.0 ? (box.xmax(0)*xB+xp)/vxr : vxl < 0.0 ? (box.xmin(0)*xB+xp)/vxl : infLength;
  double ty = vyr > 0.0 ? (box.xmax(1)*yB+yp)/vyr : vyl < 0.0 ? (box.xmin(1)*yB+yp)/vyl : infLength;
  double tz = vzr > 0.0 ? (box.xmax(2)*zB+zp)/vzr : vzl < 0.0 ? (box.xmin(2)*zB+zp)/vzl : infLength;
  double tc = std::min({tx,ty,tz}) + 1e-6;
  if(tc<= tp) {
    cout << "CascCell::wallColliisonTime tc<0? tc= "<< tc<<endl;
    cout << "r = "<< r;
    cout << "v= "<< v;
  }
  return tc;
}

InterList* CascCell::findNextCollision()
{
  InterList* inter=nullptr;
  double lastCollisionTime=1e+35;
  for(const auto& box : cascBox) {
    InterList* i = box->sortInterList();
    if(i==nullptr) continue;
    if(i->getCollisionOrderTime() < lastCollisionTime) {
      lastCollisionTime = i->getCollisionOrderTime();
      inter = i;
    }
  }
  return inter;
}

} // end namespace jam2
