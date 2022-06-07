#include <jam2/meanfield/RQMDsBox.h>

// Relativistic quantum molecular dynamics with Skyrme potential
// potentials are implemented as scalar type.

namespace jam2 {

using namespace std;

RQMDsBox::RQMDsBox(Pythia8::Settings* s, InitialCondition* initcnd) : MCell(s),MeanField(s)
{
  initCell(initcnd);
  initBox(s);
}

// This is called before starting time evolution.
void RQMDsBox::init(std::list<EventParticle*>& plist)
{
  firstEng=true;

}

void RQMDsBox::initBox(Pythia8::Settings* s)
{
  std::array<double,3> xmin, xmax;
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

    mBox.push_back(new MBoxS(xmin,xmax,pos,edge,s));

  }

  // set neighbor box for the interaction between boxes
  for(auto& box : mBox ) {
    int x = box->x(0);
    int y = box->x(1);
    int z = box->x(2);
    for(int k=max(0,-1+z);k<min(2+z,maxZ);k++)
    for(int j=max(0,-1+y);j<min(2+y,maxY);j++)
    for(int i=max(0,-1+x);i<min(2+x,maxX);i++) {
      int ib = i + maxX*(j + maxY*k);
      if(i==x && j==y && k==z) continue; // exclude myself

      // include all neighbor cell
      box->addNeighbor(mBox[ib]);

      // exclude a cell to avoid duplicate actions.
      if(mBox[ib]->haveThisSite(box)) continue;
      box->setAction();
    }
  }

  for(auto& box : mBox) {
    for(auto& ng : box->getNeighbors()) {
      ng->setNeighbor( ng->box->findNeighbor(box) );
    }
  }

}

RQMDsBox::~RQMDsBox()
{
  for(auto& box : mBox) {
    delete box;
  }
  mBox.clear();
}

void RQMDsBox::addParticle(std::list<EventParticle*>& plist,double t,int step)
{
  eFree=0.0;
  for(auto& i : plist) {
    if(i->isMeanField(t,optPotential)) {
      mBox[inside(i->getR())]->add(i);
    } else {
      eFree += i->getP();
    }
  }
  if(step==1) eFree0 = eFree;

}

void RQMDsBox::evolution(list<EventParticle*>& plist, double t,double dt, int step)
{
  for(auto& b : mBox) b->clear();
  addParticle(plist,t,step);

  for(auto& b : mBox) b->initBox();
  for(auto& b : mBox) b->qmdMatrix();
  bool optSet= optPotentialArg == 0 ? true :  false;
  for(auto& b : mBox) b->singleParticlePotential(optSet);
  for(auto& b : mBox) b->computeForce(dt);

  if(!optSet) {
    for(auto& b : mBox) b->singleParticlePotential(optSet);
  }

  // Save initial total energy-momentum.
  if(firstEng) computeEnergy(plist,step);

  if(isDebug > 1) computeEnergy(plist,step);

  //for(auto& b : mBox) b->clearMatrix();

}

// compute single particle potential energy.
Vec4 RQMDsBox::computeEnergy(list<EventParticle*>& plist, int step)
{
  pTot=0.0;
  for(auto i = plist.begin(); i != plist.end(); ++i) {
    pTot += (*i)->getP();
  }

  if(firstEng) {pTot0 = pTot; firstEng=false;}

  if(isDebug > 1) {
  double econ=abs(pTot0[0]-pTot[0])/pTot0[0]*100;
  cout << "RQMDsBox time= " << step
      << " econ= " << fixed << econ << " %"
     << scientific << setw(13) << pTot[1] 
     << scientific << setw(13) << pTot[2]
     << scientific << setw(13) << pTot[3]
     <<endl;
  }

  return pTot;

}

} // namespace jam2


