#include <array>
#include <jam2/meanfield/MCell.h>

namespace jam2 {

const double MCell::infLength=1e+6;

MCell::MCell(Pythia8::Settings* s):settings(s)
{
  optInitCell = settings->mode("MeanField:optInitializeCellParameters");
  maxX = settings->mode("MeanField:nxCell");
  maxY = settings->mode("MeanField:nyCell");
  maxZ = settings->mode("MeanField:nzCell");

  vX = settings->parm("MeanField:vxCell");
  vY = settings->parm("MeanField:vyCell");
  vZ = settings->parm("MeanField:vzCell");

  dX = settings->parm("MeanField:dxCell");
  dY = settings->parm("MeanField:dyCell");
  dZ = settings->parm("MeanField:dzCell");
  expansionStartTime = settings->parm("MeanField:cellExpansionStartTime");
  optBoxBoundary=settings->mode("MeanField:boxBoundary");  // =0: boxes at the edge are infinitely large  =1: finite box

}

void MCell::initCell(const InitialCondition* initcnd)
{
  // under construction.
  if(optInitCell==1) {
    double dx=3.0;
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

  cout << "MBox optBoxBoundary = "<<optBoxBoundary << endl;
  cout << " grid x= "<< maxX << " grid y= "<< maxY << " grid z= "<< maxZ<<endl;
  cout << " dx= "<< dX << " dy= "<< dY << " dz= "<< dZ<<endl;
  cout << " xMin= "<< xMin << " yMin= "<< yMin << " zMin= "<< zMin <<endl;
  cout << " xMax= "<< xMax << " yMax= "<< yMax << " zMax= "<< zMax <<endl;

  cout << "Meanfield Box expansion start at "<< expansionStartTime << " with the velocity Vx= "<< vX << " Vy= "<< vY << " vZ= "<< vZ << endl;

}

} // end namespace jam2
