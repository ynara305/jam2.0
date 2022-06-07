#include <jam2/initcond/Nucleus.h>

namespace jam2 {

using namespace std;
using namespace Pythia8;

Vec4 Nucleus::generateHardSphare() const
{
  double r = rHard * pow(rndPtr->flat(), 1.0/3.0);
  double costhe = 2.0*rndPtr->flat() - 1.0;
  double sinthe = sqrt(max(1.0 - costhe*costhe, 0.0));
  double phi = 2.0*M_PI*rndPtr->flat();

  return Vec4(r*sinthe*cos(phi), r*sinthe*sin(phi), r*costhe);

}

// deuteron
Vec4 Nucleus::generateDeuteron() const
{
  const double gamma = 0.2316;  // [1/fm]
  const double beta = 1.385;    // [1/fm]
  double phi_d =0.0, r=0.0;
  // max of nucleon prob density
  double rwMax = pow(beta-gamma,2);
  do { // Hulthen wave function
    r = 10.*pow(rndPtr->flat(),1.0/3.0);  // trial r, [fm]
    phi_d = (r>1e-3) ? (exp(-gamma*r) - exp(-beta*r))/r : beta-gamma;
    r /= 2.0;         // r in Hulthen wf is relative distance of n and p
  } while(rndPtr->flat()*rwMax > phi_d*phi_d);
  double cx = 1.0-2.0*rndPtr->flat();
  double sx = sqrt(1.0-cx*cx);
  double phi = 2*M_PI*rndPtr->flat();
  return Vec4(r*sx*cos(phi),r*sx*sin(phi),r*cx);
}


// This is the same as  Pythia8::GLISSANOModel::generate() except pz().
vector<Nucleon> Nucleus::generate() const
{
  //vector<Nucleon> nucl=GLISSANDOModel::generate();

  int sign = id() > 0? 1: -1;
  int pid = sign*2212;
  int nid = sign*2112;
  vector<Nucleon> nucleons;
  if ( A() == 0 ) {
    nucleons.push_back(Nucleon(id(), 0, Vec4()));
    return nucleons;
  } else if ( A() == 1 ) {
    if ( Z() == 1 ) nucleons.push_back(Nucleon(pid, 0, Vec4()));
    else  nucleons.push_back(Nucleon(nid, 0, Vec4()));
    return nucleons;
  } else if( A() == 2 ) {
    Vec4 pos = generateDeuteron();
    nucleons.push_back(Nucleon(pid,0,pos));
    nucleons.push_back(Nucleon(nid,1,Vec4(-pos.px(),-pos.py(),-pos.pz())));
    return nucleons;
  }

  Vec4 cms;
  vector<Vec4> positions;
  while ( int(positions.size()) < A() ) {
    while ( true ) {
      Vec4 pos;
      if(optSample == 1) 
        pos = generateNucleon();
      else
        pos = generateHardSphare();

      bool overlap = false;
      for ( int i = 0, N = positions.size(); i < N && !overlap; ++i )
        if ( (positions[i] - pos).pAbs() < (gaussHardCore ? RhGauss() : Rh()) )
          overlap = true;
      if ( overlap ) continue;
      positions.push_back(pos);
      cms += pos;
      break;
    }
  }

  cms /= A();
  nucleons.resize(A());
  int Np = Z();
  int Nn = A() - Z();
  for ( int i = 0, N= positions.size(); i < N; ++i ) {
    Vec4 pos(positions[i].px() - cms.px(),
             positions[i].py() - cms.py(),
             positions[i].pz() - cms.pz());//This part is the only modification.
    if ( int(rndPtr->flat()*(Np + Nn)) >= Np ) {
      --Nn;
      nucleons[i] = Nucleon(nid, i, pos);
    } else {
      --Np;
      nucleons[i] = Nucleon(pid, i, pos);
    }
  }

  return nucleons;

}


}

