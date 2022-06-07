#ifndef jam2_interaction_NuclearCluster_h
#define jam2_interaction_NuclearCluster_h

#include <jam2/collision/EventParticle.h>
//#include "MassTable.h"

#include <vector>

namespace jam2 {

using Pythia8::Vec4;

class NuclearCluster
{
private:
  double R0=3.8; // fm
  double P0=0.3; // GeV/c
  //MassTable* mastbl;
public:
  NuclearCluster(double r0=3.8, double p0=0.3):R0(r0), P0(p0) {
  cout << "nuclear cluster parameters: R0= "<< R0 << " P0= " << P0<<endl; }
  void setR0(double r) {R0=r;}
  void setP0(double p) {P0=p;}
  void setRP0(double r, double p) {R0=r,P0=p;}
  bool clust(EventParticle* p1, EventParticle* p2);
  int findCluster(std::list<EventParticle*>& plist);
  void rpCMsq(EventParticle &p1, EventParticle &p2,double& rr,double &pp);
  EventParticle* setClusterProperty(std::vector<EventParticle*>& part, std::vector<int>& num,
      int& ii, int nt,int isave, int& ibar,Pythia8::Vec4& pf);
  double getMass(int kf);

};
}
#endif
