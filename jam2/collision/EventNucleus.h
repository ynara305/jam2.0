#ifndef jam2_collision_EventNucleus_H
#define jam2_collision_EventNucleus_H

#include <jam2/collision/EventParticle.h>

namespace jam2 {

class EventNucleus : public EventParticle
{
private:
  int iz=0;          // number of proton
  int in=0;          // number of neutron
  int iy=0;          // number of Lambda
  int is[3];         // number of Sigma is[0]=Sigma-, is[1]=Sigma0, is[2]=Sigma+
  int ix[2];         // number of Xi    ix[0]=Xi-, ix[1]=Xi0
  int ig=0;          // number of Omega-
  int   l;           // angular momentum
  //double ex;       // total energy of the nuclear cluster - mass (GeV)
  //int ie;
  //int nucl;
  //std::vector<EParticle*> pc;   // hadrons in the cluster

public:
  EventNucleus(int id0=0,int iz0=0,int in0=0,int iy0=0,int ism=0, int is0=0,int isp=0,int ixm=0,int ix0=0,int ig0=0):
     //iz(iz0),in(in0),iy(iy0), //is[0](ism), is[1](is0), is[2](isp), ix[0](ixm), ix[1](ix0), 
    //ig(ig0):
    EventParticle(id0) { 
      iz=iz0;in=in0;iy=iy0;
      is[0]=ism, is[1]=is0, is[2]=isp, ix[0]=ixm, ix[1]=ix0;
      ig=ig0;
    }
  void setAngMom(double j) {l=j;}
  //void add(EParticle* p)  {pc.push_back(p);}
  //std::vector<EParticle*> getPc() {return pc;}
  //int  getSize()          {return pc.size();}
  //EParticle* getPc(int i) {return pc[i];}
};

}
#endif
