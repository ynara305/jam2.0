#ifndef jam2_meanfield_RQMDsvBox_h
#define jam2_meanfield_RQMDsvBox_h

#include <jam2/meanfield/MBoxSV.h>
#include <jam2/meanfield/MCell.h>

namespace jam2 {

class RQMDsvBox : public MCell, public MeanField
{
private:
  std::vector<MBoxSV*> mBox;
  int optBoxBoundary;
  bool firstEng;
  Vec4 eFree,eFree0;
  Vec4 pTot;
public:
  RQMDsvBox(Pythia8::Settings* s, InitialCondition* ini);
  ~RQMDsvBox();
  void evolution(std::list<EventParticle*>& plist,double t,double dt,int step);
  void init(std::list<EventParticle*>& plist);
  void singleParticlePotential();

  void initBox(Pythia8::Settings* s);
  void addParticle(std::list<EventParticle*>& plist,double dt, int step);

  Pythia8::Vec4 computeEnergy(std::list<EventParticle*>& plist, int step);

}; 
}  // namespace jam2
#endif
