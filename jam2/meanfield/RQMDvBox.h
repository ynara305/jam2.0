#ifndef jam2_meanfield_RQMDvBox_h
#define jam2_meanfield_RQMDvBox_h

#include <jam2/meanfield/MBoxV.h>
#include <jam2/meanfield/MCell.h>

namespace jam2 {

class RQMDvBox : public MCell, public MeanField
{
private:
  std::vector<MBoxV*> mBox;
  int optBoxBoundary;
  bool firstEng;
  Vec4 eFree,eFree0;
  Vec4 pTot;
public:
  RQMDvBox(Pythia8::Settings* s, InitialCondition* ini);
  ~RQMDvBox();
  void evolution(std::list<EventParticle*>& plist,double t,double dt,int step);
  void init(std::list<EventParticle*>& plist);
  void singleParticlePotential();

  void initBox(Pythia8::Settings* s);
  void addParticle(std::list<EventParticle*>& plist,double dt, int step);
  //void initMatrix(std::list<EventParticle*>& plist,double dt, int step);
  //void deleteMatrix();

  Pythia8::Vec4 computeEnergy(std::list<EventParticle*>& plist, int step);

}; 
}  // namespace jam2
#endif
