#ifndef jam2_meanfield_QMDBox_h
#define jam2_meanfield_QMDBox_h

//#include <jam2/meanfield/RQMDs.h>
#include <jam2/meanfield/MBoxN.h>
#include <jam2/meanfield/MCell.h>

namespace jam2 {

class QMDBox : public MCell, public MeanField
{
private:
  std::vector<MBoxN*> mBox;
  int optBoxBoundary;
public:
  QMDBox(Pythia8::Settings* s, InitialCondition* ini);
  ~QMDBox();
  void evolution(std::list<EventParticle*>& plist,double t,double dt,int step);
  Pythia8::Vec4 computeEnergy(std::list<EventParticle*>& plist, int step);
  void init(std::list<EventParticle*>& plist);
  void singleParticlePotential();

  void initBox(Pythia8::Settings* s);
  void addParticle(std::list<EventParticle*>& plist,double t);
  //void initMatrix(std::list<EventParticle*>& plist,double dt, int step);
  //void deleteMatrix();


}; 
}  // namespace jam2
#endif
