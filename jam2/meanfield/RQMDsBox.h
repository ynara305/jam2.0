#ifndef jam2_meanfield_RQMDsBox_h
#define jam2_meanfield_RQMDsBox_h

//#include <jam2/meanfield/RQMDs.h>
#include <jam2/meanfield/MBoxS.h>
#include <jam2/meanfield/MCell.h>

namespace jam2 {

class RQMDsBox : public MCell, public MeanField
{
private:
  std::vector<MBoxS*> mBox;
  int optBoxBoundary;
  bool firstEng;
  Vec4 eFree,eFree0;
  Vec4 pTot;
public:
  RQMDsBox(Pythia8::Settings* s, InitialCondition* ini);
  ~RQMDsBox();
  void evolution(std::list<EventParticle*>& plist,double t,double dt,int step);
  void init(std::list<EventParticle*>& plist);
  void singleParticlePotential();

  void initBox(Pythia8::Settings* s);
  void addParticle(std::list<EventParticle*>& plist,double dt, int step);
  //void initMatrix(std::list<EventParticle*>& plist,double dt, int step);
  //void deleteMatrix();

  void RecoverEnergy(std::list<EventParticle*>& plist);
  void RecoverEnergy2(std::list<EventParticle*>& plist);
  double funcEnergy(std::list<EventParticle*>& plist,double a);

  //double facSDensity(Vec4& p, double m);
  Pythia8::Vec4 computeEnergy(std::list<EventParticle*>& plist, int step);
  /*
  double facScalar(Vec4& p, double m) {
    if(optScalarDensity==0) return p.mCalc()/p[0];
    else return m/sqrt(m*m + p.pAbs2());
  }
  */

}; 
}  // namespace jam2
#endif
