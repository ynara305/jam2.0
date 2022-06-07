#ifndef jam2_hadrons_DecayWidth_h
#define jam2_hadrons_DecayWidth_h

#include <cmath>
#include <Pythia8/Basics.h>
#include <Pythia8/ParticleData.h>

namespace jam2 {

class DecayWidth
{
private:
  int optDeltaWidth, optWidth;
  int itag;
  double minKinE;
  double pWidth[70];
  int nChannel;
  Pythia8::ParticleData*   particleData;
  bool useTable;
public:
  DecayWidth(Pythia8::ParticleData* pd);
  double getTotalWidth(Pythia8::ParticleDataEntry* p, double emcm,
	    int kf1=0, int kf2=0);
  double getPartialWidth(Pythia8::ParticleDataEntry* pd,double emcm, int kf1, int kf2);
  double BlattWeisskopf(double x,int l);
  double formBlattWeisskopf(double m0, double m,double pf0, double pf, int l,int kf1,int kf2);
  double BlattWeisskopfInt(double m, double m0,int kfd1,int kfd2, int l,double R);
  double getPWidth(int i) const {return pWidth[i];}
  int getIP() const {return itag;}
  int getNChannel() {return nChannel;}
  double deltaWidth(Pythia8::ParticleDataEntry* pd, int id0, double emcm,
	    int kf1=0, int kf2=0);
  double getDeltaWidth(int iwidth, double emd);
  double lifeTime(Pythia8::ParticleDataEntry* pd,double m, double e);
  double BWMass(double emmin,double emmax,double emr,double wid);

};
}
#endif

