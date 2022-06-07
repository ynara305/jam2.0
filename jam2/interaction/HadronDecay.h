#ifndef jam2_interaction_HadronDecay_h
#define jam2_interaction_HadronDecay_h

#include <cmath>
#include <jam2/collision/EventParticle.h>
#include <jam2/collision/InterList.h>
#include <Pythia8/Pythia.h>
#include <Pythia8/Basics.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/ParticleDecays.h>
#include <Pythia8/Event.h>
#include <Pythia8/Info.h>
#include <Pythia8/Settings.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/interaction/JParticleDecays.h>

namespace jam2 {

class HadronDecay
{
private:
  int optDeltaWidth, optWidth, potentialHandling;
  int optDeltaDecayAngle,optDecayAngle, optDecayAngleS,optPosition,optPotential;
  int optVectorPotential,optVdot,optPropagate;
  int itag;
  double paramSmear;
  double minKinE;
  double pWidth[70];
  int nChannel;
  JamParticleData*         jamParticleData;
  Pythia8::ParticleData*   particleData;
  Pythia8::Pythia*         pythiaDecay;
  Pythia8::Event           pythiaEvent;
  Pythia8::Rndm*           rndm;
  Pythia8::Info*           info;
  Pythia8::Settings*       settings;
  JParticleDecays          jdecay;
  bool useTable;
public:
  HadronDecay();
  void init(Pythia8::Info *infor,Pythia8::Settings& s,JamParticleData* pd,
	    Pythia8::Pythia* d, Pythia8::Rndm* r);
  //void decay(InterList* inter, std::vector<EventParticle*>& outgoing);
  void decay(EventParticle* pa, std::vector<EventParticle*>& outgoing,
	    bool finalDecay=0);
  void findConstQuark(Event& event,int ibar,
      int* q, int* iend,int* constq);

  /*
  double deltaWidth(Pythia8::ParticleDataEntry* pd, int id0, double emcm,
	    int kf1=0, int kf2=0);
  double getTotalWidth(Pythia8::ParticleDataEntry* p, double emcm,
	    int kf1=0, int kf2=0);

  double BlattWeisskopf(double x,int l);
  double formBlattWeisskopf(double m0, double m,double pf0, double pf, int l,int kf1,int kf2);
  double BlattWeisskopfInt(double m, double m0,int kfd1,int kfd2, int l,double R);

  double getPWidth(int i) const {return pWidth[i];}
  int getIP() const {return itag;}
  int getNChannel() {return nChannel;}
  double getDeltaWidth(int iwidth, double emd);
  double lifeTime(Pythia8::ParticleDataEntry* pd,double m, double e);

  // Functions: momentum in two-particle cm.
  double  pawt(double a,double b,double c) {
	return sqrt((a*a-(b+c)*(b+c))*(a*a-(b-c)*(b-c)))/(2.0*a);}
  double BWMass(double emmin,double emmax,double emr,double wid);
  */

  /*
  static bool isDelta(int id) {
	if(id==1114 || id==2114 || id==2214 || id==2224)
	return true; else return false;}
	*/

};
}
#endif

