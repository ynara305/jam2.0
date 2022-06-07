#ifndef jam2_interaction_Scatter_h
#define jam2_interaction_Scatter_h

#include <Pythia8/Basics.h>
#include <Pythia8/ParticleData.h>
#include <jam2/collision/EventParticle.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/collision/Collision.h>
#include <jam2/collision/TwoBodyInterList.h>
#include <jam2/xsection/CollisionPair.h>
#include <jam2/xsection/CrossSection.h>
#include <jam2/interaction/SoftStrings.h>
#include <jam2/interaction/JPythia.h>

#include <vector>
#include <cmath>

namespace jam2 {

using Pythia8::Vec4;

class Scatter
{
protected:
    Pythia8::Info*  info;
    Pythia8::Settings* settings;
    Pythia8::Rndm* rndm;
    JPythia  *pythia, *pythia2;
    SoftStrings* scatt;
    int channel;
    int softModel;
    int printColl;
    bool inelOnly;
    double eCMPythia,eMinPert,eWidthPert;

public:
  enum Channel{ELASTIC, RESONANCE, SOFT, HARD, DECAY, ABSORPTION,DIFFRACTIVE, BBANIHILATION};

  Scatter(Pythia8::Info* inf,Pythia8::Settings* s,JamParticleData* jd,
	    CrossSection* xs, Pythia8::Pythia* py, 
            Pythia8::StringFlav* flav,
	    Pythia8::Rndm* r);
  virtual ~Scatter();

  virtual void scatter(InterList* inter,std::vector<EventParticle*>& outgoing,Collision* evvent);

  void setPythia(JPythia* py) {pythia=py;}
  void setPythia2(JPythia* py) {pythia2=py;}
  int usePythia(double ecm);
  int getChannel() const {return scatt->getChannel();}

};
}
#endif

