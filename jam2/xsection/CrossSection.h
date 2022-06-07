#ifndef jam2_xsection_CrossSection_h
#define jam2_xsection_CrossSection_h

#include <algorithm>
#include <jam2/xsection/SigmaBB.h>
#include <jam2/xsection/SigmaMB.h>
#include <jam2/xsection/SigmaKB.h>
#include <jam2/xsection/SigmaABB.h>
#include <jam2/xsection/SigmaMM.h>
//#include <jam2/xsection/JamParticleTable.h>
#include <jam2/collision/EventParticle.h>
#include <jam2/xsection/CollisionPair.h>
#include <Pythia8/Basics.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/xsection/SampleMass.h>


namespace jam2 {

class CrossSection
{
private:
    JamParticleData* jamParticleData;
    //CollisionPair*  cpair;
    double sigTot,sigEl;
    double eCM,pCM;
    int id1,id2;
    EventParticle *pa1,*pa2;
    Pythia8::ParticleDataEntry *pd1, *pd2;
    bool anti;
    int mSel;
    int optJAM1;
    double sigRand;
    SigmaBB*   sigBB;
    SigmaABB* sigABB;
    SigmaMB*  sigMB;
    SigmaKB*  sigKB;
    SigmaMM*  sigMM;
    SampleMass* sampleMass;
    double ecmStringBB,ecmStringMBc;
    Pythia8::Info *info;
    Pythia8::Settings* settings;
    Pythia8::Rndm* rndm;

    double SigmaBaryonBaryon(CollisionPair& cpair);
    double SigmaMesonBaryon(CollisionPair& cpair);
    double SigmaMesonMeson(CollisionPair& cpair);
    double SigmaAntiBaryonBaryon(CollisionPair& cpair);
    //double SigmaMBHeavy();
public:
    CrossSection(Pythia8::Info *inf,Pythia8::Settings* s, JamParticleData* pd, Pythia8::StringFlav* flav,Pythia8::Rndm* r);
    ~CrossSection();
    double sigmaTot() const {return sigTot;}
    double sigmaEl()  const {return sigEl;}
    double sigma(CollisionPair& cp);
    SigmaBB* getBB() {return sigBB;}
    SigmaMB* getMB() {return sigMB;}
    SigmaMM* getMM() {return sigMM;}
    SigmaKB* getKB() {return sigKB;}
    SigmaABB* getAB() {return sigABB;}
    void selectOutGoingParticle(CollisionPair& cpair);
    void sigmaTChannel(CollisionPair& cpair,double sigin);
  static int getHeavyQB(int kf) {
    int kfa=abs(kf);
    int nheavy=0;
    if((kfa/1000) %10 >=4) nheavy++;
    if((kfa/100)  %10 >=4) nheavy++;
    if((kfa/10)   %10 >=4) nheavy++;
    return nheavy;
  }
  static int getHeavyQM(int kf) {
    int kfa=abs(kf);
    int nheavy=0;
    if((kfa/100) %10 >=4) nheavy++;
    if((kfa/10)  %10 >=4) nheavy++;
    return nheavy;
  }


};
}
#endif

