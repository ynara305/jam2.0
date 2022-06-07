#ifndef Scatter_Kinh
#define Scatter_Kinh

#include <Pythia8/Basics.h>
#include <Pythia8/ParticleData.h>
#include <jam2/collision/EventParticle.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/collision/Collision.h>
#include <jam2/collision/TwoBodyInterList.h>
#include <jam2/xsection/CollisionPair.h>
#include <jam2/xsection/CrossSection.h>

//#include <jam2/interaction/pqcd/TimelikeShower.h>

#include <vector>
#include <cmath>

namespace jam2 {

using Pythia8::Vec4;

class ScatterKin
{
protected:
    Pythia8::Info*  info;
    Pythia8::Settings* settings;
    JamParticleData* jamParticleData;
    CrossSection* xSection;
    Pythia8::Rndm* rndm;

    EventParticle *inComing[2],*outGoing[2];
    Vec4 xCM;
    Vec4 pCM;
    double  eCM, sCM,pIni,pFinal;            // c.m. energy in GeV.
    double eCMcan, pInican,eCM1,eCM2;
    double tFormA,tFormB;
    Vec4 xIncoming[2];  //initial coordinate in the interaction point(cm).
    Vec4 pIncoming[2];  //  initial momenta.
    Vec4 pOut[4], xOut[4];
    bool preHadronA, preHadronB;
    int  *constQ[2];
    double qFactor[2];
    double  mIn[2],mInEff[2],mOut[4],mOutEff[4], SM[2];
    int nOperation, numCollision;
    double  phiCM, thetaCM,beX,beY,beZ;  // rotation angle and beta to c.m.
    Pythia8::RotBstMatrix MfromCM, MtoCM;
    Pythia8::RotBstMatrix MfromCMcan;
    int nOutParticle;
    Pythia8::ParticleDataEntry *pd[4];
    int     idIn[2];
    int     idOut[4],pidOut[4];
  bool   isBaryon[2];
  int    iBaryon[2];
    double sPot[2];
    Pythia8::Vec4 vPot[2];
    int channel,typeDiff;

    CollisionPair cpair;
    int printColl, optConstQuarkDiffra;
    bool inelOnly, optPreserveReactionPlane;
    int withMeanField,optVectorPotential,optVdot,optPropagate;
    int  optPotential;
    int isDebug;

public:
  enum Channel{ELASTIC, RESONANCE, SOFT, HARD, DECAY, ABSORPTION,DIFFRACTIVE, BBANIHILATION};

  ScatterKin(Pythia8::Info* inf,Pythia8::Settings* s,JamParticleData* jd,
	    CrossSection* xs, Pythia8::Rndm* r);
  virtual ~ScatterKin();

  double ecm() const {return eCM;}
  int idout(int i) const {return idOut[i];}
  const int* constq(int i) const {return constQ[i];}
  int nColl() const {return numCollision;}
  Vec4 p(int i) const {return pIncoming[i];}
  Vec4 pc(int i) const {return inComing[i]->getMomentum();}
  const int* idin() const {return idIn;}
  const double* meff() const {return mInEff;}
  const double* m() const {return mIn;}
  Vec4 xout(int i) const {return xOut[i];}
  bool prehadronA() const {return preHadronA;}
  bool prehadronB() const {return preHadronB;}
  CollisionPair& getCpair() {return cpair;}
  void checkMomentum(std::vector<EventParticle*>& outgoing);

  int typeDiffra() const {return typeDiff;}
  void setChannel(int i) {channel=i;}
  int getChannel() const {return channel;}
  void initialize(TwoBodyInterList* inter);
  int selectChannel(TwoBodyInterList* inter, int ncol);
  inline double PCM(double srt, double m1, double m2) {
        return sqrt((srt*srt-pow2(m1+m2))*(srt*srt-pow2(m1-m2)))/(2.0*srt);
  }

};
}
#endif

