#ifndef jam2_xsection_SigmaMB_h
#define jam2_xsection_SigmaMB_h

#include <jam2/xsection/CollisionPair.h>
#include <jam2/interaction/HadronDecay.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/hadrons/ParticleTable.h>
#include <jam2/xsection/SampleMass.h>

namespace jam2 {

class SigmaMB
{
protected:
    Pythia8::Settings *settings;
    Pythia8::ParticleDataEntry *pd1,*pd2;
    double sigTot,sigEl;
    double sigin[30];
    double sigRand;
    int mSel;
    int mChanel,mAbsrb;
    Pythia8::ParticleDataEntry* pRes;
    DecayWidth* decay;
    JamParticleData* jamtable;
    SampleMass* sampleMass;
    int optSChannel, optJAM1;
    double sigElBW;
    double ecmStringMB,ecmStringMBs;

public:
    inline double pow2(const double& x) {return x*x;}
    inline double pow3(const double& x) {return x*x*x;}
    SigmaMB(Pythia8::Settings *s,JamParticleData* table, SampleMass* sm);
    ~SigmaMB() { };
    void calcS0(CollisionPair& cpair);
    void calcS1(CollisionPair& cpair);
    void calcS2(CollisionPair& cpair);
    void calcKN(CollisionPair& cpair);
    double piNtchannelElastic(double srt, int izt);
    double piNtchannelElasticJAM1(double srt, int izt);
    double jamxbw1(int istr,double srt,double pr,int iz1, int iz2,
	int& kf1,int& kf2,bool anti=false);


    double BreitWigner(ParticleTable* table,int i0,
	    double srt, double pr, int iz,
	    int& kf1, int& kf2,bool anti=false);

    void jamxpin(int izpi,int izn,double srt,double pr);
    double jamxdelt(int izpi,int izn,double srt,double pr);
    void jamxkp(int kfm,int kfb,double srt,double emmes,double embar,
	double& sig,double& sigel,double& sigch,double* sigy);

};
}
#endif

