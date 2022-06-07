#ifndef jam2_xsection_SigmaMM_h
#define jam2_xsection_SigmaMM_h

#include <jam2/xsection/CollisionPair.h>
#include <jam2/xsection/SampleMass.h>
#include <jam2/xsection/HadronContent.h>
#include <jam2/collision/EventParticle.h>
#include <Pythia8/Basics.h>
#include <Pythia8/ParticleData.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/hadrons/ParticleTable.h>

namespace jam2 {

class SigmaMM
{
private:
    int mSel;
    int idRes;
    double emRf;
    double sigRand;
    EventParticle *pa1,*pa2;
    Pythia8::Settings *settings;
    Pythia8::ParticleDataEntry *pd1,*pd2;
    Pythia8::ParticleDataEntry* pRes;
    ParticleTable* mesonTable;
    DecayWidth* decay;
    JamParticleData* jamtable;
    Pythia8::ParticleData* particleData;
    SampleMass* sampleMass;
    Pythia8::Rndm* rndm;
    int mChanel, mAbsrb;
    HadronContent* hadronContent;
    int optSChannel;
    double ecmStringMM;

public:
    SigmaMM(Pythia8::Settings *s,JamParticleData* table, SampleMass* sm,Pythia8::Rndm* r);
    ~SigmaMM() {
	delete hadronContent;
    }
    void calc(CollisionPair& cp);
    double jamxbw2(double srt,double pr,int kf1,int kf2,
	                       int iz1,int iz2,int isel);
    double BreitWigner(double srt, double pr,
	              int iz, int kf1, int kf2,int isel);

    //void attflv(int kf,int& ifla,int& iflb);
    //int ifrkfc(int ia,int ib,int ic,double s);
    //double jamemjet(int kfl10,int kfl20);
};
}
#endif

