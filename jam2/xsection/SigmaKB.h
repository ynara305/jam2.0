#ifndef jam2_xsection_SigmaKB_h
#define jam2_xsection_SigmaKB_h

#include <jam2/xsection/SigmaMB.h>
#include <Pythia8/Basics.h>
#include <Pythia8/Settings.h>
#include <unordered_map>

namespace jam2 {

// Purpose: to treat Kaon-nonstrange baryon collisions.
class SigmaKB: public SigmaMB
{
protected:
    CollisionPair* cpair;
    //int id1, id2;
    int kfm, kfb;
    int izm, izb, izt;
    int kfm1,kflr1,kflb1,kflc1, kfls1;
    int ibra;
    double snew,srt,pr,plab;
    int istr;
    Pythia8::Rndm* rndm;
    std::unordered_map<int,int> deltaMto0, delta2toP;
public:
    SigmaKB(Pythia8::Settings* s,JamParticleData* table, SampleMass* sm,
	    Pythia8::Rndm* r);
    ~SigmaKB() { };
    void calcKN(CollisionPair& cpair);
    bool chargeExchange();
    void sigKNDelta();
    void sigKNKstar();
    void deltaAbsorption();
    void absorbKstar();
    void KstarD();
    void KstarDabsorb();
    double jamsigkn(int isig,double s,double plab);
};
}
#endif

