#ifndef jam2_xsection_SigmaABB_h
#define jam2_xsection_SigmaABB_h

#include <jam2/xsection/CollisionPair.h>
//#include <jam2/interaction/HadronDecay.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/hadrons/ParticleTable.h>
#include <jam2/xsection/SampleMass.h>


namespace jam2 {

class SigmaABB
{
private:
    double sigTot,sigEl,sigAnn;
    double sigin[30];
    double sigRand;
    int mSel,mChanel,mAbsrb;
    JamParticleData* jamtable;
    //HadronDecay* decay;
    SampleMass* sampleMass;
    double ecmStringABB;
public:
    inline double pow2(const double& x) {return x*x;}
    SigmaABB(Pythia8::Settings* s, JamParticleData* table, SampleMass* sm) {
	jamtable=table;sampleMass=sm;
        ecmStringABB=s->parm("Cascade:ecmStringABB");
    }
    ~SigmaABB() {;}
    void calc(CollisionPair& cpair);
    void getXsection(double srt,double plab,int iz1,int iz2);
    double sigPPAnn(double plab);
};
}
#endif

