#ifndef jam2_xsection_SigmaBB2_h
#define jam2_xsection_SigmaBB2_h

#include <jam2/xsection/SigmaBB.h>
#include <jam2/xsection/SigmaRR.h>

namespace jam2 {

using Pythia8::ParticleDataEntry;

class SigmaBB2 : public SigmaBB
{
private:
    double srt,pr;
    int    iz1,iz2;
    int optRR;
    SigmaRR *sigmaRR;
    int nError=0;
public:
    SigmaBB2(Pythia8::Info *inf,Pythia8::Settings *s,JamParticleData* jd, SampleMass *sm,
	    Pythia8::Rndm* r,std::string fname);
    ~SigmaBB2();
    void calc(CollisionPair& cpair);
    bool sampleResonanceMass(double em,int idn1, int idn2,
	                     double& m1, double& m2);
    void sigmaNuclNucl(CollisionPair& cp);
    void sigmaPP(CollisionPair& cp);
    void sigmaPN(CollisionPair& cp);
    int  sigmaNsNs(CollisionPair& cp);
    int  sigmaNsDs(CollisionPair& cp);
    int  sigmaDD(CollisionPair& cp);
};
}
#endif

