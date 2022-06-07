#ifndef jam2_xsection_SigmaBB1_h
#define jam2_xsection_SigmaBB1_h

#include <jam2/xsection/CollisionPair.h>
#include <Pythia8/Basics.h>
#include <Pythia8/ParticleData.h>
#include <jam2/collision/EventParticle.h>
#include <jam2/interaction/HadronDecay.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/xsection/SampleMass.h>
#include <jam2/xsection/SigmaNDelta.h>
#include <jam2/xsection/SigmaBB.h>

namespace jam2 {

using Pythia8::ParticleDataEntry;

class SigmaBB1 : public SigmaBB 
{
private:
  double srt,pr;
  int iz1,iz2;
    int optDeut1, optDeut2, optDetBal;
    int optSwave; // s-wave pion production is treated by N(1440) production
    int optProb, optDelta;
    double *bwR;
    SigmaNR *sigNDelta;
    static const double eps;
    static const int mxhlf;     // max. number of bin for integration

    //...Definition of Fitting Function.
    inline double  bw(double a,double b) {return b/((a*a-1)*(a*a-1)+b*b);}
public:
    SigmaBB1(Pythia8::Info *inf,Pythia8::Settings *s,JamParticleData* jd, SampleMass *sm, Pythia8::Rndm* r);
    ~SigmaBB1();
    void calc(CollisionPair& cpair);
    bool sampleResonanceMass(double em,int idn1, int idn2, double& m1, double& m2);
    void sigmaNuclNucl(CollisionPair& cp);
    void sigmaPP(CollisionPair& cp);
    void sigmaPN(CollisionPair& cp);
    int  sigmaNsNs(CollisionPair& cp);
    int  sigmaNsDs(CollisionPair& cp);
    int  sigmaDD(CollisionPair& cp);

    void jamxnnin(double x,double *sigin,int iporn,int ideut);
    double jambres2(ParticleTable* table, double mr, int id);
    double jamxdpi1(double x);
    double jamxdpi2(int iopt,double x);
    double DetBal1(int msel,Pythia8::ParticleDataEntry* p,double m1,double m2);
    double DetBal2(int msel);
    double DetBal3();
    double BW1(int idetsw,Pythia8::ParticleDataEntry* p,double em,double mr0,double m2);
    double BW2(double m1,double emin2);
    double BW3(double m1,double m2);
    Pythia8::ParticleDataEntry* selectDN(int idn, int& idr);
    Pythia8::ParticleDataEntry* selectResonance(ParticleTable* t);

};
}
#endif

