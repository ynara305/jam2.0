#ifndef jam2_xsection_SigmaBB_h
#define jam2_xsection_SigmaBB_h

#include <jam2/xsection/CollisionPair.h>
#include <Pythia8/Basics.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/SigmaTotal.h>
#include <jam2/collision/EventParticle.h>
#include <jam2/interaction/HadronDecay.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/xsection/SampleMass.h>

namespace jam2 {

using Pythia8::ParticleDataEntry;

class SigmaBB
{
protected:
    int    mMode;
    double sigTot,sigEl;
    //double srt,pr;
    //int    iz1,iz2;
    int    kf1,kf2,id1,id2;
    double em1,em2, em3,em4;
    EventParticle *pa1,*pa2;
    ParticleDataEntry *pd1,*pd2, *pout[2];
    SampleMass* sampleMass;
    Pythia8::Info *info;
    Pythia8::Settings *settings;
    Pythia8::Rndm* rndm;
    Pythia8::SigmaTotal *sigTotPy8;
    //HadronDecay* decay;
    DecayWidth* decay;
    JamParticleData* jamTable;
    ParticleTable *nStar[2], *dStar[4];
    std::vector<ParticleDataEntry*> NStar[2], DStar[4];
    std::vector<ParticleDataEntry*> protonT, neutronT;
    std::vector<ParticleDataEntry*> deltamT, delta0T, deltapT, deltappT;
    std::vector<ParticleDataEntry*> N1440T, P1440T;
    double sigin[30];
    int mChanel;
    int nMaxRes;
    bool absorptionBB;
    static const double maxRMass;
    static const double Mp,Mn,Mnucl,Mdelta,Mpion,sMin;
    static const double normD1232;
    static const int  mxhlf;
    static const int  BBnn[16][2], BBpp[16][2],BBpn[17][2];
    static const int  DD3[9][2], DD4[9][2],DDm2[9][2],DDm1[9][2];
    static const int NDm1[9][2], ND3[9][2];

public:
// 1:n* 2:p* 3:D*- 4:D*0 5:D*+ 6:D*++
    enum NDchannel {Nstar=1,Pstar=2, Dstarm=3, Dstar0=4, Dstarp=5, Dstarpp=6};
    //enum NDchannel {Nstar,Pstar, Dstarm, Dstar0, Dstarp, Dstarpp};

    Pythia8::ParticleDataEntry* outParticle(int i) {return pout[i];}

    SigmaBB(Pythia8::Info *inf,Pythia8::Settings *s,JamParticleData* jd, SampleMass *sm, Pythia8::Rndm* r);
    virtual ~SigmaBB() {delete sigTotPy8;}
    double sigmaTot() const {return sigTot;}
    double sigmaEl()  const {return sigEl;}
    virtual void calc(CollisionPair& cpair)=0;
    virtual bool sampleResonanceMass(double em,int idn1, int idn2, double& m1, double& m2)=0;

    void calcS(CollisionPair& cpair);
    void jamxnn(double ecm,int z1, int z2);
    void   makePP(double* sig1);
    void   makePN(double* sig1);
    double sigmaS1(int ipair,double ecm,int z1,int z2);
    double sigmaS2(int ipair,double ecm,int z1,int z2);
    std::vector<ParticleDataEntry*> pickRTable(int idn, int& idr);
    bool sampleMassFix(double ecm,double pf0,int iex[2],double emr[2],
	    double em0[2], double gam0[2],double mmin[2],
	    double ymin[2],double ymax[2],double& m1, double& m2);

};
}
#endif

