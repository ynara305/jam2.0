#ifndef jam2_initcond_BoostedTwoNuclei_h
#define jam2_initcond_BoostedTwoNuclei_h

#include <cmath>
#include <string>
#include <jam2/collision/Collision.h>
#include <jam2/initcond/InitialCondition.h>
#include <jam2/hadrons/JamParticleData.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/HIUserHooks.h>

namespace jam2 {

class BoostedTwoNuclei: public InitialCondition
{
protected:
    Pythia8::NucleusModel* proj;
    Pythia8::NucleusModel* targ;
    vector<Pythia8::Nucleon> projNucl, targNucl;
    vector<Pythia8::Vec4> pProj, pTarg;
    std::string compFrame;
    double gWidth;
    double numberOfParticle;
    double mP, mN;
    Pythia8::ParticleDataEntry *paA, *paB, *paP, *paN;

    int optFermiMomentum;
    bool histWS;
    const int nhist=80;
    int   nEvent;
    double rMax, pMax, dR, dP, widR, widP, facR, facP;
    double *histr, *histp, *rhor, *rhop;

public:
    BoostedTwoNuclei();
    BoostedTwoNuclei(Pythia8::Settings* s, JamParticleData* pd,
	    Pythia8::Rndm* r);
    virtual ~BoostedTwoNuclei() {
	if(histWS) printHist();
	delete proj;
	delete targ;
        delete [] histr;
        delete [] histp;
    }

    void init();
    void generate(Collision* event,int mode=0);
    vector<Pythia8::Vec4> generateFermiMomentum(vector<Pythia8::Nucleon>& nucl);

    double getNumberOfParticle() const {return numberOfParticle;}
    int    getAproj() const {return proj->A();}
    int    getZproj() const {return proj->Z();}
    int    getAtarg() const {return targ->A();}
    int    getZtarg() const {return targ->Z();}
    double getRTarg() const {return targ->R();}
    double getRProj() const {return proj->R();}
    //double getRadius(int a) { return 1.19*std::pow(a,0.333333) - 1.61*std::pow(a,-0.333333);}

    void initHist();
    void fill(vector<Pythia8::Nucleon>& nucl, std::vector<Pythia8::Vec4> pf);
    void printHist();
};

}

#endif
