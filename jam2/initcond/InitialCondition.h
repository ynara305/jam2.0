#ifndef jam2_initcond_InitialCondition_h
#define jam2_initcond_InitialCondition_h

#include <cmath>
#include <fstream>
#include <Pythia8/Basics.h>
#include <Pythia8/Settings.h>
#include <Pythia8/ParticleData.h>
#include <jam2/hadrons/JamParticleData.h>

namespace jam2 {

using namespace Pythia8;

class Collision;

class InitialCondition
{
protected:
    JamParticleData *jamParticleData;
    Pythia8::ParticleData *particleData;
    Pythia8::Settings *settings;
    Pythia8::Rndm *rndm;
    std::ofstream ofs;
    int overSample;
    int  outputPhaseSpace;
    int masA, masB;
    double pzAcm, pzBcm, eA,eB, mA, mB;
    int nCollG, nPartG;
    double yCM, gamA, gamB;
    double impactPar,bMin,bMax;
    double eCM,eLab,pLab;
    double zSep;

public:
    InitialCondition(Settings* s, JamParticleData* pd, Rndm* r) {
	settings=s;
	jamParticleData = pd;
	particleData = jamParticleData->getParticleData();
	rndm=r;
	impactPar=0.0;
	overSample=settings->mode("Cascade:overSample");
	nCollG=nPartG=0;
    }
    virtual ~InitialCondition() { }
    virtual void init()=0;
    virtual void generate(Collision* event,int mode=0)=0;

    double ecm() const {return eCM;}
    double ycm() const {return yCM;}
    double gammaA() const {return gamA;}
    double gammaB() const {return gamB;}
    double massA() const {return mA;}
    double massB() const {return mB;}
    double energyA() const {return eA;}
    double energyB() const {return eB;}
    double pzA() const {return pzAcm;}
    double pzB() const {return pzBcm;}
    int    massNumberA()  const {return masA;}
    int    massNumberB()  const {return masB;}
    int    nColl()  const {return nCollG;}
    int    nPart()  const {return nPartG;}

    int getNumberOfOverSample() const {return overSample;}
    void  setNumberOfOverSample(int n) {overSample=n;}
    double getRadius(int a) const {
       	       return 1.19*std::pow(a,0.333333) - 1.61*std::pow(a,-0.333333);}
    void setOutputPhaseSpace(char* c) {
	 ofs.open(c);
	 if(!ofs) {
	     cerr << "cannot open file " << c << std::endl;
	     std::exit(1);
	 }
	 outputPhaseSpace=1;
    }
    static int findAZ(std::string nucl);
    void  setBMin(double b)  {bMin=b;}
    void  setBMax(double b)  {bMax=b;}
    double bmin()  const {return bMin;}
    double bmax()  const {return bMax;}
    double zseparation() const {return zSep;}
    void  setImpactPar(double b) {impactPar = b;}
    double  setImpactPar() { return
             sqrt(bMin*bMin + rndm->flat()*(bMax*bMax-bMin*bMin));}
    double getImpactPar()    const {return impactPar;}


};


class UserInitialCondition: public InitialCondition
{
public:
  UserInitialCondition(Pythia8::Settings* s, JamParticleData* pd,Pythia8::Rndm* r):
    InitialCondition(s,pd,r) {}
  virtual~UserInitialCondition() {}
  virtual void init() {}
  //virtual void generate(Collision* event,int mode=0){}
  virtual void generate(Collision*, int mode=0){}

};

} //end jam

#endif
