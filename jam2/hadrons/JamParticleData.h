#ifndef jam2_hadrons_JamParticleData_h
#define jam2_hadrons_JamParticleData_h

#include <Pythia8/ParticleData.h>
#include <jam2/hadrons/ParticleTable.h>
#include <jam2/hadrons/DecayWidth.h>
#include <jam2/hadrons/DecayWidthTable.h>
#include <unordered_map>
//#include <map>

namespace jam2 {

class JamParticleData
{
private:
    Pythia8::ParticleData* table;  // particle data table for all particles.
    ParticleTable *nucleons,*deltas,*lambdas,*sigmas,*xis;
    ParticleTable *mesons;
    ParticleTable *nStar, *pStar, *dmStar,*d0Star,*dpStar,*dppStar;
    ParticleTable *smStar, *s0Star, *spStar,*x0Star,*xmStar;
    ParticleTable *light0Meson;
    ParticleTable *light1Meson0, *light1Mesonp;
    ParticleTable *strMeson0, *strMesonp;
    ParticleTable *smHadrons;
    std::unordered_map<int,int> pID;
    //std::map<std::int,std::int> pID;
    Pythia8::ParticleDataEntry *proton, *neutron;
    Pythia8::ParticleDataEntry *deltam,*delta0,*deltap,*deltapp;
    Pythia8::ParticleDataEntry *n1440, *p1440;
    std::unordered_map<int,int> nstarID;
    std::unordered_map<int,int> dstarID;
    std::unordered_map<ParticleDataEntry*,DecayWidthTable*> decayWidthTable;
    DecayWidth* decayWidth;
    Pythia8::Rndm*           rndm;
    std::ofstream out;

public:
    JamParticleData(Pythia8::ParticleData* pd, Pythia8::Rndm* r);
    ~JamParticleData();

    Pythia8::ParticleData* getParticleData() {return table;}
    Pythia8::ParticleDataEntry* find(int kf) {return table->findParticle(kf);}
    Pythia8::ParticleDataEntry* getProton() {return proton;}
    Pythia8::ParticleDataEntry* getNeutron() {return neutron;}
    Pythia8::ParticleDataEntry* getDeltam() {return deltam;}
    Pythia8::ParticleDataEntry* getDelta0() {return delta0;}
    Pythia8::ParticleDataEntry* getDeltap() {return deltap;}
    Pythia8::ParticleDataEntry* getDeltapp() {return deltapp;}
    Pythia8::ParticleDataEntry* getN1440() {return n1440;}
    Pythia8::ParticleDataEntry* getP1440() {return p1440;}
    ParticleTable* getNucl()  {return nucleons;}
    ParticleTable* getDelta()  {return deltas;}
    ParticleTable* getLambda()  {return lambdas;}
    ParticleTable* getSigma()  {return sigmas;}
    ParticleTable* getXi()     {return xis;}
    ParticleTable* getMeson()  {return mesons;}
    ParticleTable* getNstar()  {return nStar;}
    ParticleTable* getPstar()  {return pStar;}
    ParticleTable* getDmstar()  {return dmStar;}
    ParticleTable* getD0star()  {return d0Star;}
    ParticleTable* getDpstar()  {return dpStar;}
    ParticleTable* getDppstar()  {return dppStar;}
    ParticleTable* getSpstar()  {return spStar;}
    ParticleTable* getS0star()  {return s0Star;}
    ParticleTable* getSmstar()  {return smStar;}
    ParticleTable* getX0star()  {return x0Star;}
    ParticleTable* getXmstar()  {return xmStar;}
    ParticleTable* getLight0Meson()  {return light0Meson;}
    ParticleTable* getLight1Meson0()  {return light1Meson0;}
    ParticleTable* getLight1Mesonp()  {return light1Mesonp;}
    ParticleTable* getStrMeson0()  {return strMeson0;}
    ParticleTable* getStrMesonp()  {return strMesonp;}
    ParticleTable* getSMHadrons()  {return smHadrons;}
    void setPidMap();
    void addDiffractiveParticle(Pythia8::ParticleData* pd);
    int nStarID(int id) {return nstarID[id];}
    int dStarID(int id) {return dstarID[id];}
    int pid(int id) {
	std::unordered_map<int,int>::const_iterator it=pID.find(std::abs(id));
	//auto it=pID.find(id);
	if(it==pID.end() ) {
	    std::cout << "JamParticleData::pid particle not found id= " <<
		id << std::endl;
	    exit(1);
	} else {
	    return it->second;
	}
    }
  double totalWidth(Pythia8::ParticleDataEntry* pd, double m) {

    if(m < 6.0) {
      return decayWidthTable[pd]->getWidth(m);
    } else {
	std::cout << "JamParticleData::totalWidth mass exceed table limit "
	          << m << " id= "<< pd->id() << std::endl;
      return decayWidth->getTotalWidth(pd,m);
    }
  }
  DecayWidth* getDecayWidth() {return decayWidth;}
  double lifeTime(Pythia8::ParticleDataEntry* pd,double m, double e);
  double BWMass(double emmin,double emmax,double emr,double wid);
};
}
#endif
