#ifndef jam2_hadrons_Baryons_h
#define jam2_hadrons_Baryons_h

#include <Pythia8/ParticleData.h>
#include <jam2/hadrons/ParticleTable.h>

namespace jam2 {

using Pythia8::ParticleData;

class Baryons
{
public:
    //void setNucleons(ParticleData* table);
    //void setLambdas(ParticleData* table);
    //void setSigmas(ParticleData* table);
    //void setXis(ParticleData* table);
    //void addBaryons(Pythia8::ParticleData* table);

    void addNuclResonance(int i,ParticleData* table,ParticleTable* n,
	    ParticleTable* nstar,ParticleTable* pstar);
    void addDeltaResonance(int i,ParticleData* table,ParticleTable* d,
	    ParticleTable* dmstar,ParticleTable* d0star,
	    ParticleTable* dpstar,ParticleTable* dppstar);
    void addLambdaResonance(int i,ParticleData* table,ParticleTable* l);
    void addSigmaResonance(int i,ParticleData* table,ParticleTable* s,
	ParticleTable* smstar,ParticleTable* s0star,ParticleTable* spstar);
    void addXiResonance(int i,ParticleData* table,ParticleTable* x,
	    ParticleTable* xmstar,ParticleTable* x0star);

    void setNuclResonance(ParticleData* table,ParticleTable* n,
	    ParticleTable* nstar,ParticleTable* pstar);
    void setDeltaResonance(ParticleData* table,ParticleTable* d,
	    ParticleTable* dmstar,ParticleTable* d0star,
	    ParticleTable* dpstar,ParticleTable* dppstar);
    void setLambdaResonance(ParticleData* table,ParticleTable* l);
    void setSigmaResonance(ParticleData* table,ParticleTable* s,
	ParticleTable* smstar,ParticleTable* s0star,ParticleTable* spstar);
    void setXiResonance(ParticleData* table,ParticleTable* x,
	    ParticleTable* xmstar,ParticleTable* x0star);
};
}
#endif // jam2_hadrons_Baryons_h
