#ifndef jam2_hadrons_Mesons_h
#define jam2_hadrons_Mesons_h

#include <string>
#include <Pythia8/ParticleData.h>
#include <jam2/hadrons/ParticleTable.h>

namespace jam2 {

using Pythia8::ParticleData;

class Mesons
{
public:
  void addMesons(ParticleData* table,ParticleTable* m,
	    ParticleTable* l0, ParticleTable* l10, ParticleTable* l1p,
	    ParticleTable* s0, ParticleTable* sp);
  void addLightIsoSpin0Meson(ParticleData* table,ParticleTable* m,
	  ParticleTable* l0);
  void addLightIsoSpin1Meson(ParticleData* table,ParticleTable* m,
	  ParticleTable* l10,ParticleTable* l1p);
  void addStrangeMeson(ParticleData* table,ParticleTable* m,
	  ParticleTable* str0,ParticleTable* strp);
};
}
#endif // jam2_hadrons_Baryons_h
