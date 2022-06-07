#ifndef jam2_interaction_BeamHadron_h
#define jam2_interaction_BeamHadron_h

#include <Pythia8/Basics.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/BeamParticle.h>

namespace jam2 {

using Pythia8::ParticleData;
using Pythia8::Settings;
using Pythia8::Rndm;
using Pythia8::Vec4;
using Pythia8::ResolvedParton;

class ResolvedPartonQ : public ResolvedParton {
private:
    int constQ;
public:
    ResolvedPartonQ(int icq=0, int ipos=0, int id=0, double x=0.,
	    int coman=-1) : ResolvedParton(ipos,id,x,coman) {
	constQ=icq;
    }
    void setConstQ(int i) {constQ=i;}
    int getConstQ() {return constQ;}
};


class BeamHadron {
private:
    int  idBeam;
    Vec4 pBeam;
    double mBeam;
    double eCM;
    int   pickG;
    int    nValKinds, idVal[3], nVal[3];
    std::vector<ResolvedPartonQ> resolved;
    bool isBaryonBeam,isHadronBeam,isMesonBeam;
    double valencePowerMeson, valencePowerUinP, valencePowerDinP,
         valenceDiqE0, valenceDiqE,valenceDiqEnhance, eMinPert,eWidthPert;
    double companionPower, gluonPower,xGluonCutoff;
    ParticleData* particleData;
    Rndm* rndm;
    int nValenceq;

public:
  BeamHadron(Settings& settings, ParticleData* pd, Rndm* r);
  void init(double ecm,int id, int cq[2],Vec4 p, double m,int kf1, int kf2, int pickg=0);

  int id()            const {return idBeam;}
  Vec4 p()            const {return pBeam;}
  double px()         const {return pBeam.px();}
  double py()         const {return pBeam.py();}
  double pz()         const {return pBeam.pz();}
  double e()          const {return pBeam.e();}
  double m()          const {return mBeam;}
 // Overload index operator to access a resolved parton from the list.
  ResolvedPartonQ& operator[](int i) {return resolved[i];}
  const ResolvedPartonQ& operator[](int i) const {return resolved[i];}
  int size() const {return resolved.size();}
  bool isBaryon()     const {return isBaryonBeam;}
  std::vector<ResolvedPartonQ>& parton() {return resolved;}

  void initBeamKind();
  void newValenceContent(int idnew);

  // Pick unrescaled x of remnant parton (valence or sea).
  double xRemnant(int i);
  int nValence() const {return nValenceq;}

  int nValence(int idIn) const {for (int i = 0; i < nValKinds; ++i)
      if (idIn == idVal[i]) return nVal[i];
    return 0;}


};

} // end namespace jam2
#endif
