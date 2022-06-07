#ifndef jam2_interaction_SoftString_h
#define jam2_interaction_SoftString_h

#include <jam2/interaction/Scatter2.h>
#include <Pythia8/Pythia.h>
#include <Pythia8/Basics.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/SigmaTotal.h>
#include <jam2/collision/EventParticle.h>
#include <jam2/xsection/HadronContent.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/interaction/BeamHadron.h>

namespace jam2 {

using namespace Pythia8;

class SoftStrings : public Scatter2
{
private:
    double primordialKTremnant;
    bool isBaryonBeam;
    double valencePowerMeson, valencePowerUinP, valencePowerDinP;
    double companionPower, gluonPower,xGluonCutoff;
    double valenceDiqEnhance, eMinPert, eWidthPert, valenceDiqE;
    int primordialKTsChannelString,optConstQuarkDiffra;
    double probDiffra;
    double probBBAnn1, probBBAnn2, probBBAnn3;
    double probNonDiff, probSDA, probSDB;
    double probFlavorExchange;
    double mA, mB;

    int optConstQscatt,noBaBAnnihilation;
    bool optConstFormationTime;
    bool doShower;
    double aveConstFormationTime, aveFormationTime;
    int nCFormationTime, nFormationTime;

    double eCM,pCM;
    Pythia8::RotBstMatrix MfromCM;
    int idBeamA, idBeamB;
    int cqBeamA[2], cqBeamB[2];
    BeamHadron *beam1, *beam2;
    std::vector<ResolvedPartonQ> partonA, partonB;

    ParticleData* particleData;
    HadronContent* hadronContent;
    Pythia* hadronize;
    SigmaTotal sigtot;
    StringFlav* flavSel;
    std::ofstream osOut;
    int nError=0;

  // pt for diffractive process and BBar annihilation.
  double sigmaQ, sigmaQBBar, sigmaAnn;
  // Properties of the current collision. 1 or 2 is two incoming hadrons.
  // "c" or "ac" is colour or anticolour component of hadron.
  int    idc1, idac1, idc2, idac2;
  double z1, z2, mT1, mT2,
         mc1, mac1, px1, py1, pTs1, mTsc1, mTsac1, mTc1, mTac1,
         mc2, mac2, px2, py2, pTs2, mTsc2, mTsac2, mTc2, mTac2;
  double fracEtass, fracEtaPss, xPowMes, xPowBarU, xPowBarD,xPowBarS,xDiqEnhance;
  int passStringFrag;
  int allowRescatterSameString;

public:
    SoftStrings(Info* inf, Settings& settings, JamParticleData* pd,
	CrossSection* xs, Pythia* py, StringFlav* flav, Rndm* r);
    ~SoftStrings();
    bool setKinematics(int type,std::vector<EventParticle*>& outgoing);

  // Compute diffreactive and non-diffractive cross sections.
  void calcSigma();

  // Handle inelastic nondiffractive collision.
  bool nondiff();
  // Choose relative momentum of colour and anticolour constituents in hadron.
  double splitZ(int iq1, int iq2, double mRat1, double mRat2);
  // Estimate lowest possible mass state for flavour combination.
  double mThreshold( int iq1, int iq2);

  // Handle elastic and diffractive collisions.
  bool eldiff( int type, Vec4& pA, Vec4& pB);

  // Split up hadron A or B into a colour pair, with masses and pT values.
  bool splitA( double redMpT);
  bool splitB( double redMpT);

  // Split a hadron into a colour and an anticolour part.
  pair< int, int> splitFlav( int id);

  // Pick slope b of exp(b * t) for elastic and diffractive events.
  double bSlope( int type);

  EventParticle* putHadron(int id1,double m,double sm,Vec4& xout,
   double tform,Vec4& p3, int cq1[2],double spot, Pythia8::Vec4& vpot);

  bool nonDiff();

  bool hijingSoft(int type,std::vector<EventParticle*>& outgoing);
  bool partonMomentum(int sgn, int id,
	    int kfa1,int kfa2,double qm1, double qm2,
	    Vec4& p4, std::vector<ResolvedPartonQ>& q,
	    int cq[2]);

  bool absorbS(int id1, int id2,int id3, const int cq1[2], const int cq2[2],std::vector<EventParticle*>& outgoing);

  bool BaBAnnihilation(std::vector<EventParticle*>& outgoing);

  int findStringBBar(int id1, int id2,int ianti1,int ianti2,
	double ecm,
	int kfl1[3],int kfl2[3],int& kfm, double& emm);

  bool BBarQuarkMom(int ns, double pe, double pz, double emm, int jt,
	ResolvedPartonQ q[3]);

  bool stringFragment(int side,std::vector<ResolvedPartonQ>& q,
	Vec4& xout, double sm, double spot, Pythia8::Vec4& vpot,std::vector<EventParticle*>& outgoing);

  void makeHadron(int id, double sm,std::vector<ResolvedPartonQ>& q,
	Vec4& xout, double tform,double spot, Pythia8::Vec4& vpot,std::vector<EventParticle*>& outgoing);

  void makeString(int nval, std::vector<ResolvedPartonQ>& q,std::vector<pair<int,int>>& constQ, Event& event);

  void findConstQuark(Event& event,int ibar,int nconstq, 
	    std::vector<std::pair<int,int>>& q);

  void findConstQuark2(Event& event,int ibar,int nconstq, 
	    std::vector<ResolvedPartonQ>& q, int* iend,int* constq);

  void recoverEnergy(vector<EventParticle*>& hadrons, Vec4& ptot,double sm, double spot, Vec4& vpot);

    double xSample1(double srt,double xmin=0.0, double xmax=1.0);
    double xSample3(double srt,double xmin=0.0, double xmax=1.0);
    double constFormatinotime() {
	return aveConstFormationTime/std::max(1,nCFormationTime);}
    double formatinotime() {return aveFormationTime/std::max(1,nFormationTime);}
    double mJet(std::vector<ResolvedPartonQ>& parton);
    int findIdRes(std::vector<ResolvedPartonQ>& parton,double m);
    int isResonance(double m, double sm,std::vector<ResolvedPartonQ>& parton);
    double HadronMass(int id1, int id2);
    double GSHadronMass(int id1);

};
} // end namespace jam2
#endif
