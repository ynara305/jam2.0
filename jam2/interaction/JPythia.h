#ifndef jam2_interaction_JPythia_h
#define jam2_interaction_JPythia_h

#include <Pythia8/Pythia.h>
#include <Pythia8/Basics.h>
#include <Pythia8/ParticleData.h>
#include <jam2/collision/EventParticle.h>
#include <jam2/interaction/PDFHadron.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/interaction/HadronDecay.h>
#include <string>

namespace jam2 {

using Pythia8::Vec4;

class JBeamShape : public BeamShape
{
private:
    double pzIni;
    Vec4 pA, pB;
public:
    JBeamShape(double pz) : pzIni(pz) {}
    void set(Vec4& p1, Vec4& p2) { pA=p1; pB=p2;}
    void set(double pz1, double  pz2) {
	pA=Vec4(0.0,0.0,pz1); pB=Vec4(0.0,0.0,pz2);}
    virtual void pick() {
	deltaPxA = pA.px();
	deltaPyA = pA.py();
	deltaPzA = pA.pz() - pzIni;
	deltaPxB = pB.px();
	deltaPyB = pB.py();
	deltaPzB = pB.pz() + pzIni;
    }
};

class JPythia
{
public:
  Settings *settings;
  // Read in one update for a setting or particle data from a single line.
  bool readString(std::string line, bool warn = true) {
      return pythia->readString(line,warn); }

private:
  Pythia8::Pythia *pythia;
  Pythia8::ParticleData *particleData;
  JBeamShape* beamShape;
  int isBaryonA, isBaryonB;
  std::vector<std::pair<int,int> > constQuark;
  std::vector<int> beamSide, isValenceQ;
  int optDiq;
  int isDebug;
  int procCode;
  double eCM, mA, mB, eA, eB;
  int allowRescatterSameString;

  // Private UserHooks class to select a specific process.
  struct ProcessSelectorHook: public UserHooks {

    ProcessSelectorHook(): proc(0), b(-1.0) {}

    // Yes we can veto event after process-level selection.
    virtual bool canVetoProcessLevel() {
      return true;
    }

    // Veto any unwanted process.
    virtual bool doVetoProcessLevel(Event&) {
	//cout << "### proc= "<< proc << " code= "<< infoPtr->code()
	//    << " this= "<< this <<endl;
      return proc > 0 && infoPtr->code() != proc;
    }

    // Can set the overall impact parameter for the MPI treatment.
    virtual bool canSetImpactParameter() const {
      return b >= 0.0;
    }

    // Set the overall impact parameter for the MPI treatment.
    virtual double doSetImpactParameter() {
      return b;
    }

    // The wanted process;
    int proc;

    // The selected b-value.
    double b;

  };

  // Holder class to temporarily select a specific process
  struct HoldProcess {

    // Set the given process for the given hook object.
    HoldProcess(ProcessSelectorHook & hook, int proc, double b = -1.0)
      : saveHook(&hook), saveProc(hook.proc), saveB(hook.b) {
      hook.proc = proc;
      hook.b = b;
    }

    // Reset the process of the hook object given in the constructor.
    ~HoldProcess() {
      if ( saveHook ) {
        saveHook->proc = saveProc;
        saveHook->b = saveB;
      }
    }

    // The hook object.
    ProcessSelectorHook * saveHook;

    // The previous process of the hook object.
    int saveProc;

    // The previous b-value of the hook object.
    double saveB;

  };

  // The process selector for the SASD object.
  ProcessSelectorHook selectSASD;


protected:
    PDFHadron* pdfHadron;
    Pythia8::PDF *pdfProton, *pdfPion;
    Pythia8::StringFlav* flavSel;
    Pythia8::RotBstMatrix MfromCM, MtoCM;
    JamParticleData* jamParticleData;
    HadronDecay* hDecay;
    bool optPreserveReactionPlane;
    bool preHadronA, preHadronB;
    bool userHook;

    int optConstQscatt,optConstQuarkDiffra;
    double aveConstFormationTime, aveFormationTime;
    int nCFormationTime, nFormationTime;

public:
    JPythia(std::string xmlDir = "../share/Pythia8/xmldoc", 
	    bool printBanner = true);
    JPythia(Settings& settingsIn, ParticleData& particleDataIn,
    bool printBanner = true);

    ~JPythia();
    void init();
    void init0(JamParticleData* jd,HadronDecay* dec);
    //bool initBeam(int id[2],double m[2], Vec4& pA, Vec4& pB);
    void initBeam(const int id[2], const double m[2],const double m0[2],Vec4& pA, Vec4& pB,
	    bool preHA, bool preHB);

    int whichSide(int ip, int& cq);

    bool generate(int ncoll, int operation, const int cq1[2], const int cq2[2],
	Vec4& xA, Vec4& xB, vector<EventParticle*>& outgoing);

    void transportHadron(Vec4& xA, Vec4& xB, int ncoll, int operation,
	vector<EventParticle*>& outgoing);

    int convertID(int id, int& isb);
    void findValenceQ(BeamParticle& beamP, const int cq1[2], int side);
    int findValenceQ71(int ip, int cq1,int cq2,int side);
    int findValenceQ74(int ip);

    void setBeamShapePtr2( JBeamShape* beamShapePtrIn)
    { beamShape = beamShapePtrIn;}

    double constFormatinotime() {
	return aveConstFormationTime/std::max(1,nCFormationTime);}
    double formatinotime() { return aveFormationTime/std::max(1,nFormationTime);}

    //auto getIdA() { return &Pythia8::Pythia::idASave;}

};
}
#endif

