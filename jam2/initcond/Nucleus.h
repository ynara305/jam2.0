#ifndef jam2_initcond_Nucleus_h
#define jam2_initcond_Nucleus_h

#include <cmath>
#include <Pythia8/HIUserHooks.h>
#include <Pythia8/Basics.h>

namespace jam2 {

class Nucleus : public Pythia8::GLISSANDOModel
{
protected:
    //vector<Vec4> momentum;
    std::vector<Pythia8::Nucleon> generate() const;
    bool gaussHardCore;
    int  optSample;
    double R0,rHard;

    // working area for deuteron.
    double x0,y0,z0;

public:
    Nucleus() : gaussHardCore(false) {
      R0=0.03;
    }
    virtual ~Nucleus() {
    }
    void setParam() {
      optSample = settingsPtr->mode("HeavyIon:optSample");
      R0 = settingsPtr->parm("HeavyIon:WSR0");
      rHard = 1.123*pow(double(A()),1.0/3.0) - R0;
    }
    Pythia8::Vec4 generateHardSphare() const;
    Pythia8::Vec4 generateDeuteron() const;
    void initHist();
};

}

#endif // jam2_initcond_Nucleus_h

