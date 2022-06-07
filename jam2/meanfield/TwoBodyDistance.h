#ifndef jam2_meanfield_TwoBodyDistance_h
#define jam2_meanfield_TwoBodyDistance_h

#include <Pythia8/Basics.h>
#include <Pythia8/Settings.h>
#include <jam2/collision/EventParticle.h>

namespace jam2 {

class TwoBodyDistance
{
protected:
  Pythia8::Settings* settings;
  int transportModel,optP0dev,optPV,optScalarDensity;
  double gamCM,gamCM2;
  double widG, wG, facG;

  Pythia8::Vec4 r1,r2,p1,p2;
  Pythia8::Vec4 v1, v2;
  double em1,em2;

public:
  TwoBodyDistance(Pythia8::Settings* s);
  virtual ~TwoBodyDistance() {};
  virtual void density()=0;
  virtual void psq()=0;
  virtual void distanceR()=0;
  virtual void distanceP()=0;
  virtual void devV(const Pythia8::Vec4& Ai, const Pythia8::Vec4& Aj)=0;
  virtual void setRP1(Pythia8::Vec4& r,Pythia8::Vec4& p) {
    r1 = r; p1 = p; v1=p1/p1[0]; em1=p1.mCalc();
  }
  virtual void setRP2(Pythia8::Vec4& r,Pythia8::Vec4& p) {
    r2 = r; p2 = p;v2=p2/p2[0]; em2=p2.mCalc();
  }
  virtual void devGamma()=0;
  Vec4 V1() const {return v1;}
  Vec4 V2() const {return v2;}
  Pythia8::Vec4 dp2ijpi, dp2jipi, dp2ijpj, dp2jipj;
  Pythia8::Vec4 dr2ijri, dr2jiri, dr2jirj, dr2ijrj;
  Pythia8::Vec4 dr2ijpi, dr2jipi, dr2jipj, dr2ijpj;
  Pythia8::Vec4 devgam1,devgam2;
  Pythia8::Vec4 devV1,devV2;
  double density1, density2;
  double psq1, psq2;

};

class NonRelDistance : public TwoBodyDistance
{
public:
  NonRelDistance(Pythia8::Settings* s) :TwoBodyDistance(s) {};
  void density();
  void psq();
  void distanceR();
  void distanceP();
  void devGamma() {;}
  void devV(const Pythia8::Vec4& Ai, const Pythia8::Vec4& Aj);
};

class TwoBodyCM : public TwoBodyDistance
{
public:
  TwoBodyCM(Pythia8::Settings* s) :TwoBodyDistance(s) {};
  void density();
  void psq();
  void distanceR();
  void distanceP();
  void devGamma();
  void devV(const Pythia8::Vec4& Ai, const Pythia8::Vec4& Aj);
};

class RestFrame : public TwoBodyDistance
{
public:
  RestFrame(Pythia8::Settings* s) :TwoBodyDistance(s) {};
  void density();
  void psq();
  void distanceR();
  void distanceP();
  void devGamma() {;}
  void setRP1(Pythia8::Vec4& r,Pythia8::Vec4& p) {
    r1 = r; p1 = p; 
    em1 = p1.mCalc();
    v1=p1/em1;
  }
  void setRP2(Pythia8::Vec4& r,Pythia8::Vec4& p) {
    r2 = r; p2 = p;
    em2 = p2.mCalc();
    v2=p2/em2;
  }
  void devV(const Pythia8::Vec4& Ai, const Pythia8::Vec4& Aj);
};

}
#endif
