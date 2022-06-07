#ifndef jam2_interaction_Scatter2_h
#define jam2_interaction_Scatter2_h

#include <jam2/interaction/ScatterKin.h>
#include <vector>
#include <cmath>

namespace jam2 {

class Scatter2: public ScatterKin
{
protected:
  int optSuppressSoftCollision,iBack;
  double paramSoftColl;
public:
  Scatter2(Pythia8::Info* inf,Pythia8::Settings* s,JamParticleData* jd,
	    CrossSection* xs, Pythia8::Rndm* r);
  ~Scatter2();
  void absorb(std::vector<EventParticle*>& outgoing);
  void scatter2(TwoBodyInterList* inter,Collision* event,std::vector<EventParticle*>& outgoing);

  double uniformAngle() {return  2*rndm->flat()-1;}
  bool checkSoftCollisionTime(double e1, double e2, double pz1, double pz2,
       TwoBodyInterList* inter,Collision* event);

  // Functions: momentum in two-particle cm.
  static double pawt(const double a,const double b,const double c) {
	return sqrt((a*a-(b+c)*(b+c))*(a*a-(b-c)*(b-c)))/(2.0*a);}

  double getThat(const double srt, const double em1, const double em2,
                 const double em3, const double em4, const double debyemass);

  void print(ostream& ofs,std::vector<EventParticle*>& outgoing) const;
  double elasticAngle(double pr);
  double elasticSlope(int kf1,int kf2,int ibar1,int ibar2,double plab,double snew);
  double inelasticAngle(double pini,double pr);
  double angInel(double srt, double pr0, double pr);
  double angDeltaN();

};
}
#endif

