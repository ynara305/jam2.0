#ifndef jam2_hadrons_JamStdlib_h
#define jam2_hadrons_JamStdlib_h
#include <cmath>
#include <Pythia8/Basics.h>

namespace jam2 {

//const double eKinMin=0.001; // parc(41)
const double eKinMin=0.0005;  // ParticleDecays:mSafty

enum particleID {id_quark=1,id_lept=2,id_exc=3,
    id_boson=4,id_diq=5,
    id_tec=6,id_susy=7,id_special=6,
//...Light Mesons.
    id_pi=101,id_light1=100, id_light0=120,
    id_str=130,id_charm=140,id_bott=150,id_cc=160,id_bb=170,id_mdiff=199,
    id_nucl=11, id_nucls=12, id_delt=13,id_delts=14,
    id_lambda=21,id_lambdas=22, id_sigma=31,id_sigmas=32, id_xi=41,id_xis=42,
    id_omega=51, id_charmb=61,id_bottb=72, id_bdiff=99};

// Functions: momentum in two-particle cm.
inline double PCM(double a,double b,double c) {
	return sqrt((a*a-(b+c)*(b+c))*(a*a-(b-c)*(b-c)))/(2.0*a);}

inline static double PCM2(double a,double b,double c) {
	return (a*a-(b+c)*(b+c))*(a*a-(b-c)*(b-c))/(4.0*a*a);}

// Function:lab. momentum.
inline double plabsr(double a,double b,double c) {
	return sqrt((a*a-b*b-c*c)*(a*a-b*b-c*c)/(4.0*c*c)-b*b);}

inline double pow2(const double& x) {return x*x;}

inline double modulo(double x, double y) { return x - std::floor(x/y)*y; }


inline double BW(double m0, double m, double gam) {
	return m*gam /( pow2(pow2(m)-pow2(m0)) + pow2(m*gam) );
    }

inline double cmDistanceSquare(const Pythia8::Vec4& r, const Pythia8::Vec4& u)
{
  //double dot4 = r.e()* u.e() - r.px()*u.px() - r.py()*u.py() - r.pz()*u.pz();
  double dot4 = r * u;
  return r.m2Calc() - dot4*dot4;
}

inline bool isDelta(int id) {
  return (abs(id)==1114 || abs(id)==2114 || abs(id)==2214 || abs(id)==2224)
      ? true:false;
}
inline bool isPion(int id) {
  return (id==111 || abs(id)==211) ? true:false;
}
inline bool isRho(int id) {
  return (id==113 || abs(id)==213) ? true:false;
}
inline bool isN1440(int id) {
  return (abs(id)==12112 || abs(id)==12212) ? true:false;
}
inline bool isN1520(int id) {
  return (abs(id)==1214 || abs(id)==2124) ? true:false;
}

inline bool isLambda1405(int id) {
  return (abs(id)==13122) ? true:false;
}

inline bool isLambda1520(int id) {
  return (abs(id)==3124) ? true:false;
}

inline bool isSigmaMeson(int id) {
  return id==9000221 ? true:false;
}

}
#endif
