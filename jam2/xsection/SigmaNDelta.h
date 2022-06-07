#ifndef jam2_xsection_SigmaNDelta_h
#define jam2_xsection_SigmaNDelta_h

// Define conversion hbar * c = 0.2 GeV * fm = 1.
//#ifndef HBARC
//#define HBARC 0.19732698
//#endif

//inline double pow2(const double& x) {return x*x;}
//inline double pow3(const double& x) {return x*x*x;}
//inline double pow4(const double& x) {return x*x*x*x;}
// Functions: momentum in two-particle cm.
//inline double PCM(double a,double b,double c) {
//	return sqrt((a*a-(b+c)*(b+c))*(a*a-(b-c)*(b-c)))/(2.0*a);}

namespace jam2 {

class SigmaNR
{
private:
  static const double emnuc, empion, ekinmi,emr,gamr,gr;
  static const int ldec;
  int optWidth;
  double Norm;
  double prres,lam2, facF;
  double *xg, *wg;
  int NG;

public:
  SigmaNR(int opt=2);
  ~SigmaNR();
  double computeNorm();
  double deltaWidth(double emd);
  double sigma(double srt);
  double norm() {return Norm;}

};
}
#endif
