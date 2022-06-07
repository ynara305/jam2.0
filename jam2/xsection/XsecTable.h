#ifndef jam2_xsection_XsecTable_h
#define jam2_xsection_XsecTable_h
#include <cmath>

namespace jam2 {

class XsecTable
{
private:
    static const double sigfit[29][100];
    static const double sigfits1[14][100];
    static const double sigfits2[26][100];
public:
    static double jamsighh(int isig,double srt);
    static double jamchc96(int i,double plab);
    static double jamrgg96(double srt,int i);
    static double pdg2016(double srt,int i);
    static void   additiveQuarkModel(int kf1,int kf2,double& sig,double& sigel);
    static double sigBBS1(int isig,double srt);
    static double sigBBS2(int isig,double srt);

    // Functions: momentum in two-particle cm.
    static double pawt(double a,double b,double c) {
	return sqrt((a*a-(b+c)*(b+c))*(a*a-(b-c)*(b-c)))/(2.0*a);}
    static double pawt2(double a,double b,double c) {
	return (a*a-(b+c)*(b+c))*(a*a-(b-c)*(b-c))/(4.0*a*a);}

    // Function:lab. momentum.
    static double plabsr(double a,double b,double c) {
	return sqrt((a*a-b*b-c*c)*(a*a-b*b-c*c)/(4.0*c*c)-b*b);}
    static double nnElastic(double srt,int iopt);

};
}
#endif
