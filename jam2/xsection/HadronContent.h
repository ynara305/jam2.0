#ifndef jam2_xsection_HadronContent_h
#define jam2_xsection_HadronContent_h

#include <string>
#include <Pythia8/ParticleData.h>

namespace jam2 {

using Pythia8::ParticleData;

class HadronContent
{
private:
  ParticleData* particleData;  // particle data table for all particles.
  Pythia8::Rndm* rndm;

public:
  HadronContent(ParticleData* table,Pythia8::Rndm* r)
	: particleData(table), rndm(r) { }
  void findFlavor(int kf,int& ifla,int& iflb);
  int  ifrkfc(int ia,int ib,int ic,double s);
  void findFlavorMeson(int kf,int* ifla,int* iflb, double* prob,int& n);
  double jamemjet(int kfl10,int kfl20);
  int combine(int flav1, int flav2);
};
}
#endif // jam2_xsection_HadronContent_h
