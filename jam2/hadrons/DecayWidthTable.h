#ifndef jam2_hadrons_DecayWidthTable_h
#define jam2_hadrons_DecayWidthTable_h

#include <Pythia8/ParticleData.h>
#include <jam2/hadrons/ParticleTable.h>
#include <jam2/hadrons/DecayWidth.h>

namespace jam2 {

class DecayWidthTable
{
private:
  int nMass;
  double minMass,maxMass,dMass;
  double *width;
  double wid0;
  ParticleDataEntry* part;
public:
  DecayWidthTable(ParticleDataEntry* p, DecayWidth* d) {
    part = p;
    wid0=p->mWidth();
    if(wid0 < 1e-15) {
      width=0;
      nMass=1;
      minMass=p->m0();
      maxMass=p->m0();
      return;
    }
    nMass=500;
    maxMass=6.0;
    minMass=p->m0();
    if(wid0 > 1e-5) minMass=p->mMin();
    dMass = (maxMass - minMass) / (nMass - 1);
    width = new double [nMass];

    // make table loop over mass.
    for(int i=0;i<nMass;i++) {
      width[i]=d->getTotalWidth(p,minMass + i*dMass);
    }

  }
  ~DecayWidthTable() {if(width) delete [] width;}

  double getWidth(double m) {
    if(nMass==1) return wid0;
    int i=(int)((m-minMass)/dMass);
    if(i>=0 && i<nMass-1) {
      double mi=minMass+i*dMass;
      double s = (m-mi)/dMass;
      return width[i]*(1.0-s) + s*width[i+1];
    } else if(i>=maxMass) {
      std::cout << "DecayWidthTable::getWidth mass too large id= "
	  << part->id() << " m= " << m <<std::endl;
      exit(1);
      //return width[nMass-1];
    } else {
      std::cout << "DecayWidthTable::getWidth mass too small "<< m
                << " minMass= "<< minMass
		<< " name= "<< part->name()
		<< " id= "<< part->id()
                << " i= "<< i
                << " floor= " << (int)((m-minMass)/dMass)
                <<std::endl;
      return width[0];
    }
  }


};
} // end namespace jam2
#endif
