#include <jam2/interaction/Scatter.h>
#include <Pythia8/PythiaStdlib.h>
#include <jam2/hadrons/JamStdlib.h>
#include <Pythia8/PythiaStdlib.h>

namespace jam2 {

using namespace std;
using Pythia8::pow3;
using Pythia8::pow4;

Scatter::Scatter(Pythia8::Info* inf,Pythia8::Settings* s,JamParticleData* jdata,
	CrossSection* xsec, Pythia8::Pythia* hadronize,
        Pythia8::StringFlav* flavSel, Pythia8::Rndm* r)
    : info(inf), settings(s), rndm(r) 
{
  printColl=settings->mode("Cascade:PrintCollision");
  softModel = settings->mode("Collision:SoftModel");
  eMinPert   = settings->parm("Beams:eMinPert");
  eWidthPert = settings->parm("Beams:eWidthPert");

  scatt = new SoftStrings(info,*settings,jdata,xsec,hadronize,flavSel,rndm);

}

Scatter::~Scatter() 
{
  delete scatt;
}

int Scatter::usePythia(double eCM)
{
    //return 1;   //Pythia only
    //return 0;  // No Pythia

  //if(preHadronA || preHadronB) return 0;

  // currently, hard scattering of secondary interaction is BB only
  //if(abs(inComing[0]->getNColl()*inComing[1]->getNColl()) !=1
  //    && cpair.getCollType() !=1) return 0;

  //if(abs(inComing[0]->getNColl()*inComing[1]->getNColl()) !=0
  //    && cpair.getCollType() ==3) return 0;

  if(channel==DIFFRACTIVE) return 0;

  if ( exp( -(eCM - eMinPert) / eWidthPert ) > rndm->flat()  )
      return 0;

  if( eCM < 500.0) return 1;

  return 2; // use default pythia

}

void Scatter::scatter(InterList* inter0, vector<EventParticle*>& outgoing,
	Collision* event)
{
    TwoBodyInterList *inter=dynamic_cast<TwoBodyInterList*>(inter0);
    channel = scatt->selectChannel(inter,event->getNColl());
    scatt->setChannel(channel);
    CollisionPair& cpair=scatt->getCpair();
    outgoing.clear();

    int ipass=0;
    bool success=true;
    double eCM = scatt->ecm();
    int typeDiff=scatt->typeDiffra();

    const int *constQ[2];
    constQ[0]=scatt->constq(0);
    constQ[1]=scatt->constq(1);

    if(channel == SOFT || channel == DIFFRACTIVE) {
	int id1=cpair.getID(0);
	int id2=cpair.getID(1);

	// s-channel string formation.
        if(cpair.getOutGoing().id2 == 0) {
	  int id3=scatt->idout(0);

	    if(cpair.getOutGoing().id1==93) {
              scatt->setChannel(BBANIHILATION);
	      success=scatt->BaBAnnihilation(outgoing);
	    } else {
              scatt->setChannel(ABSORPTION);
	      success=scatt->absorbS(id1,id2,id3,constQ[0],constQ[1],outgoing);
	    }

	    if(!success)  ipass=1;

	// t-channel string formation.
	} else {


	    if(int usepythia = usePythia(eCM)) {
                scatt->setChannel(HARD);
                Vec4 pc1 = scatt->pc(0);
		Vec4 pc2 = scatt->pc(1);
		const int* idIn = scatt->idin();
		const double* meff = scatt->meff();
		const double* m = scatt->m();
		Vec4 xout1 = scatt->xout(0);
		Vec4 xout2 = scatt->xout(1);
		bool preHadronA = scatt->prehadronA();
		bool preHadronB = scatt->prehadronB();

		if(usepythia==1) {
  		  pythia->initBeam(idIn,meff,m,pc1,pc2,preHadronA,preHadronB);
		  success=pythia->generate(scatt->nColl(),event->getNColl(),
			constQ[0],constQ[1],xout1,xout2,outgoing);
		} else {
		  pythia2->initBeam(idIn,meff,m,pc1,pc2,preHadronA,preHadronB);
		  success=pythia2->generate(scatt->nColl(),event->getNColl(),
			constQ[0],constQ[1],xout1,xout2,outgoing);
		}

		ipass=2;
	    } else {
      
                scatt->setChannel(SOFT);
		ipass=3;
		if(softModel==2) {
		  success=scatt->setKinematics(typeDiff,outgoing);
		} else {
		  success= scatt->hijingSoft(typeDiff,outgoing);
		}

		if(!success) {
		   //cout << " Soft string excitation failed";
		   //cout << " id1= " << id1 << " id2= " << id2
		    //   << " ecm= " << eCM
		     //  << " typeDiff= "<< typeDiff
		      // <<endl;

		  if(outgoing.size()>0) {
		    cout << "Scatter::scatter error after setKinematics " 
			<< " outgoing = " << outgoing.size()
			<< endl;
                    for(int i=0; i<(int)outgoing.size();i++) delete outgoing[i];
		    outgoing.clear();
		  }
		 ipass=4;
	       }
	    }

      }

    // h1 + h2 -> h3.
    } else if(channel == ABSORPTION) {
	scatt->absorb(outgoing);
	ipass=5;

    // h1 + h2 -> h3 + h4.
    } else {
	scatt->scatter2(inter,event,outgoing);
	ipass=6;
    }

    if(!success) {
        channel=ELASTIC;
        scatt->setChannel(channel);
	scatt->scatter2(inter,event,outgoing);
    }

    if(printColl) cout << "after Scatter:: " 
	<< " ipass= "<< ipass
	<<  " outgoing size="<< outgoing.size()
	<< endl;

    if(outgoing.size()>0) scatt->checkMomentum(outgoing);

}

} // end namespace jam2.
