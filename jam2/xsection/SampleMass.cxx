#include <jam2/xsection/SampleMass.h>
#include <jam2/hadrons/JamStdlib.h>
#include <jam2/hadrons/GaussPoints.h>
#include <jam2/xsection/CrossSection.h>

namespace jam2 {

using namespace std;
using Pythia8::ParticleDataEntry;

const double SampleMass::Mnucl=0.9383, SampleMass::Mpion=0.138;
//const int SampleMass::nS=500;
const int SampleMass::nS=100;
const double SampleMass::sMinBB=2*Mnucl+Mpion;
const double SampleMass::sMinMM=3*Mpion;
const double SampleMass::sMinMB=Mnucl+2*Mpion;
const double SampleMass::sMaxBB=7.0;
const double SampleMass::sMaxMM=5.0;
const double SampleMass::sMaxMB=6.0;


SampleMass::SampleMass(Pythia8::Settings *s,string bwfile,JamParticleData* table,StringFlav* fs,Pythia8::Rndm* r)
{
  settings = s;
  jamTable=table; 
  //decay=jamTable->getDecayWidth();
  flavSel = fs;
  rndm=r;
  hadronContent = new HadronContent(jamTable->getParticleData(),rndm);

 // Mixing for eta and eta'.
  double theta    = settings->parm("StringFlav:thetaPS");
  double alpha    = (theta + 54.7) * M_PI / 180.;
  fracEtass       = pow2(sin(alpha));
  fracEtaPss      = 1. - fracEtass;

  proton  = jamTable->find(2212);
  neutron = jamTable->find(2112);
  lambda  = jamTable->find(3122);
  sigmam  = jamTable->find(3112);
  sigma0  = jamTable->find(3212);
  sigmap  = jamTable->find(3222);
  xim     = jamTable->find(3312);
  xi0     = jamTable->find(3322);
  omega   = jamTable->find(3334);

  optConstQuarkDiffra=settings->mode("Cascade:optConstQuarkDiffra");
  optQuarkExchange=settings->flag("Cascade:optQuarkExchangeProcess");
  optQuarkAnn = settings->flag("Cascade:optQuarkAnnihilationProcess");

  // Option for BW integration.
  optProb=3; // use table
  //optProb=2;// full two-dim integral
  //optProb=1;// use approximate formula
  //optProb=0;// neglect BW integration.
 
  //optWidth=0; // constant width
  optWidth=1; // momentum-dependent width

  //ParticleTable* meson=jamTable->getMeson();
 
  // sigma, eta, omega,...
  ParticleTable* meson = jamTable->getLight0Meson();
  for(int i=0;i<(int)meson->size();i++) {
    mesons.push_back(meson->getParticle(i));
    isoParticle[meson->getParticle(i)]=meson->getParticle(i);
  }

  // pion,rho...
  ParticleTable* meson1 = jamTable->getLight1Meson0();
  ParticleTable* meson1p = jamTable->getLight1Mesonp();
  for(int i=0;i<(int)meson1->size();i++) {
    mesons.push_back(meson1->getParticle(i));
    isoParticle[meson1->getParticle(i)]=meson1->getParticle(i);
    isoParticle[meson1p->getParticle(i)]=meson1->getParticle(i);
  }

  // K*
  ParticleTable *strp = jamTable->getStrMesonp();
  ParticleTable *str0 = jamTable->getStrMeson0();
  for(int i=0;i<(int)strp->size();i++) {
    mesons.push_back(strp->getParticle(i));
    isoParticle[strp->getParticle(i)]=strp->getParticle(i);
    isoParticle[str0->getParticle(i)]=strp->getParticle(i);
  }

  // N*
  baryons.push_back(proton);
  isoParticle[proton]=proton;
  isoParticle[neutron]=proton;
  ParticleTable* pstar=jamTable->getPstar();
  ParticleTable* nstar=jamTable->getNstar();
  for(int i=0;i<(int)pstar->size();i++) {
    baryons.push_back(pstar->getParticle(i));
    isoParticle[pstar->getParticle(i)]=pstar->getParticle(i);
    isoParticle[nstar->getParticle(i)]=pstar->getParticle(i);
  }

  // Delta*
  ParticleTable* dm=jamTable->getDmstar();
  ParticleTable* d0=jamTable->getD0star();
  ParticleTable* dp=jamTable->getDpstar();
  ParticleTable* dpp=jamTable->getDppstar();
  for(int i=0;i<(int)dm->size();i++) {
    baryons.push_back(dm->getParticle(i));
    isoParticle[dm->getParticle(i)]=dm->getParticle(i);
    isoParticle[d0->getParticle(i)]=dm->getParticle(i);
    isoParticle[dp->getParticle(i)]=dm->getParticle(i);
    isoParticle[dpp->getParticle(i)]=dm->getParticle(i);
  }

  // Lambda*
  baryons.push_back(lambda);
  isoParticle[lambda]=lambda;
  ParticleTable* lam=jamTable->getLambda();
  for(int i=0;i<(int)lam->size();i++) {
    baryons.push_back(lam->getParticle(i));
    isoParticle[lam->getParticle(i)]=lam->getParticle(i);
  }

  // Sigma*
  baryons.push_back(sigmam);
  isoParticle[sigmam]=sigmam;
  isoParticle[sigma0]=sigmam;
  isoParticle[sigmap]=sigmam;
  ParticleTable* sm=jamTable->getSmstar();
  ParticleTable* s0=jamTable->getS0star();
  ParticleTable* sp=jamTable->getSpstar();
  for(int i=0;i<(int)sm->size();i++) {
    baryons.push_back(sm->getParticle(i));
    isoParticle[sm->getParticle(i)]=sm->getParticle(i);
    isoParticle[s0->getParticle(i)]=sm->getParticle(i);
    isoParticle[sp->getParticle(i)]=sm->getParticle(i);
  }
 
  // Xi*
  baryons.push_back(xim);
  isoParticle[xi0]=xim;
  ParticleTable *xism=jamTable->getXmstar();
  ParticleTable *xis0=jamTable->getX0star();
  for(int i=0;i<(int)xism->size();i++) {
    baryons.push_back(xism->getParticle(i));
    isoParticle[xism->getParticle(i)]=xism->getParticle(i);
    isoParticle[xis0->getParticle(i)]=xism->getParticle(i);
  }

  baryons.push_back(omega);
  isoParticle[omega]=omega;

  if(optProb==3) if(!readBWTable(bwfile)) makeBWTable(bwfile);

}

SampleMass::~SampleMass()
{
  for(int i=0;i<(int)BWintSave.size();i++) delete BWintSave[i];
  BWintSave.clear();
  delete hadronContent;
}

bool SampleMass::readBWTable(string fname)
{
  ifstream inFile(fname.c_str(),ios::in);
  if(!inFile.is_open()) {
    cout << "Sample::readBWTable file does not exist " << fname
         << " now making table" << endl;
    return false;
  }

  std::vector<double> bwint(nS);

  vector<ParticleDataEntry*> hadrons=mesons;
  hadrons.insert(hadrons.end(),baryons.begin(),baryons.end());
  for(int i=0;i<(int)hadrons.size();i++)
  for(int j=i;j<(int)hadrons.size();j++) {
  //if(hadrons[i]->mWidth() <1e-5 && hadrons[j]->mWidth() <1e-5) continue;
    if(hadrons[i]->mWidth() <1e-5 || hadrons[j]->mWidth() <1e-5) continue;
    double sMin=sMinMB, sMax=sMaxMB;
    if(hadrons[i]->isMeson() && hadrons[j]->isMeson()) {
      sMin=sMinMM; sMax=sMaxMM;
    } else if(hadrons[i]->isBaryon() && hadrons[j]->isBaryon()) {
      sMin=sMinBB; sMax=sMaxBB;
    }
    for(int k=0;k<nS;k++) inFile >> bwint[k];
    BWintegral *bw = new BWintegral(sMin,sMax,bwint);
    BWint[make_pair(hadrons[i],hadrons[j])]=bw;
    if(i != j) BWint[make_pair(hadrons[j],hadrons[i])]=bw;
    BWintSave.push_back(bw);
  }

/*

  // baryon-baryon
  for(int i=0;i<(int)baryons.size();i++)
  for(int j=i;j<(int)baryons.size();j++) {
    for(int k=0;k<nS;k++) inFile >> bwint[k];
    BWintegral *bw = new BWintegral(nS,sMinBB,sMaxBB,bwint);
    BWint[make_pair(baryons[i],baryons[j])]=bw;
    if(i != j) BWint[make_pair(baryons[j],baryons[i])]=bw;
  }

  // meson-meson
  for(int i=0;i<(int)mesons.size();i++)
  for(int j=i;j<(int)mesons.size();j++) {
    for(int k=0;k<nS;k++) inFile >> bwint[k];
    BWintegral *bw = new BWintegral(nS,sMinMM,sMaxMM,bwint);
    BWint[make_pair(mesons[i],mesons[j])]=bw;
    if(i != j) BWint[make_pair(mesons[j],mesons[i])]=bw;
  }

  // meson-baryon
  for(int i=0;i<(int)mesons.size();i++)
  for(int j=0;j<(int)baryons.size();j++) {
    for(int k=0;k<nS;k++) inFile >> bwint[k];
    BWintegral *bw = new BWintegral(nS,sMinMB,sMaxMB,bwint);
    BWint[make_pair(mesons[i],baryons[j])]=bw;
    BWint[make_pair(baryons[j],mesons[i])]=bw;
    cout << "kf1= "<< mesons[i]->id()
	<< " kf2= "<< baryons[i]->id()
	<< " bw= "<< BWint[make_pair(mesons[i],baryons[j])]
	<<endl;
  }

*/

  inFile.close();
  return true;
}

void SampleMass::makeBWTable(string outfile)
{
  ofstream ofs(outfile.c_str());

  vector<ParticleDataEntry*> hadrons=mesons;
  hadrons.insert(hadrons.end(),baryons.begin(),baryons.end());
  std::vector<double> bwint(nS);
  for(int i=0;i<(int)hadrons.size();i++)
  for(int j=i;j<(int)hadrons.size();j++) {
    if(hadrons[i]->mWidth() <1e-5 || hadrons[j]->mWidth() <1e-5) continue;
    double sMin=sMinMB, sMax=sMaxMB;
    if(hadrons[i]->isMeson() && hadrons[j]->isMeson()) {
      sMin=sMinMM; sMax=sMaxMM;
    } else if(hadrons[i]->isBaryon() && hadrons[j]->isBaryon()) {
      sMin=sMinBB; sMax=sMaxBB;
    }
    double dS=(sMax-sMin)/nS;

    for(int k=0;k<nS;k++) {
      double srt = sMin + dS * k;
      bwint[k] = probBW2(srt,hadrons[i],hadrons[j]);
      ofs << bwint[k] << endl;
    }

    BWintegral *bw = new BWintegral(sMin,sMax,bwint);
    BWint[make_pair(hadrons[i],hadrons[j])]=bw;
    if(i != j) BWint[make_pair(hadrons[j],hadrons[i])]=bw;
    BWintSave.push_back(bw);
  }

/*
  // Baryon-Baryon
  double *bw1 = new double [nS];
  double dS=(sMaxBB-sMinBB)/nS;
  for(int i=0;i<(int)baryons.size();i++)
  for(int j=i;j<(int)baryons.size();j++) {
    for(int k=0;k<nS;k++) {
      double srt = sMinBB + dS * k;
      bw1[k] = probBW2(srt,baryons[i],baryons[j]);
      ofs << bw1[k] << endl;
    }
    BWintegral *bw = new BWintegral(nS,sMinBB,sMaxBB,bw1);
    BWint[make_pair(baryons[i],baryons[j])]=bw;
    if(i != j) BWint[make_pair(baryons[j],baryons[i])]=bw;
  }
  delete [] bw1;

  // Meson-Meson
  double *bw2 = new double [nS];
  dS=(sMaxMM-sMinMM)/nS;
  for(int i=0;i<(int)mesons.size();i++)
  for(int j=i;j<(int)mesons.size();j++) {
    for(int k=0;k<nS;k++) {
      double srt = sMinMM + dS * k;
      bw2[k] = probBW2(srt,mesons[i],mesons[j]);
      ofs << bw2[k] << endl;
    }
    BWintegral *bw = new BWintegral(nS,sMinMM,sMaxMM,bw2);
    BWint[make_pair(mesons[i],mesons[j])]=bw;
    if(i != j) BWint[make_pair(mesons[j],mesons[i])]=bw;
  }
  delete [] bw2;

  // Meson-Baryon
  double *bw3 = new double [nS];
  dS=(sMaxMB-sMinMB)/nS;
  for(int i=0;i<(int)mesons.size();i++)
  for(int j=0;j<(int)baryons.size();j++) {
    for(int k=0;k<nS;k++) {
      double srt = sMinMB + dS * k;
      bw3[k] = probBW2(srt,mesons[i],baryons[j]);
      ofs << bw3[k] << endl;
    }
    BWintegral *bw = new BWintegral(nS,sMinMB,sMaxMB,bw3);
    BWint[make_pair(mesons[i],baryons[j])]=bw;
    BWint[make_pair(baryons[j],mesons[i])]=bw;
  }
  delete [] bw3;

*/

}

//***********************************************************************
//      subroutine jamrmas2(kf1,kf2,kc1,kc2,srt,em1,em2,icon)
//...Generate masses according to the Breit-Wigner distribution.
int SampleMass::jamrmas2(int coltype,ParticleDataEntry* pd1, ParticleDataEntry* pd2,
	int kf1,int kf2,int id1,int id2, int iz1, int iz2, 
	double m1, double m2, bool preHadronA, bool preHadronB,double srt,bool isAnti)
{
  // elastic
  //return 1;

  // currently charmed hadrons are not implemented.
  if(CrossSection::getHeavyQB(kf1)+CrossSection::getHeavyQB(kf2)>0) return 1;

  eCM = srt;
  pout[0] = pd1;
  pout[1] = pd2;

  Id[0] = kf1;
  Id[1] = kf2;
  emr[0] = m1;
  emr[1] = m2;
  iZ1 = iz1;
  iZ2 = iz2;
  iD1 = id1;
  iD2 = id2;

  collPair.clear();

  if(preHadronA || preHadronB) {
    diffractive(1.0);
  } else {

    if(optQuarkExchange && optQuarkAnn) {
      diffractive(0.3);
      //quarkExchangeProcess(Id[0],Id[1],0.4);
      quarkExchangeProcess2(Id[0],Id[1],0.4);
      if(coltype!=1) quarkAnnProcess(Id[0],Id[1],0.3);
    } else if(optQuarkExchange && !optQuarkAnn) {
      diffractive(0.5);
      //quarkExchangeProcess(Id[0],Id[1],0.5);
      quarkExchangeProcess2(Id[0],Id[1],0.4);
    } else if(!optQuarkExchange && optQuarkAnn) {
      diffractive(0.5);
      if(coltype!=1) quarkAnnProcess(Id[0],Id[1],0.5);
    } else if(!optQuarkExchange && !optQuarkAnn) {
      diffractive(1.0);
    }


    /*
    if(optQuarkExchange || optQuarkAnn) {
      std::array<int,3> iq1, iq2;
      quarkContent(Id[0],iq1);
      quarkContent(Id[1],iq2);

      if(optQuarkExchange) {
        quarkExchangeProcess(iq1,iq2,0.5);
      }
      if(optQuarkAnn) {
        if(coltype!=1) quarkAnnProcess(iq1,iq2,0.3);
      }
    }
    */

    /*
    for(auto& i : collPair) {
      cout << pd1->name(kf1)  << " + " << pd2->name(kf2) << " -> "
	<< i.particleA()->name(i.particleA()->id())
	<< " + " 
	<< i.particleB()->name(i.particleB()->id())
	<<endl;
    }
    */
  }

  double totsp=0.0;
  for(auto& i : collPair) totsp += i.probability();

  if(totsp == 0.0) {
    cout << "jamrmas2 totsp=0 " << totsp 
      << " n= " << collPair.size()
      << " srt= " << srt <<endl;
    cout << " kf1= " << kf1 << " kf2= " << kf2 << endl;
    cout << " m1= " << m1 << " m2= " << m2 << endl;
    if(collPair.size()==0) return 1;

    int imin=-1;
    double mmin=10000.;
    int i=0;
    for(auto& cp : collPair) {
      ParticleDataEntry* p1=cp.particleA();
      ParticleDataEntry* p2=cp.particleB();
      if(Id[0]==p1->id() && Id[1]==p2->id()) continue;
      double mmin1= p1->mWidth()>1e-5 ? p1->mMin() : p1->m0();
      double mmin2= p2->mWidth()>1e-5 ? p2->mMin() : p2->m0();
      if(mmin1+mmin2<mmin) {
	mmin = mmin1+mmin2;
	imin=i;
      }
      i++;
      /*
      cout << " p1= "<< p1->name() 
	<< " p2= "<< p2->name() 
	<< " id1= "<< p1->id()
	<< " id2= "<< p2->id()
	<< " m1= " << mmin1
	<< " m2= " << mmin2
	<< " m1+m2= " << mmin1+mmin2
	<< " prob= " << collPair[imin].probability()
	<<endl;
	*/
    }

      ParticleDataEntry* p1=collPair[imin].particleA();
      ParticleDataEntry* p2=collPair[imin].particleB();
      double mmin1= p1->mWidth()>1e-5 ? p1->mMin() : p1->m0();
      double mmin2= p2->mWidth()>1e-5 ? p2->mMin() : p2->m0();
      cout << " p1= "<< p1->name() 
	<< " p2= "<< p2->name() 
	<< " id1= "<< p1->id()
	<< " id2= "<< p2->id()
	<< " m1= " << mmin1
	<< " m2= " << mmin2
	<< " m1+m2= " << mmin1+mmin2
	<< " prob= " << collPair[imin].probability()
	<<endl;

    return 1;
  }

  // Select outgoing particles.
  int ntry=0;
  do {
    double  xrand=totsp*rndm->flat();
    double ttp=0.0;
    for(auto& cp : collPair) {
      xrand -= cp.probability();
      ttp += cp.probability();
      if(xrand <= 0.0) {
	pout[0]=cp.particleA();
	pout[1]=cp.particleB();
	Id[0]=cp.id1();
	Id[1]=cp.id2();
	goto L100;
      }
    }
    cout << "SampleMass::  no particles? totsp= " << totsp
      << " ttp= "<< ttp << endl;
    return 1;
L100:


    // Monte Carlo sampling of particle masses.
    if(sample()) break;

  } while(ntry++<100);

  if(ntry==100) {
    cout << "SampleMass: does not converge eCM= " << srt
      << " id1= "<< Id[0]
      << " id2= "<< Id[1]
      << endl;
    return 1;
  }

  /*
  if(!pout[0]->hasAnti()) kfsign1=1;
  if(!pout[1]->hasAnti()) kfsign2=1;
  Id[0] = (pout[0]->id())*kfsign1;
  Id[1] = (pout[1]->id())*kfsign2;
  */

  //Id[0] = (pout[0]->id());
  //Id[1] = (pout[1]->id());

  if(isAnti) {
    if(pout[0]->hasAnti()) Id[0] *= -1;
    if(pout[1]->hasAnti()) Id[1] *= -1;
  }

  //cout << " id1= "<< kf1 << " kf2= "<< kf2 <<endl;
  //cout << " id3= "<< Id[0] << " id4= "<< Id[1] <<endl;

  /*
  ParticleDataEntryPtr p3 = jamTable->find(Id[0]);
  ParticleDataEntryPtr p4 = jamTable->find(Id[1]);
    if(p3==0 || p4==0) {
      cout << "p=0????"<< p3 << " p2= "<< p4 <<endl;
      exit(1);
    }
    */

  /*
  if(pd1->isBaryon() || pd2->isBaryon() || pout[0]->isBaryon() || pout[1]->isBaryon()) {
      cout << " p1= "<< pd1->name() 
	<< " p2= "<< pd2->name() 
	<< " id1= "<< pd1->id()
	<< " id2= "<< pd2->id()
	<< " id3= " << Id[0]
	<< "  " << pout[0]->name()
	<< " id4= " << Id[1]
	<< "  " << pout[1]->name()
	<<endl;
  }
  */

  return 0;

}

vector<ParticleDataEntry*> SampleMass::
findMeson(ParticleDataEntry* pa, int kf0)
{
  //....Mesons.

  ParticleTable* table=jamTable->getMeson();
  vector<ParticleDataEntry*> ncount;
  int kfsign= kf0 > 0 ? 1 : -1;
  if(pa->hasAnti()) kfsign=1;

  int kfa=abs(kf0);
  int kf1=(kfa/100) % 10;
  int kf2=(kfa/10)  % 10;
  int kfm=100*kf1+10*kf2;
  if(kfm == 220 || kfm == 330) kfm=110;
  //int kfg=(kfm+1)*kfsign;

  int itry=0, jtry=0;
  do {
    if(itry++ > 200) {
      cout << "SamplemMass::findMeson infinit loop? " << kf0 <<endl;
      exit(1);
    }

    //9000111
    //int icount=0;
    for(int i9=0;i9<=9; i9+=9) 
      for(int ir=0;ir<=2; ir++) 
	//for(int ir=0;ir<=5; ir++) 
	for(int is=1;is<=7; is+=2) {
	  //int kfmes=(ir*10000+kfm+is)*kfsign;
	  int kfmes=(i9*1000000+ir*10000+kfm+is)*kfsign;
	  //ParticleDataEntry* pr=jamTable->find(kfmes);
	  ParticleDataEntry* pr=table->find(kfmes);
	  if(pr == 0) continue;
	  ncount.push_back(pr);
	}

    jtry=0;
    if(kfm == 220) {
      kfm=330;
      jtry=1;
    } else if(kfm == 110) {
      kfm=220;
      jtry=2;
    } else if(kfm == 330) {
      kfm=110;
      jtry=3;
    }

  } while (itry <= 2 && jtry >= 1);

  return ncount;

}

double SampleMass::computeProb(ParticleDataEntry* p1, ParticleDataEntry* p2)
{
  // exclude elastic collision.
  if(Id[0]==p1->id() && Id[1]==p2->id()) return 0.0;

  // check energy.
  double mmin1= p1->mWidth()>1e-5 ? p1->mMin() : p1->m0();
  double mmin2= p2->mWidth()>1e-5 ? p2->mMin() : p2->m0();
  if(eCM < mmin1 + mmin2 + eKinMin) return 0.0;

  double pbw=0.0;
  if(optProb==3) // use table
    pbw = probBW3(eCM,p1,p2);
  else if(optProb==0 || optProb==1)  // simplified method
    pbw = probBW1(eCM,p1,p2);
  else
    pbw = probBW2(eCM,p1,p2);

  if(pbw <= 0.0) return 0.0;

  double fac= p1->id() == p2->id() ? 0.5 : 1.0;
  int ispin1=p1->spinType();
  int ispin2=p2->spinType();
  return fac * ispin1*ispin2 * pbw;

}

// Before call this, specify the incoming particles pout[2].
bool SampleMass::sample()
{    
  double em0[2],gam0[2],gamw[2];
  double mmax[2],mmin[2],bwmax[2],bw[2],ymin[2],ymax[2];
  int iex[2], id[2];

  bool ok=false;
  double wid1=pout[0]->mWidth();
  double wid2=pout[1]->mWidth();
  mmin[0]=pout[0]->m0();
  mmin[1]=pout[1]->m0();
  iex[0]=0; iex[1]=0;
  if(wid1>1e-5) {mmin[0] = pout[0]->mMin(); iex[0]=1;}
  if(wid2>1e-5) {mmin[1] = pout[1]->mMin(); iex[1]=1;}

  mmax[0]=mmin[0];
  mmax[1]=mmin[1];
  //if(iex[0]==1) mmax[0]=min(3.5,ecm-mmin[1]-eKinMin);
  //if(iex[1]==1) mmax[1]=min(3.5,ecm-mmin[0]-eKinMin);
  if(iex[0]==1) mmax[0]=max(mmin[0],eCM-mmin[1]-eKinMin);
  if(iex[1]==1) mmax[1]=max(mmin[1],eCM-mmin[0]-eKinMin);

  if(mmin[0] <= mmax[0] && mmin[1] <= mmax[1]) { ok=true; }

  if(!ok) {
    id[0]=pout[0]->id();
    id[1]=pout[1]->id();
    cout << "SampleMass::ResnanceMass too small srt? srt= " << eCM <<endl;
    cout << " ok= " << ok << endl;
    cout << " id1= " << id[0] << " min1= "<< mmin[0]<< " max= "<< mmax[0]
      << " iex= " << iex[0] <<endl;
    cout << " id2= " << id[1] << " min2= " << mmin[1]<< " max= "<< mmax[1]
      << " iex= " << iex[1] <<endl;
    exit(1);
  }

  if(eCM < mmin[0] + mmin[1] + eKinMin) {
    cout << "sampleResonanceMass too small srt? " << eCM
      << " mmin1= " << mmin[0]
      << " mmin2= " << mmin[1]
      << " id1= " << Id[0]
      << " id2= " << Id[1]
      <<endl;
    exit(1);
  }

  id[0]=pout[0]->id();
  id[1]=pout[1]->id();

  // First compute maximum value.
  for(int jt=0; jt<2;jt++) {
    bwmax[jt]=1.0;
    bw[jt]=1.0;
    ymin[jt]=0.0;
    ymax[jt]=0.0;
    em0[jt]=pout[jt]->m0();
    gam0[jt]=pout[jt]->mWidth();
    emr[jt]=em0[jt];
    if(iex[jt] == 0) continue;
    ymin[jt]=atan((pow2(mmin[jt])-pow2(em0[jt]))/(em0[jt]*gam0[jt]));
    ymax[jt]=atan((pow2(mmax[jt])-pow2(em0[jt]))/(em0[jt]*gam0[jt]));
    //gamw[jt] = decay->getTotalWidth(pout[jt],mmax[jt]);
    //bwmax[jt] = BW(em0[jt],mmax[jt],gamw[jt])/BW(em0[jt],mmax[jt],gam0[jt]);
    gamw[jt] = jamTable->totalWidth(pout[jt],em0[jt]);
    bwmax[jt] = BW(em0[jt],em0[jt],gamw[jt])/BW(em0[jt],em0[jt],gam0[jt]);

    //int pid=jamTable->pid(id[jt]);
    //if(pid == id_nucls) bwmax[jt] *= 2.0;
    //if(pid == id_delts) bwmax[jt] *= 1.5;
    //if(pid == id_delt) bwmax[jt] *= 1.5;
    //if(pid == id_str) bwmax[jt] *= 1.5;
    if(id[jt] == 333) bwmax[jt] *= 3.0;
    if(id[jt] == 313 || id[jt]==323) bwmax[jt] *= 2.0;
  }

  double pf0=PCM(eCM,mmin[0],mmin[1]);
  double bwpmax=bwmax[0]*bwmax[1]*pf0 * 5.0;

  // Generate resonance mass according to Breit-Wigner + phase space.
  int maxtry=10;
  int ntry=0;
  double pf=0.0;
  do {
    do {
      for(int jt=0;jt<2;jt++)
	if(iex[jt] == 1) {
	  // First generate mass according to B-W with constant width.
	  // A(m^2)=d(m^2)/((m^2 - m_R^2)^2 + (m*Gamma)^2)
	  double y=ymin[jt] + rndm->flat()*( ymax[jt] - ymin[jt] );
	  emr[jt] = sqrt( em0[jt]*gam0[jt]*tan(y) + em0[jt]*em0[jt] );
	}
    } while(emr[0]+emr[1]+eKinMin > eCM);

    for(int jt=0;jt<2;jt++)
      if(iex[jt] == 1) {
	double wid = jamTable->totalWidth(pout[jt],emr[jt]);
	bw[jt] = BW(em0[jt],emr[jt],wid)/BW(em0[jt],emr[jt],gam0[jt]);
	//bw[jt] = BW(em0[jt],emr[jt],gamw[jt])/BW(em0[jt],emr[jt],gam0[jt]);
      }

    // Check final phase space.
    pf=PCM(eCM,emr[0],emr[1]);
    if(bwpmax < bw[0]*bw[1]*pf) {
      cout << "SampleMass:bwm= " << bwpmax
	<< " bw= "<< bw[0]*bw[1]*pf
	<< " id1= " << id[0]
	<< " id2= " << id[1]
	<< " srt= " << eCM  << endl;
      bwpmax *= 1.05;
      continue;
    }

    if(ntry++ > maxtry) return false;

  } while (bwpmax*rndm->flat() > bw[0]*bw[1]*pf);

  return true;
}



//***********************************************************************
//      subroutine jamexpa(kf0,em,icount,ncount,ispin,igr,ig)
vector<ParticleDataEntry*> SampleMass::
jamexpa(ParticleDataEntry* pa, int kf0,int id0, int iz0)
{
  //...Purpose: to find possible hadronic excitation states.

  vector<ParticleDataEntry*> ncount;

  if(pa->isMeson()) {

    return findMeson(pa,kf0);

  } else if(pa->isBaryon()) { 

    if(kf0 < 0) iz0 *= -1;

    if(id0 == id_nucl || id0 == id_nucls
	|| id0 == id_delt || id0 == id_delts) {
      //....Find N* D* channel.

      if(iz0 == -1) {
	ParticleTable* dm=jamTable->getDmstar();
	for(int i=0;i<(int)dm->size();i++) {
	  ncount.push_back(dm->getParticle(i));
	}
      } else if(iz0 == 0) {
	ncount.push_back(neutron);
	ParticleTable* d0=jamTable->getD0star();
	for(int i=0;i<(int)d0->size();i++) {
	  ncount.push_back(d0->getParticle(i));
	}
	ParticleTable* n0=jamTable->getNstar();
	for(int i=0;i<(int)n0->size();i++) {
	  ncount.push_back(n0->getParticle(i));
	}
      } else if(iz0 == 1) {
	ncount.push_back(proton);
	ParticleTable* dp=jamTable->getDpstar();
	for(int i=0;i<(int)dp->size();i++) {
	  ncount.push_back(dp->getParticle(i));
	}
	ParticleTable* ps=jamTable->getPstar();
	for(int i=0;i<(int)ps->size();i++) {
	  ncount.push_back(ps->getParticle(i));
	}
      } else if(iz0 == 2) {
	ParticleTable* dpp=jamTable->getDppstar();
	for(int i=0;i<(int)dpp->size();i++) {
	  ncount.push_back(dpp->getParticle(i));
	}
      }

      //...Lambda/Sigma
    } else if(id0 == id_lambda || id0 == id_lambdas
	|| id0 == id_sigma || id0 == id_sigmas) {

      if(iz0 == -1) {
	ncount.push_back(jamTable->find(3112));
	ParticleTable* table=jamTable->getSmstar();
	for(int i=0;i<(int)table->size();i++) {
	  ncount.push_back(table->getParticle(i));
	}
      } else if(iz0 == 0) {
	ncount.push_back(jamTable->find(3212));
	ParticleTable* table=jamTable->getLambda();
	for(int i=0;i<(int)table->size();i++) {
	  ncount.push_back(table->getParticle(i));
	}
	table=jamTable->getS0star();
	for(int i=0;i<(int)table->size();i++) {
	  ncount.push_back(table->getParticle(i));
	}
      } else if(iz0 == 1) {
	ncount.push_back(jamTable->find(3222));
	ParticleTable* table=jamTable->getSpstar();
	for(int i=0;i<(int)table->size();i++) {
	  ncount.push_back(table->getParticle(i));
	}
      }

      //...Xi*
    } else if(id0 == id_xi || id0 == id_xis) {

      if(iz0 == -1) {
	ncount.push_back(jamTable->find(3312));
	ParticleTable* table=jamTable->getXmstar();
	for(int i=0;i<(int)table->size();i++) {
	  ncount.push_back(table->getParticle(i));
	}
      } else if(iz0 == 0) {
	ncount.push_back(jamTable->find(3322));
	ParticleTable* table=jamTable->getX0star();
	for(int i=0;i<(int)table->size();i++) {
	  ncount.push_back(table->getParticle(i));
	}
      }

    } else {
      ncount.push_back(pa);
    }

  } else {
    ncount.push_back(pa);
  }

  if(ncount.size()==0) {
    cout << "SampleMass:jamexpa no particles ? " << kf0 << " id0= " << id0
      << " iz= "<< iz0
      << endl;
    exit(1);
  }

  return ncount;
}


void SampleMass::quarkContent(int id, std::array<int,3>& iq)
{
  int idabs = abs(id);
  iq[0]=(idabs/1000)  % 10;
  iq[1]=(idabs/100)   % 10;
  iq[2]=(idabs/10)    % 10;
  //int  iq1=idabs         % 10;      // spin

  if(iq[0]==0 && iq[1] != iq[2]) {
    if (iq[1]%2 == 1) std::swap( iq[1], iq[2]);
    if (id > 0) iq[2] = -iq[2];
    else {
      iq[1] = -iq[1];
    }
  } else if(iq[0] ==0) {
    int id2=iq[1]*10 + iq[2];
    if (iq[1] < 3 || id2 == 33) {
      iq[1] = (rndm->flat() < 0.5) ? 1 : 2;
      // eta and eta' can also be s sbar.
      if (id2 == 22 && rndm->flat() < fracEtass)  {iq[1] = 3, iq[2]=-3;}
      if (id2 == 33 && rndm->flat() < fracEtaPss) {iq[1] = 3, iq[2]=-3;}
    }
    iq[2] = -iq[1];

    // anti-baryon
  } else if(id < 0) {
    iq[0]= -iq[0];
    iq[1]= -iq[1];
    iq[2]= -iq[2];
  }
}


void SampleMass::diffractive(double fac)
{
  vector<ParticleDataEntry*> pcount1 = jamexpa(pout[0],Id[0],iD1,iZ1);
  vector<ParticleDataEntry*> pcount2 = jamexpa(pout[1],Id[1],iD2,iZ2);
  int n1 = pcount1.size();
  int n2 = pcount2.size();
  int kfsign1 = Id[0] > 0 ? 1 : -1;
  int kfsign2 = Id[1] > 0 ? 1 : -1;
  for(int ik1=0;ik1<n1;ik1++) 
    for(int ik2=0;ik2<n2;ik2++) {
      double prob = computeProb(pcount1[ik1],pcount2[ik2]);
      int id1 = (pcount1[ik1]->id())*kfsign1;
      int id2 = (pcount2[ik2]->id())*kfsign2;
      collPair.push_back(CollPair(pcount1[ik1],pcount2[ik2],id1,id2,prob*fac));
    }
}

void SampleMass::quarkExchangeProcess(int id1, int id2, double fac)
{
  std::array<int,3> iq1, iq2;
  quarkContent(id1,iq1);
  quarkContent(id2,iq2);

  vector<pair<int,int> > qpair;
  qpair.clear();

  // quark exchange process.
  for(int i=0; i<3;i++) {
    if(iq1[i]==0) continue;
    for(int j=0; j<3;j++) {
      if(iq2[j]==0) continue;
      if(iq1[i]==iq2[j]) continue;
      if(iq1[i]*iq2[j] < 0) continue;
      qpair.push_back(make_pair(i,j));
    }
  }

  if(qpair.size()==0) return;

  // chose exchange quark pair randomly
  int iq=rndm->flat()*qpair.size();
  int i1 = qpair[iq].first;
  int i2 = qpair[iq].second;
  array<int,3> kq1=iq1;
  array<int,3> kq2=iq2;
  kq1[i1]=iq2[i2];
  kq2[i2]=iq1[i1];

  int q1=kq1[1];
  int qq1=kq1[2];
  int q2=kq2[1];
  int qq2=kq2[2];
  // find diquark
  if(iq1[0] !=0) {
    int j=0;
    int jq1[2]={0,0};
    for(int i=0;i<3;i++) {
      if(i != i1) jq1[j++]=iq1[i];
    }
    int spin=3;
    if(jq1[0] != jq1[1])  if(rndm->flat() <0.5) spin=1;
    qq1= 1000*max(abs(jq1[0]) ,abs(jq1[1]))+ 100*min(abs(jq1[0]),abs(jq1[1]))+spin;
    if(Id[0] <0) qq1 *= -1;
    q1=iq2[i2];
  }
  if(iq2[0] !=0) {
    int j=0;
    int jq2[2]={0,0};
    for(int i=0;i<3;i++) {
      if(i != i2) jq2[j++]=iq2[i];
    }
    int spin=3;
    if(jq2[0] != jq2[1])  if(rndm->flat() <0.5) spin=1;
    qq2= 1000*max(abs(jq2[0]) ,abs(jq2[1]))+ 100*min(abs(jq2[0]),abs(jq2[1]))+spin;
    if(Id[1] <0) qq1 *= -1;
    q2=iq1[i1];
  }

  FlavContainer flav1(q1);
  FlavContainer flav2(qq1);
  int id3=0;
  do {
    id3 = flavSel->combine(flav1, flav2);
  }while(id3==0);


  FlavContainer flav3(q2);
  FlavContainer flav4(qq2);
  int id4=0;
  do {
    id4 = flavSel->combine(flav3, flav4);
  } while(id4==0);

  ParticleDataEntry* p3 = jamTable->find(id3);
  ParticleDataEntry* p4 = jamTable->find(id4);

  bool error=false;
  if(p3==0 ) {
    error=true;
    cout << "quark exchange wrong id3 id3= "<< id3 << " p3= "<< p3<<endl;
  }
  if(p4==0 ) {
    error=true;
    cout << "quark exchange wrong id4 id4= "<< id4 << " p4= "<< p4<<endl;
  }
  if(p3 !=0 && p4 != 0) {
    int ch0=pout[0]->chargeType(Id[0]) + pout[1]->chargeType(Id[1]);
    int ch1=p3->chargeType(id3) + p4->chargeType(id4);
    if(ch0 != ch1) {
      cout << "quark exchange charge does not conserve ch0= "<< ch0 << " ch1= "<< ch1<<endl;
      error=true;
    }
  }

  if(error) {
    cout << " i1= "<< i1 << " i2= "<<i2<<endl;
    cout << " iq1= "<< iq1[0] << " " << iq1[1] << "  "<< iq1[2]  <<endl;
    cout << " iq2= "<< iq2[0] << " " << iq2[1] << "  "<< iq2[2]  <<endl;
    cout << "q1 = "<< q1 <<  "  qq1= "<< qq1 <<endl;
    cout << "q2 = "<< q2 <<  "  qq2= "<< qq2 <<endl;
    cout << "id1= "<< Id[0] << " id2=  "<< Id[1] << " id3= "<< id3 << " id4= "<< id4<<endl;
   exit(1);
  }

  int iz3=p3->chargeType(id3)/3;
  int iz4=p4->chargeType(id4)/3;
  int pd3=jamTable->pid(id3);
  int pd4=jamTable->pid(id4);

  vector<ParticleDataEntry*> pcount3=jamexpa(p3,id3,pd3,iz3);
  vector<ParticleDataEntry*> pcount4=jamexpa(p4,id4,pd4,iz4);

  int n3 = pcount3.size();
  int n4 = pcount4.size();
  int kfsign1 = id3 > 0 ? 1 : -1;
  int kfsign2 = id4 > 0 ? 1 : -1;
  for(int ik1=0;ik1<n3;ik1++) 
    for(int ik2=0;ik2<n4;ik2++) {
      double prob = computeProb(pcount3[ik1],pcount4[ik2]);
      int ido1 = (pcount3[ik1]->id())*kfsign1;
      int ido2 = (pcount4[ik2]->id())*kfsign2;
      collPair.push_back(CollPair(pcount3[ik1],pcount4[ik2],ido1,ido2,fac*prob));
  }

}

void SampleMass::quarkExchangeProcess2(int id1, int id2,double fac)
{
  vector<pair<int,int> > qpair;

  int ifla1, iflb1, ifla2,iflb2;
  hadronContent->findFlavor(id1,ifla1,iflb1);
  hadronContent->findFlavor(id2,ifla2,iflb2);
  array<int,2> iq1={ifla1,iflb1};
  array<int,2> iq2={ifla2,iflb2};
  bool isDiq1 = (abs(ifla1) < 10 && abs(iflb1) <10) ? false: true;
  bool isDiq2 = (abs(ifla2) < 10 && abs(iflb2) <10) ? false: true;

  // exclude anti-baryon-baryon collision.
  if(isDiq1 && isDiq2 && ifla1*ifla2 < 0) return;

  // quark exchange process.
  for(int i=0; i<2;i++) {
    int iqa=abs(iq1[i]);
    for(int j=0; j<2;j++) {
      int iqb=abs(iq2[j]);
      if(iq1[i]==iq2[j]) continue;
      if(iqa>10 || iqb>10) continue;
      //if(iqa>10 && iqb>10) continue;
      //if(iqa<10 && iqb>10 && iq1[i]*iq2[j] >  0) continue;
      //if(iqa>10 && iqb<10 && iq1[i]*iq2[j] >  0) continue;
      if(iqa<10 && iqb<10 && iq1[i]*iq2[j] <  0) continue;
      qpair.push_back(make_pair(i,j));
    }
  }

  int nqpair=qpair.size();
  if(nqpair==0) return;

  // chose exchange quark pair randomly
  int iq=0;
  if(nqpair==1) {
    int i1 = abs(qpair[0].first);
    int i2 = abs(qpair[0].second);
    if(i1==3 || i2 ==3 ) fac = 0.1;
    if(i1>10 || i2> 10 ) fac = 0.01;
  } else if(nqpair==2) {
    double prob[2]={1.0,1.0};
    for(int i=0;i<2;++i) {
      int i1 = abs(qpair[i].first);
      int i2 = abs(qpair[i].second);
      if(i1==3 || i2 ==3 ) prob[i] = 0.2;
      if(i1>10 || i2> 10 ) prob[i] = 0.01;
    }
    double xq=rndm->flat()*(prob[0]+prob[1]);
    if(xq-prob[0] <= 0.0) iq=0;
    else iq=1;
  } else {

    cout << "quarkExchangeProcess2 qpair.size= "<< nqpair<<endl;
    cout << " iq1= "<< iq1[0] << " " << iq1[1] <<endl;
    cout << " iq2= "<< iq2[0] << " " << iq2[1] <<endl;
    exit(1);
  }

  int i1 = qpair[iq].first;
  int i2 = qpair[iq].second;
  array<int,2> kq1=iq1;
  array<int,2> kq2=iq2;
  kq1[i1]=iq2[i2];
  kq2[i2]=iq1[i1];

  // find hadron id3.
  FlavContainer flav1(kq1[0]);
  FlavContainer flav2(kq1[1]);
  int id3=0;
  int itry=0;
  do {
    id3 = flavSel->combine(flav1, flav2);
    if(itry++>100) {
      cout << "SampleMass::quarkExchangeProcess2 infinite loop? kq1[0]= "<< kq1[0]<< " kq1[1]= "<< kq1[1]<<endl;
    cout << " iq1= "<< iq1[0] << " " << iq1[1] <<endl;
    cout << " iq2= "<< iq2[0] << " " << iq2[1] <<endl;
    cout << " kq1= "<< kq1[0] << " " << kq1[1] <<endl;
    cout << " kq2= "<< kq2[0] << " " << kq2[1] <<endl;
    cout << " id1= "<< id1 << " id2= " << id2 <<endl;
    cout << " isDiq1= "<< isDiq1 << " isDiq2= "<< isDiq2 << " id1*id2= "<< id1*id2 <<endl;
      exit(1);
    }
  }while(id3==0);


  // find hadron id4.
  FlavContainer flav3(kq2[0]);
  FlavContainer flav4(kq2[1]);
  int id4=0;
  itry=0;
  do {
    id4 = flavSel->combine(flav3, flav4);
    if(itry++>100) {
      cout << "SampleMass:quarkExchangeProcess2 infinite loop? kq2[0]= "<< kq2[0]<< " kq2[1]= "<< kq2[1]<<endl;
    cout << " iq1= "<< iq1[0] << " " << iq1[1] <<endl;
    cout << " iq2= "<< iq2[0] << " " << iq2[1] <<endl;
    cout << " kq1= "<< kq1[0] << " " << kq1[1] <<endl;
    cout << " kq2= "<< kq2[0] << " " << kq2[1] <<endl;
    cout << " id1= "<< id1 << " id2= " << id2 <<endl;
    cout << " isDiq1= "<< isDiq1 << " isDiq2= "<< isDiq2 << " id1*id2= "<< id1*id2 <<endl;
      exit(1);
    }
  } while(id4==0);

  // find hadrons.
  ParticleDataEntry* p3 = jamTable->find(id3);
  ParticleDataEntry* p4 = jamTable->find(id4);

  bool error=false;
  if(p3==0 ) {
    error=true;
    cout << "2 quark exchange wrong id3 id3= "<< id3 << " p3= "<< p3<<endl;
  }
  if(p4==0 ) {
    error=true;
    cout << "2 quark exchange wrong id4 id4= "<< id4 << " p4= "<< p4<<endl;
  }
  if(p3 !=0 && p4 != 0) {
    int ch0=pout[0]->chargeType(Id[0]) + pout[1]->chargeType(Id[1]);
    int ch1=p3->chargeType(id3) + p4->chargeType(id4);
    if(ch0 != ch1) {
      cout << "quarkExchangeProcess2::quark exchange charge does not conserve ch0= "<< ch0 << " ch1= "<< ch1<<endl;
      error=true;
    }
  }

  if(error) {
    cout << " i1= "<< i1 << " i2= "<<i2<<endl;
    cout << " iq1= "<< iq1[0] << " " << iq1[1] <<endl;
    cout << " iq2= "<< iq2[0] << " " << iq2[1] <<endl;
    cout << " kq1= "<< kq1[0] << " " << kq1[1] <<endl;
    cout << " kq2= "<< kq2[0] << " " << kq2[1] <<endl;
    cout << "id1= "<< Id[0] << " id2=  "<< Id[1] << " id3= "<< id3 << " id4= "<< id4<<endl;
   exit(1);
  }

  int iz3=p3->chargeType(id3)/3;
  int iz4=p4->chargeType(id4)/3;
  int pd3=jamTable->pid(id3);
  int pd4=jamTable->pid(id4);

  vector<ParticleDataEntry*> pcount3=jamexpa(p3,id3,pd3,iz3);
  vector<ParticleDataEntry*> pcount4=jamexpa(p4,id4,pd4,iz4);

  int n3 = pcount3.size();
  int n4 = pcount4.size();
  int kfsign1 = id3 > 0 ? 1 : -1;
  int kfsign2 = id4 > 0 ? 1 : -1;
  for(int ik1=0;ik1<n3;ik1++) 
    for(int ik2=0;ik2<n4;ik2++) {
      double prob = computeProb(pcount3[ik1],pcount4[ik2]);
      int ido1 = (pcount3[ik1]->id())*kfsign1;
      int ido2 = (pcount4[ik2]->id())*kfsign2;
      collPair.push_back(CollPair(pcount3[ik1],pcount4[ik2],ido1,ido2,fac*prob));
  }

}

void SampleMass::quarkAnnProcess(int id1,int id2,double fac)
{
  std::array<int,3> iq1, iq2;
  quarkContent(id1,iq1);
  quarkContent(id2,iq2);

  vector<pair<int,int> > qpair;
  qpair.clear();

  // find q-qbar pair.
  for(int i=0; i<3;i++) {
    if(iq1[i]==0) continue;
    for(int j=0; j<3;j++) {
      if(iq2[j]==0) continue;
      if(iq1[i]+iq2[j]==0) {
	qpair.push_back(make_pair(i,j));
      }
    }
  }

  if(qpair.size()==0) return;

  double probS=0.2;
  int i0=rndm->flat()*qpair.size();
  int i1 = qpair[i0].first;
  int i2 = qpair[i0].second;
  int qann=1;
  if(abs(iq1[i1])==1) {
    if(rndm->flat()< probS) qann=3;
    else qann=2;
  } else if(abs(iq1[i1])==2) {
    if(rndm->flat()< probS)  qann=3;
    else qann=1;
  } else if(abs(iq1[i1])==3) {
    if(rndm->flat()< 0.5) qann=1;
    else qann=2;
  } else {
    double x = rndm->flat();
    if(x< probS)  qann=3;
    else if(x< 0.5*(probS+1.0)) qann=1;
    else qann=2;
  }

  int j1=0, j2=0;
  int jq1[2]={0,0};
  int jq2[2]={0,0};
  int q1=0, q2=0;
  for(int i=0;i<3;i++) {
    if(i==i1) {
      q1 = iq1[i1] > 0 ? qann : -qann;
    } else if(iq1[i] !=0) {
      jq1[j1++]=iq1[i];
    }
    if(i==i2) {
      q2 = iq2[i2] > 0 ? qann : -qann;
    } else if(iq2[i] !=0) {
      jq2[j2++]=iq2[i];
    }
  }

  // make diquark
  int qq1=0, qq2=0;
  if(iq1[0] !=0) {
    int spin=3;
    if(jq1[0] != jq1[1])  if(rndm->flat() <0.5) spin=1;
    qq1= 1000*max(abs(jq1[0]) ,abs(jq1[1]))+ 100*min(abs(jq1[0]),abs(jq1[1]))+spin;
    if(Id[0] <0) qq1 *= -1;
  } else {
    qq1=jq1[0];
  }

  // find hadron id
  FlavContainer flav1(q1);
  FlavContainer flav2(qq1);
  int id3=0;
  do {
    id3 = flavSel->combine(flav1, flav2);
  }while(id3==0);

  // same for particle 2.
  if(iq2[0] !=0) {
    int spin=3;
    if(jq2[0] != jq2[1])  if(rndm->flat() <0.5) spin=1;
    qq2= 1000*max(abs(jq2[0]) ,abs(jq2[1]))+ 100*min(abs(jq2[0]),abs(jq2[1]))+spin;
    if(Id[1]<0) qq2 *= -1;
  } else {
    qq2=jq2[0];
  }

  FlavContainer flav3(q2);
  FlavContainer flav4(qq2);
  int id4=0;
  do {
    id4 = flavSel->combine(flav3, flav4);
  } while(id4==0);

  ParticleDataEntry* p3 = jamTable->find(id3);
  ParticleDataEntry* p4 = jamTable->find(id4);

  bool error=false;
  if(p3==0 ) {
    error=true;
    cout << " wrong id3 id3= "<< id3 << " p3= "<< p3<<endl;
  }
  if(p4==0 ) {
    error=true;
    cout << " wrong id4 id4= "<< id4 << " p4= "<< p4<<endl;
  }
  if(p3 !=0 && p4 != 0) {
    int ch0=pout[0]->chargeType(Id[0]) + pout[1]->chargeType(Id[1]);
    int ch1=p3->chargeType(id3) + p4->chargeType(id4);
    if(ch0 != ch1) {
      cout << " charge does not conserve ch0= "<< ch0 << " ch1= "<< ch1<<endl;
      error=true;
    }
  }

  if(error) {
    cout << " i1= "<< i1 << " i2= "<<i2<<endl;
    cout << " iq1= "<< iq1[0] << " " << iq1[1] << "  "<< iq1[2]  <<endl;
    cout << " iq2= "<< iq2[0] << " " << iq2[1] << "  "<< iq2[2]  <<endl;
    cout << " ann= "<< qann <<endl;;
    cout << "q1 = "<< q1 << " jq1= "<< jq1[0] << " "<< jq1[1] << "  qq1= "<< qq1 <<endl;
    cout << "q2 = "<< q2 << " jq2= "<< jq2[0] << " "<< jq2[1] << "  qq2= "<< qq2 <<endl;
    cout << "id1= "<< Id[0] << " id2=  "<< Id[1] << " id3= "<< id3 << " id4= "<< id4<<endl;
   exit(1);
  }

  int iz3=p3->chargeType(id3)/3;
  int iz4=p4->chargeType(id4)/3;
  int pd3=jamTable->pid(id3);
  int pd4=jamTable->pid(id4);

  vector<ParticleDataEntry*> pcount3=jamexpa(p3,id3,pd3,iz3);
  vector<ParticleDataEntry*> pcount4=jamexpa(p4,id4,pd4,iz4);

  int n3 = pcount3.size();
  int n4 = pcount4.size();
  int kfsign1 = id3 > 0 ? 1 : -1;
  int kfsign2 = id4 > 0 ? 1 : -1;
  for(int ik1=0;ik1<n3;ik1++) 
    for(int ik2=0;ik2<n4;ik2++) {
      double prob = computeProb(pcount3[ik1],pcount4[ik2]);
      if(prob>0.0) {
        int ido1 = (pcount3[ik1]->id())*kfsign1;
        int ido2 = (pcount4[ik2]->id())*kfsign2;
        collPair.push_back(CollPair(pcount3[ik1],pcount4[ik2],ido1,ido2,fac*prob));
      }
  }

}

double SampleMass::probBW1(double srt, ParticleDataEntry* p1, ParticleDataEntry* p2)
{
  if(optProb==0) {
    double em1=p1->m0();     // pole mass
    double em2=p2->m0();     // pole mass
    if(srt > em1 + em2 + eKinMin) return PCM(srt,em1,em2);
    return 0.0;
  }

  double em0[2], gam0[2],emin[2],emax[2], bw[2], ymin[2],ymax[2];
  int iex[2];
  ParticleDataEntry* pv[2];
  pv[0]=p1;
  pv[1]=p2;

  for(int jt=0; jt<2; jt++) {
    em0[jt]=pv[jt]->m0();     // pole mass
    gam0[jt]=pv[jt]->mWidth();
    emin[jt]=em0[jt];
    emax[jt]=emin[jt];
    iex[jt]=0;
    if(gam0[jt]> 1e-5) {
      iex[jt]=1;
      //emin[jt] -= gam0[jt];
      emin[jt]=pv[jt]->mMin();
      if(emin[jt]==0) {
	cout << "emin= "<< emin<< " id= "<< pv[jt]->id()
	    <<endl;
	exit(1);
      } 
    }
  }

  if(srt < emin[0] + emin[1] + eKinMin) return 0.0;

  for(int jt=0; jt<2; jt++)
  if(iex[jt] ==  1) {
    emax[jt]=max(emin[jt],srt-emin[1-jt]-eKinMin);
  }
  if(emin[0] > emax[0] || emin[1] > emax[1]) return 0.0;

  double pf0=0.0;
  pf0=PCM(srt,emin[0],emin[1]);
  //if(srt > em0[0]+em0[1]) pf0=PCM(srt,em0[0],em0[1]);
  //else return false;

  for(int jt=0; jt<2; jt++) {
    bw[jt]=1.0;
    if(iex[jt] >=1) {
      ymin[jt]=atan((pow2(emin[jt])-pow2(em0[jt]))/(em0[jt]*gam0[jt]));
      ymax[jt]=atan((pow2(emax[jt])-pow2(em0[jt]))/(em0[jt]*gam0[jt]));
      bw[jt]=(ymax[jt]-ymin[jt])/M_PI;
    }
  }
  if(bw[0]<=0 || bw[1]<=0) return 0.0;

  return bw[0]*bw[1]*pf0;
}

// Use table.
double SampleMass::probBW3(double srt,ParticleDataEntry* pv1,ParticleDataEntry* pv2)
{
  double gam1=pv1->mWidth();
  double gam2=pv2->mWidth();
  if(gam1<1e-5 && gam2 <1e-5) {
    double em1=pv1->m0();     // pole mass
    double em2=pv2->m0();     // pole mass
    if(srt > em1 + em2 + eKinMin) return PCM(srt,em1,em2);
    return 0.0;

  } else if (gam1 >= 1e-5 && gam2  < 1e-5) {
    return bwint(pv1, srt, pv2->m0());

  } else if (gam1 < 1e-5 && gam2  >= 1e-5) {
    return bwint(pv2, srt, pv1->m0());

  } else {
    ParticleDataEntry* p1 = isoParticle[pv1];
    ParticleDataEntry* p2 = isoParticle[pv2];
    if(p1==0 || p2==0) {
      cout << " pv1= "<< pv1->id() << " pv2= "<< pv2->id()<<endl;
      exit(1);
    }
    return BWint[make_pair(p1,p2)]->getBW(srt);
  }

}

double SampleMass::probBW2(double srt,ParticleDataEntry* pv1,ParticleDataEntry* pv2)
{
  double gam1=pv1->mWidth();
  double gam2=pv2->mWidth();
  if(gam1<1e-5 && gam2 <1e-5) {
    double em1=pv1->m0();     // pole mass
    double em2=pv2->m0();     // pole mass
    if(srt > em1 + em2 + eKinMin) return PCM(srt,em1,em2);
    return 0.0;
  } else if (gam1 >= 1e-5 && gam2  < 1e-5) {
    return bwint(pv1, srt, pv2->m0());
  } else if (gam1 < 1e-5 && gam2  >= 1e-5) {
    return bwint(pv2, srt, pv1->m0());
  } else {
    //cout << "SampleMass:probBW2 pv1= "<< pv1->name()
//	<< " pv2= "<< pv2->name() <<endl;
    return bwint2(pv1,pv2,srt);
  }

}

double SampleMass::BreitWigner(ParticleDataEntry* p,double emd)
{
  double em0=p->m0();
  double gam=p->mWidth();
  //if(optWidth >0) gam = decay->getTotalWidth(p,emd);
  if(optWidth >0) gam = jamTable->totalWidth(p,emd);
  double gamr=gam;

  return 2.0/M_PI*emd*emd*gamr/(pow2(emd*emd-em0*em0)+pow2(emd*gam));
}

// Integration of Breit-Wigner 
double SampleMass::bwint(ParticleDataEntry* p, double srt,double em2)
{
  double em0=p->m0();
  double gamr=p->mWidth();
  double emin=p->mMin();
        if(srt < emin + em2 + eKinMin) return 0.0;
  double emax=srt-em2;
  //double gr=em0*gamr;
  //double ymax=2*atan((emax*emax-em0*em0)/gr);
  //double ymin=2*atan((emin*emin-em0*em0)/gr);
  double ymax=2*atan(2*(emax-em0)/gamr);
  double ymin=2*atan(2*(emin-em0)/gamr);
  double dy=ymax-ymin;
  if(dy<=0.0) return 0.0;

  double ds=0.0;
  for(int i=0;i<NGP38;i++) {
    double y=ymin+xg38f[i]*dy;
    //double emd=sqrt(tan(0.5*y)*gr+ em0*em0);
    //double dc=gr/(1+cos(y))/(2*emd);
    double emd=tan(0.5*y)*gamr/2 + em0;
    double dc= gamr/( 0.25*gamr*gamr + pow2(emd - em0) );
    double bw=BreitWigner(p,emd);
    double pf=PCM(srt,emd,em2);
    double dw=wg38f[i]*dy;
    ds += bw*pf*dw/dc;
  }
  return ds;

}

// Integration of two-Breit-Wigner 
double SampleMass::bwint2(ParticleDataEntry* p1,ParticleDataEntry* p2,double srt)
{
  double emin1=p1->mMin();
  double emin2=p2->mMin();
        if(srt < emin1 + emin2 + eKinMin) return 0.0;

  double emax1=srt-emin2;
  double emax2=srt-emin1;
  if(srt<=emin1+emin2) return 0.0;
  if(emax1<=emin1) return 0.0;
  if(emax2<=emin2) return 0.0;

  double ds=0.0;
  double em0=p1->m0();
  double gamr=p1->mWidth();
  double emax=srt-emin2;
  //double gr=em0*gamr;
  //double ymax=2*atan((emax*emax-em0*em0)/gr);
  //double ymin=2*atan((emin1*emin1-em0*em0)/gr);
  double ymax=2*atan(2*(emax-em0)/gamr);
  double ymin=2*atan(2*(emin1-em0)/gamr);

  double dy=ymax-ymin;
  if(dy<=0.0) return 0.0;

  for(int i=0;i<NGP38;i++) {
    double y=ymin+xg38f[i]*dy;
    //double emd=sqrt(tan(0.5*y)*gr+ em0*em0);
    //double dc=gr/(1+cos(y))/(2*emd);
    double emd=tan(0.5*y)*gamr/2 + em0;
    double dc= gamr/( 0.25*gamr*gamr + pow2(emd - em0) );
    double bw=BreitWigner(p1,emd);
    double bw2=bwint(p2,srt,emd);
    double dw=wg38f[i]*dy;
    ds += bw*bw2*dw/dc;
  }
  return ds;

}


}
