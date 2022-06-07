#include <jam2/interaction/ScatterKin.h>
#include <Pythia8/PythiaStdlib.h>
#include <jam2/hadrons/JamStdlib.h>
#include <Pythia8/PythiaStdlib.h>

namespace jam2 {

using namespace std;
using Pythia8::pow3;
using Pythia8::pow4;

ScatterKin::ScatterKin(Pythia8::Info* inf,Pythia8::Settings* s,JamParticleData* jd,
	CrossSection* xs, Pythia8::Rndm* r)
     : info(inf), settings(s),jamParticleData(jd), xSection(xs),rndm(r)
{

  withMeanField = settings->mode("MeanField:mode");
  optPreserveReactionPlane = false;
  inelOnly = settings->flag("Cascade:InelasticOnly");
  printColl=settings->mode("Cascade:PrintCollision");
  optConstQuarkDiffra=settings->mode("Cascade:optConstQuarkDiffra");

  // Propagation by canonical momentum in RQMD mode
  optVectorPotential=settings->mode("MeanField:optVectorPotential");
  optVdot=settings->mode("MeanField:optVdot");
  optPropagate=0;
  if(optVectorPotential==1 && optVdot==0) optPropagate=1;

  optPotential=settings->mode("MeanField:optPotential");
  if(settings->mode("MeanField:mode")==0) optPotential=0;

  isDebug=settings->mode("Check:Debug");

}

ScatterKin::~ScatterKin() 
{
}

void ScatterKin::initialize(TwoBodyInterList* inter)
{
  // Save incoming particles. Do not change values of inComing[i] itself
  // because of the possibility of cancellation of this collision.
  inComing[0] = inter->getParticle(0);
  inComing[1] = inter->getParticle(1);
  numCollision = abs(inComing[0]->getNColl())+abs(inComing[1]->getNColl())+1;

  xCM=0.0;
  pCM=0.0;
  int nIn=2;
  Vec4 pkin[2];

  // Propagate the incoming particles to the collision time and
  // save some properties.
  for(int i=0; i<nIn; i++) {
    xIncoming[i] = inComing[i]->getR();
    pIncoming[i] = inComing[i]->getP();

    mIn[i] = inComing[i]->getMass();
    mInEff[i] = inComing[i]->getEffectiveMass();
    SM[i] = mInEff[i]/mIn[i];
    sPot[i] = inComing[i]->pots();
    vPot[i] = inComing[i]->potv();


    double meff=pIncoming[i].mCalc();
    if(abs(meff-mInEff[i])>1e-3) {
	cout << "scatter:ini meff?  meff= "<< meff
	    << " mInEff= "<< mInEff[i]
	    << " meff= "<< meff
	    << " spot= "<< inComing[i]->pots()
	    << " m0= "<< mIn[i]
	    << " m*= "<< mIn[i]+inComing[i]->pots()
	    << " id= "<< inComing[i]->getID()
	    << endl;
	exit(1);
    }
    pkin[i]=pIncoming[i];

    // set kinetic momentum p* in case of canonical momentum propagation.
    if(optPropagate==1) {
      pkin[i]=inComing[i]->getPkin();
      pIncoming[i]=pkin[i];
    }

    mOut[i] = mIn[i];
    idIn[i] = inComing[i]->getID();
    idOut[i] = idIn[i];
    pidOut[i] = inComing[i]->getPID();
    pd[i]=inComing[i]->getParticleDataEntry();
    isBaryon[i]= (abs(idIn[i])/1000)%10 > 0;
    //iBaryon[i] = cpair.baryon(i);
    iBaryon[i] = inComing[i]->baryon();

    double ctime = inter->getCollisionTime(i);
    double e = pIncoming[i][0];
    double dt = max(0.0,ctime-xIncoming[i][0]);

    //if(ctime < xIncoming[i][0]) cout << "ScatterKin::ctime= "<< ctime << " t0= "<< xIncoming[i][0]<< endl;
  
    Vec4 vel=pIncoming[i]/e;
    xOut[i][1] = xIncoming[i][1] + dt*vel[1];
    xOut[i][2] = xIncoming[i][2] + dt*vel[2];
    xOut[i][3] = xIncoming[i][3] + dt*vel[3];
    xOut[i][0] = ctime;

    //pIncoming[i][0] += vPot[i][0];
    xCM += xIncoming[i];
    pCM += pIncoming[i];

    constQ[i]=inComing[i]->constQuark();
    qFactor[i]=inComing[i]->qFactor();
    if(ctime >= inComing[i]->getTf())  {
      constQ[i][0]=1;
      constQ[i][1] = abs(inComing[i]->baryon()) == 3 ? 2: 1; 
      qFactor[i]=1.0;
    }

  }

  if(printColl) {
    cout << " before scatt n= " << nIn
	  << " constQ1= " << constQ[0][0] << " " <<constQ[0][1]
	  << " constQ2= " << constQ[1][0] << " " <<constQ[1][1]
	 <<endl;
  }

    xCM /= nIn;

    MfromCM.reset();
    MfromCM.fromCMframe(pIncoming[0],pIncoming[1]);

    preHadronA = xOut[0][0] < inComing[0]->getTf() ? true : false;
    preHadronB = xOut[1][0] < inComing[1]->getTf() ? true : false;

    //MtoCM = MfromCM;
    //MtoCM.invert();

    sCM = pCM.m2Calc();
    if(sCM<=0.0) {
	cout << "ScatterKin: s<0? " << sCM << endl;
	exit(1);
    }
    eCM = sqrt( max(0.0, sCM) );
    beX=pCM[1]/pCM[0];
    beY=pCM[2]/pCM[0];
    beZ=pCM[3]/pCM[0];

        if(eCM < mInEff[0]+mInEff[1]) {
          cout << "ScatterKin::ini eCM = " << eCM
	   << " m1+m2= "<< mInEff[0]+mInEff[1]
	   << " mInEff0= "<< mInEff[0]
	   << " mInEff1= "<< mInEff[1]
	   <<endl;
	  exit(1);
      }

    Vec4 p = pIncoming[0];
    p.bst(-beX,-beY,-beZ);
    thetaCM = p.theta();
    phiCM   = p.phi();

    // pcm in the cpair may be different because of QMD.
    //pIni = cpair.getCMmomentum();
    pIni = p.pAbs();

  MfromCMcan = MfromCM;
  pInican = pIni;
  eCMcan=eCM;
  eCM1=sqrt(mInEff[0]*mInEff[0] + pIni*pIni);
  eCM2=sqrt(mInEff[1]*mInEff[1] + pIni*pIni);

  if(optPropagate==1) {
    MfromCMcan.reset();
    MfromCMcan.fromCMframe(inComing[0]->getMomentum(),inComing[1]->getMomentum());
    eCMcan=Pythia8::m(inComing[0]->getMomentum(),inComing[1]->getMomentum());
    //pInican = PCM(eCM,mInEff[0],mInEff[1]);
    pInican = PCM(eCMcan,mInEff[0],mInEff[1]);
    eCM1=sqrt(mInEff[0]*mInEff[0] + pInican*pInican);
    eCM2=sqrt(mInEff[1]*mInEff[1] + pInican*pInican);
  }

  tFormA = inComing[0]->getTf();
  tFormB = inComing[1]->getTf();

}

int ScatterKin::selectChannel(TwoBodyInterList* inter, int ncol)
{
  cpair=inter->getCpair();
  nOperation = ncol;
  initialize(inter);

    // First check elastic scattering or not.
    double xsig;
    double sig = cpair.getSigma();
    double sigel = cpair.getSigmaElastic();
    if(inelOnly) {
	sig=sig-sigel;
	if(sig <1e-8) {
	    cout << "ScatterKin::selectChannel  There is no inelastic channel do not use InelOnly mode"
		<< endl;
	    cout << " sig= " << cpair.getSigma()
		 << " sigel= " << cpair.getSigmaElastic()
		 << " srt= "<< eCM
	<< " id1= " << cpair.getID(0)
	<< " id2= " << cpair.getID(1)
		 << endl;
	    //exit(1);
	}
	xsig=sig*rndm->flat();
    } else { 
	xsig = sig*rndm->flat() - sigel;
    }
    cpair.setXsig(xsig);
    channel=0;
    nOutParticle=2;
    if(xsig < 0.0) return ELASTIC;

    // Check inelastic channels.

    double srt = eCM - (mInEff[0] -mIn[0]) - (mInEff[1]-mIn[1]);
    cpair.setECM(srt);
    cpair.setPCM( PCM(srt,mIn[0],mIn[1]) );

    cpair.setMode(3);
    xSection->sigma(cpair);
    channel=cpair.getChannel();
    if(channel== -1) return ELASTIC;
    OutGoing out = cpair.getOutGoing();

    // unformed meson scatter only elastically.
    //if(preHadronA && inComing[0]->baryon()==0)  return ELASTIC;
    //if(preHadronB && inComing[1]->baryon()==0)  return ELASTIC;
    //if(preHadronA )  return ELASTIC;
    //if(preHadronB )  return ELASTIC;

    typeDiff=1;

    /*
    if(inComing[0]->qFactor() <1.0 && inComing[1]->qFactor() <1.0) {
	//return ELASTIC;
	if(eCM < 4.0) return ELASTIC;
	typeDiff=3;
	if(rndm->flat()<0.5) typeDiff=4;
	return DIFFRACTIVE;
    }
    if(inComing[0]->qFactor() <1.0 || inComing[1]->qFactor() <1.0) {
	if(eCM < 4.0) return ELASTIC;
	typeDiff=3;
	if(inComing[0]->qFactor() <1.0) typeDiff=4;
        return DIFFRACTIVE;
    }
    */

    if(printColl) {
	cout << "*** ScatterKin::selectChannel= " << channel << endl;
        cout << " icltyp= " << cpair.getCollType()
	    << " channel= "<< channel
	<< " ecm= " << eCM
	<< " id1= " << cpair.getID(0)
	<< " id2= " << cpair.getID(1)
	<< " out1= " << out.id1
	<< " out2= " <<  out.id2
	<< " out3= " << out.id3
	<< " out4= " <<  out.id4
	<< " mout1 = " << out.m1
	<< " mout2= " <<  out.m2
	<<endl;
    }


  if (channel== -2) { // s-channel string formation.
    cpair.setChannel(0);
    //out = cpair.getOutGoing();
    idOut[0]=cpair.getOutGoing(0).id1;
    return SOFT;
  }


    // t-channel string formation.
  if((out.id1==92 && out.id2==92) || (out.id1==-92 && out.id2==-92)) {
    channel=SOFT;


    if(optConstQuarkDiffra >0 && optConstQuarkDiffra <5) {
      if(preHadronA && preHadronB) {return channel=ELASTIC;}
      if(preHadronA || preHadronB) channel=DIFFRACTIVE;
    }

        //if(qFactor[0] <1.0 &&  qFactor[1] <1.0) return ELASTIC;

	/*
        if(qFactor[0] <1.0 ||  qFactor[1] <1.0) {
          typeDiff=5;
	  double probbb=0.0;
          if(rndm->flat() < probbb) {
            if(qFactor[0] <1.0 &&  qFactor[1] == 1.0) typeDiff=4;
	    else if(qFactor[0] == 1.0 &&  qFactor[1] < 1.0) typeDiff=3;
	    else {
              typeDiff=4;
              if(rndm->flat() <0.5) typeDiff=3;
	    }
	  }
          return DIFFRACTIVE;
        }
	*/

	// test
	//channel=ELASTIC;

    } else if(out.id1==93) {
	channel=SOFT;

    } else if (channel>= 0) {


	if(out.id2 == 0) {

	  channel=ABSORPTION;
	  idOut[0]=out.id1;
	  mOut[0]=out.m1;

	//channel=ELASTIC;

	} else {

        //if(qFactor[0] <1.0  ||  qFactor[1] <1.0) return ELASTIC;

	channel=RESONANCE;


	// Pythia only.
	//return SOFT;

	// reduce resonance production prob. 
	//if(rndm->flat() < 0.20) return SOFT;

        //pd[0]=out.pa1; pd[1]=out.pa2;
        //pd[2]=out.pa3; pd[3]=out.pa4;

	/*
	int colltype=cpair.getCollType();
	if(colltype==2 && inComing[0]->baryon() !=0) {
	    mOut[0]=out.m2;
	    mOut[1]=out.m1;
	    idOut[0]=out.id2;
	    idOut[1]=out.id1;
	} else {
	*/

        double sm1 = mInEff[0]/mIn[0];
        double sm2 = mInEff[1]/mIn[1];
	mOut[0]=out.m1;
	mOut[1]=out.m2;
	double mtot = sm1*mOut[0] + sm2*mOut[1];
	idOut[0]=out.id1;
	idOut[1]=out.id2;
	idOut[2]=out.id3;
	idOut[3]=out.id4;
	mOut[2]=0.0;
	mOut[3]=0.0;

	nOutParticle=2;
	if(out.id3 !=0) {
	    nOutParticle=3;
	    xOut[2] = xCM;
	    idOut[2]=out.id3;
	    mOut[2]=out.m3;
	    mtot += mOut[2];
	}
	if(idOut[3] !=0) {
	    nOutParticle=4;
	    xOut[3] = xCM;
	    idOut[3]=out.id4;
	    mOut[3]=out.m4;
	    mtot += mOut[3];
	}
	if(eCM < mtot) {
            mOut[0] = mIn[0];
            mOut[1] = mIn[1];
	    return ELASTIC;
	}

	}

	if(idOut[0]==0 || idOut[1]==0) {

	    cout << "ScatterKin::selectChannel id of outgoing channel? id1=" << idOut[0]
		<< " id2= " << idOut[1]
		<<endl;
	    exit(1);
	}

        if(printColl) {
	    cout << " nOutParticle= "<< nOutParticle<<endl;
	}
    }

    return channel;
}

void ScatterKin::checkMomentum(vector<EventParticle*>& outgoing)
{
    Vec4 ptot=0.0;
    int iz0 = inComing[0]->charge() + inComing[1]->charge();
    int is0 = inComing[0]->strange() + inComing[1]->strange();
    int ib0 = inComing[0]->baryon() + inComing[1]->baryon();
    int nq0 = constQ[0][0]+constQ[0][1]+ constQ[1][0]+constQ[1][1];
    if(channel==ABSORPTION) nq0 -= 2;

    //double tmin= std::min(inComing[0]->getT(), inComing[1]->getT());

    int error=0;
    int iz=0, is=0, ib=0;
    int nq=0;
    for(int i=0; i<(int)outgoing.size(); i++) {
      int id= outgoing[i]->getID();
      double mass = outgoing[i]->getMass();
	ptot += outgoing[i]->getMomentum();
	iz += outgoing[i]->charge();
	is += outgoing[i]->strange();
	ib += outgoing[i]->baryon();

	/*
	if(outgoing[i]->getT()< tmin) {
	  cout <<" time is earler ? " << id 
	    << " tmin= " << tmin
	    << outgoing[i]->getR();
	  error=31;
	}
	*/

	if(isDelta(id) && mass <1.0764)  {
	    cout << " delta mass? " << mass
		<< " id= " << id
		<<endl;
	    error=10;
	}

	if(abs(id) == 313 && mass < 0.6) {
	    cout << " K* mass strange m= "<< mass << " id= " << id<<endl;
	    exit(1);
	}


	if((id !=22 && abs(id) !=11)  && outgoing[i]->getMass() <0.1)  {
	    cout << "hadron mass ? " << outgoing[i]->getMass()
		<< " id= " << id
		<<endl;
	    error=11;
	}


        int* cq = outgoing[i]->constQuark();
	nq += cq[0]+cq[1];
	if(cq[0] !=0 && cq[0] !=1) {
	    cout << " cq1 strange id= " <<outgoing[i]->getID()
		<< " cq0= " << cq[0] << " cq1= " << cq[1]
		<<endl;
	    error=12;
	}
	if(cq[1] !=0 && cq[1] !=1 && cq[1] !=2) {
	    cout << " cq2 strange id= " << outgoing[i]->getID()
		<< " cq0= " << cq[0] << " cq1= " << cq[1]
		<<endl;
	    error=13;
	}
	if(outgoing[i]->baryon() ==0 && cq[1]==2) error=14;

	if(id > 22) {
        double meff=outgoing[i]->getMomentum().mCalc();
        if(abs(meff-outgoing[i]->getEffectiveMass())>1e-3) {
	cout << "Scatter2::checkMom meff?  meff= "<< meff
	    << " mOutEff= "<< outgoing[i]->getEffectiveMass()
	    << " id= " << id
	    << endl;
	for(int j=0; j<(int)outgoing.size(); j++) outgoing[j]->print(cout);
	exit(1);
        }
	}

	/*
	if(outgoing[i]->getT() < 0.0) {
	  cout << "t is negative?  id= "<< id << outgoing[i]->getR();
	  cout << "pcm = " << pCM << " eCM= " << eCM
	    << " channel=" << channel
	    << endl;
	cout << "incomoing kf1 = " << inComing[0]->getID()
	     << " kf2 = " << inComing[1]->getID() << endl;
	cout << " size outgoing= " <<  outgoing.size() << endl;
	  error=14;
	}
	*/

	//if(outgoing[i]->lifetime() < outgoing[i]->getTf()) error=15;

	/*
	if(outgoing[i]->getID() ==225 && outgoing[i]->getMass() < 0.5)  {
	    cout << " after scatter m = "<< outgoing[i]->getMass()
		<<endl;
	    error=20;
	}
	*/
    }

//...Check energy momentum conservation.
    Vec4 diff= ptot-pCM;
    //double ediff = abs(diff.mCalc());
    double pdiff = abs(diff.pAbs2());
    //if((optPropagate==0 && ediff/eCM > 1e-5) || (optPropagate==0 && pdiff >1e-5)) {
    if(optPropagate==0 && pdiff >1e-5) {
	cout << "(Scatter2::) Energy Momentum not conserved"
	    << " channel=" << channel
            << " id1= " << cpair.getOutGoing().id1
            << " id2= " <<  cpair.getOutGoing().id1 
           << endl;
	cout << " kf1 = " << inComing[0]->getID()
	     << " kf2 = " << inComing[1]->getID() << endl;

	cout << "pcm = " << pCM << " eCM= " << eCM << endl;
	cout << "ptot= " << ptot << endl;
	for(int i=0; i<(int)outgoing.size(); i++)
	    cout << "pf = " << outgoing[i]->getP()<< endl;
	cout << endl;
        cout << "diff= " << diff << " error = " << pdiff / eCM << endl;

     error=2;
    }
    

    if(channel!=HARD && channel != BBANIHILATION && channel != ABSORPTION  && nq0 != nq) {
	cout << " const.quark not conserved after scatter nq0= " << nq0
	    << " nq= " << nq
	    <<endl;
	//cout << "1) q1= " << constQ1[0] << " q2= " << constQ1[1]<<endl;
	//cout << "2) q1= " << constQ2[0] << " q2= " << constQ2[1]<<endl;
	error=5;
    }

    if(iz0 != iz) {
	cout << " charge not conserved after scatter iz0= " << iz0
	    << " iz= " << iz
	    <<endl;
	error=1;
    }
    if(ib0 != ib) {
	cout << "baryon number not conserved after scatter ib0= " << ib0
	    << " ib= " << ib
	    <<endl;
	error=3;
    }
    if(is0 != is) {
	cout << "strangeness not conserved after scatter is0= " << is0
	    << " is= " << is
	    <<endl;
	error=4;
    }

    if(error !=0) {
	cout << "After scatter Erorr! pcm = " << pCM << " eCM= " << eCM
	    << " channel=" << channel
	    << " error=" << error
	    << endl;
	cout << "incomoing kf1 = " << inComing[0]->getID()
	     << " kf2 = " << inComing[1]->getID() << endl;
	cout << " size outgoing= " <<  outgoing.size() << endl;
	for(int i=0; i<(int)outgoing.size(); i++) {
	    outgoing[i]->print(cout);
	}
	if( printColl ) exit(1);

	//if(error==5) ;//cin.get();
	//else
	//if(error>=20) exit(0);
    }


    /*
  int pid1=inComing[0]->getPID();
  int pid2=inComing[1]->getPID();
  if(pid1 != id_nucl || pid2 != id_nucl) {
    if(outgoing.size()==2) {
      int pid3=outgoing[0]->getPID();
      int pid4=outgoing[1]->getPID();
      if(pid3==id_nucl && pid4==id_nucl) {
	  cout << " collision of RR -> NN! "
	      << "id1= "<< inComing[0]->getID()
	      << " id2= "<< inComing[1]->getID()
	      << " -> id3= "<< outgoing[0]->getID()
	      << " id4= "<< outgoing[1]->getID()
	      <<endl;
	      cin.get();
      }
    }
  }
  */


}

} // end namespace jam2.
