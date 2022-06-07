#include <jam2/interaction/JPythia.h>
#include <jam2/hadrons/JamStdlib.h>

namespace jam2 {

using namespace std;

JPythia::JPythia(string xmlDir, bool printBanner)
{
  pythia = new Pythia(xmlDir,printBanner);
  particleData = &pythia->particleData;
  settings = &pythia->settings;
  init();
}

JPythia::JPythia(Settings& settingsIn, ParticleData& particleDataIn,
  bool printBanner)
{
  pythia = new Pythia(settingsIn,particleDataIn,printBanner);
  particleData = &particleDataIn;
  settings = &pythia->settings;
  init();
}

void JPythia::init()
{
  flavSel=0;
  pdfHadron=0;
  aveConstFormationTime=0;
  aveFormationTime=0;
  nCFormationTime=0;
  nFormationTime=0;
  allowRescatterSameString=settings->mode("Cascade:allowRescatterSameString");
}

void JPythia::init0(JamParticleData* jd, HadronDecay* dec)
{
  isDebug=pythia->settings.mode("Check:Debug");
  jamParticleData=jd;
  hDecay=dec;

  // need for beamA and beamB
  flavSel = new StringFlav();
  flavSel->init(pythia->settings,  particleData, &pythia->rndm, &pythia->info);

  // making PDF any hadrons
  pdfProton = pythia->getPDFPtr(2212);
  pdfPion  = pythia->getPDFPtr(211);
  pdfHadron = new PDFHadron();
  pdfHadron->setPDF( pdfProton, pdfPion);
  pythia->setPDFPtr( pdfHadron, pdfHadron);

  optConstQscatt=pythia->settings.mode("Cascade:optConstQuarkScattHard");
  optConstQuarkDiffra=pythia->settings.mode("Cascade:optConstQuarkDiffra");
  optDiq=0;

  userHook=settings->flag("Cascade:constQuarkDiffractive");
  //if(userhook) addUserHooksPtr(&selectSASD);
  if(userHook) pythia->setUserHooksPtr(&selectSASD);

  pythia->init();

}


JPythia::~JPythia()
{
  cout << "JPythia::ave const. formation time= "<< constFormatinotime()
       << " formation time= " << formatinotime() 
       <<endl;

  if(flavSel !=0) delete flavSel;
  if(pdfHadron !=0) delete pdfHadron;
  delete pythia;
}

void JPythia::initBeam(const int id[2], const double m[2], const double m0[2],Vec4& pA, Vec4& pB,
	bool prehadron1, bool prehadron2)
{
    preHadronA=prehadron1;
    preHadronB=prehadron2;
    int idA=convertID(id[0],isBaryonA);
    int idB=convertID(id[1],isBaryonB);

    // mA and mB are the effective mass in case potential is on.
    // m0[2] are the rest mass.
    mA = m[0];
    mB = m[1];
    eA = sqrt( pA.pAbs2() + mA*mA);
    eB = sqrt( pB.pAbs2() + mB*mB);
    eCM = (pA + pB).mCalc();

    if(mA > 4.0 || mB > 4.0) {
	cout << "JPythia::initBeam mass of hadrons too large ? mA= "<< mA
	    << " mB= "<< mB
	    <<endl;
    }


    // Find boost+rotation needed to move from/to CM frame.
    MfromCM.reset();
    MfromCM.fromCMframe( pA, pB);
    MtoCM = MfromCM;
    MtoCM.invert();

    double pzAcm    = 0.5 * sqrt( (eCM + mA + mB) * (eCM - mA - mB)
           * (eCM - mA + mB) * (eCM + mA - mB) ) / eCM;
    double pzBcm    = -pzAcm;

   // This is forbidden because JPythia is not a friend of Infor....
  //info.setBeamA( idA, pzAcm, eA, m[0]);
  //info.setBeamB( idB, pzBcm, eB, m[1]);
  //info.setECM( eCM);

    // So this function is added in Pythia.h
    //setBeam(eCM,idA,idB,mA,mB);
    pythia->setBeam(eCM,idA,idB,m0[0],m0[1]);

    // and call nextKinematics() to set new eCM.
    //beamShape->set(pzAcm,pzBcm);

    /*  We do not treat lepton beam.
    int idAabs = abs(idA);
    int idBabs = abs(idB);
    bool isLeptonA    = (idAabs > 10 && idAabs < 17);
    bool isLeptonB    = (idBabs > 10 && idBabs < 17);
    bool isUnresLep   = !settings.flag("PDF:lepton");
    bool isUnresolvedA = ( isLeptonA && (idAabs%2 == 0 || isUnresLep) );
    bool isUnresolvedB = ( isLeptonB && (idBabs%2 == 0 || isUnresLep) );
    */
    bool isUnresolvedA = false;
    bool isUnresolvedB = false;

    /*
    PDF* pdfAPtr=pdfBaryon;
    PDF* pdfBPtr=pdfBaryon;
    if ((abs(idA)/1000)%10 == 0) pdfAPtr = pdfMeson;
    if ((abs(idB)/1000)%10 == 0) pdfBPtr = pdfMeson;
    */

    PDF* pdfAPtr=pdfHadron;
    PDF* pdfBPtr=pdfHadron;
    Pythia8::PDF  *pdfHardAPtr = pdfAPtr;
    Pythia8::PDF  *pdfHardBPtr = pdfBPtr;

    pythia->beamA.init(idA, pzAcm, eA, m0[0], &pythia->info,pythia->settings,particleData,
	    &pythia->rndm,pdfAPtr,pdfHardAPtr,isUnresolvedA, flavSel);
    pythia->beamB.init(idB, pzBcm, eB, m0[1], &pythia->info,pythia->settings,particleData,
	    &pythia->rndm,pdfBPtr,pdfHardBPtr,isUnresolvedB, flavSel);

    if(isDebug) {
    cout << " JPythia::initBeam before MfromCM  id1= " << id[0]
	<< " id2= " << id[1]
	<<endl;

    cout << " p1= " << pA
	<<endl
	<< " p2= " << pB
	<< " pCM= " << pA+pB
	<< " eCM= " << (pA+pB).mCalc()
	<<endl;

    cout << "idA= " << idA
	<< " idB= " << idB
	<< " ecm= " << eCM
	<<endl;
    cout << "After init beamA= "<< pythia->beamA.size()
         << " beamB= "<< pythia->beamB.size()
	 <<endl;
    }

}

// some part of Pythia does not work for ID like 12212 N(1440),
// so we convert e.g. 12212 -> 2212.
int JPythia::convertID(int id, int& isBaryon)
{
    int idBeamAbs=abs(id);
    int idnew=id;
    int id0 = (idBeamAbs)%10;
    if ((idBeamAbs/1000)%10 == 0) {
	int id1 = (idBeamAbs/100)%10;
	int id2 = (idBeamAbs/10)%10;
	idnew=id1*100+id2*10+id0;
	isBaryon=0;

    } else {
	isBaryon=1;
	int id1 = (idBeamAbs/1000)%10;
	int id2 = (idBeamAbs/100)%10;
	int id3 = (idBeamAbs/10)%10;
	int id4 = id1*1000+id2*100+id3*10;
	// PDG ID of 1214 does not work for diffractive scattering in pythia.
	if(id4 == 1210) id4=2110;
	if(id4 == 2120) id4=2210;
	idnew=id4+id0;
	//idnew=id1*1000+id2*100+id3*10+id0;
    }
    idnew = id > 0 ? idnew : -idnew;
    return idnew;
}


// find two quarks from 74=combination of two junction quarks (+ nearby gluons)
//to a diquark recombinations)
int JPythia::findValenceQ74(int ip)
{
  int nq=0;  // number of const. quarks.
  for(int i=0;i<2;i++) {  // loop over  first and second quark.

  int ipos=ip;
  int ntry=0;
  do {
    int mo1=pythia->event[ipos].mother1();
    int mo2=pythia->event[ipos].mother2();
    int id0=pythia->event[ipos].id();
    if(mo1 < 3 && mo2 == 0) {
      nq += constQuark[ipos].first + constQuark[ipos].second;
      break;
    }

    //int id1=event[mo1]->id();
    int id2=pythia->event[mo2].id();
    if(mo1==mo2)
      ipos=mo1;
    else if(i==0 && ntry==0)
      ipos=mo1;
    else if(i==1 && ntry==0)
      ipos=mo2;
    else if(abs(id2) >10) {ipos=mo1;}
    else if(id0*id2 < 0) {ipos=mo1;}
    else if(mo1 !=0 && mo2 !=0) {

	if(i==0) {
	    ipos=mo1;
	    if(pythia->event[ipos].id()==21) ipos=mo2;
	} else if(i==1) {
	    ipos=mo2;
	    if(pythia->event[ipos].id()==21) ipos=mo1;
	}

    }

    } while(ntry++<100);

  }

  return nq;

}

// Set constituent quarks in the quark whose status is -71.
int JPythia::findValenceQ71(int ip, int cq1, int cq2, int side)
{
  Event &event = pythia->event;

  int id0 = pythia->event[ip].id();
  bool isq = pythia->event[ip].isQuark();
  int ipos=ip;
  int ntry=0;
  do {
    int da1 = pythia->event[ipos].daughter1();
    int da2 = pythia->event[ipos].daughter2();
    //int id0=pythia->event[ipos].id();
    if(da1==0 && da2==0) return 0;

    // this quark is not spectator: hard scattered.
    if(da1< ipos && da2 < ipos) return 0;


    if(pythia->event[da1].status()>0 && pythia->event[da2].status()>0) {
      beamSide[ipos]=side;
      constQuark[ipos]=make_pair(cq1,cq2);
      isValenceQ[ipos]=cq1+cq2;
      return ipos;
    } else if(pythia->event[da1].status()>0) {
      beamSide[da1]=side;
      constQuark[da1]=make_pair(cq1,cq2);
      isValenceQ[da1]=cq1+cq2;
      return da1;
    }

    constQuark[da1]=make_pair(cq1,cq2);
    isValenceQ[da1]=cq1+cq2;
    //cout << " da1= "<< da1 << " const q = "<< constQuark[da1].first << " q2= "<< constQuark[da1].second<<endl;

    int ipos0=ipos;
    ipos=da1;
    //if(pythia->event[ipos].id()==21) ipos=da2;
    if(side==2) ipos=da2;
    int id1=pythia->event[da1].id();
    int id2=pythia->event[da2].id();
    if(id1==21 && id2 !=21) ipos=da2;
    else if(id1!=21 && id2 ==21) ipos=da1;

    if(id1==id0 && id2 != id0) ipos=da1;
    else if(id1!=id0 && id2 == id0) ipos=da2;

    if(da1==0 && da2 >0) ipos=da2;
    else if(da1>0 && da2==0) ipos=da1;

    //cout << " ipos= "<< ipos << " isq= "<< isq << " isqq= "<< event[ipos].isDiquark() <<endl;

    // two quarks was combined into diquark.
    if(isq && event[ipos].isDiquark()) {
      int mo1=event[ipos].mother1();
      int mo2=event[ipos].mother2();
      if(mo1!=mo2) {
      int mo=mo1;
      if(mo1 != ipos0 && mo2 == ipos0) {
	mo=mo1;
      } else if(mo1 == ipos0 && mo2 != ipos0) {
	mo=mo2;
      } else {
	cout << "findValenceQ71 mo1= "<< mo1 << " mo2= "<< mo2 << " ipos0= "<<ipos0
	  << " ipos= "<< ipos <<endl;
	exit(1);
      }
      int itry=0;
      int idm = event[mo].id();
      do {
        int mo3 = event[mo].mother1();
        int mo4 = event[mo].mother2();
        if((mo3==1 || mo3==2) && mo4==0) break;
        mo=mo3;
	if(event[mo].id() != idm) mo=mo4;
      } while(++itry<100);
      int nq=constQuark[mo].first + constQuark[ipos].first;
      cq2 = nq; cq1=0;
    }
    }

    if(pythia->event[ipos].id()==21) {
      pythia->event.list();
      pythia->beamA.list();
      pythia->beamB.list();
      cout << "gluon has const.q? ipos= "<< ipos
	<< " ip= "<< ip
	<< " da1= "<< da1 << " da2= "<< da2
	<< " side= "<< side
	<<endl;
      exit(1);
    }

    /*
    // what is this???
    if(da1 != da2 && da2 !=0) {
      if(side==2) {
        ipos=da2;
        if(pythia->event[da2].id()==21) ipos=da1;
      }
    }
    */

  } while(ntry++<100);

  cout << "JPythia::findValcenceQ71 infinite loop? ip= " << ip 
      << " id= " << pythia->event[ip].id() <<endl;
  pythia->event.list();
  pythia->beamA.list();
  pythia->beamB.list();
  exit(1);
}

// determine valence quarks in preparation of hadronization process
// i.e. status -71-79
// cq1[2] constituent quark of beamA cq1[0]: q cq1[1]:qq
void JPythia::findValenceQ(BeamParticle& beamP, const int cq1[2], int side)
{
  // side = 1: beamA,  side=2:beamB
  int nq=cq1[0]+cq1[1]; // total number of constituent quarks.
  bool antiB=false;
  bool isBaryon=beamP.isBaryon();
  if(isBaryon && beamP.id()<0) antiB=true;

  vector<int> ipos;
  Event &event = pythia->event;

  // find the position of valence quark at the side 1.
  // outgoing diffractively scattered
  //if(event[3]->status() == -15 && beamP.size()<3) {
  int da1=event[side+2].daughter1();
  int da2=event[side+2].daughter2();
  bool noPomeron = event[da1].id()!=990 && event[da2].id()!=990;
  if(event[side+2].status() == -15 && noPomeron) {
    ipos.push_back(event[side+2].daughter1());
    ipos.push_back(event[side+2].daughter2());
    // gluon pick up.
    if(event[ipos[0]].id()==21) ipos[0]++;
    constQuark[ipos[0]]=make_pair(cq1[0],0);
    constQuark[ipos[1]]=make_pair(0,cq1[1]);
    isValenceQ[ipos[0]]=cq1[0];
    isValenceQ[ipos[1]]=cq1[1];

  // outgoing elastically scattered
  } else if (event[side+2].status()==14) {
    beamSide[side+2]=1;
    constQuark[side+2]=make_pair(cq1[0],cq1[1]);
    isValenceQ[side+2]=cq1[0]+cq1[1];

  } else if(beamP.size()>0) {

    for(int i=0;i<beamP.size();i++) {
      //bool isConstQ = beamP[i].isValence();
      //if(optConstQscatt==2 && beamP[i].id() !=21) isConstQ=true;

      // Valence quark
      if(beamP[i].isValence()) {
	if(beamP[i].id()==21) {
	  cout << " gluon is valence? "<< beamP[i].id() <<endl;
	  beamP.list();
	  exit(1);
	}
        ipos.push_back(beamP[i].iPos());
	isValenceQ[ipos.back()]=cq1[0]+cq1[1];
	if(!antiB) {
	  if(beamP[i].col() >0) constQuark[ipos.back()]=make_pair(cq1[0],0);
	  if(beamP[i].acol()>0) constQuark[ipos.back()]=make_pair(0,cq1[1]);
	} else {
	  if(beamP[i].acol()>0) constQuark[ipos.back()]=make_pair(cq1[0],0);
	  if(beamP[i].col()>0)  constQuark[ipos.back()]=make_pair(0,cq1[1]);
	}

      // sea quarks
      //} else if((eCM > 500.0) & (optConstQscatt>=3) && (beamP[i].id() !=21)) {
      } else if((optConstQscatt>=3) && (beamP[i].id() !=21)) {
        ipos.push_back(beamP[i].iPos());
	if(beamP[i].col() >0) constQuark[ipos.back()]=make_pair(1,0);
	if(beamP[i].acol()>0) constQuark[ipos.back()]=make_pair(0,1);
	isValenceQ[ipos.back()]=0;
      }
    }

    // junction
    if(ipos.size()==3) {
      if(nq==3) {
	constQuark[ipos[0]]=make_pair(1,0);
	constQuark[ipos[1]]=make_pair(1,0);
	constQuark[ipos[2]]=make_pair(1,0);
        isValenceQ[ipos[0]]=1;
        isValenceQ[ipos[1]]=1;
        isValenceQ[ipos[2]]=1;
      } else if (nq==2) {
	constQuark[ipos[0]]=make_pair(1,0);
	constQuark[ipos[1]]=make_pair(1,0);
        isValenceQ[ipos[0]]=1;
        isValenceQ[ipos[1]]=1;
      } else if(nq==1) {
	constQuark[ipos[0]]=make_pair(1,0);
        isValenceQ[ipos[0]]=1;
      }

    } 

  }

  if(isDebug)cout << " ------ side " << side << " --- npos =" << ipos.size() <<endl;

  for(const auto &i: ipos) {

     if(isDebug)  cout << i << " ipos1 " << i
	 << " q1= " << constQuark[i].first
	 << " q2= " << constQuark[i].second
	    <<endl;

	int pos = findValenceQ71(i,constQuark[i].first,
		constQuark[i].second,side);

         if(isDebug) {
	   cout << "ipos1 " << i << " ipos= " << pos
	        << " q1= " << constQuark[pos].first
	        << " q2= " << constQuark[pos].second
	       << endl;
	 }
  }

}

// determine whether quark ip come from side A or B, 
// and it is valence quark or not.
int JPythia::whichSide(int ip, int& cq)
{
  Event &event = pythia->event;
  BeamParticle &beamA = pythia->beamA;
  BeamParticle &beamB = pythia->beamB;
    cq=0;
    int ipos=ip;
    int ntry=0;
    do {
	//int mo=event[ipos].mother1();
	int mo1=event[ipos].mother1();
	int mo2=event[ipos].mother2();
	if((mo1==3 || mo1==4) && mo2==0) {
	      if(event[ipos].col()>0) cq=1;
	      else if(event[ipos].acol()>0) {
		  if(abs(event[ipos].id())>100) cq=2;
		  else cq=1;
	      }

	      if(mo1==3) return 1;
	      if(mo1==4) return 2;

	} else if(mo1==1 && mo2==0) {
	    if(beamA.size()==0) {
		cq=100;
	    } else if(beamA.size()>0) {
		int iq=-1;
		for(int i=0;i<pythia->beamA.size();i++) {
		    if(beamA[i].iPos()==ipos) {
			iq=i; break;
		    }
		}
		if(iq==-1) {
		    cout << " ipos? " << ipos << " iq= " << iq << endl;
		    for(int i=0;i<beamA.size();i++) {
		      cout << "ipos= "<< beamA[i].iPos()<< endl;
		    }

		    beamA.list(); exit(1); 
		}

	    if(beamA[iq].isValence()) {
	      if(beamA[iq].col()>0) cq=1;
	      else if(beamA[iq].acol()>0) {
		  if(abs(beamA[iq].id())>100) cq=2;
		  else cq=1;
	      }
	    }
	    } else {
	      cout << "JPythia::whichSide side 1 ipos= " << ipos << " cq= " << cq <<endl;
	      exit(1);
	    }


	    return 1;

	} else if(mo1==2 && mo2==0) {
	    if(beamB.size()==0) {
		cq=100;
		//if(abs(event[ipos].id())<10) cq=1;
		//if(abs(event[ipos].id())>10) cq=2;
	    } else if(beamB.size()>0) {
		int i1=-1;
		for(int i=0;i<beamB.size();i++) {
		    if(beamB[i].iPos()==ipos) {i1=i;break;}
		}
	    if(beamB[i1].isValence()) {
	      if(beamB[i1].col()>0) cq=1;
	      else if(beamB[i1].acol()>0) {
		  if(abs(beamB[i1].id())>100) cq=2;
		  else cq=1;
	      }
	    }
	    } else {
	      cout << "JPythia::wichSide side 2 ipos= " << ipos << " cq= " << cq <<endl;
	      exit(1);
	    }
	    return 2;
	}

	if(mo2==0) {
	    ipos=mo1;
	} else {
	int id0=event[ipos].id();
	int id1=event[mo1].id();
	int id2=event[mo2].id();
	if(id0==id1) ipos=mo1;
	else if(id0==id2) ipos=mo2;
	else if(id1==21 && id2==21) ipos=mo1;
	else {
	    // quark fusion to diquark  quarks came from side 1 and 2....
	    if(abs(id0)>100) {
		ipos=mo1; // temporal 
	    } else {
	    cout << "JPythia::whichSide strange? ip= " << ip << endl;
	    cout << "ipos= " << ipos << " id0= " << id0 <<endl;
	    cout << "mo1= " << mo1 << " id1= " << id1 <<endl;
	    cout << "mo2= " << mo2 << " id2= " << id2 <<endl;
	    //event.list();

	    // temporal check later.
	            int ixxx=1;
		    if(pythia->rndm.flat() < 0.5) ixxx=2;
		    return ixxx;

	    //exit(1);
	    }
	}
	}

    } while(ntry++<100);

    cout << "JPythia::wihchSide infinite loop? ip= " << ip <<endl;
    event.list();
    exit(1);
}
    

bool JPythia::generate(int ncoll, int operation, const int cq1[2], const int cq2[2],
	Vec4& xA, Vec4& xB, vector<EventParticle*>& outgoing)
{
  Event &event = pythia->event;
  int procid=0;

  if(userHook && optConstQuarkDiffra==5) {
    if(preHadronA && preHadronB) procid=105;
    else if(preHadronA  &&  !preHadronB)  procid=104;
    else if(!preHadronA  &&  preHadronB)  procid=103;
  }

    /*
    HoldProcess hold(selectSASD, procid, -1);
    if(procid>0) {
       cout << " procid= "<< procid
	   << " slectSASD.proc= "<< selectSASD.proc
	   <<endl;
    }
    */

  selectSASD.proc = procid;

  // Generate pythia event.
  if(!pythia->next()) {
    cout << " jpythia failed code="<< pythia->info.code()
	 << " beamA.size= " << pythia->beamA.size()
	 << " beamB.size= " << pythia->beamB.size()
	 <<endl;
    event.list(0);
    return false;
  }
  procCode=pythia->info.code();

  if(procid>0 && procid != procCode) {

  cout << " selectSASD= "<< &selectSASD << " .proc= "<< selectSASD.proc <<endl;
    cout << " phthia generatene procid= "<< procid
          << " code= "<< pythia->info.code()
	<< " preHA= " << preHadronA
	<< " preHB= " << preHadronB
	<< " idA= " << pythia->beamA.id()
	<< " idB= " << pythia->beamB.id()
	<< " Ecm= "<< eCM
	<<endl;
    exit(1);
  }


  if(isDebug) {
    event.list(0);
    pythia->beamA.list();
    pythia->beamB.list();

  // 101 nondiffractive
  // 102 elastic
  // 103,104 single diffractive
  // 105 double diffractive
  // 106 central diffractive
  cout << " code = " << procCode
       << " beamA.size= " << pythia->beamA.size()
       << " beamB.size= " << pythia->beamB.size()
       <<endl;
  }
    

  constQuark.resize(event.size(),make_pair(0,0));
  beamSide.resize(event.size(),0);
  isValenceQ.resize(event.size(),0);

  // find the constituent quarks and put them into the array constQuark.
  //findValenceQ(cq1,cq2);
  findValenceQ(pythia->beamA,cq1,1);
  findValenceQ(pythia->beamB,cq2,2);

  int ipos=0;
  do {
    if(!event[ipos].isFinal()) continue;
    //int status=event[ipos].status();
    //if(status==81 || status==82 || hadrons from ministring into one hadrons.
    //   status==83 || status==84 || hadrons from normal string.
    //   status==85 ||  // hadron from a junction
    //   status==86 ||  // hadron from a junction
    //   status==87 ||  // baryon from a junction
    //   status==89 ||  // baryon from a junction from ministring
    //   status==14 ||
    //   status==62)  

    int mo1 = event[ipos].mother1();  // parton1
    int mo2 = event[ipos].mother2();  // parton2
    int da1 = event[mo1].daughter1(); // parton1's daughter
    int da2 = event[mo1].daughter2(); // parton2's daughter
    int side1=0, side2=0;
    int cnstq1=0, cnstq2=0;
    // directory produced from diffractive scattering.
    // status must be 14
    if(mo2==0 || mo1==mo2) {
      da1=ipos;
      da2=ipos;
      //da1 = event[ipos].daughter1();
      //da2 = event[ipos].daughter2();
      side1 = mo1;
      side2 = mo1;
      beamSide[ipos]=mo1;
      if(mo1==1) {
	// exclude leptons.
        if(event[ipos].isHadron()) {
	  constQuark[ipos]=make_pair(cq1[0],cq1[1]);
	  isValenceQ[ipos]=cq1[0]+cq1[1];
	}
      } else if(mo1==2) {
        if(event[ipos].isHadron()) {
	  constQuark[ipos]=make_pair(cq2[0],cq2[1]);
	  isValenceQ[ipos]=cq1[0]+cq1[1];
	}
      }
    } else {
      side1 = whichSide(mo1,cnstq1);
      side2 = whichSide(mo2,cnstq2);
    }

    int status1=event[da1].status();
    int status2=event[da2].status();

    if(isDebug) {
      cout << " status1= " << status1 << " status2= " << status2<<endl;
      /*
      for(int i=mo1;i<=mo2;i++) {
	cout << i << " side= " << beamSide[i]
	     << setw(4) << event[i].status()
	    << setw(6) << " id= " << event[i].id()
	    << setw(5) << "  " << event[i].name()
	    << setw(12) << "  q1= " << constQuark[i].first
	    << setw(7) << " q2= " << constQuark[i].second
	    <<endl;
      }
      */
    }

    // ministring into one hadron.
    if(status1==81 || status2==81) {

      if(da1 != da2) {
        cout << "ministring  da1 = " << da1 << " da2= " << da2 <<endl;
      }

      if(abs(event[mo1].id())>10) 
        beamSide[da1]=side1;
      else
        beamSide[da1]=side2;

 // form ministring into one hadron.
    int nq = constQuark[mo1].first+constQuark[mo1].second
           + constQuark[mo2].first+constQuark[mo2].second;
    constQuark[da1]=make_pair(0,0);
    if(event[da1].particleDataEntry().isMeson()) {
      if(nq==2) {
        constQuark[da1]=make_pair(1,1);
      } else if(nq==1){
        constQuark[da1]=make_pair(1,0);
      }
    } else {
      if(nq==3) {
        constQuark[da1]=make_pair(1,2);
      } else if(nq==2){
        constQuark[da1]=make_pair(1,1);
      } else if(nq==1){
	constQuark[da1]=make_pair(1,0);
      }
    }

    if(isDebug) {
	cout << "* 81 " << da1<<endl;
	cout << " constq1 da1= " << constQuark[da1].first
	     << " constq2 da1= " << constQuark[da1].second
	     <<endl;
	cout << " constq1 mo1= " << constQuark[mo1].first
	     << " constq2 mo1= " << constQuark[mo1].second
	     <<endl;
	cout << " constq1 mo2= " << constQuark[mo2].first
	     << " constq2 mo2= " << constQuark[mo2].second
	     <<endl;
    }

    // 82 ministring into two hadrons. 89: from junction.
    } else if((status1==82 || status1==89) && (status2==82 || status2==89)) { 
      int idq1=event[mo1].id();
      int idq2=event[mo2].id();
      if(particleData->isBaryon(event[da1].id())) {
	    if(abs(idq1)>10) {
		beamSide[da1]=side1;
		beamSide[da2]=side2;
	    } else if(abs(idq2)>10) {
		beamSide[da1]=side2;
		beamSide[da2]=side1;
	    }
      } else if(particleData->isBaryon(event[da2].id())) {
	    if(abs(idq1)>10) {
		beamSide[da2]=side1;
		beamSide[da1]=side2;
	    } else if(abs(idq2)>10) {
		beamSide[da2]=side2;
		beamSide[da1]=side1;
	    }
      }

    } else if (status1 >=83 && status1 <=88) {
	for(int i=da1;i<=da2;i++) {
	    if(event[i].status()==83 ||
	       event[i].status()==85 ||
	       event[i].status()==87) beamSide[i]=1;
	    if(event[i].status()==84 ||
	       event[i].status()==86 ||
	       event[i].status()==88) beamSide[i]=2;
	}

    } else if (status1 ==14) { // outgoing elastically scattered

    } else if (status1 ==91) { // normal decay product

    } else if (status1 ==61) {
	beamSide[da1]=1;
    } else if (status1 ==62) {
	beamSide[da1]=2;
    } else if (status1 ==63) { // outgoing beam remnant
    } else {

	cout << "JPythia::generate ? " <<endl;
        cout << " status1= " << status1
	    << " status2= " << status2<<endl;
	event.list();
	exit(1);

    }

    if(da2-da1>0) {

      /*  This is not needed anymore because it is taken account in findValenceQuark71()
      //check 74=combination of two junction quarks (+ nearby gluons)
      for(int ip=mo1;ip<=mo2;ip++) {
        if(event[ip].status()==-74) {
          if(event[ip].id() > 0) {
            constQuark[ip]=make_pair(0,findValenceQ74(ip));
          } else {
            constQuark[ip]=make_pair(findValenceQ74(ip),0);
	  }
        }
      }
      */


      int posq[3]={0,0,0};
      int nquark=0;
      int istart=0;
      int nq1=0, nq2=0;
      for(int i=mo1;i<=mo2;i++) {
	if(event[i].id() == 21) continue;
	if(istart==0) istart=i;

	nq1 += constQuark[i].first;
	nq2 += constQuark[i].second;
	if(event[i].isQuark()) posq[nquark++]=i;
      }

      if(isDebug) {
      cout << "mo1= " << mo1 << " mo2= " << mo2<<endl;
      cout << "da1= " << da1 << " da2= " << da2<<endl;
      cout << " istart= " << istart<<endl;
      cout << "da1= " << da1 << " da2= " << da2<< " da2-da1= " <<da2-da1<<endl;
      cout << "da1 isBaryon= " <<event[da1].particleDataEntry().isBaryon()
           << " da2 isBaryon= " <<event[da2].particleDataEntry().isBaryon()
	   <<endl;
      cout << " nquark= " << nquark << " posq= " << posq[0]
	  << " 1= " << posq[1]
	  << " 2= " << posq[2]
	  <<endl;
      }

  // junction.  
  if(nquark==3) {
    constQuark[da1]   = constQuark[posq[0]];
    constQuark[da2-1] = constQuark[posq[1]];
    constQuark[da2]   = constQuark[posq[2]];

  } else if(da2-da1==1) {
    if(event[da1].particleDataEntry().isBaryon()) {
      if(event[mo1].isDiquark()) { 
        constQuark[da1]=constQuark[mo1];
        constQuark[da2]=constQuark[mo2];
      } else {
        constQuark[da1]=constQuark[mo2];
        constQuark[da2]=constQuark[mo1];
      }

      /*
      // di-quark breaking.
      if(constQuark1[da2]==0 && constQuark2[da1]==2) {
	  constQuark1[da2]=1;
	  constQuark2[da1]=1;
      }
      */

    } else {
      if(event[mo1].isDiquark()) { 
        constQuark[da1]=constQuark[mo2];
        constQuark[da2]=constQuark[mo1];
      } else {
        constQuark[da1]=constQuark[mo1];
        constQuark[da2]=constQuark[mo2];
      }

      /*
      // di-quark breaking.
      if(constQuark1[da1]==0 && constQuark2[da2]==2) {
	  constQuark1[da1]=1;
	  constQuark2[da2]=1;
      }
      */

    }


  } else {

    if(constQuark[istart].first==1) {
	constQuark[da1].first=1; 
    }
    if(mo2-istart > 2 && constQuark[mo1+1].first==1) {
	constQuark[da1+1].first=1; 
    }

    if(constQuark[mo2].second==2) {
	if(particleData->isBaryon(event[da2].id())) {
	  if(optDiq==0) {
	    constQuark[da2].second=2;
	  } else {
	    constQuark[da2].second=1;
	    constQuark[da2-1].second=1;
	  }
	} else {
	    constQuark[da2].first=1;
	    constQuark[da2-1].first=1;
	}
    } else if(constQuark[mo2].second==1) {
	constQuark[da2].second=1; 
    }

    // opposite side.
    if(constQuark[mo2].first==1) {
	  constQuark[da2].first=1; 
    }
    if(mo2-istart ==2  && constQuark[mo2-1].first==1) {
	constQuark[da2-1].first=1; 
    }

    if(constQuark[istart].second==2) {
	if(particleData->isBaryon(event[da1].id())) {
	  if(optDiq==0) {
	    constQuark[da1].second=2;
	  }else {
	    constQuark[da1].second=1;
	    constQuark[da1+1].second=1;
	  }
	} else {
	    constQuark[da1].first=1;
	    constQuark[da1+1].first=1;
	}
    } else if(constQuark[istart].second==1) {
	constQuark[da1].second=1; 
    }

  }

  /*
  // baryon junction.
  for(int i=da1;i<=da2;i++)
    if(event[i].status()==87 || event[i].status()==88) {
	if(constQuark1[i]==0) constQuarkk1[i]=1;
	else if(constQuark1[i]==1) 
	    constQuark2[i]=1;
  }
  */

    /*
    // baryon junction.
    for(int i=da1;i<=da2;i++) {
	if(event[i].status()==87 && constQuark1[i]==0) constQuark1[i] +=1;
	if(event[i].status()==88 && constQuark2[i]==0) constQuark2[i] +=1;
    }
    */


    }

    if(isDebug) {
    cout << "side1= " << side1 << " cq1= " << cnstq1 <<endl;
    cout << "side2= " << side2 << " cq2= " << cnstq2  <<endl;
    cout << "==" << " mo1= "<< mo1 << " mo2 "<< mo2
         << " ============================================"<<endl;
    for(int i=mo1;i<=mo2;i++) {
	cout << i << " side= " << beamSide[i]
	     << setw(4) << event[i].status()
	    << setw(6) << " id= " << event[i].id()
	    << setw(5) << "  " << event[i].name()
	    << setw(12) << "  q1= " << constQuark[i].first
	    << setw(7) << " q2= " << constQuark[i].second
	    <<endl;
    }
    cout << "==" << " da1= "<< da1 << " da2 "<< da2 << 
	" =========================================="<<endl;
    for(int i=da1;i<=da2;i++) {
	cout << i << " side= " << beamSide[i]
	     << setw(4) << event[i].status()
	    << setw(6) << " id= " << event[i].id()
	    << setw(5) << "  " << event[i].name()
	    << setw(12) << "  q1= " << constQuark[i].first
	    << setw(7) << " q2= " << constQuark[i].second
	    <<endl;
    }
    cout << "============================================"<<endl;
    }


    ipos=da2;
  } while(++ipos < event.size()); // end parton loop.

  int nq=0; 
  for(int i=0;i<(int)constQuark.size();i++) {
    if(event[i].isFinal()) nq += constQuark[i].first+constQuark[i].second;
  }
  int nq0 = cq1[0]+cq1[1]+cq2[0]+cq2[1];

    /*
    if(nq0 != nq) {
     cout << "*** nq0= " << nq0 << " nq= " << nq <<endl;
	    event.list();
       for(int i=0;i<(int)event.size();i++) {
	 cout << i << setw(6) << event[i].id()
	   << setw(15) << event[i].name()
	   << setw(5) << event[i].status()
	   << setw(6) << " q1= "<< constQuark[i].first
	   << setw(6) << " q2= "<< constQuark[i].second
	   <<endl;
       }
       exit(1);
    }
    */

     if(isDebug) {
     cout << "*** nq= " << nq << " cq1+cq2= " << nq0 <<endl;
    //if(nq0 != nq)  cin.get();
       //if(code == 103 || code == 104 || code == 105) cin.get();
       for(int i=0;i<(int)event.size();i++) {
	 if(event[i].isFinal()) {
	 cout << i << setw(6) << event[i].id()
	   << setw(15) << event[i].name()
	   << setw(5) << event[i].status()
	   << setw(6) << " isValence= "<< isValenceQ[i]
	   << setw(6) << " q1= "<< constQuark[i].first
	   << setw(6) << " q2= "<< constQuark[i].second
	   << setw(6) << " pz= "<< event[i].pz()
	   << setw(6) << " pt= "<< event[i].pT()
	   <<endl;
	 }
       }
     }

    // Boost event from CM frame to lab frame.
    event.rotbst(MfromCM);

    transportHadron(xA, xB, ncoll,operation,outgoing);

    constQuark.clear();
    beamSide.clear();
    isValenceQ.clear();

    return true;

}

// transfer pythia hadron into jam array.
void JPythia::transportHadron(Vec4& xA, Vec4& xB, int ncoll, int operation,
	vector<EventParticle*>& outgoing)
{
  Event &event = pythia->event;
  int kaon=0;
  for(int i=1; i< event.size();i++)
  if(event[i].isFinal()) {
    int id=event[i].id();
    if(id==130 || id==310) {
      if(kaon==0) {
        if(pythia->rndm.flat() < 0.5) {
	  id=311; kaon=1;
        } else {
          id=-311; kaon=-1;
        }
      } else if(kaon==1) {
        id=-311; kaon=0;
      } else if(kaon== -1) {
        id=311; kaon=0;
      }
    }
    ParticleDataEntry* pd =particleData->findParticle(id);
    //ParticleDataEntry* pd =&(event[i].particleDataEntry());

    int pid=jamParticleData->pid(id);
    EventParticle* pa = new EventParticle(id,pd);
    pa->setPID(pid);
    pa->setMass(event[i].m());
    Vec4 p=event[i].p();
    pa->setMomentum(p);
    pa->setOnShell();
    pa->setNumberOfColl(ncoll);
      if(allowRescatterSameString) pa->lastColl(-1);
      else pa->lastColl(operation);
    pa->setBoostMatrix(MfromCM);

    // Find leading hadron.
    bool isq = constQuark[i].first>0 || constQuark[i].second>0;
    int nq=constQuark[i].first + constQuark[i].second;
    if(isq && nq==0) {
      cout << "isq= "<<isq << " nq= "<< nq<<endl;
      exit(1);
    }

    Vec4 v=xA;
    if(beamSide[i]==2) v=xB;
    if(isq==0) v = 0.5*(xA+xB);
    Vec4 xprod = event[i].vProd()*MM2FM;
    double tform=event[i].tProd()*MM2FM;
    double tf = v[0]+tform;

    // no constituent quark scattering.
    if(optConstQscatt==0) v +=  xprod;

    else if(optConstQscatt==1 || optConstQscatt==3) {
      if(isq==0 ) v += xprod;

    // no formation time
    } else if(optConstQscatt==2) {
      constQuark[i].first=1;
      constQuark[i].second= pa->baryon() == 0 ? 1 : 2;
      tform=0.0;
      tf=v[0];

    } else if(optConstQscatt==4) { // nondiffractive scattering
      if(isq==0 ) v += xprod;
      else if(procCode==101) {
	if(nq<2 && isValenceQ[i]==0) v += xprod;  // hadrons with one quark do not interact.
      }
    } else if(optConstQscatt==5) { // nondiffractive scattering
      if(isq==0 ) v += xprod;
      else if(procCode==101) {
        if(pa->baryon()==0 && isValenceQ[i]==0) v += xprod;
      } else {
	v += xprod;
      }
    } else if(optConstQscatt==6) {
      if(isq==0 ) v += xprod;
      else if(isValenceQ[i]==0) {
	v += xprod;
      }
    }

    pa->setConstQuark(0,constQuark[i].first);
    pa->setConstQuark(1,constQuark[i].second);
    pa->setCoordinate(v);
    pa->setVertex(v);
    pa->setFormationTime(tf);

    if(beamSide[i]==1) {
      pa->setParent(100);
    } else {
      pa->setParent(200);
    }

    // status=14 is a outgoing elastically scattered, thus
    // no need to set formation time.
    //if(event[i].status() ==14) pa->setFormationTime(v[0]);

    // compute decay time.
    //double m=p.mCalc();
    double m=event[i].m();
    double e=event[i].e();
    double dect = jamParticleData->lifeTime(pd,m,e);
    //pa->setLifeTime(v[0] + dect);
    pa->setLifeTime(tf + dect);

    //double tform=0.8*e/m;
    //pa->setFormationTime(t+dect+tform);

    outgoing.push_back(pa);

    // accumulate the statistics of the formation time.
    double gam=e/m;
    if(isq==1) {
      aveConstFormationTime += tform/gam;
      nCFormationTime++;
    } else {
      aveFormationTime += tform/gam;
      nFormationTime++;
    }

  } // end particle loop.

}

} // end namespace jam2
