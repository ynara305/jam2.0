#include <jam2/interaction/HadronDecay.h>
#include <jam2/xsection/XsecTable.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/hadrons/JamStdlib.h>
#include <jam2/hadrons/GaussPoints.h>
#include <Pythia8/PythiaStdlib.h>

using namespace std;

namespace jam2 {

using Pythia8::pow3;
using Pythia8::pow4;
using Pythia8::Vec4;
using Pythia8::ParticleDataEntry;

HadronDecay::HadronDecay()
{
  optDeltaWidth=4; // Randrap mstc(63)
  //optDeltaWidth=1; // Frankfurt mstc(63)
  //optDeltaWidth=3; // GiBUU mstc(63)

  optWidth=3; //  BlattWeisskopf mstc(65) 
  //optWidth=1; //mstc(65) Frankfurt

  minKinE=0.003; // parj(64)

  useTable=true;
}

void HadronDecay::init(Pythia8::Info *infor,Pythia8::Settings& s,
	JamParticleData* pd, Pythia8::Pythia* d, Pythia8::Rndm* r) 
{
  settings = &s;
  info=infor;
  jamParticleData=pd;
  particleData=jamParticleData->getParticleData();
  pythiaDecay=d;
  rndm = r;
  pythiaEvent.init("(hadron decay)",particleData);
  jdecay.init(info,s,particleData,rndm);

  potentialHandling=settings->mode("ParticleDecays:potentialHandling");
  optDeltaDecayAngle=settings->mode("ParticleDecays:optDeltaDecayAngle");
  optDecayAngle=settings->mode("ParticleDecays:optDecayAngle");
  optDecayAngleS=settings->mode("ParticleDecays:optDecayAngleSchannel");
  optPotential=settings->mode("MeanField:optPotential");
  if(settings->mode("MeanField:mode")==0) optPotential=0;
  
  optPosition=settings->mode("Cascade:optDecayPosition");
  paramSmear=settings->parm("Cascade:decayPositionSmearParam");

  // Propagation by canonical momentum in RQMD mode
  optVectorPotential=settings->mode("MeanField:optVectorPotential");
  optVdot=settings->mode("MeanField:optVdot");
  optPropagate=0;
  if(optVectorPotential==1 && optVdot==0) optPropagate=1;

}

void HadronDecay::decay(EventParticle* pa0, vector<EventParticle*>& outgoing,
	bool finalDecay)
{
  double dectime = pa0->lifetime();
  int id0=pa0->getID();
  int baryon0 = pa0->baryon();
  Vec4 r0 = pa0->getR();
  Vec4 p0 = pa0->getP();
  if(!finalDecay) r0 = pa0->propagate(dectime,optPropagate);
  int* cq=pa0->constQuark();
  double ftime=pa0->getTf();
  Vec4 potv0 = pa0->potv();
  double pots0 = pa0->pots();

  //double meff=pa0->getEffectiveMass();
  double meff0=p0.mCalc();
  double m0 = pa0->getMass();
  //double SM = meff0 / m0;
  double SM=1.0;

  //if(pa0->baryon()==0 && abs(pa0->getMass() - p0.mCalc())>1e-8) {
  if(abs(pa0->getEffectiveMass() - p0.mCalc())>1e-5) {
	cout << " decay mass ? m= " << pa0->getMass()
	    << " meff = "<< meff0
	    << " pots = "<< pa0->pots()
	    << " m0 = "<< p0.mCalc()
	    << " id= "<< id0
	    << " ncol= "<<pa0->getNColl()
	    << " last= "<<pa0->lastColl()
	    << " sta= "<<pa0->getStatus()
	    << endl;
  }

  if(optPropagate==1) {
    p0 -= potv0;
    p0[0] = sqrt(meff0*meff0+p0.pAbs2());
  }

  if(dectime >= ftime) {
    cq[0]=1;
    cq[1]=abs(baryon0) == 3 ? 2 : 1;
  }

  // go to the previous collision C.M.frame.
  Vec4 pcm = p0;
  if(SM==1.0) pcm[0]=sqrt(m0*m0+pcm.pAbs2());

  RotBstMatrix MfromCM = pa0->getMfromCM();
  RotBstMatrix MtoCM = MfromCM;
  MtoCM.invert();
  pcm.rotbst(MtoCM);

  // This resonance was created by s-channel reaction?
  bool schannel = (pcm.pAbs2() < 1e-9) ? true : false;

  pythiaEvent.clear();
  pythiaEvent.append(id0, 1,0,0,pcm/SM, m0);
  //pythiaEvent.append(id0, 1,0,0,pcm/SM, meff0);

  //pythiaEvent[0].setPDEPtr(pa0->getParticleDataEntry());

  int pid0=pa0->getPID();
  if( pid0 == id_delt) {
    jdecay.setOptAngle(optDeltaDecayAngle);
  //} else if(schannel) {
  //  jdecay.setOptAngle(optDecayAngleS);
  //2020/1/3
  } else if( pa0->baryon() !=0) {
  //} else if( pid0 == id_delts || pid0 == id_nucls) {
    jdecay.setOptAngle(optDecayAngle);
  } else {
    jdecay.setOptAngle(0);
  }

  bool idec=false;
  switch (pid0) {
    case id_pi:
    case id_light0:
    case id_light1:
    case id_str:
    case id_nucl:
    case id_nucls:
    case id_delt:
    case id_delts:
    case id_lambda:
    case id_lambdas:
    case id_sigma:
    case id_sigmas:
    case id_xi:
    case id_xis:
    case id_omega : idec = true; break;
    default: idec=false;
  }

  bool icon;
  // Perform decay of resonance.
  if(idec) {

    int optDec=1;
    if(optDec==1) {
	icon=jdecay.decay(0,schannel,pythiaEvent);
    } else {
      ParticleDataEntry* pd=pa0->getParticleDataEntry();
      double totwid=jamParticleData->totalWidth(pd,m0);
      if(totwid > 0.0) {
        icon=jdecay.decay2(0,schannel,pythiaEvent,totwid,pWidth,nChannel);
      } else {
        if(finalDecay) {
          icon=jdecay.decay(0,schannel,pythiaEvent);
        } else {
          cout <<" HadronDecay::decay totwid=0.0 "<< totwid
	     << " id= "<< id0 << " m= "<< m0 
	     << " wid= "<< pd->mWidth()
	     << " dectime= "<< dectime
	     << " finaldecay= "<< finalDecay
	     <<endl;
          exit(1);
        }
      }
    }

  } else {
      pythiaDecay->event.reset();
      pythiaDecay->event.append(pythiaEvent[0]);
      icon=pythiaDecay->moreDecays();
      //icon=pythiaDecay->forceHadronLevel();
      pythiaEvent = pythiaDecay->event;
  }

  pythiaEvent.rotbst(MfromCM);

  if(!icon) {
    pythiaEvent.list();
    ParticleDataEntry* pd=pa0->getParticleDataEntry();
    cout << "Decay failed decay= " << id0 << " name= " << pd->name()
	<< " m= "<< pa0->getMass()
	<< " t= "<< dectime
	<< " SM= "<< SM
	<< " idec= "<< idec
	<< " schannel= "<< schannel
	<< endl;
    cout << " width= " << pd->mWidth()
	  << " useBreitWigner= " << pd->useBreitWigner() <<endl;

    if(abs(id0)<22) {
	cout << " lepton phton decay??? " << id0 <<endl;
    }
    int nd=pd->sizeChannels();

    if(finalDecay && (pd->mWidth()< 1e-9 || nd==0)) {
       ParticleDataEntry* pd0=pa0->getParticleDataEntry();

	cout << "This particle is stable! nd= "<< nd<<endl;
	cout << " canDecay= "<< pd0->canDecay()
	     << " mayDecay= "<< pd0->mayDecay()
	     <<endl;

	outgoing.push_back(new EventParticle(*pa0));
	outgoing[0]->setLifeTime(1e+35);
	return;
    }

    double brsum=0.0;
    for(int idl=0;idl<nd;idl++) {
	DecayChannel d=pd->channel(idl);
	double mtot=0.0;
	for(int ibr=0;ibr<d.multiplicity();ibr++) {
	  int id=d.product(ibr);
          double mr=particleData->mWidth(id) > 1e-7 ?
                    particleData->mMin(id) : particleData->m0(id);
          mtot += mr;
	  cout << id << " m= "<< mr << endl;
	}
	cout << endl;
	if(pa0->getMass() >= mtot + 0.0005) brsum+=d.bRatio();
	cout << " r= "<< d.bRatio() << " onmode= "<< d.onMode()<<endl;
    }
    cout << "brsum= "<< brsum<<endl;
    exit(1);
  }

  // Find a leading hadron.
  int iend[3]={0,0,0}, constq[3]={0,0,0};

  // If hadrons decay before their formation time.
  //findConstQuark(pythiaEvent,isBaryon,cq,iend,constq);

  int numCollision=abs(pa0->getNColl());
  int lastcl=pa0->lastColl();
  int kaon=0;
  Vec4 ptot=0.0;
  // Loop over decayed particles.
  for(int i=0;i<(int)pythiaEvent.size(); i++) {
	if(pythiaEvent[i].status() < 0) continue;

	int id=pythiaEvent[i].id();

	if(!finalDecay) {
	    if(id==130 || id==310) {
		if(kaon==0) {
		    if(rndm->flat() < 0.5) {
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
	}

	if(id==83) {
	    cout << " particle decay id? " << id<<endl;
	    cout << " m= "<< pythiaEvent[i].p()
		<<endl;
	    exit(1);
	}
	Vec4 p=pythiaEvent[i].p();

	double m=pythiaEvent[i].m();
	double e=p[0];
	p *= SM;

	Vec4 rp = r0;
	if(optPosition && i>0) {
	  if(optPosition==1) {
          double del=paramSmear*sqrt(rndm->flat());
	  double cos1=1.0 - 2.0*rndm->flat();
	  double sin1=sqrt(1.0-cos1*cos1);
	  double phi1=2*M_PI*rndm->flat();
	  Vec4 dr(del*sin1*cos(phi1),del*sin1*sin(phi1),del*cos1,0.0);
          dr.bst( p0, meff0 );
	  rp += dr;
	  if(rp[0]< r0[0]) {
	    cout << " after decay time is strange r0= "<< r0[0]
	      << " t= "<< rp[0]
	      << " dt= "<< dr[0]<<endl;
	    exit(1);
	  }
	  } else {
            double delt=paramSmear*sqrt(rndm->flat());
	    rp[1] += delt*p[1]/p[0];
	    rp[2] += delt*p[2]/p[0];
	    rp[3] += delt*p[3]/p[0];
	    rp[0] += delt;
	  }

	}

	//ParticleDataEntry* pd=&pythiaEvent[i].particleDataEntry();
	ParticleDataEntry* pd = particleData->findParticle(id);

	EventParticle* ep = new EventParticle(id,m,rp,p,pd,numCollision);
	int pid=jamParticleData->pid(id);
	ep->setPID(pid);
        ep->setParent(id0);

/*
//check
  if(ep->baryon()==0 && abs(m - p.mCalc())>1e-3) {
	cout << " decay mass ? m= " << ep->getMass()
	  << " m= "<< m
	    << " meff = "<< p.mCalc()
	    << " p = "<< p
	    << " pots = "<< ep->pots()
	    << endl;
	pythiaEvent.list();
  }
*/

	double mef2 = p.m2Calc();
	// photon can have very small negative mass....
	if(id > 22 && mef2 < 0.0) {
	    cout << "HadronDecay::decay after decay meff < 0 ? "<< mef2
		<< " id0= " << id0
		<< " id= " << id
		<< " pd= " << p
		<< " m= " << m
		<< " SM= " << SM
		<< " S= " << ep->pots()
		 <<endl;
	    pythiaEvent.list();
	    exit(1);
	} else {
	  if(mef2>0.0) {
	    if(SM !=1.0) ep->setPotS( sqrt(mef2)-m );
	  }
	}

	//ep->setOnShell();
	double dect = jamParticleData->lifeTime(pd,m,e);
	ep->setLifeTime(rp[0]+dect);
	ep->lastColl(lastcl);
	//ep->setNumberOfColl(numCollision);
	ep->setBoostMatrix(MfromCM);

	// Find virtual hadrons.
	if(iend[0]==i) {
	  ep->setConstQuark(0,constq[0]);
          ep->setFormationTime(ftime);
        } else if(iend[1]==i) {
	  ep->setConstQuark(1,constq[1]);
          ep->setFormationTime(ftime);
        }

        if(ep->isMeanField(rp[0],optPotential)) {
          ep->setPotS( pots0 );
          ep->setPotV( potv0 );
        }

	outgoing.push_back(ep);
	ptot += ep->getP();

  } // End Loop over decayed particles

  // Option that only the leading hadron has an effective mass after decay.
  // Other particles are put into free.
  bool optPotH=false;
  if(potentialHandling==1 && abs(SM-1.0) > 1e-10) {
    optPotH=true;
    int leading=-1;
    double ene = 0.0;
    for(int i=0; i<(int)outgoing.size(); i++) {
      if(leading==-1 && outgoing[i]->baryon()==baryon0) {
        leading=i;
      } else {
	outgoing[i]->setPotS(0.0);
        outgoing[i]->setOnShell();
	ene += outgoing[i]->getPe();
      }
    }
    if(leading == -1) {
	cout << " decay leading = "<< leading << " baryon0= " << baryon0 <<endl;
	exit(1);
    }
    double el = p0[0] - ene;
    outgoing[leading]->setE(el);
    double meff= el*el - outgoing[leading]->pAbs2();
    if(meff < 0.0) {
	cout << setprecision(8)<< "after decay meff<0? " << meff
             << " SM= " << SM
             << " id0= " << id0
             << " p0= " <<  p0[0]
             << " id= " << outgoing[leading]->getID()
              << " m= " << outgoing[leading]->getMass()
		 <<endl;
	exit(1);
    }
    if(optPropagate==1)  outgoing[leading]->addP(potv0);
      outgoing[leading]->setPotS(sqrt(meff) - outgoing[leading]->getMass());
      double th=outgoing[leading]->getT();
      if(outgoing[leading]->isMeanField(th,optPotential)) {
	outgoing[leading]->setPotV( pots0 );
	outgoing[leading]->setPotV( potv0 );
      }
  }

  if(!optPotH && optPropagate==1) {
    int n=outgoing.size();
    for(int i=0; i < n; i++) {
      outgoing[i]->addP(potv0/n);
      outgoing[i]->setOnShell();
    }
  }

//...Check energy momentum conservation.
    Vec4 diff= ptot-p0;
    //double pdiff = abs(diff.mCalc());
    double pdiff = diff.pAbs();
    if(pdiff > 1e-5) {
      // more detailed condition using relative pdiff
      double const ptot_euclid_norm2 = ptot.e() * ptot.e() + ptot.pAbs2();
      double const p0_euclid_norm2 = p0.e() * p0.e() + p0.pAbs2();
      double const relative_pdiff = pdiff / std::sqrt(std::max(ptot_euclid_norm2, p0_euclid_norm2));
      if (relative_pdiff > 1e-5) {
	std::ios_base::fmtflags saved_flags(std::cerr.flags());
	cerr << std::scientific << std::setprecision(8);
	cerr << "(HadronDecay:) Energy Momentum not conserved id= "
	     << id0
	     << " m= " << pa0->getMass()
	     << " outgoing = " << outgoing.size()
	     << endl;
	cerr << "p0  = " << p0 << endl;
	cerr << "ptot= " << ptot << endl;
	for(int i=0; i<(int)outgoing.size(); i++)
	    cerr << outgoing[i]->getID() << " pf = " << outgoing[i]->getP()<< endl;
	cerr << endl;
        cerr << "diff= " << diff << " error = " << pdiff  << " relative_error = " << relative_pdiff << endl;
	cerr.flags(saved_flags); // restore settings
	pythiaEvent.list();
	exit(1);
      }
    }

    /*
    // check constituent quark.
    if(nconstq != nq ) {
       pythiaEvent.list();
	cout << "HadronDecay nconstq= " << nconstq
	    << " nq= " << nq <<endl;
    cout << " isbaryon= "<< isBaryon<<endl;
    cout << " q0= "<< cq[0] 
         << " q1= "<< cq[2]
	 <<endl;

    cout << " constq1=" << constq[0]
         << " constq2=" << constq[1]
         << " constq3=" << constq[2]
	 <<endl;
    cout << " iend1  =" << iend[0]
         << " iend2  =" << iend[1]
         << " iend3  =" << iend[2]
	 <<endl;
	exit(1);
    }
    */

}

void HadronDecay::findConstQuark(Event& event,int ibar,
      int* q, int* iend,int* constq)
{
// q[0],q[1] const. quark of mother hadron.
// ibar:  baryon number of mother hadron.

  //const int optDiq1=1;
  //const int optDiq2=1;

  int nconstq = q[0]+q[1];

  if(ibar == 0) {
    if(nconstq==2) return; // no vertical hadron.
    if(q[0]==0 && q[1]==1) {
      constq[0]=0;
      for(int i=1; i<event.size();i++) 
      if(event[i].isFinal()) {iend[0]=i; break;}
    } else if(q[1]==0 && q[0]==1) {
      iend[1]=event.size()-1;
      constq[1]=0;
    } else if(q[1]==0 && q[0]==0) {
      constq[0]=0;
      constq[1]=0;
      for(int i=1; i<event.size();i++) 
      if(event[i].isFinal()) {iend[0]=i; break;}
      iend[1]=event.size()-1;
    }
    return;
  }

  // Baryon mother start here.
  if(nconstq==3) return; // no vertical hadron.

  // Count number of hadrons.
  int nhadron=0;
  int istart=0;
  for(int i=1; i<event.size();i++) {
    if(event[i].isFinal()) {
	if(istart==0)  istart = i;
	nhadron++;
    }
  }
  int ilast= event.size()-1;
  int ibar1 = event[istart].particleDataEntry().isBaryon();
  int ibar2 = event[ilast].particleDataEntry().isBaryon();

  /*
  cout << " q1= "<< q[0] << " q2= " << q[1] <<endl;
  cout << " istart= "<< istart << " ilast= " << ilast<<endl;
  cout << " ibar1= "<< ibar1 << " ibar2= " << ibar2 <<endl;
  */

  if(q[0]==1) {
    iend[1]= ibar1 ==0 ? ilast : istart;
    constq[1]=q[1];

  } else if(q[0]==0) {
    if(ibar1==0 && ibar2 !=0) {
      iend[0]=istart;
      constq[0]=0;
      if(q[1]==1) {
        iend[1]=ilast;
        constq[1]=q[1];
      }
    } else if(ibar2 == 0 && ibar1 !=0) {
      constq[0]=0;
      iend[0]=ilast;
      if(q[1]==1) {
        iend[1]=istart;
        constq[1]=q[1];
      }
    } else {
      cout << "HadronDecay::findCostQuark strange q[0]= "<< q[0]
	   << " q1= "<< q[1]
	   << " ibar1= "<< ibar1 << " ibar2= "<< ibar2
	   <<endl;
      exit(1);
    }
  } else {
    cout << "HadronDecay::findCostQuark strange q[0]= "<< q[0]<<endl;
    exit(1);
  }

  if(ibar<0) {
      std::swap(constq[0],constq[1]);
      std::swap(iend[0],iend[1]);
  }

  return;


  if(nhadron==2) {
    iend[0]=istart;
    iend[1]=event.size()-1;
    int mo1 = event[istart].mother1();
    if(event[istart].particleDataEntry().isBaryon()) {
      iend[0]=event.size()-1;
      iend[1]=istart;

      constq[0]=q[0];
      constq[1]=q[1];

    } else {
      if(event[mo1].isDiquark()) { 
        constq[1]=q[0];
        constq[0]=q[1];
      } else {
        constq[0]=q[0];
        constq[1]=q[1];
      }
    }

    if(constq[0]==0 && constq[1]==2) {
      constq[0]=1;
      constq[1]=1;
    }
    if(constq[0]==0) iend[0]=0;
    if(constq[1]==0) iend[1]=0;

    return;
  }


  if(q[0]==1) {
	constq[0]=1;
	for(int i=1; i<event.size();i++) 
	if(event[i].isFinal()) {iend[0]=i; break;}
  } else if(q[0]==0) {
	for(int i=1; i<event.size();i++) 
	if(event[i].isFinal()) {iend[0]=i; break;}
	if(event[iend[0]].particleDataEntry().isBaryon()) {
	    constq[0]=0;
	} else if(event[iend[0]+1].particleDataEntry().isBaryon()) {
	    constq[0]=0;
	} else {
	    cout << "1 HadronDecay::findConstQuark strange constq1 "<<endl;
	    event.list();
	    exit(1);
	}
  }

  if(q[1]==1) {
	iend[1]=event.size()-1;
	if(event[iend[1]].particleDataEntry().isBaryon()) {
	    constq[1]=1;
	} else if(event[iend[1]-1].particleDataEntry().isBaryon()) {
	    constq[1]=1;
	} else {
	    cout << "HadronDecay strange constq1 " << endl;
	    cout << " q0= "<< q[0] << " q1= "<< q[1] << endl;
	    event.list();
	    exit(1);
	}

  } else if(q[1]==0) {
	iend[1]=event.size()-1;
	constq[1]=0;
  }

  if(ibar<0) {
      std::swap(constq[0],constq[1]);
      std::swap(iend[0],iend[1]);
  }

  /*
    if(constq[0]==0 && constq[1]==2) {
      constq[0]=1;
      constq[1]=1;
    }
    if(constq[0]==0) iend[0]=0;
    if(constq[1]==0) iend[1]=0;
    */

  bool ok=true;
  if(nconstq != constq[0]+constq[1]+constq[2]) ok=false;
  //if(iend[0]==iend[1]) ok=false;
  //if(iend[0]==iend[2]) ok=false;
  //if(iend[1]==iend[2]) ok=false;

  if(ok==false) {
    cout << "HadronDecay const quark wrong number!" << endl;
    cout << " q0= "<< q[0]
         << " q1= "<< q[1]
	 <<endl;
    cout << " nconstq= " << nconstq
         << " q1= " << constq[0] << " iend= "<< iend[0]
	 << " q2= " << constq[1] << " iend= "<< iend[1]
	 << " q3= " << constq[2] << " iend= "<< iend[2]
	 <<endl;
     event.list();
     exit(1);
  }


}

} // end namespace jam2
