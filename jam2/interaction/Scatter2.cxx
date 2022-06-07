#include <jam2/interaction/Scatter2.h>
#include <Pythia8/PythiaStdlib.h>
#include <jam2/hadrons/JamStdlib.h>
#include <Pythia8/PythiaStdlib.h>

namespace jam2 {

using namespace std;
using Pythia8::pow3;
using Pythia8::pow4;

Scatter2::Scatter2(Pythia8::Info* inf,Pythia8::Settings* s,JamParticleData* jd,
	CrossSection* xs, Pythia8::Rndm* r): ScatterKin(inf,s,jd,xs,r)
{
  optSuppressSoftCollision= settings->mode("Cascade:optSuppressSoftCollision");
  paramSoftColl= settings->parm("Cascade:paramSoftCollision");
  iBack=settings->mode("Cascade:optElasticBackWardScattering");
}

Scatter2::~Scatter2() 
{
}

void Scatter2::absorb(vector<EventParticle*>& outgoing)
{
  nOutParticle=1;

  pd[0]=jamParticleData->find(idOut[0]);
  pidOut[0]=jamParticleData->pid(idOut[0]);
  outGoing[0] = new EventParticle(idOut[0],pd[0]);
  outGoing[0]->setPID(pidOut[0]);

  Vec4 pcmf = pIncoming[0]/SM[0] + pIncoming[1]/SM[1];
  //Vec4 pcmf = pIncoming[0] + pIncoming[1];
  mOut[0] = pcmf.mCalc();
  pOut[0] = pIncoming[0] + pIncoming[1];
  mOutEff[0] = pOut[0].mCalc();

  // converted into a canonical momentum.
  if(optPropagate==1) {
    pOut[0] += inComing[0]->potv() + inComing[1]->potv();
    pOut[0].e( sqrt(mOutEff[0]*mOutEff[0]+pOut[0].pAbs2()) );
  }

  // Check minimum resonance mass in case potential is on.
  //if(inComing[0]->meanFieldOn() || inComing[1]->meanFieldOn())
  if(abs(mOut[0]-mOutEff[0]) > 1e-5)
    mOut[0] = max(pd[0]->m0Min(),  Pythia8::m(pIncoming[0],pIncoming[1]) );

  outGoing[0]->setMass(mOut[0]);
  outGoing[0]->setMomentum(pOut[0]);
  outGoing[0]->setOnShell();

  int qmax=2;
  if(outGoing[0]->baryon()==0) qmax=1;
  int q[2];
  q[0]= min(1,constQ[0][0]+constQ[1][0]);
  q[1]= min(qmax,constQ[0][1]+constQ[1][1]);
  outGoing[0]->setConstQuark(q);
  bool formed = false;
  if(outGoing[0]->baryon() == 0) {
    if(q[0]+q[1]==2)  formed = true;
  } else {
    if(q[0]+q[1]==3)  formed = true;
  }

  int jt=0;
  double ftime = max(inComing[0]->getTf(),inComing[1]->getTf());
  if(inComing[0]->baryon() !=0) {
    outGoing[0]->setCoordinate(xOut[0]);
    outGoing[0]->setVertex(xOut[0]);
    if(formed) ftime = xOut[0].e();
    //ftime = inComing[0]->getTf();
  } else {
    jt=1;
    outGoing[0]->setCoordinate(xOut[1]);
    outGoing[0]->setVertex(xOut[1]);
    if(formed) ftime = xOut[1].e();
    //ftime = inComing[1]->getTf();
  }

 // if(inComing[jt]->meanFieldOn())  {
    //outGoing[0]->setPotS(sPot[jt]);
    outGoing[0]->setPotV(vPot[jt]);
    outGoing[0]->setPotS(mOutEff[0]-mOut[0]);
  //}

  // Store outgoing particles.

  outGoing[0]->setNumberOfColl(numCollision);
  outGoing[0]->lastColl(nOperation);
  outGoing[0]->setBoostMatrix(MfromCM);
  outGoing[0]->setParent(ABSORPTION);


    double meff=outGoing[0]->mCalc();
    if(abs(meff-mOutEff[0])>1e-9) {
	cout << "absorb meff?  meff= "<< meff
	    << " mOutEff= "<< mOutEff[0]
	    << " spot1 "<< sPot[0]
	    << " spot2 "<< sPot[1]
	    << " id1 "<< idIn[0]
	    << " +id2 -> "<< idIn[1]
	    << " id3 "<< idOut[0]
	    << endl;
	exit(1);
    }

  outGoing[0]->setFormationTime(ftime);
  double m=outGoing[0]->getMass();
  double e=outGoing[0]->getE0();
  double dect = jamParticleData->lifeTime(pd[0],m,e);
  outGoing[0]->setLifeTime(xOut[0].e()+dect);
  //outGoing[0]->setLifeTime(ftime+dect);
  outgoing.push_back(outGoing[0]);

}

void Scatter2::scatter2(TwoBodyInterList* inter,Collision* event,vector<EventParticle*>& outgoing)
{
  double phi0=0.0;
  if(!optPreserveReactionPlane) phi0=2*M_PI*rndm->flat();

  // Find scattering angle.
  double cth=0.0;
  if(channel==ELASTIC) {
    outGoing[0] = new EventParticle(idIn[0],inComing[0]->getParticleDataEntry());
    outGoing[1] = new EventParticle(idIn[1],inComing[1]->getParticleDataEntry());

    outGoing[0]->setPID(inComing[0]->getPID());
    outGoing[1]->setPID(inComing[1]->getPID());

    mOutEff[0]=mInEff[0];
    mOutEff[1]=mInEff[1];
    //mOutEff[0]=pIncoming[0].mCalc();
    //mOutEff[1]=pIncoming[1].mCalc();

    pFinal=PCM(eCM,mOutEff[0],mOutEff[1]);
    cth = elasticAngle(pIni);
    outGoing[0]->setConstQuark(constQ[0]);
    outGoing[1]->setConstQuark(constQ[1]);

  } else {

    pd[0]=jamParticleData->find(idOut[0]);
    pd[1]=jamParticleData->find(idOut[1]);

    outGoing[0] = new EventParticle(idOut[0],pd[0]);
    outGoing[1] = new EventParticle(idOut[1],pd[1]);
    pidOut[0]=jamParticleData->pid(idOut[0]);
    pidOut[1]=jamParticleData->pid(idOut[1]);
    outGoing[0]->setPID(pidOut[0]);
    outGoing[1]->setPID(pidOut[1]);
    mOutEff[0]=mOut[0]*SM[0];
    mOutEff[1]=mOut[1]*SM[1];

    if(eCM < mOutEff[0]+mOutEff[1]) {
          cout << " eCM = " << eCM
   	   << " moutEff1= "<< mOutEff[0]
	   << " moutEff2= "<< mOutEff[1]
	   << " m1+m2= "<< mOutEff[0]+mOutEff[1]
	   << " mInEff0= "<< mInEff[0]
	   << " mInEff1= "<< mInEff[1]
	   <<endl;
	  exit(1);
    }

    pFinal=PCM(eCM,mOutEff[0],mOutEff[1]);
    cth = inelasticAngle(pIni,pFinal);

    if(inComing[0]->baryon() == outGoing[0]->baryon() ) {
      outGoing[0]->setConstQuark(constQ[0]);
      outGoing[1]->setConstQuark(constQ[1]);
    } else {
      outGoing[0]->setConstQuark(constQ[1]);
      outGoing[1]->setConstQuark(constQ[0]);
    }
  }


  pOut[0]=0.0;
  pOut[1]=0.0;

  double s1 = mOutEff[0]*mOutEff[0];
  double s2 = mOutEff[1]*mOutEff[1];

  double e1 = sqrt(s1 + pFinal*pFinal);
  double e2 = sqrt(s2 + pFinal*pFinal);
  double pz1 =  pFinal;
  double pz2 = -pFinal;

  /*
  pOut[0].e( 0.5 * (eCM + ( mOutEff[0]*mOutEff[0] - mOutEff[1]*mOutEff[1] ) / eCM ) );
  pOut[0].pz(sqrt(max(0.0,pOut[0][0]*pOut[0][0]-mOutEff[0]*mOutEff[0])));
  pOut[1].e(eCM - pOut[0][0]);
  pOut[1].pz(-pOut[0][3]);

  pOut[0].e(sqrt(mOutEff[0]*mOutEff[0]+pFinal*pFinal));
  pOut[0].pz(pFinal);
  pOut[1].e(sqrt(mOutEff[1]*mOutEff[1]+pFinal*pFinal));
  pOut[1].pz(-pFinal);
  */

  pOut[0].e( e1 );
  pOut[0].pz( pz1 );
  pOut[1].e( e2 );
  pOut[1].pz( pz2 );


//...Rotate outgoing partons using cos(theta)=(th-uh)/lam(sh,sqm3,sqm4)
  pOut[0].rot(acos(cth),phi0);
  pOut[1].rot(acos(cth),phi0);
  pz1 = pOut[0].pz();
  pz2 = pOut[1].pz();

  // back to the original frame.
  for(int i=0; i<nOutParticle; i++) {
    //pOut[i].rot(thetaCM,phiCM); 
    //pOut[i].bst(beX,beY,beZ);
    pOut[i].rotbst(MfromCM);
  }

  // Quantum mechanical suppression of very soft multiple collision.
  if(optSuppressSoftCollision) {
    //if(channel==ELASTIC) {
    if(checkSoftCollisionTime(e1,e2,pz1,pz2,inter,event)) return;
  }

  // convert to cannonical momentum.
  if(optPropagate==1) {
    pOut[0] += inComing[0]->potv();
    pOut[0].e( sqrt(mOutEff[0]*mOutEff[0]+pOut[0].pAbs2()) );
    pOut[1] += inComing[1]->potv();
    pOut[1].e( sqrt(mOutEff[1]*mOutEff[1]+pOut[1].pAbs2()) );
  }

  // Store outgoing particles.
  for (int jt=0; jt<nOutParticle; jt++) {
    //pOut[jt][0] -= vPot[jt][0];
    outGoing[jt]->setCoordinate(xOut[jt]);
    outGoing[jt]->setVertex(xOut[jt]);
    outGoing[jt]->setMass(mOut[jt]);
    outGoing[jt]->setMomentum(pOut[jt]);

    //outGoing[jt]->setPotS(mOutEff[jt]-mOut[jt]);
    outGoing[jt]->setPotS(sPot[jt]);
    outGoing[jt]->setPotV(vPot[jt]);

    outGoing[jt]->setNumberOfColl(numCollision);
    outGoing[jt]->lastColl(nOperation);
    outGoing[jt]->setBoostMatrix(MfromCM);
    outGoing[jt]->setParent(channel);

    double ftime=inComing[jt]->getTf();
    outGoing[jt]->setFormationTime(ftime);

    if(channel !=ELASTIC) {
      double m=outGoing[jt]->getMass();
      double e=outGoing[jt]->getE0();
      double dect = jamParticleData->lifeTime(pd[jt],m,e);
      outGoing[jt]->setLifeTime(xOut[jt].e()+dect);
      //outGoing[jt]->setLifeTime(ftime+dect);
    }
    outgoing.push_back(outGoing[jt]);
  }

}

bool Scatter2::checkSoftCollisionTime(double e1, double e2, double pz1, double pz2,
    TwoBodyInterList* inter,Collision* event)
{
    /*
    Vec4 x1=xOut[0], x2=xOut[1];
    x1.bst(-beX,-beY,-beZ); x1.rot(0.0,-phiCM); x1.rot(-thetaCM,0.0); 
    x2.bst(-beX,-beY,-beZ); x2.rot(0.0,-phiCM); x2.rot(-thetaCM,0.0); 
    */

    /*
    MtoCM = MfromCM;
    MtoCM.invert();
    Vec4 x1=xOut[0], x2=xOut[1];
    x1.rotbst(MtoCM);
    x2.rotbst(MtoCM);
    cout << " x1= "<< x1;
    cout << " x2= "<< x2;
    Vec4 p1=pIncoming[0], p2=pIncoming[1];
    p1.rotbst(MtoCM);
    p2.rotbst(MtoCM);
    cout << " p1= "<< p1;
    cout << " p2= "<< p2;
    */


  // p^- and p^+ before collision.
  double pm01=sqrt(mInEff[0]*mInEff[0]+pIni*pIni) - pIni;  // p-
  double pp02=sqrt(mInEff[1]*mInEff[1]+pIni*pIni) - pIni;  // p+

  // light-cone momentum after scatter.
  double pm1 = e1 - pz1;
  double pp2 = e2 + pz2;

  // collision time in the rest frame of particle.
  //double gam1=pIncoming[0][4]/mInEff[0];
  //double gam2=pIncoming[1][4]/mInEff[1];
  double gamCM=pCM[0]/eCM;
  //double dtau1=paramSoftColl*HBARC/max(1e-8,abs(pm01-pm1))/gam1;
  //double dtau2=paramSoftColl*HBARC/max(1e-8,abs(pp02-pp2))/gam2;
  double dtau1=paramSoftColl*HBARC/max(1e-8,abs(pm01-pm1))*gamCM;
  double dtau2=paramSoftColl*HBARC/max(1e-8,abs(pp02-pp2))*gamCM;
  double tc1= xOut[0][0] + dtau1;
  double tc2= xOut[1][0] + dtau2;

  // collision time in the computational frame.
  //double gam3=pOut[0].e()/mOutEff[0];
  //double gam4=pOut[1].e()/mOutEff[1];
  //double tc1= xOut[0][0] + dtau1*gam3;
  //double tc2= xOut[1][0] + dtau2*gam4;
  //double tcnext1= (inComing[0]->nextCollisionTime()-xOut[0][0])/gam1;
  //double tcnext2= (inComing[1]->nextCollisionTime()-xOut[1][0])/gam2;

  //double gg=1.0/sqrt(1.0 - beX*beX-beY*beY-beZ*beZ);
  /*
  if( (inComing[0]->nextCollisionTime() < tc1)
    || ( inComing[1]->nextCollisionTime() < tc2)) {
    cout << " t01= "<< xIncoming[0][0] << " t02= "<< xIncoming[1][0] <<endl;
    cout << " t 1= "<< xOut[0][0] << " t 2= "<< xOut[1][0]<<endl;
    cout << " tc1= "<< tc1 << " nextcol= " << inComing[0]->nextCollisionTime() <<endl;
    cout << " dtau1= " << dtau1 << " gam1= " << gam1 << " gam3= "<< gam3 << " gamcm= "<< gamCM 
       << " t1= "<< xOut[0][0] <<endl;
    cout << "coll time1= "<< endl;
    inComing[0]->printCollisionTime();


    cout << " id1= " << idIn[0] << " id2= " << idIn[1] << " id3= "<< idOut[0] << " id4= "<< idOut[1]<<endl;
    cout << " tc2= "<< tc2 << " nextcol= " << inComing[1]->nextCollisionTime() <<endl;
    cout << " dtau2= " << dtau2 << " gam2= " << gam2 << " gam4= "<< gam4 << " gamcm= "<< gamCM
      << " t2= "<< xOut[1][0]  <<endl;
    cout << "coll time2= "<< endl;
    inComing[1]->printCollisionTime();
    cout << " id1= " << idIn[0] << " id2= " << idIn[1] << " id3= "<< idOut[0] << " id4= "<< idOut[1]<<endl;
    cin.get();
  }
  */

  /*
  if(!preHadronB && inComing[1]->nextCollisionTime() < tc2) {
    cout << " tc2= "<< tc2 << " nextcol= " << inComing[1]->nextCollisionTime() <<endl;
    cout << " dtau2= " << dtau2 << " gam2= " << gam2 << " gam4= "<< gam4 << " gamcm= "<< gamCM
      << " t2= "<< xOut[1][0]  <<endl;
    cout << " id1= " << idIn[0] << " id2= " << idIn[1] << " id3= "<< idOut[0] << " id4= "<< idOut[1]<<endl;
    cin.get();
  }
  */

  if(event->checkNextCollisionTime(inter,tc1,tc2,preHadronA,preHadronB)) return true;
  //if(!preHadronA && inComing[0]->nextCollisionTime() < tc1) return;
  //if(!preHadronB && inComing[1]->nextCollisionTime() < tc2) return;

  return false;

}

/*
void Scatter2::scatter2to2(double phi, double cth)
{
    if(mOut[0] + mOut[1] >= eCM) {
	cerr << "Scatter2::scatter2to2 error emf1 + emf2 > srt " << eCM
	     << " m1= " << mOut[0]
	     << " m2= " << mOut[1]
	     << " m1+m2= " << mOut[0]+mOut[1]
	     << endl;
	exit(1);
    }

    pOut[0]=0.0;
    pOut[1]=0.0;
    pOut[0].e(0.5*(eCM+(mOut[0]*mOut[0] - mOut[1]*mOut[1])/eCM));
    pOut[0].pz(sqrt(max(0.0,pOut[0][0]*pOut[0][0]-mOut[0]*mOut[0])));
    //pFinal=pOut[0].pz();

    pOut[1].e(eCM - pOut[0][0]);
    pOut[1].pz(-pOut[0][3]);

//...Rotate outgoing partons using cos(theta)=(th-uh)/lam(sh,sqm3,sqm4)
    pOut[0].rot(acos(cth),phi);
    pOut[1].rot(acos(cth),phi);

}

// back to the original frame.
void Scatter2::boostBack(vector<EventParticle*>& outgoing)
{
    double tm0 = min(xOut[0].e(),xOut[1].e());
    Vec4 ptot=0.0;
    int nop = outgoing.size();
    for(int i=0; i<(int)outgoing.size(); i++) {
	Vec4 p = outgoing[i]->getMomentum();
	p.rot(thetaCM,phiCM); p.bst(beX,beY,beZ);
	outgoing[i]->setMomentum(p);

	if(nop > 2) {
	Vec4 x = outgoing[i]->getR();
	x.rot(thetaCM,phiCM); x.bst(beX,beY,beZ);
	if(x.e() < 0 )
	    cerr << "after boost t<0 " << i << " t = " << x.e() << endl;
	x += xCM;
	outgoing[i]->setCoordinate(x);
	//} else {
	 //   outgoing[i]->setCoordinate(xIncoming[i]);
	}
	if(outgoing[i]->getT() < tm0) {
	    cerr << "backforwd collision??? tm0="<< tm0
		 << " tf= " << outgoing[i]->getT()
		 << " tcm= " << xCM.e()
		 << endl;
	}
	ptot += p;
    }
}
*/

// gives cos(th) for the distribution 1/(t-debyemass**2)**2
double Scatter2::getThat(const double srt, const double em1,
       	const double em2, const double em3, const double em4,
       	const double debyemass)
{
    double pp1 = pawt(srt,em1,em2);
    double pp2 = pawt(srt,em3,em4);
    double x = (em1*em1 - em3*em3 - em2*em2 + em4*em4)/(2*srt);
    double ppx1 = (pp1 + pp2)*(pp1 + pp2);
    double ppx2 = (pp1 - pp2)*(pp1 - pp2);
    double t0 = x*x - ppx1 - debyemass*debyemass;
    double t1 = x*x - ppx2 - debyemass*debyemass;

    double t = 1.0/(rndm->flat()*(1.0/t1 -1.0/t0) + 1.0/t0 )
	+debyemass*debyemass;

    double e1  = sqrt(em1*em1+pp1*pp1);
    double e3  = sqrt(em3*em3+pp2*pp2);
    double cth = (t-em1*em1 - em3*em3 +2*e1*e3)/(2*pp1*pp2);
    cth = rndm->flat() < 0.5 ? cth : -cth;
    return cth;
}

void Scatter2::print(ostream& ofs,vector<EventParticle*>& outgoing) const
{

    ofs << inComing[0]->getParticleDataEntry()->name()
	<< " + "
        << inComing[1]->getParticleDataEntry()->name()
	<< " -> ";

	ofs << outgoing[0]->getParticleDataEntry()->name()
	    <<  " + " 
	    << outgoing[1]->getParticleDataEntry()->name()
	    <<endl;

    for(int i=0; i< (int)outgoing.size(); i++) {
	outgoing[i]->print(ofs);
    }
}


//******************************************************************

double Scatter2::elasticAngle(double pr)
{
//...Purpose: to determine elastic scattering angle.
//=================================================================*
//...determine the new momentum direction direction
//...(relative to the original direction).
//...calculate the cosine of the polar angle, c1.
//
//        pr     : magnitude of rel. momentum.      (input)
//        sig    : total cross section.             (input)
//        t1     : azimuthal angle                  (output)
//        c1     : cosine(polar angle)              (output)
//
//=================================================================*

    static const double emnuc=0.938;

//...Elastic scattering of identical particles
    double factor=1;
    double c1=0.0;

//...Elastic scattering of non-identical particles
    int kf1=cpair.getID(0);
    int kf2=cpair.getID(1);
    if(kf1 != kf2) factor=2;

    //double pr=cpair.getCMmomentum();
    double sh=4.0*(emnuc*emnuc+pr*pr);
    double snew=sqrt(sh);
    double plab=sqrt(sh*(sh-4.0*emnuc*emnuc))/(2.0*emnuc);

//....Get slope parameter
    int ibar1=cpair.baryon(0);
    int ibar2=cpair.baryon(1);
    double a=elasticSlope(kf1,kf2,ibar1,ibar2,plab,snew);
    double sig=cpair.getSigma();
    double fac=sig/40.0;
    if((kf1==2112 || kf1==2212) && (kf2 == 2112 || kf2 == 2212)) fac=1.0;
    a *=fac;
    double  ta=2.0*pr*pr;
    double  ata=ta*a;
    if(ata < 1.e-7) {
        return 1.0-2.0*rndm->flat();
    } else {
        double y3=rndm->flat();
        double tt1=log((1.0-y3)*exp(-min(50.0,factor*ata))+y3)/ata;
        c1=1.0+tt1;
    }

//...Backward scattering prob.
  if(iBack >= 1) {
    if(cpair.getCollType() ==  1 && kf1 != kf2 && snew <= 3.0) {
      double bprob=1.0;
      if(plab > 0.8) {
	bprob=0.8/plab; // Cugnon
	if(iBack==2) bprob = pow2(bprob); // Niita
        bprob=bprob/(1.0+bprob);
      }
      if(iBack < 3 && cpair.getTotalStrangeness()==0) {
        if(rndm->flat() <= bprob) c1 *= -1;
      }
    }
  } else if(iBack==-1) {
      c1 *= -1; 
  }
 
  if(abs(c1) > 1.0) c1= c1 > 0 ? 1.0 : -1.0; 
  return c1;
}

//***********************************************************************

double Scatter2::elasticSlope(int kf1,int kf2,int ibar1,int ibar2,double plab,
	double snew)
{
    if(snew >= 10.0) {
        double b1=2.3;
        double b2=2.3;
        if(ibar1==0) {
          b1=1.4;
          int kfl1=(abs(kf1)/100) % 10;
          int kfl2=(abs(kf1)/10)  % 10;
          if(10*kfl1+kfl2 == 44) b1=0.23;  // J/psi
	}
        if(ibar2==0) {
          b2=1.4;
          int kfl1=(abs(kf2)/100) % 10;
          int kfl2=(abs(kf2)/10)  % 10;
          if(10*kfl1+kfl2 == 44) b2=0.23;  // J/psi
	}
        return 2.0*(b1+b2) + 4 * pow(snew*snew,0.0808)-4.2;
    }

//...Anti-p p
    if(ibar1*ibar2 == -9) {
//...M.R.Clover,et al., Phys. Rev. C26 (1982) 2138.
//...This fit is valid from 100MeV to 2GeV lab. eng.
//       a=12.94+39.03*exp(-2.075*plab)
//...J. Cugnon, et al., Phys. Rev. C41 (1990) 1701.
	return 3.34;
    }

    double a;
//...p+n
    if(kf1 != kf2) {

	if( plab <=  0.6 )
            a = 6.2 * ( plab - 0.225 )/ 0.375;
        else if( plab <= 1.6 )
            a =  - 1.63 * plab + 7.16;
        else if( plab <= 2.0 )
            a = 5.5 * pow(plab,8) / ( 7.7 + pow(plab,8) );
        else if(plab <= 4.2)
            a = 5.34 + 0.67 * ( plab - 2.0 );
        else
//          a = 5.656 + 0.863*log(plab)
	    a = 5.656 + 1.300*log(plab);

    } else {

          if( plab <= 2.0 )
            a = 5.5 * pow(plab,8) / ( 7.7 + pow(plab,8) );
          else if(plab <= 4.2)
            a = 5.34 + 0.67 * ( plab - 2.0 );
          else
//           a = 5.656 + 0.863*log(plab);
//......srt=19GeV
//           a = 5.656 + 1.100*log(plab);
//......srt=10GeV
	    a = 5.656 + 1.300*log(plab);

    }

    return a;
}

//******************************************************************

double Scatter2::inelasticAngle(double pr0, double pr)
{
//...Purpose: to determine inelastic scattering angle.

//        srt    : c.m. energy                              (input)*
//        pr0    : magnitude of initial rel. momentum       (input)*
//        pr     : magnitude of final rel. momentum         (input)*
//        em1    : mass of the particle 1                   (input)*
//        em2    : mass of the particle 2                   (input)*
//        ibar1  : baryon number of the particle 1          (input)*
//        ibar2  : baryon number of the particle 2          (input)*
//        kf01   : KF code  of the initial  particle 1      (input)*
//        kf02   : KF code  of the initial  particle 2      (input)*
//        kc01   : KC code  of the initial  particle 1      (input)*
//        kc02   : KC code  of the initial  particle 2      (input)*
//        t1     : azimuthal angle                         (output)*
//        c1     : cosine(polar angle)                     (output)*

    //const int mstc67=2, mstc68=2, mstc69=2, mstc79=2;

    int ibar1=cpair.baryon(0);
    int ibar2=cpair.baryon(1);
    double srt=eCM;
    //double pr0=cpair.getCMmomentum();
    //double pr=pFinal;
    //double em1=mOut[0];
    //double em2=mOut[1];

//...BB collisions.
    if(ibar1*ibar2 == 9) {
 
//.....No strange BB collisions.
        int str1=cpair.getStrange(0);
        int str2=cpair.getStrange(1);
        if(abs(str1)+abs(str2) == 0 && eCM < 2.17) {
             return angDeltaN();
	}

        return angInel(srt,pr0,pr);

//...MM
    } else if(ibar1 == 0 && ibar2 == 0) {
           return angInel(srt,pr0,pr);
//...MB
    } else if(ibar1*ibar2 == 0) {
           return angInel(srt,pr0,pr);
//...Ani-B B
    } else if(ibar1*ibar2 == -9) {
           return angInel(srt,pr0,pr);
    } else {
        cout << "(jamangin:) unrecognize baryon number " <<endl;
	exit(1);
    }

}

//***********************************************************************
// part of subroutine jamangrr(pr0,pr,srt,c1,em1,em2,iang)
double Scatter2::angInel(double srt, double pr0, double pr)
{
    double a = 2.5 + max(0.0,0.7*(log(srt*srt)-log(2.0)));
    double ata=4.0*pr0*pr*a;
    if(ata < 1.e-7) return 1.0;
    double xran=rndm->flat();
    double c1=1.0+log((1.0-xran)*exp(-min(50.0,ata))+xran)/ata;
    if(abs(c1) > 1.0) c1= c1 > 0 ? 1.0 : -1.0;
    //sign(1.0d0,c1)

    return c1;
}

/*
      subroutine jamangrr(pr0,pr,srt,c1,em1,em2,iang)

//...Purpse: to generate angular distribution for hh interactions.
//     input: pr0 : c.m. mom. before collision
//            pr  : c.m. mom. after  collision
//            srt : inv. mass of collision
//     output: c1 : scattering angle in c.m.
//======================================================================c
//     Method:
// (1) Follwoing HIJING, 
//         prob(pt^2)=1.0/(x**2+hip0**2)/(x**2+hic**2)
//    &                  /(1+exp((x-hip0)/0.4))
//     with hic=0.1, hip0=1.0
//     Since it is difficult to generate this pt distribution, 
//     by using Monte-Carlo Method (without Table), 
//     we fit the above function as follows,
//       prob(pt^2)=a*(
//    &             exp(-x**2/b**2)/b/b		! Gauss
//    &           c*exp(-x**2/d**2)/d/d		! Gauss
//    &          +e*(exp(-x/f)-exp(-2*x/f))/x/f)	! 2-Exp.
//     Then, each part can be generated by 
//       pt=b*sqrt(-log(1.e0-rnp*(1.0-exp(-x0**2/b**2)))) ! Gauss
//     or
//       pt=-f*log(1.0-sqrt(rn(0))*(1.0-exp(-x0/f))) ! 2-Exp.
//     Selection of each part is done by Monte-Carlo, therefore,
//     This procedure uses two random numbers
//
// (2) Consideration in small pr
//     The above pt^2 distribution of HIJING can be used around 24 GeV/c.
//     At lower energies, we considered as follows.
//         d(sigma)/d(p_t) = d(sigma)/d(theta)/pr
//     where theta is the CM scattering angle (< pi/2).
//     This is true for large pr.
//     Normalization change due to the finite range of pt is ignored.
//
//     Fitting part is done by A.Ohnishi
//======================================================================c

//...Original Parameter Set
    const double ag1=1.0,ag2=1.60207,ag3=5.63708;
    const double bg1=0.0944812,bg2=0.196705,bg3=0.288274;
//...Parameter set II:
//     parameter(ag1=1.0d0,ag2=0.78d0,ag3=2.78d0)
//     parameter(bg1=0.11d0,bg2=0.28d0,bg3=0.29d0)

//...Ohnishi
//     parameter(ag1=1.0d0,ag2=1.0d0,ag3=2.5d0)
//     parameter(bg1=0.15d0,bg2=0.32d0,bg3=0.35d0)
//...Artificial Parameter Set: 
//     parameter(ag1=0.0d0,ag2=0.0d0,ag3=1.00d0)
//     parameter(bg1=0.11d0,bg2=0.30d0,bg3=0.32d0)

    if(iang == 1) {

//...Select the Functional Form of Transverse Momentum Distribution
//...pt=pr*theta, and theta is limited to be < pi/2
        ptx=pr*paru(1)/2
        rnx=rn(0)*(ag1+ag2+ag3)
        rnp=rn(0)
                if(rnx.le.ag3) then
        expf3=1.0d0-exp(-ptx/bg3)
        pt=-bg3*log(1.0d0-sqrt(rnp)*expf3)
                elseif(rnx.le.ag3+ag2) then
        expf2=1.0d0-exp(-ptx**2/bg2**2)
        pt=bg2*sqrt(-log(1.d0-rnp*expf2))
                else
        expf1=1.0d0-exp(-ptx**2/bg1**2)
        pt=bg1*sqrt(-log(1.d0-rnp*expf1))
                endif

//...In the original (HIJING) treatment, pt is real transverse momentum.
//...However, at lower energies, it easily becomes larger than pr,
//...and the measure d(pt^2) does not coincide with the normal measure
//...pr^2 d(cos(theta)), we re-interpret d(sigma)/d(pt)
//...as d(sigma)/d(theta)/pr at high energies and forward angles,
//...these become the same.
        c1=cos(pt/pr)

    } else if(iang == 2)  {

//...See for example, P.T.P.suppl. 41 and 42(1967) p291
//       a=parc(44)
//...Diffraction dissosiation
//       a=7.72+max(0.0d0,0.6d0*(log(srt*srt)-log(2.d0)))
//...binary diagram
        a = 2.5 + max(0.0,0.7*(log(srt*srt)-log(2.0)))

        ata=4.0*pr0*pr*a
        if(ata < 1.e-7) return 1.0;

        double xran=rndm->flat();
        c1=1.0+log((1.0-xran)*exp(-min(50.0,ata))+xran)/ata;
        if(abs(c1) > 1.0) c1= c1 > 0 ? 1.0 : -1.0;
	//sign(1.0d0,c1)

    } else if(iang == 3) {

        b=2.5;  // parc(44)
        c=0.7;  //parc(45)
        sqm1=pcp(5,1)**2
        sqm2=pcp(5,2)**2
        sqm3=em1**2
        sqm4=em2**2
        sh=srt*srt

c...Determine maximum possible t range and coefficients of generation.
        sqla12=(sh-sqm1-sqm2)**2-4.d0*sqm1*sqm2
        sqla34=(sh-sqm3-sqm4)**2-4.d0*sqm3*sqm4
        tha=sh-(sqm1+sqm2+sqm3+sqm4)+(sqm1-sqm2)*(sqm3-sqm4)/sh
        thb=sqrt(max(0d0,sqla12))*sqrt(max(0d0,sqla34))/sh
        thc=(sqm3-sqm1)*(sqm4-sqm2)+(sqm1+sqm4-sqm2-sqm3)*
     &  (sqm1*sqm4-sqm2*sqm3)/sh
        thl=-0.5d0*(tha+thb)
        thu=thc/thl
        if(c.gt.0.0d0) then
          icpos=1
          thl=max(thl,-0.5d0*b/c)
          thrnd=exp(max(-50.d0,0.5d0*b*(thl-thu)))-1.d0
          thwmx=exp(max(-50.d0,0.5d0*b*thu+c*thu**2))
        else
          icpos=0
          thrnd=exp(max(-50.d0,b*(thl-thu)))-1.d0
          thwmx=exp(max(-50.d0,c*thu**2))
        endif

C...Select t according to exp(B*t + C*t^2).
        itry=0
  140   if(icpos.eq.1) then
          itry1=0
  144     th=thu+(2.d0/b)*log(1.+thrnd*rn(0))
          itry1=itry1+1
          if(exp(max(-50.d0,0.5d0*b*th+c*th**2)).lt.rn(0)*thwmx)
     $          goto 144
          if(itry1.gt.100) goto 147
        else
          itry2=0
  146     th=thu+(1.d0/b)*log(1.d0+thrnd*rn(0))
          itry2=itry2+1
          if(exp(max(-50.d0,c*th**2)).lt.rn(0)*thwmx) goto 146
          if(itry2.gt.100) goto 147
        endif

  147    itry=itry+1
         if(itry.ge.200) then
           ata=4.d0*pr0*pr*a
           if(ata.lt.1.d-7) then
             c1=1.d0
             return
           else
             xran=rn(0)
             c1=1.d0+log((1.d0-xran)*exp(-min(50.0d0,ata))+xran)/ata
           end if
           if(abs(c1) .gt. 1.0d0) c1=sign(1.0d0,c1)
           goto 148
         endif

        if(th.lt.thl.or.th.gt.thu) goto 140
  148   c1=(th-sqm1-sqm3+2.d0*sqrt(sqm1+pr0**2)*sqrt(sqm3+pr**2))/
     $        (2.d0*pr0*pr)

//...Gaussian dist.
    } else if(iang == 4) {
        ptsq=-parc(47)**2*log(1.d0-rn(0)*(1.0d0-exp(-(pr/parc(47))**2)))
        c1=cos(sqrt(ptsq)/pr)

//...Gaussian + exponential dist.
    }  else if(iang == 5) {

        ptx=pr*paru(1)/2
        wg1=parc(47)
        wg2=parc(48)
        rnx=rn(0)*(1.0d0+parc(49))
        if(rnx.le.parc(49)) then
          pt=-wg2*log(1.0d0-rn(0)*(1.0d0-exp(-ptx/wg2)))
        else
          pt=wg1*sqrt(-log(1.d0-rn(0)*(1.0d0-exp(-(ptx/wg1)**2))))
        endif
        c1=cos(pt/pr)

    }  else if(iang == 6) {
        ptmx=pr*paru(1)/2
        ptmx2=ptmx*ptmx
        pt=sqrt(hirnd2(8,0.0d0,ptmx2))
//       pt=jamrnd2(8,0.0d0,min(6.0d0,ptmx))
        if(pt.gt.parc(70)) then
          expf=exp(-(parc(70)/parc(66))**2)
          pt=parc(66)*sqrt(-log(expf-rn(0)
     $             *(expf-exp(-ptmx2/parc(66)**2))))
        endif
        c1=cos(pt/pr)
    }

}

*/

//******************************************************************

double Scatter2::angDeltaN()
{
//...Purpose: to determine inelastic scattering angle for nn->nd.
//-----------------------------------------------------------------*
//   determine the new momentum direction                          *
//   for inelastic scattering (relative to the original direction).*
//   return the azimuthal angle, t1,                               *
//   and the cosine of the polar angle, c1.                        *
//                                                                 *
//        srt    : c.m. energy                              (input)*
//        pr     : magnitude of rel. momentum               (input)*
//        em1    : mass of the particle 1                   (input)*
//        em2    : mass of the particle 2                   (input)*
//        kc1    : KC code of the outgoing particle 1       (input)*
//        kc2    : KC code of the outgoing particle 2       (input)*
//        c1     : cosine(polar angle)                     (output)*
//-----------------------------------------------------------------*

    double em1=cpair.M(0);
    double em2=cpair.M(1);
    double srt=eCM;
    double pr=pFinal;

    double srte=srt-em1-em2;
    double srtd=srt;
    double srtde=srte;
    double srtne=srte;
    double prad=pr;

    int id1=cpair.getPID(0);
    int id2=cpair.getPID(1);
//...Rescale except for the delta-n final state.
    if(CollisionPair::Pair(id1,id2) != 
	    CollisionPair::Pair(id_delt,id_nucl)) {
        const double em_nuc= 0.938, em_del=1.232;
        double dmasp=pd[0]->m0() + pd[1]->m0() -em_nuc - em_del;
        srtd=srt-dmasp;
	srtne=srte-dmasp;
	if(srtne > 0.0) {
          srtde=srtne;
          prad=sqrt(pow2(srtd*srtd-em1*em1-em2*em2)
                    -4.0*pow2(em1*em2))/(2.0*srtd);
	}
    }

    double c1=0.0;
    if(rndm->flat() < 0.5) {
        if(srtne > 0.0) {
          double as=pow(3.65*srtde,6);
          double a=as/(1.0+as)*pow4(srtd)*0.14;
          double ta=-2.0*prad*prad;
          double x=rndm->flat();
          double t1=log((1-x)*exp(max(-50.0,2.0*a*ta))+x)/a;
          c1=1.0-t1/ta;
          if(abs(c1) > 1.0) c1=2.0*x-1.0;
	} else {
          c1 =2.0*rndm->flat()-1.0;
	}
    } else {
        if( srtd < 2.14 ) {
           c1=2.0*rndm->flat()-1.0;
	} else {
	    double b1,b3;
            if(srtd > 2.4) {
		b1=0.06;
		b3=0.4;
	    } else {
		b1= 29.0286 - 23.749  * srtd + 4.86549 * pow2(srtd);
		b3=-30.3283 + 25.5257 * srtd - 5.30129 * pow2(srtd);
	    }
	    double pp3=b1/(3.0*b3);
	    double qq3=0.5 * (0.5 - rndm->flat()) / b3;
	    double pq3=sqrt(pow2(qq3) + pow3(pp3));
	    double uu =pow(-qq3 + pq3 ,1.0/3.0);
	    double vv =pow( qq3 + pq3 ,1.0/3.0);
	    c1 =uu-vv;
	    if(abs(c1) > 1.0) c1=c1/abs(c1);
        }
    }

    return c1;
}

} // end namespace jam2.
