#include <vector>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <iostream>

#include <jam2/collision/Collision.h>
#include <jam2/xsection/CollisionPair.h>
#include <jam2/xsection/XsecTable.h>

#include <jam2/hadrons/JamStdlib.h>
#include <Pythia8/PythiaStdlib.h>

namespace jam2 {

using namespace std;

void Collision::clearPlist()
{
  if(plist.size()>0) {
  for(auto& p : plist) delete p;
  plist.clear();
  }

  if(pnew.size()>0) {
  for(auto& p : pnew) delete p;
  pnew.clear();
  }
}


Collision::Collision(Pythia8::Settings* set, JamParticleData* jp, CrossSection* in,
	Pythia8::Rndm* r):
    settings(set), jamParticleData(jp),xsection(in), rndm(r)
{
    lastCollisionTime  = -1000;
    //maxCell     = 27; // box simulation
    numberOfCollision = 0;
    numberOfDecay = 0;
    numberOfWallCollision=0;
    minKinEnergy  = 0.001;
    predictedColl = 0;
    overSample = settings->mode("Cascade:overSample");
    optCollisionOrdering = settings->mode("Cascade:optCollisionOrder");
    optTau=false;
    if((optCollisionOrdering > 10 && optCollisionOrdering < 20) || optCollisionOrdering == 21) {
      optTau=true;
    }
    optCollisionTimeLimit = settings->mode("Cascade:optCollisionTimeLimit");
    optFluctuation=settings->mode("Cascade:optFluctuation");
    ggcf=0;
    omegaSigma=0.2;
    if(optFluctuation) ggcf = new GGCF(rndm);

    decayOn = settings->flag("Cascade:Decay");
    //finalTime = settings->parm("Cascade:FinalTime");
    double dt = settings->parm("Cascade:TimeStepSize");
    int nstep = settings->mode("Cascade:TimeStep");
    finalTime = dt*nstep;

    int  optSoftCol = settings->mode("Cascade:optSuppressSoftCollision");
    if(optSoftCol) finalTime = max(5.0,finalTime);

    constQCollisionOnly = settings->flag("Cascade:constQCollisionOnly");
    bbCollisionOnly = settings->flag("Cascade:BBCollisionOnly");
    noMMCollision = settings->flag("Cascade:noMMCollision");
    withBox=settings->mode("Cascade:box");
    xBox=settings->parm("Cascade:boxLx");
    yBox=settings->parm("Cascade:boxLy");
    zBox=settings->parm("Cascade:boxLz");

    gWidth = settings->parm("Cascade:GaussianWidth");
    int optPauli=settings->mode("Cascade:PauliBlocking");
    if(optPauli==1) {
      pauliR = 1.0/(2*gWidth);
      pauliP = 2.0*gWidth/(HBARC*HBARC);
      pauliC = 4.0;
    } else {
    // Fusimi function.
      pauliR = 1.0/(4*gWidth);
      pauliP = gWidth/(HBARC*HBARC);
      pauliC = 0.5;
    }

    passingTime= 0.0;

    rMaxCutSqNN   = 0.1*55/M_PI;
    rMaxCutSqBB   = 0.1*200/M_PI;
    rMaxCutSqMB   = 0.1*200/M_PI;
    rMaxCutSqMM   = 0.1*150/M_PI;
    //rMaxCutSqBBar = 0.1*350/M_PI;
    rMaxCutSqBBar = 0.1*200/M_PI;
    bMax=0.0;

  // Propagation by canonical momentum in RQMD mode
  optVectorPotential=settings->mode("MeanField:optVectorPotential");
  optVdot=settings->mode("MeanField:optVdot");
  optPropagate=0;
  if(optVectorPotential==1 && optVdot==0) optPropagate=1;

  cascadeMethod = settings->mode("Cascade:model");
  removeSpectator = settings->flag("Cascade:removeSpectator");

}

Collision::~Collision()
{
    list<EventParticle*>::iterator cp;
    for(cp = plist.begin(); cp != plist.end(); ++cp) {
	delete *cp;
    }
    plist.clear();

    for(cp = pnew.begin(); cp != pnew.end(); ++cp) {
	delete *cp;
    }
    pnew.clear();

    /*
    InterIt it;
    for(it = interList.begin(); it != interList.end(); ++it) {
	delete *it;
    }
    interList.clear();
    */

    delete ggcf;

    //cout << "Collision max. impact. par= "<< bMax<<endl;
}

TwoBodyInterList* Collision::hit(EventParticle* i1,EventParticle* i2)
{
  // Avoid dead particles.
  //if((*i1)->getStatus() < 0) return 0;
  //if((*i2)->getStatus() < 0) return 0;

  // Avoid first collisions within the same nucleus
  hitPath=1;
  if(i1->getNColl()*i2->getNColl() == 1 ) return nullptr;

  hitPath=2;
  // Avoid second collisions for the same pairs
  if((i1->lastColl() == i2->lastColl()) && (i2->lastColl() != -1))
    return nullptr;

  // Get collision type.
  int icltyp = collisionType(i1, i2);
  if(icltyp==0) return nullptr;

  // currently parton-hadron collision has not been implemented.
  //if(icltyp > 5 ) return nullptr;

  hitPath=3;
  // BB collision only.
  if(bbCollisionOnly && icltyp !=1) return nullptr;

  hitPath=4;
  // No meson-meson collision.
  if(noMMCollision && icltyp == 3) return nullptr;

  double ctime=0.0;
  double tcol1=0.0;
  double tcol2=0.0;
  double brel=0.0;
  double sig=0.0;

  double t1 = i1->getT();
  double t2 = i2->getT();
  Vec4 dr = i1->getR()- i2->getR();

  if(withBox) {
    dr[1] = modulo(dr.px() + xBox/2, xBox) - xBox/2;
    dr[2] = modulo(dr.py() + yBox/2, yBox) - yBox/2;
    dr[3] = modulo(dr.pz() + zBox/2, zBox) - zBox/2;
  }

  // Determine max. cross section and max. impact par.
  // as well as low energy cutoff
  Vec4 p1 = optPropagate==1 ? i1->getPkin() : i1->getP(); 
  Vec4 p2 = optPropagate==1 ? i2->getPkin() : i2->getP(); 

  double m1sq=p1.m2Calc();
  double m2sq=p2.m2Calc();
  double m1 = sqrt(m1sq);
  double m2 = sqrt(m2sq);
  double s = (p1+p2).m2Calc();
  double srt = sqrt(s);

  hitPath=5;
  // Low energy cut off.
  if(srt < m1 + m2 + minKinEnergy) return nullptr;

  hitPath=6;
  double prsq=(s-(m1+m2)*(m1+m2))*(s-(m1-m2)*(m1-m2))/(4*s);
  // Too low relative momentum.
  if(prsq < 0.000001) return nullptr;

  double drp1 = dr*p1;
  double drp2 = dr*p2;
  double dp12 = p1*p2;
  double dn = dp12*dp12 - m1sq*m2sq;
  if(dn < 1e-5) return nullptr;
  double  b12sq = drp1*drp1*m2sq + drp2*drp2*m1sq -2*drp1*drp2*dp12;
  // impact parameter squared.
  double bsq = -dr*dr - b12sq/dn;

  hitPath=9;
  if(srt < 10.0 && bsq > rMaxCutSq) return nullptr;

  // collision times in the computational frame.
  if(optCollisionOrdering < 10) {

      // Will particles get closest point in this time interval ?
      tcol1 = t1 + p1.e()*(drp1*m2sq-drp2*dp12)/dn;
      tcol2 = t2 - p2.e()*(drp2*m1sq-drp1*dp12)/dn;

      // Avoid backward collision.
      double tlast1 = i1->TimeLastColl();
      double tlast2 = i2->TimeLastColl();
      if(tcol1 <= tlast1 || tcol2 <= tlast2) return nullptr;

      // Define collision ordering time.
      switch (optCollisionOrdering) {
        case 1: ctime=0.5*(tcol1+tcol2);break;
        case 2: ctime=min(tcol1,tcol2); break;
        case 3: ctime=max(tcol1,tcol2); break;
        case 4: ctime=0.5*(tcol1+tcol2);
	      tcol1=ctime; tcol2=ctime;
              //if(tcol1 <= t1 || tcol2 <= t2) return nullptr;
              if(tcol1 <= tlast1 || tcol2 <= tlast2) return nullptr;
	      break;
        case 5: ctime=min(tcol1,tcol2);
	      tcol1=ctime; tcol2=ctime;
              //if(tcol1 <= t1 || tcol2 <= t2) return nullptr;
              if(tcol1 <= tlast1 || tcol2 <= tlast2) return nullptr;
	      break;
        default: ctime=min(tcol1,tcol2); break;
      }

      if(tcol1 > finalTime && tcol2 > finalTime) return nullptr;
      if(tcol1 < initialTime || tcol2 < initialTime) return nullptr;

      //if(optCollisionTimeLimit==1 && (ctime < t1 || ctime < t2)) return nullptr;
     // avoid collision such as t1 < t1col < t2  < t2ccol
     //if(tcol1 <= t2 || tcol2 <= t1) return nullptr;

  // Bjorken coordinate.
  } else if(optCollisionOrdering < 20) {

      double dt1 =   p1.e()*(drp1*m2sq-drp2*dp12)/dn;
      double dt2 = - p2.e()*(drp2*m1sq-drp1*dp12)/dn;
      tcol1 = t1 + dt1;
      tcol2 = t2 + dt2;

      // Avoid backward collision.
      double tlast1 = i1->TimeLastColl();
      double tlast2 = i2->TimeLastColl();
      if(tcol1 <= tlast1 || tcol2 <= tlast2) return nullptr;

      double zcol1 = i1->getZ() + p1.pz()/p1.e()*dt1;
      double zcol2 = i2->getZ() + p2.pz()/p2.e()*dt2;
      if((abs(tcol1)<abs(zcol1)) || (abs(tcol2)<abs(zcol2))) {
	cout <<"(hit) tcol1= "<< tcol1 << " zcol1= "<< zcol1
	     <<" tcol2= "<< tcol2 << " zcol2= "<< zcol2
	     <<endl;
	cout <<"t1= "<< t1 << " z1= "<< i1->getZ()
	     <<" t2= "<< t2 << " z2= "<< i2->getZ()
	     <<endl;
	cout << " srt= "<< srt << " b= "<< sqrt(bsq)
	  << " ncol1= " << i1->getNColl()
	  << " ncol2= "<< i2->getNColl()
	  <<endl;

	return nullptr;
      }

      double tau1 = sqrt(tcol1*tcol1 - zcol1*zcol1);
      double tau2 = sqrt(tcol2*tcol2 - zcol2*zcol2);

      switch (optCollisionOrdering) {
        case 11: ctime=0.5*(tau1+tau2);break;
        case 12: ctime=min(tau1,tau2); break;
        case 13: ctime=max(tau1,tau2); break;
        case 14: ctime=0.5*(tau1+tau2);
	      tcol1=0.5*(tcol1+tcol2); tcol2=tcol1;
              if(tcol1 <= tlast1 || tcol2 <= tlast2) return nullptr;
	      break;
        case 15: ctime=min(tau1,tau2);
	      tcol1=min(tcol1,tcol2); tcol2=tcol1;
              if(tcol1 <= tlast1 || tcol2 <= tlast2) return nullptr;
	      break;
        default: ctime=min(tau1,tau2); break;
      }

      if(tcol1 > finalTime && tcol2 > finalTime) return nullptr;
      if(tcol1 < initialTime || tcol2 < initialTime) return nullptr;


  // Time of collision is defined by the time of closet approach
  // in the observational (computational) frame.
  } else {

      Vec4 v1 = p1/p1.e();
      Vec4 v2 = p2/p2.e();
      Vec4 dv = v1 -v2;
      double dvsq= dv.pAbs2();
      if(dvsq < 1e-8) return nullptr;
      double dt1 =  (-dot3(dv,dr)+ dr.e()*dot3(dr,v2))/dvsq;
      // Avoid backward collision.
      if(dt1 < 0.0) return nullptr;
      double dt2 =  (-dot3(dv,dr)+ dr.e()*dot3(dr,v1))/dvsq;
      if(dt2 < 0.0) return nullptr;

      tcol1=t1 + dt1;
      tcol2=t2 + dt2;
      ctime = tcol1;
      // note that tcol1=tcol2
      if(ctime > finalTime) return nullptr;

      if(optCollisionOrdering == 21) {
        double zcol1 = i1->getZ() + p1.pz()/p1.e()*dt1;
        double zcol2 = i2->getZ() + p2.pz()/p2.e()*dt2;
        if((tcol1>zcol1) || (tcol2>zcol2)) {
  	  cout <<" tcol1= "<< tcol1 << " zcol1= "<< zcol1
	       <<" tcol2= "<< tcol2 << " zcol2= "<< zcol2
	       <<endl;
	  exit(1);
        }
        double tau1 = sqrt(tcol1*tcol1 - zcol1*zcol1);
        double tau2 = sqrt(tcol2*tcol2 - zcol2*zcol2);
	ctime = 0.5*(tau1+tau2);
      }

  }

   // Avoid collision that will happen after wall collision.
   // This predicts larger collision number for optCollisionOrdering=2,
   // and smaller collision number for optCollisionOrdering=3.
   // Others are fine.
   if(tcol1 > i1->tWall() || tcol2 > i2->tWall()) return nullptr;
   //if(ctime > i1->tWall() || ctime > i2->tWall()) return nullptr;
        
    //if(ctime < lastCollisionTime) return nullptr;
    //if(ctime < initialTime) return nullptr;

    // Check max. time.
    //if(ctime > finalTime) return nullptr;
    //if(tcol1 > finalTime && tcol2 > finalTime) return nullptr;

     // Avoid collision that will happen after decay.
    if(tcol1 > i1->lifetime() || tcol2 > i2->lifetime()) return nullptr;


    // Can const. quark interact within a formation time?
    double qfac1 = i1->getTf() > tcol1 ? i1->qFactor() : 1.0;
    double qfac2 = i2->getTf() > tcol2 ? i2->qFactor() : 1.0;

    if(constQCollisionOnly) {
      if(qfac1 == 1.0 && qfac2 == 1.0) {
      if(icltyp == 2) return nullptr; // exclude MM collision
      if(icltyp == 3) return nullptr; // exclude MM collision
      if(icltyp == 4) return nullptr; // exclude BBar collision
      if(icltyp == 5) return nullptr; // exclude BbarBar collision
      }
    }

    //...Get total cross section.
    //double pr=sqrt(prsq);
    double m01 = i1->getMass();
    double m02 = i2->getMass();
    srt -=  (m1-m01) + (m2-m02);
    if(srt<=m01+m02) {
      cout << "Collision::hit srt < m01+m02 srt= "<<srt << " m01= "<<m01
	<< " m02= "<< m02<<endl;
      exit(1);
    }
    double pr=PCM(srt,m01,m02);
    //double srt0=sqrt(m01*m01+pr*pr)+sqrt(m02*m02+pr*pr);

    CollisionPair cpair = CollisionPair(icltyp,i1,i2,srt,pr);
    cpair.qFactor(qfac1,qfac2);

    // Compute total cross section.
    if(qfac1 > 0.0 && qfac2 > 0.0) {
      sig = xsection->sigma(cpair)/overSample;
    } else {
      double sigel=0.0;
      XsecTable::additiveQuarkModel(i1->getID(),i2->getID(),sig,sigel);
      //cpair.setXS(sig,sigel);
      //cpair.setXS(sig,sig);
      cpair.setXS(sigel,sigel);
    }

    // Glauber-Gribov Color fluctuation by Strickman
    if(ggcf && (sig > 30.0 && sig < rMaxCutSq*M_PI*10)) {
      ggcf->setParam(sig,omegaSigma);
      //ggcf->computeParameter2(sig,1.0);
      if(ggcf->computeParameter(sig)) {
        sig = 0.5*(ggcf->sample() + ggcf->sample());
      }
    }

    // XXX test for Y-pi collision.
    //if(icltyp ==1) {
      //double sig0=sig;
      //if(i1->baryon() !=0 && i1->strange() < 0) sig = 0.0;
      //if(i2->baryon() !=0 && i2->strange() < 0) sig = 0.0;
      //if((i1->baryon() !=0 && i1->strange() < 0) ||
      //  (i2->baryon() !=0 && i2->strange() < 0)) {
	//cout << "srt= " << srt << "bsq= "<< bsq << " cut= "<< rMaxCutSq
	//  <<" sig= "<< sig << " i1= "<< i1->getID() << " i2= " << i2->getID()
	//   << " " << i1->getParticleDataEntry()->name()
	//   << " " << i2->getParticleDataEntry()->name()
	//  <<endl;
    //}

    hitPath=14;
    // Is their impact parameter small enough?
    if(bsq*M_PI > 0.1*sig*qfac1*qfac2) return nullptr;

    brel = sqrt(max(0.0,bsq));
    bMax=max(brel,bMax);

    TwoBodyInterList* it = new TwoBodyInterList(cpair,i1,i2, ctime,tcol1,tcol2,brel);

    hitPath=0;
    return it;

}

void Collision::PrintCollision(ostream& os) const
{
    os << "** after scatter pnew size= " << pnew.size() << endl;
    if(pnew.size()==1) os << "<< Absorption >>" << endl;

    list<EventParticle*>::const_iterator i;
    for(i = pnew.begin(); i != pnew.end(); ++i) {
	(*i)->print(os);
    }
}

//***********************************************************************

int Collision::collisionType(EventParticle* i1, EventParticle* i2)
{
//...Purpose: to determine type of collision and maximum distance.

    rMaxCutSq = 5 * 5;
    int kf1 = i1->getID();
    int kf2 = i2->getID();

    // skip if lepton, gamma,etc
    if(abs(kf1)>10 && abs(kf1)<100 && kf1 != 21) return 0;
    if(abs(kf2)>10 && abs(kf2)<100 && kf2 != 21) return 0;

    int ibar1 = i1->baryon();
    int ibar2 = i2->baryon();

    int kfl1 = (abs(kf1)/10) % 10;
    int kfl2 = (abs(kf2)/10) % 10;
    int iq1=0;
    int iq2=0;
    if(kfl1 == 0 || kf1 == 21 || abs(kf1) <= 10) iq1=1;
    if(kfl2 == 0 || kf2 == 21 || abs(kf2) <= 10) iq2=1;

//...Hadron-hadron collisions.
      if(iq1 == 0 && iq2 == 0) {


      if(ibar1*ibar2 ==  9) {               // B-B or antiB-antiB
	rMaxCutSq = rMaxCutSqBB;
	//int kfa1=abs(kf1);
	//int kfa2=abs(kf2);
	//int inucl1=0, inucl2=0;
        //if(kfa1==2112 || kfa1==2212  || kfa1 ==3122)inucl1=1;
        //if(kfa2==2112 || kfa2==2212  || kfa2 ==3122)inucl2=1;
        //if(inucl1*inucl2 == 1) rMaxCutSq = rMaxCutSqNN;
        if(i1->lifetime() >1e+20 && i2->lifetime() >1e+20) rMaxCutSq = rMaxCutSqNN;
        if(kf1 < 0) return 5; // Bbar-Bbar
        return 1;
      } else if(ibar1 == 0 &&  ibar2 == 0) { // M-M
          rMaxCutSq = rMaxCutSqMM;
          return 3;
      } else if(ibar1*ibar2 == 0) {          // M-B
          rMaxCutSq = rMaxCutSqMB;
          return 2;
      } else if(ibar1*ibar2 == -9) {        // AntiB-B
          rMaxCutSq = rMaxCutSqBBar;
          return 4;
      } else {
          cout << " Collision::collisionType invalid kf "
              <<  " kf1= " << kf1
              <<  " iba= " << ibar1
              <<  " kf2= " << kf2
              <<  " iba= " << ibar2
              << endl;
	  return 0;
      }


//...Parton-parton collisions
      } else if(iq1 != 0 && iq2 != 0) {
         return 7;
//...Parton-hadron collisions.
      } else {
         return 6;
      }

}

void Collision::deleteParticle()
{
  int idecp=0;
  list<EventParticle*>::iterator first = plist.begin();
  while(first != plist.end()) {
    list<EventParticle*>::iterator next =  first;
    ++next;
    if((*first)->getStatus()<0) {
      idecp++;
      delete *first;
      plist.erase(first);
    }
    first = next;
  }

}

void Collision::deleteParticle(int istat)
{
  int idecp=0;
  list<EventParticle*>::iterator first = plist.begin();
  while(first != plist.end()) {
    list<EventParticle*>::iterator next =  first;
    ++next;
    if((*first)->getStatus()==istat) {
      idecp++;
      delete *first;
      plist.erase(first);
    }
    first = next;
  }

}

// make collision list into collision predictor arrays
// for boosted two nuclei.
// This is used to estimate initial BB collisions.
void Collision::makeCollisionListTwoNuclei(double itime,double gtime)
{
  std::multiset<InterList*,TimeOrder> interlist;

  initialTime=itime;
  finalTime=gtime;
  for(const auto& i1 : pnewA)
  for(const auto& i2 : pnewB) {
    TwoBodyInterList* it = hit(i1,i2);
    if(it !=nullptr) interlist.insert(it);
  }

  predictedColl = interlist.size();
  numberOfCollision=0;
  if(predictedColl==0) {
    lastCollisionTime = 0.0;
    passingTime=0.0;
    numberOfParticipants =0;
    //return;

  } else {
    lastCollisionTime = (*interlist.begin())->getCollisionOrderTime();
    //passingTime = lastCollisionTime;
    passingTime = (*(--interlist.end()))->getCollisionOrderTime();
    numberOfParticipants = computeNumberOfParticipants(interlist);
  }

  for(auto p=pnewA.begin();p != pnewA.end();p++) {
    plist.push_front(*p);
  }
  for(auto p=pnewB.begin();p != pnewB.end();p++) {
    plist.push_front(*p);
  }

  pnewA.clear();
  pnewB.clear();

  for(auto& it: interlist) delete it;
  interlist.clear();

}

int Collision::computeNumberOfParticipants(std::multiset<InterList*,TimeOrder>& interlist)
{
  for(const auto& it: interlist) {
    EventParticle* i1 = it->getParticle(0);
    i1->setParent(1);
    i1 = it->getParticle(1);
    i1->setParent(1);
  }

  int npart=0;
  for(auto& p : pnewA) {
    if(p->getParent()==1) {
      npart++; p->setParent(0);
    }else if(removeSpectator) {
      p->updateR(finalTime,0);
    }
  }
  for(auto& p : pnewB) {
    if(p->getParent()==1) {
      npart++; p->setParent(0);
    }else if(removeSpectator) {
      p->updateR(finalTime,0);
    }
  }

  return npart;
}

bool Collision::doPauliBlocking(InterList* inter,
	vector<EventParticle*>& outgoing,int opt)
{
  int np = inter->getNumberOfInComing();
  EventParticle* i1 = inter->getParticle(0);
  EventParticle* i2 = (np == 2) ? inter->getParticle(1):i1;

  // Loop over newly produced particles. 
  for(const auto& ip : outgoing) {
    if(ip->getStatus() < 0) continue;
      int idp = ip->getID();
      if(idp != 2212 && idp != 2112) continue;
      double ctime=ip->getT();
      //if(ctime < ip->getT()) continue; // not formed
      if(ctime < ip->getTf()) continue; // not formed
      Vec4 r = ip->getR();
      Vec4 p = ip->getP();
      double phase = 0.0;

    // Loop over all particles. 
    for(const auto& i3 : plist) {
      if(i1 == i3) continue; // exclude incoming particle 1.
      if(i2 == i3) continue; // exclude incoming particle 2.
      if(idp != i3->getID()) continue;  // not a nucleon.
      if(ctime < i3->getT()) continue; // not formed
      if(ctime < i3->getTf()) continue; // not formed
      Vec4 r3 = i3->propagate(ctime,opt);
      Vec4 p3 = i3->getP();
      Vec4 dr = r - r3;
      Vec4 dp = p - p3;
      dr[0]=0.0;

      if(withBox) {
        dr[1] = modulo(dr[1] + xBox/2, xBox) - xBox/2;
        dr[2] = modulo(dr[2] + yBox/2, yBox) - yBox/2;
        dr[3] = modulo(dr[3] + zBox/2, zBox) - zBox/2;
      }
//    dr[1] = modulo(dr[1] + Lbox(1)/2, Lbox(1)) - Lbox(1)/2
//    dr[2] = modulo(dr[2] + Lbox(2)/2, Lbox(2)) - Lbox(2)/2
//    dr[3] = modulo(dr[3] + Lbox(3)/2, Lbox(3)) - Lbox(3)/2

      double s = m2(p,p3);
      Vec4 P = p+p3;
      //Vec4 u  = (p+p3)/m(p,p3);
      //double dr2 = cmDistanceSquare(dr,u);
      //double dp2 = cmDistanceSquare(dp,u);
      double dr2 = dr.m2Calc() - pow2(dr*P)/s;
      double dp2 = dp.m2Calc() - pow2(dp*P)/s;
      phase += exp(pauliR*dr2 + pauliP*dp2);
    }

    //cout << "phase= "<< pauliC*phase <<endl;
    //cin.get();

    if(pauliC*phase/overSample > rndm->flat()) return true; // Pauli-blocked.

  }

  // No Pauli blocked.
  return false;

}

void Collision::propagate(double ctime, int opt, int step)
{
  if(step==1) return;
  for(auto i = plist.begin(); i != plist.end(); ++i) {
    if(ctime > (*i)->getT()) (*i)->updateR(ctime,opt);
  }
}


} // end namespace jam2


