#include <jam2/initcond/BoostedTwoNuclei.h>
#include <jam2/initcond/Nucleus.h>
#include <jam2/collision/Collision.h>
#include <cstdlib>

using namespace std;
using namespace Pythia8;

namespace jam2 {

BoostedTwoNuclei::BoostedTwoNuclei(Settings* s, JamParticleData* p, Rndm* r)
    : InitialCondition(s,p,r)
{
    eCM=settings->parm("Beams:eCM");
    pLab=settings->parm("Beams:pLab");
    eLab=settings->parm("Beams:eLab");
    bMin=settings->parm("Beams:bmin");
    bMax=settings->parm("Beams:bmax");
    zSep=settings->parm("Beams:zseparation");
    compFrame=settings->word("Beams:compFrame");
    optFermiMomentum=settings->mode("Cascade:optFermiMomentum");
    gWidth=1.0;

    paP = particleData->particleDataEntryPtr(2212);
    paN = particleData->particleDataEntryPtr(2112);
    mP=paP->m0();
    mN=paN->m0();

    histWS = settings->flag("HeavyIon:histWS");
    histr = 0; histp=0;
    if(histWS) initHist();

}

void BoostedTwoNuclei::init()
{
    string nuclA = settings->word("Beams:beamA");
    string nuclB = settings->word("Beams:beamB");
    int idProj=findAZ(nuclA);
    int idTarg=findAZ(nuclB);

    //proj = new GLISSANDOModel();
    proj = new Nucleus();
    proj->initPtr(idProj, *settings,*particleData,*rndm);
    proj->init();
    dynamic_cast<Nucleus*>(proj)->setParam();
    //targ = new GLISSANDOModel();
    targ = new Nucleus();
    targ->initPtr(idTarg, *settings,*particleData,*rndm);
    targ->init();
    dynamic_cast<Nucleus*>(targ)->setParam();

    mA=paP->m0();
    mB=paP->m0();

    /*
    if(eLab > 0.0) {
      pLab=sqrt(eLab*(2*mA+eLab));
      eCM = sqrt(pow2(eLab+mA+mB)-pLab*pLab);
      settings->parm("Beams:eCM",eCM);
    } else if(pLab > 0.0) {
      eLab=sqrt(mA*mA+pLab*pLab)-mA;
      eCM = sqrt(pow2(eLab+mA+mB)-pLab*pLab);
      settings->parm("Beams:eCM",eCM);
    } else {
      eLab=(eCM*eCM -mA*mA - mB*mB)/(2*mB)-mA;
      pLab=sqrt(eLab*(2*mA+eLab));
    }
    */

    pzAcm  = 0.5 * sqrtpos( (eCM + mA + mB) * (eCM - mA - mB)
           * (eCM - mA + mB) * (eCM + mA - mB) ) / eCM;
    eA = sqrt(mA*mA + pzAcm*pzAcm);
    eB = sqrt(mB*mB + pzAcm*pzAcm);
    yCM = 0.5*log((eA + pzAcm)/(eA - pzAcm));
    gamA= eA/mA;
    gamB= eB/mB;
    masA = atoi(nuclA.c_str());
    masB = atoi(nuclB.c_str());

    //cout << "nuclA=" << nuclA << " nuclB= " << nuclB  <<endl;
    //cout << "massA= " << masA << " masB=  " << masB  <<endl;
    //cout << "idA=" << idProj << " idTarg= " << idTarg  <<endl;
    //cout << " ecm= "<<eCM << " pLab= "<< pLab << " eLab= "<< eLab <<endl;
}

void BoostedTwoNuclei::generate(Collision* event,int mode)
{
  // Generate impact parameter.
  impactPar = setImpactPar();
  double tStart = settings->parm("Cascade:timeStart");

  for(int iev=0;iev<overSample;iev++) {

    // Generate the position of nucleons inside nucleus.
    projNucl = proj->generate();
    targNucl = targ->generate();

    // Generate Fermi momentum.
    pProj=generateFermiMomentum(projNucl);
    pTarg=generateFermiMomentum(targNucl);
    if(histWS) {
      fill(projNucl,pProj);
      fill(targNucl,pTarg);
    }

    paA = particleData->particleDataEntryPtr(projNucl[0].id());
    paB = particleData->particleDataEntryPtr(targNucl[0].id());
    mA=paA->m0();
    mB=paB->m0();
    int nA = projNucl.size();
    int nB = targNucl.size();
    if(projNucl.size()>1) mA=(mP+mN)/2.0;
    if(targNucl.size()>1) mB=(mP+mN)/2.0;

    // compFrame=nn
    pzAcm  = 0.5 * sqrtpos( (eCM + mA + mB) * (eCM - mA - mB)
           * (eCM - mA + mB) * (eCM + mA - mB) ) / eCM;
    pzBcm  = -pzAcm;

    if(compFrame=="cm") {
      double m1 = mA*nA;
      double m2 = mB*nB;
      double etotal=eLab*nA + m1 + m2;
      double plab = pLab*nA;
      double ecm = sqrt(etotal*etotal - plab*plab);
      double  pz = 0.5 * sqrtpos( (ecm + m1 + m2) * (ecm - m1 - m2)
           * (ecm - m1 + m2) * (ecm + m1 - m2) ) / ecm;
      pzAcm =   pz / nA;
      pzBcm = - pz / nB;
    } else if (compFrame=="lab") {
      pzAcm = pLab;
      pzBcm = 0.0;
    }

    eA = sqrt(mA*mA + pzAcm*pzAcm);
    eB = sqrt(mB*mB + pzBcm*pzBcm);
    gamA= eA/mA;
    gamB= eB/mB;
    yCM = 0.5*log((eA + pzAcm)/(eA - pzAcm));

    double zmaxA=-100;
    for ( int i = 0, N = projNucl.size(); i < N; ++i )
	zmaxA=max(zmaxA,projNucl[i].nPos().pz());
    double rzA = (proj->R()+1.0);
    if(zmaxA> rzA + gamA*zSep) {
      rzA=zmaxA;
      //cout << "rzA= "<< zmaxA << " r= "<< proj->R()+1.0 <<endl;
    }


    double zminB=100;
    for ( int i = 0, N = targNucl.size(); i < N; ++i )
	zminB=min(zminB,targNucl[i].nPos().pz()/gamB);
    double rzB = (targ->R()+1.0);
    if(zminB < -rzB - gamB*zSep) {
      rzB=zminB;
      //cout << "rzB= "<< rzB << " r= "<< targ->R()+1.0 <<endl;
    }


    double xA,zA,xB,zB;
    if(compFrame != "lab"){
       xA= impactPar/2.0;
       zA=  -abs(rzA)/gamA - zSep;
       xB = -impactPar/2.0;
       zB =  abs(rzB)/gamB + zSep;
    } else {
       xA= impactPar;
       zA=  -abs(rzA)/gamA - zSep -abs(rzB);
       xB = 0.0;
       zB = 0.0;
    }


    //Vec4 bvecA( impactPar/2.0, 0.0, -abs(rzA)/gamA - zSep, 0.0);
    //Vec4 bvecB(-impactPar/2.0, 0.0,  abs(rzB)/gamB + zSep, 0.0);
    Vec4 bvecA( xA, 0.0, zA, 0.0);
    Vec4 bvecB( xB, 0.0, zB, 0.0);


    Vec4 ptot=0.0;
    for ( int i = 0, N = projNucl.size(); i < N; ++i ) {
	//projNucl[i].reset();
	projNucl[i].bShift(bvecA);
	Vec4 r=projNucl[i].nPos();
	r[3] /= gamA;
	r += bvecA;
	r[0]= tStart;

	pProj[i][3] = gamA*pProj[i][3] + pzAcm;
        int id=projNucl[i].id();
	double m=mA;
        Pythia8::ParticleDataEntry* pa = paA;
	if(id==2212) {m=mP; pa=paP;}
	if(id==2112) {m=mN; pa=paN;}
	EventParticle* cp = new EventParticle(id,m,r,pProj[i],pa,1);
	cp->setPID(jamParticleData->pid(abs(id)));
	cp->setOnShell();
	cp->setNumberOfColl(1);
	//cp->setConstQuark();
	//event->setPnewList(cp);
	if(mode==0) {
	  event->setPList(cp);
	} else {
	  event->setPnewAList(cp);
	}
        ptot +=cp->getP();
    }

    for ( int i = 0, N = targNucl.size(); i < N; ++i ) {
	//targNucl[i].reset();
	targNucl[i].bShift(bvecB);
	Vec4 r=targNucl[i].nPos();
	r[3] /= gamB;
	r += bvecB;
	r[0]= tStart;

	pTarg[i][3] = gamB*pTarg[i][3] + pzBcm;
        int id=targNucl[i].id();
	double m=mB;
        Pythia8::ParticleDataEntry* pa = paB;
	if(id==2212) {m=mP; pa=paP;}
	if(id==2112) {m=mN; pa=paN;}
	EventParticle* cp = new EventParticle(id,m,r,pTarg[i],pa,-1);
	cp->setPID(jamParticleData->pid(id));
	cp->setOnShell();
	cp->setNumberOfColl(-1);
	//cp->setConstQuark();
	//event->setPnewList(cp);
	if(mode==0) {
	  event->setPList(cp);
	} else {
	  event->setPnewBList(cp);
	}
        ptot +=cp->getP();
    }

  }

    if(mode !=0 ) {
      event->makeCollisionListTwoNuclei(tStart,1000000.0);
      //nCollG = event->getInterListSize(); // collision number.
      nCollG = event->getNumberOfGlauberCollision(); // collision number.
      nPartG = event->getNumberOfParticipants();
    }

    //event->makeCollisionList();

}

vector<Vec4> BoostedTwoNuclei::generateFermiMomentum(vector<Pythia8::Nucleon>& nucl)
{
  int N=nucl.size();
  vector<Vec4> pv(N,0.0);
  if(N ==1 || optFermiMomentum==0) return pv;

  double fac=pow(4*M_PI*gWidth,1.5);
  std::vector<double> rhon(N,0.0);

  for(int i=0;i<N; i++) {
	for(int j=i+1; j<N; j++) {
	    Vec4 dr=nucl[i].nPos() - nucl[j].nPos();
	    double den=(dr.px()*dr.px() + dr.py()*dr.py()+dr.pz()*dr.pz());
	    den = exp(-den/4/gWidth)/fac;
	    rhon[i] += den;
	    rhon[j] += den;
	}
  }

    Vec4 ptot=0.0;
    for(int i=0;i<N; i++) {
	double pfermi=HBARC*pow(1.5*M_PI*M_PI*rhon[i],1.0/3.0);
	double pf=pfermi*pow(rndm->flat(),1.0/3.0);
	double cth=1.0 - 2.0*rndm->flat();
	double sth=sqrt(1.0-cth*cth);
	double phi=2*M_PI*rndm->flat();
	//Vec4 p;
	pv[i][3]=pf*cth;
	pv[i][1]=pf*sth*cos(phi);
	pv[i][2]=pf*sth*sin(phi);
	ptot += pv[i];
	//pv.push_back(p);
    }

    ptot /= N;
    //Vec4 ptotal = 0.0;
    for(int i=0;i<N; i++) {
	int id=nucl[i].id();
	double m=particleData->m0(id);
	pv[i] -= ptot;
	pv[i][0]=sqrt(m*m + pv[i].pAbs2());
	//ptotal += pv[i];
    }

    // check total momentum if it is zero.
    //cout << N << " generateFermi p= "<< ptotal << endl;
    
    //delete [] rhon;

    return pv;
}

void BoostedTwoNuclei::printHist()
{
  ofstream ofst("ws.hist");
  ofst << "# width= "<< widR/2.0 << endl;
  for( int i=0; i< nhist; i++) {
    ofst << setw(7) << i*dR 
 	<< setw(12) << histr[i]/nEvent/overSample
	<< setw(12) << rhor[i]/nEvent/overSample
        << setw(7) << i*dP
	<< setw(12) << histp[i]/nEvent/overSample
	<< setw(12) << rhop[i]/nEvent/overSample
	<<endl;
  }

  ofst.close();
}

void BoostedTwoNuclei::fill(vector<Pythia8::Nucleon>& nucl, vector<Vec4> pf)
{
  nEvent++;
  int N = nucl.size();
  for ( int i = 0;  i < N; ++i ) {
    Vec4 r = nucl[i].nPos();
    Vec4 p = pf[i]/HBARC;
    double rr = r.pAbs2();
    double pp = p.pAbs2();
    int ir = int(sqrt(rr)/dR+0.5);
    int ip = min(nhist,int(sqrt(pp)/dP+0.5));
    histr[ir] += 1.0/(4*M_PI*rr*dR);
    histp[ip] += 1.0/(4*M_PI*pp*dP);
    for (int l=0;l<nhist;l++) {
      Vec4 x(l*dR, 0.0, 0.0, 0.0);
      Vec4 px(l*dP/HBARC,0.0, 0.0, 0.0);
      rhor[l] += exp(-(r-x).pAbs2()/widR)*facR;
      rhop[l] += exp(-(p-px).pAbs2()/widP)*facP;
    }
  }
}

void BoostedTwoNuclei::initHist()
{
  nEvent = 0;
  rMax=15.0; pMax=2.0;
  dR=rMax / nhist;
  dP=pMax / nhist;
  histr = new double [nhist];
  histp = new double [nhist];
  rhor = new double [nhist];
  rhop = new double [nhist];
  for(int i=0;i<nhist;i++) {
    histr[i]=histp[i]=0.0;
    rhor[i]=rhop[i]=0.0;
  }
  widR = 2 * settings->parm("MeanField:gaussWidth");
  widP = 1.0/widR;
  facR=1.0/pow(M_PI*widR,1.5);
  facP=1.0/pow(M_PI*widP,1.5);
}

}
