#include <jam2/fluid/Fluid2Particle.h>
#include <jam2/hadrons/GaussPoints.h>

namespace jam2 {

using namespace std;

// pariticle freezeout multiplicity table.
int Fluid2Particle::nT=170;
int Fluid2Particle::nBmu=60;
int Fluid2Particle::nSmu=40;
double Fluid2Particle::tMin=0.003;
double Fluid2Particle::bMuMin=0.0;
double Fluid2Particle::sMuMin=0.0;
double Fluid2Particle::dTmp=0.001;
double Fluid2Particle::dbMu=0.02;
double Fluid2Particle::dsMu=0.02;

Fluid2Particle::Fluid2Particle(Pythia8::Settings* s,JamParticleData* pd,
	FreezeOut* fr, Fluid* fl,Pythia8::Rndm* r)
{
  settings = s;
  particleTable = pd;
  freezeout = fr;
  fluid = fl;
  rndm = r;
  nTest=1;
  nFreezeout=0;

  // Option for sampling. =1: only total momentum 
  // =2: energy and momentum conservation are recovered
  optEcon=1;

  mstc144= settings->mode("Hydro:OptFreezeOut");
  isDebug=settings->mode("Check:Debug");

  // Some constants to be used.
  // 4pi/(2pi*hbar)^3
  facStat=1.0/(2*M_PI*M_PI*HBARC*HBARC*HBARC);
  volF    = fluid->getVolH();
  volume  = fluid->getVol();
  tFreezeCut = freezeout->getTCut();

  makeHadronTable();

  useMTable=0;
  // compute total multiplicity at freezeout.
  if(useMTable) makeMultiplicityTable();
}

Fluid2Particle::~Fluid2Particle()
{
  cout << "total freeze-out surface "<< nFreezeout<<endl;
}

void Fluid2Particle::makeHadronTable()
{
  smHadrons.clear();

  ParticleTable* hadrons = particleTable->getSMHadrons();

  // Loop over particles.
  int n=hadrons->size();
  for(int i=0;i<n;i++) {
    ParticleDataEntry* p=hadrons->getParticle(i);
    int kf = p->id();
    if(kf == 130 || kf == 310) continue;    // K_L0 or K_S0

    int charge=p->chargeType(kf)/3;
    int ispin=kf % 10;      // spin 2J+1
    int ibary = p->baryonNumberType(kf)/3;
    int strangeType=p->nQuarksInCode(3);
    if(ibary==0) {
      strangeType = kf > 0 ? strangeType : -strangeType;
      // eta' phi f0(1710)...
      if((abs(kf)/100)%10 ==3 && (abs(kf)/10)%10==3) strangeType=0;
     } else {
       strangeType = kf < 0 ? strangeType : -strangeType;
     }
     int istat=1;
     if(ibary == 0) istat=-1;
     //int cmu = strangeType*smu + ibary*bmu;
     double m=p->m0();                 // particle pole mass
     double w=p->mWidth();

     smHadrons.push_back(
	     SMHadron(kf,ibary,charge,strangeType,ispin,istat,m,w,p));

    // Anti-particle if it exists.
    if(p->hasAnti()) {
      smHadrons.push_back(
	      SMHadron(-kf,-ibary,-charge,-strangeType,ispin,istat,m,w,p));
    }

  }
}

// The information of this table is already contains the EoS table, which is saved in eosTabD[]
// if you use other EoS table, you may need to construct the table.
void Fluid2Particle::makeMultiplicityTable()
{
  cout << "Fluid2Particle::making multiplicity table "<<endl;
  for(int k=0;k<nSmu;k++) {
    double s=sMuMin+k*dsMu;  
    for(int j=0;j<nBmu;j++) {
      double b=bMuMin+j*dbMu;  
      for(int i=0;i<nT;i++) {
        double t=tMin+i*dTmp;
	//if((b+s)/t <40.0) {
        pFreezeoutTable.push_back(totalMultiplicity(t,b,s));
	cout << " T= "<< t << " B= "<< b << " s= "<< s<< " dens= "<< pFreezeoutTable.back()<< endl;
      }
    }
  }

}


//...Poisson distribution for n>=0
int Fluid2Particle::generate_poisson(double mean)
{
  //std::poisson_distribution<int> distribution(mean);
  //return distribution(generator);

  if(mean < 100) {
    double expmean=exp(-mean);   // mean of Poisson
    double pir = 1;
    int N = -1;
    do {
      N++;
      pir *= rndm->flat();
    } while(pir > expmean);
    return N;

  }

  return int(mean+0.5)+sqrt(mean)*rndm->gauss();

}

//--------------------------------------------------------------------------
// Convert fluid element to jam particle.   hydro2jamall
void Fluid2Particle::convert(double gtime,int icheck)
{
//...gtime: current time for fluid simulation.
//...icheck=0: final call, i.e. all fluid elements are below
//...particlization energy density e_p.
//...icheck>0: some of fluid elements are below e_p.

  outgoing.clear();

  // Compute total energy-momentum that will be converted into particles.
  fConv=0.0;
  fConvB=0.0;
  int nfreezeout = freezeout->getN();
  nFreezeout += nfreezeout;
  Vec4 ptot=0.0;
  int ntot=0;
  double dtot=0.0;
  int nbar=0;
  int nch=0;
  int nstr=0;

  if(isDebug) cout << "Fluid2Partile::convert time= "<< gtime << " nfreezeout= "<< nfreezeout << endl;

  if(nfreezeout <= 0) {

    return;

    /*
    if(icheck != 0) return; // not final.
    freezeout->clear();
    nfreezeout = freezeout->isochronousFreezeout(fluid->getElement(),gtime,1);
    if(nfreezeout == 0) {
        cout << "(Fluid2Particle::) still nfreezeout=0?"<<endl;
        return;
    }
    ntot=sampling1p(icheck,nbar,nch,nstr,ptot,dtot);
    */

  // sample particles.
  } else {
    int ntry=0;
    do {
      ntot = sampling4c(icheck,nbar,nch,nstr,ptot,dtot);
      if(ntot==0 && dtot<1.0e-3) break;
      if(ntry++ == 20) break;
    } while(ntot ==0 );


    /*
    if(isDebug) {
      cout << "1 Fluid2Particle::convert time= "<< gtime << " ntot= "<< ntot
        << " ntry= "<< ntry
	<< " ntot= "<<  dtot
	<< " nbar= "<<  nbar
        << " nfreezeout= "<< nfreezeout <<endl;
      for(auto& p: outgoing) {
        cout << p->getT() << " id= " << p->getID() << " p= " << p->getX() << " " << p->getY() << " " << p->getZ()<<endl;
      }
  }
  */

    // in this time-step, there is no particle in dynamical freezeout.
    if(ntot==0 && icheck !=0) {
      //freezeout->clear();
      return;
    }

  }

  if(ntot == 0) {
    ntot=sampling1p(icheck,nbar,nch,nstr,ptot,dtot);
    if(ntot <= 0) {
      if(ntot == 0 && icheck != 0) return;
      freezeout->clear();
      freezeout->isochronousFreezeout(fluid->getElement(),gtime,1);
      if(nfreezeout == 0) {
          cout << "Fluid2Particle::still nfreezeout=0 at the end of hydro?"<<endl;
          return;
      }
      ntot=sampling1p(icheck,nbar,nch,nstr,ptot,dtot);
      if(ntot == 0) {
           cout << "no particle after hydro? ntot= "<< ntot <<endl;
           return;
      }
    }
  }


  Vec4 ptotal = fConv;
  nbar = fConvB;
  if(icheck == 0) {
    ptotal = fluid->pFluidTotal();
    fluid->resetNC();
    fluid->reset();

  // Real time particlization from isochronos hyper surface.
  //} else if(mstc144 == 11) {
  } else {
    fConv=0.0;
    for(int i=0;i<nfreezeout;i++) {
      int ix=freezeout->r(1,i);
      int iy=freezeout->r(2,i);
      int iz=freezeout->r(3,i);
      fConv += fluid->u(ix,iy,iz)*volF;
      fluid->reset_fluid_element(ix,iy,iz);
    }
    ptotal = fConv;
  //} else {
  //  ptotal = fConv;
    nbar = fConvB;
  }


  //if(icheck==0 && optEcon >= 1) recover_mom(ptot, ptotal);
  if(optEcon >= 1) recover_mom(ptot, ptotal);

  // Reset fluid information: Dynamical freeze-out.
  fluid->removeTotalU(fConv);
  fluid->removeTotalBaryon((double)nbar);
  fluid->removeTotalCharge((double)nch);
  fluid->removeTotalStrangeness((double)nstr);

  // Reset freeze-out vectors.

  if(mstc144>=1) {
    if(icheck==0) freezeout->reset(fluid->getElement());
    else freezeout->clear();
  }

//...Output energy-momentum, baryon number, charge, strangeness.
  if(isDebug) {
    cout << "Fluid2Particle::convert time= "<< gtime << " ntot= "<< ntot << "  ptotal= "<< ptotal
      << " ptot_p= "<< ptot;
    for(auto& p: outgoing) {
      cout << "fluid2p t= "<< p->getT() << " id= " << p->getID() << " p= " << p->getX() << " " << p->getY() << " " << p->getZ()<<endl;
    }
    if(abs(ptot.pAbs()-ptotal.pAbs())>1e-8) {
      cout << "Fluid2Particle momentum not conserverd "<<endl;
      exit(1);
    }
  }

}

//------------------------------------------------------------------------
//  MC sampling of particles from freeze-out surface by the method of
//  Scott Pratt III.
int Fluid2Particle::sampling4c(int icheck,int& nbar,int& nch,int& nstr,
	Vec4& ptot, double& nftot)
{
  nbar=0;
  nch=0;
  nstr=0;
  ptot=0.0;
  nftot=0.0;
  int ntot=0;

  // Loop over freeze-out hyper-surface.
  int nfreezeout = freezeout->getN();
  for(int n=0;n<nfreezeout;n++) {

    if(freezeout->temperature(n) < tFreezeCut) continue;

    double gam=freezeout->v(0,n);
    double v1=freezeout->v(1,n);
    double v2=freezeout->v(2,n);
    double v3=freezeout->v(3,n);
    double dst=freezeout->dsigma(0,n);
    double dsx=freezeout->dsigma(1,n);
    double dsy=freezeout->dsigma(2,n);
    double dsz=freezeout->dsigma(3,n);
    double vv=v1*v1 + v2*v2 + v3*v3;
    if(vv > 1.0) {
      cout << n << " sampling4c v>1 ? " << v1
   	<< " v2= " << v2
	<< " v3= " << v3
	<< " vv= " << vv
	<< " T= " << freezeout->temperature(n)*HBARC
	<<endl;
      continue;
    }
    //double gam=1.0/sqrt(1.0-v1*v1-v2*v2-v3*v3);
    double vol=gam*(dst + dsx*v1 + dsy*v2 + dsz*v3);
    double vopt=vol+sqrt(vol*vol - (dst*dst - dsx*dsx-dsy*dsy-dsz*dsz));
    double tf=freezeout->temperature(n)*HBARC;
    double bmu=freezeout->mub(n)*HBARC;
    double smu=freezeout->mus(n)*HBARC;

    // Determine number of particle to be generated.
    double pdns=freezeout->num(n);
    nftot += pdns*vol*nTest;

    int npar=generate_poisson(pdns*vopt)*nTest;
    if(npar == 0) continue;

    // Calculate partial multiplicities.
    //pdns=totalMultiplicity(tf,bmu,smu);
 

    for(int i=0;i<npar;i++) {

      //SMHadron* hc = selectParticle(pdns);

      SMHadron* hc = selectParticle(pdns,tf,bmu,smu);
    
      int ibary = hc->baryonType();
      int str= hc->strangeType();
      int cha= hc->chargeType();
      int istat=hc->stat();
      double cmu = str*smu + ibary*bmu;
      //double mh=hc->m0(); // pole mass
      double mh=hc->particle()->mSel();

      double pp=thermaldist3(mh,tf,cmu,istat);
      double cost=2.0*rndm->flat()-1.0;
      double sint=sqrt(1.0-cost*cost);
      double phi=2*M_PI*rndm->flat();
      //pc[1]=pp*sint*cos(phi); pc[2]=pp*sint*sin(phi); pc[3]=pp*cost; pc[0]=sqrt(pp*pp+mh*mh);
      Vec4 pc(pp*sint*cos(phi), pp*sint*sin(phi), pp*cost,sqrt(pp*pp+mh*mh));
      double er=pc[0];
      pc.bst(v1,v2,v3,gam);
      double pds = pc[0]*dst+pc[1]*dsx+pc[2]*dsy+pc[3]*dsz;
      if(pds <= 0.0) continue;
      double prob=pds/(er*vopt);
      if(prob < 0.0 || prob > 1.0) {
        cout << "prob? " << prob << " vopt= "<< vopt << endl;
	exit(1);
      }
      if(rndm->flat() > prob) continue;

      //...Set coordinate of particle.
      double r[4];
      r[1]=fluid->getX( freezeout->r(1,n), rndm->flat() );
      r[2]=fluid->getY( freezeout->r(2,n), rndm->flat() );
      r[3]=fluid->getZ( freezeout->r(3,n), rndm->flat() );

      //...Set time and formation time of particle.
      r[0]=freezeout->t(n);
      Vec4 rc(r[1],r[2],r[3],r[0]);

      /*
      // Real-time freeze-out by the Cooper-Frye formula:remove
      // energy-momentum of particle from fluid element.
      if(mstc144 == 12 && icheck != 0) {
        int ix=freezeout->r(1,n);
        int iy=freezeout->r(2,n);
        int iz=freezeout->r(3,n);
	Vec4 utmp = fluid->u(ix,iy,iz) - pc/volF;;
        double det=utmp.m2Calc();
        //if(utmp[0] <= 0.0 || det <= 0.0) {
        if(det <= 0.0) {
	  cout << " det= "<< det << " utmp[0]= "<< utmp[0]<<endl;
            continue; // Skip this particle because of too much energy.

//...Nevertheless let's include this particle setting fluid element to be zero
//            fconv(:) = fconv(:) + u(:,ix,iy,iz)*volf
//            call reset_fluid_element(frr(1,n),frr(2,n),frr(3,n))
//            frtmp(n)=0.0;

//...No problem. Subtract particle energy-momentum from the fluid element. 
        } else {
          //fConv += fluid->u(ix,iy,iz) - utmp;
          fConv += pc/volF;
          fluid->setU(ix,iy,iz,utmp);
          fluid->setB(ix,iy,iz,ibary/volume);
          fluid->setLocalQuantity(ix,iy,iz);
	  //fluid->removeTotalU(pc/volF);
	  //fluid->removeTotalBaryon(ibary);
	  //fluid->removeTotalCharge(cha);
	  //fluid->removeTotalStrangeness(str);
	}
        fConvB += ibary;
      }
      */

      /*
      if(mstc144 == 12 && icheck != 0) {
        int ix=freezeout->r(1,n);
        int iy=freezeout->r(2,n);
        int iz=freezeout->r(3,n);
        fConv += fluid->u(ix,iy,iz)*volF;
        fConvB += ibary;
        fluid->reset_fluid_element(ix,iy,iz);
      }
      */

      p2jam(hc,rc,pc,mh,v1,v2,v3);

      ntot++;
      ptot += pc/nTest;
      nbar += ibary/nTest;
      nstr += str/nTest;
      nch  += cha/nTest;
    } // end loop over npar

  } // end loop freezeout surface.

  return ntot;

}


int Fluid2Particle::sampling1p(int icheck,int& nbar,int& nch,int& nstr,
	Vec4& ptot, double& nftot)
{
  nbar=0;
  nch=0;
  nstr=0;
  ptot=0.0;
  int ntot=0;
  nftot=0.0;

  int nfreezeout = freezeout->getN();

  cout << "sampling1p nfreezeout= " << nfreezeout<<endl;

  // Loop over freeze-out hyper-surface.
  for(int n=0;n<nfreezeout;n++) {

    if(freezeout->temperature(n) < tFreezeCut) continue;

    double v1=freezeout->v(1,n);
    double v2=freezeout->v(2,n);
    double v3=freezeout->v(3,n);
    double dst=freezeout->dsigma(0,n);
    double dsx=freezeout->dsigma(1,n);
    double dsy=freezeout->dsigma(2,n);
    double dsz=freezeout->dsigma(3,n);
    double gam=1.0/sqrt(1.0-v1*v1-v2*v2-v3*v3);
    double vol=gam*(dst + dsx*v1 + dsy*v2 + dsz*v3);
    //double vopt=vol+sqrt(vol*vol - (dst*dst - dsx*dsx-dsy*dsy-dsz*dsz));
    double tf=freezeout->temperature(n)*HBARC;
    double bmu=freezeout->mub(n)*HBARC;
    double smu=freezeout->mus(n)*HBARC;

    nftot += freezeout->num(n)*vol*nTest;
    // Determine number of particle to be generated.
    double pdns=freezeout->num(n);
    if(pdns <= 0.0) continue;

    int npar=1;

    // Calculate partial multiplicities.
    //pdns=totalMultiplicity(tf,bmu,smu);

    for(int i=0;i<npar;i++) {

      double pds=0.0;
      double mh;
      Vec4 pc;
      int itry=0;
      SMHadron* hc;
      int ibary, str;
      do {
        //hc = selectParticle(pdns);
        hc = selectParticle(pdns,tf,bmu,smu);

        ibary = hc->baryonType();
        str= hc->strangeType();
        int istat=hc->stat();
        double cmu = str*smu + ibary*bmu;
        mh=hc->particle()->mSel();
        double pp=thermaldist3(mh,tf,cmu,istat);
        double cost=2.0*rndm->flat()-1.0;
        double sint=sqrt(1.0-cost*cost);
        double phi=2*M_PI*rndm->flat();
        pc[1]=pp*sint*cos(phi);
        pc[2]=pp*sint*sin(phi);
        pc[3]=pp*cost;
        pc[0]=sqrt(pp*pp + mh*mh);
        //double er=pc[0];
        pc.bst(v1,v2,v3,gam);
        pds = pc[0]*dst+pc[1]*dsx+pc[2]*dsy+pc[3]*dsz;
	if(itry <= 20) {
	    cout << itry << " samplingp1 infinite loop? T= " << tf
		<< " mu= " << bmu
		<< " m= " << mh
		<< " id= " << hc->id()
		<< " pds= " << pds
		<< " vol= " << vol
		<< endl;
	}
      } while (pds <= 0.0 && ++itry < 20);

//      cout << " pds= " << pds << " n= "<<n<<endl;

      if(pds <=0.0) continue;

      //double prob=pds/(er*vopt);
      //if(rndm->flat() > prob) goto L100;

      //...Set coordinate of particle.
      double r[4];
      r[1]=fluid->getX( freezeout->r(1,n), rndm->flat() );
      r[2]=fluid->getY( freezeout->r(2,n), rndm->flat() );
      r[3]=fluid->getZ( freezeout->r(3,n), rndm->flat() );

      //...Set time and formation time of particle.
      r[0]=freezeout->t(n);
      Vec4 rc(r[1],r[2],r[3],r[0]);

      // Real-time freeze-out by the Cooper-Frye formula:remove
      // energy-momentum of particle from fluid element.
      if(mstc144 == 12 && icheck != 0) {
        int ix=freezeout->r(1,n);
        int iy=freezeout->r(2,n);
        int iz=freezeout->r(3,n);

	Vec4 utmp = fluid->u(ix,iy,iz) - pc/volF;;
        double det=utmp.m2Calc();
        if(utmp[0] <= 0.0 || det <= 0.0) {
            fConv += fluid->u(ix,iy,iz)*volF;
            fluid->reset_fluid_element(ix,iy,iz);
            freezeout->temperature(n,0.0);

       // No problem. Subtract particle energy-momentum from the fluid element. 
        } else {
          //fConv += fluid->u(ix,iy,iz) - utmp;
          fConv += pc/volF;
          fluid->setU(ix,iy,iz,utmp);
          fluid->setB(ix,iy,iz,ibary/volume);
          fluid->setLocalQuantity(ix,iy,iz);
	  //fluid->removeTotalU(pc/volF);
	  //fluid->removeTotalBaryon(ibary);
	  //fluid->removeTotalCharge(cha);
	  //fluid->removeTotalStrangeness(str);
	}
        //fconv(4) = fconv(4) + ibary
        fConvB += ibary;
      }

      p2jam(hc,rc,pc,mh,v1,v2,v3);

      int cha= hc->chargeType();
      ntot++;
      ptot += pc/nTest;
      nbar += ibary/nTest;
      nstr += str/nTest;
      nch  += cha/nTest;

      // opne particle is produced, it is enough.
      return ntot;

    }

  } // end loop freezeout surface.

  if(ntot==0) cout << " no particles are created in sampling1p"<<endl;
  return ntot;

}

// Put particle into JAM.
void Fluid2Particle::p2jam(SMHadron* hc,Vec4& rc,Vec4& pc,double mh,
	double v1,double v2, double v3)
{
  //ParticleDataEntry* pd = jamParticleData->find(id1);
  ParticleDataEntry* pd = hc->particle();
  int id=hc->id();
  int pid=particleTable->pid(id);
  EventParticle* pa = new EventParticle(id,pd);
  pa->setPID(pid);
  pa->setMass(mh);
  pa->setCoordinate(rc);
  pa->setVertex(rc);
  RotBstMatrix  MfromCM;
  MfromCM.bst(v1,v2,v3);
  //p3.rotbst(MfromCM);
  pa->setBoostMatrix(MfromCM);
  pa->setMomentum(pc);
  pa->setOnShell();
  pa->setNumberOfColl(1000);
  pa->setParent(1000);
  //pa->lastColl(-1);
  pa->setConstQuark();
  double m=pc.mCalc();
  double e=pc.e();
  double t = rc[0];
  double dect = particleTable->lifeTime(pd,m,e);
  pa->setLifeTime(t+dect);
  outgoing.push_back(pa);

 if(m<0.0 || abs(m-mh) > 1e-8) {
  cout << "p2jam m<0? " << pc.mCalc() << " id= "<< hc->id()<<endl;
  cout << "pc = " << pc  <<endl;
  cout << "mh = " << mh  <<endl;
  exit(1);
  }

}

//--------------------------------------------------------------------------
SMHadron* Fluid2Particle::selectParticle(double& tdns,double tch,double bmu,double smu)
{
  int ntry=0;
  double xpart=tdns*rndm->flat();
  double tdns2=0.0;
  do {
    for(int i=0;i<(int)smHadrons.size();i++) {
      double dens;
      if(ntry==0) {
        double m=smHadrons[i].m0();
        int ibary=smHadrons[i].baryonType();
        int str=smHadrons[i].strangeType();
        int ispin=smHadrons[i].spinType();
        int istat=smHadrons[i].stat();
        double cmu = str*smu + ibary*bmu;
        dens=pdensity(m,tch,cmu,istat)*ispin*facStat;
        smHadrons[i].setMult(dens);
      } else {
        dens = smHadrons[i].multiplicity();
      }
      xpart -= dens;
      if(xpart <= 0.0) return &smHadrons[i];
      tdns2 += dens;
    }
    xpart=tdns2*rndm->flat();
    //cout << ntry << " tdns= "<< tdns << " tdns2= "<< tdns2<<endl;
    tdns=tdns2;
  } while(ntry++<3);

  cout << " selectParticle infinite loop? " << tdns << endl;
  exit(1);
}

//--------------------------------------------------------------------------
SMHadron* Fluid2Particle::selectParticle(double& tdns)
{
  int ntry=0;
  double xpart=tdns*rndm->flat();
  double tdns2=0.0;
  do {
    for(int i=0;i<(int)smHadrons.size();i++) {
      xpart -= smHadrons[i].multiplicity();
      tdns2 += smHadrons[i].multiplicity();
      if(xpart <= 0.0) return &smHadrons[i];
    }
    xpart=tdns2*rndm->flat();
    tdns=tdns2;
  } while(ntry++<3);

  cout << " selectParticle infinite loop? " <<endl;
  exit(1);
}

//--------------------------------------------------------------------------
//...Initialize particle multiplicities.
double Fluid2Particle::totalMultiplicity(double tch,double bmu,double smu)
{
  double pmult=0.0;    // hadron density.

  // Loop over all hadrons.
  for(int i=0;i<(int)smHadrons.size();i++) {
    double m=smHadrons[i].m0();
    int ibary=smHadrons[i].baryonType();
    int str=smHadrons[i].strangeType();
    int ispin=smHadrons[i].spinType();
    int istat=smHadrons[i].stat();
    double cmu = str*smu + ibary*bmu;
    double dens=pdensity(m,tch,cmu,istat)*ispin*facStat;
    smHadrons[i].setMult(dens);
    pmult += dens;
  }

  if(pmult <= 0.0) {
    cout<< "pmult <0 ? " << pmult
        << " T= " << tch
        << " muB= " << bmu
        << " muS= " << smu
        << endl;
    exit(1);
  }

  return pmult;

}

//--------------------------------------------------------------------------
double Fluid2Particle::pdensity(double pm,double tf,double mu,int istat)
{
  // istat = -1: Bosons.
  // istat = 1:  Fermions.

  // Boltzmann distribution
  //if(istat==0) return pm*pm*tf*std::cyl_bessel_k(2.0, pm/tf)*exp(mu/tf);

  double pdensity=0.0;

  /*
    // use modified Bessel function of the second kind.  -std=c++17
   double stat=-istat;
   double pden=0.0;
   for(int j=1;j<=10;j++) {
     pden += (pow(stat,(j+1))/j)*std::cyl_bessel_k(2.0, j*pm/tf)*exp(j*mu/tf);
   }
   double pdensity2 = pm*pm*tf*pden;
   */

  // Numerical integration.

  /*
  for(int i=0;i<NGP38;i++) {
    double p=xg38[i];
    double e=sqrt(p*p+pm*pm);
    double dis=exp(-(e-mu)/tf);
    pdensity += p*p*wg38[i]*dis/(1.0+istat*dis);
  }
  */

  for(int i=0;i<NGP2;i++) {
    double p=xg12[i];
    double e=sqrt(p*p+pm*pm);
    double dis=exp(-(e-mu)/tf);
    pdensity += p*p*wg12[i]*dis/(1.0+istat*dis);
  }

  /*
  if(abs(pdensity-pdensity1)>1e-8)  {
  cout << "pm= "<< pm << " tf= "<<tf << " mu= "<< mu << " istat= "<< istat
    << " pden= "<< pdensity << " pden1= "<< pdensity << " eps= "<< pdensity-pdensity1
    <<endl;
  cin.get();
  }
  */

  bool err=false;
  //if(istat==-1 && abs(pdensity-pdensity2)>1e-3)  err=true;

  if(pdensity < 0.0 || err) {
    cout<< setprecision(16) << "pdensity <0 ? " << pdensity
        << " T= " << tf
        << " mu= " << mu
        << " m= " << pm
        << " istat " << istat
        << endl;

    /*
   double stat=-istat;
   double pden=0.0;
   for(int j=1;j<=10;j++) {
     pden += (pow(stat,(j+1))/j)*std::cyl_bessel_k(2.0, j*pm/tf)*exp(j*mu/tf);
   }
   pden = pm*pm*tf*pden;
   cout << " pden2= "<< pden <<endl;
   */


  for(int i=0;i<NGP38;i++) {
    double p=xg38[i];
    double e=sqrt(p*p+pm*pm);
    double dis=exp(-(e-mu)/tf);
    cout << " i= " << i<< " dis= "<< p*p*wg38[i]*dis/(1.0+istat*dis)
      << " dis= "<< dis
      <<endl;
  }
  }

  return pdensity;
}

//--------------------------------------------------------------------------
 // Generate the momentum from the thermal distribution.
double Fluid2Particle::thermaldist3(double pm,double tf,double mu,int a)
{
  if(pm == 0.0) {
    cout << "thermaldist m=" << pm <<endl;
    exit(1);
  }

  double m=pm/tf;
  double xm=mu/tf;
  double  w0=1.0;
  if(a == -1) w0 = 1.0 - exp(-m+xm);
  double pmaxsq = 2 + 2*sqrt(1.0+m*m);      // Peak mass squared
  double emax=sqrt(pmaxsq+m*m);
  double pmax=sqrt(pmaxsq);

  // Log-Logistic with alpha=2: g(x)=2*(x/b)/(b*(1+(x/b)**2)**2)
  double beta=pmax*sqrt(3.);          // Peak is chosen to be the same as MB
  double gmax=9./(beta*8*sqrt(3.0));
  double gfac=1.0/gmax*2/(beta*beta);    // Normalization factor to one

  double p, e;
  do {
    double g, ex;
    do {
      p=rndm->flat();
      p=beta*pow(p/(1.0-p),0.5);
      g=gfac*p/pow2(1.0+pow2(p/beta));
      e=sqrt(p*p+m*m);
      ex=2*log(p/pmax)-e+emax;
    } while(g*rndm->flat() > exp(max(-50.0,min(50.0,ex))));

  } while (rndm->flat() > w0/( 1.0 + a*exp(-e + xm)) );

  return p*tf;

}

void Fluid2Particle::recover_mom(Vec4& ptot,Vec4& pconv)
{
// ntot: total number of particles produced from fluid.
// ptot: total momentum of particles produced from fluid.
// pconv: total momentum of fluid elements converted into particles.

  // Shift momenta to match the required total momentum.
  Vec4 pf=0.0;
  int ntot=outgoing.size();
  Vec4 psub = -(ptot - pconv)/ntot;
  for(int i=0;i<ntot;i++) {
    outgoing[i]->addP(psub);
    outgoing[i]->setOnShell();
    pf += outgoing[i]->getP();
  }

  // Momentum conservation is recovered now.
  ptot=pf;
  if(optEcon == 1) return;

  // Go to the C.M. frame of all particles.
  for(int i=0;i<ntot;i++) {
    outgoing[i]->bstback(pconv);
  }

  // Find a factor to correct total energy.
  //double e=sqrt(pe_tot**2-px_tot**2-py_tot**2-pz_tot**2)
  double pe = pconv.mCalc();
  //double pe = pconv[0];
  double a=1.0;
  int ntry=0;
  double f;
  do {
    f=-pe;
    double  df=0.0;
    for(int i=0;i<ntot;i++) {
      Vec4 p = outgoing[i]->getP();
      double m = outgoing[i]->getMass();
      double e = sqrt(m*m + a*a*p.pAbs2());
      f +=  e;
      df += a/e;
    }
    a -= f / df;

    if(ntry++>50) {
      cout << "recover_mom iter " << ntry
        << " f= "<<  f
        << "f/df= " << a
        << endl;
      exit(1);
    }
  } while(abs(f) > 1e-8);


//...Go back to the computational frame after scaling the momenta.
  pf=0.0;
  for(int i=0;i<ntot;i++) {
    Vec4 p = outgoing[i]->getP();
    double m = outgoing[i]->getMass();
    p *= a;
    p[0] = sqrt(m*m + p.pAbs2());
    p.bst(pconv);
    pf += p;
    outgoing[i]->setP(p);
  }

  // This is not really needed.
  psub = -(pf - pconv)/ntot;
  pf=0.0;
  for(int i=0;i<ntot;i++) {
    outgoing[i]->addP(psub);
    outgoing[i]->setOnShell();
    pf += outgoing[i]->getP();
  }

}

double Fluid2Particle::BWmass(double m0, double wid, double mmin, double mmax)
{
//...Breit Wigner distribution.
  double cons=2.0*(mmin-m0)/wid;
  double const1=atan(cons);
  double const2=M_PI/2.0-const1;
  double xmax=(atan(2.0*(mmax-m0)/wid)-const1)/const2;
  xmax=min(1.0,xmax);
  double x=xmax*rndm->flat();
  double t=tan(x*const2);
  return m0+0.5*wid*(cons+t)/(1.0-cons*t);

  //pjmass= m0 + 0.5*wid*tan( (2.0*pjr(0)-1.0)*atan(2.0*pmas(kc,3)/wid) );

}

} // namespace jam2

