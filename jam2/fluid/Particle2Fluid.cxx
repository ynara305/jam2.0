#include <jam2/fluid/Particle2Fluid.h>

namespace jam2 {

//list<EventParticle*>& plist,

Particle2Fluid::Particle2Fluid(Pythia8::Settings* s,
	Fluid* f, JamParticleData* jpd)
{
  fluid = f;
  particleTable = jpd;
  settings = s;
  isDebug=settings->mode("Check:Debug");
  optGauss = settings->mode("Hydro:optGaussSmear");  // mstc(146)
  //optGauss = 3;  // mstc(146)
  n_gaussp = 1; // mstc(138);
  widG = settings->parm("Hydro:gaussWidth");
  //widG = 0.5;
  widG2=2*widG*widG; // Gaussian width
  widCof=pow(1.0/(M_PI*widG2),1.5);
  gVolume=4.0/3.0*M_PI*widG*widG*widG;
  dCut=2.0;        // parc(147)
  optDens = 0; //mstc137=0;
  optConv = settings->mode("Hydro:optFluidConversion"); // mstc142=3;
  optPreHadron=1;// mstc89
  optEdmin=1;
  optCoreCorona = settings->mode("Hydro:optCoreCorona"); // =mstc(147)

  // GeV/fm^3
  minEnergyDensity = settings->parm("Hydro:FluidizationEnergyDensity");

  optTauEta = fluid->getOptTauEta();

  volume=fluid->getVol();

  //...Cut off parameter for Gaussian smearing.
  dX=fluid->dx();
  dY=fluid->dy();
  dZ=fluid->dz();
  iP=int(dCut/dX+0.5);
  iQ=int(dCut/dY+0.5);
  iR=int(dCut/dZ+0.5);
  maxX=fluid->xDim();
  maxY=fluid->yDim();
  maxZ=fluid->zDim();
 
  // Propagation by canonical momentum in RQMD mode
  optVectorPotential=settings->mode("MeanField:optVectorPotential");
  optVdot=settings->mode("MeanField:optVdot");
  if(optVectorPotential==1 && optVdot==0) optPropagate=1;
}

// Check fluid conversion after string fragmentation when
// strings are fragmented immediately after the collision
// with their formation time.
//***************************************************************************
void Particle2Fluid::checkFluidConversion(InterList* inter,
	vector<EventParticle*>& outgoing,double gtime)
{
  // all fluids were already converted into particles.
  if(fluid->isFluid()== -1)  return;

  vector<EventParticle*> outA, outB;
  for(int i=0; i<(int)outgoing.size(); i++) {
    if(outgoing[i]->getParent()==100) outA.push_back(outgoing[i]);
    else if(outgoing[i]->getParent()==200) outB.push_back(outgoing[i]);
  }

  if(outA.size()==0 && outB.size()==0) return;

  //TwoBodyInterList *inter2=dynamic_cast<TwoBodyInterList*>(inter);
  if(outA.size() > 0) {
    EventParticle* p = inter->getParticle(0);
    checkFluidConversion(p, outA, gtime);
  }
  if(outB.size() > 0) {
    EventParticle* p = inter->getParticle(1);
    checkFluidConversion(p, outB, gtime);
  }

  return;

  for(int i=0; i<(int)outgoing.size(); i++) {
     // if(outgoing[i]->getStatus() < 0 ) 
    cout << i << " status= "<<  outgoing[i]->getStatus()
	<< " pa= " << outgoing[i]->getParent()
	<< " id= " << outgoing[i]->getID()
	<< " t= " << outgoing[i]->getT()
	<< " tf= " << outgoing[i]->getTf()
	<< " x= " << outgoing[i]->getR()
	<<endl;
  }
  for(int i=0; i<(int)outA.size(); i++) {
     // if(outgoing[i]->getStatus() < 0 ) 
    cout << i << " statusA= "<<  outA[i]->getStatus()
	<< " paA= " << outA[i]->getParent()
	<< " id= " << outA[i]->getID()
	<< " t= " << outA[i]->getT()
	<< " tf= " << outA[i]->getTf()
	<<endl;
  }
  for(int i=0; i<(int)outB.size(); i++) {
     // if(outgoing[i]->getStatus() < 0 ) 
    cout << i << " statusB= "<<  outB[i]->getStatus()
	<< " paB= " << outB[i]->getParent()
	<< " id= " << outB[i]->getID()
	<< " t= " << outB[i]->getT()
	<< " tf= " << outB[i]->getTf()
	<<endl;
  }


}

//...Convert hadrons from string or resonance decays into fluid element.
//***************************************************************************
int Particle2Fluid::checkFluidConversion(EventParticle* ip,
	vector<EventParticle*>& outgoing,double gtime)
{
  // all fluids were already converted into particles.
  if(fluid->isFluid()== -1)  return -1;

  // First check energy density of this point. 
  int x,y,z;
  Vec4 r=0.0;
  if(optCoreCorona) {
    if(!checkEnergyDensity(ip,gtime,r,x,y,z)) return 0;
  }

  int iconv=0;
  // Loop over decaying particles.
  for(int i=0; i<(int)outgoing.size(); i++) {
    bool conv=outgoing[i]->isconv(optConv);

    // This particle will be converted into fluid after formation time.
    if(conv) {
      double ptime=outgoing[i]->getT();
      double dtp=gtime - ptime;
      if(dtp <  0.0) {
        outgoing[i]->setStatus(-1200);  // flag for fluization after formation time
        outgoing[i]->setLifeTime(ptime);
	iconv++;
	continue;
      } else {
        Vec4 r0 = outgoing[i]->propagate(gtime,optPropagate);
        conv = fluid->inside(r0,x,y,z);
      }
    }

    // This particle will be converted into fluid now.
    if(conv) {
      outgoing[i]->setStatus(-1000);
      iconv++;
    }

  } // end loop over outgoing particle.

  return iconv;

  // remove absorbed particles.
  /*
  auto result=std::remove_if(outgoing.begin(), outgoing.end(),
    [](EventParticle* p) { return p->getStatus() == -1;});
  for_each(result,outgoing.end(),[](EventParticle* p){delete p;});
  outgoing.erase(result,outgoing.end());
  */

}

// compute minimum energy density for the fluidzation.
//----------------------------------------------------------------------------
double Particle2Fluid::edminfluid(double binv)
{
  if(optEdmin==1) return minEnergyDensity;

  const double nq=2.5;
  const double pisq=9.86960440108936;
  const double hc=0.1973269788;
  const double Bagconst=pow(0.235,4)/(hc*hc*hc);

  // Compute minimum energy density (at T=0) at a baryon density with
  // massless QGP bag model.
  double x=abs(binv);
  double edminQGP=nq/(108*pisq)*pow(81*pisq/nq*x,4.0/3.0)*hc + Bagconst;
  return minEnergyDensity + edminQGP;
}

//----------------------------------------------------------------------------
// Check if this particle can converted into fluid or not.
bool Particle2Fluid::checkEnergyDensity(EventParticle* ip,double ctime,
	Vec4& r, int& x, int& y, int& z)
{
  // Check density.
  r = ip->propagate(ctime,optPropagate);
  if(fluid->inside(r,x,y,z) == 0) return false;

  double efluid=fluid->energyDensity(x,y,z)*HBARC;
  double bfluid=fluid->baryonDensity(x,y,z);
  double rho=0.0, rhob=0.0, einv=0.0;
  int icon=dens(ip,ctime,rho,rhob,einv,1);
  double bden=bfluid+rhob;
  if(icon == 0 && einv+efluid >= edminfluid(bden)) return true;
  return false;
}

// convert newly produced particles from string.
void Particle2Fluid::convert(Collision* event,EventParticle* ip,double ctime)
{
  if(isDebug) {
    cout << " p2fluid convert "<< ip->getID()
      << " name= "<< ip->getParticleDataEntry()->name() 
      <<endl;
  }

  Vec4 r = ip->propagate(ctime,optPropagate);
  int x=0,y=0,z=0;
  fluid->inside(r,x,y,z);

  p2fluid(ip,r,x,y,z);
  //delete ip;
  //event->removeInteraction(ip);

}

// Convert particle into fluid element after formation time.
//----------------------------------------------------------------------------
void Particle2Fluid::convertAform(Collision* event,EventParticle* ip,
	vector<EventParticle*>& outgoing,double ctime)
{
  bool convert=false;
  // Too late; All fluid was already converted to particle.
  //  if(mste(41).eq.-1) goto 1000
  if(fluid->isFluid()== -1)  convert=false;

  int x,y,z;
  Vec4 r=0.0;
  // If density is high enogh, this particle is converted into fluid.
  if(optCoreCorona) {
      convert = checkEnergyDensity(ip,ctime,r,x,y,z);
  } else {
    convert = fluid->inside(r,x,y,z);
  }

  if(convert) {
    //cout << " p2fluid from convertAform "<< ip->getStatus() << endl;
    p2fluid(ip,r,x,y,z);
    //event->removeInteraction(ip);

  // Cancel conversion.
  } else {

    ip->updateR(ctime);
    outgoing.push_back(new EventParticle(*ip));
    // Set resonance decay time.
    double m=ip->mCalc();
    double e=ip->getPe();
    double tdec = particleTable->lifeTime(ip->getParticleDataEntry(),m,e);
    outgoing.back()->setLifeTime(ctime+tdec);
  }

}

// Put particle into fluid elements.
//----------------------------------------------------------------------------
void Particle2Fluid::p2fluid(EventParticle* ep, Vec4& r, int ix, int iy, int iz,
    int opt)
{
  //double tau=1.0;
// if(opt_taueta.eq.1) tau=tau0 + it*dtt/3.0

  ep->setStatus(-1); // flag for dead particle.
  Vec4 pc = ep->getP();
  double bar = ep->baryon()/3.;
  double ch  = ep->charge()/3.;
  double str = ep->strange();
  //if(bar<0) { cout << "p2fluid bar<0 "<< bar << " id= "<< ep->getID() <<endl; }

  // Update total fluid energy-momentum and baryon number.
  fluid->updatePTotal(pc);
  fluid->updateBTotal(bar);
  fluid->updateCTotal(ch);
  fluid->updateSTotal(str);
  fluid->updateNC();

  if(optGauss == 0) {
    double vol=fluid->getVol();
    fluid->addU(ix,iy,iz, pc/vol/HBARC );
    fluid->addB(ix,iy,iz, bar/vol );
    if(opt==1) fluid->setLocalQuantity(ix,iy,iz);
    return;
  }

    // Gauss smearing.
    Vec4 u=0.0;
    if(optGauss == 3)  {
      u=pc/pc.mCalc();
    } else if(optGauss == 2)  {
      double emt=sqrt(pc[0]*pc[0] -pc[3]*pc[3]);
      u[3]=pc[3]/emt;
    }

    int jx1=max(0,ix-iP); int jx2=min(maxX-1,ix+iP);
    int jy1=max(0,iy-iQ); int jy2=min(maxY-1,iy+iQ);
    int jz1=max(0,iz-iR); int jz2=min(maxZ-1,iz+iR);
    double  enorm=0.0;
    for(int jx=jx1;jx<=jx2;jx++)
    for(int jy=jy1;jy<=jy2;jy++)
    for(int jz=jz1;jz<=jz2;jz++) {
      enorm += fluid->gauss(r, u, widG2, jx,jy,jz);
    }

    if(enorm == 0.0) {
      cout << "p2hydro::enorm? " <<  enorm 
          << " name= "<< ep->getParticleDataEntry()->name() 
	  << " x= "<< ep->getR()
	  << " r= "<< r
	  << " u= " << u 
	  << " status= " << ep->getStatus()
	  << " parent= " << ep->getParent()
	  <<endl;
      cout << " maxX= " << maxX << " ny "<< maxY<< " nz= " << maxZ <<endl;
      cout << " ip= " << iP << " iq= "<< iQ << " ir= " << iR <<endl;
      cout << " ix= " << ix << " iy= "<< iy << " iz= " << iz <<endl;
      cout << " jx1= " << jx1 << " jx2= "<< jx2 <<endl;
      cout << " jy1= " << jy1 << " jy2= "<< jy2 <<endl;
      cout << " jz1= " << jz1 << " jz2= "<< jz2 <<endl;
      cout << "vol= " << fluid->getVol()<<endl;
      exit(1);
    }

    // Put energy-momentum and baryon number into the fluid cells.
    for(int jx=jx1;jx<=jx2;jx++)
    for(int jy=jy1;jy<=jy2;jy++)
    for(int jz=jz1;jz<=jz2;jz++) {
      double ex=fluid->gauss(r, u, widG2,jx,jy,jz)/enorm/volume;
      if(ex < 1e-12) continue;
      fluid->addU(jx,jy,jz, pc/HBARC*ex );
      fluid->addB(jx,jy,jz, bar*ex );
      if(opt==1) fluid->setLocalQuantity(jx,jy,jz);
    }
}

//...Purpose: compute energy-momentum tensor and hydrodynamics velocity.
//***********************************************************************
//int  Particle2Fluid::dens(list<EventParticle*>& plist,EventParticle* i1,
int  Particle2Fluid::dens(EventParticle* i1,
	double ctime,double& rho,double& rhob,double& einv,int iopt)
{
  rho=0.0;
  rhob=0.0;
  einv=0.0;
  Vec4 u=0.0;
  Vec4 r0 = i1->propagate(ctime,optPropagate);
  int ix0, iy0,iz0;
  bool isfluid=fluid->inside(r0,ix0,iy0,iz0);
  if(!isfluid) {
      cout << "Particle2Fluid::dens outside the fluid"<< r0 << endl;
      return 1;
  }

  double cur[4]={0.0}, curb[4]={0.0}, tens[4][4]={0.0};
  int ncount=0;

  // copy particle list.
  list<EventParticle*> partlist;
  CascBox* box=i1->box();
  if(box) {
    for(auto& b : box->getNeighbors2()) {
      if(b->getParticles().size()>0) partlist.insert(partlist.begin(),b->getParticles().begin(),b->getParticles().end());
    }
  } else {
    partlist.insert(partlist.begin(),particleList->begin(),particleList->end());
  }

//....Loop over all particles
  for(auto ip : partlist) {
    if(ip ==  i1) continue;
    if(ip->getMass()  < 1e-5) continue;

    // not yet collide
    if(abs(ip->getNColl()) <=  1) continue;  

    //pre-formed hadrons.
    if(optPreHadron == 0 && ip->getTf() > ctime) continue;

    double dt=ctime - ip->getT();
    if(iopt == 1 && dt <  0.0) continue;

    Vec4 r1 = ip->propagate(ctime,optPropagate);
    int ix,iy,iz;
    if(!fluid->inside(r1,ix,iy,iz)) continue;

    double bar=ip->baryon()/3.0;
    double facq=1.0;
    double facb=bar;
    if(iopt == 1 && ip->getTf() > ctime) {
	facq = ip->qFactor();
	facb = facq*bar;
    }

    Vec4 pv = (ip)->getP();
    /*
    double den=0.0;
    // Gauss smearing.
    if(optGauss == 0) {
	den = gaussSmear(pv,r1-r0);
    } else {
      if(ix == ix0 && iy ==  iy0 && iz == iz0) den += 1.0/volume;
      else continue;
    }
    */

    double den = gaussSmear(pv,r1-r0);
    double bden = den*facb;

    //...Compute current and energy-momentum tensor.
    for(int im=0; im<4; im++) {
      cur[im]  += pv[im]/pv[0]*den;
      curb[im] += pv[im]/pv[0]*bden;
      for(int ik=0; ik<4;ik++) {
        tens[im][ik] += pv[im]*pv[ik]/pv[0]*den;
      }
    }

    ncount++;

  } // End loop over particles.


  if(ncount ==  0)   return 1;
  if(cur[0] <  1e-7) return 1;

  double  cc=cur[0]*cur[0] - cur[1]*cur[1]-cur[2]*cur[2]-cur[3]*cur[3];
  if(cc < 1e-7) return 1;

  // Compute local energy density and baryon density by using EoS.
  if(optDens == 1) {
    Vec4 uu(tens[1][0], tens[2][0], tens[3][0], tens[0][0]);
    double pre, sva, vel;
    fluid->thermal(uu,curb[0],einv,rhob,pre,sva,vel);
    double tau=1.0;
    double fg = uu[0] + tau*pre;
    u[1] = uu[1]/fg;
    u[2] = uu[2]/fg;
    u[3] = uu[3]/fg;
    u[0]=sqrt(1.0-u[1]*u[1]-u[2]*u[2]-u[3]*u[3]);
    return 0;
  }

  cc=sqrt(cc);
  u[0]=cur[0]/cc;
  u[1]=cur[1]/cc;
  u[2]=cur[2]/cc;
  u[3]=cur[3]/cc;

//     if(mstc(47).eq.2) call landauframe(u,tens)

//...Lorentz invariant Scalar Number Density
  rho=cur[0]*u[0]-cur[1]*u[1]-cur[2]*u[2]-cur[3]*u[3];

//...Lorentz invariant Baryon density
  rhob=curb[0]*u[0]-curb[1]*u[1]-curb[2]*u[2]-curb[3]*u[3];

//...Lorentz invariant pressure and energy density
  int g[4][4]={{-1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,1}};
  int gg[4][4]={{1,1,1,-1},  {1,1,1,-1}, {1,1,1,-1}, {-1,-1,-1,1}};
  int dl[4][4]={{1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}};

  double gam=u[0];
  double vz=u[3]/u[0];
  double gamz=1.0/sqrt(1.0-vz*vz);
  double z[4];
  z[0]=gamz*vz;
  z[3]=gamz;
  z[2]=0.0;
  z[1]=0.0;
  double pfree=0.0;
  double pxx=0.0;
  double pyy=0.0;
  double pzz=0.0;
  double pperp=0.0;
  double plong=0.0;
  for(int i=0;i<4;i++)
  for(int j=0;j<4;j++) {
    double tmp=gg[i][j]*u[i]*u[j];
    einv += tens[i][j]*tmp;                     // energy density
    pfree += - 1.0/3.0*tens[i][j]*(g[i][j]-tmp); // pressure
    double pl=gg[i][j]*z[i]*z[j];
    plong += tens[i][j]*pl;
    pperp -= 0.5*tens[i][j]*(g[i][j]-tmp+pl);

//...Txx,Tyy,Tzz at the local rest frame.
    double ci=-1.0;
    double cj=-1.0;
    if(i > 0) ci=u[i]/(1+gam);
    if(j > 0) cj=u[j]/(1+gam);
    pxx += tens[i][j]*(dl[1][i]+u[1]*ci)*(dl[1][j]+u[1]*cj);
    pyy += tens[i][j]*(dl[2][i]+u[2]*ci)*(dl[2][j]+u[2]*cj);
    pzz += tens[i][j]*(dl[3][i]+u[3]*ci)*(dl[3][j]+u[3]*cj);
  }

  //cout << "dens einv= "<< einv << " T00= "<< tens[0][0] << " d= "<< rho/0.168 << " nv= "<< partlist.size()<<endl;
  //partlist.clear();
  //cin.get();

  return 0;

}

//------------------------------------------------------------
double Particle2Fluid::gaussSmear(const Vec4& pv, const Vec4& dr)
{
  const double xg3[]={0.0, -.77459666924148337703,.77459666924148337703};
  const double wg3[]={.8888888888888888,.55555555555555,.555555555555555};
  const int optg=2;

  Vec4 u=0.0;
  double gam=1.0;
  double pe=pv.e();
  if(optGauss == 0 || optGauss == 3)  {
    double  emd= pv.mCalc();
    u=pv/emd;
    gam=pe/emd;
    if(optGauss == 0) {
	return cmDistanceSquare(dr,u) < widG ? 1.0/gVolume: 0.0;
    }
  } else if(optGauss == 2) {
    double pz=pv.pz();
    double emt=sqrt(pe*pe -pz*pz);
    u[3]=pz/emt;
    gam=pe/emt;
  }

  if(optg == 1) {
    double xtra = cmDistanceSquare(dr,u);
    return widCof*gam*exp(xtra/widG2);

  } else {

    Vec4 dr2 = Vec4(dX/2, dY/2, dZ/2, 0.0);
    double den=0.0;
    for(int ig=0;ig<3;ig++)
    for(int jg=0;jg<3;jg++)
    for(int kg=0;kg<3;kg++) {
      Vec4 r3 = dr - xg3[ig]*dr2;
      double xtra=cmDistanceSquare(r3,u);
      den += exp(xtra/widG2)*wg3[ig]*wg3[jg]*wg3[kg]/8;
    }
    return widCof*gam*den;

  }


}

//...Convert jam particles into fluid elements when a particle is in the dense region.
void Particle2Fluid::convertAtDenseRegion(double ctime,Collision* event)
{
  auto plist = &event->plist;
  if(plist->size() == 0) return;

  // First find the particles which will be converted into fluid.
  // loop over all particles.
  for(auto& ip : *plist) {
    if(abs(ip->getNColl()) <=  1) continue;  // not yet interact
    if((ip)->getMass()  < 1e-5) continue;
    Vec4 r1 = ip->propagate(ctime,optPropagate);
    int x,y,z;
    if(!fluid->inside(r1,x,y,z)) continue;
    if(fluid->energyDensity(x,y,z)*HBARC <1e-10) continue;

    if(optCoreCorona == 2) {
        if((ip)->getTf() > ctime) continue;  //...pre-formed hadrons.
    } else if(optCoreCorona == 3) {
        if((ip)->getT() > ctime) continue;  // not formed hadrons.
        if((ip)->getTf() > ctime) continue;  // pre-formed hadrons.
    }

    if(checkEnergyDensity(ip,ctime,r1,x,y,z)) ip->setStatus(-1500);

    }

  int nv=0;
  // Convert particles into fluid.
  /*
  for(auto ip=plist->begin();ip !=plist->end();ip++) {
    if((*ip)->getStatus() != -1500) continue;
    Vec4 rp = (*ip)->propagate(ctime,optPropagate);
    fluid->inside(rp,x,y,z);
    p2fluid(*ip,rp,x, y, z, 1);
    event->removeInteraction(*ip);
    nv++;
  }
  */

  // Convert particles into fluid.
  list<EventParticle*>::iterator first=plist->begin();
  while(first != plist->end()) {
    list<EventParticle*>::iterator next=first;
    ++next;
    if((*first)->getStatus() == -1500) {
      Vec4 rp = (*first)->propagate(ctime,optPropagate);
      int x=0,y=0,z=0;
      fluid->inside(rp,x,y,z);
      p2fluid(*first,rp,x, y, z, 1);
      if(isDebug) cout << "p2fluid at dense region t= "<< 
	rp[0] << " id= "<< (*first)->getID()<<endl;
      //event->removeInteraction(*first);
      //event->eraseParticle(*first);
      event->erasePList(first);
      nv++;
    }
    first=next;
  }

  // delete particles.
  //event->deleteParticle(-1500);

  // Update collision list.
  //event->makeCollisionList(ctime,ftime);

} 

//...Convert all jam particles into fluid elements. jam2hydroall
void Particle2Fluid::convertAll(double ctime,double ftime,Collision* event)
{
  //list<EventParticle*> plist = event->plist;
  auto plist = &event->plist;

  if(plist->size() == 0) {
    cout << "(Particle2Fluid::convertAll: nv= " << plist->size() <<endl;
    exit(1);
  }

  int isel=0;
  if(optCoreCorona == 2) isel=1;

  double tau0=ctime;
  double tau=1.0;
  if(optTauEta == 1) tau=tau0;
  int x=0, y=0, z=0;

  // First find the particles which will be converted into fluid.
  // loop over all particles.
  for(auto& ip : *plist) {
    //if(abs((ip)->getNColl()) <=  1) continue;  // not yet interact
    if((ip)->getMass()  < 1e-5) continue;
    Vec4 r1 = (ip)->propagate(ctime,optPropagate);
    if(!fluid->inside(r1,x,y,z)) continue;

    if(optCoreCorona >= 1) {
      if(abs((ip)->getNColl()) <=  1) continue;  // Exclude spectator
      if(optCoreCorona == 2) {
        if((ip)->getTf() > ctime) continue;  //...pre-formed hadrons.
      } else if(optCoreCorona == 3) {
        if((ip)->getT() > ctime) continue;  // not formed hadrons.
        if((ip)->getTf() > ctime) continue;  // pre-formed hadrons.
      }
      double rho=0.0, rhob=0.0, einv=0.0;
      if(dens(ip,ctime,rho,rhob,einv,isel) !=0) continue;
      if(einv >= edminfluid(rhob)) (ip)->setStatus(-1500);
    } else {
      (ip)->setStatus(-1500);
    }

  }

  fluid->reset();

  int nv=0;
  // Convert particles into fluid.
  /*
  for(auto ip=plist->begin();ip !=plist->end();ip++) {
    if((*ip)->getStatus() != -1500) continue;
    Vec4 rp = (*ip)->propagate(ctime,optPropagate);
    fluid->inside(rp,x,y,z);
    p2fluid(*ip,rp,x, y, z, 0);
    event->removeInteraction(*ip);
    nv++;
  }
  */

  list<EventParticle*>::iterator first=plist->begin();
  while(first != plist->end()) {
    list<EventParticle*>::iterator next=first;
    ++next;
    if((*first)->getStatus() == -1500) {
      Vec4 rp = (*first)->propagate(ctime,optPropagate);
      fluid->inside(rp,x,y,z);
      p2fluid(*first,rp,x, y, z, 0);
      //event->removeInteraction(*first);
      //event->eraseParticle(*first);
      event->erasePList(first);
      nv++;
    }
    first=next;
  }

  /*
  auto& it = plist->begin();
  while (it !=plist->end()) {
    if((*it)->getStatus() == -1500) {
      Vec4 rp = (*it)->propagate(ctime,optPropagate);
      fluid->inside(rp,x,y,z);
      p2fluid(*it,rp,x, y, z, 0);
      it=event->removeInteraction(*it);
      nv++;
    } else ++it;
  }
  */

  if(nv == 0) {
    cout << "Particle2Fluid::convertAll no particle is converted into fluid " << nv <<endl;
    exit(1);
  }

  if(isDebug)
  cout << nv << "  particles are converted" <<endl;

  // Compute local valuables at each fluid element.
  for(int ix=0;ix<maxX;ix++)
  for(int iy=0;iy<maxY;iy++)
  for(int iz=0;iz<maxZ;iz++) {
    if(fluid->U(ix,iy,iz,0) < 1e-12) {
      fluid->reset_fluid_element(ix,iy,iz);
      continue;
    }
    fluid->setLocalQuantity(ix,iy,iz,tau);
  }

  if(isDebug)
  cout << "before particle size = "<< plist->size() << endl;

  // delete particles.
  //event->deleteParticle();

  if(isDebug)
  cout << "after particle size = "<< plist->size() << endl;

  // Update collision list.
  event->makeCollisionList(ctime,ftime);

  if(isDebug) {
    double etot,pxtot,pytot,pztot,btot;
    fluid->hydroTotalEnergy(0,ctime, etot,pxtot,pytot,pztot,btot);
  }

} 

} // end namespace jam2
