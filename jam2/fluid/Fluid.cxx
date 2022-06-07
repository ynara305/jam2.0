#include <jam2/fluid/Fluid.h>

namespace jam2 {

using namespace std;

Fluid::Fluid(Pythia8::Settings* s,EoS* eos0, FreezeOut* f)
{
  settings = s;
  eos = eos0;
  freezeout=f;
  eFreezeOut=freezeout->getFreezeOutEnergyDensity(); // 1/fm
  optFreezeOut= settings->mode("Hydro:optFreezeOut");

  maxX= settings->mode("Hydro:nx");
  maxY= settings->mode("Hydro:ny");
  maxZ= settings->mode("Hydro:nz");
  origX=maxX/2;
  origY=maxY/2;
  origZ=maxZ/2;
  volume = maxX*maxY*maxZ;
  dT= settings->parm("Hydro:dt");
  dX= settings->parm("Hydro:dx");
  dY= settings->parm("Hydro:dy");
  dZ= settings->parm("Hydro:dz");
  vol = dX*dY*dZ;
  volH=vol*HBARC;

  passTime=5.0;
  secondOrder = 0;
  tau0=0.0;
  optTimeLike = settings->mode("Hydro:optTimeLike");
  optRescale=0;
  optTauEta=0;
  jobMode=settings->mode("Check:Debug");

  for(int z=0;z<maxZ;z++)
  for(int y=0;y<maxY;y++)
  for(int x=0;x<maxX;x++) {
    fluid.push_back(new FluidElement());
  }

  fluidx.resize(maxX);
  fluidy.resize(maxY);
  fluidz.resize(maxZ);

  solverx = new RSolverHLLE(maxX, dX, eos, optTimeLike,optTauEta);
  solvery = new RSolverHLLE(maxY, dY, eos, optTimeLike,optTauEta);
  solverz = new RSolverHLLE(maxZ, dZ, eos, optTimeLike,optTauEta);

  reset();

  // Propagation by canonical momentum in RQMD mode
  optVectorPotential=settings->mode("MeanField:optVectorPotential");
  optVdot=settings->mode("MeanField:optVdot");
  if(optVectorPotential==1 && optVdot==0) optPropagate=1;

  particleDens = new ParticleDensity(settings);

}

Fluid::~Fluid()
{
  for(auto& f:fluid) delete f;
  delete solverx;
  delete solvery;
  delete solverz;
  delete particleDens;
}

int Fluid::evolution(int step,double gtime, Collision* event)
{
  //double htime = step* dT + tau0;

  int it=fluidTimeStep;

  // 3D evolution by the Operator splitting.
  if(!secondOrder) {
    evolution3d(it,gtime);
  } else {
    evolution3d2(it,gtime);
  }

  if(jobMode >= 1) {
    double etot,pxtot,pytot,pztot,btot;
    hydroTotalEnergy(step,gtime, etot,pxtot,pytot,pztot,btot);
  }

  // isoenergy density.
  if(optFreezeOut % 10 == 2) {
    freezeout->isoenergyFreezeout(fluid,it,gtime);
  } else if(optFreezeOut == 11) {
    freezeout->isochronousFreezeout(fluid,gtime,0);
  }

  // check: number of surface which is greater than the particlization energy density.
  int check = checkout(gtime);

  if(jobMode >= 2) {
    //if(it % 5 == 0) outputhzx(step);
    outputhzx(step);
  }

  fluidTimeStep++;

  return check;

}

//************************************************************
void Fluid::reset()
{
  for(int i=0;i<volume;i++) {
      fluid[i]->reset();
  }
  uzero=0.0; uout=0.0;
  bzero=0.0;
  bout=0.0;
  baryonTot=0.0;
  chargeTot=0.0;
  strangeTot=0.0;
  entroTot=0.0;
  pTotal=0.0;
  fluidTimeStep=0; 
  nConversion=0;
  freezeout->clear();
  freezeout->reset(fluid);
}

//************************************************************
//....hydrodynamical evolution in full 3-dimension by the operator
//...splitting method.
void Fluid::evolution3d2(int it,double tau)
{
  int is = it % 6;
  int iz=0;
  if(optTauEta >= 1) iz=1;
  double dt0=dT;
  double dt;

  switch(is) {
    case 0:
      dt=0.5*dt0;                 update_z(tau,dt,0);
      dt=0.5*dt0; tau=tau+dt/3.0; update_x(tau,dt,0);
      dt=dt0;   ; tau=tau+dt/3.0; update_y(tau,dt,0);
      dt=0.5*dt0; tau=tau+dt/3.0; update_x(tau,dt,0);
      dt=0.5*dt0; tau=tau+dt/3.0; update_z(tau,dt,iz);
      break;
    case 1:
      dt=0.5*dt0;                 update_x(tau,dt,0);
      dt=0.5*dt0; tau=tau+dt/3.0; update_z(tau,dt,0);
      dt=    dt0; tau=tau+dt/3.0; update_y(tau,dt,0);
      dt=0.5*dt0; tau=tau+dt/3.0; update_z(tau,dt,0);
      dt=0.5*dt0; tau=tau+dt/3.0; update_x(tau,dt,iz);
      break;
    case 2:
      dt=0.5*dt0;                 update_y(tau,dt,0);
      dt=0.5*dt0; tau=tau+dt/3.0; update_x(tau,dt,0);
      dt=dt0;     tau=tau+dt/3.0; update_z(tau,dt,0);
      dt=0.5*dt0; tau=tau+dt/3.0; update_x(tau,dt,0);
      dt=0.5*dt0; tau=tau+dt/3.0; update_y(tau,dt,iz);
      break;
    case 3:
      dt=0.5*dt0;                 update_z(tau,dt,0);
      dt=0.5*dt0;  tau=tau+dt/3.0;update_y(tau,dt,0);
      dt=dt0;      tau=tau+dt/3.0;update_x(tau,dt,0);
      dt=0.5*dt0;  tau=tau+dt/3.0;update_y(tau,dt,0);
      dt=0.5*dt0;  tau=tau+dt/3.0;update_z(tau,dt,iz);
      break;
    case 4:
      dt=0.5*dt0;                update_x(tau,dt,0);
      dt=0.5*dt0; tau=tau+dt/3.0;update_y(tau,dt,0);
      dt=dt0;     tau=tau+dt/3.0;update_z(tau,dt,0);
      dt=0.5*dt0; tau=tau+dt/3.0;update_y(tau,dt,0);
      dt=0.5*dt0; tau=tau+dt/3.0;update_x(tau,dt,iz);
      break;
    case 5:
      dt=0.5*dt0;                 update_y(tau,dt,0);
      dt=0.5*dt0; tau=tau+dt/3.0; update_z(tau,dt,0);
      dt=dt0;     tau=tau+dt/3.0; update_x(tau,dt,0);
      dt=0.5*dt0; tau=tau+dt/3.0; update_z(tau,dt,0);
      dt=0.5*dt0; tau=tau+dt/3.0; update_y(tau,dt,iz);
      break;
   default:
      cout << " error in Fluid::evolution3d2 is = " << is << endl; 
      exit(1);
  }

}

//************************************************************
//....hydrodynamical evolution in full 3-dimension by the operator
//...splitting method.
void Fluid::evolution3d(int it,double tau)
{
  int s1=0, s2=0, s3=0;
  double tau1=1.0;
  double tau2=1.0;
  double tau3=1.0;

  double dt = dT;
  if(optTauEta >= 1) {
//  s1=3; s2=3; s3=1  // Include source term for all direction
    s1=0; s2=0; s3=1; // Include source term only for z-direction
    tau1=tau;
    tau2=tau+dt/3.0;
    tau3=tau+2.0*dt/3.0;
  }

  int is = it % 6;
  switch(is) {
    case 0:
      update_z(tau1,dt,s3);
      update_x(tau2,dt,s1);
      update_y(tau3,dt,s2);
      break;
    case 1:
      update_x(tau1,dt,s1);
      update_z(tau2,dt,s3);
      update_y(tau3,dt,s2);
      break;
    case 2:
      update_y(tau1,dt,s2);
      update_x(tau2,dt,s1);
      update_z(tau3,dt,s3);
      break;
    case 3:
      update_z(tau1,dt,s3);
      update_y(tau2,dt,s2);
      update_x(tau3,dt,s1);
      break;
    case 4:
      update_x(tau1,dt,s1);
      update_y(tau2,dt,s2);
      update_z(tau3,dt,s3);
      break;
    case 5:
      update_y(tau1,dt,s2);
      update_z(tau2,dt,s3);
      update_x(tau3,dt,s1);
      break;
   default:
      cout << " error in Fluid::evolution3d is = " << is << endl; 
      exit(1);
  }

}

//******************************************************************
void Fluid::update_x(double tau,double dt,int is)
{
  double tau1=1.0;
  double tau3=1.0;
  if(optTauEta > 0) {
    tau1=tau;
    tau3=tau+dt/3.0;
  }

  vector<FluidElement*> fluidx0;
  for(int iy=N; iy< maxY-N;iy++)
  for(int iz=N; iz< maxZ-N;iz++) {

    fluidx.clear();
    fluidx0.clear();
    for(int ix=0; ix< maxX;ix++) {
      fluidx.push_back(site(ix,iy,iz));
      fluidx0.push_back(fluidx[ix]);
    }

    /*
    for(int ix=0; ix< maxX;ix++) {
      if(fluidx[ix]->u(0)>0.0)
      //if(site(ix,iy,iz)->u(0)>0.0)
      //cout << " f= "<< site(ix,iy,iz)->u() << " fx= "<< fluidx[ix]->u() <<endl;
      cout <<ix << scientific << setprecision(16) 
	<< " before f= "<< site(ix,iy,iz)->u(1)
	<< " fx= "<< fluidx[ix]->u(1)
	<< " fx0= "<< fluidx0[ix]->u(1)
	<<endl;
    }
    */

    solverx->solve(fluidx,dt,tau1,1,is);
    update_local_val(fluidx,tau3,maxX);

    /*
    for(int ix=0; ix< maxX;ix++) {
      if(fluidx[ix]->u(0)>0.0)
      //if(site(ix,iy,iz)->u(0)>0.0)
      //cout << " f= "<< site(ix,iy,iz)->u() << " fx= "<< fluidx[ix]->u() <<endl;
      cout <<ix << scientific << setprecision(8) 
	<< " before f= "<< site(ix,iy,iz)->u(0)
	<< " fx= "<< fluidx[ix]->u(0)
	<< " fx0= "<< fluidx0[ix]->u(0)
	<<endl;
    }
    */


  }
    //cout << "size= "<< fluidx.size()<<endl;
    //cin.get();

  checkboundary_x();
}

//******************************************************************
void Fluid::update_y(double tau,double dt, int is)
{
  double tau1=1.0;
  double tau3=1.0;
  if(optTauEta > 0) {
    tau1=tau;
    tau3=tau+dt/3.0;
  }

  for(int ix=N; ix< maxX-N;ix++)
  for(int iz=N; iz< maxZ-N;iz++) {

    fluidy.clear();
    for(int iy=0; iy< maxY;iy++)
      fluidy.push_back(site(ix,iy,iz));
    solvery->solve(fluidy,dt,tau1,2,is);
    update_local_val(fluidy,tau3,maxY);

  }

  checkboundary_y();

}

//******************************************************************
void Fluid::update_z(double tau,double dt, int is)
{
  double tau1=1.0;
  double tau3=1.0;

  if(optTauEta >0 ) {
    tau1=tau;
    tau3=tau+dt/3.0;
  }

  for(int ix=N;ix<maxX-N;ix++)
  for(int iy=N;iy<maxY-N;iy++) {
    fluidz.clear();
    for(int iz=0; iz< maxZ;iz++)
      fluidz.push_back(site(ix,iy,iz));
    solverz->solve(fluidz,dt,tau1,3,is);
    update_local_val(fluidz,tau3,maxZ);
  }

  checkboundary_z();
}

void Fluid::checkboundary_x()
{
  double ch=1.0; double sh=0.0;
  int n=2;
  for(int i=0;i<2;i++) {
    int ix=(maxX-1)*i;
    for(int iz=n;iz<maxZ-n;iz++) {
      if(optTauEta != 0) {
        double h=dZ*(iz-origZ);
	ch = cosh(h); sh=sinh(h);
      }
      for(int iy=n;iy<maxY-n;iy++) {
        FluidElement *f = site(ix,iy,iz);
        bout +=  f->U(4);
        uout[1] +=  f->u(1);
        uout[2] +=  f->u(2);
        uout[0] +=  ch*f->u(0) + sh*f->u(3);
        uout[3] +=  sh*f->u(0) + ch*f->u(3);
        f->resetU();
      }
    }
  }

}

void Fluid::checkboundary_y()
{
  double ch=1.0; double sh=0.0;
  int n=2;
  for(int i=0;i<2;i++) {
    int iy=(maxY-1)*i;
    for(int iz=n;iz<maxZ-n;iz++) {
      if(optTauEta != 0) {
        double h=dZ*(iz-origZ);
	ch = cosh(h); sh=sinh(h);
      }
      for(int ix=n;ix<maxX-n;ix++) {
        FluidElement *f = site(ix,iy,iz);
        bout +=  f->U(4);
        uout[1] +=  f->u(1);
        uout[2] +=  f->u(2);
        uout[0] +=  ch*f->u(0) + sh*f->u(3);
        uout[3] +=  sh*f->u(0) + ch*f->u(3);
        f->resetU();
      }
    }
  }

}

void Fluid::checkboundary_z()
{
  double ch=1.0; double sh=0.0;
  int n=2;
  for(int i=0;i<2;i++) {
    int iz=(maxZ-1)*i;
    if(optTauEta != 0) {
      double h=dZ*(iz-origZ);
      ch = cosh(h); sh=sinh(h);
    }
    for(int ix=n;ix<maxX-n;ix++) {
      for(int iy=n;iy<maxY-n;iy++) {
        FluidElement *f = site(ix,iy,iz);
        bout +=  f->U(4);
        uout[1] +=  f->u(1);
        uout[2] +=  f->u(2);
        uout[0] +=  ch*f->u(0) + sh*f->u(3);
        uout[3] +=  sh*f->u(0) + ch*f->u(3);
        f->resetU();
      }
    }
  }

}

//-----------------------------------------------------------------------
//...Update thermo dynamical quantities from uu
void Fluid::update_local_val(vector<FluidElement*>& fluid0,double tau, int maxx)
{
  //const double mim = 0.00001;

//     uo(:) = uo(:) + uu(:,1)
//     uo(:) = uo(:) + uu(:,mxx-1)
//     uo(:) = uo(:) + uu(:,0)
//     uo(:) = uo(:) + uu(:,mxx)
//     uu(:,0)=0.0; uu(:,1)=0.0; uu(:,mxx)=0.0; uu(:,mxx-1)=0.0
//-----------------------------------------------------------------------      
  if(optRescale) {
    double b0=0.0;
    for(int i=1;i<maxx-1;i++) {
      uzero += fluid0[i]->scaleU(b0,optRescale);
      bzero += b0;
    }
  }

  for(int ix=1;ix<maxx-1;ix++) {

    double u4 = fluid0[ix]->U(4)/tau; // baryon density
    Vec4 uh = fluid0[ix]->u()/tau;

    if(uh[0] <= 1e-15) { uh=0.0; u4=0.0; }

    // baryon density is negative  <--- SHOULD BE CHANGED!!!
    if(u4 < 0.0) { 
      //cout << "baryon density is negative "<< u4 << " ix= "<< ix
      //   << " det= "<< uh.m2Calc() << " uh0= "<< uh[0] <<endl;
      fluid0[ix]->setB(0.0); u4=0.0; 
    }
  
    double det = uh.m2Calc();
    if (uh[0] < 0.0 || det <= 0.0) {

      uzero += fluid0[ix]->u();
      bzero += fluid0[ix]->b();
      fluid0[ix]->reset(); 

    } else {
      double ene,nba,pre,sva,vel;
      thermal(uh,u4,ene,nba,pre,sva,vel);
      fluid0[ix]->ed(ene);
      fluid0[ix]->bd(nba);
      fluid0[ix]->pr(pre);
      double fg = uh[0] + tau*pre;

      /*
      double gam = 1.0-uh.pAbs2()/fg/fg;
      if((fg > 0.0) && (gam > 0.0))
          gam = sqrt(1.0/gam);
        else
          gam = 1.0;
      fluid0[ix]->bd(u4/gam/tau);
      cout << "nba= "<< nba << " u4/g= "<< u4/gam<<endl;
      */


      if(fg != 0.0) {
        fluid0[ix]->vx(uh[1]/fg);
        fluid0[ix]->vy(uh[2]/fg);
        fluid0[ix]->vz(uh[3]/fg);
        fluid0[ix]->cs(sva);
      } else {
        fluid0[ix]->vx(0.0);
        fluid0[ix]->vy(0.0);
        fluid0[ix]->vz(0.0);
        fluid0[ix]->cs(0.0);
      }

      double vx=fluid0[ix]->vx();
      double vy=fluid0[ix]->vy();
      double vz=fluid0[ix]->vz();
      double v2 = vx*vx + vy*vy + vz*vz;
      if(v2 > 1.0) {
	  cout << " v> 1? " << scientific << v2
	      << " vx= " << vx
	      << " vy= " << vy
	      << " vz= " << vz
	      <<endl;
	  cout << " fg = " << fg
	      << " pre= " << pre
	      << " u0= " << uh[0]
	      << " vx= " << uh[1]/fg 
	      << " vy= " <<  uh[2]/fg
	      << " vz= " <<  uh[3]/fg
	      <<endl;
	  exit(1);
      }

    }


  } // end loop over grid.

}

//-----------------------------------------------------------------------      
void Fluid:: setLocalQuantity(int ix,int iy,int iz,double tau)
{
  FluidElement *f = site(ix,iy,iz);
  Vec4 u = f->u()/tau;
  double u4 = f->U(4)/tau;
  double ene,nba,pre,sva,vel;

  thermal(u,u4,ene,nba,pre,sva,vel);

  double fg = u[0] + tau*pre;
  f->vx( u[1]/fg );
  f->vy( u[2]/fg );
  f->vz( u[3]/fg );

  double va=u.pAbs()/fg;
  if(va > 1.0) {
    cout << "setlocalq v>=1? " << va << " fg= "<< fg << " vel=" << vel <<endl;
  }
  f->ed(ene);
  f->bd(nba);
  f->pr(pre);
  f->cs(sva);

}

//-----------------------------------------------------------------------      

//...Compute thermodynamical values from uu()
void Fluid::thermal(Vec4& uu,double u4,double& elc,double& nblc,double& plc,
	double& cslc, double& v)
{
  const int mxit=30;
  const double eps=1e-8;

  int icon=0;
  //double vv1 = uu[1];
  //double vv2 = uu[2];
  //double vv3 = uu[3];
  //double m = sqrt(vv1*vv1+vv2*vv2+vv3*vv3);
  double vv4 = uu[0];
  double vv5 = max(0.0,u4);
  double m = uu.pAbs();

  if(m == 0.0) {
    nblc=vv5;
    elc=vv4;
    plc=eos->getPressure(nblc,elc);
    cslc=eos->getCs(nblc,elc);
    return;
  }

  if(vv4 < 1e-12) {
    elc=0.0;
    nblc=0.0;
    plc=0.0;
    cslc=sqrt(1.0/3.0);
    return;
  }

//------------------------------------------------------
  icon = iteratev(vv4,vv5,m,elc,nblc,plc,v);

  cslc=eos->getCs(nblc,elc);
  double fe=vv4+plc;
  double g1=1.0-m*m/fe/fe;
  if(g1 <= 0.0) return;
  if(icon == 0) return;

//------------------------------------------------------
  double vp = 1.0;
  double vm = 0.0;
 
  for(int iv=0;iv< mxit;iv++) {
    v = (vp+vm)/2.;
    elc = vv4-m*v;
    nblc = vv5*sqrt(1.0-v*v);
    if(elc < 0.0) {
        cout << "thermal: elc<0" << elc <<endl;
	exit(1);
    }
    plc = eos->getPressure(nblc,elc);
    cslc= eos->getCs(nblc,elc);
    double fv = v-m/(vv4+plc);
    if(fv > 0.0) {
      vp = v;      
    } else {
      vm = v;
    }

    if(abs(vp-vm) < eps) return;

  }

  cout << "not converge? v=" << vp << endl;
  exit(1);

}

//***********************************************************************
int Fluid::iteratev(double u0,double u4,double m,
	double& elc,double& nblc,double& plc,double& vf)
{
  const double eps=1e-7;
  const int mxit=20;

  vf = 0.0;
  int iv = 0;
  do {
    double v = vf;
    elc = u0-m*v;
    if(elc < 0.0)  return 2;
    nblc = u4*sqrt(max(0.0, 1.0 - v*v) );
    plc = eos->getPressure(nblc,elc);
    vf = m / (u0  + plc );
    if (abs((v-vf)/vf) < eps) return 0;

  } while(++iv < mxit);

  return 1;

}

//***********************************************************************
void Fluid::fluidFraction(const double tau,Collision* event,ostream& out,int opt)
{
  double zmin=1.0;

  auto plist = event->plist;
  double epart=0.0;
  for(auto& p : plist) {
    if(p->getStatus() <  0) continue;
    if(abs(p->getNColl()) <=  1) continue;  // not yet interact
    if(p->getMass()  < 1e-5) continue;
    if(tau - p->getT() < 0.0) continue;
    Vec4 r = p->propagate(tau,optPropagate);
    //if(inside(r,x,y,z)) continue;
    if(abs(r[3]) < zmin) epart += p->getPe();
  }

  double efluid = 0.0, bfluid=0.0;
  double ch = 1.0, sh = 0.0;
  for(int ix=0;ix<maxX;ix++)
  for(int iy=0;iy<maxY;iy++)
  for(int iz=0;iz<maxZ;iz++) {
    double h=dZ*(iz-origZ);
    if(abs(h)>zmin) continue;
    if(optTauEta != 0) {
      ch=cosh(h); sh=sinh(h);
    }
    FluidElement *f = site(ix,iy,iz);
    //if(f->u(0)*HBARC < eFreezeOut) continue;
    //cout << " e= "<< f->u(0)*HBARC << " h= "<< h << " u0= "<< f->u(0) <<endl;
 
    efluid  += volH*(ch*f->u(0)+sh*f->u(3) );
    //pztot += volH*(sh*f->u(0)+ch*f->u(3) );
    //pxtot += f->u(1)*volH;
    //pytot += f->u(2)*volH;
    bfluid  += f->U(4)*vol;

  }

  if(opt==0) {
    particleDens->computeEnergyMomentum(plist,0,tau,0);
    double efrac= epart+efluid > 1e-20 ? efluid/(efluid+epart) : 0.0;
    out << setw(7) << tau
        << setw(15) << epart
        << setw(15) << efluid
        << setw(15) << efrac
        << setw(15) << particleDens->getEdens()
        << setw(15) << site(maxX/2,maxY/2,maxZ/2)->ed()*HBARC
        << setw(15) << particleDens->getRhob()
        << setw(15) << site(maxX/2,maxY/2,maxZ/2)->bd()
        <<endl;

  // after all fluid elements are converted into particles at once, 
  // we should check the fluid freezeout at early times.
  } else {
    int n=tau/dT;
    for(int i=1;i<n;i++) {
      particleDens->computeEnergyMomentum(plist,0,i*dT,1);
      out << setw(7) << i*dT
          << setw(15) << 0.0
          << setw(15) << 0.0
          << setw(15) << 0.0
          << setw(15) << particleDens->getEdens()
          << setw(15) << 0.0
          << setw(15) << particleDens->getRhob()
          << setw(15) << 0.0
          <<endl;
    }

  }

}

//***********************************************************************
void Fluid::hydroTotalEnergy(const int it,const double tau,
	double& etot,double& pxtot,double& pytot,double& pztot,double& btot)
{
  double etot0 = pTotal[0];
  double px0   = pTotal[1];
  double py0   = pTotal[2];
  double pz0   = pTotal[3];
  double btot0 = baryonTot;
  double entro0 = entroTot;

  etot = 0.0;
  pxtot = 0.0;
  pytot = 0.0;
  pztot = 0.0;
  btot  = 0.0;
  double entro=0.0;

  for(int ix=0;ix<maxX;ix++)
  for(int iy=0;iy<maxY;iy++)
  for(int iz=0;iz<maxZ;iz++) {
    FluidElement *f = site(ix,iy,iz);
    pxtot += f->u(1)*volH;
    pytot += f->u(2)*volH;
    btot  += f->U(4)*vol;
    double ch = 1.0, sh = 0.0;
    if(optTauEta != 0) {
      double h=dZ*(iz-origZ);
      ch=cosh(h); sh=sinh(h);
    }
    etot  += volH*(ch*f->u(0)+sh*f->u(3) );
    pztot += volH*(sh*f->u(0)+ch*f->u(3) );

    if(f->u(0) > 1e-12) {
      double v2=f->vsq();
      if(v2 < 1.0) {
        double gam=1.0/sqrt(1.0-v2);
        entro += eos->getEntropy(f->bd(),f->ed())*vol*gam;
      } else {
        cout << "hydroTotalenergy v>1? " << v2 <<endl;
      }
    }
  }

  double peerr = abs((etot-etot0)/etot0)*100;
  //double pzerr = abs((pztot-pz0)/max(1.0,pz0))*100;
  double berr  = abs((btot-btot0)/max(1.0,btot0))*100;
  double berr2  = abs((btot-btot0+bzero*vol+bout*vol)/max(1.0,btot0))*100;
  double serr  = abs((entro-entro0)/max(1.0,entro0))*100;

  cout << "======================================================" <<endl;
  cout << setprecision(6) << "btot0= "<< btot0 << " btot= "<< btot<< " err= "<< btot0-btot<<endl;

  cout << "*h it="<< it
      << " time= " << tau
      << " E= " << etot
      << " B= " << btot
      << " bout= " << bout
      << " B (%)= " << berr
      << " bzero= " << bzero*vol
      << " B+bzero (%)= " << berr2
      << " S (%)= " << serr
      <<endl;
  cout << "etot0= " << etot0
      << setw(12) << px0
      << setw(12) <<  py0
      << setw(12) <<  pz0
      << setw(12) <<  btot0
      <<endl;
  cout << "pzero= " << scientific << uzero[0]*volH
      << setw(15) << scientific << uzero[1]*volH
      << setw(15) << scientific << uzero[2]*volH
      << setw(15) << scientific << uzero[3]*volH
      <<endl;
  cout << "p_out=   " << uout[0]*volH
      << setw(15) << uout[1]*volH
      << setw(15) << uout[2]*volH
      << setw(15) << uout[3]*volH
      <<endl;

  cout << "p_total= " << etot
      << setw(15) << pxtot
      << setw(15) << pytot
      << setw(15) << pztot
      << " econ(%)= " << peerr
      <<endl;

  peerr = abs((etot+uzero[0]*volH-etot0)/etot0)*100;

  cout << "p_tot+pzero= " << px0-pxtot+uzero[1]*volH
       << setw(15) << py0-pytot+uzero[2]*volH
       << setw(15) << pz0-pztot+uzero[3]*volH
       << " econ(%)= "<<peerr
       <<endl;

  peerr = abs((etot+uout[0]*volH-etot0)/etot0)*100;
  cout << "p_tot+p_out= " << px0-pxtot+uout[1]*volH
       << setw(15) << py0-pytot+uout[2]*volH
       << setw(15) << pz0-pztot+uout[3]*volH
       << " econ(%)= "<<peerr
       <<endl;

  etot  += (uout[0]+uzero[0])*volH;
  pxtot += (uout[1]+uzero[1])*volH;
  pytot += (uout[2]+uzero[2])*volH;
  pztot += (uout[3]+uzero[3])*volH;
  peerr = abs((etot-etot0)/etot0)*100;

  cout << "p_tot+p_out+pzero= "
      << etot0-etot
      << setw(15) << px0-pxtot
      << setw(15) << py0-pytot
      << setw(15) << pz0-pztot
      << " econ(%)= "<<  peerr
      <<endl;

  if(abs(peerr)> 1e-5) {
    cout << "fluid energy does not converve"<<endl;
  }
  if((pxtot-px0)*(pxtot-px0)+(pytot-py0)*(pytot-py0)+(pztot-pz0)*(pztot-pz0)> 1e-5) {
    cout << "fluid momentum does not converve"<<endl;
    exit(1);
  }

}

//************************************************************************

int Fluid::checkout(double t)
{
  int check=0;

  double emax=-100.0;
  double pmax=-100.0;
  double tmax=-100.0;
  double bmax=-100.0;
  double r=0.0;
  int ir = 0;
  double aveb=0.0;
  double umax=0.0;

  for(int ix = N; ix < maxX-N; ix++)
  for(int iy = N; iy < maxY-N; iy++)
  for(int iz = N; iz < maxZ-N; iz++) {
    FluidElement *f = site(ix,iy,iz);
    umax=max(umax,f->u(0));
    tmax=max(tmax,eos->getT(f->bd(),f->ed()));
    pmax=max(pmax,f->pr());
    emax=max(emax,f->ed());
    bmax=max(bmax,f->bd());
    double b = f->bd();
    r += b;
    if(b > 0) ir++;
    if(f->ed() >= eFreezeOut) check++;
  }

  if(t <= passTime) check=1;

  if(r > 0.0) aveb = r / ir;
  double edens0=site(maxX/2,maxY/2,maxZ/2)->ed()*HBARC;

  if(jobMode > 0) {

    cout << " emax= "<< emax*HBARC
         << " bmax/rho0= "<< bmax/0.168
 	 << " pmax= "<< pmax*HBARC
 	 << " Tmax= "<< tmax*HBARC*1000
 	 << " check= "<< check
 	 << endl;

    cout << " e(0,0,0)= "<<  edens0
	 << " <rho_B/rho0>= "<< aveb/0.168
	 //<< " <T>= "<< avet
	 //<< " <mu>= " << avemu
	 << " passing time= " << passTime
	 << " umax = " << umax*HBARC
	 <<endl;

  }

  return check;
}



//***************************************************************************
void Fluid::outputhzx(int it)
{
  ofstream out("dens.dat");
  //out.open("dens.dat");

  cout << "#time = " << it*dT <<endl;
  int iy = origY;        
  for(int ix = 0; ix<maxX;ix++) {
  for(int iz = 0; iz<maxZ;iz++) {
    double x = dX*(ix-origX);
    double z = dZ*(iz-origZ);
    FluidElement* f = site(ix,iy,iz);
    out  << setw(6) << z
	 << setw(6) << x
	 << setw(14) << f->U(4)
	 << setw(14) << f->u(0)*HBARC
	 << setw(14) << f->bd()
	 << setw(14) << f->ed()
	 <<endl;
  }
    out << endl;
  }

  out.close();

}

double Fluid::gaussInter(Vec4& r, Vec4& u, double wid, int jx, int jy, int jz,
	int n_gaussp)
{
  double x[5],w[5];

  if(n_gaussp ==  1) {
    x[0]=0.0;
    w[0]=2.0;
  } else if(n_gaussp == 2) {
    x[0]= .57735026918962576451; // 1/sqrt(3)
    x[1]=-x[0];
    w[0]=1.0;
    w[1]=1.0;
  } else if(n_gaussp == 3) {
    x[0]=0;
    x[1]=-sqrt(3.0/5.0);
    x[2]=sqrt(3.0/5.0);
    w[0]=8.0/9.0;
    w[1]=5.0/9.0;
    w[2]=5.0/9.0;
  } else if(n_gaussp == 4) {
    x[0]= 0.33998104358485626480;
    x[1]= -x[0];
    x[2]=0.86113631159405257521;
    x[3]= -x[2];
    w[0]= 0.65214515486254614262;
    w[1]= w[0];
    w[2]= 0.34785484513745385737;
    w[3]= w[2];
  } else if(n_gaussp == 5) {
    x[0]=0;
    x[1]=0.53846931010568309103;
    x[2]=-x[1];
    x[3]=0.90617984593866399278;
    x[4]=-x[3];
    w[0]=0.56888888888888888888;
    w[1]=0.47862867049936646804;
    w[2]=w[1];
    w[3]=0.23692688505618908751;
    w[4]=w[3];
  } else {
    cout << "gaussInter wrong n_gaussp=" << n_gaussp <<endl;
    exit(1);
  }

  double eg=0.0;
  for(int i=0; i<n_gaussp;i++)
  for(int j=0; j<n_gaussp;j++)
  for(int k=0; k<n_gaussp;k++) {
   Vec4 r1(dX*(jx-origX)+x[i]*dX/2, 
	  dY*(jy-origY)+x[j]*dY/2,
	  dZ*(jz-origZ)+x[k]*dZ/2,
	  r[0]);
   double xtra = cmDistanceSquare(r-r1,u);
   eg += w[i]*w[j]*w[k]*exp( xtra/ wid );
  }
  return eg*vol/8;
}

// A, B: mass number of nucleus, b: impact parameter
// opt=1: Woods-Saxon  =2: Hard sphere nucleus
double Fluid::makeTwoNuclei(double ecm,double A, double B, double b, int opt)
{
  // Reset fluid variables.
  reset();
  updateNC();

  const int vacuum = 1;
  double em1=0.938;
  double em2=0.938;
  double s2=ecm*ecm;
  double p01=sqrt((s2-pow2(em1+em2))*(s2-pow2(em1-em2)))/(2*ecm);
  double e01=sqrt(em1*em1+p01*p01);
  //double ycm=0.5*log((e01+p01)/(e01-p01));
  double gam=e01/em1;
  double rad1=1.17*pow(A,0.333333) - 1.61*pow(A,-0.333333);
  double rad2=1.17*pow(B,0.333333) - 1.61*pow(B,-0.333333);
  double zini1=1.0+rad1/gam;
  double zini2=1.0+rad2/gam;
  if(opt !=  1) {
    zini1=0.01+rad1/gam;
    zini2=0.01+rad2/gam;
  }
  double t_pass=(zini1+zini2)/p01*e01;

  double rho0=0.17;
  //double rho0=0.15891;
  double  Be=0.016;

  /*
  int ivol=0;
  for(int iy = 0; iy<maxY;iy++)
  for(int ix = 0; ix<maxX;ix++)
  for(int iz = 0; iz<maxZ;iz++) {
    double x = dX*(ix-origX) - dX/2;
    double y = dY*(iy-origY) - dY/2;
    double z = dZ*(iz-origZ) - dZ/2;
    double pz = p01;
    if(z <= 0.0) {
      x -= b/2;
      z = (z + zini)*gam;
    } else {
      x += b/2;
      z = (z - zini)*gam
      pz=-p01;
    }

    if(opt ==  1) {
      nba=Ws(x,y,z,rad0,rho0)
    } else {
      nba=hard_sphere(x,y,z,rad0,rho0)
    }
    if(nba >= 1e-10) ivol++;
  }
  */

  double petot0=0.0;
  double pztot0=0.0;
  double bar=0.0;
  int ip=0;
  int im=0;
  //int N=2;
  for(int iy = N; iy<maxY-N;iy++)
  for(int ix = N; ix<maxX-N;ix++)
  for(int iz = N; iz<maxZ-N;iz++) {
    double x = dX*(ix-origX);
    double y = dY*(iy-origY);
    double z = dZ*(iz-origZ);
    double pz=p01;
    double rad0=rad1;
    if(z < 0.0) {
      x -= b/2;
      z = (z + zini1)*gam;
      ip++;
    } else {
      x += b/2;
      z = (z - zini2)*gam;
      pz=-p01;
      rad0=rad2;
      im++;
    }

    if(iz == origZ) pz=0.0;

    double nba=0.0;
    if(opt == 1) {
      nba=WS(x,y,z,rad0,rho0);
    } else {
      nba=hard_sphere(x,y,z,rad0,rho0);
    }

    if(vacuum == 1 && nba < 1e-10) continue;
    if(vacuum == 0) nba=max(0.3*rho0,nba);

    double pf=pow(1.5*M_PI*nba,1.0/3.0)*HBARC;
    double ene = (sqrt(em1*em1+pf*pf)-Be)*nba/HBARC;

    double pre=eos->getPressure(nba,ene);
    double sva=eos->getCs(nba,ene);

    FluidElement* f=site(ix,iy,iz);
    double t00 = (pre+ene)*gam*gam - pre;
    double uz =(pre+ene)*gam*gam*pz/e01;
    f->U(0,t00);
    f->U(1,0.0);
    f->U(2,0.0);
    f->U(3, uz);
    f->U(4,nba*gam);
    f->ed(ene);
    f->bd(nba);
    f->pr(pre);
    f->cs(sva);
    f->vx(0.0);
    f->vy(0.0);
    f->vz(pz/e01);
    if(t00*t00 - uz*uz <=0) {
	cout << " T00= " << t00 <<endl;
        cout << "ene= " <<  ene << " nba= " <<  nba << endl;
	exit(1);
    }

    pztot0 += uz;
    petot0 += t00*volH;
    bar += nba*gam;

  }

  Vec4 pc(0.0,0.0,pztot0,petot0);
  updatePTotal(pc);
  updateBTotal(bar);

  cout << "etot= " << petot0 << " pztot= " << pztot0 <<endl;
  double etot,pxtot,pytot,pztot,btot;
  hydroTotalEnergy(0,0.0, etot,pxtot,pytot,pztot,btot);

  return t_pass;

}

//***********************************************************************
double Fluid::WS(double x,double y,double z,double rad,double rho0)
{
  double rr=sqrt(x*x+y*y+z*z);
  return rho0/(1.0+exp((rr-rad)/0.53));
}

//***********************************************************************
double Fluid::hard_sphere(double x,double y,double z,double rad,double rho0)
{
  double rr=sqrt(x*x+y*y+z*z);
  return rr <= rad ? rho0 : 0.0;
}


} // namespace jam2


