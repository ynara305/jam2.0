#include <jam2/fluid/FreezeOut.h>

namespace jam2 {

FreezeOut::FreezeOut(Pythia8::Settings* s,EoS* eos0)
{
  settings = s;
  eos=eos0;
  eFreezeOut = settings->parm("Hydro:ParticlizationEnergyDensity");
  eFreezeOut /= HBARC;  // 1/fm^4

  maxX=settings->mode("Hydro:nx");
  maxY=settings->mode("Hydro:ny");
  maxZ=settings->mode("Hydro:nz");
  dT=settings->parm("Hydro:dt");
  dX=settings->parm("Hydro:dx");
  dY=settings->parm("Hydro:dy");
  dZ=settings->parm("Hydro:dz");
  dV=dX*dY*dZ;

  FT=1; FX=1; FY=1; FZ=1;
  //FT=2;FX=2;FY=2;FZ=2;
  //FT=1;FX=2;FY=2;FZ=2;

  dfx[0]=dT*FT;
  dfx[1]=dX*FX;
  dfx[2]=dY*FY;
  dfx[3]=dZ*FZ;

  // minium temperature
  tFreezeOutCut = 0.003 / HBARC;

 int idm=4;
 cornelius.init(idm,eFreezeOut,dfx); 

  HyperCube = new double ***[2];
  dsite = new double ***[2];
  vxsite = new double ***[2];
  vysite = new double ***[2];
  vzsite = new double ***[2];
  for(int i=0;i<2;i++) {
    HyperCube[i] = new double **[2];
    dsite[i] = new double **[2];
    vxsite[i] = new double **[2];
    vysite[i] = new double **[2];
    vzsite[i] = new double **[2];
    for(int j=0;j<2;j++) {
      HyperCube[i][j] = new double *[2];
      dsite[i][j] = new double *[2];
      vxsite[i][j] = new double *[2];
      vysite[i][j] = new double *[2];
      vzsite[i][j] = new double *[2];
      for(int k=0;k<2;k++) {
        HyperCube[i][j][k] = new double [2];
        dsite[i][j][k] = new double [2];
        vxsite[i][j][k] = new double [2];
        vysite[i][j][k] = new double [2];
        vzsite[i][j][k] = new double [2];
      }
    }
  }

}

FreezeOut::~FreezeOut()
{
  for(int i=0;i<2;i++) {
    for(int j=0;j<2;j++) {
      for(int k=0;k<2;k++) {
        delete [] HyperCube[i][j][k];
        delete [] dsite[i][j][k];
        delete [] vxsite[i][j][k];
        delete [] vysite[i][j][k];
        delete [] vzsite[i][j][k];
      }
      delete [] HyperCube[i][j];
      delete [] dsite[i][j];
      delete [] vxsite[i][j];
      delete [] vysite[i][j];
      delete [] vzsite[i][j];
    }
    delete [] HyperCube[i];
    delete [] dsite[i];
    delete [] vxsite[i];
    delete [] vysite[i];
    delete [] vzsite[i];
  }
 
  delete [] HyperCube;
  delete [] dsite;
  delete [] vxsite;
  delete [] vysite;
  delete [] vzsite;

}

void FreezeOut::clear()
{
  freezeOutTime.clear();
  freezeOutR[0].clear();
  freezeOutR[1].clear();
  freezeOutR[2].clear();
  freezeOutV[0].clear();
  freezeOutV[1].clear();
  freezeOutV[2].clear();
  freezeOutV[3].clear();
  freezeOutT.clear();
  freezeOutMuB.clear();
  freezeOutMuS.clear();
  freezeOutDsgma[0].clear();
  freezeOutDsgma[1].clear();
  freezeOutDsgma[2].clear();
  freezeOutDsgma[3].clear();
  freezeOutBden.clear();
  freezeOutNum.clear();
  nFreezeout=0;
}

void FreezeOut::reset(vector<FluidElement*>& fluid)
{
  for(int iz=N;iz<=maxZ-N-FZ;iz++)
  for(int iy=N;iy<=maxY-N-FY;iy++)
  for(int ix=N;ix<=maxX-N-FX;ix++) {
    //FluidElement *f = fluid[site(ix,iy,iz)];
    fluid[site(ix,iy,iz)]->savePreviousLocalVal();
    /*
    fE[0][ix][iy][iz]  = f->ed();
    fD[0][ix][iy][iz]  = f->bd();
    fVx[0][ix][iy][iz] = f->vx();
    fVy[0][ix][iy][iz] = f->vy();
    fVz[0][ix][iy][iz] = f->vz();
    */
  }
  clear();
}

void FreezeOut::isoenergyFreezeout(vector<FluidElement*>& fluid,int it,double htime) 
{
  if (it % FT  != 0) return;

  /*
  for(int iz=N;iz<=maxZ-N-FZ;iz++)
  for(int iy=N;iy<=maxY-N-FY;iy++)
  for(int ix=N;ix<=maxX-N-FX;ix++) {
    FluidElement *f = fluid[ix + maxX* ( iy + maxY*iz)];
    fE[1][ix][iy][iz]  = f->ed();
    fD[1][ix][iy][iz]  = f->bd();
    fVx[1][ix][iy][iz] = f->vx();
    fVy[1][ix][iy][iz] = f->vy();
    fVz[1][ix][iy][iz] = f->vz();
  }
  */

  freezeout(fluid,htime);

  // Save the values for next freezeout.
  for(int iz=N;iz<=maxZ-N-FZ;iz++)
  for(int iy=N;iy<=maxY-N-FY;iy++)
  for(int ix=N;ix<=maxX-N-FX;ix++) {
    fluid[site(ix,iy,iz)]->savePreviousLocalVal();
    /*
    fE[0][ix][iy][iz]  = fE[1][ix][iy][iz];
    fD[0][ix][iy][iz]  = fD[1][ix][iy][iz];
    fVx[0][ix][iy][iz] = fVx[1][ix][iy][iz];
    fVy[0][ix][iy][iz] = fVy[1][ix][iy][iz];
    fVz[0][ix][iy][iz] = fVz[1][ix][iy][iz];
    */
  }

}

//--------------------------------------------------------
void FreezeOut::freezeout(vector<FluidElement*>& fluid,double htime) 
{
  for(int iz=N;iz<=maxZ-N-FZ;iz++)
  for(int iy=N;iy<=maxY-N-FY;iy++)
  for(int ix=N;ix<=maxX-N-FX;ix++)
  if ((ix%FX  == 0) && (iy%FY == 0) && (iz%FZ == 0 )) {
      int Ngp = 0;
      for(int l=0;l<2;l++) // z
      for(int k=0;k<2;k++) // y
      for(int j=0;j<2;j++) // x
      for(int i=0;i<2;i++) // t
      if( fluid[site(ix+j*FX,iy+k*FY,iz+l*FZ)]->ed(i) >= eFreezeOut) Ngp++;

      //if( fE[i][ix+j*FX][iy+k*FY][iz+l*FZ] >= eFreezeOut) Ngp++;

      if ((Ngp > 0) && (Ngp < 16)) {
        for(int l=0;l<2;l++) // z
        for(int k=0;k<2;k++) // y
        for(int j=0;j<2;j++) // x
        for(int i=0;i<2;i++)  {// t
          FluidElement *f = fluid[site(ix+j*FX,iy+k*FY,iz+l*FZ)];
          HyperCube[i][j][k][l] = f->ed(i);
          dsite[i][j][k][l]     = f->bd(i);
          vxsite[i][j][k][l]    = f->vx(i);
          vysite[i][j][k][l]    = f->vy(i);
          vzsite[i][j][k][l]    = f->vz(i);

	  /*
	  cout << " e= "<< f->ed(i)*HBARC
	       << " d= "<< f->bd(i)/0.168
	    << " T= "<< eos->getT(f->bd(i),f->ed(i))*HBARC<<endl;
	    */


	}

	cornelius.find_surface_4d(HyperCube);
	int Nsurf = cornelius.get_Nelements();

	for(int k=0;k<Nsurf;k++) {
	  double rt = cornelius.get_centroid_elem(k,0)/dfx[0];
	  double rx = cornelius.get_centroid_elem(k,1)/dfx[1];
	  double ry = cornelius.get_centroid_elem(k,2)/dfx[2];
	  double rz = cornelius.get_centroid_elem(k,3)/dfx[3];
          double tmid = htime -FT*dT + rt*dfx[0];

          //double xmid = dX*(ix-origX);
          //double ymid = dY*(iy-origY);
          //double zmid = dZ*(iz-origZ);

          double emid=interpolation4(HyperCube, rt, rx, ry, rz);
          double dmid=interpolation4(dsite, rt, rx, ry, rz);
          double vxmid=interpolation4(vxsite, rt, rx, ry, rz);
          double vymid=interpolation4(vysite, rt, rx, ry, rz);
          double vzmid=interpolation4(vzsite, rt, rx, ry, rz);
          double dst=cornelius.get_normal_elem(k,0);
          double dsx=cornelius.get_normal_elem(k,1);
          double dsy=cornelius.get_normal_elem(k,2);
          double dsz=cornelius.get_normal_elem(k,3);

          double vv=vxmid*vxmid + vymid*vymid + vzmid*vzmid;
          if(vv > 1.0) {
            cout << " v>1 ? vxmid= " << vxmid
	     << " vymid= " << vymid
     	     << " vzmid= " << vzmid
             << " v^2= " << vxmid*vxmid + vymid*vymid + vzmid*vzmid
	     <<endl;
           continue;
          }
	  double gam=1.0/sqrt(1.0 - vv);
          double vol=gam*(dst + dsx*vxmid + dsy*vymid + dsz*vzmid);
          double vopt=vol+sqrt(vol*vol - (dst*dst - dsx*dsx-dsy*dsy-dsz*dsz));
	  if(vopt<=0.0) continue;


          //eos->geteos2(emid,dmid,Tempmid,mumid,smumid,smid)
	  double Temp = eos->getT(dmid,emid);
	  double mub = eos->getMuB(dmid,emid);
	  double mus = eos->getMuS(dmid,emid);
          nFreezeout++;
          freezeOutTime.push_back(tmid);
          freezeOutR[0].push_back(ix);
          freezeOutR[1].push_back(iy);
          freezeOutR[2].push_back(iz);
          freezeOutV[0].push_back(gam);
          freezeOutV[1].push_back(vxmid);
          freezeOutV[2].push_back(vymid);
          freezeOutV[3].push_back(vzmid);
          freezeOutT.push_back(max(Temp,tFreezeOutCut));
          freezeOutMuB.push_back(mub);
          freezeOutMuS.push_back(mus);
          freezeOutDsgma[0].push_back(cornelius.get_normal_elem(k,0));
          freezeOutDsgma[1].push_back(cornelius.get_normal_elem(k,1));
          freezeOutDsgma[2].push_back(cornelius.get_normal_elem(k,2));
          freezeOutDsgma[3].push_back(cornelius.get_normal_elem(k,3));
          freezeOutBden.push_back(dmid);
          freezeOutNum.push_back(eos->getMult(dmid,emid));

	  /*
          FluidElement *f = fluid[site(ix,iy,iz)];
          double em=interpolation(HyperCube, rt, rx, ry, rz)*HBARC;
	  double eave=0.0;
	  for(int i=0; i<2;i++)
	  for(int j=0; j<2;j++)
	  for(int k=0; k<2;k++)
	  for(int l=0; l<2;l++) {
	    eave += HyperCube[i][j][k][l]*HBARC/16;
	    cout << " Hypercube= "<< HyperCube[i][j][k][l]*HBARC<<endl;
	  }

	  cout << " eave= "<< eave<<endl;
    cout <<" t= "<< htime << " em= "<< em << " e0= "<< f->ed(0)*HBARC << " e= "<< f->ed(1)*HBARC
      << " t0= "<< HBARC*eos->getT(f->bd(0),f->ed(0))
      << " t= "<< HBARC*eos->getT(f->bd(),f->ed())
      << " tmid= "<< Temp*HBARC
      << " nf= "<< nFreezeout
      << " emid= "<< emid*HBARC
      << " dmid= "<< dmid/0.168
      <<endl;
    cin.get();
    */

	  /*
	  cout << " freezeout T= "<< tmid*HBARC
	    << " htime= "<< htime
            << " rt= "<< rt << " dfx0= " <<dfx[0]
	     << " ft*dt= "<< FT*dT 
	      << " emid= "<< emid*HBARC
	      << " mu= "<< mub
	      << " vx= "<< vxmid
	      << " vy= "<< vymid
	      << " vz= "<< vzmid
	      <<endl;
	  */

	} //  end Nsurf

    }
  } // end loop over space.

}

//--------------------------------------------------------
double FreezeOut::interpolation(double ****qsite, double rt, double rx,
	double ry, double rz)
{
  double st = 1.0 - rt;
  double sx = 1.0 - rx;
  double sy = 1.0 - ry;
  double sz = 1.0 - rz;

  return  st*sx*sy*sz*qsite[0][0][0][0] 
        + st*sx*sy*rz*qsite[0][0][0][1]
        + st*sx*ry*sz*qsite[0][0][1][0] 
        + st*sx*ry*rz*qsite[0][0][1][1] 
        + st*rx*sy*sz*qsite[0][1][0][0] 
        + st*rx*sy*rz*qsite[0][1][0][1] 
        + st*rx*ry*sz*qsite[0][1][1][0] 
        + st*rx*ry*rz*qsite[0][1][1][1] 
        + rt*sx*sy*sz*qsite[1][0][0][0] 
        + rt*sx*sy*rz*qsite[1][0][0][1] 
        + rt*sx*ry*sz*qsite[1][0][1][0] 
        + rt*sx*ry*rz*qsite[1][0][1][1] 
        + rt*rx*sy*sz*qsite[1][1][0][0] 
        + rt*rx*sy*rz*qsite[1][1][0][1] 
        + rt*rx*ry*sz*qsite[1][1][1][0] 
        + rt*rx*ry*rz*qsite[1][1][1][1];
}

//--------------------------------------------------------
double FreezeOut::interpolation4(double ****qsite, double rt, double rx,
	double ry, double rz)
{
  double qmid=0.0;
  for(int j4=0;j4<2;j4++)
  for(int j3=0;j3<2;j3++)
  for(int j2=0;j2<2;j2++)
  for(int j1=0;j1<2;j1++) {
    qmid +=  (j1*rt+(1-j1)*(1-rt))*(j2*rx+(1-j2)*(1-rx))
             *(j3*ry+(1-j3)*(1-ry))*(j4*rz+(1-j4)*(1-rz))
             *qsite[j1][j2][j3][j4];
  }
  return qmid;

}

//--------------------------------------------------------
int FreezeOut::isochronousFreezeout(vector<FluidElement*>& fluid,double htime,int opt) 
{
  for(int iz=N;iz<=maxZ-N-FZ;iz++)
  for(int iy=N;iy<=maxY-N-FY;iy++)
  for(int ix=N;ix<=maxX-N-FX;ix++) {
    FluidElement *f = fluid[site(ix,iy,iz)];
    if(f->u(0) <= 1e-7) continue;
    double e=f->ed();
    //if(opt == 0 &&  e > eFreezeOut) continue;

    double d=f->bd();
    double temp = eos->getT(d,e);
    double mub = eos->getMuB(d,e);
    double mus = eos->getMuS(d,e);


    if(opt==0) {
      if(e > eFreezeOut) {
        f->savePreviousLocalVal();
	continue;
      }
      if(f->ed(0) < 1e-9 ) {
        f->savePreviousLocalVal();
	continue;
      }
      /*
      if(f->ed(0) < eFreezeOut ) {
        f->savePreviousLocalVal();
	continue;
      }
      */

      int jz1=max(0,iz-1), jz2=min(iz+2,maxZ);
      int jy1=max(0,iy-1), jy2=min(iy+2,maxY);
      int jx1=max(0,ix-1), jx2=min(ix+2,maxX);
      bool vacuum=false;

      /*
      for(int jz=jz1;jz<jz2;jz++)
      for(int jy=jy1;jy<jy2;jy++)
      for(int jx=jx1;jx<jx2;jx++) {
        if(fluid[site(jx,jy,jz)]->u(0) < 1e-9) vacuum=true;
      }
      */

      double dex=fluid[site(jx2,iy,iz)]->u(0)-fluid[site(jx1,iy,iz)]->u(0);
      double dey=fluid[site(ix,jy2,iz)]->u(0)-fluid[site(ix,jy1,iz)]->u(0);
      double dez=fluid[site(ix,iy,jz2)]->u(0)-fluid[site(ix,iy,jz1)]->u(0);
      if(f->u(1)*dex + f->u(2)*dey + f->u(3)*dez > 0) vacuum=true;

      if(vacuum==false) {
        f->savePreviousLocalVal();
	continue;
      }

    }

    cout <<" t= "<< htime << " e0= "<< f->ed(0)*HBARC << " e= "<< e*HBARC
      << " T= "<< temp*HBARC
      << " nf= "<< nFreezeout
      <<endl;

    freezeOutTime.push_back(htime);
    freezeOutR[0].push_back(ix);
    freezeOutR[1].push_back(iy);
    freezeOutR[2].push_back(iz);
    double vv= f->vx()*f->vx() + f->vy()*f->vy() + f->vz()*f->vz();
    freezeOutV[0].push_back(1.0/sqrt(1.0-vv));
    freezeOutV[1].push_back(f->vx());
    freezeOutV[2].push_back(f->vy());
    freezeOutV[3].push_back(f->vz());
    freezeOutT.push_back(max(temp,tFreezeOutCut));
    freezeOutMuB.push_back(mub);
    freezeOutMuS.push_back(mus);
    freezeOutDsgma[0].push_back(dV);
    freezeOutDsgma[1].push_back(0.0);
    freezeOutDsgma[2].push_back(0.0);
    freezeOutDsgma[3].push_back(0.0);
    freezeOutBden.push_back(d);
    freezeOutNum.push_back(eos->getMult(d,e));
    nFreezeout++;
    f->savePreviousLocalVal();

    /*
    double v1= f->vx();
    double v2= f->vy();
    double v3= f->vz();
    double p=eos->getPressure(d,e);
    double g2=1./(1-vv);
    cout << "E= "<< f->u(0) << " T00= "<< (e+p)*g2-p <<endl;
    cout << "Mx= "<< f->u(1) << " T0x= "<< (e+p)*g2*v1 << endl;
    cout << "My= "<< f->u(2) << " T0y= "<< (e+p)*g2*v2 << endl;
    cout << "Mz= "<< f->u(3) << " T0z= "<< (e+p)*g2*v3 << endl;
    */

  } // end loop over fluid element

    //cin.get();

  return nFreezeout;
}

} // end namespace jam2
