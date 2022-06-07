#include <jam2/hadrons/DecayWidth.h>
//#include <jam2/xsection/XsecTable.h>
#include <jam2/hadrons/JamStdlib.h>
#include <jam2/hadrons/GaussPoints.h>
#include <Pythia8/PythiaStdlib.h>

using namespace std;

namespace jam2 {

using Pythia8::pow3;
using Pythia8::pow4;
using Pythia8::pow6;
using Pythia8::pow8;
using Pythia8::Vec4;
using Pythia8::ParticleDataEntry;

DecayWidth::DecayWidth(Pythia8::ParticleData* pd)
{
  particleData=pd;

  //optDeltaWidth=1; // Frankfurt mstc(63)
  optDeltaWidth=3; // GiBUU mstc(63)
  //optDeltaWidth=4; // Randrap mstc(63)

  optWidth=3; //  BlattWeisskopf mstc(65) 
  //optWidth=1; //mstc(65) Frankfurt

  minKinE=0.003; // parj(64)

  useTable=true;
}

double DecayWidth::getPartialWidth(Pythia8::ParticleDataEntry* pd,double emcm, int kf1, int kf2)
{
//...Purpose: to calculate momentum dependent partial decay width
//...for resonances.
//======================================================================*
//     kc     : Compressed particle code                                *
//     kf1 kf2: ingoing particle KF code.                               *
//     emcm   :  Mass of the particle (GeV)                             *
//     pwid   :  partial decay width  (GeV) (output)                    *
//     totwid :  total width                                            *
//     itag   : branch which is identical to the kf1 and kf2.           *
//======================================================================*

  static double delta=0.1;
  itag=-1;
  double widr=pd->mWidth();
  if(widr < 1e-15) return widr;
  int id0=pd->id();

  int nd=pd->sizeChannels();
  nChannel=nd;
  for(int ibra=0;ibra<nd;ibra++) {
    pWidth[ibra]=0.0;
    Pythia8::DecayChannel d=pd->channel(ibra);
    int kfd1=d.product(0);
    int kfd2=d.product(1);
    if(kfd1 == 0 || kfd2 ==0) continue;
    if(id0 < 0 && particleData->hasAnti(kfd1)) kfd1=-kfd1;
    if(id0 < 0 && particleData->hasAnti(kfd2)) kfd2=-kfd2;
    if( (kfd1 == kf2 && kfd2==kf1) || (kfd1 == kf1 && kfd2==kf2)) {
      //int  ldec = d.meMode();
      if(d.onMode()==0) continue;
      if(d.bRatio() <= 1.e-10) continue;
      // Check threshold. 
      double totmas = 0.0;
      int nprod=d.multiplicity();
      for (int i = 0; i < nprod; ++i) {
      int idNow = d.product(i);
      double mr = particleData->mWidth(idNow) > 1e-7 ?
                  particleData->mMin(idNow) : particleData->m0(idNow);
        totmas += mr;
      }
      if(emcm < totmas + eKinMin) continue;
      itag=ibra;

      // constant width.
      if(optWidth==0) {
        pWidth[ibra]=d.bRatio()*widr;
        return pWidth[ibra];
      }
    }
  }
  if(itag==-1) return 0.0;

  double emres=pd->m0();
  int ibra = itag;
  Pythia8::DecayChannel d=pd->channel(ibra);
  int kfd1=d.product(0);
  int kfd2=d.product(1);

  double m1=particleData->m0(kfd1);
  double m2=particleData->m0(kfd2);
  double mth = m1 + m2 + minKinE;
  if(emcm < mth || emres < mth) return 0.0;

  double widp=d.bRatio();
  // Special for Delta(1232).
  if(isDelta(id0)) {
    pWidth[ibra]=widp*getDeltaWidth(optDeltaWidth,emcm);
    return pWidth[ibra];
  }

  int ldec;
  if(d.meMode() >=3 && d.meMode() <=7) ldec = d.meMode() - 3;
  else if(d.meMode() ==2) ldec=3;
  else ldec=1;
  double prres=PCM(emres,m1,m2);
  double pr=PCM(emcm,m1,m2);
  double form;
  // Frankfurt
  if(optWidth == 1) {
    form=1.2/(1.0+0.2*pow(pr/prres,2*ldec));
    pWidth[ibra]=widp*pow(pr/prres,2*ldec+1)*(emres/emcm)*form*widr;
  } else if (optWidth==2) {
    //form=pow((prres*prres+delta)/(pr*pr+delta),ldec+1);
    form=pow((prres*prres+delta)/(pr*pr+delta),ldec);
    pWidth[ibra]=widp*pow(pr/prres,2*ldec+1)*(emres/emcm)*form*widr;
  } else {
    pWidth[ibra]=widp*widr*formBlattWeisskopf(emres,emcm,prres,pr,ldec,kfd1,kfd2);
  }
  return pWidth[ibra];
}

double DecayWidth::getTotalWidth(Pythia8::ParticleDataEntry* pd,double emcm, int kf1, int kf2)
{
//...Purpose: to calculate momentum dependent partial decay width
//...for resonances.
//======================================================================*
//     kc     : Compressed particle code                                *
//     kf1 kf2: ingoing particle KF code.                               *
//     emcm   :  Mass of the particle (GeV)                             *
//     pwid   :  partial decay width  (GeV) (output)                    *
//     totwid :  total width                                            *
//     itag   : branch which is identical to the kf1 and kf2.           *
//======================================================================*

  //static double delta=0.09;
  //static double delta=HBARC*HBARC;
  static double delta=0.1;

  double totwid=0.0;
  itag=-1;

  //if(emcm < pd->mMin() + eKinMin) return 0.0;

  double widr=pd->mWidth();
  if(widr < 1e-15) return widr;

  int id0=pd->id();
  // Special for Delta(1232).
  if(isDelta(id0)) return deltaWidth(pd,id0,emcm,kf1,kf2);

  double emres=pd->m0();
  int nd=pd->sizeChannels();
  nChannel=nd;

  for(int ibra=0;ibra<nd;ibra++) {
    pWidth[ibra]=0.0;
    Pythia8::DecayChannel d=pd->channel(ibra);
    int isd=d.onMode();
    if(isd==0) continue;
    int ldec;
    if(d.meMode() >=3 && d.meMode() <=7) ldec = d.meMode() - 3;
    else if(d.meMode() ==2) ldec=3;
    else ldec=1;
    double widp=d.bRatio();
    if(widp <= 1.e-10) continue;

    // Check threshold. 
    double totmas = 0.0;
    int nprod=d.multiplicity();
    for (int i = 0; i < nprod; ++i) {
      int idNow = d.product(i);
      double mr = particleData->mWidth(idNow) > 1e-7 ?
                  particleData->mMin(idNow) : particleData->m0(idNow);
      totmas += mr;
    }
    if(emcm < totmas + eKinMin) continue;

    // constant width or more than two products.
    if(optWidth ==0  || nprod != 2) {
      pWidth[ibra]=widp*widr;
      totwid += pWidth[ibra];
      continue;
    }
      
    //.....Lepton, etc.
    //if(abs(kfd1) <= 100 || abs(kfd2) <= 100) continue;

    int kfd1=d.product(0);
    int kfd2=d.product(1);
    if(kfd1 == 0 || kfd2 ==0) continue;
    if(id0 < 0 && particleData->hasAnti(kfd1)) kfd1=-kfd1;
    if(id0 < 0 && particleData->hasAnti(kfd2)) kfd2=-kfd2;
    if( (kfd1 == kf2 && kfd2==kf1) ||
        (kfd1 == kf1 && kfd2==kf2)) itag=ibra;

    // Mom. dep. width
    double m1=particleData->m0(kfd1);
    double m2=particleData->m0(kfd2);
    double mth = m1 + m2 + minKinE;
    if(emcm >= mth && emres >= mth) {
      double prres=PCM(emres,m1,m2);
      double pr=PCM(emcm,m1,m2);
      double form;
      // Frankfurt
      if(optWidth == 1) {
        form=1.2/(1.0+0.2*pow(pr/prres,2*ldec));
        pWidth[ibra]=widp*pow(pr/prres,2*ldec+1)*(emres/emcm)*form*widr;
      } else if (optWidth==2) {
        //form=pow((prres*prres+delta)/(pr*pr+delta),ldec+1);
        form=pow((prres*prres+delta)/(pr*pr+delta),ldec);
        pWidth[ibra]=widp*pow(pr/prres,2*ldec+1)*(emres/emcm)*form*widr;
      } else {
        pWidth[ibra]=widp*widr*formBlattWeisskopf(emres,emcm,prres,pr,ldec,kfd1,kfd2);
      }
    } else {
      //pWidth[ibra]=widp*widr;
      pWidth[ibra]=0.0;
    }
    totwid += pWidth[ibra];

  //cout << id0 << " kf1= "<< kfd1 << " kf2= "<<kfd2 << " itag= "<< itag <<endl;
  //cout << " g1= "<< g1 << " g2= "<< g2 << " l= "<< ldec <<endl;

  } // loop over products.

  return totwid;
}

double DecayWidth::formBlattWeisskopf(double m0, double m,
	double pf0, double pf, int l, int kfd1, int kfd2)
{
  //double mr = particleData->mWidth(idNow) > 1e-7 ?
  //                particleData->mMin(idNow) : particleData->m0(idNow);

  bool optInt=true;
  optInt=false;
  bool ist1 = particleData->mWidth(kfd1) < 1e-7 ?  true : false;
  bool ist2 = particleData->mWidth(kfd2) < 1e-7 ?  true : false;
  double R=1.0/HBARC;

  //cout << " kfd1= "<< kfd1 << " ist1= " << ist1 << " wid= "<< particleData->mWidth(kfd1)<<endl;
  //cout << " kfd2= "<< kfd2 << " ist2= " << ist2 << " wid= "<< particleData->mWidth(kfd2)<<endl;

  if(!optInt || (ist1 && ist2)) {
    double bl  = BlattWeisskopf(pf*R,l);
    double bl0 = BlattWeisskopf(pf0*R,l);
    //cout << " stable= "<<  pf/pf0 * m0/m * (bl*bl)/(bl0*bl0);
    //cin.get();
    return pf/pf0 * m0/m * (bl*bl)/(bl0*bl0);
  }
  if(ist1 && !ist2) {
    double bl0 = BlattWeisskopfInt(m0,m0,kfd2,kfd1,l,R);
    double bl  = BlattWeisskopfInt(m,m0, kfd2,kfd1,l,R);
    //cout << " p1 is stable "<<  bl/bl0 << " b= "<< bl << " b0= "<<bl0 <<endl;
    //cin.get();
    return bl/bl0;
  } else if(!ist1 && ist2) {
    double bl0 = BlattWeisskopfInt(m0,m0,kfd1,kfd2,l,R);
    double bl  = BlattWeisskopfInt(m, m0,kfd1,kfd2,l,R);
    //cout << " p2 is stable "<<  bl/bl0 << " b= "<< bl << " b0= "<< bl0 <<endl;
    //cin.get();
    return bl/bl0;
  } else {
    cout << "DecayWidth::formBlattWeisskopf both are unstable particles?"
	<<endl;
    cout << " kf1= "<< kfd1 << " kfd2= "<<kfd2<<endl;
    exit(1);
  }

}

double DecayWidth::BlattWeisskopfInt(double m, double m0,int kfd1,int kfd2, int l,double R)
{
  ParticleDataEntry* pd=particleData->findParticle(kfd1);
  double widr=pd->mWidth();
  double mr=pd->m0();
  double mmin= pd->mMin();
  int ich= (isN1440(kfd1)||isDelta(kfd1)||isLambda1405(kfd1)||isLambda1520(kfd1)) ? 1:0;
  Pythia8::DecayChannel d=pd->channel(ich);
  int kfd3=d.product(0);
  int kfd4=d.product(1);
  double m3=particleData->m0(kfd3);
  double m4=particleData->m0(kfd4);
  double prres=PCM(mr,m3,m4);

    if(mmin < m3+m4) {
	cout << "kf= "<< kfd1;
	cout << " m1= "<< mmin << " m3= "<< m3 << " m4= " << m4<<endl;
	exit(1);
    }

  //cout << " kfd3= "<< kfd3 << " kfd4= "<< kfd4 <<endl;
  //cout << " m3= "<< m3 << " m4= "<< m4 <<endl;
  //cin.get();

  //double widr = particleData->mWidth(kfd1);
  //double mr= particleData->m0(kfd1);
  //double mmin= particleData->mMin(kfd1);

  double m2 = particleData->m0(kfd2);
  double mmax=m - m2;
  double s0=mmin + m2;
  //ParticleDataEntry* pd2=jamParticleData->find(kfd2);
  //double lam=pow4(0.8);
  double lam=0.4096; // 0.8^4  e.g. pi + rho
  if(pd->isBaryon()) lam=16.0; // =2^4 unstable baryon pi+Delta/pi+N(1440)
  else if(particleData->isBaryon(kfd2)) lam=6.5536; // =1.6^4 unstable meson rho+N

  /*
  if(isDelta(kfd1)) { // R -> D + pi
      lam=16;  // 2^4
        //mr=emdelt;
        //m3=emnuc;
        //m4=empion;
        //widr=0.117;
        //emin=emnuc+empion;
  } else if(isN1440(kfd1)) { // D* -> N*(1440) + pi
      lam=16;  // 2^4
        //mr=1.44;
        //m3=emnuc;
        //m4=empion;
        //widr=0.35;
        //emin=emnuc+empion;
  } else if (isRho(kfd1)) { // R -> N + rho
      lam=6.5536; // 1.6^4
        //mr=0.77;
        //m3=empion;
        //m4=empion;
        //emin=2*empion;
        //widr=0.15;
  } else if(isSigmaMeson(kfd1)) { // N* -> N + sigma meson
      lam=6.5536; // 1.6^4
        //mr=0.55;
        //m3=empion;
        //m4=empion;
        //emin=2*empion;
        //widr=0.5;
  } else {
      cout << " DecayWidth::BlattWeisskopfInt kfd1= "<< kfd1 <<endl;
      exit(1);
  }
  */

  // integral over mass
  double ds=0.0;
  for(int i=0;i<38;i++) {
    double m1=mmin+xg38f[i]*(mmax-mmin);
    double dw=wg38f[i]*(mmax-mmin);
    double pr=PCM(m1,m3,m4);
    double form=(prres*prres+R*R)/(pr*pr+R*R);
    double gam=widr*pow(pr/prres,3)*(mr/m1)*form;
    double bw=2.0/M_PI*m1*m1*gam/(pow2(m1*m1-mr*mr)+pow2(m1*gam));
    double cut=(lam+0.25*pow2(s0-m0*m0))/(lam+pow2(m*m-0.5*(s0+m0*m0)));
    double pf=PCM(m,m1,m2);
    double bl = BlattWeisskopf(pf*R,l);
    ds += cut*cut*bw*pf/m*bl*bl*dw;

    /*
    cout << "cut= "<< cut
	<< " bw= "<< bw
	<< " m= " << m
	<< " bl= "<< bl
	<< " l= "<< l
	<<endl;
	*/
  }
  return ds;

}


double DecayWidth::BlattWeisskopf(double x,int l)
{
  switch (l) {
    case 0: return 1;
    case 1: return x/sqrt(1+x*x);
    case 2: return x*x/sqrt(9+3*x*x+pow4(x));
    case 3: return pow3(x)/sqrt(225+45*x*x+6*pow4(x)+pow6(x));
    case 4: return pow4(x)/sqrt(11025+1575*x*x+135*pow4(x)+10*pow6(x)+pow8(x));
    default:
      cout << "DecayWidth::BlatWeisskopf error l= "<< l<<endl;
      exit(1);
  }
}

double DecayWidth::deltaWidth(Pythia8::ParticleDataEntry* pd, int id0, double emcm, int kf1, int kf2)
{
  double totwd=getDeltaWidth(optDeltaWidth,emcm);
  //double totwd=getDeltaWidth(1,emcm);
  if(totwd < 0.0001) {
    cout << " DecayWidth::getTotalWidth delta width? id0= " << id0
	<< " m= " << emcm
	<< " w= " << totwd
	<< endl;
    exit(1);
  }

  double totwid=0.0;
  int nd=pd->sizeChannels();
  nChannel=nd;
  for(int ibra=0;ibra<nd;ibra++) {
    pWidth[ibra]=0.0;
    Pythia8::DecayChannel d=pd->channel(ibra);
    int kfd1=d.product(0);
    int kfd2=d.product(1);
    if(id0 < 0 && particleData->hasAnti(kfd1)) kfd1=-kfd1;
    if(id0 < 0 && particleData->hasAnti(kfd2)) kfd2=-kfd2;
    if( (kfd1 == kf2 && kfd2==kf1) ||
	(kfd1 == kf1 && kfd2==kf2)) itag=ibra;
    double widp=d.bRatio();
    pWidth[ibra]=totwd*widp;
    totwid += pWidth[ibra];

    //cout << emcm << " kf= " << id0 << " kf1= "<< kfd1 << " kf2= "<< kfd2<< endl;
    //cout << ibra << " pwid= "<< pWidth[ibra] << " tot= "<< totwid << " itag= "<< itag<<endl;

  }

  /*
  double w=0.0;
  for(int ibra=0;ibra<nd;ibra++) {
      w+= pWidth[ibra];
      cout << scientific << " ** pwid= "<< ibra << " " <<  pWidth[ibra] <<endl;
  }
      cout << scientific << " delta width strange " << abs(w-totwid)
	  << " totwid = "<< totwid
	  << " w= " << w
	  <<endl;
	  */

  return totwid;
}

//***********************************************************************

double DecayWidth::getDeltaWidth(int iwidth, double emd)
{
//....Purpose: to calculate momentum dependent total decay width
//...  of delta(1232).
//
//optDeltaWidth: no decay(frozen delta)
//...         1: Frankfurt
//...         2: Giessen
//...         3: GiBUU
//...         4: Randrup
//...           Ref: Randrup, NP A314 (1979) 429.
//...                Rittenberg, REV.MOD.PHYS. 43 (1971) S1.
//...         5: Kitazoe
//...         6: Barz/Iwe  NP A453(1986)728
//...   other: constant width 0.12GeV

    static const double emnuc=0.9383,empion=0.138,ekinmi=0.0001;
    static const double emin=emnuc+empion+ekinmi;
    static const double mpion2=empion*empion;
    static const double emdelt=1.232,widdlt=0.12;
    static const double bet2=0.3*0.3;
    static const double bet3=HBARC*HBARC, gamr3=0.117,qqr=0.227894932;
    static const double qqr2=0.051936,gamr=0.11;
    static const double pscal1= 0.238, pscal2= 0.318, p0ref=0.227;

    if(emd < emin) {
         cout << scientific << setprecision(8) << 
	     "(jamdlwid:)invalid d(1232) mass " << emd 
	     << " min= " << emin << endl;
	 return 1e-35;
    }

    double pp=PCM(emd,emnuc,empion);
    double pp2=pp*pp;
    switch (iwidth) {
	case 1:
           return 0.12*emdelt/emd*pow(sqrt(pp2/qqr2),3)*1.2/(1+0.2*pp2/qqr2);
	case 2:
	   {
            double form= (1.0+qqr2/bet2)/(1.0+pp2/bet2);
            return pow(sqrt(pp2/qqr2),3)*emdelt/emd*gamr*form*form;
	   }
	case 3:
            return gamr3*pow3(pp/qqr)*emdelt/emd*(bet3+qqr2)/(bet3+pp2);
	case 4:
            return widdlt*(pp*pp*pp/(1.0+pow2(pp/pscal1)+pow4(pp/pscal2)))
          /(pow3(p0ref)/(1.0+pow2(p0ref/pscal1)+pow4(p0ref/pscal2)));
	case 5:
	    return 0.47/(1.0+0.6*pp2/mpion2*pp2)/mpion2*pp;
	case 6:
	    return 29.0*pow3(pp)/(1.0+40.0*pp2);
	case 0:
	    return 1e+35;
	default:
	    return widdlt;
    }


}

}
