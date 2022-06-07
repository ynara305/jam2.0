#include <jam2/meanfield/VectorPotential.h>
#include <jam2/hadrons/JamStdlib.h>

namespace jam2 {

using namespace std;

VectorPotential::VectorPotential(Pythia8::Settings* s) : PotentialType(s)
{
  optPotentialType=2;

  double C1=0.0, C2=0.0, mu1=1.0, mu2=1.0;
  t0=0.0;
  t5=1.0;
  t2f=0.0;

  withMomDep=1;
  if(eosType==1) {     // Skyrme Hard Nara2019 K=380 
    alpha = -0.124838736062805;
    beta  =  0.0707373967100532;
    gam   =  2.00261915156202;
    withMomDep=0;

  } else if(eosType==2) { // Skyrme Soft Nara2019 K=210
    alpha = -0.310514700691548;
    beta  =  0.256413361338797;
    gam   =  1.20292928384297;
    withMomDep=0;

  } else if(eosType==3) {  // MH2 Skyrme Hard Nara2019 K=380 
    // U( Pinf= 1.7 )= 0.06  Einf= 1.0036086114353737 Pinf= 1.7
    // U( 0.65 )= 0.0 Elab= 0.20320287416392357

    alpha = -0.00595607567906286;
    beta = 0.08157802422694654;
    gam = 1.718563134604995;
    C1 = -0.38616561269531274;
    mu1 = 2.02*HBARC;
    C2= 0.34267680471490275;
    mu2= 1.0*HBARC;

    /*
    alpha = 0.00877512904142221;
    beta = 0.07048243772253593;
    gam = 1.8907272542475995;
    C1 = -0.3907221709140509;
    C2= 0.3341868667966711;
    mu1 = 2.02*HBARC;
    mu2= 1.0*HBARC;
    */

  } else if(eosType==4) {  // MS2 Skyrme Hard Nara2019 K=210 
    // U( Pinf= 1.7 )= 0.06  Einf= 1.0036086114353737 Pinf= 1.7
    // U( 0.65 )= 0.0 Elab= 0.20320287416392357
    alpha = -0.3085250472810656;
    beta = 0.3804296791388347;
    gam = 1.1142288606816626;
    C1 = -0.38616561269531274;
    C2= 0.34267680471490275;
    mu1 = 2.02*HBARC;
    mu2= 1.0*HBARC;

  } else if(eosType==5) {  // MH3 Skyrme Hard Nara2019 K=380 

    alpha = 0.031029748068250693;
    beta = 0.04755488958383087;
    gam = 2.1354953844773155;
    C1 = -0.20192630016292473;
    C2= 0.05694208245922472;
    mu1 = 2.8*HBARC;
    mu2= 1.0*HBARC;

  // 2019/1/14
  } else if(eosType==6) {  // MS3 Skyrme Soft Nara2019 K=210 
    //alpha = -0.6285006244768873;
    //beta = 0.7070852697952635;
    //gam = 1.0496537865787445;
    //C1 = -0.20192630016292473;
    //C2= 0.05694208245922472;
    //mu1 = 2.8*HBARC;
    //mu2= 1.0*HBARC;

    alpha = -0.8290165881158311;
    beta = 0.907601225767913;
    gam = 1.0386838048982703;
    C1 = -0.20192630016293525;
    C2= 0.05694208245927023;
    mu1 = 2.8*HBARC;
    mu2= 1.0*HBARC;



  } else if(eosType==7) {  // MH4 Skyrme Hard Nara2019 K=380 
    alpha = 0.03953021450755792;
    beta = 0.040627668934799493;
    gam = 2.2928923888385695;
    C1 = -0.17086291974074935;
    C2 = 0.0;
    mu1 = 3.1466990715061636*HBARC;
    mu2 = 1.0*HBARC;

  } else if(eosType==8) {  // MS4 Skyrme Soft Nara2019 K=210 
    alpha = -0.22912691249560302;
    beta = 0.30928479593796043;
    gam = 1.1087616187247877;
    C1 = -0.17086291974074935;
    C2 = 0.0;
    mu1 = 3.1466990715061636*HBARC;
    mu2 = 1.0*HBARC;

  } else if(eosType==9) {  // MM4 Skyrme medium Nara2019 K=265
    // U( Pinf= 1.7 )= 0.06  Einf= 1.0036086114353737
    // U( 0.65 )= 0.0 Elab= 0.20320287416392357

    alpha = -0.0006561772779394997;
    beta = 0.08081406072029691;
    gam = 1.4918627502376933;
    C1 = -0.17086291974074935;
    C2= 0.0;
    mu1 = 3.1466990715061636*HBARC;
    mu2= 5.067730758770427;

  // optical potential is defined by sqrt{(m_N+S)^2+p^2} - sqrt{m_N^2 + p^2}
  } else if(eosType==11) {  // MH2 Skyrme Nara2021 K=380
    // U( Pinf= 1.7 )= 0.06  Einf= 1.0036086114353737 Pinf= 1.7
    // U( 0.65 )= 0.0 Elab= 0.20320287416392357
    alpha = -0.013119515535259911;
    beta = 0.08885442751779084;
    gam = 1.674140687543709;
    C1 = -0.3989613178044121;
    C2= 0.36728513480692454;
    mu1 = 2.02*HBARC;
    mu2= 1.0*HBARC;

  } else if(eosType==12) {  // MS2 Skyrme Nara2021 K=210
    // U( Pinf= 1.7 )= 0.06  Einf= 1.0036086114353737 Pinf= 1.7
    // U( 0.65 )= 0.0 Elab= 0.20320287416392357
    alpha = -0.5157475041588349;
    beta = 0.5906455475692723;
    gam = 1.0708570434690778;
    C1 = -0.3989613178044121;
    C2= 0.36728513480692454;
    mu1 = 2.02*HBARC;
    mu2= 1.0*HBARC;

  } else if(eosType==13) {  // MH4 Skyrme Nara2021 K=380
    // U( Pinf= 1.7 )= 0.06  Einf= 1.0036086114353737 Pinf= 1.7
    // U( 0.65 )= 0.0 Elab= 0.20320287416392357
    alpha = 0.03894524543017354;
    beta = 0.04170816920028206;
    gam = 2.2732990944321685;
    C1 = -0.16975960664857545;
    mu1 = 3.2272077792680443*HBARC;
  } else if(eosType==14) {  // MS4 Skyrme Nara2021 K=210
    // U( Pinf= 1.7 )= 0.06  Einf= 1.0036086114353737 Pinf= 1.7
    // U( 0.65 )= 0.0 Elab= 0.20320287416392357
    alpha = -0.23308866850927748;
    beta = 0.31374208313973306;
    gam = 1.1090643781089917;
    C1 = -0.16975960664857545;
    mu1 = 3.2272077792680443*HBARC;

  } else if(eosType==15) {  // MH2 Skyrme Nara2021 K=380 Hama2
    // U( Pinf= 1.6 )= 0.052  Einf= 0.916681643840797 Pinf= 1.6
    // U( 0.65 )= 0.0 Elab= 0.20320287416392357
    alpha = -0.00207556129465497;
    beta = 0.06962453317870083;
    gam = 1.8217624153926744;
    C1 = -0.34090163562062514;
    C2= 0.2800674442822501;
    mu1 = 2.02*HBARC;
    mu2= 1.0*HBARC;

  } else if(eosType==16) {  // MS2 Skyrme Nara2021 K=210 Hama2
    // U( Pinf= 1.6 )= 0.052  Einf= 0.916681643840797 Pinf= 1.6
    // U( 0.65 )= 0.0 Elab= 0.20320287416392357
    alpha = -0.9436482354010339;
    beta = 1.011197330566653;
    gam = 1.0379014197593113;
    C1 = -0.34090163562062514;
    C2= 0.2800674442822501;
    mu1 = 2.02*HBARC;
    mu2= 1.0*HBARC;

  } else if(eosType==17) {  // MH4 Skyrme Nara2021 K=380 Hama2
    //U( Pinf= 1.6 )= 0.052  Einf= 0.916681643840797 Pinf= 1.6
    // U( 0.65 )= 0.0 Elab= 0.20320287416392357
    alpha = 0.03117613882368128;
    beta = 0.03994951439261929;
    gam = 2.31156555163795;
    C1 = -0.1631814665874966;
    mu1 = 2.9741420289664107*HBARC;

  } else if(eosType==18) {  // MS4 Skyrme Nara2021 K=210 Hama2
    // U( Pinf= 1.6 )= 0.052  Einf= 0.916681643840797 Pinf= 1.6
    // U( 0.65 )= 0.0 Elab= 0.20320287416392357
    alpha = -0.21344609876974371;
    beta = 0.2845717519860443;
    gam = 1.1177471683656501;
    C1 = -0.1631814665874966;
    mu1 = 2.9741420289664107*HBARC;

  //  U=a*ρ/ρ0+b*(ρ/ρ0)^(4/3)+c*(ρ/ρ0)^(5/3)
  } else if(eosType==20) {  // Ohnishi 2022 M1 K=250MeV
    alpha = -202.8*1e-3;
    beta =  132.4*1e-3;
    gam = 4./3.;
    t5=18.3*1e-3/(1.+5./3.)/pow(rho0,5./3.);
    C1 = 0.0;
    mu1 = 3.2272077792680443*HBARC;
    withMomDep=0;

  //  U=a*ρ/ρ0+b*(ρ/ρ0)^(4/3)+c*(ρ/ρ0)^(5/3)+U_mom
  } else if(eosType==21) {  // Ohnishi 2022 MM1 K=250MeV
    alpha = -21.4*1e-3;
    beta =  92.1*1e-3;
    gam = 4./3.;
    t5=12.6*1e-3/(1.+5./3.)/pow(rho0,5./3.);
    C1 = -0.16975960664857545;
    mu1 = 3.2272077792680443*HBARC;

//Vector implementation (2022/05/19)	a(MeV)	b(MeV)	c(MeV)	C1(MeV)	C2(MeV)	μ1(/fm)	μ2(/fm)
//MH2 (K=380MeV)	-12.35778635	-3.515757907	91.60845624	-398.9613178	367.2851348	2.02	1
//MM2 (K=295MeV)	-97.35778635	194.8175754	-21.72487709	-398.9613178	367.2851348	2.02	1
//MS2 (K=210MeV)	-182.3577864	393.1509088	-135.0582104	-398.9613178	367.2851348	2.02	1
  //  U=a*ρ/ρ0+b*(ρ/ρ0)^(4/3)+c*(ρ/ρ0)^(5/3)+U_mom
  } else if(eosType==22) {  // Jinno 2022 MH2 K=380MeV
    alpha = -12.3577863512256*1e-3;
    beta  = -3.51575790693152*1e-3;
    gam = 4./3.;
    t5 = 91.6084562406881*1e-3/(1.+5./3.)/pow(rho0,5./3.);
    C1 = -398.961317804411*1e-3;
    C2 =  367.285134806923*1e-3;
    mu1 = 2.02*HBARC;
    mu2 = 1.00*HBARC;

  } else if(eosType==23) {  // Jinno 2022 MM2 K=295Mev
    alpha = -97.3577863512375*1e-3;
    beta  = 194.817575426425*1e-3;
    gam = 4./3.;
    t5 = -21.7248770926572*1e-3/(1.+5./3.)/pow(rho0,5./3.);
    C1 = -398.961317804411*1e-3;
    C2 =  367.285134806923*1e-3;
    mu1 = 2.02*HBARC;
    mu2 = 1.00*HBARC;

  } else if(eosType==24) {  // Jinno 2022 MS2 K=210Mev
    alpha = -182.357786351229*1e-3;
    beta  =  393.150908759742*1e-3;
    gam = 4./3.;
    t5 = -135.058210425982*1e-3/(1.+5./3.)/pow(rho0,5./3.);
    C1 = -398.961317804411*1e-3;
    C2 =  367.285134806923*1e-3;
    mu1 = 2.02*HBARC;
    mu2 = 1.00*HBARC;

  } else {
    cout << "VecgtorPotential:EoStype not implemented " << eosType << endl;
    exit(1);
  }

  pmu1= mu1*mu1;
  pmu2= mu2*mu2;

  if(transportModel==1) {
    t1 = 0.5*alpha/rho0;
    t3 = beta/(gam+1.0)/pow(rho0,gam);
    t3f = gam*t3;
    vex1=C1/(2*rho0);
    vex2=C2/(2*rho0);
  } else {
    t1 = alpha/rho0;
    t3 = beta/pow(rho0,gam);
    t3f = gam*t3;
    vex1=C1/rho0;
    vex2=C2/rho0;
  }


  if(firstCall) {
  cout << "# RQMDv mode eosType= "<< eosType  << " rho0= " << rho0
      << " alpha= " << alpha
      << " beta= "<< beta << " gamma= "<< gam
      << " t1= "<< t1 << " t3= "<< t3
      << " vex1= "<< vex1 << " pmu1= "<< pmu1
      << " vex2= "<< vex2 << " pmu2= "<< pmu2
      <<endl;
    firstCall=false;
  }
}

void VectorPotential::setPotentialParam(EventParticle* p)
{
  // Lambda potential U_Lambda(rho)= a*rho/rho0 + b*(rho/rho0)^4/3 + c*(rho/rho0)^5/3
 
  int islam = isHyperon(p);

  //            wid,alpha,beta,gam,  c1,  gam2,  cs, cv   mus muv
  double fp[10]={1.0, 1.0, 1.0, gam, 0.0, 5./3., 1.0,1.0, 1.0,1.0};
  int idp=1;

  if(islam==0) {
    if(t5 != 1.0) fp[4]=1.0; // three-range Skyrme potential
    p->setPotentialParam(fp);
    p->setPotentialId(idp);
    return;
  }
  idp=3;

  // Sigma potential
  if(optStrangeBaryonPot==1 && islam==2) {
    for(int i=0;i<10;i++) {
      fp[i]=facPotS[i];
    }
    if(t5 == 1.0) fp[4]=0.0; // two-range Skyrme potential
    fp[3]=gam;
    fp[5] = 5.0/3.0;
    p->setPotentialParam(fp);
    p->setPotentialId(idp);
    return;
  }


  // Hyperon potential is the same as nucleon potential multiplied by some factor.
  if(optLambdaPot==0) {
    for(int i=0;i<10;i++) {
      fp[i]=facPotL[i];
    }
    if(t5 == 1.0) fp[4]=0.0; // two-range Skyrme potential
    fp[3]=gam;
    fp[5] = 5.0/3.0;
    p->setPotentialParam(fp);
    p->setPotentialId(idp);
    return;
  }

  setLambdaPotential(fp);
  p->setPotentialParam(fp);
  p->setPotentialId(idp);

}

void VectorPotential::setSigmaPotential(double* fp)
{
  // Kohno-Kohno 2 body
  double gamS1=4.0/3.0;
  double gamS2=5.0/3.0;
  double aSigma= -17.38274141729;  // in MeV
  double bSigma= 118.398500622047; // in MeV
  double cSigma= -23.9710881588274; // in MeV
  double c1= -77.6961685779991; // in MeV
  double c2=0.0;
  //double mu1=3.23*HBARC; // in GeV
  double mu1 = 3.2272077792680443*HBARC;
  double mu2=1.0*HBARC; // in GeV

  fp[0]= facPotS[0];
  fp[3]= gamS1;
  fp[5]= gamS2;

  fp[1]= fp[0]*aSigma*1e-3/rho0/t1;
  fp[2]= fp[0]*bSigma*1e-3/pow(rho0,gamS1)/t3;
  fp[4]= fp[0]*cSigma*1e-3/pow(rho0,gamS2)/t5;

  fp[6]=0.0;
  fp[7]=0.0;
  if(vex1 !=0.0) fp[6] = c1*1e-3/(rho0)/vex1;
  if(vex2 !=0.0) fp[7] = c2*1e-3/(rho0)/vex2;

  fp[8] = mu1*mu1/pmu1;
  fp[9] = mu2*mu2/pmu2;

  // option for vector density dependent potential (not implemented in the vector potential!)
  if(optPotentialDensity>0) {
    fp[1] /= 2;
    fp[2] /= (gamS1+1.0);
    fp[4] /= (gamS2+1.0);
    fp[6] /=2;
    fp[7] /=2;
  }

}

void VectorPotential::setLambdaPotential(double* fp)
{
// Jinno 2022/1/5
//=1: 2+3BF upper
//=2: 2+3BF lower
//=3: 2BF   upper
//=4: 2BF   lower
//=5: GKW3 (2+3-body) average
//=6: GKW3 (2-body) average
//=7: GKW3 + Kohno (2+3-body) average mu1 fixed
//=8: GKW3 + Kohno (2-body) average mu1 fixed
//=9: Kohno + Kohno (2+3body) average mu1 fixed
//=10:Kohno + Kohno (2body) average mu1 fixed
//=11:GKW3 + Kohno (2+3-body) average
//=12:GKW3 + Kohno (2-body) average
//=13:Kohno + Kohno (2+3body) average
//=14:Kohno + Kohno (2body) average
//=15:GKW3 + Kohno (2+3-body) U(ρ0, p_inf=1.7 GeV) = 20 MeV, U(ρ0, p=0) = -30MeV
//=16:GKW3 + Kohno (2+3-body) U(ρ0, p_inf=1.7 GeV) = 40 MeV, U(ρ0, p=0) = -30MeV
//=17:GKW3 + Kohno (2+3-body) U(ρ0, p_inf=1.7 GeV) = 60 MeV, U(ρ0, p=0) = -30MeV
//=18:GKW3 + Kohno (2+3-body) U(ρ0, p_inf=1.7 GeV) = 10 MeV, U(ρ0, p=0) = -30MeV,  U(ρ0, p=2.5 /fm ) = 4.1MeV
//=19:GKW3 + Kohno (2+3-body) U(ρ0, p_inf=1.7 GeV) =  0 MeV, U(ρ0, p=0) = -30MeV,  U(ρ0, p=2.5 /fm ) = -23.4 MeV

  // Kohno-Kohno 2+3 body
  //optLambdaPot=9
  double gamL1=4.0/3.0;
  double gamL2=5.0/3.0;
  double aLambda= 12.3616606066095;  // in MeV
  double bLambda= 57.3358896065954; // in MeV
  double cLambda=  6.24830460415641;// in MeV
  double c1= -116.208988050382; // in MeV
  double c2=0.0;
  //double mu1=3.23*HBARC; // in GeV
  double mu1 = 3.2272077792680443*HBARC;
  double mu2=1.0*HBARC; // in GeV

  if(optLambdaPot<=4) {
    // Jinno 2022/1/5
    //=1: 2+3BF upper
    //=2: 2+3BF lower
    //=3: 2BF   upper
    //=4: 2BF   lower
    double J[4]={-29.15, -29.95, -33.44, -34.23};
    double L[4]={17.05, 7.64, 3.01, -6.66};
    double K[4]={556.59, 453.04, 364.43, 347.58};

    double JL=0.0,LL=0.0,KL=0.0;
    JL=J[optLambdaPot-1];
    LL=L[optLambdaPot-1];
    KL=K[optLambdaPot-1];

    aLambda= 10*JL-3*LL+KL/2;
    bLambda=-15*JL+5*LL-KL;
    cLambda=  6*JL-2*LL+KL/2;
    c1=0.0; c2=0.0;

  // GKW3 (2+3-body) average
  } else if(optLambdaPot==5) {

    aLambda=-80.1275;
    bLambda=0.16;
    cLambda=50.4175;
    c1=0.0; c2=0.0;

    // GKW2 (2body) average
  } else if(optLambdaPot==6) {
      aLambda=-154.8725;
      bLambda= 142.395;
      cLambda= -21.3575;
      c1=0.0; c2=0.0;

  // GKW3+Kohno (2+3body) average
  } else if(optLambdaPot==7) {
    //aLambda= 37.18423415;
    //bLambda= -10.76572129;
    //cLambda=  49.4892164;
    //c1= -116.1517136;

    //               lower              upper
    aLambda= 0.5*(51.5414237388091 + 66.1280377032306);
    bLambda= 0.5*(-21.722020209495 + -63.471866244167);
    cLambda= 0.5*(46.1264512880473 + 73.2896833582977);
    c1= -116.208988050382;
    c2= 0;

  // GKW2+Kohno (2body) average
  } else if(optLambdaPot==8) {
    //aLambda= -15.5759449;
    //bLambda= 89.42358328;
    //cLambda= -9.44069475;
    //c1= -104.2440349;

    //               lower             upper
    aLambda= 0.5*(13.7698072798742 - 6.68237330154632);
    bLambda= 0.5*(43.9660801950354 + 74.0106840000922);
    cLambda= 0.5*(6.70801810185758 - 2.88440512177866);
    c1= -104.298927774409;
    c2= 0;

  // Kohno+Kohno (3body) average
  } else if(optLambdaPot==9) {
    aLambda= 12.3616606066095;  // in MeV
    bLambda= 57.3358896065954; // in MeV
    cLambda=  6.24830460415641;// in MeV
    c1= -116.208988050382; // in MeV
    c2=0.0;

  // Kohno+Kohno (2body) average
  } else if(optLambdaPot==10) {
    aLambda= -5.18489127246517;
    bLambda= 72.7792470056277;
    cLambda= -3.15045015639529;
    c1= -104.298927774409;
    c2= 0;

  // GKW3+Kohno2 (2+3body) average
  } else if(optLambdaPot==11) {
    //               lower              upper
    aLambda= 0.5*(24.1707057443248 + 38.0042879575274);
    bLambda= 0.5*(-25.5959778717981 -66.0459593007329);
    cLambda= 0.5*(47.1067647690446 + 73.7231639847769);
    c1= -90.2339057147568;
    c2= 0;
    mu1 = 2.34783557120811*HBARC;

  // GKW2+Kohno2 (2body) average
  } else if(optLambdaPot==12) {
    //               lower             upper
    aLambda= 0.5*(-11.4394497501377 - 31.587841756003);
    bLambda= 0.5*( 39.606270943554  + 69.1243384432747);
    cLambda= 0.5*(7.89155855018705  - 1.47811694366828);
    c1= -79.739484950681;
    c2= 0;
    mu1 = 2.31923073877857*HBARC;

  // Kohno+Kohno2 (3body) average
  } else if(optLambdaPot==13) {
    aLambda= 24.1707057443248;
    bLambda=-25.5959778717981;
    cLambda= 47.1067647690446;
    c1= -90.2339057147568;
    c2= 0;
    mu1=2.34783557120811*HBARC;

  // Kohno+Kohno2 (2body) average
  } else if(optLambdaPot==14) {
    aLambda= -31.587841756003;
    bLambda= 69.1243384432747;
    cLambda= -1.47811694366828;
    c1= -79.739484950681;
    c2= 0;
    mu1=2.31923073877857*HBARC;

  // 2022/6/2 Jinnio parametrization
  // Weise+Kohno2 (3body) average
  // U(ρ0, p_inf=1.7 GeV) = 20 MeV ,  U(ρ0, p=0) = -30MeV
  } else if(optLambdaPot==15) {
    aLambda=  0.5*(6.33487624292818-7.36605124098376);  // -0.515587499;
    bLambda=  0.5*(43.9660801950354 + 74.0106840000922);//-31.20714418;
    cLambda=  0.5*(6.70801810185758 - 2.88440512177866);//59.57021749;
    c1=  -63.8049363598459;
    c2= 0;
    mu1=3.22720777926804*HBARC;

  // 2022/6/2 Jinnio parametrization
  // Weise+Kohno2 (3body) average
  //  U(ρ0, p_inf=1.7 GeV) = 40 MeV ,  U(ρ0, p=0) = -30MeV
  } else if(optLambdaPot==16) {
    aLambda=  0.5*(32.835815920346 + 19.0973970947713);
    bLambda=  0.5*(-53.7793093675166 -13.4923083796226);
    cLambda=  0.5*(71.9299650288679 + 45.3813828665485);
    c1=  -89.3269014736294;
    c2= 0;
    mu1=3.22720777926804*HBARC;

  // 2022/6/2 Jinnio parametrization
  // Weise+Kohno2 (3body) average
  //  U(ρ0, p_inf=1.7 GeV) = 60 MeV ,  U(ρ0, p=0) = -30MeV
  } else if(optLambdaPot==17) {
    aLambda=  0.5*(59.336750092793 + 45.5608398706086);
    bLambda=  0.5*(-56.240095705508 - 15.8888506249168);
    cLambda=  0.5*(71.0287980772872 + 44.4534632188805);
    c1=  -114.848861193033;
    c2= 0;
    mu1=3.22720777926804*HBARC;

  // 2022/6/3 Jinnio parametrization
  // Weise+Kohno2 (3body) average
  } else if(optLambdaPot==18) {
    aLambda=  0.5*(93.5756275550927 + 83.6253501899192);
    bLambda=  0.5*(-208.847972954215-175.067718159051);
    cLambda=  0.5*(125.799173793322 + 101.969196363332);
    c1= -260.083277587013;
    c2= 0;
    mu1=0.387779653268034*HBARC;

  // 2022/6/3 Jinnio parametrization
  // Weise+Kohno2 (3body) average
  } else if(optLambdaPot==19) {
    aLambda=  0.5*(7.9028792996856+5.75819167596766);
    bLambda=  0.5*(-49.502478157297 -89.6584713681133);
    cLambda=  0.5*(58.3276669138862 + 84.8225891490491);
    c1= -54.9929051628254;
    c2= 0;
    mu1=1.12431370250245*HBARC;

  } else if(optLambdaPot>100 && optLambdaPot<110){
    double J2[9]={-30.0,-30.0,-30.0, -30.0,-30.0,-30.0, -30.0,-30.0,-30.0};
    double L2[9]={-90.0,-90.0,-90.0,   0.0,  0.0,  0.0,  15.0, 15.0, 15.0};
    double K2[9]={  0.0,350.0,500.0,   0.0,350.0,500.0,   0.0,350.0,500.0};
    double JL2=0.0,LL2=0.0,KL2=0.0;
    JL2=J2[optLambdaPot-101];
    LL2=L2[optLambdaPot-101];
    KL2=K2[optLambdaPot-101];

    aLambda= 10*JL2-3*LL2+KL2/2;
    bLambda=-15*JL2+5*LL2-KL2;
    cLambda=  6*JL2-2*LL2+KL2/2;
    c1=0.0; c2=0.0;
  } else if(optLambdaPot>=110 && optLambdaPot<=115){
  //                20     21    22    23    24    25
    double J3[6]={-50.0,-30.0, 10.0, -50.0,-30.0, 10.0};
    double L3[6]={ 15.0, 15.0, 15.0,  -2.0, -2.0,  -20};
    double K3[6]={500.0,500.0,500.0, 350.0,350.0,350.0};
    double JL3=0.0,LL3=0.0,KL3=0.0;
    JL3=J3[optLambdaPot-110];
    LL3=L3[optLambdaPot-110];
    KL3=K3[optLambdaPot-110];

    aLambda= 10*JL3-3*LL3+KL3/2;
    bLambda=-15*JL3+5*LL3-KL3;
    cLambda=  6*JL3-2*LL3+KL3/2;
    c1=0.0; c2=0.0;
  } else {
    cout << "VectorPotential::setLambdaPotential wrong optLmambdaPot = "<< optLambdaPot<<endl;
    exit(1);
  }

  fp[0]= facPotL[0];
  fp[3]= gamL1;
  fp[5]= gamL2;

  fp[1]= fp[0]*aLambda*1e-3/rho0/t1;
  fp[2]= fp[0]*bLambda*1e-3/pow(rho0,gamL1)/t3;
  fp[4]= fp[0]*cLambda*1e-3/pow(rho0,gamL2)/t5;

  fp[6]=0.0;
  fp[7]=0.0;
  if(vex1 !=0.0) fp[6] = c1*1e-3/(rho0)/vex1;
  if(vex2 !=0.0) fp[7] = c2*1e-3/(rho0)/vex2;

  fp[8] = mu1*mu1/pmu1;
  fp[9] = mu2*mu2/pmu2;

  // option for vector density dependent potential (not implemented in the vector potential!)
  if(optPotentialDensity>0) {
    fp[1] /= 2;
    fp[2] /= (gamL1+1.0);
    fp[4] /= (gamL2+1.0);
    fp[6] /=2;
    fp[7] /=2;
  }

}


// vector potential --------------------------------------------------------------------------
void VectorPotential::Vmd(double psq1,double psq2,double vfac1,double vfac2,double pfi1,double pfi2,double pfj1,double pfj2)
{
  vmomsi=0.0;
  vmomsj=0.0;
  vmom4i = vfac1*vex1/(1.0-psq1/(pfi1*pmu1)) + vfac2*vex2/(1.0-psq1/(pfi2*pmu2));
  vmom4j = vfac1*vex1/(1.0-psq2/(pfj1*pmu1)) + vfac2*vex2/(1.0-psq2/(pfj2*pmu2));
}

void VectorPotential::dVdns(double rhij,double rhji)
{
  rhomij=rhij;
  rhomji=rhji;
  gam1=gama1;
  gam2=gamb1;

  //pFac1=pfac1*pfac2;
  //pFac2=pfac3*pfac4;
  //pFac3=pfac5*pfac6;

  pFac1=pfac1;
  pFac2=pfac3;
  pFac3=pfac5;

  Vec4 Ai = facV(JBi,rhoi,rhogi,rhog2i,v2);

  pFac1=pfac2;
  pFac2=pfac4;
  pFac3=pfac6;

  gam1=gama2;
  gam2=gamb2;

  Vec4 Aj = facV(JBj,rhoj,rhogj,rhog2j,v1);
  fskyi = -wG*(Ai * vk1)*bar1*bar2*rhomij;
  fskyj = -wG*(Aj * vk2)*bar1*bar2*rhomji;
}

void VectorPotential::dVdmd(double psq1,double psq2,double fs,double fv,
    double pmi1,double pmj1,double pmi2,double pmj2)
{
  vf1=1.0;
  vf2=1.0;
  if(optVectorPotential==1) {
    vf1 = vk1 * v2;
    vf2 = vk2 * v1;
  }
  fmomdi = -wG*devVmd2(psq2,fs,fv,pmi1,pmi2)*rhomij*vf1;
  fmomdj = -wG*devVmd2(psq1,fs,fv,pmj1,pmj2)*rhomji*vf2;
  fmomei =     devVme2(psq2,fs,fv,pmi1,pmi2)*rhomij*vf1;
  fmomej =     devVme2(psq1,fs,fv,pmj1,pmj2)*rhomji*vf2;

  //fmomdi = -wG*devVmd2(psq2,pmi1,pmi2)*rhomij*vf1*fv1*fv2;
  //fmomdj = -wG*devVmd2(psq1,pmj1,pmj2)*rhomji*vf2*fv1*fv2;
  //fmomei =     devVme2(psq2,pmi1,pmi2)*rhomij*vf1*fv1*fv2;
  //fmomej =     devVme2(psq1,pmj1,pmj2)*rhomji*vf2*fv1*fv2;
}

void VectorPotential::dVdp()
{
  // Vector Derivative of p^\mu/p0_j term in the vector potential.

  //pFac1=pfac1*pfac2;
  //pFac2=pfac3*pfac4;
  //pFac3=pfac5*pfac6;

  pFac1=pfac1;
  pFac2=pfac3;
  pFac3=pfac5;

  gam1=gama1;
  gam2=gamb1;
  Vec4 A1 = facV(JBi,rhoi,rhogi,rhog2i,vk1);
  pFac1=pfac2;
  pFac2=pfac4;
  pFac3=pfac6;

  gam1=gama2;
  gam2=gamb2;
  Vec4 A2 = facV(JBj,rhoj,rhogj,rhog2j,vk2);
  devV(A1,A2);
  forceri = devV1*rhomji*bar1*bar2;
  forcerj = devV2*rhomij*bar1*bar2;
  facsk= -(fskyi+fskyj)/wG;
}

  // Derivative of p_j / p_0j part in the momentum-dependent potential.
void VectorPotential::dVdpm()
{
  double facm1 = -fmomdj/(wG*vf2);
  double facm2 = -fmomdi/(wG*vf1);
  devV(vk1,vk2);
  forceri = facm1*devV1;
  forcerj = facm2*devV2;
  //forceri += facm1*(optP0dev*v1 - vk2)/p01;
  //forcerj += facm2*(optP0dev*v2 - vk1)/p02;
  facmom = -(fmomdi + fmomdj)/wG;
}

} // namespace jam2
