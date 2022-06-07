#include <jam2/meanfield/NonRelPotential.h>
#include <jam2/hadrons/JamStdlib.h>

namespace jam2 {

using namespace std;

NonRelPotential::NonRelPotential(Pythia8::Settings* s) : PotentialType(s)
{
  withMomDep=1;

  double C1=0.0, C2=0.0;
  double mu1 = 2.02*HBARC;
  double mu2 = 1.0*HBARC;

  t0=0.0;
  t5=1.0;
  t2f=0.0;

  if(eosType==1) {     // Skyrme Hard  K=380 MeV
    alpha = -0.12483873606280496;
    beta = 0.07073739671005325;
    gam = 2.0026191515620178;
    //alpha = -124.0e-3;
    //beta = 70.5e-3;
    //gam = 2.0;
    withMomDep=0;
  } else if(eosType==2) { // Skyrme Soft K=210 MeV
    alpha = -0.3105147006915484;
    beta = 0.2564133613387967;
    gam = 1.2029292838429695;
    //alpha = -356.0e-3;
    //beta = 303.0e-3;
    //gam = 1.16667;
    withMomDep=0;

//    K= 0.2 alpha= -0.3807617551886489 beta= 0.3266604158358972 gam= 1.1558887033889078

  } else if(eosType==3) { // Skyrme + mom.dep Isse2005 K=448
    alpha = -33.e-3;
    beta  = 110.e-3;
    gam = 1.6666;  // 5.0/3.0
    mu1 = 2.35*HBARC;
    mu2 = 0.4*HBARC;
    C1  =  -0.277;
    C2 =    0.663;

  } else if(eosType==4) { // Skyrme + mom.dep Medium Isse2005 K=?
    alpha = -116.e-3;
    beta  = 193.e-3;
    gam = 1.33333;
    mu1 = 2.35*HBARC;
    mu2 = 0.4*HBARC;
    C1 =  -0.277;
    C2 =   0.663;
  } else if(eosType==5) { // Skyrme + mom.dep Medium Isse2005 K=314
   alpha = -268.e-3;
   beta = 345.e-3;
   gam = 1.16667; // 7/6
   mu1 = 2.35*HBARC;
   mu2 = 0.4*HBARC;
   C1 =  -0.277;
   C2 =   0.663;

  } else if(eosType==6) { // Skyrme + mom.dep Ohnishi2015 K=370
    alpha = -0.0122455748584756;
    beta  =  0.0873961711644606;
    gam   =  5.0/3.0;
    mu1   =  2.02*HBARC;
    mu2   =  1.0*HBARC;
    C1    = -0.38314;
    C2    =  0.33741;

  } else if(eosType==7) { // Skyrme + mom.dep Ohnishi2015 K=272.59

    alpha = -0.208886959978512;
    beta  = 0.284037556284497;
    gam   = 7.0/6.0;
    C1  =  -0.38314;
    C2  =   0.33741;
    mu1 = 2.02*HBARC;
    mu2 = 1.0*HBARC;

  } else if(eosType==8) { // MH2 Skyrme + mom.dep Nara2019 K=380 m*/m=0.848
    alpha = -0.0065641639085685915;
    beta = 0.08219105302398698;
    gam = 1.7228705539199192;
    C1= -0.386579535843693;
    C2= 0.3439279165932326;
    mu1= 2.02*HBARC;
    mu2= 1.0*HBARC;

//  } else if(eosType==9) { // MS2 Skyrme + mom.dep Nara2019 K=245 m*/m=0.8768
//    alpha = -1.1778187011310486;
//    beta = 1.253445590246467;
//    gam = 1.035433123202534;
//    C1= -0.38657953584369276;
//    C2= 0.34392791659323196;
//    mu1= 2.02*HBARC;
//    mu2= 1.0*HBARC;

  } else if(eosType==9) { // MS2 Skyrme + mom.dep Nara2019 K=210
    alpha = -0.31501938078689967;
    beta = 0.3874745511416311;
    gam = 1.1127892793311491;
    C1= -0.38657953584369287;
    C2= 0.34392791659323213;
    mu1= 2.02*HBARC;
    mu2= 1.0*HBARC;

  } else if(eosType==10) { // MH3 Skyrme + mom.dep Nara2019 K=380 m*/m=0.72
    alpha = 0.029643038556651874;
    beta = 0.04896806260897878;
    gam = 2.121401412903513;
    C1= -0.20263202109735332;
    mu1= 2.8*HBARC;
    C2= 0.06002279810547833;
    mu2= 1.0*HBARC;

  } else if(eosType==11) { // MS3 Skyrme + mom.dep Nara2019 K=210 m*/m= 0.7394
    //alpha = -7.019484656837315;
    //beta = 6.992149609691952;
    //gam = 1.0109840903345495;
    //C1= -0.1697658694585711;
    //C2= 0.33741;
    //mu1= 3.172915754445973*HBARC;
    //mu2= 1.0*HBARC;

    alpha = -0.6765491747331438;
    beta = 0.7551602758987743;
    gam = 1.047703735031756;
    C1= -0.20263202109735337;
    C2= 0.06002279810547853;
    mu1= 2.8*HBARC;
    mu2= 1.0 *HBARC;

  } else if(eosType==12) { // MH4 Skyrme + mom.dep Nara2019 K=380 m*/m=0.651
     alpha = 0.03862935265120862;
     beta = 0.04169348484309908;
     gam = 2.280403522511386;
     C1= -0.16976586945727515;
     C2= 0.0;
     mu1= 3.172915754474445*HBARC;
     mu2= 1.0*HBARC;

  } else if(eosType==13) { // MS4 Skyrme + mom.dep Nara2019 K=210 m*/m=0.705395
    alpha = -0.2068210445187974;
    beta = 0.28714388661643747;
    gam = 1.120133469619578;
    C1= -0.16976586945727515;
    C2= 0.0;
    mu1= 3.172915754474445*HBARC;
    mu2= 1.0*HBARC;

    //alpha = -0.20784369609866815;
    //beta = 0.2881665335936613;
    //gam = 1.1197071552352167;
    //C1= -0.16976586945850788;
    //mu1= 3.1729157544574456*HBARC;
    //C2= 0.0;  mu2= 1.0*HBARC;

  } else if(eosType==14) { // MH2 Skyrme + mom.dep Nara2021 K=380 Hama2
    alpha = -0.036917343467595724;
    beta = 0.10277983945195997;
    gam = 1.5851487200433838;
    C1= -0.7163341712031479;
    C2= 0.7526821171198502;
    mu1= 1.55 *HBARC;
    mu2= 1.0 *HBARC;

  } else if(eosType==15) { // MS2 Skyrme + mom.dep Nara2021 K=210 Hama2
    alpha = -0.8885528341665909;
    beta = 0.9544152142736043;
    gam = 1.0432230334789112;
    C1= -0.7163341712031479;
    C2= 0.7526821171198502;
    mu1= 1.55 *HBARC;
    mu2= 1.0 *HBARC;

  } else if(eosType==20) {  // Ohinishi 2022 M1 K=250MeV
    alpha = -202.8*1e-3;
    beta =  132.4*1e-3;
    gam = 4./3.;
    t5=18.3*1e-3/(1.+5./3.)/pow(rho0,5./3.);
    C1 = 0.0;
    mu1 = 3.2272077792680443*HBARC;
    withMomDep=0;

  } else if(eosType==21) {  // Ohinishi 2022 MM1 K=250MeV
    alpha = -21.4*1e-3;
    beta =  92.1*1e-3;
    gam = 4./3.;
    t5=12.6*1e-3/(1.+5./3.)/pow(rho0,5./3.);
    C1 = -0.16975960664857545;
    mu1 = 3.2272077792680443*HBARC;

  } else if(eosType==100) {
    beta = 125.93*1e-3;
    gam = 1.676;
    alpha = -0.93*2.0*87.67*1e-3;
    withMomDep=0;

  } else {
    cout << "NonRelPotential::EoStype not implemented " << eosType << endl;
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
  cout << "# eosType= "<< eosType  << " rho0= " << rho0
      << " alpha= " << alpha
      << " beta= "<< beta << " gamma= "<< gam
      << " t5 = "<< t5
      << " withmomdep = "<< withMomDep
      <<endl;
    firstCall=false;
  }
}

void NonRelPotential::setPotentialParam(EventParticle* p)
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

void NonRelPotential::setSigmaPotential(double* fp)
{
  // Kohno-Kohno 2 body
  double gamS1=4.0/3.0;
  double gamS2=5.0/3.0;
  double aSigma= -17.0883105169774;  // in MeV
  double bSigma= 118.512834395095; // in MeV
  double cSigma= -24.0500980143978; // in MeV
  double c1= -78.041869780504; // in MeV
  double c2=0.0;
  double mu1=3.23*HBARC; // in GeV
  //double mu1 = 3.2272077792680443*HBARC;
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

void NonRelPotential::setLambdaPotential(double* fp)
{
// Jinno 2022/1/5
//=1: 2+3BF upper
//=2: 2+3BF lower
//=3: 2BF   upper
//=4: 2BF   lower
//=5: GKW3 (2+3-body) average
//=6: GKW3 (2-body) average
//=7: GKW3 + Kohno (2+3-body) average
//=8: GKW3 + Kohno (2-body) average
//=9: Kohno+Kohno (2+3body) average
//=10:Kohno+Kohno (2body) average

  // Kohno-Kohno 2+3 body
  double gamL1=4.0/3.0;
  double gamL2=5.0/3.0;
  double aLambda= 15.441112151603;  // in MeV
  double bLambda= 51.7156158787796; // in MeV
  double cLambda= 9.30908485961732; // in MeV
  double c1=-116.7513727; // in MeV
  double c2=0.0;
  double mu1=3.23*HBARC; // in GeV
  //double mu1 = 3.2272077792680443*HBARC;
  double mu2=1.0*HBARC; // in GeV

  if(optLambdaPot<=4) {
    // Ohnishi 2021/12/31
    //double J[4]={-29.28, -29.61, -32.78, -33.61};
    //double L[4]={16.11, 6.19, -0.43, -7.};
    //double K[4]={556.28, 447.66, 360.65, 345.1};

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
    aLambda=40.59694627;
    bLambda= -10.42255776;
    cLambda=  46.29142439;
    c1= -116.7513727;
    c2= 0;
  /*
  2+3body force	lower 28.90026782 18.24560978 29.3199353  -116.7513727	  3.23
   lower -7.327857338	21.38895389 27.48728513	-102.9620846  20.29065006 2.02
   upper 52.29362471	-39.09072529  63.26291347 -116.7513727	  3.23
   upper 15.67772195 -35.28818956	61.1588493  -102.9620846  20.29065006 2.02
 */

  // GKW2+Kohno (2body) average
  } else if(optLambdaPot==8) {
    aLambda=-10.32799492;
    bLambda= 83.55785844;
    cLambda=-8.307349732;
    c1= -104.8008152;
    c2= 0;

  // Kohno+Kohno (2+3body) average
  } else if(optLambdaPot==9) {
    aLambda= 15.441112151603;  // in MeV
    bLambda= 51.7156158787796; // in MeV
    cLambda= 9.30908485961732; // in MeV
    c1=-116.7513727; // in MeV
    c2=0.0;

  // Kohno+Kohno (2body) average
  } else if(optLambdaPot==10) {
    aLambda= -4.770236003;
    bLambda= 72.97351005;
    cLambda= -3.280760248;
    c1= -104.8008152;
    c2= 0;


  /*
                   a(MeV)           b(MeV)          c(MeV)          C1(MeV)      C2(MeV)  μ1 mu2
2BF	lower	0.5212105714	67.4130712	-3.01176798	-104.8008152		3.23
	lower	-31.71118355	68.41058175	-3.973140057	-90.09173048	16.41984089	2.02	1
	upper	-21.1772004	99.70264568	-13.60293148	-104.8008152		3.23
	upper	-53.27114434	100.4676676	-14.47026513	-90.09173048	16.41984089	2.02	1
   */
  } else if(optLambdaPot>10 && optLambdaPot<20){
    double J2[9]={-30.0,-30.0,-30.0, -30.0,-30.0,-30.0, -30.0,-30.0,-30.0};
    double L2[9]={-90.0,-90.0,-90.0,   0.0,  0.0,  0.0,  15.0, 15.0, 15.0};
    double K2[9]={  0.0,350.0,500.0,   0.0,350.0,500.0,   0.0,350.0,500.0};
    double JL2=0.0,LL2=0.0,KL2=0.0;
    JL2=J2[optLambdaPot-11];
    LL2=L2[optLambdaPot-11];
    KL2=K2[optLambdaPot-11];

    aLambda= 10*JL2-3*LL2+KL2/2;
    bLambda=-15*JL2+5*LL2-KL2;
    cLambda=  6*JL2-2*LL2+KL2/2;
    c1=0.0; c2=0.0;
  } else if(optLambdaPot>=20 && optLambdaPot<=25){
  //                19    20     21    22    23    24
    double J3[6]={-50.0,-30.0, 10.0, -50.0,-30.0, 10.0};
    double L3[6]={ 15.0, 15.0, 15.0,  -2.0, -2.0,  -20};
    double K3[6]={500.0,500.0,500.0, 350.0,350.0,350.0};
    double JL3=0.0,LL3=0.0,KL3=0.0;
    JL3=J3[optLambdaPot-20];
    LL3=L3[optLambdaPot-20];
    KL3=K3[optLambdaPot-20];

    aLambda= 10*JL3-3*LL3+KL3/2;
    bLambda=-15*JL3+5*LL3-KL3;
    cLambda=  6*JL3-2*LL3+KL3/2;
    c1=0.0; c2=0.0;
    // cout << "# a=" << aLambda << ", b=" << bLambda << ", c=" << cLambda << endl;
  } else if(optLambdaPot==26){
    // GKW3+Kohno/2
    mu1=3.23*HBARC; // in GeV
    mu2=1.0*HBARC; // in GeV
    aLambda=-19.47;bLambda=-5.569;cLambda=48.515;c1=-116.75/2;c2=0.0;
  } else if(optLambdaPot==27){
    // GKW3+Kohno/3
    mu1=3.23*HBARC; // in GeV
    mu2=1.0*HBARC; // in GeV
    aLambda=-39.69;bLambda=-3.660;cLambda=49.149;c1=-116.75/3;c2=0.0;
  } else if(optLambdaPot==28){
    // GKW3+Kohno/5
    mu1=3.23*HBARC; // in GeV
    mu2=1.0*HBARC; // in GeV
    aLambda=-55.86;bLambda=-2.13;cLambda=49.66;c1=-23.35;c2=0.0; // GKW3+Kohno/5 (2+3-body)
  } else {
    cout << "NonRelPotential::setLambdaPotential wrong optLmambdaPot = "<< optLambdaPot<<endl;
    exit(1);
  }

  fp[0]= facPotL[0];
  fp[3]= gamL1;
  fp[5]= gamL2;

  fp[1]= fp[0]*aLambda*1e-3/rho0/t1;
  fp[2]= fp[0]*bLambda*1e-3/pow(rho0,gamL1)/t3;

  //double f3= t5 == 0.0 ? 1.0 : t5;
  fp[4]= fp[0]*cLambda*1e-3/pow(rho0,gamL2)/t5;

  fp[6]=0.0;
  fp[7]=0.0;
  if(vex1 !=0.0) fp[6] = c1*1e-3/(rho0)/vex1;
  if(vex2 !=0.0) fp[7] = c2*1e-3/(rho0)/vex2;

  fp[8] = mu1*mu1/pmu1;
  fp[9] = mu2*mu2/pmu2;
  if(optPotentialDensity>0) {
    fp[1] /= 2;
    fp[2] /= (gamL1+1.0);
    fp[4] /= (gamL2+1.0);
    fp[6] /=2;
    fp[7] /=2;
  }

}

} // namespace jam2
