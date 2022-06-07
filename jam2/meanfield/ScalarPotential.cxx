#include <jam2/meanfield/ScalarPotential.h>
#include <jam2/hadrons/JamStdlib.h>

namespace jam2 {

using namespace std;

// scalar potential --------------------------------------------------------------------------
ScalarPotential::ScalarPotential(Pythia8::Settings* s) : PotentialType(s)
{
  optPotentialType=1;
  //potParam = new PotentialParam(settings,optPotentialType);

  withMomDep=1;
  double C1=0.0, C2=0.0, mu1=1.0, mu2=1.0;

  if(eosType==1) {     // Skyrme Hard Nara2019 K=380 
    //alpha = -124.0e-3;
    //beta = 70.5e-3;
    //gam = 2.0;

    alpha = -0.12384212774318967;
    beta = 0.06875425035612921;
    gam = 2.123648826081624;
    withMomDep=0;

  } else if(eosType==2) { // Skyrme Soft Nara2019 K=210
    //alpha = -356.0e-3;
    //beta = 303.0e-3;
    //gam = 1.16667;

    alpha = -0.2621174235726939;
    beta = 0.20645616373308995;
    gam = 1.2652357572720947;
    withMomDep=0;

  } else if(eosType==3) { // MH2 Skyrme + mom.dep Nara2019 K=380
    alpha = 0.00877512904142221;
    beta = 0.07048243772253593;
    gam = 1.8907272542475995;
    C1 = -0.3907221709140509;
    C2= 0.3341868667966711;
    mu1 = 2.02*HBARC;
    mu2= 1.0*HBARC;

  } else if(eosType==4) { // MS2 Skyrme + mom.dep Nara2019 K=210

    // U( Pinf= 1.7 )= 0.06  Einf= 1.0036086114353737 Pinf= 1.7
    // U( 0.65 )= 0.0 Elab= 0.20320287416392357
    alpha = -0.07785801989125717;
    beta = 0.1565487207865633;
    gam = 1.261622285760052;
    C1 = -0.3907221709140509;
    C2= 0.3341868667966711;
    mu1 = 2.02*HBARC;
    mu2= 1.0*HBARC;

  } else if(eosType==5) { // MH3 Skyrme + mom.dep Nara2019 K=380

    // U( Pinf= 1.7 )= 0.06  Einf= 1.0036086114353737 Pinf= 1.7
    // U( 0.65 )= 0.0 Elab= 0.20320287416392357
    // meff=  0.8823864683629444  M*/M= 0.9407105206427979 rho0= 0.168
    alpha = 0.03911636591740648;
    beta = 0.043064133086746524;
    gam = 2.36458860928333;
    C1 = -0.2048087947472817;
    mu1 = 2.8*HBARC;
    C2= 0.04713369448491099;
    mu2= 1.0*HBARC;

  } else if(eosType==6) { // MS3 Skyrme + mom.dep Nara2019 K=210

    // U( Pinf= 1.7 )= 0.06  Einf= 1.0036086114353737 Pinf= 1.7
    // U( 0.65 )= 0.0 Elab= 0.20320287416392357
    // meff=  0.8823864683629444  M*/M= 0.9407105206427979 rho0= 0.168
    alpha = -0.04155376493236816;
    beta = 0.1231543745372214;
    gam = 1.294657727096535;
    C1 = -0.2048087947472817;
    mu1 = 2.8*HBARC;
    C2= 0.04713369448491099;
    mu2= 1.0*HBARC;

    /*
    alpha = -0.6285006244768873;
    beta = 0.7070852697952635;
    gam = 1.0496537865787445;
    C1 = -0.20192630016292473;
    C2= 0.05694208245922472;
    mu1 = 2.8*HBARC;
    mu2= 1.0*HBARC;
    */

  } else if(eosType==7) { // MH4 Skyrme + mom.dep Nara2019 K=380
    alpha = 0.045166262659641;
    beta = 0.03824971687635796;
    gam = 2.509835588271163;
    C1 = -0.1787443736909779;
    mu1 = 3.0754425376340975*HBARC;

  } else if(eosType==8) { // MS4 Skyrme + mom.dep Nara2019 K=210
    alpha = -0.0896591573281629;
    beta = 0.17249044531651028;
    gam = 1.202859798739601;
    C1 = -0.1787443736909779;
    mu1 = 3.0754425376340975*HBARC;

  //--------------------------------------------------------------------------
  // optical potential is defined by sqrt{(m_N+S)^2+p^2} - sqrt{m_N^2 + p^2}
  // m*/m=0.88
  } else if(eosType==11) { // NH2 Skyrme + mom.dep Nara2021 K=380

    // U( Pinf= 1.7 )= 0.05  Einf= 1.0036086114353737 Pinf= 1.7
    // U( 0.685 )= 0.0 Elab= 0.22349429615474214
    // meff=  0.8823864683629444  M*/M= 0.9407105206427979 rho0= 0.168
    alpha = -0.14746618076876553;
    beta = 0.28231864726575806;
    gam = 1.3086087855271107;
    C1 = -0.8500057544566262;
    mu1 = 2.02*HBARC;
    C2= 1.0508877826863514;
    mu2= 1.0*HBARC;

    /*
    alpha = 0.1436791875727167;
    beta = 0.04604753299756481;
    gam = 2.3300535324726885;
    C1 = -0.272173218332104;
    mu1 = 5.077642517924502;
    */

  } else if(eosType==12) { // MS2 Skyrme + mom.dep Nara2021 K=210
    // U( Pinf= 1.7 )= 0.05  Einf= 1.0036086114353737 Pinf= 1.7
    // U( 0.685 )= 0.0 Elab= 0.22349429615474214
    // meff=  0.8823864683629444  M*/M= 0.9407105206427979 rho0= 0.168
    alpha = -1.7403675255660511;
    beta = 1.8746769579585103;
    gam = 1.0351617222744887;
    C1 = -0.8500057544566262;
    C2= 1.0508877826863514;
    mu1 = 2.02*HBARC;
    mu2= 1.0*HBARC;

  } else if(eosType==13) { // NH4 Skyrme + mom.dep Nara2021 K=380 Hama1
    // U( Pinf= 1.7 )= 0.06  Einf= 1.0036086114353737 Pinf= 1.7
    // U( 0.685 )= 0.0 Elab= 0.22349429615474214
    // meff=  0.8823864683629444  M*/M= 0.9407105206427979 rho0= 0.168
    alpha = 0.16795158543549493;
    beta = 0.04930245416029725;
    gam = 2.286443474993692;
    C1 = -0.2952228390578033;
    mu1 = 5.880351961800975*HBARC;

  } else if(eosType==14) { // NH4 Skyrme + mom.dep Nara2021 K=210 Hama1
    // U( Pinf= 1.7 )= 0.06  Einf= 1.0036086114353737 Pinf= 1.7
    // U( 0.685 )= 0.0 Elab= 0.22349429615474214
    // meff=  0.8823864683629444  M*/M= 0.9407105206427979 rho0= 0.168
    alpha = -0.014163922715488364;
    beta = 0.23083815155546436;
    gam = 1.176883586803349;
    C1 = -0.2952228390578033;
    mu1 = 5.880351961800975*HBARC;

  } else if(eosType==15) { // NH4-2 Skyrme + mom.dep Nara2021 K=380 Hama2
    // U( Pinf= 1.7 )= 0.052  Einf= 1.0036086114353737 Pinf= 1.7
    // U( 0.65 )= 0.0 Elab= 0.20320287416392357
    // meff=  0.8823864683629444  M*/M= 0.9407105206427979 rho0= 0.168
    alpha = 0.11908916113172166;
    beta = 0.04506950566822351;
    gam = 2.3485654474084523;
    C1 = -0.2480442014823663;
    mu1 = 4.670304914324092*HBARC;

  } else if(eosType==16) { // NS4-2 Skyrme + mom.dep Nara2021 K=210 Hama2
    // U( Pinf= 1.7 )= 0.052  Einf= 1.0036086114353737 Pinf= 1.7
    // U( 0.65 )= 0.0 Elab= 0.20320287416392357
    // meff=  0.8823864683629444  M*/M= 0.9407105206427979 rho0= 0.168
    alpha = -0.061588753664328644;
    beta = 0.22516601093851452;
    gam = 1.1694228159044495;
    C1 = -0.2480442014823663;
    mu1 = 4.670304914324092*HBARC;



  } else if(eosType==21) { // Skyrme + mom.dep Ohnishi2015 K=272.59
    alpha = -0.208886959978512;
    beta  = 0.284037556284497;
    gam   = 7.0/6.0;
    mu1 = 2.02*HBARC;
    mu2 = 1.0*HBARC;
    C1  =  -0.38314;
    C2  =   0.33741;

  } else if(eosType==22) { // Skyrme  K=550
    alpha = -0.11096186950393147;
    beta = 0.05645200738039593;
    gam = 2.774173096033283;
    withMomDep=0;

  } else if(eosType==23) { // Skyrme + mom.dep K=450
    alpha = 0.050658464813240926;
    beta = 0.03300384907908334;
    gam = 3.0470954207158307;
    C1 = -0.17874437370225157;
    mu1 = 3.0754425374321817;



  } else {
    cout << "RQMDs mode MeanField:EoStype not implemented " << eosType << endl;
    exit(1);
  }

  pmu1= mu1*mu1;
  pmu2= mu2*mu2;

  t5=1.0;
  if(transportModel==1) {
  t0 = 0.5*alpha/rho0;
  t2 = beta/(gam+1.0)/pow(rho0,gam);
  t2f = gam*t2;
  vex1=C1/(2*rho0);
  vex2=C2/(2*rho0);

  } else {
    t0 = alpha/rho0;
    t2 = beta/pow(rho0,gam);
    t2f = gam*t2;
    vex1=C1/rho0;
    vex2=C2/rho0;
  }

  t1=0.0;
  t3=0.0;
  t3f=0.0;

  if(firstCall) {
  cout << "# RQMDs mode eosType= "<< eosType  << " rho0= " << rho0
      << " alpha= " << alpha
      << " beta= "<< beta << " gamma= "<< gam
      <<endl;
    firstCall=false;
  }

}


// scalar potential -------------------------------------------------------------------------------
void ScalarPotential::Vmd(double psq1,double psq2,double vfac1,double vfac2,double pfi1,double pfi2,double pfj1,double pfj2)
{
  vmomsi = vfac1*vex1/(1.0-psq1/(pfi1*pmu1)) + vfac2*vex2/(1.0-psq1/(pfi2*pmu2));
  vmomsj = vfac1*vex1/(1.0-psq2/(pfj1*pmu1)) + vfac2*vex2/(1.0-psq2/(pfj2*pmu2));
  vmom4i=0.0;
  vmom4j=0.0;
}

void ScalarPotential::dVdns(double rhij,double rhji)
{
  rhomij=rhij;
  rhomji=rhji;

  s1 = pfac1*pfac2*t0 + gama1*t2*pfac3*pfac4*rhosgi + gamb1*pfac5*pfac6*pow(rhosi,gamb1-1);
  s2 = pfac1*pfac2*t0 + gama2*t2*pfac4*pfac3*rhosgj + gamb2*pfac5*pfac6*pow(rhosj,gamb2-1);

  //s1 = pfac1*t0*pfac2 + gama1*t2*pfac3*rhosgi + t5*gamb1*pfac5*pow(rhosi,gamb1-1);
  //s2 = pfac2*t0*pfac1 + gama2*t2*pfac4*rhosgj + t5*gamb2*pfac6*pow(rhosj,gamb2-1);

  fskyi = -wG*fengi*s1*fj*rhomij;
  fskyj = -wG*fengj*s2*fi*rhomji;
}

void ScalarPotential::dVdmd(double psq1,double psq2,double fs,double fv,
    double pmi1,double pmj1,double pmi2,double pmj2)
{
  fmomdi = -wG*fengi*devVmd2(psq2,fs,fv,pmi1,pmi2)*rhomij*fj;
  fmomdj = -wG*fengj*devVmd2(psq1,fs,fv,pmj1,pmj2)*rhomji*fi;
  fmomei =     fengi*devVme2(psq2,fs,fv,pmi1,pmi2)*rhomij*fj;
  fmomej =     fengj*devVme2(psq1,fs,fv,pmj1,pmj2)*rhomji*fi;
}

// Derivative of p^\mu/m_j term in the vector potential.
void ScalarPotential::dVdp()
{
  if(optTwoBodyDistance != 3) {
    forceri = fskyj*v1/(wG*p01)*optP0dev;
    forcerj = fskyi*v2/(wG*p02)*optP0dev;
  }else {
    forceri=0.0;
    forcerj=0.0;
  }
  facsk= -(fskyi+fskyj)/wG;
}

void ScalarPotential::dVdpm()
{
  if(optTwoBodyDistance != 3) {
    forceri = fmomdj*v1/(wG*p01)*optP0dev;
    forcerj = fmomdi*v2/(wG*p02)*optP0dev;
  } else {
    forceri=0.0;
    forcerj=0.0;
  }
  facmom = -(fmomdi + fmomdj)/wG;
}

} // namespace jam2
