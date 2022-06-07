#include <jam2/meanfield/ScalarVectorPotential.h>
#include <jam2/hadrons/JamStdlib.h>

namespace jam2 {

using namespace std;
ScalarVectorPotential::ScalarVectorPotential(Pythia8::Settings* s) : PotentialType(s)
{
  if(optLambdaPot>0) {
    cout << "MeanField:optLambdaPotential > 0 cannot be used in RMF mode" <<endl;
    exit(1);
  }

  optPotentialType=3;

  double C1=0.0, C2=0.0, mu1=1.0, mu2=1.0;
  withMomDep=1;
  double Cv,Cs,Cd;
  g4 = 0.0;

  // non-linear vector potential
   // K=380MeV meff=  0.5658193115240782  M*/M= 0.6032188822218317 rho0= 0.168
  if(eosType==1) {

    Cs= 306.72053219515175; // GeV^-2
    Cv= 239.85427189549577;  // GeV^-2
    Cd= 0.123197568454716;
    withMomDep=0;

  // K=210 MeV meff=  0.6483817101103596  M*/M= 0.69123849691936 rho0= 0.168
  } else if(eosType==2) {

    Cs= 235.3864964552273;
    Cv= 192.23574594599532;
    Cd= 0.25240953167965297;
    withMomDep=0;

  //K=170 meff=  0.6772920823786632  M*/M= 0.7220597893162721 rho0= 0.168
  } else if(eosType==3) {

    Cs= 211.08806268066527;
    Cv= 174.24894924265055;
    Cd= 0.28540760488970235;
    withMomDep=0;

  // K=100 MeV meff=  0.7487859111178186  M*/M= 0.7982792229401051 rho0= 0.168
  } else if(eosType==4) {

    Cs= 152.0308580296217;
    Cv= 127.68894231421507;
    Cd= 0.34816130194083716;
    withMomDep=0;

    // MV1
  } else if(eosType==11) {  //MD1 like  K=272MeV  m*/m=0.65
   //U( Pinf= 1.7 )= 0.058, U( 0.65 )= 0.0, U( 3.0 )= 0.065

    Cs=  241.97268739051128;
    Cv=  90.45600092469897;
    Cd=  0.19387957509551007;
    gam=  0.3333333333333333;
    //C1 = 33.547026068603934; // GeV^-2
    //C2= 129.06726683835495;
    C1= 0.04330347090157576;  // GeV/rho0
    C2= 0.16660375862990676;
    mu1 = 0.6411295603162918;
    mu2 = 1.8413488870476014; // in GeV

    // MV2
  } else if(eosType==12) {  //MD2 like  K=272 MeV m*/m=0.65
  // U( Pinf= 1.7 )= 0.056  U( 0.65 )= 0.0  U( 6.0 )= 0.025

    Cs=  254.1301519559795;
    Cv=  41.015359694798384;
    Cd=  0.19387464043557337;
    gam=  0.3333333333333333;
    // C1 = 20.69381998040592 C2= 177.50261883100623 in GeV^-2
    C1 = 0.026712180970440987;
    C2= 0.2291255109704515;
    mu1 = 0.4896732480502826;
    mu2 = 2.4886149663285533; //  in GeV

  } else if(eosType==13) {  //MD4 like  K= 272 GeV M*/M= 0.65
  // U( Pinf= 1.7 )= 0.06  U( 0.65 )= 0.0  U( 5.0 )= 0.07 Elab= 4.148408082707736
   
    Cs=  248.5126151508354;
    Cv=  90.45669623967359;
    Cd=  0.1938829214447074;
    gam=  0.3333333333333333;
    // C1 = 26.702866396661847 C2= 128.62465992254627
    C1 = 0.03446883177163635;
    C2 = 0.16603242882975386; //  GeV
    mu1 = 0.5652108383665899;
    mu2 = 1.9441947692051738; //  in GeV

    // MV3
  } else if(eosType==14) {  //MD2 like  K=163MeV m*/m=0.75
   //U( Pinf= 1.7 )= 0.056, U( 0.65 )= 0.0, U( 6.0 )= 0.025

    Cs=  132.88357280361947;
    Cv=  33.611953709186515;
    Cd=  0.2561535832076794;
    gam=  0.3333333333333333;
    //C1 = 70.05807993746818;
    //C2= 120.91980188731048;
    C1 = 0.09043299456101145;
    C2= 0.15608677537454843;
    mu1 = 0.6699871063319663;
    mu2 = 3.1108637632415297;

  } else if(eosType==15) {  // MD4 like K= 130 GeV M*/M= 0.8
  // U( Pinf= 1.7 )= 0.06  U( 0.65 )= 0.0  U( 5.0 )= 0.07 Elab= 4.148408082707736
   // momentum-dep. pot is the same as MD4
 
    Cs=  77.17220896448848;
    Cv=  82.06638123140485;
    Cd=   0.2739146565236662;
    gam=  0.3333333333333333;
    // C1 = 90.1797480774058 C2= 39.004302594584075
    C1 = 0.11640662539798463;
    C2= 0.05034788117993155; //  GeV
    mu1 = 0.6919208986787576;
    mu2 = 2.6668979434247753; //  in GeV

  } else if(eosType==16) {  // MD4 like K= 121 GeV M*/M= 0.83
  // U( Pinf= 1.7 )= 0.06  U( 0.65 )= 0.0  U( 5.0 )= 0.07 Elab= 4.148408082707736
   // momentum-dep. pot is the same as MD4
  
    Cs=  44.410898422960116;
    Cv=  75.49691074521036;
    Cd=  0.28157627938305196;
    gam=  0.3333333333333333;
    // C1 = 101.60435633640503 C2= 25.148757563277442
    C1 = 0.13115383995864768;
    C2 = 0.03246274316912446; //  GeV
    mu1 = 0.7035464013689996;
    mu2 = 4.252421078959087; //  in GeV

 

//--  gamma not equal to 1/3  ----------------------------------------------
  } else if(eosType==17) {  // M4 like K= 210 GeV M*/M= 0.8
// U( Pinf= 1.7 )= 0.06  Einf= 1.0030919133758833 Pinf= 1.7
// U( 0.65 )= 0.0 Elab= 0.20302495594448378
// U( 5.0 )= 0.07 Elab= 4.148408082707736

    Cs=  77.17220896448848;
    Cv=  289.4932142949801;
    Cd=  122.1765864997752;
    gam=  0.9045586477893912;
    //C1 = 90.1797480774058 C2= 39.004302594584075
    C1 = 0.11640662539798463;
    C2= 0.05034788117993155; //  GeV
    mu1 = 0.6919208986787576;
    mu2 = 2.6668979434247753; //  in GeV
      

  } else if(eosType==18) {  // M4 like K= 210 GeV M*/M= 0.81
  // U( Pinf= 1.7 )= 0.06  Einf= 1.0030919133758833 Pinf= 1.7
  // U( 0.65 )= 0.0 Elab= 0.20302495594448378
  // U( 5.0 )= 0.07 Elab= 4.148408082707736

   Cs=  66.32741826135909;
   Cv=  388.06042038941143;
   Cd=  210.0253015399415;
   gam=  0.9318087212484244;
   // C1 = 93.89977948877997 C2= 33.463607863999194
   C1 = 0.12120854946857387;
   C2 = 0.043195792271963344;  //  GeV
   mu1 = 0.695537617417831;
   mu2 = 2.8944623003059498; //  in GeV

  } else if(eosType==19) {  // M4 like K= 343 GeV M*/M= 0.61
  // U( Pinf= 1.7 )= 0.06  Einf= 1.0030919133758833 Pinf= 1.7
  // U( 0.65 )= 0.0 Elab= 0.20302495594448378
  // U( 5.0 )= 0.07 Elab= 4.148408082707736
  // initial Usk=  0.13968202562672366

    Cs=  296.679048289649;
    Cv=  90.53454462473317;
    Cd=  0.15439847979476506;
    gam= 0.3333333333333333;
    //# C1 = 8.09280463137856 C2= 152.8444774753584
    C1 = 0.010446426134783;
    C2= 0.197296069383038; //  GeV
    mu1 = 0.35998936003490106;
    mu2 = 1.8696352493989794; //  in GeV




// non-linear scalar potential ----------------------------------------------------------------
  // K=380 MeV  M*/M= 0.597 rho0= 0.168
  } else if(eosType==21) {
    Cs= 300.8290744958282;
    Cv= 233.2877062500283;
    Cd= 0.1221503440591248;
    gam=  0.3333333333333333;
    withMomDep=0;
    optPotentialType=4;

  } else if(eosType==22) {  // K= 252 GeV M*/M= 0.65  MD1 like
  // U( Pinf= 1.7 )= 0.058 U( 0.65 )= 0.0 U( 3.0 )= 0.065 Elab= 2.204520478699002

    Cs=  223.76408272256893;
    Cv=  74.10216618086122;
    Cd=  0.20827481975232537;
    gam=  0.3333333333333333;
    // C1 = 33.54702606857792 C2= 129.06726683835328
    C1 = 0.04330347090157576;
    C2 = 0.16660375862990676;
    mu1 = 0.6411295603159707;
    mu2 = 1.841348887047491; //  in GeV
    optPotentialType=4;

  } else if(eosType==23) {  // K= 252 GeV M*/M= 0.65 MD2 like
  // U( Pinf= 1.7 )= 0.056  U( 0.65 )= 0.0 U( 6.0 )= 0.025

    Cs=  235.92176813150945;
    Cv=  24.661941191880853;
    Cd=  0.20827136868731128;
    gam=  0.3333333333333333;
    // C1 = 20.69381998040592 C2= 177.50261883100623
    C1 = 0.026712180970440987;
    C2= 0.2291255109704515;
    mu1 = 0.4896732480502826;
    mu2 = 2.4886149663285533; //  in GeV
    optPotentialType=4;

  } else if(eosType==24) {  // K= 346 GeV M*/M= 0.6 MD4 like
  // U( Pinf= 1.7 )= 0.06 U( 0.65 )= 0.0 U( 5.0 )= 0.07 Elab= 4.148408082707736

    Cs=  294.30006535947695;
    Cv=  78.32240888710614;
    Cd=  0.15504973979210573;
    gam=  0.3333333333333333;
    // C1 = 5.308698422700437 C2= 158.9441887837254
    C1 = 0.0068526213680673515;
    C2 = 0.20516975304757276;
    mu1 = -0.19612260187380054;
    mu2 = 1.8469370105201242;  //  in GeV
    optPotentialType=4;


// non-linear sigma-field ---------------------------------------------------------------------
     // mstc(106).eq.153
  } else if(eosType==31) {   // New NS1 Hard K=380 m*/m=0.83
    gs=  6.447696268747421;
    g2= -38.001619264755874;
    g3= 339.59550852599887;
    gv= 6.858858568503943;
    gsp=0.0; gvp=0.0;
    mu1 = 1.0; mu2 = 1.0;
    withMomDep=0;
    optPotentialType=5;

  } else if(eosType==32) {   // New NS2 Soft K=210 m*/m=0.83

    gs=  7.901789322846869;
    g2= 44.31272165001529;
    g3= 21.98880833814041;
    gv= 6.858858568503943;
    gsp=0.0; gvp=0.0;
    mu1 = 1.0; mu2 = 1.0;
    withMomDep=0;
    optPotentialType=5;

  } else if(eosType==33) {  //MD1  K=380 MeV m*/m=0.65 Hama1 U=0.065 at plab=3
    gs=  9.02956219322301;
    g2= 4.218109600022776;
    g3= 6.667054503467467;
    gv= 6.7402687603431675;
    gsp =  3.1855887031669394;
    gvp =  8.89548883191133;
    mu1 = 0.6411295603159707;
    mu2 = 1.841348887047491;
    optPotentialType=5;
    gam = 3.0;


  } else if(eosType==34) {  // MD2   = JAM1 204 
    // K=380 MeV m*/m=0.65 U(1.7)=0.056  U(6)=0.025
    gs=  9.232564820463068;
    g2= 4.012495169512666;
    g3= 5.520251163159458;
    gv= 3.8884399009614437;
    gsp =  2.5019753284300768;
    gvp =  10.431917516759935;
    mu1 = 0.4896732480502826;
    mu2 = 2.4886149663285533;
    gam = 3.0;
    optPotentialType=5;

  // MD4   210 in JAM1
  } else if(eosType==35) { // K=210 MeV m*/m=0.83 U(1.7)=0.06 U(5)=0.075

    gs=  4.0586271572892025;
    g2= -160.31303134630116;
    g3= 2684.378830396555;
    gv= 5.63247647754081;
    gsp =  5.543944245008469;
    gvp =  3.9266304417034466;
    mu1 = 0.7035464013689996;
    mu2 = 4.252421078959087;
    gam = 3.0;
    optPotentialType=5;


  // MD5
  } else if(eosType==36) { // K=380 MeV m*/m=0.65 U(1.7)=0.06 U(4)=0.06
// U( Pinf= 1.7 )= 0.06  Einf= 1.0030919133758833 Pinf= 1.7
// U( 0.65 )= 0.0 Elab= 0.20302495594448378
// U( 4.0 )= 0.06 Elab= 3.169737153919681
     gs=  9.218555727507903;
     g2= 4.017566480669699;
     g3= 5.625329209008978;
     gv= 6.11270077474097;
     gsp =  2.5605990856596312;
     gvp =  9.309827349691702;
     mu1 = 0.4982896543033027;
     mu2 = 2.163386179379825;
     gam = 3.0;
     optPotentialType=5;

  } else if(eosType==37) { // modified MD1 K=380 MeV m*/m=0.65
  // m*/m =  0.65  M_sigma= 0.55 M_omega= 0.783
// This results are obtained by  /export/ynara/tex/paper/2021/qmd/sv/rmfp.py
// optMom =  0 # of Gauss point= 12
// U( Pinf= 1.7 )= 0.06  Einf= 1.0030919133758833 Pinf= 1.7
// U( 0.65 )= 0.0 Elab= 0.20302495594448378
// U( 3.0 )= 0.075 Elab= 2.204520478699002

     gs=  8.737517140523442;
     g2= 4.4948912175680595;
     g3= 8.915928956362198;
     gv= 7.477535505594842;
     gsp =  3.9518787360160923;
     gvp =  8.330662131818302;
     mu1 = 0.7675482504131891;
     mu2 = 1.506156163664052;
    gam = 3.0;
    optPotentialType=5;

  } else if(eosType==38) { // K=210 MeV m*/m=0.65 U(1.7)=0.06 U(5)=0.075
    gsp =  3.1855887031669394;
    gvp =  8.89548883191133;
    mu1 = 0.6411295603159707;
    mu2 = 1.841348887047491;
    gs=  9.531539908087288;
    g2= 21.884205176899933;
    g3= -63.40629377388759;
    gv= 6.7402687603431675;
    gam = 3.0;
    optPotentialType=5;

  } else {
    cout << "RQMDsv mode: MeanField:EoStype not implemented " << eosType << endl;
    exit(1);
  }

  t5=0.0;
  double CONV=HBARC*HBARC*HBARC;
  // non-linear vector potential.
  if(optPotentialType==3) {
    t0 = -0.5*Cs*CONV;
    t2 =  0.0;
    t2f =  0.0;
    t1 =  0.5*Cv*CONV;
    t3 = -Cd/(gam+1.0)*HBARC;
    t3f = gam*t3;
    vex1 = -C1/(2*rho0);
    vex2 =  C2/(2*rho0);

    if(transportModel==2) {
      t0 = -Cs*CONV;
      t1 =  Cv*CONV;
      t3 = -Cd*HBARC;
      t3f = gam*t3;
      vex1 *= 2;
      vex2 *= 2;
    }

  // non-linear scalar potential.
  } else if(optPotentialType==4) {
    t0 = -0.5*Cs*CONV;
    t2 = -Cd/(gam+1.0)*HBARC;
    t2f = t2*gam;
    t1 =  0.5*Cv*CONV;
    t3 =  0.0;
    t3f = 0.0;
    vex1 = -C1/(2*rho0);
    vex2 =  C2/(2*rho0);

    if(transportModel==2) {
      t0 = -Cs*CONV;
      t2 = -Cd*HBARC;
      t2f = t2*gam;
      t1 =  Cv*CONV;
      vex1 *= 2;
      vex2 *= 2;
    }

  // non-linear sigma field.
  } else {

    mSigma = 0.55;      // sigma mass (GeV)
    mOmega = 0.783;     // omega mass (GeV)
    mSigma2 = mSigma*mSigma;  // sigma mass square in  in GeV^2
    mSigmaFM2 = mSigma2/(HBARC*HBARC);  // sigma mass square in 1/fm^2

    t0  = -0.5*gs;
    //t1f = -0.5*pow2(gs/mSigma)*CONV; // This is for linear-sigma model that is not used now.
    t2 = 0.0;
    t2f = 0.0;
    t1  =  0.5*pow2(gv/mOmega)*CONV;  // Gevfm^3
    t3=0.0;
    t3f=0.0;
    gam=1.0;

    C1 = -pow2(gsp/mSigma)*CONV;
    C2 =  pow2(gvp/mOmega)*CONV;
    vex1 = C1/2;
    vex2 = C2/2;
    GS = gs*CONV;
    G2 = 2*g2*HBARC;
    G3 = 3*g3;

    G24=2.0/3.0*g2*g4/HBARC;
    G34=5.0/4.0*g3*g4/HBARC/HBARC;
    GS4=g4*gs/HBARC/HBARC;
    G4=0.5*g4/HBARC/HBARC;

    if(transportModel==2) {
      t0  = -gs;
      t1  =  pow2(gv/mOmega)*CONV;  // Gevfm^3
      vex1 = C1;
      vex2 = C2;
    }
  }

  pmu1 = mu1*mu1;
  pmu2 = mu2*mu2;

  if(firstCall) {
  cout << "# RQMDsv mode eosType= "<< eosType
      << " potentialType= "<< optPotentialType
      << " rho0= " << rho0
      << " Cs^2= " << Cs
      << " Cv^2= "<< Cv
      << " Cd^2= "<< Cd
      << " gamma= "<< gam
      << " t0= "<< t0
      << " t2= "<< t2
      << " t1= "<< t1
      << " t3= "<< t3
      << " vex1=  "<< vex1
      << " vex2=  "<< vex2
      << " mu1=  "<< mu1
      << " mu2=  "<< mu2
      <<endl;
    firstCall=false;
  }

}

void ScalarVectorPotential::setPotentialParam(EventParticle* p)
{
  //double f3=t1, f5=0.0;
  // non-linear vector potential
  //if(optPotentialType==3) {f3=t1, f5=t3;}
  // non-linear scalar potential
  //else if(optPotentialType==4) {f3=t2;f5=t1;}

  //PotentialType::setPotentialParam(p,t0,f3,f5,gam,vex1,vex2,pmu1,pmu2);

  int islam = isHyperon(p);

  //            wid,alpha,beta,gam,  c1,  gam2,  cs, cv   mus muv
  double fp[10]={1.0, 1.0, 1.0, gam, 0.0, 5./3., 1.0,1.0, 1.0,1.0};
  int idp=1;

  if(islam==0) {
    //if(t5 != 1.0) fp[4]=1.0; // three-range Skyrme potential
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
    //if(t5 == 1.0) fp[4]=0.0; // two-range Skyrme potential
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
    //if(t5 == 1.0) fp[4]=0.0; // two-range Skyrme potential
    fp[3]=gam;
    fp[5] = 5.0/3.0;
    p->setPotentialParam(fp);
    p->setPotentialId(idp);
    return;
  }

  cout << "ScalarVectorPotential::setPotentialParam ";
  cout << " something is wrong optLambaPot= "<< optLambdaPot << " islam= "<< islam<< endl;
  exit(1);

}


// scalar-vector potential --------------------------------------------------------------------------

// spotv = spot - vmoms
double ScalarVectorPotential::dSigma(double spotv, double rhos)
{
  //double sigma=(spot-vmoms)/t0;
  double sigma=spotv/t0;
  return t0*GS*pow2(1.0+G4*sigma*sigma)/
        ( mSigma2+G2*sigma+G3*sigma*sigma
          +G24*pow3(sigma)+G34*pow4(sigma)
         -2*(1.0+G4*pow2(sigma))*GS4*sigma*rhos ) ;
}

void ScalarVectorPotential::Vmd(double psq1,double psq2,double vfac1,double vfac2,double pfi1,double pfi2,double pfj1,double pfj2)
{
  // momentum-dependent scalar potential
  vmomsi = vfac1*vex1/(1.0-psq1/(pfi1*pmu1));
  vmomsj = vfac1*vex1/(1.0-psq2/(pfj1*pmu1));

  // momentum-dependent vector potential
  vmom4i = vfac2*vex2/(1.0-psq1/(pfi2*pmu2));
  vmom4j = vfac2*vex2/(1.0-psq2/(pfj2*pmu2));

}

void ScalarVectorPotential::dVdns(double rhij, double rhji)
{
  rhomij=rhij;
  rhomji=rhji;

  double pfac=pfac1*pfac2;
  s1 = optPotentialType==5 ? pfac*dSigma(spotvi,rhosi) : pfac*t0+gama1*t2*rhosgi;
  s2 = optPotentialType==5 ? pfac*dSigma(spotvj,rhosj) : pfac*t0+gama2*t2*rhosgj;

  fskys1 = fengi*s1*fj*rhomij;
  fskys2 = fengj*s2*fi*rhomji;

  pFac1=pfac3*pfac4;
  pFac2=0.0;
  pFac3=0.0;
  gam1=1.0;
  gam2=1.0;
  Vec4 Ai = facV(JBi,rhoi,rhogi,rhog2i,v2);
  Vec4 Aj = facV(JBj,rhoj,rhogj,rhog2j,v1);
  fskyv1 =(Ai * vk1)*bar1*bar2*rhomij;
  fskyv2 =(Aj * vk2)*bar1*bar2*rhomji;

  fskyi = -wG*(fskys1 + fskyv1 );
  fskyj = -wG*(fskys2 + fskyv2 );

}

void ScalarVectorPotential::dVdmd(double psq1,double psq2,double fs,double fv,
    double pmi1,double pmj1,double pmi2,double pmj2)
{
  vf1=1.0;
  vf2=1.0;
  if(optVectorPotential==1) {
    vf1 = vk1 * v2;
    vf2 = vk2 * v1;
  }
  fmomds1 = fengi*devVmd(psq2,pmi1*pmu1,vex1)*fj*rhomij*fs;
  fmomdv1 = devVmd(psq2,pmi1*pmu2,vex2)*rhomij*fv;
  fmomdi  = -wG*(fmomds1 + fmomdv1*vf1);

  fmomds2 = fengj*devVmd(psq1,pmj1*pmu1,vex1)*fi*rhomji*fs;
  fmomdv2 = devVmd(psq1,pmj2*pmu2,vex2)*rhomji*fv;
  fmomdj = -wG*(fmomds2 + fmomdv2*vf2);

  fmomei = (fengi*devVme(psq2,pmi1*pmu1,vex1)*fj*fs + fv*devVme(psq2,pmi2*pmu2,vex2)*vf1)*rhomij;
  fmomej = (fengj*devVme(psq1,pmj1*pmu1,vex1)*fi*fs + fv*devVme(psq1,pmj2*pmu2,vex2)*vf2)*rhomji;
}

void ScalarVectorPotential::dVdp()
{
  // Scalar Derivative of e^\mu/m_j term in the scalar potential.
  if(optTwoBodyDistance != 3) {
    forceri =  -fskys2*v1/p01;
    forcerj =  -fskys1*v2/p02;
  } else {
    forceri= 0.0;
    forcerj= 0.0;
  }

  // Vector Derivative of p^\mu/p0_j term in the vector potential.
  pFac1=pfac1*pfac2;
  pFac2=pfac3*pfac4;
  pFac3=pfac5;
  Vec4 A1 = facV(JBi,rhoi,rhogi,rhog2i,vk1);
  pFac3=pfac6;
  Vec4 A2 = facV(JBj,rhoj,rhogj,rhog2j,vk2);
  devV(A1,A2);
  forceri += devV1*rhomji*bar1*bar2;
  forcerj += devV2*rhomij*bar1*bar2;
  facsk= -(fskyi+fskyj)/wG;
}

void ScalarVectorPotential::dVdpm()
{
  if(optTwoBodyDistance != 3) {
    forceri = -fmomds2*v1/p01;
    forcerj = -fmomds1*v2/p02;
  } else {
    forceri=0.0;
    forcerj=0.0;
  }

  // Derivative of p_j / p_0j part.
  devV(vk1,vk2);
  forceri = fmomdv2*devV1;
  forcerj = fmomdv1*devV2;

  //forceri += fmomdv2*(optP0dev*v1 - vk2)/p01;
  //forcerj += fmomdv1*(optP0dev*v2 - vk1)/p02;

  //forceri += fmomdv2*(optP0dev*vk1 - vk2)/p01;
  //forcerj += fmomdv1*(optP0dev*vk2 - vk1)/p02;

  facmom = -(fmomdi + fmomdj)/wG;
}

} // namespace jam2
