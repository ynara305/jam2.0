#include <jam2/meanfield/RQMDw.h>
#include <jam2/meanfield/ScalarVectorPotential.h>

// Relativistic quantum molecular dynamics with sigma-omega interactions
// with momentum-dependent potentials

namespace jam2 {

using namespace std;

RQMDw::RQMDw(Pythia8::Settings* s) : MeanField(s)
{
  double widG = settings->parm("MeanField:gaussWidth");

  // interaction density 
  //facG = 1.0/pow(4.0*M_PI*widG, 1.5);
  //wG = 1.0/(4*widG);

  facG = 1.0/pow(2.0*M_PI*widG, 1.5);
  wG = 1.0/(2*widG);

  // accuracy of energy conservation
  epsEnergy = 1e-3;  // 0.1 %

  optVdot=settings->mode("MeanField:optVdot");

  double C1,C2,mu1,mu2;

  // In case only time component of vector potential is included, we do not
  // need to distinguish between kinetic momenta and canonical momenta
  if(optVectorPotential !=1) optVdot=0;

  nonlinearSigma=1;  //mstd(102)=1  ! nonlinear sigma-field is used.
  withMomDep=1;
  mSigma = 0.55;      // sigma mass (GeV)
  mOmega = 0.783;     // omega mass (GeV)
  g4 = 0.0;

//..2020/1/12 NS1'change Mnucl=0.939 and m_sigma, m_omega from 126
  //    else if (mstc(106).eq.150) then !NS1' K=210 MeV m*/m=0.8 rho0=0.168
  if(eosType==1) {     // NS1  Soft K=210 
    gs=  8.176015492139769;
    g2= 31.48609618904411;
    g3= -3.9926278758544633;
    gv= 7.724939376757135;
    gsp=0.0; gvp=0.0;
    mu1 = 1.0; mu2 = 1.0;
    withMomDep=0;

  } else if(eosType==2) {     // NS2  Hard K=380 
//..2020/1/12 NS2'change Mnucl=0.939 and m_sigma, m_omega from 117
//      else if (mstc(106).eq.151) then !NS2' K=380 MeV m*/m=0.8 rho0=0.168
    gs=  7.208208087191382;
    g2= -17.740099004710263;
    g3= 196.10588999399084;
    gv= 7.724939376757135;
    mu1 = 1.0; mu2 = 1.0;
    withMomDep=0;

    // New NS3 m*/m=0.7
  } else if(eosType==3) {     // NS3  Hard K=380 
//..2020/1/12 NS3'change Mnucl=0.939 and m_sigma, m_omega from 124
//      else if (mstc(106).eq.152) then !N23' K=380 MeV m*/m=0.7 rho0=0.168
    gs=  8.863811172925521;
    g2= 2.190629372074654;
    g3= 27.070648189159165;
    gv= 10.067747992111274;
    mu1 = 1.0; mu2 = 1.0;
    withMomDep=0;


  //if (mstc(106).eq.117) then   ! K=380 MeV m*/m=0.8 rho0=0.168
  } else if(eosType==4) {     // NS2  Hard K=380 
    gs =  7.2110738469868885;  // Scalar part g_sigma   pard(101)
    g2 = -17.888646553233265;  // Scalar part g2 (1/fm) pard(102)
    g3 =  197.63846964767876;  // Scalar part g3        pard(103)
    gv =  7.721187454057742;   // Vector part g_omega   pard(111)
    gsp=0.0; gvp=0.0;
    mu1 = 1.0; mu2 = 1.0;
    mSigma = 2.79*HBARC;     // sigma mass (GeV)
    mOmega = 3.97*HBARC;     // omega mass (GeV)
    withMomDep=0;

  //else if (mstc(106).eq.126) then   ! K=230 MeV m*/m=0.8 rho0=0.168
  } else if(eosType==5) {     // NS1  Hard K=230 
    gs=   8.18183693589196;
    g2= 31.622719664403174;
    g3= -3.797739825658857;
    gv= 7.721187454057742;
    gsp=0.0; gvp=0.0;
    mu1 = 1.0; mu2 = 1.0;
    mSigma = 2.79*HBARC;     // sigma mass (GeV)
    mOmega = 3.97*HBARC;     // omega mass (GeV)
    withMomDep=0;

  //else if (mstc(106).eq.124) then   ! K=380 MeV m*/m=0.722 rho0=0.168
  } else if(eosType==6) {     // NS3  Hard K=380 
    gs =  8.562176910455312;
    g2 = 0.4429217593591168;
    g3 = 44.70405920222624;
    gv = 9.601382237013418;
    gsp=0.0; gvp=0.0;
    mu1 = 1.0; mu2 = 1.0;
    mSigma = 2.79*HBARC;     // sigma mass (GeV)
    mOmega = 3.97*HBARC;     // omega mass (GeV)
    withMomDep=0;

  // mstc(106).eq.153
  } else if(eosType==7) {   // New NS1 Hard K=380 m*/m=0.83
//    gs=  6.448930391435675;
//    g2= -38.28483359158677;
//    g3= 342.08655094514916;
//    gv= 6.8548080859628415;
//    mSigma=2.79*HBARC;
//    mOmega=3.97*HBARC;

      gs=  6.447696268747421;
      g2= -38.001619264755874;
      g3= 339.59550852599887;
      gv= 6.858858568503943;

    gsp=0.0; gvp=0.0;
    mu1 = 1.0; mu2 = 1.0;
    withMomDep=0;

  } else if(eosType==8) {   // New NS2 Soft K=210 m*/m=0.83

    gs=  7.901789322846869;
    g2= 44.31272165001529;
    g3= 21.98880833814041;
    gv= 6.858858568503943;
    gsp=0.0; gvp=0.0;
    mu1 = 1.0; mu2 = 1.0;
    withMomDep=0;


  } else if(eosType==11) {  // K=380 MeV m*/m=0.65 SEP fit
    gs= 8.446240226735835;
    g2= 3.448041549527271;
    g3= 12.667394877335811;
    gv= 4.325553092079072;

    gsp =  4.669818923094996;
    gvp =  10.400937587189995;
    mu1 = 0.6051656607160294; // GeV
    mu2 = 1.5582862936704895; // GeV

  } else if(eosType==12) {  // K=380 MeV m*/m=0.65 Uopt->0 SEP fit
    gs=  8.69894808411858;
    g2= 3.3123747312338723;
    g3= 10.151504907242417;
    gv= 3.0515019650861466;
    gsp =  4.147458199005445;
    gvp =  10.805320816731728;
    mu1 = 0.5377037526882253;
    mu2 = 1.784448816788669;

  } else if(eosType==13) {  // K=380 MeV m*/m=0.65 gv = 0.0 SEP fit

    gs=  4.239088246904197;
    g2= -83.38119751325448;
    g3= 1039.0072461869995;
    gv= 0.0;
    gsp =  8.299983097235295;
    gvp =  11.415550399348907;
    mu1 = 1.2;
    mu2 = 1.3;

  // MD1
  } else if(eosType==14) {  // K=380 MeV m*/m=0.65 Hama1 U=0.065 at plab=3
    gs=  9.02956219322301;
    g2= 4.218109600022776;
    g3= 6.667054503467467;
    gv= 6.7402687603431675;
    gsp =  3.1855887031669394;
    gvp =  8.89548883191133;
    mu1 = 0.6411295603159707;
    mu2 = 1.841348887047491;

  // 
  } else if(eosType==15) {  // K=210 MeV m*/m=0.8 U(1.7)=0.058 U(3)=0.065
    gs=  6.551389416446063;
    g2= 90.3843215960433;
    g3= 150.1299910657155;
    g4= 0.0;
    gv= 5.757865467744954;
    gsp =  5.324936778292298;
    gvp =  5.198211638397418;
    mu1 = 0.7053607242477306;
    mu2 = 2.4912130056684383;

  } else if(eosType==16) {  // K=210 MeV m*/m=0.8 U(1.7)=0.06  U(3)=0.072
    gs=  6.566268863085936;
    g2= 90.1744321156568;
    g3= 135.0994461030604;
    g4= 0.0;
    gv= 6.462757384929801;
    gsp =  5.309336339274904;
    gvp =  4.291158097321677;
    mu1 = 0.7020113593277154;
    mu2 = 2.0301445602745294;

    // MD2   = JAM1 204 
  } else if(eosType==17) {  

    // K=380 MeV m*/m=0.65 U(1.7)=0.06  U(6)=0.04
    //gs=  9.283794468334284;
    //g2= 3.9415832307299548;
    //g3= 5.341587689277067;
    //g4= 0.0;
    //gv= 4.892556211086802;
    //gsp =  2.306806519850028;
    //gvp =  9.993473536946235;
    //mu1 = 0.4268055932276202;
    //mu2 = 2.4686246424999148;

    // K=380 MeV m*/m=0.65 U(1.7)=0.056  U(6)=0.025
      gs=  9.232564820463068;
      g2= 4.012495169512666;
      g3= 5.520251163159458;
      g4= 0.0;
      gv= 3.8884399009614437;
      gsp =  2.5019753284300768;
      gvp =  10.431917516759935;
      mu1 = 0.4896732480502826;
      mu2 = 2.4886149663285533;


  // MD3 smallest 205 in JAM1
  } else if(eosType==18) { // Hama2 K=380 MeV m*/m=0.65 U(1.6)=0.052  U(3)=0.04

    gs=  9.152032162363772;
    g2= 4.1091984804037205;
    g3= 5.886441070518481;
    g4= 0.0;
    gv= 3.568789983273673;
    gsp =  2.7878500218565083;
    gvp =  10.555177917527832;
    mu1 = 0.5668789032308825;
    mu2 = 2.3990906874902547;

  // close to NS1 parametrization
  } else if(eosType==19) { // Hama2 K=380 MeV m*/m=0.65 U(1.7)=0.062 U(4)=0.085
    //gs=  8.917093241304265;
    //g2= 4.310332275162393;
    //g3= 7.507858107772047;
    //g4= 0.0;
    //gv= 7.640246919750964;
    //gsp =  3.5101362450650186;
    //gvp =  8.159416185913237;
    //mu1 = 0.6890923623901617;
    //mu2 = 1.564516307198077;

    // K=380 MeV m*/m=0.65 U(1.7)=0.063 U(8)=0.105
    gs=  8.855455035851648;
    g2= 4.363525432945469;
    g3= 8.000210299348272;
    g4= 0.0;
    gv= 7.818179738285339;
    gsp =  3.6730776666171483;
    gvp =  8.000962678688673;
    mu1 = 0.7124877369137642;
    mu2 = 1.4835205859122296;

  // md4  Upot(infty) = 125 MeV
  } else if(eosType==20) { // K=380 MeV m*/m=0.65 U(1.7)=0.072 U(8)=0.13
    gs=  8.722412706239542;
    g2= 4.496378775852544;
    g3= 9.10233237324333;
    gv= 7.742536986570211;
    gsp =  3.991493206821452;
    gvp =  8.091790166485383;
    mu1 = 0.7654663164433807;
    mu2 = 1.435226370278356;


  // md4  Upot(infty) = 150 MeV  208 in JAM1
  } else if(eosType==21) { // K=380 MeV m*/m=0.65 U(1.7)=0.072 U(8)=0.13
    gs=  9.15499996067294;
    g2= 4.0336963331563656;
    g3= 6.143870343679022;
    g4= 0.0;
    gv= 8.620286411780615;
    gsp =  2.81531192821113;
    gvp =  7.084464641127519;
    mu1 = 0.5269385823575612;
    mu2 = 1.4655608399069255;

  // md4  Upot(infty) = 175 MeV 207 in JAM1
  } else if(eosType==22) { // K=380 MeV m*/m=0.65 U(1.7)=0.08 U(8)=0.15
    gs=  9.304844831497922;
    g2= 3.835117791722945;
    g3= 5.5579843719039905;
    g4= 0.0;
    gv= 9.199390296353101;
    gsp =  2.277148182541522;
    gvp =  6.277561735662023;

    mu1 = 0.3652112126911749;
    mu2 = 1.4301429116848454;

  //  largest Upot(infty) = 200 MeV
  } else if(eosType==23) { // K=380 MeV m*/m=0.65 U(1.7)=0.09 U(9)=0.18
    gs=  8.943540882257514;
    g2= 3.917672746640286;
    g3= 8.728327369798889;
    g4= 0.0;
    gv= 9.87112877312085;
    gsp =  3.636579180635797;
    gvp =  5.481831782950353;
    mu1 = 0.5039127940027004;
    mu2 = 0.6998101782971953;

  // MD5   209 in JAM1
  } else if(eosType==24) { // K=210 MeV m*/m=0.8 U(1.7)=0.06 U(5)=0.075
    gs=  6.594093254363855;
    g2= 89.86702225927918;
    g3= 108.51032370440886;
    g4= 0.0;
    gv= 6.3580918751932325;
    gsp =  5.277392175886514;
    gvp =  4.4396533278183234;
    mu1 = 0.6983616678464525;
    mu2 = 2.207541181176542;

  // MD4   210 in JAM1
  } else if(eosType==25) { // K=210 MeV m*/m=0.83 U(1.7)=0.06 U(5)=0.075

    gs=  4.0586271572892025;
    g2= -160.31303134630116;
    g3= 2684.378830396555;
    g4= 0.0;
    gv= 5.63247647754081;
    gsp =  5.543944245008469;
    gvp =  3.9266304417034466;
    mu1 = 0.7035464013689996;
    mu2 = 4.252421078959087;


  // MD3 211 in JAM1  extremely small omega=0.0!
  } else if(eosType==26) { // K=380 MeV m*/m=0.65 U(1.6)=0.05 U(10)=0.0
    gs=  5.4390747082116775;
    g2= -15.587197606469752;
    g3= 391.86527371750225;
    g4= 0.0;
    gv= 0.0;
    gsp =  7.711348424066859;
    gvp =  11.22012602741299;
    mu1 = 1.7024475104128964;
    mu2 = 1.898341169747188;


  } else if(eosType==27) { // K=380 MeV m*/m=0.65 U(1.7)=0.055 U(3)=0.055
       gs=  9.048354922793935;
       g2= 4.20958829512232;
       g3= 6.5093840573668444;
       g4= 0.0;
       gv= 5.900891294185006;
       gsp =  3.1245160106443817;
       gvp =  9.467708263559453;
       mu1 = 0.6357372647745847;
       mu2 = 2.005796795262418;




  } else if(eosType==30) {  // K=210 MeV m*/m=0.65 Hama1 U=0.065 at plab=3

//      cout << " RQMD.RMF eoType 15 not implemented yet"<<endl;
//      exit(1);
    gs=  9.534423675188515;
    g2= 21.44720328182629;
    g3= 0.0;
    g4= 6.379088354903728;
    gv= 6.7402687603431675;
    gsp =  3.1855887031669394;
    gvp =  8.89548883191133;
    mu1 = 0.6411295603159707;
    mu2 = 1.841348887047491;


  } else if(eosType==31) {  // K=380 MeV m*/m=0.65 Hama2
    gs=  8.947935883831653;
    g2= 4.318268151695466;
    g3= 7.139928458733808;
    g4= 0.0;
    gv= 5.782871932979855;
    gsp =  3.410696145313405;
    gvp =  9.550962907207156;
    mu1 = 0.6962348654589589;
    mu2 = 1.9446923764954853;

  } else if(eosType==32) {  // K=380 MeV m*/m=0.7 Hama1
    gs=  7.907719118613026;
    g2= -0.20607388656104766;
    g3= 50.211648401219904;
    g4= 0.0;
    gv= 6.464411434678581;
    gsp =  4.133793490330552;
    gvp =  7.8385846277021;
    mu1 = 0.6792515472633837;
    mu2 = 1.93223351514093;

  } else if(eosType==33) {  // K=380 MeV m*/m=0.65 Hama1 U=0.07 at plab=3
    gs=  8.550244792065492;
    g2= 4.690697469746118;
    g3= 10.686819254196784;
    g4= 0.0;
    gv= 7.2498054890753645;
    gsp =  4.354282269065862;
    gvp =  8.547781489920993;
    mu1 = 0.835226553464443;
    mu2 = 1.4776246824790673;

  } else if(eosType==34) {  // K=380 MeV m*/m=0.65 optmom=2

    gs=  8.677350755431574;
    g2= -1.224855619153818;
    g3= 17.343422138548316;
    g4= 0.0;
    gv= 6.6742719282765774;
    gsp =  5.17632980740873;
    gvp =  9.0911380776431;
    mu1 = 0.2;
    mu2 = 1.3;


  } else {
    cout << "MeanField:EoStype not implemented " << eosType << endl;
    exit(1);
  }

  mSigma2 = mSigma*mSigma;  // sigma mass square in  in GeV^2
  mSigmaFM2 = mSigma2/(HBARC*HBARC);  // sigma mass square in 1/fm^2

  double CONV=HBARC*HBARC*HBARC;
  C1 = -pow2(gsp/mSigma)*CONV;
  C2 =  pow2(gvp/mOmega)*CONV;

  t1  = -0.5*gs;
  t1f = -0.5*pow2(gs/mSigma)*CONV;
  t2  =  0.5*pow2(gv/mOmega)*CONV;  // Gevfm^3

  pmu1 = mu1*mu1;
  pmu2 = mu2*mu2;
  vex1 = C1/2;
  vex2 = C2/2;
  GS = gs*CONV;
  G2 = 2*g2*HBARC;
  G3 = 3*g3;

  G24=2.0/3.0*g2*g4/HBARC;
  G34=5.0/4.0*g3*g4/HBARC/HBARC;
  GS4=g4*gs/HBARC/HBARC;
  G4=0.5*g4/HBARC/HBARC;


  optMethod=2;
  if(optMethod==2) {
  // use PotentialType ===================================================
  potential = new ScalarVectorPotential(settings);
  //t0 = potential->getT0();
  //t2 = potential->getT2();
  t1 = potential->getT0();
  t2 = potential->getT3();
  //gam = potential->getGam();
  withMomDep = potential->isMomDep();
  // use PotentialType end ===================================================
  }



  cout << "# RQMD.RMF mode eosType= "<< eosType
      << " gs = " <<  gs
      << " gv = "<< gv
      << " t1 = "<< t1
      << " t1f = "<< t1f
      << " t2 = "<< t2
      << " vex1 = "<< vex1
      << " vex2 = "<< vex2
      << " pmu1 = "<< pmu1
      << " pmu2 = "<< pmu2
      <<endl;
}

void RQMDw::evolution(list<EventParticle*>& plist,double t, double dt, int step)
{
  // pre-hadrons feel potential when formation time is less than next time
  //double cutime = step * dt;

  globalTime = t;
  part.clear();
  pFree=0.0;
  for(auto i = plist.begin(); i != plist.end(); ++i) {
    if((*i)->isMeanField(t,optPotential)) {
	part.push_back(*i);
	//if(optVdot==1) part.back()->addP(part.back()->potv());
    } else {
      //eFree += (*i)->getPe() + (*i)->potv(0);
      pFree += (*i)->getP();
      pFree[0] += (*i)->potv(0);
    }
  }

  NV = part.size();
  force = new Vec4 [NV];
  forcer = new Vec4 [NV];
  rho = new double [NV];   // invariant baryon density
  rhos = new double [NV];   // scalar density
  //rhog = new double [NV];  // rho^{gam-1}
  rhom = new double *[NV]; // rho_{ij}
  vmoms = new double [NV]; // rho_{ij}
  vmom4 = new Vec4 [NV];   // momentum-dependent potential
  JB = new Vec4 [NV];      // baryon current
  Vdot = new Vec4 [NV];    // time derivative of vector potential
  vpot = new Vec4 [NV];
  for(int i=0;i<NV;i++) rhom[i] = new double [NV];

  for(int i=0;i<NV;i++) {
    force[i]= 0.0;
    forcer[i]= 0.0;
    //rhos[i]=0.0;
    //vmoms[i]=0.0;
    //vmom4[i]=0.0;
    //JB[i]=0.0;
    //Vdot[i]=0.0;
  }

  // change from kinetic to canonical momenta.
  //if(optPotentialArg>=1 || optVdot==1) {
  if(optVdot==1) {
    for(int i=0;i<NV;i++) part[i]->setFreeMass();
    //if(optPotentialArg ==1) for(int i=0;i<NV;i++) part[i]->setFreeMass();
    //else for(int i=0;i<NV;i++) part[i]->setFree();
  }

  // save vector potentials.
  for(int i=0;i<NV;i++) vpot[i] = part[i]->potv();
  qmdMatrix();
  singleParticlePotential();

  if(optMethod==1) {
    computeForce();
  } else {
    computeForceP();
  }

  /*
  if(optTwoBodyDistance == 1) {
    computeForce1();
  } else if(optTwoBodyDistance==2) {
    computeForce2();
  } else if(optTwoBodyDistance==3) {
    computeForce3();
  } else {
    computeForce4();
  }
  */


  // Note that optVdot=0,1: all momenta are canonical.
  // Canonical momenta update.
  if(optVdot == 0) {
    for(int i=0; i< NV; i++) {
      part[i]->updateByForce(force[i], forcer[i],dt);
    }

  // canonical momenta are used for the evaluation of forces,
  // after update canonical momenta, set kinetic momenta.
  } else if(optVdot == 1) {

    for(int i=0; i< NV; i++) {
      // canonical momenta update.
      part[i]->updateByForce(force[i], forcer[i],dt);
      // change to kinetic momenta.
      //part[i]->setKinetic();
    }

  // compute time derivatives of vector potential V_i^\mu numerically.
  } else if(optVdot == 2) {

    // Vdot[] is already computed in singleParticlePotential().
    for(int i=0; i< NV; i++) {
      part[i]->updateByForce(force[i] - Vdot[i]/dt, forcer[i],dt);
    }

  // compute time derivatives of vector potential analytically.
  } else if(optVdot == 3) {

    computeVdot();
    for(int i=0; i< NV; i++) {
      Vec4 f = (force[i] - Vdot[i]);
      part[i]->updateByForce(f, forcer[i],dt);
    }

  // compute time derivatives of V^\mu numerically; (V(t+dt) - V(t))/dt
  // (under construction)
  } else if(optVdot == 4 ) {

    for(int i=0;i<NV;i++) {
      //Vec4 f = force[i]*dt + part[i]->potv();
      part[i]->updateByForce(force[i], forcer[i],dt);
      //part[i]->updateByForce(force[i]*dt/2, forcer[i]*dt/2);
    }
    for(int i=0;i<NV;i++) vpot[i] = part[i]->potv();
    qmdMatrix();
    singleParticlePotential();
    for(int i=0;i<NV;i++) {
      //part[i]->updateByForce(force[i]*dt/2-Vdot[i], forcer[i]*dt/2);
      part[i]->addP(-Vdot[i]);
      part[i]->setOnShell();
      //cout << " Vdot= "<< Vdot[i];
    }

    //for(int i=0;i<NV;i++) part[i]->setKinetic();

  } else {
    cout << "RQMDw wrong option optVdot = "<< optVdot<<endl;
    exit(1);
  }

  qmdMatrix();
  singleParticlePotential();

  // Save initial total energy-momentum.
  if(step==1) computeEnergy(plist,step);

  if(optRecoverEnergy==1) RecoverEnergy(plist);
  else if(optRecoverEnergy==2) RecoverEnergy2();
  else if(optRecoverEnergy==3) RecoverEnergy3();

  if(isDebug > 1) computeEnergy(plist,step);


  /*
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  for(int i=0;i<NV;i++) {
    double meff=part[i]->getMomentum().mCalc();
    if(abs(meff-part[i]->getEffectiveMass())>1e-9) {
        cout << "RQMDw meff?  meff= "<< meff
            << " mOutEff= "<< part[i]->getEffectiveMass()
            << endl;
        exit(1);
        }
  }
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  */

  delete [] rho;
  delete [] rhos;
  delete [] vmoms;
  delete [] vmom4;
  delete [] JB;
  for(int i=0;i<NV;i++) delete [] rhom[i];
  delete [] rhom;
  delete [] force;
  delete [] forcer;
  delete [] vpot;
}

// Compute single particle potential.
void RQMDw::singleParticlePotential()
{
  // Compute sigma field.
  sigmaField();

  Vec4 vTotal=0.0;
  Vec4 vc = 0.0;
  Vec4 vmtot=0.0;
  bool optV=false;
  //bool optV=true;
  //if(optVdot==1) optV=true;
  //if(optVdot==2 || optVdot==4) optV=true;
  if(optV && optVectorPotential==1) {
    for(int i=0; i< NV; i++) {
      int ibary = part[i]->baryon()/3;
      vTotal += ibary*t2*JB[i] + vmom4[i];
      vmtot += vmom4[i];
    }
    vc = vTotal / NV;
    vc[0]=0.0;
    vmtot /= NV;
  }

  Vec4 vdott = 0.0;
  for(int i=0; i< NV; i++) {
    //if(rho[i] < 1e-15) continue;
    //Vec4 vpot0 = part[i]->potv();
    int ibar = part[i]->baryon()/3;

    // four-components of vector potential are fully included.
    if(optVectorPotential==1) {
      vmom4[i] -= vmtot;
      part[i]->setPotV( ibar*t2*JB[i] + vmom4[i] -vc);

      // compute time derivatives of vector potential V_i^\mu numerically.
      if(optVdot == 2 || optVdot==4) {
        //Vdot[i] = part[i]->potv() - vpot0;
        Vdot[i] = part[i]->potv() - vpot[i];
        vdott += Vdot[i];
      }

    // only time component of vector potential is included.
    } else if(optVectorPotential==2) {
      part[i]->setPotV(0,ibar*t2*JB[i][0] + vmom4[i][0] );

    // only time component of vector potential is included in the form
    // of V(rho_B) where rho_B is an invariant baryon density.
    } else {
      part[i]->setPotV(0,ibar*t2*rho[i] + vmom4[i][0] );
    }


  }

  if(optVdot == 2 || optVdot==4) {
    vdott /= NV;
    for(int i=0; i< NV; i++) Vdot[i] -= vdott;

  } 


}

// compute single particle potential energy.
Vec4 RQMDw::computeEnergy(list<EventParticle*>& plist,int step)
{
  pTot=0.0;   
  for(auto i = plist.begin(); i != plist.end(); ++i) {
    if(optVectorPotential==1 && optVdot==0) {
      Vec4 pk= (*i)->getP() - (*i)->potv();
      double m = (*i)->getEffectiveMass();
       pTot[0] += sqrt( m*m + pk.pAbs2());
       //pTot[1] += pk[1];
       //pTot[2] += pk[2];
       //pTot[3] += pk[3];
      pTot[1] += (*i)->getP(1);
      pTot[2] += (*i)->getP(2);
      pTot[3] += (*i)->getP(3);
    } else {
      pTot += (*i)->getP();
    }
    pTot[0] += (*i)->potv(0);
  }

  if(step==1) pTot0 = pTot;

  if(isDebug > 1) {
  double econ=abs(pTot0[0]-pTot[0])/pTot0[0]*100;
  cout << "RQMDw: step= " << step
     << " econ= " << fixed << econ << " %"
     << scientific << setw(13) << pTot[1] 
     << scientific << setw(13) << pTot[2]
     << scientific << setw(13) << pTot[3]
     <<endl;
  }

  return pTot;

}

double RQMDw::funcEnergy(list<EventParticle*>& plist,double a)
{
  qmdMatrix(a);
  singleParticlePotential();
  double etot=0.0;
  for(auto i = plist.begin(); i != plist.end(); ++i) {
    double m = (*i)->getEffectiveMass();
    double pa2 = (*i)->pAbs2();
    double e = sqrt(m*m + a*a*pa2) + (*i)->potv(0);
    etot += e;
  }
  return etot;
}

// Momenta of all particles are changed 
// to recover the total energy of the system.
void RQMDw::RecoverEnergy(list<EventParticle*>& plist)
{
  // total energy of the system 
  Vec4 ptot=pTot0;
  double etot0 = ptot[0];

  Vec4 ptotcm=0.0;
  bool boost = false;
  if(ptot.pAbs2() > 1e-9)  {
    boost = true;
    etot0 = ptot.mCalc();
    for(auto i = plist.begin(); i != plist.end(); ++i) ptotcm += (*i)->getP();
    Vec4 ptot4 = 0.0;
    for(auto i = plist.begin(); i != plist.end(); ++i) {
	(*i)->bstback(ptotcm);
	ptot4 += (*i)->getP();
    }
    cout << "RQMDw::RecoverEnergy not implemented yet  ptot4= "<< ptot4 ;
    cout  << " e= "<< etot0 << endl;
    exit(1);
  }

  double a=1.0;
  double etot = funcEnergy(plist,a);
  int itry=0;
  while(abs(etot-etot0) > epsEnergy*etot0) {
    a = etot0/etot*a;
    etot = funcEnergy(plist,a);
    if(++itry>10) break;
  }
  
  double efinal=0.0;
  for(auto i = plist.begin(); i != plist.end(); ++i) {
    (*i)->multP(a);
    (*i)->setOnShell();
    efinal += (*i)->getPe();
  }

  if(boost) {
    for(auto i = plist.begin(); i != plist.end(); ++i) (*i)->bst(ptotcm);
  }

}

double RQMDw::funcEnergy4(double a)
{
  for (int i=0;i<NV; i++) {
      double meff=part[i]->getEffectiveMass();
      double m0=part[i]->getMass();
      part[i]->setPotS(a*meff - m0);
  }

  qmdMatrix();
  singleParticlePotential();
  double etot=0.0;
  for (int i=0;i<NV; i++) {
    double m = part[i]->getEffectiveMass();
    double pa2 = part[i]->pAbs2();
    double e = sqrt(m*m + pa2) + part[i]->potv(0);
    etot += e;
  }
  return etot;
}

void RQMDw::RecoverEnergy4()
{
  // total energy of the system 
  double etot0 = pTot0[0]-pFree[0];
  cout << " pTot= "<< pTot0[0] << " efree= "<< pFree[0] <<endl;

  double a=1.0;
  int itry=0;
  double etot=0.0;
  do {
    etot = funcEnergy4(a);
    a *= etot0/etot;
    if(++itry>10) break;

  } while(abs(etot-etot0) > epsEnergy*etot0);

  //cout << "itry= "<< itry << " a= "<< a
  //    << " etot= "<< (etot-etot0)/etot0 <<endl;
  
  /*
  double efinal=0.0;
  for(auto i = plist.begin(); i != plist.end(); ++i) {
    (*i)->multP(a);
    (*i)->setOnShell();
    efinal += (*i)->getPe();
  }
  */


}

double RQMDw::funcEnergy2(double a)
{
  for (int i=0;i<NV; i++) {
    part[i]->multP(a);
    part[i]->setOnShell();
  }

  qmdMatrix();
  singleParticlePotential();
  double etot=0.0;
  for(int i = 0; i < NV; ++i) {
    //double m = part[i]->getEffectiveMass();
    //double pa2 = part[i]->pAbs2();
    //double e1 = sqrt(m*m + pa2) + part[i]->potv(0);
    double e = part[i]->getPe() + part[i]->potv(0);
    //cout << " e= " << e << " e1= "<< e1 <<endl;
    etot += e;
  }
  return etot;
}

// Momenta of particles interacting with potentials are changed 
// to recover the total energy of the system.
void RQMDw::RecoverEnergy2()
{
  // total energy of the system 
  Vec4 ptot=pTot0 - pFree;
  double etot0 = pTot[0] - pFree[0];

  Vec4 ptotcm=0.0;
  bool boost = false;
  if(ptot.pAbs2() > 1e-9)  {
    boost = true;
    // total energy in the c.m. frame.
    etot0 = ptot.mCalc();
    for(int i = 0; i < NV; ++i) ptotcm += part[i]->getP();
    for(int i = 0; i < NV; ++i) part[i]->bstback(ptotcm);
  }

  int itry=0;
  double a=1.0;
  double etot = funcEnergy2(a);
  while(abs(etot-etot0) > epsEnergy*etot0) {
    a = etot0/etot*a;
    etot = funcEnergy2(a);
    if(++itry>10) break;
  }

  if(abs(etot-etot0) > epsEnergy*etot0) {
      cout << itry << " abs(etot-etot0) = "<< abs(etot-etot0)/etot0 <<endl;
  }
  
  double efinal=0.0;
  for(int i = 0; i < NV; ++i) {
    efinal += part[i]->getPe();
  }

  if(boost) {
    ptotcm[0]=efinal;
    Vec4 pfinal=0.0;
    for(int i = 0; i < NV; ++i) {
	part[i]->bst(ptotcm,ptotcm[0]); 
	part[i]->setOnShell();
	pfinal +=part[i]->getP();
    }
    qmdMatrix();  // recalculate potentials in the computational frame.
    singleParticlePotential();
  }

}

double RQMDw::funcEnergy3(double a)
{
  qmdMatrix(a);
  singleParticlePotential();
  double etot=0.0;
  for(int i = 0; i < NV; ++i) {
    double m = part[i]->getEffectiveMass();
    double pa2 = part[i]->pAbs2();
    double e = sqrt(m*m + a*a*pa2) + part[i]->potv(0);
    etot += e;
  }
  return etot;
}

// Momenta of particles interacting with potentials are changed 
// to recover the total energy of the system.
void RQMDw::RecoverEnergy3()
{
  // total energy of the system 
  Vec4 ptot=pTot0 - pFree;
  double etot0 = pTot[0] - pFree[0];

  Vec4 ptotcm=0.0;
  bool boost = false;
  if(ptot.pAbs2() > 1e-9)  {
    boost = true;
    // total energy in the c.m. frame.
    etot0 = ptot.mCalc();
    for(int i = 0; i < NV; ++i) ptotcm += part[i]->getP();
    Vec4 ptot4 = 0.0;
    for(int i = 0; i < NV; ++i) {
	//part[i]->bstback(ptot);
	//part[i]->bstback(ptotcm,etot0);
	part[i]->bstback(ptotcm);
	ptot4 += part[i]->getP();
    }

    /*
    cout << "RQMDw::RecoverEnergy3 ptot4= "<< ptot4 ;
    cout  << " e= "<< scientific << setw(13) << etot0
	<<  " e0= "<< scientific << setw(13) << etot0<<endl;

    cout << " ptotcm= "<< scientific << setw(13) << ptotcm;
    cout << "   ptot= "<< scientific << setw(13) << ptot;
    */
  }

  int itry=0;
  double a=1.0;
  double etot = funcEnergy3(a);
  while(abs(etot-etot0) > epsEnergy*etot0) {
    a = etot0/etot*a;
    etot = funcEnergy3(a);
    if(++itry>10) break;
  }
  
  //Vec4 pfinal = 0.0;
  double efinal=0.0;
  for(int i = 0; i < NV; ++i) {
    part[i]->multP(a);
    part[i]->setOnShell();
    efinal += part[i]->getPe();
    //pfinal += part[i]->getP();
  }

  //cout << itry << " edif= "<< etot-etot0 << " pfinal = "<< pfinal;

  if(boost) {
    ptotcm[0]=efinal;
    Vec4 pfinal=0.0;
    for(int i = 0; i < NV; ++i) {
	//part[i]->bst(ptotcm,etot0); 
	//part[i]->bst(ptot); 
	part[i]->bst(ptotcm,ptotcm[0]); 
	part[i]->setOnShell();
	pfinal +=part[i]->getP();
    }

  //cout << " ptotcm= "<< ptotcm
  //     << " pfinal= " << pfinal;
  //cin.get();
  }

}

// Pre-factor from the baryon current.
Vec4 RQMDw::facV(int i, Vec4& p)
{
  Vec4 vi=0.0;
  if(rho[i]<1e-8) return vi;

  // 4-components of the vector potential is fully included.
  if(optVectorPotential==1) {
    vi = p/p[0];

  // Only time-component V^0 term is included.
  } else if(optVectorPotential==2) {
    vi[0]=1.0;

  // Only time-component V^0 term is included with the form of V(rho_B).
  } else if (optVectorPotential==3) {
    return t2/rho[i] * JB[i];
    //return (t2 + t3f*rhog[i] + t4f*pow(rho[i],gam))/rho[i] * JB[i];
  }

  return t2*vi;

  double vj = JB[i] * vi;
  //double vv = t2 + t3*rhog[i] + t4*pow(rho[i],gam-1);  // V/rho_B
  //double dv = (gam-1.0)*t3*pow(rho[i],gam-3);  // del(V/rho)/rho
  double vv = t2;  // V/rho_B
  double dv = 0.0;
  return vj*dv*JB[i] + vv*vi;

}

// non-relativistic distance.
void RQMDw::computeForce1()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = p1;
    double emi = part[i]->getEffectiveMass();
    if(optVdot <= 1) pk1 = part[i]->getPkin();
    Vec4 A1 = facV(i,pk1);
    double feng1=emi/pk1[0];
    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    Vec4 v1 = p1/p1[0];
    double s1 = dSigma(i);
    double fi = p1.mCalc()/p1[0];
    int  bar1 = part[i]->baryon()/3;

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      Vec4 pk2 = p2;
      double emj = part[j]->getEffectiveMass();
      if(optVdot <= 1) pk2 = part[j]->getPkin();
      if(optPotentialArg==1) p2 = part[j]->getPcan(optVdot);
      Vec4 v2 = p2/p2[0];
      Vec4 A2 = facV(j,pk2);
      double feng2=emj/pk2[0];
      double s2 = dSigma(j);
      double fj = p2.mCalc()/p2[0];
      int  bar2 = part[j]->baryon()/3;

      // Scalar part.
      double fsky1 = feng1 * fj * s1 * rhom[i][j];
      double fsky2 = feng2 * fi * s2 * rhom[j][i];

      // Vector part.
      double fsky3 = (A1 * v2)*rhom[i][j]*bar1*bar2;
      double fsky4 = (A2 * v1)*rhom[j][i]*bar1*bar2;
      double fsky = -wG*(fsky1 + fsky2 + fsky3 + fsky4);

      force[i] += -2*fsky*(r1 - r2);
      force[j] += -2*fsky*(r2 - r1);

      if(optDerivative) {
      // Derivative of m_j/p^0_j in the scalar density.
      forcer[i] += -fsky2 / p1[0] * v1 * optP0dev;
      forcer[j] += -fsky1 / p2[0] * v2 * optP0dev;

      // Derivative of p^\mu_j/p^0_j in the vector density.
        double fsgam2i = optP0dev*rhom[j][i]*dot3(A2,v1)/p1[0];
        double fsgam2j = optP0dev*rhom[i][j]*dot3(A1,v2)/p2[0];
        forcer[i] += fsgam2i*v1 - A2/p1[0]*rhom[j][i];
        forcer[i] += fsgam2j*v2 - A1/p2[0]*rhom[i][j];
      }

      if(!withMomDep) continue;

      double vf=2.0;
      double vf1 = 1.0;
      double vf2 = 1.0;
      if(optVectorPotential==1) {
        vf = (pk1 * v2)/pk1[0] + (pk2 * v1)/pk2[0];
        vf1 = (pk1 * v2)/pk1[0];
	vf2 = (pk2 * v1)/pk2[0];
      }

      Vec4 dp = p1 - p2;
      double psq = dp.pAbs2() - (1-optPV) * dp[0]*dp[0];
      /*
      //double v01 = vmom4[i][0]-vmom4[j][0];
      //double psq = dp.pAbs2() - (1-optPV) * (dp[0]+v01)*(dp[0]+v01);
      if(psq < 0.0) {
	  cout << " psq < 0 " << psq <<endl;
	  exit(1);
      }
      */

      //psq = max(0.0,psq);

      double fac1 = 1.0 + psq/pmu1;
      double fac2 = 1.0 + psq/pmu2;
      double vm1 = vex1/fac1;
      double vm2 = vex2/fac2;
      double dvm1 = -vex1/(pmu1*fac1*fac1);
      double dvm2 = -vex2/(pmu2*fac2*fac2);

     //cout << " dp= "<< dp.pAbs2() << " psq0= "<< psq<<endl;

      /*
      // Scalar part
      double feng12 = feng1*fj*rhom[i][j] + feng2*fi*rhom[j][i];
      // Vector part
      double fengij = vf1*rhom[i][j] + vf2*rhom[j][i];
      // r derivative term.
      double fmomd0 = -wG*feng12*vex1/fac1 - wG*fengij*vex2/fac2;
      // p derivative term.
      double fmome0  = -feng12*vex1/(pmu1*fac1*fac1)
	              -fengij*vex2/(pmu2*fac2*fac2);
      */
		      

      double sf = feng1*fj + feng2*fi;
      double fmomd = -wG*( sf * vm1  + vf * vm2  )*rhom[i][j];
      double fmome =     ( sf * dvm1 + vf * dvm2 )*rhom[i][j];

      Vec4 dp2pi =  dp -optP0dev*(1-optPV)*dp[0]*v1;
      Vec4 dp2pj = -dp +optP0dev*(1-optPV)*dp[0]*v2;

      force[i]  += -2*fmomd*(r1 - r2);
      force[j]  += -2*fmomd*(r2 - r1);
      forcer[i] +=  2*fmome*dp2pi;
      forcer[j] +=  2*fmome*dp2pj;

      if(optDerivative) {

      // Derivative of m_j/p^0_j in the scalar density.
      double fs1 = -vm1*feng2 * fi / p1[0] * rhom[j][i] * optP0dev;
      double fs2 = -vm1*feng1 * fj / p2[0] * rhom[i][j] * optP0dev;
      forcer[i] += fs1 * v1;
      forcer[j] += fs2 * v2;

      // Derivative of p_j/p^0_j in the vector potential.
      double fm1 = vm2*rhom[j][i]/p1[0];
      double fm2 = vm2*rhom[i][j]/p2[0];
      vf1 = dot3(pk1, v2)/pk1[0];
      vf2 = dot3(pk2, v1)/pk2[0];
      forcer[i] += fm1 * ( -pk2/pk2[0] + vf2 * v1 * optP0dev );
      forcer[j] += fm2 * ( -pk1/pk1[0] + vf1 * v2 * optP0dev );
      }

    }
  }

}

// two-body distance is defined by the c.m. frame of two-particles.
void RQMDw::computeForce2()
{
  for(int i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double em1 = part[i]->getEffectiveMass();
    Vec4 pk1 = p1;
    if(optVdot <= 1) pk1 = part[i]->getPkin();
    double feng1=em1/pk1[0];
    Vec4 A1 = facV(i,pk1);

    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    double emi = p1.mCalc();
    Vec4 v1 = p1/p1[0];
    double fi = emi/p1[0];

    double s1 = dSigma(i);
    int  b1 = part[i]->baryon()/3;

    for(int j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double em2 = part[j]->getEffectiveMass();
      Vec4 pk2 = p2;
      if(optVdot <= 1) pk2 = part[j]->getPkin();
      double feng2=em2/pk2[0];
      Vec4 A2 = facV(j,pk2);

     if(optPotentialArg==1) p2 = part[j]->getPcan(optVdot);
      double emj = p2.mCalc();
      Vec4 v2 = p2/p2[0];
      double fj = emj/p2[0];
      double s2=dSigma(j);
      int  b2 = part[j]->baryon()/3;

      // Scalar part:
      double fsky1 = feng1*fj*s1*rhom[i][j];
      double fsky2 = feng2*fi*s2*rhom[j][i];

      // Vector part:
      double fsky3 = (A1 * v2)*b1*b2*rhom[i][j];
      double fsky4 = (A2 * v1)*b2*b1*rhom[j][i];
      double fsky = -wG*(fsky1 + fsky2 + fsky3 + fsky4 );

      Vec4 pcm = p1 + p2;
      Vec4 bbi = pcm/pcm[0] - optP0dev*v1;
      Vec4 bbj = pcm/pcm[0] - optP0dev*v2;

      Vec4 dr  = r1 - r2;dr[0]=0.0;
      double s = pcm.m2Calc();
      double rbij = -dr * pcm/s;
      Vec4 dr2ri = dr + rbij*pcm;
      Vec4 dr2rj = -dr2ri;
      Vec4 dr2pi = (dr + pcm[0]*rbij*bbi)*rbij;
      Vec4 dr2pj = (dr + pcm[0]*rbij*bbj)*rbij;

      force[i]  += -2*fsky*dr2ri;
      force[j]  += -2*fsky*dr2rj;
      forcer[i] +=  2*fsky*dr2pi;
      forcer[j] +=  2*fsky*dr2pj;

      if(optDerivative) {
      // Derivative of gamma_ij in front of Gaussian.
        double facsk = fsky1 + fsky2 + fsky3 + fsky4;
	double fsgam1 = facsk / s;
	double fsgam2 = optP0dev*facsk*(1.0/pcm[0] - pcm[0]/s);
	forcer[i] += fsgam1*pcm + fsgam2*v1;
	forcer[j] += fsgam1*pcm + fsgam2*v2;

        // Derivative of m_j/p^0_j in the scalar density.
        double fs1 = -fsky2 / p1[0] * optP0dev;
        double fs2 = -fsky1 / p2[0] * optP0dev;
        forcer[i] += fs1 * v1;
        forcer[j] += fs2 * v2;

      // Derivative of p^\mu_j/p^0_j in the vector density.
      double fsgam2i = optP0dev*rhom[j][i]*dot3(A2,v1)/p1[0];
      double fsgam3i = -rhom[j][i]/p1[0];
      double fsgam2j = optP0dev*rhom[i][j]*dot3(A1,v2)/p2[0];
      double fsgam3j = -rhom[i][j]/p2[0];
      forcer[i] += fsgam2i*v1 + A2*fsgam3i;
      forcer[i] += fsgam2j*v2 + A1*fsgam3j;
      }

      if(!withMomDep) continue;

      double vf1=1.0;
      double vf2=1.0;
      if(optVectorPotential==1) {
        vf1 = (pk1 * v2)/pk1[0];
        vf2 = (pk2 * v1)/pk2[0];
      }

      Vec4 dp = p1 - p2;
      double psq = dp.m2Calc() - optPV * pow2(dp * pcm)/s;
      double fac1 = 1.0 - psq/pmu1;
      double fac2 = 1.0 - psq/pmu2;
      double pma = pow2(emi*emi - emj*emj)/s;

      // Scalar part.
      double fengs = feng1*fj*rhom[i][j] + feng2*fi*rhom[j][i];

      // Vector part.
      double fengv = vf1*rhom[i][j] + vf2*rhom[j][i];

      // p derivative term.
      double fmome = -fengs*vex1/(pmu1*fac1*fac1)-fengv*vex2/(pmu2*fac2*fac2);

      // r derivative term.
      double fmomd = -wG*fengs*vex1/fac1 -wG*fengv*vex2/fac2;

      Vec4 dp2pi =  dp - optP0dev*dp[0]*v1 + optPV*pcm[0]/s*pma*bbi;
      Vec4 dp2pj = -dp + optP0dev*dp[0]*v2 + optPV*pcm[0]/s*pma*bbj;

      force[i]  += -2*fmomd*dr2ri;
      force[j]  += -2*fmomd*dr2rj;
      forcer[i] +=  2*fmomd*dr2pi + 2*fmome*dp2pi;
      forcer[j] +=  2*fmomd*dr2pj + 2*fmome*dp2pj;

      if(optDerivative) {
	// derivative of gamma_ij
	double facij=fengs*vex1/fac1 + fengv*vex2/fac2;
	double fmgam1 = facij/s;
	double fmgam2 = optP0dev*facij*(1.0/pcm[0] - pcm[0]/s);
        forcer[i] += pcm*fmgam1 + v1*fmgam2;
        forcer[j] += pcm*fmgam1 + v2*fmgam2;

        // Derivative of m_j/p^0_j in the scalar density.
        double fn1 = -vex1/fac1*feng2 * fi / p1[0] * rhom[j][i] * optP0dev;
        double fn2 = -vex1/fac1*feng1 * fj / p2[0] * rhom[i][j] * optP0dev;
        forcer[i] += fn1 * v1;
        forcer[j] += fn2 * v2;

        // Derivative of p_j/p^0_j in the vector density.
        double fm1 = vex2/fac2*rhom[j][i]/p1[0];
        double fm2 = vex2/fac2*rhom[i][j]/p2[0];
        vf1 = dot3(pk1, v2)/pk1[0];
        vf2 = dot3(pk2, v1)/pk2[0];
        forcer[i] += -fm1 * pk2/pk2[0] + fm1*vf2 * v1 * optP0dev;
        forcer[j] += -fm2 * pk1/pk1[0] + fm2*vf1 * v2 * optP0dev;

      }

    }
  }

}

// two-body distance is defined by the rest frame of particle.
void RQMDw::computeForce3()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = p1;
    if(optVdot <= 1) pk1 = part[i]->getPkin();
    double feng1=part[i]->getEffectiveMass()/pk1[0];
    Vec4 Ai = facV(i,pk1);

    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    double emi = p1.mCalc();
    Vec4 vi = p1/emi;
    int  bar1= part[i]->baryon()/3;
    double s1=dSigma(i);

    for(auto j=i+1; j< NV; j++) {

      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      Vec4 pk2 = p2;
      if(optVdot <= 1) pk2 = part[j]->getPkin();
      Vec4 Aj = facV(j,pk2);
      double feng2=part[j]->getEffectiveMass()/pk2[0];

      if(optPotentialArg==1) p2 = part[j]->getPcan(optVdot);
      double emj = p2.mCalc();
      Vec4 vj = p2/emj;

      int  bar2 = part[j]->baryon()/3;
      double s2=dSigma(j);

      double fskyi = -wG*(feng1*s1 + bar1*bar2*(Ai * vj))*rhom[i][j]; 
      double fskyj = -wG*(feng2*s2 + bar1*bar2*(Aj * vi))*rhom[j][i]; 

      Vec4 dr = r1 - r2;dr[0]=0.0;
      double rbi = dot3(dr,vi);
      double rbj = dot3(dr,vj);
      Vec4 dr2ri =  2*(dr + rbj*vj);  // R^2_{ij}/dp_i
      Vec4 dr2rj =  2*(dr + rbi*vi);  // R^2_{ji}/dp_i
      Vec4 dr2pi =  2*dr*rbi/emi;     // R^2_{ji}/dp_i
      Vec4 dr2pj =  2*dr*rbj/emj;     // R^2_{ij}/dp_j

      force[i]  += -fskyi*dr2ri - fskyj*dr2rj;
      force[j]  +=  fskyi*dr2ri + fskyj*dr2rj;
      forcer[i] +=  fskyj*dr2pi;
      forcer[j] +=  fskyi*dr2pj;

      if(optDerivative) {
      // Derivative of p^\mu/m_j term in the vector potential.
      forcer[i] += (optP0dev*Aj[0]*vi/p1[0] - Aj/emi)*rhom[j][i];
      forcer[j] += (optP0dev*Ai[0]*vj/p2[0] - Ai/emj)*rhom[i][j];
      }

      if(!withMomDep) continue;

      double vf1=vj[0];
      double vf2=vi[0];
      if(optVectorPotential==1) {
        vf1 = (pk1 * vj) / pk1[0];
        vf2 = (pk2 * vi) / pk2[0];
      }

      if(optMDarg3==1) distanceP1(p1,p2);
      else if(optMDarg3==2) distanceP2(p1,p2);
      else distanceP3(p1,p2);

      double fac1=1.0-psq2/pmu1;
      double fac2=1.0-psq2/pmu2;
      double fmomdi = -wG*( feng1*vex1/fac1 + vf1*vex2/fac2 )*rhom[i][j];
      double fmomei =    -( feng1*vex1/(pmu1*fac1*fac1)
	                    + vf1*vex2/(pmu2*fac2*fac2) )*rhom[i][j];

      fac1=1.0-psq1/pmu1;
      fac2=1.0-psq1/pmu2;
      double fmomdj = -wG*( feng2*vex1/fac1 + vf2*vex2/fac2 )*rhom[j][i];
      double fmomej =    -( feng2*vex1/(pmu1*fac1*fac1)
	                    + vf2*vex2/(pmu2*fac2*fac2) )*rhom[j][i];

      /*
      double fmomdi = -wG*( feng1*devVmd(psq2,pmu1,vex1)
    	                    + vf1*devVmd(psq2,pmu2,vex2) )*rhom[i][j];
      double fmomei =     ( feng1*devVme(psq2,pmu1,vex1)
	                    + vf1*devVme(psq2,pmu2,vex2) )*rhom[i][j];
      double fmomdj = -wG*( feng2*devVmd(psq1,pmu1,vex1)
 	                    + vf2*devVmd(psq1,pmu2,vex2) )*rhom[j][i];
      double fmomej =     ( feng2*devVme(psq1,pmu1,vex1)
	                    + vf2*devVme(psq1,pmu2,vex2) )*rhom[j][i];
      */

      force[i]  += -fmomdi*dr2ri   - fmomdj*dr2rj;
      force[j]  +=  fmomdi*dr2ri   + fmomdj*dr2rj;
      forcer[i] +=  fmomei*dp2ijpi + fmomej*dp2jipi + fmomdj*dr2pi;
      forcer[j] +=  fmomej*dp2jipj + fmomei*dp2ijpj + fmomdi*dr2pj;

      if(optDerivative) {
      // Derivative of p_j / m_j part.
      double facm1 = vex2/(1.0-psq1/pmu2)*rhom[j][i]/emi;
      double facm2 = vex2/(1.0-psq2/pmu2)*rhom[i][j]/emj;
      forcer[i] += facm1*( optP0dev*p1/p1[0] - pk2/pk2[0] );
      forcer[j] += facm2*( optP0dev*p2/p2[0] - pk1/pk1[0] );
      }

    }
  }

}

// Test version.
// two-body distance is defined by the c.m. frame of two-particles.
// gamma factor and m_j/p_j^0 is ignored in the scalar density.
void RQMDw::computeForce4()
{
  for(int i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double em1 = part[i]->getEffectiveMass();
    Vec4 pk1 = p1;
    if(optVdot <= 1) pk1 = part[i]->getPkin();
    double feng1=em1/pk1[0];
    Vec4 A1 = facV(i,pk1);

    if(optPotentialArg==1) p1 = part[i]->getPcan(optVdot);
    double emi = p1.mCalc();
    Vec4 v1 = p1/emi;
    //double fi = emi/p1[0];
    double fi=1.0;

    double s1 = dSigma(i);
    int  b1 = part[i]->baryon()/3;

    for(int j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double em2 = part[j]->getEffectiveMass();
      Vec4 pk2 = p2;
      if(optVdot <= 1) pk2 = part[j]->getPkin();
      double feng2=em2/pk2[0];
      Vec4 A2 = facV(j,pk2);

     if(optPotentialArg==1) p2 = part[j]->getPcan(optVdot);
      double emj = p2.mCalc();
      Vec4 v2 = p2/emj;
      double s2=dSigma(j);
      int  b2 = part[j]->baryon()/3;
      //double fj = emj/p2[0];
      double fj=1.0;

      // Scalar part:
      double fsky1 = feng1*fj*s1*rhom[i][j];
      double fsky2 = feng2*fi*s2*rhom[j][i];

      // Vector part:
      double fengi = (A1 * v2)*b1*b2;
      double fengj = (A2 * v1)*b2*b1;
      double fsky3 = fengi*rhom[i][j];
      double fsky4 = fengj*rhom[j][i];
      double fsky = -wG*(fsky1 + fsky2 + fsky3 + fsky4 );

      Vec4 pcm = p1 + p2;
      Vec4 bbi = pcm/pcm[0] - optP0dev*v1;
      Vec4 bbj = pcm/pcm[0] - optP0dev*v2;

      Vec4 dr  = r1 - r2;dr[0]=0.0;
      double s = pcm.m2Calc();
      double rbij = -dr * pcm/s;
      Vec4 dr2ri = dr + rbij*pcm;
      Vec4 dr2rj = -dr2ri;
      Vec4 dr2pi = (dr + pcm[0]*rbij*bbi)*rbij;
      Vec4 dr2pj = (dr + pcm[0]*rbij*bbj)*rbij;

      force[i]  += -2*fsky*dr2ri;
      force[j]  += -2*fsky*dr2rj;
      forcer[i] +=  2*fsky*dr2pi;
      forcer[j] +=  2*fsky*dr2pj;

      if(optDerivative) {

      // Derivative of p^\mu_j/m_j in the vector density.
      forcer[i] += (optP0dev*A2[0]*v1/p1[0] - A2/emi)*rhom[j][i];
      forcer[j] += (optP0dev*A1[0]*v2/p2[0] - A1/emj)*rhom[i][j];

      // Derivative of p^\mu_j/p^0_j in the vector density.
      //double fsgam2i = optP0dev*rhom[j][i]*dot3(A2,v1)/p1[0];
      //double fsgam3i = -rhom[j][i]/p1[0];
      //double fsgam2j = optP0dev*rhom[i][j]*dot3(A1,v2)/p2[0];
      //double fsgam3j = -rhom[i][j]/p2[0];
      //forcer[i] += fsgam2i*v1 + A2*fsgam3i;
      //forcer[i] += fsgam2j*v2 + A1*fsgam3j;
      }

      if(!withMomDep) continue;

      double vf1=v2[0];
      double vf2=v1[0];
      if(optVectorPotential==1) {
        vf1 = (pk1 * v2)/pk1[0];
        vf2 = (pk2 * v1)/pk2[0];
      }

      Vec4 dp = p1 - p2;
      double psq = dp.m2Calc() - optPV * pow2(dp * pcm)/s;
      double fac1 = 1.0 - psq/pmu1;
      double fac2 = 1.0 - psq/pmu2;
      double pma = pow2(emi*emi - emj*emj)/s;

      // Scalar part.
      double fengs = feng1*fj*rhom[i][j] + feng2*fi*rhom[j][i];

      // Vector part.
      double fengv = vf1*rhom[i][j] + vf2*rhom[j][i];

      // p derivative term.
      double fmome = -fengs*vex1/(pmu1*fac1*fac1)-fengv*vex2/(pmu2*fac2*fac2);

      // r derivative term.
      double fmomd = -wG*fengs*vex1/fac1 -wG*fengv*vex2/fac2;

      Vec4 dp2pi =  dp - optP0dev*dp[0]*v1 + optPV*pcm[0]/s*pma*bbi;
      Vec4 dp2pj = -dp + optP0dev*dp[0]*v2 + optPV*pcm[0]/s*pma*bbj;

      force[i]  += -2*fmomd*dr2ri;
      force[j]  += -2*fmomd*dr2rj;
      forcer[i] +=  2*fmomd*dr2pi + 2*fmome*dp2pi;
      forcer[j] +=  2*fmomd*dr2pj + 2*fmome*dp2pj;

      if(optDerivative) {
      // Derivative of p_j / m_j part.
      double facm1 = vex2/(1.0-psq1/pmu2)*rhom[j][i]/emi;
      double facm2 = vex2/(1.0-psq2/pmu2)*rhom[i][j]/emj;
      forcer[i] += facm1*( optP0dev*p1/p1[0] - pk2/pk2[0] );
      forcer[j] += facm2*( optP0dev*p2/p2[0] - pk1/pk1[0] );
      }

    }
  }

}


// non-relativistic two-body distance
// (distance at the global cm frame where the total momentum is zero)
void RQMDw::qmdMatrix1(double a)
{
  for(int i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double m1 = part[i]->getEffectiveMass();
    if(optPotentialArg==1) {
      p1 = part[i]->getPcan(optVdot);
      m1 = part[i]->getMass();
    }
    p1 *= a; p1[0]= sqrt(m1*m1+p1.pAbs2());
    double fi = m1/p1[0];
    int bi = part[i]->baryon()/3;
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;

    for(int j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double m2 = part[j]->getEffectiveMass();
      if(optPotentialArg==1) {
        p2 = part[j]->getPcan(optVdot);
        m2 = part[j]->getMass();
      }
      p2 *=a; p2[0]= sqrt(m2*m2+p2.pAbs2());
      double fj=m2/p2[0];
      int bj = part[j]->baryon()/3;
      double qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;

      double drsq = (r1 - r2).pAbs2();
      rhom[i][j] = facG * exp(-drsq*wG)*qfac1*qfac2;
      rhom[j][i] = rhom[i][j];
      rhos[i] += rhom[i][j]*fj;
      rhos[j] += rhom[j][i]*fi;

      JB[i] += rhom[i][j]*p2/p2[0]*bj;
      JB[j] += rhom[j][i]*p1/p1[0]*bi;

      if(!withMomDep) continue;
      double ps = (p1 - p2).pAbs2() - (1-optPV) * pow2((p1-p2)[0]);
      ps = max(0.0,ps);

      double pmom2s = vex1/(1.0+ps/pmu1);
      double pmom2v = vex2/(1.0+ps/pmu2);
      vmoms[i] += pmom2s*rhom[i][j]*fj;
      vmoms[j] += pmom2s*rhom[j][i]*fi;
      vmom4[i] += pmom2v*rhom[i][j]*p2/p2[0];
      vmom4[j] += pmom2v*rhom[j][i]*p1/p1[0];
    }
  }

  for(int i=0;i<NV;i++) {
    rho[i] = sqrt(max(0.0, JB[i].m2Calc()));
    //rhog[i] = pow(rho[i],gam-1);
  }

}

// two-body distance at their center-of-momentum frame.
void RQMDw::qmdMatrix2(double a)
{
  for(int i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    int bi = part[i]->baryon()/3;
    double m1 = part[i]->getEffectiveMass();
    if(optPotentialArg==1) {
      p1 = part[i]->getPcan(optVdot);
      m1 = part[i]->getMass();
    }
    p1 *= a; p1[0]= sqrt(m1*m1+p1.pAbs2());
    double fi = m1/p1[0];
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;

    for(int j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      int bj = part[j]->baryon()/3;
      double m2 = part[j]->getEffectiveMass();
      if(optPotentialArg==1) {
	p2 = part[j]->getPcan(optVdot);
        m2 = part[j]->getMass();
      }
      p2 *=a; p2[0]= sqrt(m2*m2+p2.pAbs2());
      double fj=m2/p2[0];
      double qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;

      Vec4 dr = r1 - r2; dr[0] = 0.0;  // just in case;
      Vec4 P  = p1 + p2;
      double s = P.m2Calc();
      double drcmsq = dr.m2Calc() - pow2(dr * P)/s;
      rhom[i][j] = P[0]/sqrt(s) * facG * exp(drcmsq*wG)*qfac1*qfac2;
      rhom[j][i] = rhom[i][j];
      rhos[i] += rhom[i][j]*fj;
      rhos[j] += rhom[j][i]*fi;
      JB[i] += rhom[i][j]*p2/p2[0]*bj;
      JB[j] += rhom[j][i]*p1/p1[0]*bi;

      if(!withMomDep) continue;
      Vec4 dp = p1 - p2;
      double ps = dp.m2Calc() - optPV * pow2(dp * P)/s;
      double pmom2s = vex1/(1.0-ps/pmu1);
      double pmom2v = vex2/(1.0-ps/pmu2);
      vmoms[i] += pmom2s*rhom[i][j]*fj;
      vmoms[j] += pmom2s*rhom[j][i]*fi;
      vmom4[i] += pmom2v*rhom[i][j]*p2/p2[0];
      vmom4[j] += pmom2v*rhom[j][i]*p1/p1[0];
    }
  }

  for(int i=0;i<NV;i++) {
    rho[i] = sqrt(max(0.0, JB[i].m2Calc()));
    //rhog[i] = pow(rho[i],gam-1);
  }

}

// two-body distance is computed from the rest-frame of particle j.
void RQMDw::qmdMatrix3(double a)
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    int bi = part[i]->baryon()/3;
    double m1 = part[i]->getEffectiveMass();
    if(optPotentialArg==1) {
      p1 = part[i]->getPcan(optVdot);
      m1 = part[i]->getMass();
    }

    p1 *= a; p1[0]= sqrt(m1*m1+p1.pAbs2());
    //double m1 = p1.mCalc();
    Vec4 v1 = p1/m1;
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      int bj = part[j]->baryon()/3;
      double m2 = part[j]->getEffectiveMass();
      if(optPotentialArg==1) {
	 p2 = part[j]->getPcan(optVdot);
         m2 = part[j]->getMass();
      }

      p2 *=a; p2[0]= sqrt(m2*m2+p2.pAbs2());
      //double m2 = p2.mCalc();
      Vec4 v2 = p2/m2;
      double qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;

      Vec4 dr = r1 - r2; dr[0] = 0.0;  // just in case;
      double drsq1 = dr.m2Calc() - pow2(dr * v2);
      rhom[i][j] = facG * exp(drsq1*wG)*qfac1*qfac2;

      double drsq2 = dr.m2Calc() - pow2(dr * v1);
      rhom[j][i] = facG * exp(drsq2*wG)*qfac1*qfac2;

      rhos[i] += rhom[i][j];
      rhos[j] += rhom[j][i];

      JB[i] += rhom[i][j]*v2*bj;
      JB[j] += rhom[j][i]*v1*bi;

      if(!withMomDep) continue;

      Vec4 dp = p1 - p2;
      double dpsq = dp.m2Calc();
      double ps1 = dpsq;
      double ps2 = dpsq;

      if(optPV !=0 ) {
      if(optMDarg3==1)  {
        dpsq = dp.pAbs2();
	ps1 = -dpsq; ps2= -dpsq;
      } else if(optMDarg3==2) {
        double s = (p1+p2).m2Calc();
        dpsq += - optPV * pow2(dp * (p1+p2))/s;
	ps1 = dpsq; ps2= dpsq;
      } else {
        ps1 = dpsq - optPV * pow2(dp * v1);
        ps2 = dpsq - optPV * pow2(dp * v2);
      }
      }

      vmoms[i] += vex1/(1.0-ps2/pmu1)*rhom[i][j];
      vmoms[j] += vex1/(1.0-ps1/pmu1)*rhom[j][i];
      vmom4[i] += vex2/(1.0-ps2/pmu2)*rhom[i][j]*v2;
      vmom4[j] += vex2/(1.0-ps1/pmu2)*rhom[j][i]*v1;
    }
  }

  for(int i=0;i<NV;i++) {
    rho[i] = sqrt(max(0.0, JB[i].m2Calc()));
    //rhog[i] = pow(rho[i],gam-1);
  }

}

// two-body distance at their center-of-momentum frame.
void RQMDw::qmdMatrix4(double a)
{
  for(int i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    int bi = part[i]->baryon()/3;
    double m1 = part[i]->getEffectiveMass();
    if(optPotentialArg==1) {
      p1 = part[i]->getPcan(optVdot);
      m1 = part[i]->getMass();
    }
    p1 *= a; p1[0]= sqrt(m1*m1+p1.pAbs2());
    for(int j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      int bj = part[j]->baryon()/3;
      double m2 = part[j]->getEffectiveMass();
      if(optPotentialArg==1) {
	p2 = part[j]->getPcan(optVdot);
        m2 = part[j]->getMass();
      }
      p2 *=a; p2[0]= sqrt(m2*m2+p2.pAbs2());

      Vec4 dr = r1 - r2; dr[0] = 0.0;  // just in case;
      Vec4 P  = p1 + p2;
      double s = P.m2Calc();
      double drcmsq = dr.m2Calc() - pow2(dr * P)/s;
      rhom[i][j] = facG * exp(drcmsq*wG);
      rhom[j][i] = rhom[i][j];
      rhos[i] += rhom[i][j];
      rhos[j] += rhom[j][i];
      JB[i] += rhom[i][j]*p2/m2*bj;
      JB[j] += rhom[j][i]*p1/m1*bi;

      if(!withMomDep) continue;
      Vec4 dp = p1 - p2;
      double ps = dp.m2Calc() - optPV * pow2(dp * P)/s;
      double pmom2s = vex1/(1.0-ps/pmu1);
      double pmom2v = vex2/(1.0-ps/pmu2);
      vmoms[i] += pmom2s*rhom[i][j];
      vmoms[j] += pmom2s*rhom[j][i];
      vmom4[i] += pmom2v*rhom[i][j]*p2/m2;
      vmom4[j] += pmom2v*rhom[j][i]*p1/m1;
    }
  }

  for(int i=0;i<NV;i++) {
    rho[i] = sqrt(max(0.0, JB[i].m2Calc()));
    //rhog[i] = pow(rho[i],gam-1);
  }

}


//*********************************************************************
//...Compute time-derivatives of the vector potential.
void RQMDw::computeVdot()
{
  bool opt=true;
  //opt=false;

  for(int i=0;i<NV; i++) Vdot[i]=0.0;

  for(int i=0;i<NV; i++) {
    double vvi = t2; // + t3*rhog[i];  // V_i/rho_i
    double dvi=0.0;
    // del(V_i/rho_i)/rho_i
    //if(abs(rho[i]) > 1e-7) dvi = (gam-1.0)*t3*pow(rho[i],gam-3.0);

    Vec4 Bi = dvi*JB[i];
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double mi = part[i]->getEffectiveMass();
    double e1 = sqrt(mi*mi + p1.pAbs2());
    Vec4 vk1 = p1/e1;
    double m1 = p1.mCalc();
    if(optPotentialArg==1) {
      p1 = part[i]->getPcan(optVdot);
      m1 = part[i]->getMass();
    }

    double emi=p1[0];
    if(optTwoBodyDistance==3)  emi = m1;
    Vec4 vi = p1 / emi;
    Vec4 b1 = p1 / p1[0];

    for(int j=i+1;j<NV; j++) {

      double vvj = t2; // + t3*rhog[j];
      double dvj=0.0;
      //if(abs(rho[j]) > 1e-7) dvj = (gam-1.0)*t3*pow(rho[j],gam-3.0);
      Vec4 Bj = dvj*JB[j];
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double mj = part[i]->getEffectiveMass();
      double e2 = sqrt(mj*mj + p2.pAbs2());
      Vec4 vk2 = p2/e2;
      double m2 = p2.mCalc();
      if(optPotentialArg==1) {
        p2 = part[j]->getPcan(optVdot);
	m2 = part[j]->getMass();
      }

      double emj=p2[0];
      if(optTwoBodyDistance==3)  emj = m2;
      Vec4 vj = p2 / emj;
      Vec4 b2 = p2 / p2[0];
   
      Vec4 Aj = (vj * JB[i]) * Bi;
      Vec4 Ai = (vi * JB[j]) * Bj;

      // Compute derivatives

      // Non-relativistic.
      Vec4 rk = r1 - r2;
      Vec4 dgami=0.0;
      Vec4 dgamj=0.0;
      Vec4 dr2ri   = 2 * rk;
      Vec4 dr2rj   = -dr2ri;
      Vec4 drji2ri = dr2ri;
      Vec4 drji2rj = dr2rj;
      Vec4 dr2pi = 0.0;
      Vec4 dr2pj = 0.0;
      Vec4 drji2pi=0.0;
      Vec4 drji2pj=0.0;

      // two-body c.m. frame.
      if(optTwoBodyDistance==2) {

      Vec4 pcm = p1 + p2;
      Vec4 bbi = pcm/pcm[0] - optP0dev*b1;
      Vec4 bbj = pcm/pcm[0] - optP0dev*b2;
      double s = pcm.m2Calc();
      double rbij = dot3(rk,pcm)/s;
      dr2ri =  2*(rk + rbij*pcm);
      dr2rj = -dr2ri;
      dr2pi = 2*(rk + pcm[0]*rbij*bbi)*rbij;
      dr2pj = 2*(rk + pcm[0]*rbij*bbj)*rbij;
      drji2ri =  dr2ri;
      drji2rj = -drji2ri;
      drji2pi =  dr2pi;
      drji2pj =  dr2pj;

      // derivatives from gamma_{ij}.
      dgami = optP0dev*(1.0/pcm[0]-pcm[0]/s)*b1+pcm/s;
      dgamj = optP0dev*(1.0/pcm[0]-pcm[0]/s)*b2+pcm/s;

      // rest frame of particle i or j.
      } else if(optTwoBodyDistance == 3) {

        double rbi=dot3(rk,vi);
        double rbj=dot3(rk,vj);
        dr2ri =  2*(rk+rbj*vj);    // dR~^2_ij/dR_i
        dr2rj =  -dr2ri;           // dR~^2_ij/dR_j
        dr2pi =  0.0;              // dR~^2_ij/dP_i
        dr2pj =  2*rk*rbj/m2;     // dR~^2_ij/dP_j

        drji2rj = 2*(rk+rbi*vi);
        drji2rj =  -drji2ri;
        drji2pi =  2*rk*rbi/m1;
        drji2pj =  0.0;
      }

      double xdoti=dot3(vk1+forcer[i],dr2ri)   + dot3(vk2+forcer[j],dr2rj);
      double xdotj=dot3(vk2+forcer[j],drji2rj) + dot3(vk1+forcer[i],drji2ri);

      double pdoti=dot3(force[i],dr2pi   - dgami/wG)
	         + dot3(force[j],dr2pj   - dgamj/wG);
      double pdotj=dot3(force[j],drji2pj - dgamj/wG)
	         + dot3(force[i],drji2pi - dgami/wG);

      double doti=-wG*(xdoti+pdoti)*rhom[i][j];
      double dotj=-wG*(xdotj+pdotj)*rhom[j][i];

      Vdot[i] += doti*(Aj + vvi*vj);
      Vdot[j] += dotj*(Ai + vvj*vi);

      // Momentum dependent potential.
      double fmomdi=0.0, fmomdj=0.0;
      if(withMomDep) {

      if(optTwoBodyDistance==1) {
        distanceP1(p1,p2);
      } else if(optTwoBodyDistance==2) {
        distanceP2(p1,p2);
      } else if(optTwoBodyDistance==3) {
        if(optMDarg3==1) distanceP1(p1,p2);
        else if(optMDarg3==2) distanceP2(p1,p2);
        else distanceP3(p1,p2);
      }

      // derivative of rho_{ij} and gamma_{ij}
      fmomdi = devVmd(psq2,pmu2,vex2);
      fmomdj = devVmd(psq1,pmu2,vex2);
      Vdot[i] += doti*fmomdi*vj;
      Vdot[j] += dotj*fmomdj*vi;

      // derivative of D(p_{ij}) part.
      pdoti=dot3(force[i],dp2ijpi) + dot3(force[j],dp2ijpj);
      pdotj=dot3(force[j],dp2jipj) + dot3(force[i],dp2jipi);

      Vdot[i] += pdoti*devVme(psq2,pmu2,vex2)*rhom[i][j]*vj;
      Vdot[j] += pdotj*devVme(psq1,pmu2,vex2)*rhom[j][i]*vi;

      }


      if(opt) continue;

      //Compute derivatives of pre-factors v_j^\mu in the vector potential.
      double fai = - dot3(force[j],Bi/emj) * rhom[i][j];
      double faj = - dot3(force[i],Bj/emi) * rhom[j][i];
      double fai2 = dot3(force[j],vj);
      double faj2 = dot3(force[i],vi); 
      double fai3 = dot3(vj,Bi/p2[0]);
      double faj3 = dot3(vi,Bj/p1[0]);

      double  fvi2=0.0;
      double  fvj2=0.0;
      if(optTwoBodyDistance==3) {
        fai3=Bi[0]/p2[0];
        faj3=Bj[0]/p1[0];
      } else {
        fvi2=-optP0dev*fai2*rhom[i][j]*vvi/p2[0];
        fvj2=-optP0dev*faj2*rhom[j][i]*vvj/p1[0];
      }

      fai += optP0dev*fai2*fai3*rhom[i][j];
      faj += optP0dev*faj2*faj3*rhom[j][i];

      double fvi = vvi*rhom[i][j]/emj;
      double fvj = vvj*rhom[j][i]/emi;

      Vdot[i] += fai*JB[i] + fvi*force[j] + fvi2*vj;
      Vdot[j] += faj*JB[j] + fvj*force[i] + fvj2*vi;

      // Momentum dependent potential.
      if(withMomDep) {
        Vdot[i] += fmomdi*rhom[i][j]*force[j]/emj;
        Vdot[j] += fmomdj*rhom[j][i]*force[i]/emi;
        if(optTwoBodyDistance != 3) {
          Vdot[i] -= optP0dev*fai2*rhom[i][j]*fmomdi/p2[0]*vj;
          Vdot[i] -= optP0dev*faj2*rhom[j][i]*fmomdj/p1[0]*vi;
	}
      }


    } // end j loop
  }   // end i loop

  //return;

  Vec4 vdott=0.0;
  for(int i=0;i<NV;i++) {
    vdott += Vdot[i];
  }

  Vec4 dv = vdott/NV;
  for(int i=0;i<NV;i++) {
    Vdot[i] -= dv;
  }

  /*
  cout << " vdott = "<< scientific
       << setw(13) << setprecision(6) << vdott[1]
       << setw(13) << setprecision(6) << vdott[2]
       << setw(13) << setprecision(6) << vdott[3]
       <<endl;
  cin.get();
  */

}

// non-relativistic two-body distance
// (distance at the global cm frame where the total momentum is zero)
void RQMDw::scalarDensity1(double a)
{
  for(auto i=0; i< NV; i++) rhos[i]=0.0;

  for(int i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double m1 = part[i]->getEffectiveMass();
    p1 *= a; p1[0]= sqrt(m1*m1+p1.pAbs2());
    double fi = m1/p1[0];
    for(int j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double m2 = part[j]->getEffectiveMass();
      p2 *=a; p2[0]= sqrt(m2*m2+p2.pAbs2());
      double fj = m2/p2[0];
      double dens = facG * exp(-(r1 - r2).pAbs2()*wG);
      rhos[i] += dens*fj;
      rhos[j] += dens*fi;
    }
  }

}

// two-body distance at their center-of-momentum frame.
void RQMDw::scalarDensity2(double a)
{
  for(auto i=0; i< NV; i++) rhos[i]=0.0;

  for(int i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double m1 = part[i]->getEffectiveMass();
    p1 *= a; p1[0]= sqrt(m1*m1+p1.pAbs2());
    double fi = m1/p1[0];
    for(int j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double m2 = part[j]->getEffectiveMass();
      p2 *=a; p2[0]= sqrt(m2*m2+p2.pAbs2());
      double fj = m2/p2[0];
      Vec4 dr = r1 - r2; dr[0] = 0.0;  // just in case;
      Vec4 P  = p1 + p2;
      double s = P.m2Calc();
      double drcmsq = dr.m2Calc() - pow2(dr*P)/s;
      double den = P[0]/sqrt(s) * facG * exp(drcmsq*wG);
      rhos[i] += den*fj;
      rhos[j] += den*fi;
    }
  }

}

// two-body distance is computed from the rest-frame of particle j.
void RQMDw::scalarDensity3(double a)
{
  for(auto i=0; i< NV; i++) rhos[i]=0.0;

  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    double m1 = part[i]->getEffectiveMass();
    p1 *= a; p1[0]= sqrt(m1*m1+p1.pAbs2());
    Vec4 v1 = p1/m1;
    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      double m2 = part[j]->getEffectiveMass();
      p2 *=a; p2[0]= sqrt(m2*m2+p2.pAbs2());
      Vec4 v2 = p2/m2;
      Vec4 dr = r1 - r2; dr[0] = 0.0;  // just in case;
      rhos[i] += facG * exp(wG*( dr.m2Calc() - pow2(dr*v2) ) );
      rhos[j] += facG * exp(wG*( dr.m2Calc() - pow2(dr*v1) ) );
    }
  }

}

void RQMDw::qmdMatrix(double a)
{
  for(int i=0;i<NV;i++) {
    rhos[i]=0.0;
    vmoms[i]=0.0;
    vmom4[i]=0.0;
    JB[i]=0.0;
  }

  qmdMatrix0(a);

  /*
  if(optTwoBodyDistance == 1)      qmdMatrix1(a);
  else if(optTwoBodyDistance == 2) qmdMatrix2(a);
  else if(optTwoBodyDistance == 3) qmdMatrix3(a);
  else                             qmdMatrix4(a);
  */

//  if(optPotentialArg == 2) for(int i=0;i<NV;i++) part[i]->setKinetic();

}

void RQMDw::qmdMatrix0(double a)
{
  for(int i=0; i< NV; i++) {
    Vec4 r1  = part[i]->getR();
    Vec4 p1 = optPotentialArg>=1? part[i]->getPcan(optVdot) : part[i]->getP();
    double m1 = part[i]->getEffectiveMass();
    p1 *= a; p1[0]= sqrt(m1*m1+p1.pAbs2());
    distance->setRP1(r1,p1);

    Vec4 v1 = p1/p1[0];
    double fi=p1.mCalc()/p1[0];

    Vec4 fri = optBaryonCurrent ? part[i]->forceR() : 0.0;
    int bi   = part[i]->baryon()/3;
    double qfac1 = part[i]->getTf() > globalTime ? part[i]->qFactor() : 1.0;

    for(int j=i+1; j< NV; j++) {
      Vec4 r2  = part[j]->getR();
      Vec4 p2 = optPotentialArg>=1 ? part[j]->getPcan(optVdot) : part[j]->getP();
      double m2 = part[j]->getEffectiveMass();
      p2 *= a; p2[0]= sqrt(m2*m2+p2.pAbs2());
      distance->setRP2(r2,p2);
      Vec4 v2 = p2/p2[0];
      double fj=p2.mCalc()/p2[0];

      Vec4 frj = optBaryonCurrent ? part[j]->forceR() : 0.0;
      int bj = part[j]->baryon()/3;
      double qfac2 = part[j]->getTf() > globalTime ? part[j]->qFactor() : 1.0;
      
      distance->density();
      rhom[i][j] = distance->density1*qfac1*qfac2;
      rhom[j][i] = distance->density2*qfac1*qfac2;

      // scalar density
      rhos[i] += rhom[i][j]*fj;
      rhos[j] += rhom[j][i]*fi;

      // vector current
      JB[i] += rhom[i][j]*(v2+frj)*bj;
      JB[j] += rhom[j][i]*(v1+fri)*bi;

      if(!withMomDep) continue;
      distance->psq();

      // momentum-dependent scalar potential
      vmoms[i] += vex1/(1.0-distance->psq1/pmu1)*rhom[i][j]*fj;
      vmoms[j] += vex1/(1.0-distance->psq2/pmu1)*rhom[j][i]*fi;

      // momentum-dependent vector potential
      vmom4[i] += vex2/(1.0-distance->psq1/pmu2)*rhom[i][j]*v2;
      vmom4[j] += vex2/(1.0-distance->psq2/pmu2)*rhom[j][i]*v1;


    }
  }

}

void RQMDw::computeForce()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = optVdot <= 1 ? part[i]->getPkin() : p1;
    if(optPotentialArg>=1) p1 = part[i]->getPcan(optVdot);
    int bar1= part[i]->baryon()/3;

    distance->setRP1(r1,p1);
    Vec4 v1 = p1/p1[0];
    double fi = p1.mCalc()/p1[0];
    double fengi=fi;
    if(optPotentialArg==3) {
      double meff1 = part[i]->getEffectiveMass();
      fengi = meff1/sqrt(meff1*meff1 + p1.pAbs2());
    }
    double s1 = dSigma(i);

    for(auto j=i+1; j< NV; j++) {
      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      Vec4 pk2 = optVdot <= 1 ? part[j]->getPkin() : p2;
      if(optPotentialArg>=1) p2 = part[j]->getPcan(optVdot);
      int bar2 = part[j]->baryon()/3;

      distance->setRP2(r2,p2);
      Vec4 v2 = p2/p2[0];
      double fj=p2.mCalc()/p2[0];
      double fengj=fj;
      if(optPotentialArg==3) {
        double meff2 = part[j]->getEffectiveMass();
        fengj = meff2/sqrt(meff2*meff2 + p2.pAbs2());
      }

      distance->distanceR();
      double s2 = dSigma(j);

      // scalar part
      double fsky1 = -wG*fengi*s1*rhom[i][j]*fj; 
      double fsky2 = -wG*fengj*s2*rhom[j][i]*fi; 
      force[i]  += -fsky1*distance->dr2ijri - fsky2*distance->dr2jiri;
      force[j]  += -fsky1*distance->dr2ijrj - fsky2*distance->dr2jirj;
      forcer[i] +=  fsky1*distance->dr2ijpi + fsky2*distance->dr2jipi;
      forcer[j] +=  fsky1*distance->dr2ijpj + fsky2*distance->dr2jipj;

      // vector part
      Vec4 Ai = facV(i,p2);
      Vec4 Aj = facV(j,p1);
      double fskyi = -wG*(Ai * pk1)/pk1[0]*rhom[i][j]*bar1*bar2; 
      double fskyj = -wG*(Aj * pk2)/pk2[0]*rhom[j][i]*bar1*bar2; 
      force[i]  += -fskyi*distance->dr2ijri - fskyj*distance->dr2jirj;
      force[j]  += -fskyi*distance->dr2ijri - fskyj*distance->dr2jirj;
      forcer[i] +=  fskyi*distance->dr2ijpi + fskyj*distance->dr2ijpi;
      forcer[j] +=  fskyi*distance->dr2ijpj + fskyj*distance->dr2jipj;

      // Derivative of p^\mu/m_j term in the vector potential.
      if(optDerivative) {

        if(optTwoBodyDistance!=3) {
        forcer[i] +=  fsky2*p1/(wG*p1[0]*p1[0]);
        forcer[j] +=  fsky1*p2/(wG*p2[0]*p2[0]);
        // gamma derivative.
        distance->devGamma();
        double facsk= -(fsky1+fsky2)/wG;
        forcer[i] += distance->devgam1*facsk;
        forcer[j] += distance->devgam2*facsk;
	}

        Vec4 A1 = facV(i,pk1);
        Vec4 A2 = facV(j,pk2);
	distance->devV(A1,A2);
        forcer[i] += distance->devV1*rhom[j][i]*bar1*bar2;
        forcer[j] += distance->devV2*rhom[i][j]*bar1*bar2;

	double facsk= -(fskyi+fskyj)/wG;
        forcer[i] += distance->devgam1*facsk;
        forcer[j] += distance->devgam2*facsk;
      }

      if(!withMomDep) continue;

      distance->distanceP();
      psq1=distance->psq1;
      psq2=distance->psq2;

      // scalar part
      double fmomd1 = -wG*fengi*devVmd(psq2,pmu1,vex1)*rhom[i][j]*fj;
      double fmomd2 = -wG*fengj*devVmd(psq1,pmu1,vex1)*rhom[j][i]*fi;
      double fmome1 =     fengi*devVme(psq2,pmu1,vex1)*rhom[i][j]*fj;
      double fmome2 =     fengj*devVme(psq1,pmu1,vex1)*rhom[j][i]*fi;

      force[i]  += -fmomd1*distance->dr2ijri - fmomd2*distance->dr2jiri;
      force[j]  += -fmomd1*distance->dr2ijrj - fmomd2*distance->dr2jirj;
      forcer[i] +=  fmome1*distance->dp2ijpi + fmome2*distance->dp2jipi
	          + fmomd1*distance->dr2ijpi + fmomd2*distance->dr2jipi;
      forcer[j] +=  fmome1*distance->dp2ijpj + fmome2*distance->dp2jipj
	          + fmomd1*distance->dr2ijpj + fmomd2*distance->dr2jipj;


      // vector part
      double vf1=1.0;
      double vf2=1.0;
      if(optVectorPotential==1) {
        vf1 = pk1/pk1[0] * v2;
        vf2 = pk2/pk2[0] * v1;
      }

      double fmomdi = -wG*devVmd(psq2,pmu2,vex2)*rhom[i][j]*vf1;
      double fmomdj = -wG*devVmd(psq1,pmu2,vex2)*rhom[j][i]*vf2;
      double fmomei =     devVme(psq2,pmu2,vex2)*rhom[i][j]*vf1;
      double fmomej =     devVme(psq1,pmu2,vex2)*rhom[j][i]*vf2;

      force[i]  += -fmomdi*distance->dr2ijri  - fmomdj*distance->dr2jirj;
      force[j]  += -fmomdi*distance->dr2ijrj  + fmomdj*distance->dr2jirj;
      forcer[i] +=  fmomei*distance->dp2ijpi + fmomej*distance->dp2jipi
	          + fmomdi*distance->dr2ijpi + fmomdj*distance->dr2jipi;
      forcer[j] +=  fmomej*distance->dp2jipj + fmomei*distance->dp2ijpj
	          + fmomdi*distance->dr2ijpj + fmomdj*distance->dr2jipj;

      if(optDerivative) {

        if(optTwoBodyDistance!=3 && optP0dev) {
          forcer[i] +=  fmomd2*p1/(wG*p1[0]*p1[0]);
          forcer[j] +=  fmomd1*p2/(wG*p2[0]*p2[0]);
	}
        // gamma derivative.
        double facmom = -(fmomd1 + fmomd2)/wG;
        forcer[i] += distance->devgam1*facmom;
        forcer[j] += distance->devgam2*facmom;


        // Derivative of p_j / p_0j part.
        double facm1 = devVmd(psq1,pmu2,vex2)*rhom[j][i];
        double facm2 = devVmd(psq2,pmu2,vex2)*rhom[i][j];
	distance->devV(pk1/pk1[0],pk2/pk2[0]);
        forcer[i] += facm1*distance->devV1;
        forcer[j] += facm2*distance->devV2;

	// gamma derivative.
        facmom = facm1*vf1 + facm2*vf2;
        forcer[i] += distance->devgam1*facmom;
        forcer[j] += distance->devgam2*facmom;
      }

    } // end loop over j
  } // end loop over i

}

void RQMDw::computeForceP()
{
  for(auto i=0; i< NV; i++) {
    Vec4 r1 = part[i]->getR();
    Vec4 p1 = part[i]->getP();
    Vec4 pk1 = optVdot <= 1 ? part[i]->getPkin() : p1;
    if(optPotentialArg>=1) p1 = part[i]->getPcan(optVdot);
    int bar1= part[i]->baryon()/3;

    double fi = p1.mCalc()/p1[0];
    double fengi=fi;
    if(optPotentialArg==3) {
      double meff1 = part[i]->getEffectiveMass();
      fengi = meff1/sqrt(meff1*meff1 + pk1.pAbs2());
    }

    distance->setRP1(r1,p1);
    distanceP->setRP1(r1,p1);
    double* gf1=part[i]->facPotential();
    //potential->set1(p1,pk1,JB[i],fi,fengi,rho[i],0.0,
    //                part[i]->pots()-vmoms[i],0.0,0.0,bar1);
    potential->set1(p1,pk1,JB[i],fi,fengi,rho[i],0.0,0.0,
                    part[i]->pots()-vmoms[i],rhos[i],0.0,bar1,gf1[1],gf1[2],gf1[4],gf1[3],gf1[5]);

    for(auto j=i+1; j< NV; j++) {

      Vec4 r2 = part[j]->getR();
      Vec4 p2 = part[j]->getP();
      Vec4 pk2 = optVdot <= 1 ? part[j]->getPkin() : p2;
      if(optPotentialArg>=1) p2 = part[j]->getPcan(optVdot);
      int bar2 = part[j]->baryon()/3;

      double fj=p2.mCalc()/p2[0];
      double fengj=fj;
      if(optPotentialArg==3) {
        double meff2 = part[j]->getEffectiveMass();
        fengj = meff2/sqrt(meff2*meff2 + pk2.pAbs2());
      }

      distance->setRP2(r2,p2);
      distanceP->setRP2(r2,p2);
      double* gf2=part[j]->facPotential();
      //potential->set2(p2,pk2,JB[j],fj,fengj,rho[j],0.0,
      // part[j]->pots()-vmoms[j],0.0,0.0,bar2);
      potential->set2(p2,pk2,JB[j],fj,fengj,rho[j],0.0,0.0,
	  part[j]->pots()-vmoms[j],rhos[j],0.0,bar2,gf2[1],gf2[2],gf2[3],gf2[4],gf2[5]);

      distance->distanceR();
      potential->dVdns(rhom[i][j],rhom[j][i]);

      double fsky1=potential->fskyi;
      double fsky2=potential->fskyj;

      force[i]  += -fsky1*distance->dr2ijri - fsky2*distance->dr2jiri;
      force[j]  += -fsky1*distance->dr2ijrj - fsky2*distance->dr2jirj;
      forcer[i] +=  fsky1*distance->dr2ijpi + fsky2*distance->dr2jipi;
      forcer[j] +=  fsky1*distance->dr2ijpj + fsky2*distance->dr2jipj;

      // Derivative of p0/m and p/p0 term.
      if(optDerivative) {

        potential->dVdp();
        forcer[i] +=  potential->forceri;
        forcer[j] +=  potential->forcerj;

        // gamma derivative.
        distance->devGamma();
        forcer[i] += distance->devgam1*potential->facsk;
        forcer[j] += distance->devgam2*potential->facsk;
      }

      if(!withMomDep) continue;

      distanceP->distanceP();

      potential->dVdmd(distanceP->psq1,distanceP->psq2,gf1[6]*gf2[6],gf1[7]*gf2[7],gf1[8],gf2[8],gf1[9],gf2[9]);
      double fmomdi=potential->fmomdi;
      double fmomdj=potential->fmomdj;
      double fmomei=potential->fmomei;
      double fmomej=potential->fmomej;

      force[i]  += -fmomdi*distance->dr2ijri   - fmomdj*distance->dr2jiri;
      force[j]  += -fmomdi*distance->dr2ijrj   - fmomdj*distance->dr2jirj;
      forcer[i] +=  fmomei*distanceP->dp2ijpi + fmomej*distanceP->dp2jipi
	          + fmomdi*distance->dr2ijpi + fmomdj*distance->dr2ijpi;
      forcer[j] +=  fmomej*distanceP->dp2jipj + fmomei*distanceP->dp2ijpj
	          + fmomdi*distance->dr2ijpj + fmomdj*distance->dr2jipj;

      if(optDerivative) {
        potential->dVdpm();
        forcer[i] +=  potential->forceri;
        forcer[j] +=  potential->forcerj;
        // gamma derivative.
        //distanceP->devGamma();
        forcer[i] += distance->devgam1*potential->facmom;
        forcer[j] += distance->devgam2*potential->facmom;
      }

    } // end loop over j
  } // end loop over i

}

void RQMDw::sigmaField(double a)
{
  int maxit=5;
  const double  eps = 1e-5;
  bool optSet=true;
  if(optPotentialArg > 0 || optVdot <= 1) maxit=1; 

  //if(optTwoBodyDistance == 1)    scalarDensity1(a);
  //else if(optTwoBodyDistance==2) scalarDensity2(a);
  //else scalarDensity3(a);

  double diff=0.0;
  for(int itry=0;itry<maxit;itry++) {
    diff=0.0;
    for(int i=0; i< NV; i++) {
      double vpots0 = part[i]->pots();
      double gfac=1.0;
      scalarDens=gs*rhos[i]/gfac;// factor gr is already contained in rhos(i)
      double gss=t1*gfac*HBARC;

      double sigma1=abs((part[i]->getMass() + vmoms[i])/gss);
      //double sigma1=abs((part[i]->getMass())/gss);
      double vpots=gss*sigma_bisec(sigma1, i ) + vmoms[i];

      //double vpots=gss*sigma_iterate() + vmoms[i];
      //double vpots=gss*sigma_newton() + vmoms[i];

      part[i]->setPotS(vpots,optSet);
      diff += abs(vpots-vpots0);
    }

    //diff=abs(pot0-pot);
    if(diff/NV < eps) return;
      //cout << " diff = "<< diff << " itry= "<< itry << " mxit= "<< maxit << endl;

    if(maxit==1) return;
    qmdMatrix();
    //if(optTwoBodyDistance == 1)    scalarDensity1(a);
    //else if(optTwoBodyDistance==2) scalarDensity2(a);
    //else scalarDensity3(a);

  }

    cout << "sigma does not converge diff= " << diff 
	<< " itry= "<< maxit << endl;
    //exit(1);

}

// Find sigma field in 1/fm
double RQMDw::sigma_bisec(double sigma1, int i)
{
  double sigma0=0.0;
  double f0 = funcSigma(sigma0);
  double f1 = funcSigma(sigma1);
  double sig0=sigma1;
  // seeking a new initial condition.
  if(f0*f1 > 0.0) {
    do {
      sigma1 -= 0.03;
      if(sigma1 < sigma0) {
	  cout << "sigma_bisec no solution  sigma1 = "<< sigma1
	      << " rhos= "<< scalarDens
	      << " sigma0= "<< sig0
	      << " id= "<< part[i]->getID()
	      << " mass= "<< part[i]->getEffectiveMass()
	      <<endl;
	  return sigma1-1e-5;
      }
      f1=funcSigma(sigma1);
      //cout << " sigma1 = " << sigma1  << " f1= "<< f1 <<endl;
      if(f0*f1 < 0.0) break;
    }while(sigma1 > sigma0);
  }

  // check if sigma field is positive.
  if(sigma1 <0.0) {
      cout << " sigma < 0 ? "<< sigma1 << " sig0= "<< sig0
	   <<endl;
      exit(1);
  }

  // check if initial condition is ok.
  if(f0*f1 > 0.0) {
    cout << "bisec no solution sigma= " << sigma1
	 << " dens= "<< scalarDens <<endl;
    cout << "func1= " << funcSigma(sigma0) << " f2= "<< funcSigma(sigma1)
	<<endl;
    return 0.0;
  }

  // bisection method start.
  double f=1.0;
  double sigma=0.0;
  int itry=0;
  do {
    sigma=0.5*(sigma0+sigma1);
    f =  funcSigma(sigma);
    if(f*f1 > 0.0) sigma1=sigma;
    else sigma0=sigma;
    if(++itry > 50) {
      cout << "does not converge sigma_bisec " << sigma << endl;
      exit(1);
    }
  }while (abs(f) > 1e-5);

  return sigma;
}

// Find sigma field in 1/fm
double RQMDw::sigma_iterate()
{
  double sigma=0.0;
  int itry=0;
  double f=1.0;
  do {
    double sigma0=sigma;
    sigma=fSigma(sigma0);
    f =  abs(sigma0-sigma);
    if(++itry > 50) {
      cout << "does not converge sigma_iterate " << sigma << endl;
      exit(1);
    }
  }while (abs(f) > 1e-5);
  cout << " sigma= "<< sigma << " f= "<< f <<" ity= "<< itry<<endl;

  return sigma;
}

// Find sigma field in 1/fm
double RQMDw::sigma_newton(double s0)
{
  const double fac=1.0;
  double sigma=s0;
  int itry=0;
  double f=1.0;
  double df0 = funcdSigma(sigma);  // for simplified Newton method.
  do {
    double sigma0=sigma;

    // Newton method.
    //sigma -= fac*funcSigma(sigma0)/funcdSigma(sigma0);

    // Simplified Newton method.
    sigma -= fac*funcSigma(sigma0)/df0;

    // Modified Newton method.
    //double fx = funcSigma(sigma0);
    //double df = funcdSigma(sigma0);
    //double f0 = funcSigma(sigma0 - fx/df);
    //sigma -= fx*fx / ( df * (fx - f0) );
   
    f =  abs(sigma0-sigma);
    if(++itry > 100) {
      cout << "does not converge sigma_newton " << sigma
	  << " eps= "<< f
	   << " rhos= "<< scalarDens << endl;
      exit(1);
    }
  }while (abs(f) > 1e-5);
  cout << " sigma= "<< sigma << " f= "<< f <<" ity= "<< itry
       << " rhos= "<< scalarDens <<endl;

  return sigma;
}

// relative distance and its derivatives between p1 and p2
void RQMDw::distanceP1(const Vec4& p1,const Vec4& p2)
{
  Vec4 dp = p1 - p2;
  psq1 = -dp.pAbs2() + (1-optPV) * dp[0]*dp[0];
  psq2 = psq1;
  dp2ijpi = 2*( dp -optP0dev*(1-optPV)*dp[0]*p1/p1.e());
  dp2ijpj = 2*(-dp -optP0dev*(1-optPV)*dp[0]*p2/p2.e());
  dp2jipi = dp2ijpi;
  dp2jipj = dp2ijpj;
}

// relative distance and its derivatives between p1 and p2
// in the two body CM frame.
void RQMDw::distanceP2(const Vec4& p1,const Vec4& p2)
{
  Vec4 dp  = p1 - p2;
  Vec4 pcm = p1 + p2;
  double s = pcm.m2Calc();
  //double pma = pow2(m1*m1 - m2*m2)/s;
  double pma = pow2(dp * pcm)/s;
  psq1 = dp.m2Calc() - optPV * pma;
  psq2 = psq1;
  Vec4 b1 = p1/p1.e();
  Vec4 b2 = p2/p2.e();
  Vec4 bbi = pcm/pcm[0] - optP0dev*b1;
  Vec4 bbj = pcm/pcm[0] - optP0dev*b2;
  dp2ijpi = 2*(dp - optP0dev*dp[0]*b1 + optPV * pcm[0]/s*pma*bbi);
  dp2ijpj = 2*(-dp + optP0dev*dp[0]*b2 + optPV * pcm[0]/s*pma*bbj);
  dp2jipi = dp2ijpi;
  dp2jipj = dp2ijpj;
}

// relative distance squared and its derivatives between p1 and p2
// in the rest frame of p1 or p2.
void RQMDw::distanceP3(const Vec4& p1,const Vec4& p2)
{
  Vec4 dp = p1 - p2;
  double psq = dp.m2Calc();

  // distance squared in the rest frame of particle 2
  double m2 = p2.mCalc();
  double dot4j = dp * p2 / m2;
  psq2 = psq - optPV*dot4j*dot4j;

  // distance squared in the rest frame of particle 1
  double m1 = p1.mCalc();
  double dot4i = (dp * p1)/m1;
  psq1 = psq - optPV*dot4i*dot4i;

  // derivatives
  Vec4 bi = p1/p1.e();
  Vec4 bb = optP0dev * p2.e()*bi - p2;
  dp2ijpi = 2*(dp - optP0dev*dp[0]*bi + optPV*bb*dot4j/m2);
  dp2jipi = 2*(dp - optP0dev*dp[0]*bi - optPV*bb*dot4i/m1);

  Vec4 bj = p2/p2.e();
  Vec4 bb2 = optP0dev * p1.e()*bj - p1;
  dp2ijpj = 2*(-dp + optP0dev*dp[0]*bj + optPV*bb2*dot4j/m2);
  dp2jipj = 2*(-dp + optP0dev*dp[0]*bj - optPV*bb2*dot4i/m1);
}


} // namespace jam2


