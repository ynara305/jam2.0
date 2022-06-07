#include <jam2/interaction/BeamHadron.h>
#include <jam2/hadrons/JamStdlib.h>

namespace jam2 {


using namespace std;
//using namespace Pythia8;

BeamHadron::BeamHadron(Settings& settings, ParticleData* pd, Rndm* r) 
{
    rndm=r;
    particleData = pd;

  // Enhancement factor of x of diquark.
  valenceDiqE0 = settings.parm("BeamRemnants:valenceDiqEnhanceSoft");

  // Power of (1-x)^power/sqrt(x) for remnant valence quark distribution.
  valencePowerMeson = settings.parm("BeamRemnants:valencePowerMeson");
  valencePowerUinP  = settings.parm("BeamRemnants:valencePowerUinP");
  valencePowerDinP  = settings.parm("BeamRemnants:valencePowerDinP");

  // Assume g(x) ~ (1-x)^power/x to constrain companion to sea quark.
  companionPower    = settings.mode("BeamRemnants:companionPower");

  // Assume g(x) ~ (1-x)^power/x with a cut-off for low x.
  gluonPower        = settings.parm("BeamRemnants:gluonPower");
  xGluonCutoff      = settings.parm("BeamRemnants:xGluonCutoff");

  xGluonCutoff      = 1e-4;

  eMinPert=13.0;
  eWidthPert=1.0;
  valenceDiqE=6.0;

}

void BeamHadron::init(double ecm,int id, int cq[2], Vec4 p, double m,int kf1, int kf2, int pickg)
{ 
  eCM=ecm;
  valenceDiqEnhance=valenceDiqE0 + valenceDiqE/(1.0 + exp(-(eCM-eMinPert)/eWidthPert));

  resolved.clear();
  idBeam=id;
  pBeam=p;
  mBeam=m;
  nValenceq=2;
  pickG=pickg;

  initBeamKind();
  int junction=0;


  bool anti=false;
  if(abs(kf1) < 10) {
      if(kf1<0) anti=true;
  } else if(abs(kf1)>1000) {
      if(kf1>0) anti=true;
  } else {
      cout << " BeamHadron::init kf1= " << kf1 << endl;
      exit(1);
  }

  // check constituent quarks are virtual or not.
  int cq1=cq[0];
  int cq2=cq[1];
  int cq3=0;
  if(id<0 && cq2==2) {
      cq1=cq[1];
      cq2=cq[0];
  }
  if(junction && (cq2==2 || cq1==2)) {
      cq2=1;
      cq3=1;
  }

 // cout << "id= "<< id << " cq1= "<< cq1 << " cq2= "<< cq2 << " cq3= "<< cq3 <<endl;

  resolved.push_back(ResolvedPartonQ(cq1,0,kf1) );
  if(!anti) 
      resolved[0].cols(101,0);
  else
      resolved[0].cols(0,101);

  if(pickG) {
      for(int i=0;i<pickG;i++)
	 resolved.push_back(ResolvedPartonQ(0,0,21) );
  }

    resolved.push_back(ResolvedPartonQ(cq2,0,kf2) );
    if(anti) resolved.back().col(101);
    else  resolved.back().acol(101);


      if(junction && isBaryonBeam) {
          int id2 = (kf2 > 0) ? (kf2/100)%10 : -(((-kf2)/100)%10);
          int id1 = (kf2 > 0) ?  kf2/1000 : -((-kf2)/1000);
	  nValenceq=3;

	  resolved[resolved.size()-1].id(id1);
	  if(id1>0) resolved.back().col(102);
	  else  resolved.back().acol(102);

	  resolved.push_back(ResolvedPartonQ(cq3,0,id2));
	  if(id2>0) resolved.back().col(103);
	  else  resolved.back().acol(103);
      }

      for(int i=0;i<(int)resolved.size();i++) {
	resolved[i].m( particleData->m0( resolved[i].id() ) );
      }

}

// Pick unrescaled x values for beam remnant sharing.

double BeamHadron::xRemnant( int i) {

  double x = 0.;

    // Resolve diquark into sum of two quarks.
    int id1 = resolved[i].id();
    int id2 = 0;
    if (abs(id1) > 1000) {
      id2 = (id1 > 0) ? (id1/100)%10 : -(((-id1)/100)%10);
      id1 = (id1 > 0) ? id1/1000 : -((-id1)/1000);
    }

    if(id1==21) {
	do x = pow(xGluonCutoff, 1 - rndm->flat());
	while ( pow(1. - x, gluonPower) < rndm->flat() );

    } else {

    // Loop over (up to) two quarks; add their contributions.
    for (int iId = 0; iId < 2; ++iId) {
      int idNow = (iId == 0) ? id1 : id2;
      if (idNow == 0) break;
      double xPart = 0.;

      // Assume form (1-x)^a / sqrt(x).
      double xPow = valencePowerMeson;
      if (isBaryonBeam) {
        if (nValKinds == 3 || nValKinds == 1)
          xPow = (3. * rndm->flat() < 2.)
            ? valencePowerUinP : valencePowerDinP ;
        else if (nValence(idNow) == 2) xPow = valencePowerUinP;
        else xPow = valencePowerDinP;
      }
      do {
	  xPart = pow2( rndm->flat() );
      } while ( pow(1. - xPart, xPow) < rndm->flat() );

     //cout << iId << " xpert= " << xPart << " xPow= " << xPow <<endl;

      // End loop over (up to) two quarks. Possibly enhancement for diquarks.
      x += xPart;
    }
   if (id2 != 0) x *= valenceDiqEnhance;
   //if (id2 != 0) x *= valenceDiqEnhance + 6.0/(1.0 + exp(-(eCM-eMinPert)/eWidthPert));

    }

  return x;

}

void BeamHadron::initBeamKind() {

  // Reset.
  int idBeamAbs         = abs(idBeam);
  isHadronBeam      = false;
  isMesonBeam       = false;
  isBaryonBeam      = false;

  nValKinds         = 0;

  if(idBeamAbs<100) {
  cout << "BeamHadron:: idbeam? " << idBeam <<endl;
  exit(1);
  };

  // Resolve valence content for assumed meson. Flunk unallowed codes.
  //if (idBeamAbs < 1000) {
  if ((idBeamAbs/1000)%10 == 0) {
    int id1 = (idBeamAbs/100)%10;
    int id2 = (idBeamAbs/10)%10;
    int idnew=id1*10+id2;

    isMesonBeam = true;

    // Store valence content of a confirmed meson.
    nValKinds = 2;
    nVal[0]   = 1 ;
    nVal[1]   = 1;
    if (id1%2 == 0) {
      idVal[0] = id1;
      idVal[1] = -id2;
    } else {
      idVal[0] = id2;
      idVal[1] = -id1;
    }
    newValenceContent(idnew);

  // Resolve valence content for assumed baryon. Flunk unallowed codes.
  } else {
    int id1 = (idBeamAbs/1000)%10;
    int id2 = (idBeamAbs/100)%10;
    int id3 = (idBeamAbs/10)%10;

    //int idnew=id1*100+id2*10+id3;
    isBaryonBeam = true;

    // Store valence content of a confirmed baryon.
    nValKinds = 1; idVal[0] = id1; nVal[0] = 1;
    if (id2 == id1) ++nVal[0];
    else {
      nValKinds = 2;
      idVal[1]  = id2;
      nVal[1]   = 1;
    }
    if (id3 == id1) ++nVal[0];
    else if (id3 == id2) ++nVal[1];
    else {
      idVal[nValKinds] = id3;
      nVal[nValKinds]  = 1;
      ++nValKinds;
    }
  }

  // Flip flavours for antimeson or antibaryon, and then done.
  if (idBeam < 0) for (int i = 0; i < nValKinds; ++i) idVal[i] = -idVal[i];
  isHadronBeam = true;
  //Q2ValFracSav = -1.;

}

void BeamHadron::newValenceContent(int idnew) {

  // A pi0, rho and omega oscillates between d dbar and u ubar.
  //if (idBeam == 111 || idBeam == 113 || idBeam == 223) {
  if (idnew == 11 || idnew == 22) {
    idVal[0] = (rndm->flat() < 0.5) ? 1 : 2;
    idVal[1] = -idVal[0];

  // A K0S or K0L oscillates between d sbar and s dbar.
  } else if (idBeam == 130 || idBeam == 310) {
    idVal[0] = (rndm->flat() < 0.5) ?  1 :  3;
    idVal[1] = (idVal[0] == 1)      ? -3 : -1;

  // If phi meson set content to s sbar.
  } else if (idBeam == 333) {
    idVal[0] = 3;
    idVal[1] = -idVal[0];

  // Other hadrons so far do not require any event-by-event change.
  } else return;

}

} // end namespace jam2
