// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// ParticleDecays.cc is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// ParticleDecays class.

#include <jam2/interaction/JParticleDecays.h>
#include <jam2/hadrons/JamStdlib.h>

namespace jam2 {

using namespace std;
//using Pythia8::pow2;
using Pythia8::pow3;
using Pythia8::pow4;
using Pythia8::sqrtpos;
using Pythia8::Event;
using Pythia8::Particle;
using Pythia8::Vec4;

//==========================================================================

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of times one tries to let decay happen (for 2 nested loops).
//const int    JParticleDecays::NTRYDECAY   = 20;
const int    JParticleDecays::NTRYDECAY   = 50;

// Number of times one tries to pick valid hadronic content in decay.
const int    JParticleDecays::NTRYPICK    = 100;

// Number of times one tries to pick decay topology.
const int    JParticleDecays::NTRYMEWT    = 1000;

// Maximal loop count in Dalitz decay treatment.
const int    JParticleDecays::NTRYDALITZ  = 1000;

// Minimal Dalitz pair mass this factor above threshold.
const double JParticleDecays::MSAFEDALITZ = 1.000001;

// These numbers are hardwired empirical parameters,
// intended to speed up the M-generator.
const double JParticleDecays::WTCORRECTION[11] = { 1., 1., 1.,
  2., 5., 15., 60., 250., 1250., 7000., 50000. };

//--------------------------------------------------------------------------

// Initialize and save pointers.

void JParticleDecays::init(Pythia8::Info* infoPtrIn, 
	Pythia8::Settings& settings,
  Pythia8::ParticleData* particleDataPtrIn, Pythia8::Rndm* rndmPtrIn)
{
    optAngle=0;

  // Save pointers to error messages handling and flavour generation.
  infoPtr         = infoPtrIn;
  particleDataPtr = particleDataPtrIn;
  rndmPtr         = rndmPtrIn;

  // Safety margin in mass to avoid troubles.
  mSafety       = settings.parm("ParticleDecays:mSafety");

  // Lifetime and vertex rules for determining whether decay allowed.
  limitTau0     = settings.flag("ParticleDecays:limitTau0");
  tau0Max       = settings.parm("ParticleDecays:tau0Max");
  limitTau      = settings.flag("ParticleDecays:limitTau");
  tauMax        = settings.parm("ParticleDecays:tauMax");
  limitRadius   = settings.flag("ParticleDecays:limitRadius");
  rMax          = settings.parm("ParticleDecays:rMax");
  limitCylinder = settings.flag("ParticleDecays:limitCylinder");
  xyMax         = settings.parm("ParticleDecays:xyMax");
  zMax          = settings.parm("ParticleDecays:zMax");
  limitDecay    = limitTau0 || limitTau || limitRadius || limitCylinder;

  // B-Bbar mixing parameters.
  mixB          = settings.flag("ParticleDecays:mixB");
  xBdMix        = settings.parm("ParticleDecays:xBdMix");
  xBsMix        = settings.parm("ParticleDecays:xBsMix");

  // Suppression of extra-hadron momenta in semileptonic decays.
  sigmaSoft     = settings.parm("ParticleDecays:sigmaSoft");

  // Selection of multiplicity and colours in "phase space" model.
  multIncrease     = settings.parm("ParticleDecays:multIncrease");
  multIncreaseWeak = settings.parm("ParticleDecays:multIncreaseWeak");
  multRefMass      = settings.parm("ParticleDecays:multRefMass");
  multGoffset      = settings.parm("ParticleDecays:multGoffset");
  colRearrange     = settings.parm("ParticleDecays:colRearrange");

  // Minimum energy in system (+ m_q) from StringFragmentation.
  stopMass      = settings.parm("StringFragmentation:stopMass");

  // Parameters for Dalitz decay virtual gamma mass spectrum.
  sRhoDal       = pow2(particleDataPtr->m0(113));
  wRhoDal       = pow2(particleDataPtr->mWidth(113));

  // Allow showers in decays to qqbar/gg/ggg/gammagg.
  doFSRinDecays = settings.flag("ParticleDecays:FSRinDecays");
  doGammaRad    = settings.flag("ParticleDecays:allowPhotonRadiation");

  // u-t ratio in decay
  utRatiot       = settings.parm("ParticleDecays:utRatio");
  utRatios       = settings.parm("ParticleDecays:utRatioS");


  // Gaussian width for decay
  gWidth        = settings.parm("ParticleDecays:gWidth");
}

//--------------------------------------------------------------------------
//
bool JParticleDecays::decay2( int iDec, bool schan, Event& event,
	double totwid, double* pWidth, int nChannel)
{
  sChannel = schan;
  utRatio = sChannel ? utRatios : utRatiot;

  // Check whether a decay is allowed, given the upcoming decay vertex.
  Particle& decayer = event[iDec];
  if (limitDecay && !checkVertex(decayer)) return true;

  // Do not allow resonance decays (beyond handling capability).
  if (decayer.isResonance()) {
    cout << "Warning in JParticleDecays::decay: resonance left undecayed"<<endl;
    return true;
  }

  // Fill the decaying particle in slot 0 of arrays.
  idDec = decayer.id();
  iProd.resize(0);
  idProd.resize(0);
  mProd.resize(0);
  iProd.push_back( iDec );
  idProd.push_back( idDec );
  mProd.push_back( decayer.m() );

  // Particle data for decaying particle.
  decDataPtr = &decayer.particleDataEntry();

    // Allow up to ten tries to pick a channel.
    //if (!decDataPtr->preparePick(idDec, decayer.m())) return false;

    bool foundChannel = false;
    bool hasStored    = false;
    for (int iTryChannel = 0; iTryChannel < NTRYDECAY; ++iTryChannel) {

      // Remove previous failed channel.
      if (hasStored) event.popBack(mult);
      hasStored = false;

      // Pick new channel. Read out basics.
      //Pythia8::DecayChannel& channel = decDataPtr->pickChannel();
    double xwid = totwid*rndmPtr->flat();
    int ibra=-1;
    for(int i=0;i<nChannel;i++) {
      xwid -= pWidth[i];
      if(xwid <= 0.0) {
        ibra=i; break;
      }
    }

    if(ibra==-1) {
	cout << " no channel ? idDec= "<< idDec
	     << scientific << " totwid= "<< totwid
	     << " nch= "<< nChannel
	     <<endl;
	double ttw=0.0;
        for(int i=0;i<nChannel;i++) {
	    ttw += pWidth[i];
	    cout << scientific << " pWidth= "<< pWidth[i]<<endl;
	}
	cout << " ttw = "<< ttw 
             << " xwid= "<< xwid
	     <<endl;
	exit(1);
    }

    //cout << "*** decay = "<< idDec << " m= "<< decayer.m()<<endl;
    //cout << " ibra= "<< ibra << " pwid= "<< pWidth[ibra]<<endl;

      Pythia8::DecayChannel& channel = decDataPtr->channel(ibra);
      meMode = channel.meMode();
      mult = channel.multiplicity();

      // Allow up to ten tries for each channel (e.g with different masses).
      bool foundMode = false;
      iProd.resize(1);
      for (int iTryMode = 0; iTryMode < NTRYDECAY; ++iTryMode) {
        idProd.resize(1);
        mProd.resize(1);
        scale = 0.;

        // Extract and store the decay products in local arrays.
	double totmas = 0.0;
        for (int i = 0; i < mult; ++i) {
          int idNow = channel.product(i);
	  double mr = particleDataPtr->mWidth(idNow) > 1e-7 ?
	     particleDataPtr->mMin(idNow) : particleDataPtr->m0(idNow);
	  totmas += mr;
	}

        for (int i = 0; i < mult; ++i) {
          int idNow = channel.product(i);
          if (idDec < 0 && particleDataPtr->hasAnti(idNow)) idNow = -idNow;
          //double mNow = particleDataPtr->mSel(idNow);
	  double mNow = BWmass(idNow, totmas);
	  if(mNow < 0.0) continue;
	  //cout << "idNow= "<< idNow << " mNow= "<< mNow <<endl;

          idProd.push_back( idNow);
          mProd.push_back( mNow);
        }

	if(totmas > mProd[0]) {
	    cout << " totmas ? "<<totmas << " mProd0= "<< mProd[0]
		<< " id= " << idDec
                << " M= "<<  decayer.m()
		<< " br= "<< pWidth[ibra]
		<<endl;
           for (int i = 0; i < mult; ++i) {
            int idNow = channel.product(i);
	    cout << " id= "<< idNow  << " wid= " << particleDataPtr->mWidth(idNow)
	     << " mmin= "<<  particleDataPtr->mMin(idNow)
	     << " m0= "<< particleDataPtr->m0(idNow)
	    <<endl;
	}
	    exit(1);
	}

        // Need to set colour flow if explicit decay to partons.
        cols.resize(0);
        acols.resize(0);
        for (int i = 0; i <= mult; ++i) {
          cols.push_back(0);
          acols.push_back(0);
        }

        // Check that enough phase space for decay.
        if (mult > 1) {
          double mDiff = mProd[0];
          for (int i = 1; i <= mult; ++i) mDiff -= mProd[i];
          if (mDiff < mSafety) continue;
        }

        // End of inner trial loops. Check if succeeded or not.
        foundMode = true;
        break;
      }
      if (!foundMode) continue;


      // Store decay products in the event record.
      for (int i = 1; i <= mult; ++i) {
        int iPos = event.append( idProd[i], 91, iDec, 0, 0, 0,
          cols[i], acols[i], Vec4(0., 0., 0., 0.), mProd[i], scale);
        iProd.push_back( iPos);
	 //cout << " decay id= "<< idProd[i] << " m= "<< mProd[i]<<endl;
      }
      hasStored = true;

      // Pick mass of Dalitz decay. Temporarily change multiplicity.
      if ( (meMode == 11 || meMode == 12 || meMode == 13)
        && !dalitzMass() ) continue;

      // Do a decay, split by multiplicity.
      bool decayed = false;
      if      (mult == 1) decayed = oneBody(event);
      else if (mult == 2) decayed = twoBody(event);
      else if (mult == 3) decayed = threeBody(event);
      else                decayed = mGenerator(event);
      if (!decayed) continue;

      // Kinematics of gamma* -> l- l+ in Dalitz decay. Restore multiplicity.
      if (meMode == 11 || meMode == 12 || meMode == 13)
        dalitzKinematics(event);

      // End of outer trial loops.
      foundChannel = true;
      break;
    }

    // If the decay worked, then mark mother decayed and store daughters.
    if (foundChannel) {
      event[iDec].statusNeg();
      event[iDec].daughters( iProd[1], iProd[mult]);

    // Else remove unused daughters and return failure.
    } else {
      if (hasStored) event.popBack(mult);
      cout << "Error in JParticleDecays::decay: "
        "failed to find workable decay channel" << endl;
      cout << " id= " << idDec << "  " <<
	    decayer.name() << " m= " << decayer.m()<<endl;
      return false;
    }

  // Set decay vertex when this is displaced.
  if (event[iDec].hasVertex() || event[iDec].tau() > 0.) {
    Vec4 vDec = event[iDec].vDec();
    for (int i = event[iDec].daughter1(); i <= event[iDec].daughter2(); ++i)
      event[i].vProd( vDec );
  }

  // Set lifetime of daughters. Check if they decayed in their turn.
  for (int i = iProd[1]; i <= iProd[mult]; ++i) {
    event[i].tau( event[i].tau0() * rndmPtr->exp() );
    if (!event[i].isFinal() && (event[i].hasVertex() || event[i].tau() > 0.)) {
      Vec4 vDecR = event[i].vDec();
      for (int iR = event[i].daughter1(); iR <= event[i].daughter2(); ++iR)
        event[iR].vProd( vDecR );
    }
  }

  // Done.
  return true;

}

double JParticleDecays::BWmass(int idNow,double totmas)
{
//...Purpose: to generate mass according to the B.W. distribution.
//==================================================================*
//  emmin : minimam mass           (input)
//  emmax : max. mass              (input)
//  emr   : resonance peak mass    (input)
//  wid   : resonance full width   (input)
//  em    : resonance mass         (output)
//==================================================================*

  double wid=particleDataPtr->mWidth(idNow);
  double emr=particleDataPtr->m0(idNow);
  if(wid< 1e-7) return emr;
  double emmin=particleDataPtr->mMin(idNow);
  double emmax=mProd[0] - totmas + emmin - eKinMin;
  if(abs(emmax-emmin)<1e-8) return emmin;

  // Check boundary.
  if(emmax < emmin) {
      cout << " idDec= "<< idDec 
	  << " BWmass id = " << idNow 
	  << scientific << " mprod0 = " << mProd[0]
	  << " totmas = " << totmas
	  << " emr= "<<  emr
	  << " mmin= "<< emmin
	  << " mmax "<< emmax
	  << " def= "<< abs(emmax-emmin)
	  <<endl;
      return -1.0;
      //return emmax;
  }

//...Breit Wigner distribution.
  double cnst=2.0*(emmin-emr)/wid;
  double const1=atan(cnst);
  double const2=M_PI/2.0-const1;
  double xmax=(atan(2.0*(emmax-emr)/wid)-const1)/const2;
  xmax=min(1.0,xmax);
  double x=xmax*rndmPtr->flat();
  double t=tan(x*const2);
  double em=emr+0.5*wid*(cnst+t)/(1.0-cnst*t);
  
  return em;

}

// Decay a particle; main method.

bool JParticleDecays::decay( int iDec, bool schan, Event& event) {

  sChannel = schan;
  utRatio = sChannel ? utRatios : utRatiot;

  // Check whether a decay is allowed, given the upcoming decay vertex.
  Particle& decayer = event[iDec];
  if (limitDecay && !checkVertex(decayer)) return true;

  // Do not allow resonance decays (beyond handling capability).
  if (decayer.isResonance()) {
    cout << "Warning in JParticleDecays::decay: resonance left undecayed"<<endl;
    return true;
  }

  // Fill the decaying particle in slot 0 of arrays.
  idDec = decayer.id();
  iProd.resize(0);
  idProd.resize(0);
  mProd.resize(0);
  iProd.push_back( iDec );
  idProd.push_back( idDec );
  mProd.push_back( decayer.m() );

  // Particle data for decaying particle.
  decDataPtr = &decayer.particleDataEntry();

    // Allow up to ten tries to pick a channel.
    //if (!decDataPtr->preparePick(idDec, decayer.m())) return false;

    if (!decDataPtr->preparePick(idDec, decayer.m())) {
      cout << " allow up to ten tries?  idDec= " << idDec 
	  << " name= "<< decayer.name()
	   << " m= " << decayer.m() <<endl;
    	return false;
    }
    bool foundChannel = false;
    bool hasStored    = false;
    for (int iTryChannel = 0; iTryChannel < NTRYDECAY; ++iTryChannel) {

      // Remove previous failed channel.
      if (hasStored) event.popBack(mult);
      hasStored = false;

      // Pick new channel. Read out basics.
      Pythia8::DecayChannel& channel = decDataPtr->pickChannel();
      meMode = channel.meMode();
      mult = channel.multiplicity();

      //cout << "iTryCh= "<< iTryChannel << " mult= "<< mult<<endl;

      // Allow up to ten tries for each channel (e.g with different masses).
      bool foundMode = false;
      iProd.resize(1);
      for (int iTryMode = 0; iTryMode < NTRYDECAY; ++iTryMode) {
        idProd.resize(1);
        mProd.resize(1);
        scale = 0.;

        // Compute total minium mass.
	double totmas = 0.0;
        for (int i = 0; i < mult; ++i) {
          int idNow = channel.product(i);
	  double mr = particleDataPtr->mWidth(idNow) > 1e-7 ?
	     particleDataPtr->mMin(idNow) : particleDataPtr->m0(idNow);
	  totmas += mr;
	}

        // Extract and store the decay products in local arrays.
        for (int i = 0; i < mult; ++i) {
          int idNow = channel.product(i);
          if (idDec < 0 && particleDataPtr->hasAnti(idNow)) idNow = -idNow;

          //double mNow = particleDataPtr->mSel(idNow);
	  double mNow = BWmass(idNow, totmas);

          idProd.push_back( idNow);
          mProd.push_back( mNow);
        }

        // Need to set colour flow if explicit decay to partons.
        cols.resize(0);
        acols.resize(0);
        for (int i = 0; i <= mult; ++i) {
          cols.push_back(0);
          acols.push_back(0);
        }

        // Check that enough phase space for decay.
        if (mult > 1) {
          double mDiff = mProd[0];
          for (int i = 1; i <= mult; ++i) mDiff -= mProd[i];
          if (mDiff < mSafety) continue;
        }

	/*
          for (int i = 1; i <= mult; ++i) {
                 double emmin=particleDataPtr->mMin(idProd[i]);
                 double emmax=mProd[0] - totmas + emmin - eKinMin;
	      cout << " idprod= "<< idProd[i]
		  << " m= " << mProd[i]
		  << " emin= "<< emmin
		  << " emax= "<< emmax
		  <<endl;
	  }
	  */


        // End of inner trial loops. Check if succeeded or not.
        foundMode = true;
        break;
      }
      if (!foundMode) continue;


      // Store decay products in the event record.
      for (int i = 1; i <= mult; ++i) {
        int iPos = event.append( idProd[i], 91, iDec, 0, 0, 0,
          cols[i], acols[i], Vec4(0., 0., 0., 0.), mProd[i], scale);
        iProd.push_back( iPos);
	 //cout << " decay id= "<< idProd[i] << " m= "<< mProd[i]<<endl;
      }
      hasStored = true;

      // Pick mass of Dalitz decay. Temporarily change multiplicity.
      if ( (meMode == 11 || meMode == 12 || meMode == 13)
        && !dalitzMass() ) continue;

      // Do a decay, split by multiplicity.
      bool decayed = false;
      if      (mult == 1) decayed = oneBody(event);
      else if (mult == 2) decayed = twoBody(event);
      else if (mult == 3) decayed = threeBody(event);
      else                decayed = mGenerator(event);
      if (!decayed) continue;

      // Kinematics of gamma* -> l- l+ in Dalitz decay. Restore multiplicity.
      if (meMode == 11 || meMode == 12 || meMode == 13)
        dalitzKinematics(event);

      // End of outer trial loops.
      foundChannel = true;
      break;
    }

    // If the decay worked, then mark mother decayed and store daughters.
    if (foundChannel) {
      event[iDec].statusNeg();
      event[iDec].daughters( iProd[1], iProd[mult]);

    // Else remove unused daughters and return failure.
    } else {
      if (hasStored) event.popBack(mult);
      cout << "Error in JParticleDecays::decay: "
        "failed to find workable decay channel" << endl;
      cout << " id= " << idDec << "  " <<
	    decayer.name() << " m= " << decayer.m()<<endl;
      return false;
    }

  // Set decay vertex when this is displaced.
  if (event[iDec].hasVertex() || event[iDec].tau() > 0.) {
    Vec4 vDec = event[iDec].vDec();
    for (int i = event[iDec].daughter1(); i <= event[iDec].daughter2(); ++i)
      event[i].vProd( vDec );
  }

  // Set lifetime of daughters. Check if they decayed in their turn.
  for (int i = iProd[1]; i <= iProd[mult]; ++i) {
    event[i].tau( event[i].tau0() * rndmPtr->exp() );
    if (!event[i].isFinal() && (event[i].hasVertex() || event[i].tau() > 0.)) {
      Vec4 vDecR = event[i].vDec();
      for (int iR = event[i].daughter1(); iR <= event[i].daughter2(); ++iR)
        event[iR].vProd( vDecR );
    }
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Check whether a decay is allowed, given the upcoming decay vertex.

bool JParticleDecays::checkVertex(Particle& decayer) {

  // Check whether any of the conditions are not fulfilled.
  if (limitTau0 && decayer.tau0() > tau0Max) return false;
  if (limitTau && decayer.tau() > tauMax) return false;
  if (limitRadius && pow2(decayer.xDec()) + pow2(decayer.yDec())
    + pow2(decayer.zDec()) > pow2(rMax)) return false;
  if (limitCylinder && (pow2(decayer.xDec()) + pow2(decayer.yDec())
    > pow2(xyMax) || abs(decayer.zDec()) > zMax) ) return false;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Do a one-body decay. (Rare; e.g. for K0 -> K0_short.)

bool JParticleDecays::oneBody(Event& event) {

  // References to the particles involved.
  Particle& decayer = event[iProd[0]];
  Particle& prod    = event[iProd[1]];

  // Set momentum and expand mother information.
  prod.p( decayer.p() );
  prod.m( decayer.m() );
  prod.mother2( iProd[0] );

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Do a two-body decay.

bool JParticleDecays::twoBody(Event& event) {

  // References to the particles involved.
  Particle& decayer = event[iProd[0]];
  Particle& prod1   = event[iProd[1]];
  Particle& prod2   = event[iProd[2]];

  // Masses.
  double m0   = mProd[0];
  double m1   = mProd[1];
  double m2   = mProd[2];

  // Energies and absolute momentum in the rest frame.
  if (m1 + m2 + mSafety > m0) return false;
  double e1   = 0.5 * (m0*m0 + m1*m1 - m2*m2) / m0;
  double e2   = 0.5 * (m0*m0 + m2*m2 - m1*m1) / m0;
  double pAbs = 0.5 * sqrtpos( (m0 - m1 - m2) * (m0 + m1 + m2)
    * (m0 + m1 - m2) * (m0 - m1 + m2) ) / m0;


  // When meMode = 2, for V -> PS2 + PS3 (V = vector, pseudoscalar),
  // need to check if production is PS0 -> PS1/gamma + V.
  int iMother = event[iProd[0]].mother1();
  int idSister = 0;
  if (meMode == 2) {
    if (iMother <= 0 || iMother >= iProd[0]) meMode = 0;
    else {
      int iDaughter1 = event[iMother].daughter1();
      int iDaughter2 = event[iMother].daughter2();
      if (iDaughter2 != iDaughter1 + 1) meMode = 0;
      else {
        int idMother = abs( event[iMother].id() );
        if (idMother <= 100 || idMother%10 !=1
          || (idMother/1000)%10 != 0) meMode = 0;
        else {
          int iSister = (iProd[0] == iDaughter1) ? iDaughter2 : iDaughter1;
          idSister = abs( event[iSister].id() );
          if ( (idSister <= 100 || idSister%10 !=1
            || (idSister/1000)%10 != 0) && idSister != 22) meMode = 0;
        }
      }
    }
  }

  // Begin loop over matrix-element corrections.
  double wtME, wtMEmax;
  int loop = 0;
  do {
    wtME = 1.;
    wtMEmax = 1.;
    ++loop;

    double pT, pZ, cosTheta;
    if(optAngle==0) {

    // Isotropic angles give three-momentum.
    cosTheta = 2. * rndmPtr->flat() - 1.;
    double sinTheta = sqrt(1. - cosTheta*cosTheta);
    pT = pAbs * sinTheta;
    pZ = pAbs * cosTheta;

    // Anisotropic decay
    } else {

     pair<double, double> pr = decayAngle(pAbs);
     pT = pr.first;
     pZ = pr.second;
     cosTheta = pZ/pAbs;

    }

    if(decayer.pz() < 0.0) pZ *= -1;
    if(rndmPtr->flat() < utRatio)  pZ *= -1;

    double phi      = 2. * M_PI * rndmPtr->flat();
    double pX       = pT * cos(phi);
    double pY       = pT * sin(phi);

    // Fill four-momenta and boost them away from mother rest frame.
    prod1.p(  pX,  pY,  pZ, e1);
    prod2.p( -pX, -pY, -pZ, e2);
    prod1.bst( decayer.p(), decayer.m() );
    prod2.bst( decayer.p(), decayer.m() );

    // Matrix element for PS0 -> PS1 + V1 -> PS1 + PS2 + PS3 of form
    // cos**2(theta02) in V1 rest frame, and for PS0 -> gamma + V1
    // -> gamma + PS2 + PS3 of form sin**2(theta02).
    if (meMode == 2) {
      double p10 = decayer.p() * event[iMother].p();
      double p12 = decayer.p() * prod1.p();
      double p02 = event[iMother].p() * prod1.p();
      double s0  = pow2(event[iMother].m());
      double s1  = pow2(decayer.m());
      double s2  =  pow2(prod1.m());
      if (idSister != 22) wtME = pow2(p10 * p12 - s1 * p02);
      else wtME = s1 * (2. * p10 * p12 * p02 - s1 * p02*p02
        - s0 * p12*p12 - s2 * p10*p10 + s1 * s0 * s2);
      wtME = max( wtME, 1e-6 * s1*s1 * s0 * s2);
      wtMEmax = (p10*p10 - s1 * s0) * (p12*p12 - s1 * s2);
    } else if(meMode == 3) {
	wtMEmax=2.0;
	wtME = 1.0 + cosTheta*cosTheta;
    }

    // Break out of loop if no sensible ME weight.
    if(loop > NTRYMEWT) {
      infoPtr->errorMsg("JParticleDecays::twoBody: "
        "caught in infinite ME weight loop");
      wtME = abs(wtMEmax);
    }

  // If rejected, try again with new invariant masses.
  } while ( wtME < rndmPtr->flat() * wtMEmax );

  // Done.
  return true;

}

pair<double, double> JParticleDecays::decayAngle(double pAbs)
{
  double ptx0=gWidth;
  double pr2=pAbs*pAbs;
  double pT, pZ;
 
  // Gauss
  if(optAngle==1) {
    pT=ptx0*sqrt(-log(1.0-rndmPtr->flat()*(1.0-exp(-pr2/ptx0/ptx0))));
    pZ=sqrt(pr2-pT*pT);

  // Gauss + Isotropic
  } else if(optAngle==2) {
    //double ptx0=min(pAbs,ptx0);
    pT=ptx0*sqrt(-log(max(1.e-10,rndmPtr->flat())));
    if(pT > pAbs-0.001) {
      pZ=pAbs*(1.0-2*rndmPtr->flat());
      pT=sqrt(pr2-pZ*pZ);
    } else {
      pZ=sqrt(pr2-pT*pT);
    }

  } else if(optAngle==3) {
    double ptx=pAbs*M_PI/2;
    double expf1=1.0-exp(-ptx*ptx/ptx0/ptx0);
    pT=ptx0*sqrt(-log(1.0-rndmPtr->flat()*expf1));
    double theta1=pT/pAbs;
    pT=pAbs*sin(theta1);
    pZ=pAbs*cos(theta1);

 } else {

   double cos1=0.0;
   do {
     cos1=-1.0+2.0*rndmPtr->flat();
   } while(rndmPtr->flat() > (1.0+cos1*cos1)/2.0 );

   double sin1=sqrt(1.0-cos1*cos1);
   pT=pAbs*sin1;
   pZ=pAbs*cos1;

  }

  return make_pair(pT,pZ);
}
//--------------------------------------------------------------------------

// Do a three-body decay (except Dalitz decays).

bool JParticleDecays::threeBody(Event& event) {

  // References to the particles involved.
  Particle& decayer = event[iProd[0]];
  Particle& prod1   = event[iProd[1]];
  Particle& prod2   = event[iProd[2]];
  Particle& prod3   = event[iProd[3]];

  // Mother and sum daughter masses. Fail if too close.
  double m0      = mProd[0];
  double m1      = mProd[1];
  double m2      = mProd[2];
  double m3      = mProd[3];
  double mSum    = m1 + m2 + m3;
  double mDiff   = m0 - mSum;
  if (mDiff < mSafety) return false;

  // Kinematical limits for 2+3 mass. Maximum phase-space weight.
  double m23Min  = m2 + m3;
  double m23Max  = m0 - m1;
  double p1Max   = 0.5 * sqrtpos( (m0 - m1 - m23Min) * (m0 + m1 + m23Min)
    * (m0 + m1 - m23Min) * (m0 - m1 + m23Min) ) / m0;
  double p23Max  = 0.5 * sqrtpos( (m23Max - m2 - m3) * (m23Max + m2 + m3)
    * (m23Max + m2 - m3) * (m23Max - m2 + m3) ) / m23Max;
  double wtPSmax = 0.5 * p1Max * p23Max;

  // Begin loop over matrix-element corrections.
  double wtME, wtMEmax, wtPS, m23, p1Abs, p23Abs;
  do {
    wtME     = 1.;
    wtMEmax  = 1.;

    // Pick an intermediate mass m23 flat in the allowed range.
    do {
      m23    = m23Min + rndmPtr->flat() * mDiff;

      // Translate into relative momenta and find phase-space weight.
      p1Abs  = 0.5 * sqrtpos( (m0 - m1 - m23) * (m0 + m1 + m23)
        * (m0 + m1 - m23) * (m0 - m1 + m23) ) / m0;
      p23Abs = 0.5 * sqrtpos( (m23 - m2 - m3) * (m23 + m2 + m3)
        * (m23 + m2 - m3) * (m23 - m2 + m3) ) / m23;
      wtPS   = p1Abs * p23Abs;

    // If rejected, try again with new invariant masses.
    } while ( wtPS < rndmPtr->flat() * wtPSmax );

    // Set up m23 -> m2 + m3 isotropic in its rest frame.
    double cosTheta = 2. * rndmPtr->flat() - 1.;
    double sinTheta = sqrt(1. - cosTheta*cosTheta);
    double phi      = 2. * M_PI * rndmPtr->flat();
    double pX       = p23Abs * sinTheta * cos(phi);
    double pY       = p23Abs * sinTheta * sin(phi);
    double pZ       = p23Abs * cosTheta;
    double e2       = sqrt( m2*m2 + p23Abs*p23Abs);
    double e3       = sqrt( m3*m3 + p23Abs*p23Abs);
    prod2.p(  pX,  pY,  pZ, e2);
    prod3.p( -pX, -pY, -pZ, e3);

    // Set up m0 -> m1 + m23 isotropic in its rest frame.
    cosTheta        = 2. * rndmPtr->flat() - 1.;
    sinTheta        = sqrt(1. - cosTheta*cosTheta);
    phi             = 2. * M_PI * rndmPtr->flat();
    pX              = p1Abs * sinTheta * cos(phi);
    pY              = p1Abs * sinTheta * sin(phi);
    pZ              = p1Abs * cosTheta;
    double e1       = sqrt( m1*m1 + p1Abs*p1Abs);
    double e23      = sqrt( m23*m23 + p1Abs*p1Abs);
    prod1.p( pX, pY, pZ, e1);

    // Boost 2 + 3 to the 0 rest frame.
    Vec4 p23( -pX, -pY, -pZ, e23);
    prod2.bst( p23, m23 );
    prod3.bst( p23, m23 );

    // Matrix-element weight for omega/phi -> pi+ pi- pi0.
    if (meMode == 1) {
      double p1p2 = prod1.p() * prod2.p();
      double p1p3 = prod1.p() * prod3.p();
      double p2p3 = prod2.p() * prod3.p();
      wtME = pow2(m1 * m2 * m3) - pow2(m1 * p2p3) - pow2(m2 * p1p3)
        - pow2(m3 * p1p2) + 2. * p1p2 * p1p3 * p2p3;
      wtMEmax = pow3(m0 * m0) / 150.;

    // Effective matrix element for nu spectrum in tau -> nu + hadrons.
    } else if (meMode == 21) {
      double x1 = 2. *  prod1.e() / m0;
      wtME = x1 * (3. - 2. * x1);
      double xMax = min( 0.75, 2. * (1. - mSum / m0) );
      wtMEmax = xMax * (3. - 2. * xMax);

    // Matrix element for weak decay (only semileptonic for c and b).
    } else if ((meMode == 22 || meMode == 23) && prod1.isLepton()) {
      wtME = m0 * prod1.e() * (prod2.p() * prod3.p());
      wtMEmax = min( pow4(m0) / 16., m0 * (m0 - m1 - m2) * (m0 - m1 - m3)
        * (m0 - m2 - m3) );

    // Effective matrix element for weak decay to hadrons (B -> D, D -> K).
    } else if (meMode == 22 || meMode == 23) {
      double x1 = 2. * prod1.pAbs() / m0;
      wtME = x1 * (3. - 2. * x1);
      double xMax = min( 0.75, 2. * (1. - mSum / m0) );
      wtMEmax = xMax * (3. - 2. * xMax);

    // Effective matrix element for gamma spectrum in B -> gamma + hadrons.
    } else if (meMode == 31) {
      double x1 = 2. * prod1.e() / m0;
      wtME = pow3(x1);
      double x1Max = 1. - pow2(mSum / m0);
      wtMEmax = pow3(x1Max);

    // Matrix-element weight for "onium" -> g + g + g or gamma + g + g.
    } else if (meMode == 92) {
      double x1 = 2. * prod1.e() / m0;
      double x2 = 2. * prod2.e() / m0;
      double x3 = 2. * prod3.e() / m0;
      wtME = pow2( (1. - x1) / (x2 * x3) ) + pow2( (1. - x2) / (x1 * x3) )
        + pow2( (1. - x3) / (x1 * x2) );
      wtMEmax = 2.;
      // For gamma + g + g require minimum mass for g + g system.
      if (prod1.id() == 22 && sqrt(1. - x1) * m0 < 2. * stopMass) wtME = 0.;
      if (prod2.id() == 22 && sqrt(1. - x2) * m0 < 2. * stopMass) wtME = 0.;
      if (prod3.id() == 22 && sqrt(1. - x3) * m0 < 2. * stopMass) wtME = 0.;
    }

  // If rejected, try again with new invariant masses.
  } while ( wtME < rndmPtr->flat() * wtMEmax );

  // Boost 1 + 2 + 3 to the current frame.
  prod1.bst( decayer.p(), decayer.m() );
  prod2.bst( decayer.p(), decayer.m() );
  prod3.bst( decayer.p(), decayer.m() );

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Do a multibody decay using the M-generator algorithm.

bool JParticleDecays::mGenerator(Event& event) {

  // Mother and sum daughter masses. Fail if too close or inconsistent.
  double m0      = mProd[0];
  double mSum    = mProd[1];
  for (int i = 2; i <= mult; ++i) mSum += mProd[i];
  double mDiff   = m0 - mSum;
  if (mDiff < mSafety) return false;

  // Begin setup of intermediate invariant masses.
  mInv.resize(0);
  for (int i = 0; i <= mult; ++i) mInv.push_back( mProd[i]);

  // Calculate the maximum weight in the decay.
  double wtPS, wtME, wtMEmax;
  double wtPSmax = 1. / WTCORRECTION[mult];
  double mMax    = mDiff + mProd[mult];
  double mMin    = 0.;
  for (int i = mult - 1; i > 0; --i) {
    mMax        += mProd[i];
    mMin        += mProd[i+1];
    double mNow  = mProd[i];
    wtPSmax     *= 0.5 * sqrtpos( (mMax - mMin - mNow) * (mMax + mMin + mNow)
                 * (mMax + mMin - mNow) * (mMax - mMin + mNow) ) / mMax;
  }

  // Begin loop over matrix-element corrections.
  do {
    wtME    = 1.;
    wtMEmax = 1.;

    // Begin loop to find the set of intermediate invariant masses.
    do {
      wtPS  = 1.;

      // Find and order random numbers in descending order.
      rndmOrd.resize(0);
      rndmOrd.push_back(1.);
      for (int i = 1; i < mult - 1; ++i) {
        double rndm = rndmPtr->flat();
        rndmOrd.push_back(rndm);
        for (int j = i - 1; j > 0; --j) {
          if (rndm > rndmOrd[j]) swap( rndmOrd[j], rndmOrd[j+1] );
          else break;
        }
      }
      rndmOrd.push_back(0.);

      // Translate into intermediate masses and find weight.
      for (int i = mult - 1; i > 0; --i) {
        mInv[i] = mInv[i+1] + mProd[i] + (rndmOrd[i-1] - rndmOrd[i]) * mDiff;
        wtPS   *= 0.5 * sqrtpos( (mInv[i] - mInv[i+1] - mProd[i])
          * (mInv[i] + mInv[i+1] + mProd[i]) * (mInv[i] + mInv[i+1] - mProd[i])
          * (mInv[i] - mInv[i+1] + mProd[i]) ) / mInv[i];
      }

    // If rejected, try again with new invariant masses.
    } while ( wtPS < rndmPtr->flat() * wtPSmax );

    // Perform two-particle decays in the respective rest frame.
    pInv.resize(mult + 1);
    for (int i = 1; i < mult; ++i) {
      double pAbs = 0.5 * sqrtpos( (mInv[i] - mInv[i+1] - mProd[i])
        * (mInv[i] + mInv[i+1] + mProd[i]) * (mInv[i] + mInv[i+1] - mProd[i])
        * (mInv[i] - mInv[i+1] + mProd[i]) ) / mInv[i];

      // Isotropic angles give three-momentum.
      double cosTheta = 2. * rndmPtr->flat() - 1.;
      double sinTheta = sqrt(1. - cosTheta*cosTheta);
      double phi      = 2. * M_PI * rndmPtr->flat();
      double pX       = pAbs * sinTheta * cos(phi);
      double pY       = pAbs * sinTheta * sin(phi);
      double pZ       = pAbs * cosTheta;

      // Calculate energies, fill four-momenta.
      double eHad     = sqrt( mProd[i]*mProd[i] + pAbs*pAbs);
      double eInv     = sqrt( mInv[i+1]*mInv[i+1] + pAbs*pAbs);
      event[iProd[i]].p( pX, pY, pZ, eHad);
      pInv[i+1].p( -pX, -pY, -pZ, eInv);
    }

    // Boost decay products to the mother rest frame.
    event[iProd[mult]].p( pInv[mult] );
    for (int iFrame = mult - 1; iFrame > 1; --iFrame)
      for (int i = iFrame; i <= mult; ++i)
        event[iProd[i]].bst( pInv[iFrame], mInv[iFrame]);

    // Effective matrix element for nu spectrum in tau -> nu + hadrons.
    if (meMode == 21 && event[iProd[1]].isLepton()) {
      double x1 = 2. * event[iProd[1]].e() / m0;
      wtME = x1 * (3. - 2. * x1);
      double xMax = min( 0.75, 2. * (1. - mSum / m0) );
      wtMEmax = xMax * (3. - 2. * xMax);

    // Effective matrix element for weak decay (only semileptonic for c and b).
    // Particles 4 onwards should be made softer explicitly?
    } else if ((meMode == 22 || meMode == 23) && event[iProd[1]].isLepton()) {
      Vec4 pRest = event[iProd[3]].p();
      for (int i = 4; i <= mult; ++i) pRest += event[iProd[i]].p();
      wtME = m0 * event[iProd[1]].e() * (event[iProd[2]].p() * pRest);
      for (int i = 4; i <= mult; ++i) wtME
        *= exp(- event[iProd[i]].pAbs2() / pow2(sigmaSoft) );
      wtMEmax = pow4(m0) / 16.;

    // Effective matrix element for weak decay to hadrons (B -> D, D -> K).
    } else if (meMode == 22 || meMode == 23) {
      double x1 = 2. * event[iProd[1]].pAbs() / m0;
      wtME = x1 * (3. - 2. * x1);
      double xMax = min( 0.75, 2. * (1. - mSum / m0) );
      wtMEmax = xMax * (3. - 2. * xMax);

    // Effective matrix element for gamma spectrum in B -> gamma + hadrons.
    } else if (meMode == 31) {
      double x1 = 2. * event[iProd[1]].e() / m0;
      wtME = pow3(x1);
      double x1Max = 1. - pow2(mSum / m0);
      wtMEmax = pow3(x1Max);
    }

  // If rejected, try again with new invariant masses.
  } while ( wtME < rndmPtr->flat() * wtMEmax );

  // Boost decay products to the current frame.
  pInv[1].p( event[iProd[0]].p() );
  for (int i = 1; i <= mult; ++i) event[iProd[i]].bst( pInv[1], mInv[1] );

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Select mass of lepton pair in a Dalitz decay.

bool JParticleDecays::dalitzMass() {

  // Mother and sum daughter masses.
  double mSum1 = 0;
  for (int i = 1; i <= mult - 2; ++i) mSum1 += mProd[i];
  if (meMode == 13) mSum1 *= MSAFEDALITZ;
  double mSum2 = MSAFEDALITZ * (mProd[mult -1] + mProd[mult]);
  double mDiff = mProd[0] - mSum1 - mSum2;

  // Fail if too close or inconsistent.
  if (mDiff < mSafety) return false;
  if (idProd[mult - 1] + idProd[mult] != 0
    || mProd[mult - 1] != mProd[mult]) {
    infoPtr->errorMsg("Error in JParticleDecays::dalitzMass:"
    " inconsistent flavour/mass assignments");
    return false;
  }
  if ( meMode == 13 && (idProd[1] + idProd[2] != 0
    || mProd[1] != mProd[2]) ) {
    infoPtr->errorMsg("Error in JParticleDecays::dalitzMass:"
    " inconsistent flavour/mass assignments");
    return false;
  }

  // Case 1: one Dalitz pair.
  if (meMode == 11 || meMode == 12) {

    // Kinematical limits for gamma* squared mass.
    double sGamMin = pow2(mSum2);
    double sGamMax = pow2(mProd[0] - mSum1);
    // Select virtual gamma squared mass. Guessed form for meMode == 12.
    double sGam, wtGam;
    int loop = 0;
    do {
      if (++loop > NTRYDALITZ) return false;
      sGam = sGamMin * pow( sGamMax / sGamMin, rndmPtr->flat() );
      wtGam = (1. + 0.5 * sGamMin / sGam) *  sqrt(1. - sGamMin / sGam)
        * pow3(1. - sGam / sGamMax) * sRhoDal * (sRhoDal + wRhoDal)
        / ( pow2(sGam - sRhoDal) + sRhoDal * wRhoDal );
    } while ( wtGam < rndmPtr->flat() );

    // Store results in preparation for doing a one-less-body decay.
    --mult;
    mProd[mult] = sqrt(sGam);

  // Case 2: two Dalitz pairs.
  } else {

    // Kinematical limits for 1 + 2 and 3 + 4 gamma* masses.
    double s0 = pow2(mProd[0]);
    double s12Min = pow2(mSum1);
    double s12Max = pow2(mProd[0] - mSum2);
    double s34Min = pow2(mSum2);
    double s34Max = pow2(mProd[0] - mSum1);

    // Select virtual gamma squared masses. Guessed form for meMode == 13.
    double s12, s34, wt12, wt34, wtPAbs, wtAll;
    int loop = 0;
    do {
      if (++loop > NTRYDALITZ) return false;
      s12 = s12Min * pow( s12Max / s12Min, rndmPtr->flat() );
      wt12 = (1. + 0.5 * s12Min / s12) *  sqrt(1. - s12Min / s12)
        * sRhoDal * (sRhoDal + wRhoDal)
        / ( pow2(s12 - sRhoDal) + sRhoDal * wRhoDal );
      s34 = s34Min * pow( s34Max / s34Min, rndmPtr->flat() );
      wt34 = (1. + 0.5 * s34Min / s34) *  sqrt(1. - s34Min / s34)
        * sRhoDal * (sRhoDal + wRhoDal)
        / ( pow2(s34 - sRhoDal) + sRhoDal * wRhoDal );
      wtPAbs = sqrtpos( pow2(1. - (s12 + s34)/ s0)
        - 4. * s12 * s34 / (s0 * s0) );
      wtAll = wt12 * wt34 * pow3(wtPAbs);
      if (wtAll > 1.) infoPtr->errorMsg(
        "Error in JParticleDecays::dalitzMass: weight > 1");
    } while (wtAll < rndmPtr->flat());

    // Store results in preparation for doing a two-body decay.
    mult = 2;
    mProd[1] = sqrt(s12);
    mProd[2] = sqrt(s34);
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Do kinematics of gamma* -> l- l+ in Dalitz decay.

bool JParticleDecays::dalitzKinematics(Event& event) {

  // Restore multiplicity.
  int nDal = (meMode < 13) ? 1 : 2;
  mult += nDal;

  // Loop over one or two lepton pairs.
  for (int iDal = 0; iDal < nDal; ++iDal) {

    // References to the particles involved.
    Particle& decayer = event[iProd[0]];
    Particle& prodA = (iDal == 0) ? event[iProd[mult - 1]]
      : event[iProd[1]];
    Particle& prodB = (iDal == 0) ? event[iProd[mult]]
      : event[iProd[2]];

    // Reconstruct required rotations and boosts backwards.
    Vec4 pDec    = decayer.p();
    int  iGam    = (meMode < 13) ? mult - 1 : 2 - iDal;
    Vec4 pGam    = event[iProd[iGam]].p();
    pGam.bstback( pDec, decayer.m() );
    double phiGam = pGam.phi();
    pGam.rot( 0., -phiGam);
    double thetaGam = pGam.theta();
    pGam.rot( -thetaGam, 0.);

    // Masses and phase space in gamma* rest frame.
    double mGam     = (meMode < 13) ? mProd[mult - 1] : mProd[2 - iDal];
    double mA       = prodA.m();
    double mB       = prodB.m();
    double mGamMin  = MSAFEDALITZ * (mA + mB);
    double mGamRat  = pow2(mGamMin / mGam);
    double pGamAbs  = 0.5 * sqrtpos( (mGam - mA - mB) * (mGam + mA + mB) );

    // Set up decay in gamma* rest frame, reference along +z axis.
    double cosTheta, cos2Theta;
    do {
      cosTheta      = 2. * rndmPtr->flat() - 1.;
      cos2Theta     = cosTheta * cosTheta;
    } while ( 1. + cos2Theta + mGamRat * (1. - cos2Theta)
      < 2. * rndmPtr->flat() );
    double sinTheta = sqrt(1. - cosTheta*cosTheta);
    double phi      = 2. * M_PI * rndmPtr->flat();
    double pX       = pGamAbs * sinTheta * cos(phi);
    double pY       = pGamAbs * sinTheta * sin(phi);
    double pZ       = pGamAbs * cosTheta;
    double eA       = sqrt( mA*mA + pGamAbs*pGamAbs);
    double eB       = sqrt( mB*mB + pGamAbs*pGamAbs);
    prodA.p(  pX,  pY,  pZ, eA);
    prodB.p( -pX, -pY, -pZ, eB);

    // Boost to lab frame.
    prodA.bst( pGam, mGam);
    prodB.bst( pGam, mGam);
    prodA.rot( thetaGam, phiGam);
    prodB.rot( thetaGam, phiGam);
    prodA.bst( pDec, decayer.m() );
    prodB.bst( pDec, decayer.m() );
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

} // end namespace Pythia8
