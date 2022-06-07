// ParticleDecays.h is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the classes to perform a particle decay.
// DecayHandler: base class for external handling of decays.
// ParticleDecays: decay a particle.

#ifndef jam2_interaction_JParticleDecays_H
#define jam2_interaction_JParticleDecays_H

#include <Pythia8/Basics.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/Info.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PythiaStdlib.h>
#include <Pythia8/Settings.h>
#include <Pythia8/TimeShower.h>
#include <Pythia8/TauDecays.h>
//#include <jam2/interaction/HadronDecay.h>

namespace jam2 {

//==========================================================================

// The ParticleDecays class contains the routines to decay a particle.

class JParticleDecays {

public:

  // Constructor.
  JParticleDecays() {}

  // Initialize: store pointers and find settings
  void init(Pythia8::Info* infoPtrIn, Pythia8::Settings& settings,
    Pythia8::ParticleData* particleDataPtrIn, Pythia8::Rndm* rndmPtrIn);

  // Perform a decay of a single particle.
  bool decay(int iDec,bool schannel, Pythia8::Event& event);
  bool decay2( int iDec, bool schan, Pythia8::Event& event,
	double totwid, double* pWidth, int nChannel);
  double BWmass(int id, double mtot);
  void setOptAngle(int i) {optAngle=i;}

private:

  // Constants: could only be changed in the code itself.
  static const int    NTRYDECAY, NTRYPICK, NTRYMEWT, NTRYDALITZ;
  static const double MSAFEDALITZ, WTCORRECTION[11];

  bool sChannel; // resonance were prodeced by s-channel reaction.

  // Pointer to various information on the generation.
  Pythia8::Info*         infoPtr;

  // Pointer to the particle data table.
  Pythia8::ParticleData* particleDataPtr;

  // Pointer to the random number generator.
  Pythia8::Rndm*         rndmPtr;

  // Initialization data, read from Settings.
  bool   limitTau0, limitTau, limitRadius, limitCylinder, limitDecay,
         mixB, doFSRinDecays, doGammaRad;
  double mSafety, tau0Max, tauMax, rMax, xyMax, zMax, xBdMix, xBsMix,
         sigmaSoft, multIncrease, multIncreaseWeak, multRefMass, multGoffset,
         colRearrange, stopMass, sRhoDal, wRhoDal;

  // Multiplicity. Decay products positions and masses.
  int    idDec, meMode, mult;
  double scale;
  double utRatio,utRatiot, utRatios, gWidth;
  int optAngle;
  std::vector<int>    iProd, idProd, motherProd, cols, acols, idPartons;
  std::vector<double> mProd, mInv, rndmOrd;
  std::vector<Pythia8::Vec4>   pInv, pProd;

  // Pointer to particle data for currently decaying particle
  Pythia8::ParticleDataEntry* decDataPtr;

  // Check whether a decay is allowed, given the upcoming decay vertex.
  bool checkVertex(Pythia8::Particle& decayer);

  // Do a one-body decay.
  bool oneBody(Pythia8::Event& event);

  // Do a two-body decay;
  bool twoBody(Pythia8::Event& event);

  // Do a three-body decay;
  bool threeBody(Pythia8::Event& event);

  // Do a multibody decay using the M-generator algorithm.
  bool mGenerator(Pythia8::Event& event);

  // Select mass of lepton pair in a Dalitz decay.
  bool dalitzMass();

  // Do kinematics of gamma* -> l- l+ in Dalitz decay.
  bool dalitzKinematics(Pythia8::Event& event);

  // Anisotropic resonance decay angle.
  std::pair<double, double> decayAngle(double pAbs);

};

//==========================================================================

} // end namespace jam2

#endif // ParticleDecays_H
