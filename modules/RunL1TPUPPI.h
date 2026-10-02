/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//------------------------------------------------------------------------------

#ifndef RunL1TPUPPI_h
#define RunL1TPUPPI_h

/** \class RunL1TPUPPI
 *
 *  Phase-2 Level-1 PUPPI, following the CMS correlator emulator
 *  (L1Trigger/Phase2L1ParticleFlow, LinPuppiEmulator::linpuppi_ref,
 *  fwdlinpuppi_ref and sum2puppiPt_ref).
 *
 *  Unlike RunPUPPI, which is the offline PUPPI algorithm, this module does not
 *  estimate a per-event median and RMS of the pile-up alpha distribution. At
 *  L1 the charged reference is only populated by tracks above the ~2 GeV track
 *  finder threshold, which leaves the median sample far too sparse to be
 *  estimated event by event. The firmware instead compares alpha to tuned
 *  constants:
 *
 *    alpha = log( sum_j min(pt_j, ptMax)^2 / max(dR_ij, drMin)^2 )   (0 if empty)
 *    x2a   = clamp( alphaSlope * (alpha - alphaZero), +-alphaCrop )
 *    x2pt  = ptSlope * (pt - ptZero)
 *    w     = 1 / ( 1 + exp( -(x2a + x2pt - prior) ) )
 *
 *  The sum runs over primary-vertex tracks in regions with UseTracks true, and
 *  over the other neutral candidates (self excluded) where it is false, which
 *  is how the no-tracker and HF regions are handled in the firmware.
 *
 *  Because the alpha term is cropped and a pt term is added, a hard neutral
 *  with no nearby primary-vertex track is suppressed but not set to exactly
 *  zero, which is the main behavioural difference with respect to RunPUPPI.
 *
 *  Charged candidates follow linpuppi_chs_ref: weight 1 from the primary
 *  vertex, 0 otherwise, with no pt threshold. The firmware exempts muons from
 *  the vertex requirement; here that is handled the way the existing cards
 *  already do it, by filtering leptons out before this module and merging them
 *  back afterwards.
 *
 *  As in RunPUPPI the output four-vectors are not scaled by the weight; the
 *  weight is stored in Candidate::puppiW. Candidates whose weighted pt falls
 *  below PtCutBin are emitted with puppiW = 0.
 *
 */

#include "classes/DelphesModule.h"
#include <vector>

class TObjArray;
class TIterator;

class RunL1TPUPPI: public DelphesModule
{

public:
  RunL1TPUPPI();
  ~RunL1TPUPPI();

  void Init();
  void Process();
  void Finish();

  // one entry per |eta| region
  struct Region
  {
    double etaMin, etaMax;
    bool useTracks;
    double coneSize, coneSizeMin, ptMax;
    double ptCut;
    double ptSlope, ptSlopePhoton;
    double ptZero, ptZeroPhoton;
    double alphaSlope, alphaZero, alphaCrop;
    double prior, priorPhoton;
  };

  // light-weight stand-in for a candidate in the alpha sums. idx identifies the
  // candidate within its own collection so that the no-tracker regions can drop
  // the self term; it is never negative, so a skip of -1 excludes nothing.
  struct Seed
  {
    double pt, eta, phi;
    Int_t idx;
  };

private:
  Int_t GetRegion(Double_t eta) const;
  Double_t Alpha(const std::vector<Seed> &seeds, Double_t eta, Double_t phi,
    const Region &region, Int_t skip) const;
  Double_t Weight(Double_t alpha, Double_t pt, Bool_t isPhoton, const Region &region) const;

  TIterator *fItTrackInputArray; //!
  TIterator *fItNeutralInputArray; //!
  TIterator *fItPVInputArray; //!

  const TObjArray *fTrackInputArray; //!
  const TObjArray *fNeutralInputArray; //!
  const TObjArray *fPVInputArray; //!

  Double_t fDeltaZMax; //!

  std::vector<Region> fRegions; //!

  TObjArray *fOutputArray; //!
  TObjArray *fOutputTrackArray; //!
  TObjArray *fOutputNeutralArray; //!

  ClassDef(RunL1TPUPPI, 1)
};

#endif
