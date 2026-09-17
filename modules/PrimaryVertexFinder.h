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

#ifndef PrimaryVertexFinder_h
#define PrimaryVertexFinder_h

/** \class PrimaryVertexFinder
 *
 *  Chooses which vertex the rest of the chain should treat as the primary one.
 *
 *  PileUpMerger emits the generator signal vertex first and flags every pile-up
 *  vertex with IsPU = 1, so downstream modules that look for the vertex with
 *  IsPU == 0 (TrackPileUpSubtractor) or simply take the first entry (RunPUPPI)
 *  always get the signal vertex, however soft it happens to be. In data the
 *  primary vertex is instead whichever vertex the vertex finder likes best,
 *  which in a minimum-bias event is the hardest of many collisions and not a
 *  typical one. This module reproduces that choice.
 *
 *  Method:
 *    FastHisto  L1-style, following L1Trigger/VertexFinder (fastHisto): fill a
 *               z0 histogram with track pt (saturated at MaxTrackPt) for tracks
 *               above MinTrackPt and take the sliding window of WindowSize bins
 *               with the largest sum. With the CMS defaults the window is 4.8 mm
 *               wide, so at high pile-up it usually spans several interactions.
 *               The primary vertex is the interaction contributing the most
 *               saturated track pt inside the window (reco-to-truth vertex
 *               matching through generator ancestry), falling back to the vertex
 *               nearest the window's pt-weighted mean z. Delphes tracks carry
 *               the unsmeared production z of their generator particle, so the
 *               vertex z resolution is optimistic; the selection bias, and the
 *               merging of nearby interactions, are what this reproduces.
 *    SumPT2     the vertex with the largest SumPT2, as computed by PileUpMerger
 *               over its charged generator particles. A perfect offline-style
 *               vertex finder, with no track thresholds or acceptance.
 *    Signal     the generator signal vertex, i.e. what Delphes does today.
 *               Provided so the module can be shown to be a no-op.
 *
 *  Output is the full vertex collection, cloned, with the chosen vertex placed
 *  first and flagged IsPU = 0 and every other vertex flagged IsPU = 1, so
 *  TrackPileUpSubtractor, RunPUPPI and RunL1TPUPPI all pick it up with no
 *  configuration beyond their vertex input array. The number of vertices is
 *  preserved, which matters because RunPUPPI uses it as the vertex multiplicity
 *  behind its neutral pt threshold. The chosen vertex is always one of the
 *  inputs and keeps its constituents, which TrackPileUpSubtractor uses to tell
 *  the primary interaction's tracks from pile-up once the signal vertex is no
 *  longer the primary one.
 *
 *  Distances are in mm, the Delphes convention. The CMS defaults quoted in
 *  L1Trigger/VertexFinder are in cm and have been converted.
 *
 */

#include "classes/DelphesModule.h"
#include <string>

class TObjArray;
class TIterator;
class Candidate;

class PrimaryVertexFinder: public DelphesModule
{

public:
  PrimaryVertexFinder();
  ~PrimaryVertexFinder();

  void Init();
  void Process();
  void Finish();

private:
  Candidate *FindFastHisto();
  Candidate *FindBySumPT2();
  Candidate *FindSignal();

  std::string fMethod; //!

  Double_t fMinTrackPt; //!
  Double_t fMaxTrackPt; //!
  Double_t fHistogramMin; //!
  Double_t fHistogramMax; //!
  Double_t fBinWidth; //!
  Int_t fWindowSize; //!

  TIterator *fItVertexInputArray; //!
  TIterator *fItTrackInputArray; //!

  const TObjArray *fVertexInputArray; //!
  const TObjArray *fTrackInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(PrimaryVertexFinder, 1)
};

#endif
