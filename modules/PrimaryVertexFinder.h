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
 *               above MinTrackPt, take the sliding window of WindowSize bins
 *               with the largest sum, and return the pt-weighted mean z of the
 *               tracks in it. Note that Delphes tracks carry the unsmeared
 *               production z of their generator particle, so the binning does
 *               almost no work here and the result is close to "the vertex with
 *               the largest track-pt sum above MinTrackPt". That is the
 *               selection bias we are after; the z resolution is optimistic.
 *    SumPT2     the vertex with the largest SumPT2, as computed by PileUpMerger
 *               over its charged generator particles. A perfect offline-style
 *               vertex finder, with no track thresholds or acceptance.
 *    Signal     the generator signal vertex, i.e. what Delphes does today.
 *               Provided so the module can be shown to be a no-op.
 *
 *  Output is the full vertex collection with exactly one vertex flagged
 *  IsPU = 0 and placed first, so TrackPileUpSubtractor, RunPUPPI and
 *  RunL1TPUPPI all pick it up with no changes. RunPUPPI also uses
 *  GetEntries() of this array as the vertex multiplicity behind its neutral pt
 *  threshold; that count is preserved for SumPT2 and Signal, and is one larger
 *  for FastHisto, whose fitted vertex is a new object rather than one of the
 *  inputs. The shift is 0.7 of one vertex out of ~200 and is negligible there.
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
