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

/** \class PrimaryVertexFinder
 *
 *  \see modules/PrimaryVertexFinder.h
 *
 */

#include "modules/PrimaryVertexFinder.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

#include "TMath.h"
#include "TObjArray.h"

#include <stdexcept>
#include <vector>

using namespace std;

//------------------------------------------------------------------------------

PrimaryVertexFinder::PrimaryVertexFinder() :
  fItVertexInputArray(0),
  fItTrackInputArray(0),
  fTrackInputArray(0)
{
}

//------------------------------------------------------------------------------

PrimaryVertexFinder::~PrimaryVertexFinder() {}

//------------------------------------------------------------------------------

void PrimaryVertexFinder::Init()
{
  fMethod = GetString("Method", "FastHisto");
  if(fMethod != "FastHisto" && fMethod != "SumPT2" && fMethod != "Signal")
  {
    throw runtime_error("PrimaryVertexFinder: Method must be FastHisto, SumPT2 or Signal");
  }

  fVertexInputArray = ImportArray(GetString("VertexInputArray", "PileUpMerger/vertices"));
  fItVertexInputArray = fVertexInputArray->MakeIterator();

  if(fMethod == "FastHisto")
  {
    fTrackInputArray = ImportArray(GetString("TrackInputArray", "HCal/eflowTracks"));
    fItTrackInputArray = fTrackInputArray->MakeIterator();
  }

  // L1Trigger/VertexFinder defaults, converted from cm to mm
  fMinTrackPt = GetDouble("MinTrackPt", 2.0);
  fMaxTrackPt = GetDouble("MaxTrackPt", 127.0);
  fHistogramMin = GetDouble("HistogramMin", -204.6912512);
  fHistogramMax = GetDouble("HistogramMax", 204.6912512);
  fBinWidth = GetDouble("BinWidth", 1.5991504);
  fWindowSize = GetInt("WindowSize", 3);

  if(fBinWidth <= 0.0 || fHistogramMax <= fHistogramMin || fWindowSize < 1)
  {
    throw runtime_error("PrimaryVertexFinder: invalid histogram configuration");
  }

  fOutputArray = ExportArray(GetString("OutputArray", "vertices"));
}

//------------------------------------------------------------------------------

void PrimaryVertexFinder::Finish()
{
  if(fItVertexInputArray) delete fItVertexInputArray;
  if(fItTrackInputArray) delete fItTrackInputArray;
}

//------------------------------------------------------------------------------

Candidate *PrimaryVertexFinder::FindSignal()
{
  Candidate *vertex;
  fItVertexInputArray->Reset();
  while((vertex = static_cast<Candidate *>(fItVertexInputArray->Next())))
  {
    if(!vertex->IsPU) return vertex;
  }
  return 0;
}

//------------------------------------------------------------------------------

Candidate *PrimaryVertexFinder::FindBySumPT2()
{
  Candidate *vertex, *best = 0;
  fItVertexInputArray->Reset();
  while((vertex = static_cast<Candidate *>(fItVertexInputArray->Next())))
  {
    if(!best || vertex->SumPT2 > best->SumPT2) best = vertex;
  }
  return best;
}

//------------------------------------------------------------------------------

// fastHisto, following L1Trigger/VertexFinder/src/VertexFinder.cc
Candidate *PrimaryVertexFinder::FindFastHisto()
{
  Int_t nBins = Int_t((fHistogramMax - fHistogramMin) / fBinWidth);
  if(nBins < fWindowSize) return 0;

  vector<Double_t> hist(nBins, 0.0);
  Candidate *candidate, *particle;

  fItTrackInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItTrackInputArray->Next())))
  {
    Double_t pt = candidate->Momentum.Pt();
    if(pt < fMinTrackPt) continue;
    if(pt > fMaxTrackPt) pt = fMaxTrackPt; // VxMaxTrackPtBehavior = 1, saturate

    particle = static_cast<Candidate *>(candidate->GetCandidates()->At(0));
    Double_t z = particle->Position.Z();
    if(z < fHistogramMin || z >= fHistogramMax) continue;

    hist[Int_t((z - fHistogramMin) / fBinWidth)] += pt;
  }

  // sliding window of fWindowSize bins with the largest pt sum
  Double_t sum = 0.0;
  for(Int_t i = 0; i < fWindowSize; ++i) sum += hist[i];

  Double_t bestSum = sum;
  Int_t bestBin = 0;
  for(Int_t i = fWindowSize; i < nBins; ++i)
  {
    sum += hist[i] - hist[i - fWindowSize];
    if(sum > bestSum)
    {
      bestSum = sum;
      bestBin = i - fWindowSize + 1;
    }
  }
  if(bestSum <= 0.0) return 0;

  // pt-weighted mean z of the tracks inside the winning window
  Double_t zLow = fHistogramMin + bestBin * fBinWidth;
  Double_t zHigh = zLow + fWindowSize * fBinWidth;
  Double_t sumZ = 0.0, sumPt = 0.0;

  fItTrackInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItTrackInputArray->Next())))
  {
    Double_t pt = candidate->Momentum.Pt();
    if(pt < fMinTrackPt) continue;
    if(pt > fMaxTrackPt) pt = fMaxTrackPt;

    particle = static_cast<Candidate *>(candidate->GetCandidates()->At(0));
    Double_t z = particle->Position.Z();
    if(z < zLow || z >= zHigh) continue;

    sumZ += pt * z;
    sumPt += pt;
  }
  if(sumPt <= 0.0) return 0;

  Candidate *vertex = GetFactory()->NewCandidate();
  vertex->Position.SetXYZT(0.0, 0.0, sumZ / sumPt, 0.0);
  vertex->SumPT2 = bestSum;
  return vertex;
}

//------------------------------------------------------------------------------

void PrimaryVertexFinder::Process()
{
  Candidate *chosen = 0;
  if(fMethod == "FastHisto")
    chosen = FindFastHisto();
  else if(fMethod == "SumPT2")
    chosen = FindBySumPT2();
  else
    chosen = FindSignal();

  // No vertex could be reconstructed: fall back to the signal vertex rather
  // than leaving the chain without a primary vertex at all.
  Bool_t synthetic = (fMethod == "FastHisto" && chosen);
  if(!chosen) chosen = FindSignal();
  if(!chosen) return;

  Candidate *output = static_cast<Candidate *>(chosen->Clone());
  output->IsPU = 0;
  fOutputArray->Add(output);

  Candidate *vertex;
  fItVertexInputArray->Reset();
  while((vertex = static_cast<Candidate *>(fItVertexInputArray->Next())))
  {
    if(!synthetic && vertex == chosen) continue; // already emitted as the PV
    Candidate *other = static_cast<Candidate *>(vertex->Clone());
    other->IsPU = 1;
    fOutputArray->Add(other);
  }
}

//------------------------------------------------------------------------------
