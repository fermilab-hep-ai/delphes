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

/** \class RunL1TPUPPI
 *
 *  Phase-2 Level-1 PUPPI, following the CMS correlator emulator
 *  (L1Trigger/Phase2L1ParticleFlow, LinPuppiEmulator).
 *
 *  \see modules/RunL1TPUPPI.h for the algorithm and for how it differs from
 *  the offline PUPPI implemented in RunPUPPI.
 *
 */

#include "modules/RunL1TPUPPI.h"

#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootConfReader.h"

#include "TMath.h"
#include "TObjArray.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

using namespace std;

namespace
{

bool CompareSeedEta(const RunL1TPUPPI::Seed &a, const RunL1TPUPPI::Seed &b)
{
  return a.eta < b.eta;
}

// read a "add <Name> v0 v1 ..." list of doubles
void ReadDoubleList(DelphesModule *module, const char *name, vector<double> &out)
{
  ExRootConfParam param = module->GetParam(name);
  out.clear();
  for(Int_t i = 0; i < param.GetSize(); ++i) out.push_back(param[i].GetDouble());
}

void ReadBoolList(DelphesModule *module, const char *name, vector<bool> &out)
{
  ExRootConfParam param = module->GetParam(name);
  out.clear();
  for(Int_t i = 0; i < param.GetSize(); ++i) out.push_back(param[i].GetBool());
}

} // namespace

//------------------------------------------------------------------------------

RunL1TPUPPI::RunL1TPUPPI() :
  fItTrackInputArray(0),
  fItNeutralInputArray(0),
  fItPVInputArray(0)
{
}

//------------------------------------------------------------------------------

RunL1TPUPPI::~RunL1TPUPPI() {}

//------------------------------------------------------------------------------

void RunL1TPUPPI::Init()
{
  fTrackInputArray = ImportArray(GetString("TrackInputArray", "Calorimeter/eflowTracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();
  fNeutralInputArray = ImportArray(GetString("NeutralInputArray", "Calorimeter/eflowTowers"));
  fItNeutralInputArray = fNeutralInputArray->MakeIterator();
  fPVInputArray = ImportArray(GetString("PVInputArray", "PileUpMerger/vertices"));
  fItPVInputArray = fPVInputArray->MakeIterator();

  // track to primary vertex association: an explicit |dz| cut in mm, as in the
  // firmware, or the IsRecoPU flag already set by TrackPileUpSubtractor
  fDeltaZMax = GetDouble("DeltaZMax", -1.0);

  vector<double> etaMin, etaMax, coneSize, coneSizeMin, ptMax, ptCut;
  vector<double> ptSlope, ptSlopePhoton, ptZero, ptZeroPhoton;
  vector<double> alphaSlope, alphaZero, alphaCrop, prior, priorPhoton;
  vector<bool> useTracks;

  ReadDoubleList(this, "EtaMinBin", etaMin);
  ReadDoubleList(this, "EtaMaxBin", etaMax);
  ReadBoolList(this, "UseTracks", useTracks);
  ReadDoubleList(this, "ConeSizeBin", coneSize);
  ReadDoubleList(this, "ConeSizeMinBin", coneSizeMin);
  ReadDoubleList(this, "PtMaxBin", ptMax);
  ReadDoubleList(this, "PtCutBin", ptCut);
  ReadDoubleList(this, "PtSlopeBin", ptSlope);
  ReadDoubleList(this, "PtSlopePhotonBin", ptSlopePhoton);
  ReadDoubleList(this, "PtZeroBin", ptZero);
  ReadDoubleList(this, "PtZeroPhotonBin", ptZeroPhoton);
  ReadDoubleList(this, "AlphaSlopeBin", alphaSlope);
  ReadDoubleList(this, "AlphaZeroBin", alphaZero);
  ReadDoubleList(this, "AlphaCropBin", alphaCrop);
  ReadDoubleList(this, "PriorBin", prior);
  ReadDoubleList(this, "PriorPhotonBin", priorPhoton);

  size_t n = etaMin.size();
  if(n == 0 || etaMax.size() != n || useTracks.size() != n || coneSize.size() != n
    || coneSizeMin.size() != n || ptMax.size() != n || ptCut.size() != n
    || ptSlope.size() != n || ptSlopePhoton.size() != n || ptZero.size() != n
    || ptZeroPhoton.size() != n || alphaSlope.size() != n || alphaZero.size() != n
    || alphaCrop.size() != n || prior.size() != n || priorPhoton.size() != n)
  {
    throw runtime_error("RunL1TPUPPI: every region list must be non-empty and of the same length");
  }

  fRegions.clear();
  for(size_t i = 0; i < n; ++i)
  {
    Region region;
    region.etaMin = etaMin[i];
    region.etaMax = etaMax[i];
    region.useTracks = useTracks[i];
    region.coneSize = coneSize[i];
    region.coneSizeMin = coneSizeMin[i];
    region.ptMax = ptMax[i];
    region.ptCut = ptCut[i];
    region.ptSlope = ptSlope[i];
    region.ptSlopePhoton = ptSlopePhoton[i];
    region.ptZero = ptZero[i];
    region.ptZeroPhoton = ptZeroPhoton[i];
    region.alphaSlope = alphaSlope[i];
    region.alphaZero = alphaZero[i];
    region.alphaCrop = alphaCrop[i];
    region.prior = prior[i];
    region.priorPhoton = priorPhoton[i];
    fRegions.push_back(region);
  }

  fOutputArray = ExportArray(GetString("OutputArray", "PuppiParticles"));
  fOutputTrackArray = ExportArray(GetString("OutputArrayTracks", "puppiTracks"));
  fOutputNeutralArray = ExportArray(GetString("OutputArrayNeutrals", "puppiNeutrals"));
}

//------------------------------------------------------------------------------

void RunL1TPUPPI::Finish()
{
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItNeutralInputArray) delete fItNeutralInputArray;
  if(fItPVInputArray) delete fItPVInputArray;
}

//------------------------------------------------------------------------------

Int_t RunL1TPUPPI::GetRegion(Double_t eta) const
{
  Double_t absEta = TMath::Abs(eta);
  for(size_t i = 0; i < fRegions.size(); ++i)
  {
    if(absEta >= fRegions[i].etaMin && absEta < fRegions[i].etaMax) return Int_t(i);
  }
  return -1;
}

//------------------------------------------------------------------------------

// alpha = log( sum min(pt, ptMax)^2 / max(dR, drMin)^2 ), 0 if the cone is empty.
// The firmware evaluates log2 and folds a log(2) into alphaSlope; that is the
// same number as the natural log used here.
Double_t RunL1TPUPPI::Alpha(const vector<Seed> &seeds, Double_t eta, Double_t phi,
  const Region &region, Int_t skip) const
{
  Double_t cone2 = region.coneSize * region.coneSize;
  Double_t coneMin2 = region.coneSizeMin * region.coneSizeMin;
  Double_t sum = 0.0;

  // seeds are sorted in eta, so only the |deta| < coneSize window is scanned
  Seed low;
  low.pt = 0.0;
  low.phi = 0.0;
  low.idx = -1;
  low.eta = eta - region.coneSize;

  vector<Seed>::const_iterator it = lower_bound(seeds.begin(), seeds.end(), low, CompareSeedEta);
  for(; it != seeds.end() && it->eta < eta + region.coneSize; ++it)
  {
    if(it->idx == skip) continue;

    Double_t dEta = it->eta - eta;
    Double_t dPhi = TMath::Abs(it->phi - phi);
    if(dPhi > TMath::Pi()) dPhi = 2.0 * TMath::Pi() - dPhi;

    Double_t dr2 = dEta * dEta + dPhi * dPhi;
    if(dr2 > cone2) continue;
    if(dr2 < coneMin2) dr2 = coneMin2;

    Double_t pt = it->pt < region.ptMax ? it->pt : region.ptMax;
    sum += pt * pt / dr2;
  }

  return sum > 0.0 ? TMath::Log(sum) : 0.0;
}

//------------------------------------------------------------------------------

Double_t RunL1TPUPPI::Weight(Double_t alpha, Double_t pt, Bool_t isPhoton, const Region &region) const
{
  Double_t x2a = region.alphaSlope * (alpha - region.alphaZero);
  if(x2a > region.alphaCrop) x2a = region.alphaCrop;
  if(x2a < -region.alphaCrop) x2a = -region.alphaCrop;

  Double_t ptSlope = isPhoton ? region.ptSlopePhoton : region.ptSlope;
  Double_t ptZero = isPhoton ? region.ptZeroPhoton : region.ptZero;
  Double_t prior = isPhoton ? region.priorPhoton : region.prior;

  Double_t x2 = x2a + ptSlope * (pt - ptZero) - prior;

  if(x2 > 30.0) return 1.0;
  if(x2 < -30.0) return 0.0;
  return 1.0 / (1.0 + TMath::Exp(-x2));
}

//------------------------------------------------------------------------------

void RunL1TPUPPI::Process()
{
  Candidate *candidate, *particle, *output;

  // primary vertex, same convention as TrackPileUpSubtractor
  Double_t pvz = 0.0;
  fItPVInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItPVInputArray->Next())))
  {
    if(!candidate->IsPU) pvz = candidate->Position.Z();
  }

  // ---- charged candidates: CHS, weight 1 from the primary vertex ----------

  vector<Seed> pvTracks;

  fItTrackInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItTrackInputArray->Next())))
  {
    Bool_t fromPV;
    if(fDeltaZMax > 0.0)
    {
      particle = static_cast<Candidate *>(candidate->GetCandidates()->At(0));
      fromPV = TMath::Abs(particle->Position.Z() - pvz) < fDeltaZMax;
    }
    else
    {
      fromPV = (candidate->IsRecoPU == 0);
    }

    if(fromPV)
    {
      Seed seed;
      seed.pt = candidate->Momentum.Pt();
      seed.eta = candidate->Momentum.Eta();
      seed.phi = candidate->Momentum.Phi();
      seed.idx = Int_t(pvTracks.size());
      pvTracks.push_back(seed);
    }

    output = static_cast<Candidate *>(candidate->Clone());
    output->puppiW = fromPV ? 1.0 : 0.0;
    fOutputArray->Add(output);
    fOutputTrackArray->Add(output);
  }

  sort(pvTracks.begin(), pvTracks.end(), CompareSeedEta);

  // ---- neutral candidates ------------------------------------------------

  // the no-tracker regions build alpha out of the neutrals themselves, so the
  // whole collection has to be read before any weight can be computed
  vector<Candidate *> neutrals;
  vector<Seed> neutralSeeds;

  fItNeutralInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItNeutralInputArray->Next())))
  {
    Seed seed;
    seed.pt = candidate->Momentum.Pt();
    seed.eta = candidate->Momentum.Eta();
    seed.phi = candidate->Momentum.Phi();
    seed.idx = Int_t(neutrals.size());
    neutralSeeds.push_back(seed);
    neutrals.push_back(candidate);
  }

  vector<Seed> sortedNeutrals(neutralSeeds);
  sort(sortedNeutrals.begin(), sortedNeutrals.end(), CompareSeedEta);

  for(size_t i = 0; i < neutrals.size(); ++i)
  {
    candidate = neutrals[i];
    const Seed &seed = neutralSeeds[i];

    Double_t weight = 0.0;
    Int_t iRegion = GetRegion(seed.eta);
    if(iRegion >= 0)
    {
      const Region &region = fRegions[iRegion];
      Bool_t isPhoton = (TMath::Abs(candidate->PID) == 22);

      Double_t alpha = region.useTracks
        ? Alpha(pvTracks, seed.eta, seed.phi, region, -1)
        : Alpha(sortedNeutrals, seed.eta, seed.phi, region, seed.idx);

      weight = Weight(alpha, seed.pt, isPhoton, region);
      if(weight * seed.pt < region.ptCut) weight = 0.0;
    }

    output = static_cast<Candidate *>(candidate->Clone());
    output->puppiW = weight;
    fOutputArray->Add(output);
    fOutputNeutralArray->Add(output);
  }
}

//------------------------------------------------------------------------------
