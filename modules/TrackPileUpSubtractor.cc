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

/** \class TrackPileUpSubtractor
 *
 *  Subtract pile-up contribution from tracks.
 *
 *  A track is flagged as pile-up (IsRecoPU = 1) when it does not come from the
 *  primary interaction and its production z lies outside the z resolution
 *  around the primary vertex. Tracks from the primary interaction are never
 *  flagged, however displaced.
 *
 *  The primary vertex is the last vertex in VertexInputArray with IsPU == 0.
 *  With the PileUpMerger collection that is the generator signal vertex, the
 *  primary interaction is the set of particles with IsPU == 0, and the truth
 *  flag is used directly. A vertex finder such as PrimaryVertexFinder may
 *  instead promote a pile-up vertex; the primary interaction is then taken to
 *  be that vertex's constituents, which PileUpMerger fills with the charged
 *  particles of each interaction, and a track belongs to it when its generator
 *  ancestor is one of them. The signal interaction's tracks are then subject
 *  to the z resolution like any other pile-up.
 *
 *  Only charged particles are vertex constituents, so with a promoted vertex
 *  a track whose ancestor is neutral (a photon conversion, for example) is
 *  treated as coming from another interaction.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/TrackPileUpSubtractor.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace
{

// the generator particle a reconstructed candidate descends from
Candidate *GetGeneratorAncestor(Candidate *candidate)
{
  while(candidate->GetCandidates()->GetEntriesFast() > 0)
  {
    candidate = static_cast<Candidate *>(candidate->GetCandidates()->At(0));
  }
  return candidate;
}

} // namespace

//------------------------------------------------------------------------------

TrackPileUpSubtractor::TrackPileUpSubtractor() :
  fFormula(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

TrackPileUpSubtractor::~TrackPileUpSubtractor()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void TrackPileUpSubtractor::Init()
{
  // import input array

  fVertexInputArray = ImportArray(GetString("VertexInputArray", "PileUpMerger/vertices"));
  fItVertexInputArray = fVertexInputArray->MakeIterator();

  // read resolution formula in m
  fFormula->Compile(GetString("ZVertexResolution", "0.001"));

  fPTMin = GetDouble("PTMin", 0.);

  // import arrays with output from other modules

  ExRootConfParam param = GetParam("InputArray");
  Long_t i, size;
  const TObjArray *array;
  TIterator *iterator;

  size = param.GetSize();
  for(i = 0; i < size / 2; ++i)
  {
    array = ImportArray(param[i * 2].GetString());
    iterator = array->MakeIterator();

    fInputMap[iterator] = ExportArray(param[i * 2 + 1].GetString());
  }
}

//------------------------------------------------------------------------------

void TrackPileUpSubtractor::Finish()
{
  map<TIterator *, TObjArray *>::iterator itInputMap;
  TIterator *iterator;

  for(itInputMap = fInputMap.begin(); itInputMap != fInputMap.end(); ++itInputMap)
  {
    iterator = itInputMap->first;

    if(iterator) delete iterator;
  }

  if(fItVertexInputArray) delete fItVertexInputArray;
}

//------------------------------------------------------------------------------

void TrackPileUpSubtractor::Process()
{
  Candidate *candidate, *particle, *vertex = 0;
  map<TIterator *, TObjArray *>::iterator itInputMap;
  TIterator *iterator;
  TObjArray *array;
  Double_t z, zvtx = 0;
  Double_t pt, eta, phi, e;
  Bool_t fromPrimary;

  // find z position of primary vertex

  fItVertexInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItVertexInputArray->Next())))
  {
    if(!candidate->IsPU)
    {
      zvtx = candidate->Position.Z();
      vertex = candidate;
      // break;
    }
  }

  // Identify the primary interaction. Its constituents carry the truth IsPU
  // flag of the interaction they came from, so a promoted pile-up vertex is
  // recognisable even though its own IsPU has been reset to 0. A vertex with
  // no constituents is treated as the signal vertex, as before.
  Bool_t primaryIsSignal = kTRUE;
  set<UInt_t> primaryParticles;
  if(vertex && vertex->GetCandidates()->GetEntriesFast() > 0)
  {
    primaryIsSignal = !static_cast<Candidate *>(vertex->GetCandidates()->At(0))->IsPU;
    if(!primaryIsSignal)
    {
      TIter itConstituents(vertex->GetCandidates());
      while((particle = static_cast<Candidate *>(itConstituents.Next())))
      {
        primaryParticles.insert(particle->GetUniqueID());
      }
    }
  }

  // loop over all input arrays
  for(itInputMap = fInputMap.begin(); itInputMap != fInputMap.end(); ++itInputMap)
  {
    iterator = itInputMap->first;
    array = itInputMap->second;

    // loop over all candidates
    iterator->Reset();
    while((candidate = static_cast<Candidate *>(iterator->Next())))
    {
      particle = static_cast<Candidate *>(candidate->GetCandidates()->At(0));
      const TLorentzVector &candidateMomentum = particle->Momentum;

      eta = candidateMomentum.Eta();
      pt = candidateMomentum.Pt();
      phi = candidateMomentum.Phi();
      e = candidateMomentum.E();

      z = particle->Position.Z();

      // apply pile-up subtraction
      // assume perfect pile-up subtraction for tracks outside fZVertexResolution

      if(primaryIsSignal)
        fromPrimary = !candidate->IsPU;
      else
        fromPrimary = primaryParticles.count(GetGeneratorAncestor(candidate)->GetUniqueID()) > 0;

      if(candidate->Charge != 0 && !fromPrimary && TMath::Abs(z - zvtx) > fFormula->Eval(pt, eta, phi, e) * 1.0e3)
      {
        candidate->IsRecoPU = 1;
      }
      else
      {
        candidate->IsRecoPU = 0;
        if(candidate->Momentum.Pt() > fPTMin) array->Add(candidate);
      }
    }
  }
}
