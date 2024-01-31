/*
 * Copyright (c) 2014-2023 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// Define BaseClass_t
#include "k4FWCore/BaseClass.h"

#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCRecoTrackerHitPlaneAssociationCollection.h"
#include "edm4hep/MutableTrack.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

#include "GaudiAlg/Transformer.h"

#include "IActsGeoSvc.h"
#include "k4Interface/IGeoSvc.h"

#include "TrackFitterFunction.h"

using edm4hep::MCParticleCollection;
using edm4hep::MCRecoTrackerHitPlaneAssociationCollection;
using edm4hep::TrackCollection;
using edm4hep::TrackerHitPlaneCollection;
using edm4hep::EventHeaderCollection;

    class TruthTracking final
    : public Gaudi::Functional::MultiTransformer<std::tuple<TrackCollection>(
                                                     const EventHeaderCollection&, const MCParticleCollection&,
                                                     const TrackerHitPlaneCollection&,
                                                     const MCRecoTrackerHitPlaneAssociationCollection&),
                                                 BaseClass_t> {
public:
  TruthTracking(const std::string& name, ISvcLocator* svcLoc);
  std::tuple<TrackCollection> operator()(const EventHeaderCollection& evtHeader, const MCParticleCollection& mcp,
                                         const TrackerHitPlaneCollection&                  hits,
                                         const MCRecoTrackerHitPlaneAssociationCollection& assoc) const override;
  StatusCode                  initialize() override;

private:
  SmartIF<IGeoSvc>     m_geoSvc;
  SmartIF<IActsGeoSvc> m_actsGeoSvc;

  Acts::GeometryContext      m_geoCtx;
  Acts::MagneticFieldContext m_fieldCtx;
  Acts::CalibrationContext   m_calibCtx;

  // FIXME: I guess unique would be good enough
  std::shared_ptr<TrackFitterFunction::GeneralFitterOptions> m_options;

  std::shared_ptr<TrackFitterFunction> m_fitter;
};