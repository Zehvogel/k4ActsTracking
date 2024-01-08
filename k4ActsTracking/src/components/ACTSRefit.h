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

#include "Acts/Geometry/TrackingGeometry.hpp"

#include "edm4hep/MutableTrack.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

#include "Gaudi/Property.h"
#include "GaudiAlg/Transformer.h"

// Define BaseClass_t
#include "k4FWCore/BaseClass.h"

#include "IActsGeoSvc.h"
#include "k4Interface/IGeoSvc.h"

#include "DDRec/SurfaceManager.h"

#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "TrackFitterFunction.h"

// Which type of collection we are reading and writing
// using colltype_in  = edm4hep::MCParticleCollection;
// using colltype_out = edm4hep::MCParticleCollection;
using edm4hep::TrackCollection;
using edm4hep::TrackerHitPlaneCollection;

class ACTSRefit final : public Gaudi::Functional::Transformer<TrackCollection(const TrackerHitPlaneCollection&), BaseClass_t> {
public:
  ACTSRefit(const std::string& name, ISvcLocator* svcLoc);
  TrackCollection operator()(const TrackerHitPlaneCollection& input) const override;
  StatusCode initialize() override;
private:
  SmartIF<IGeoSvc> m_geoSvc;
  SmartIF<IActsGeoSvc> m_actsGeoSvc;
  const dd4hep::rec::SurfaceMap* m_surfMap;
  const Acts::MagneticFieldProvider& m_field;

  using traj_t = Acts::VectorMultiTrajectory;
  using kfOptions_t = Acts::KalmanFitterOptions<traj_t>;

  std::shared_ptr<TrackFitterFunction> m_fitter;

};
