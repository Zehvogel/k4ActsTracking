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

#include "TruthTracking.h"
#include "IndexSourceLink.h"
#include "MeasurementCalibration.h"
#include "edm4hep/utils/vector_utils.h"

// needed for lookup, might be able to remove it later
#include "DD4hep/Detector.h"

#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"

DECLARE_COMPONENT(TruthTracking)

TruthTracking::TruthTracking(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(
          name, svcLoc,
          {KeyValue("EventHeaderCollection", "EventHeader"), KeyValue("MCParticleCollection", "MCParticles"), KeyValue("TrackerHitCollection", "PixelBarrelTrackerHits"),
           KeyValue("TrackerHitAssociations", "PixelBarrelTrackerHitRelations")},
          {KeyValue("TrackCollection", "TruthTracks")}) {}

StatusCode TruthTracking::initialize() {
  m_geoSvc     = service("GeoSvc");
  m_actsGeoSvc = service("ActsGeoSvc");

  auto geo = std::make_shared<Acts::TrackingGeometry>(m_actsGeoSvc->trackingGeometry());
  // FIXME: use dd4hep supplied field
  // https://github.com/acts-project/acts/pull/2593
  // and distribute it via the ActsGeoSvc
  auto field = std::make_shared<Acts::ConstantBField>(Acts::Vector3(0., 0., 2.0 * Acts::UnitConstants::T));

  // verbose dies somewhere because of a surface nullpointer but is still useful if the track fit crashes earlier e.g. during smoothing
  // m_fitter = makeKalmanFitterFunction(geo, field, true, true, 0.0, Acts::FreeToBoundCorrection(), *Acts::getDefaultLogger("Kalman", Acts::Logging::Level::VERBOSE));
  m_fitter = makeKalmanFitterFunction(geo, field, true, true, 0.0, Acts::FreeToBoundCorrection(), *Acts::getDefaultLogger("Kalman", Acts::Logging::Level::DEBUG));
  // m_fitter = makeKalmanFitterFunction(geo, field, true, true, 0.0, Acts::FreeToBoundCorrection(), *Acts::getDefaultLogger("Kalman", Acts::Logging::Level::INFO));

  // FIXME: use dd4hep supplied surface?
  auto ipSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3{0., 0., 0.});

  m_geoCtx  = m_actsGeoSvc->geometryContext();
  m_options = std::make_shared<TrackFitterFunction::GeneralFitterOptions>(
      m_geoCtx, m_fieldCtx, m_calibCtx, ipSurface.get(), Acts::PropagatorPlainOptions());
  // TODO:: also do calibrator setup here

  return StatusCode::SUCCESS;
}

std::tuple<TrackCollection> TruthTracking::operator()(
    const EventHeaderCollection& evtHeader, const MCParticleCollection& MCParticles,
    const TrackerHitPlaneCollection& hits, const MCRecoTrackerHitPlaneAssociationCollection& associations) const {
  using edm4hep::MCParticle;
  using edm4hep::TrackerHitPlane;

  info() << "Hello from event: " << evtHeader[0].getEventNumber() << endmsg;

  std::map<const int, std::vector<TrackerHitPlane>> hitMap;
  // There could be more than one SimTrackerHit per TrackerHit caused by different particles
  // but ignoring this should not be so bad as this way we will just use the hit in the track
  // reconstruction of both particles, which is technically correct?
  // However, this means that we might produce more tracks than the TruthTrackFinder!
  // TODO: add another loop here to loop over the vector of association collections once we have such a thing...
  for (const auto& assoc : associations) {
    const auto& simHit = assoc.getSim();
    const auto& recHit = assoc.getRec();
    // If the hit was produced by a secondary which was not saved to the MCParticle collection
    if (simHit.isProducedBySecondary()) {
      continue;
    }
    const auto& mcp = simHit.getMCParticle();
    // the empty vector should be default constructible so that this works?
    hitMap[mcp.id().index].push_back(recHit);
  }

  auto                                trackContainer      = std::make_shared<Acts::VectorTrackContainer>();
  auto                                trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackFitterFunction::TrackContainer tracks(trackContainer, trackStateContainer);

  for (auto&& [mcp_idx, recHits] : hitMap) {
    // sort hits?
    // FIXME: what about spiraling tracks?? TruthTrackFinder only sorts hits if first fit attempt fails
    // std::sort(recHits.begin(), recHits.end(), [](const TrackerHitPlane& hit1, const TrackerHitPlane& hit2) {
    //   const auto& pos1 = hit1.getPosition();
    //   const auto& pos2 = hit2.getPosition();
    //   using edm4hep::utils::magnitudeTransverse;
    //   return magnitudeTransverse(pos1) < magnitudeTransverse(pos2);
    // });

    // do all the ugly conversion and setup stuff
    // iterate over recHits and create measurements and sourceLinks
    MeasurementContainer measurements;
    measurements.reserve(recHits.size());
    std::vector<IndexSourceLink> sourceLinks;
    sourceLinks.reserve(recHits.size());

    for (const auto& hit : recHits) {
      const auto& pos = hit.getPosition();
      // auto        z   = pos.z;
      // auto        r   = sqrt(pos.x * pos.x + pos.y * pos.y);
      // info() << "hello world from hit with z, r: " << z << ", " << r << endmsg;

      // TODO: ask Andre if this is the way to do it or if something easier is clean enough
      // in principle this gives you always the same as just taking the first 32 bit of the cellID
      auto       detElement = m_geoSvc->getDetector()->volumeManager().lookupDetElement(hit.getCellID());
      auto       volID      = detElement.volumeID();
      const auto is         = m_actsGeoSvc->surfaceMap().find(volID);
      if (is == m_actsGeoSvc->surfaceMap().end()) {
        error() << " volID (" << volID << ")  not found in acts surfaces" << endmsg;
        continue;
      }
      const Acts::Surface* surf = is->second;

      // point to the measurement that we are adding now which index will be current last measurement index plus one i.e. the size
      IndexSourceLink sourceLink{surf->geometryId(), measurements.size()};
      sourceLinks.push_back(std::move(sourceLink));

      // FIXME: figure out what to do with the momentum direction (the three 0's)
      // FIXME: check the Acts::Result before accessing value()
      auto localPos = surf->globalToLocal(m_geoCtx, {pos.x, pos.y, pos.z}, {0, 0, 0}).value();

      Acts::Vector2 loc     = Acts::Vector2::Zero();
      loc[Acts::eBoundLoc0] = localPos.x();
      loc[Acts::eBoundLoc1] = localPos.y();

      Acts::SquareMatrix2 cov = Acts::SquareMatrix2::Zero();
      cov(0, 0)               = std::pow(hit.getDu(), 2) * std::pow(Acts::UnitConstants::mm, 2);
      cov(1, 1)               = std::pow(hit.getDv(), 2) * std::pow(Acts::UnitConstants::mm, 2);

      auto measurement =
          Acts::makeMeasurement(Acts::SourceLink{sourceLink}, loc, cov, Acts::eBoundLoc0, Acts::eBoundLoc1);
      measurements.push_back(measurement);
    }

    const auto& mcp = MCParticles[mcp_idx];
    const auto& p   = mcp.getMomentum();
    using Acts::UnitConstants::GeV;
    using Acts::UnitConstants::um;

    Acts::BoundVector params;
    params(Acts::eBoundLoc0)   = 0.0;  // 0 because we put the surface for this state directly
    params(Acts::eBoundLoc1)   = 0.0;  // at the position of the particle
    params(Acts::eBoundPhi)    = edm4hep::utils::angleAzimuthal(p);
    params(Acts::eBoundTheta)  = edm4hep::utils::anglePolar(p);
    params(Acts::eBoundQOverP) = mcp.getCharge() / (edm4hep::utils::magnitude(p) * GeV);

    Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Zero();
    // FIXME: values copy pasted from EIC, check what would be physical for a HF detector
    cov(Acts::eBoundLoc0, Acts::eBoundLoc0)     = 1000 * um * 1000 * um;
    cov(Acts::eBoundLoc1, Acts::eBoundLoc1)     = 1000 * um * 1000 * um;
    cov(Acts::eBoundPhi, Acts::eBoundPhi)       = 0.05 * 0.05;
    cov(Acts::eBoundTheta, Acts::eBoundTheta)   = 0.01 * 0.01;
    cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = (0.1 * 0.1) / (GeV * GeV);
    cov(Acts::eBoundTime, Acts::eBoundTime)     = 1000 * Acts::UnitConstants::us;

    const auto& pos      = mcp.getVertex();
    auto        pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3{pos.x, pos.y, pos.z});

    TrackFitterFunction::TrackParameters initialParams{pSurface, params, cov, Acts::ParticleHypothesis::pion()};

    PassThroughCalibrator        calib;
    MeasurementCalibratorAdapter calibrator(calib, measurements);

    // FIXME: I think I tried to get rid of them before but ran into trouble
    std::vector<Acts::SourceLink> trackSourceLinks;
    trackSourceLinks.reserve(sourceLinks.size());

    for (const auto& sl : sourceLinks) {
      trackSourceLinks.push_back(Acts::SourceLink{sl});
    }

    auto result = (*m_fitter)(trackSourceLinks, initialParams, *m_options, calibrator, tracks);
    if (result.ok()) {
      info() << "track fit ok :)" << endmsg;
    } else {
      warning() << "track fit failed :(" << endmsg;
      continue;
    }

    // TODO: take fitted track parameters and convert them to an edm4hep::Track!!
    auto track = result.value();
    info() << "track momentum: " << track.absoluteMomentum() << endmsg;
    // fitter.fit(hits, etc. pp.)
    // TODO: if fit fails sort hits by r and fit again
    // TODO: what about sorting by abs(z)? I did not see this but wouldn't it make sense??
    // TODO: as I have the SimHits I could also sort by perfectly resolved time, which should be the most correct...
    // TODO: don't forget about atCalo track state
    // convert result back to edm4hep::Track(States)
    // create association between hits used for the fit and the track
  }

  // just a dummy so it compiles
  auto coll_out = TrackCollection();
  return std::make_tuple(std::move(coll_out));
}