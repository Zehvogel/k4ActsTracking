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

#include "ACTSRefit.h"
#include "GaudiKernel/StatusCode.h"
#include "DD4hep/Detector.h"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"

#include "TrackFitterFunction.h"
#include "IndexSourceLink.h"

DECLARE_COMPONENT(ACTSRefit)

ACTSRefit::ACTSRefit(const std::string& name, ISvcLocator* svcLoc)
    : Transformer(name, svcLoc, KeyValue("InputCollection", "MCParticles"),
                  KeyValue("OutputCollection", "NewMCParticles")), m_field{Acts::ConstantBField({0., 0., 2.0 * Acts::UnitConstants::T})}
                   {
}

StatusCode ACTSRefit::initialize() {
  m_geoSvc = service("GeoSvc");
  m_actsGeoSvc = service("ActsGeoSvc");

  m_surfMap = m_geoSvc->getDetector()->extension<dd4hep::rec::SurfaceManager>()->map("tracker");

  auto geo = std::make_shared<Acts::TrackingGeometry>(m_actsGeoSvc->trackingGeometry());
  auto field = std::make_shared<Acts::ConstantBField>(Acts::Vector3(0., 0., 2.0 * Acts::UnitConstants::T));

  m_fitter = makeKalmanFitterFunction(geo, field, true, true, 0.0, Acts::FreeToBoundCorrection(), *Acts::getDefaultLogger("Kalman", Acts::Logging::Level::VERBOSE));

  return StatusCode::SUCCESS;
}

// This is the function that will be called to transform the data
// Note that the function has to be const, as well as all pointers to collections
// we get from the input
TrackCollection ACTSRefit::operator()(const TrackerHitPlaneCollection& input) const {
  // first convert tracker hit planes to acts measurements
  // TODO: refactor
  MeasurementContainer measurements;
  measurements.reserve(input.size());
  std::vector<IndexSourceLink> sourceLinks;
  sourceLinks.reserve(input.size());
  for (const auto& hit : input) {
    const auto& pos = hit.getPosition();
    auto z = pos.z;
    auto r = sqrt(pos.x * pos.x + pos.y * pos.y);
    info() << "hello world from hit with z, r: " << z << ", " << r << endmsg;
    // basically the same as taking the first 32 bits but I guess the cleaner way to do it
    auto detElement = m_geoSvc->getDetector()->volumeManager().lookupDetElement(hit.getCellID());
    auto volID = detElement.volumeID();
    const auto is = m_actsGeoSvc->surfaceMap().find(volID);
    if (is == m_actsGeoSvc->surfaceMap().end()) {
      error() << " volID (" << volID << ")  not found in acts surfaces" << endmsg;
      continue;
    }
    auto geo_ctx = m_actsGeoSvc->geometryContext();
    const Acts::Surface* surf = is->second;
    info() << "got surface " << surf->toString(m_actsGeoSvc->geometryContext()) << endmsg;
    Acts::Vector3 actsPos(pos.x, pos.y, pos.z);
    auto localPos = surf->globalToLocal(m_actsGeoSvc->geometryContext(), {pos.x, pos.y, pos.z}, {0, 0, 0}).value();
    info() << "localPos: (" << localPos.x() << ", " << localPos.y() << ")" << endmsg;
    auto center = surf->center(geo_ctx);
    info() << "center: (" << center.x() << ", " << center.y() << ", " << center.z() << ")" << endmsg;
    auto localCenter = surf->globalToLocal(geo_ctx, center, {0, 0, 0}).value();
    info() << "localCenter: (" << localCenter.x() << ", " << localCenter.y() << ")" << endmsg;

    const auto iter = m_surfMap->find(volID);
    if (iter == m_surfMap->end()) {
      error() << " volID (" << volID << ")  not found in dd4hep surfaces" << endmsg;
      continue;
    }
    const auto* dd_surf = iter->second;

    dd4hep::rec::Vector3D dd_pos(pos.x, pos.y, pos.z);
    dd_pos = dd4hep::mm * dd_pos;
    auto dd_localPos = dd_surf->globalToLocal(dd_pos);
    auto dd_origin = dd_surf->globalToLocal(dd_surf->origin());
    auto origin = dd_surf->origin();
    info() << "origin: (" << origin.x() << ", " << origin.y() << ", " << origin.z() << ")" << endmsg;
    info() << "localOrigin: (" << dd_origin.u() << ", " << dd_origin.v() << ")" << endmsg;
    info() << "dd_localPos: (" << dd_localPos.u() << ", " << dd_localPos.v() << ")" << endmsg;
    info() << "lengths: (" << dd_surf->length_along_u() << ", " << dd_surf->length_along_v() << ")" << endmsg;
    auto inbounds = dd_surf->insideBounds(dd_pos) ? "yes" : "no";
    info() << "dd_localPos inbounds? " << inbounds << endmsg;

    // same results as for acts but factor 10 from unit conversion which I need to be careful with
    // especially nice, ddrec u,v lines up with acts 2d x,y therefore cov is just diag(squared errors)
    // hopefully this will be the case for all detectors :^), will need a change for stereo angle strips where <u, v> != 0
    // TODO: create cov matrix -> need to match first and second dim of dd4hep surface with those of acts surface...
      Acts::SquareMatrix2 cov = Acts::SquareMatrix2::Zero();
      cov(0, 0)            = std::pow(hit.getDu() * Acts::UnitConstants::mm, 2);
      cov(1, 1)            = std::pow(hit.getDv() * Acts::UnitConstants::mm, 2);

      Acts::Vector2 loc     = Acts::Vector2::Zero();
      loc[Acts::eBoundLoc0] = localPos.x();
      loc[Acts::eBoundLoc1] = localPos.y();

      // point to the measurement that we are adding now which index will be current last measurement index plus one i.e. the size
      IndexSourceLink sourceLink{surf->geometryId(), measurements.size()};
      // Acts::SourceLink sourceLink{surf->geometryId()};
      // just do the simplest thing and use the geometry id as the source link as we will not care about calibration for the next 15-20 years
      auto measurement = Acts::makeMeasurement(Acts::SourceLink{sourceLink}, loc, cov, Acts::eBoundLoc0, Acts::eBoundLoc1);
      measurements.push_back(measurement);
      sourceLinks.push_back(sourceLink);
  }

  // TODO: use dd4hep supplied IP surface
  auto ipSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3{0., 0., 0.});

  // TODO: calibrator could be instantiated here

  auto geoContext = m_actsGeoSvc->geometryContext();
  Acts::MagneticFieldContext fieldContext;
  Acts::CalibrationContext calibContext;

  TrackFitterFunction::GeneralFitterOptions options{
    geoContext, fieldContext, calibContext, ipSurface.get(), Acts::PropagatorPlainOptions()
  };

  // TODO: only one track at the moment
  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackFitterFunction::TrackContainer tracks(trackContainer, trackStateContainer);
  // TODO: this is where the loop over tracks would start

  Acts::BoundVector params;
  params(Acts::eBoundLoc0) = 0.0;
  params(Acts::eBoundLoc1) = 0.0;
  auto firstHit = input.at(0);
  auto pos = firstHit.getPosition();
  Acts::Vector3 a_pos(pos.x, pos.y, pos.z);
  params(Acts::eBoundPhi) = Acts::VectorHelpers::phi(a_pos);
  params(Acts::eBoundTheta) = Acts::VectorHelpers::theta(a_pos);
  params(Acts::eBoundQOverP) = -1.0 / 1 * Acts::UnitConstants::GeV;

  TrackFitterFunction::TrackParameters initialParams{ipSurface, params, {}, Acts::ParticleHypothesis::pion()};
  std::vector<Acts::SourceLink> trackSourceLinks;
  trackSourceLinks.reserve(sourceLinks.size());

  for (const auto& sl : sourceLinks) {
    trackSourceLinks.push_back(Acts::SourceLink{sl});
  }

  PassThroughCalibrator calib;
  MeasurementCalibratorAdapter calibrator(calib, measurements);

  // FIXME: can not use voidFitterCalibrator urgh
  auto result = (*m_fitter)(trackSourceLinks, initialParams, options, calibrator, tracks);
  if (result.ok()) {
    info() << "track fit ok :o" << endmsg;
  } else {
    error() << "track fit failed :(" << endmsg;
  }

  auto coll_out = TrackCollection();
  return coll_out;
}
