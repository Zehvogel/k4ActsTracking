#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

#include <cassert>

class MeasurementCalibrator {
 public:
  virtual void calibrate(const Acts::GeometryContext& gctx,
                 const Acts::CalibrationContext& cctx,
                 const Acts::SourceLink& sourceLink,
                 Acts::VectorMultiTrajectory::TrackStateProxy trackState) const = 0;

    virtual ~MeasurementCalibrator() = default;

};

// Calibrator to convert an index source link to a measurement as-is
class PassThroughCalibrator : public MeasurementCalibrator {
 public:
  void calibrate(const Acts::GeometryContext& gctx,
                 const Acts::CalibrationContext& cctx,
                 const Acts::SourceLink& sourceLink,
                 Acts::VectorMultiTrajectory::TrackStateProxy trackState) const override {
                 }

};