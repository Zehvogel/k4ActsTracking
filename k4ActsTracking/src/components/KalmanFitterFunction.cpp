#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
// for v32
#include "Acts/Propagator/Navigator.hpp"

#include "TrackFitterFunction.h"
// #include "Calibrator.h"
#include "MeasurementCalibration.h"
#include "IndexSourceLink.h"

namespace Acts {
class MagneticFieldProvider;
class SourceLink;
class Surface;
class TrackingGeometry;
}  // namespace Acts

namespace {

using Stepper = Acts::EigenStepper<>;
using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
using Fitter = Acts::KalmanFitter<Propagator, Acts::VectorMultiTrajectory>;
using DirectPropagator = Acts::Propagator<Stepper, Acts::DirectNavigator>;
using DirectFitter =
    Acts::KalmanFitter<DirectPropagator, Acts::VectorMultiTrajectory>;

using TrackContainer =
    Acts::TrackContainer<Acts::VectorTrackContainer,
                         Acts::VectorMultiTrajectory, std::shared_ptr>;

struct SimpleReverseFilteringLogic {
  double momentumThreshold = 0;

  bool doBackwardFiltering(
      Acts::VectorMultiTrajectory::ConstTrackStateProxy trackState) const {
    auto momentum = fabs(1 / trackState.filtered()[Acts::eBoundQOverP]);
    return (momentum <= momentumThreshold);
  }
};

struct KalmanFitterFunctionImpl final : public TrackFitterFunction {
  Fitter fitter;
  DirectFitter directFitter;

  Acts::GainMatrixUpdater kfUpdater;
  Acts::GainMatrixSmoother kfSmoother;
  SimpleReverseFilteringLogic reverseFilteringLogic;

  bool multipleScattering = false;
  bool energyLoss = false;
  Acts::FreeToBoundCorrection freeToBoundCorrection;

  IndexSourceLink::SurfaceAccessor slSurfaceAccessor;

  KalmanFitterFunctionImpl(Fitter&& f, DirectFitter&& df,
                           const Acts::TrackingGeometry& trkGeo)
      : fitter(std::move(f)),
        directFitter(std::move(df)),
        slSurfaceAccessor{trkGeo} {}

  template <typename calibrator_t>
  auto makeKfOptions(const GeneralFitterOptions& options,
                     const calibrator_t& calibrator) const {
    Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions;
    extensions.updater.connect<
        &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
        &kfUpdater);
    extensions.smoother.connect<
        &Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(
        &kfSmoother);
    extensions.reverseFilteringLogic
        .connect<&SimpleReverseFilteringLogic::doBackwardFiltering>(
            &reverseFilteringLogic);

    Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory> kfOptions(
        options.geoContext, options.magFieldContext, options.calibrationContext,
        extensions, options.propOptions, &(*options.referenceSurface));

    kfOptions.referenceSurfaceStrategy =
        Acts::KalmanFitterTargetSurfaceStrategy::first;
    kfOptions.multipleScattering = multipleScattering;
    kfOptions.energyLoss = energyLoss;
    kfOptions.freeToBoundCorrection = freeToBoundCorrection;
    kfOptions.extensions.calibrator.connect<&calibrator_t::calibrate>(
        &calibrator);
    // not good enough :(
    // auto calibrator = [](Acts::GeometryContext& gtx, Acts::CalibrationContext& ctx, Acts::SourceLink& sl, Acts::VectorMultiTrajectory::TrackStateProxy ts){};
    // kfOptions.extensions.calibrator.connect(calibrator);
    kfOptions.extensions.surfaceAccessor
        .connect<&IndexSourceLink::SurfaceAccessor::operator()>(
            &slSurfaceAccessor);

    return kfOptions;
  }

  TrackFitterResult operator()(const std::vector<Acts::SourceLink>& sourceLinks,
                               const TrackParameters& initialParameters,
                               const GeneralFitterOptions& options,
                               const MeasurementCalibratorAdapter& calibrator,
                               TrackContainer& tracks) const override {
    const auto kfOptions = makeKfOptions(options, calibrator);
    // const auto kfOptions = makeKfOptions(options);
    return fitter.fit(sourceLinks.begin(), sourceLinks.end(), initialParameters,
                      kfOptions, tracks);
  }

//   TrackFitterResult operator()(
//       const std::vector<Acts::SourceLink>& sourceLinks,
//       const TrackParameters& initialParameters,
//       const GeneralFitterOptions& options,
//       const RefittingCalibrator& calibrator,
//       const std::vector<const Acts::Surface*>& surfaceSequence,
//       TrackContainer& tracks) const override {
//     const auto kfOptions = makeKfOptions(options, calibrator);
//     return directFitter.fit(sourceLinks.begin(), sourceLinks.end(),
//                             initialParameters, kfOptions, surfaceSequence,
//                             tracks);
//   }
};

}  // namespace

std::shared_ptr<TrackFitterFunction>
makeKalmanFitterFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    bool multipleScattering, bool energyLoss,
    double reverseFilteringMomThreshold,
    Acts::FreeToBoundCorrection freeToBoundCorrection,
    const Acts::Logger& logger) {
  // Stepper should be copied into the fitters
  const Stepper stepper(std::move(magneticField));

  // Standard fitter
  const auto& geo = *trackingGeometry;
  Acts::Navigator::Config cfg{std::move(trackingGeometry)};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg, logger.cloneWithSuffix("Navigator"));
  Propagator propagator(stepper, std::move(navigator),
                        logger.cloneWithSuffix("Propagator"));
  Fitter trackFitter(std::move(propagator), logger.cloneWithSuffix("Fitter"));

  // Direct fitter
  Acts::DirectNavigator directNavigator{
      logger.cloneWithSuffix("DirectNavigator")};
  DirectPropagator directPropagator(stepper, std::move(directNavigator),
                                    logger.cloneWithSuffix("DirectPropagator"));
  DirectFitter directTrackFitter(std::move(directPropagator),
                                 logger.cloneWithSuffix("DirectFitter"));

  // build the fitter function. owns the fitter object.
  auto fitterFunction = std::make_shared<KalmanFitterFunctionImpl>(
      std::move(trackFitter), std::move(directTrackFitter), geo);
  fitterFunction->multipleScattering = multipleScattering;
  fitterFunction->energyLoss = energyLoss;
  fitterFunction->reverseFilteringLogic.momentumThreshold =
      reverseFilteringMomThreshold;
  fitterFunction->freeToBoundCorrection = freeToBoundCorrection;

  return fitterFunction;
}