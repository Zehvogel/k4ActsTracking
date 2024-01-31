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

#ifndef IACTSGEOSVC_H
#define IACTSGEOSVC_H

#include <GaudiKernel/IService.h>
#include <unordered_map>

namespace dd4hep {
  namespace rec {
    class Surface;
  }
  class Detector;
}  // namespace dd4hep

namespace Acts {
  class TrackingGeometry;
  class Surface;
  class ContextType;
  // class GeometryContext;
  // XXX: there _must_ be a better way to do this
  using GeometryContext = ContextType;
}  // namespace Acts

class GAUDI_API IActsGeoSvc : virtual public IService {
public:
  using VolumeSurfaceMap = std::map<uint64_t, const Acts::Surface*>;

public:
  DeclareInterfaceID(IActsGeoSvc, 1, 0);

  virtual const Acts::TrackingGeometry& trackingGeometry() const = 0;
  virtual const VolumeSurfaceMap& surfaceMap() const = 0;

  virtual const Acts::GeometryContext& geometryContext() const = 0;

  virtual ~IActsGeoSvc() {}
};

#endif  // IACTSGEOSVC_H
