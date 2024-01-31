#
# Copyright (c) 2014-2023 Key4hep-Project.
#
# This file is part of Key4hep.
# See https://key4hep.github.io/key4hep-doc/ for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import os

from Gaudi.Configuration import INFO, DEBUG
from Configurables import ACTSRefit
from Configurables import ActsGeoSvc, GeoSvc
from Configurables import ApplicationMgr
from Configurables import k4DataSvc
from Configurables import PodioOutput
from Configurables import PodioInput

event_data_svc = k4DataSvc("EventDataSvc")
# TODO change to something publicly accessible
# event_data_svc.input = "track_test5_edm4hep.root"
event_data_svc.input = "single_mu-_1GeV_80deg_fix-offset_REC_edm4hep.root"
# event_data_svc.input = "/eos/experiment/clicdp/grid/ilc/user/L/LReichenbac/resolutions/rec_e4h/CLD_o2_v05/REC_e-_80deg_10GeV_1000evt_edm4hep.root"
event_data_svc.FirstEventEntry = 0

dd4hep_geo = GeoSvc("GeoSvc")
dd4hep_geo.detectors = [f"{os.environ['OPENDATADETECTOR']}/xml/OpenDataDetector.xml"]
dd4hep_geo.EnableGeant4Geo = False

acts_geo = ActsGeoSvc("ActsGeoSvc")
acts_geo.GeoSvcName = dd4hep_geo.name()
acts_geo.debugGeometry = True
acts_geo.outputFileName = "MyObjFile"

out = PodioOutput("out")
out.filename = "test.root"
# The collections that we don't drop will also be present in the output file
out.outputCommands = ["keep EmptyGarbage"]

inp = PodioInput("InputReader")
inp.collections = [
    "PixelBarrelTrackerHits"
]

transformer = ACTSRefit("ACTSRefit",
                        InputCollection="PixelBarrelTrackerHits",
                        OutputCollection="EmptyGarbage")
transformer.OutputLevel = DEBUG

ApplicationMgr(TopAlg=[inp, transformer, out],
               EvtSel="NONE",
               EvtMax=5,
               ExtSvc=[event_data_svc, dd4hep_geo, acts_geo],
               OutputLevel=INFO,
               )
