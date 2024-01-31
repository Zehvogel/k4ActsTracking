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
import sys
from Gaudi.Configuration import *

from Configurables import k4DataSvc, MarlinProcessorWrapper
algList = []
svcList = []


evtsvc = k4DataSvc("EventDataSvc")
svcList.append(evtsvc)

CONFIG = {
             "Tracking": "Truth",
            #  "Tracking": "Conformal",
             "TrackingChoices": ["Truth", "Conformal"],
             "InputMode": "LCIO",
             "InputModeChoices": ["LCIO", "EDM4hep"], # could also mix inputs but then things get ugly
             "OutputMode": "EDM4Hep",
             "OutputModeChoices": ["LCIO", "EDM4hep"] #, "both"] FIXME: both is not implemented yet
}

from Configurables import GeoSvc, TrackingCellIDEncodingSvc
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [os.environ["OPENDATADETECTOR"]+"/xml/OpenDataDetector.xml"]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False
svcList.append(geoservice)

# I could probably write some kind of CellIDTranslator for the ODD that unifies the mapping...

# FIXME: I forgot about this! Somewhat correct encoding is needed for truth-tracking
# cellIDSvc = TrackingCellIDEncodingSvc("CellIDSvc")
# cellIDSvc.EncodingStringParameterName = "GlobalTrackerReadoutID"
# cellIDSvc.GeoSvcName = geoservice.name()
# cellIDSvc.OutputLevel = INFO
# svcList.append(cellIDSvc)

output_basename = "output"

from k4FWCore.parseArgs import parser
parser.add_argument("--inputFiles", action="extend", nargs="+", metavar=("file1", "file2"), help="One or multiple input files")
parser.add_argument("--outputBasename", help="Basename of the output file(s)", default=output_basename)
my_opts = parser.parse_known_args()[0]

output_basename = my_opts.outputBasename

# Set input files here or via --inputFiles
input_files = []

if my_opts.inputFiles is not None:
    input_files = my_opts.inputFiles
print(f"opts: {my_opts}")
print(f"input_files: {input_files}")

if not input_files:
    print("Error: missing input files, set them via --inputFiles")
    sys.exit(1)

if input_files[0].endswith(".slcio"):
    CONFIG["InputMode"] = "LCIO"
elif input_files[0].endswith(".root"):
    CONFIG["InputMode"] = "EDM4hep"

if CONFIG["InputMode"] == "LCIO":
    from Configurables import LcioEvent
    read = LcioEvent()
    read.OutputLevel = WARNING
    read.Files = input_files
    algList.append(read)
elif CONFIG["InputMode"] == "EDM4hep":
    evtsvc.inputs = input_files
    from Configurables import PodioInput
    inp = PodioInput('InputReader')
    inp.collections = [
      'EventHeader',
      'MCParticles',
    #   Conveniently also named this way in ODD but no reco atm
    #   'ECalEndcapCollection',
    #   'ECalEndcapCollectionContributions',
    #   'ECalBarrelCollection',
    #   'ECalBarrelCollectionContributions',
    #   'HCalBarrelCollection',
    #   'HCalBarrelCollectionContributions',
    #   'HCalEndcapCollection',
    #   'HCalEndcapCollectionContributions',
    # Only use barrel for now...
      'LongStripBarrelReadout',
      'LongStripEndcapReadout',
      'PixelBarrelReadout',
      'PixelEndcapReadout',
      'ShortStripBarrelReadout',
      'ShortStripEndcapReadout',
    ]
    inp.OutputLevel = WARNING
    algList.append(inp)

MyAIDAProcessor = MarlinProcessorWrapper("MyAIDAProcessor")
MyAIDAProcessor.OutputLevel = WARNING
MyAIDAProcessor.ProcessorType = "AIDAProcessor"
MyAIDAProcessor.Parameters = {
                              "Compress": ["1"],
                              "FileName": [f"{output_basename}_aida"],
                              "FileType": ["root"]
                              }

if CONFIG["InputMode"] == "EDM4hep":
    from Configurables import EDM4hep2LcioTool
    EDM4hep2Lcio = EDM4hep2LcioTool("EDM4hep2Lcio")
    EDM4hep2Lcio.convertAll = False
    EDM4hep2Lcio.collNameMapping = {
        'EventHeader':                     'EventHeader',
        'MCParticles':                     'MCParticle',
        'LongStripBarrelReadout':          'LongStripBarrelReadout',
        'LongStripEndcapReadout':          'LongStripEndcapReadout',
        'PixelBarrelReadout':              'PixelBarrelReadout',
        'PixelEndcapReadout':              'PixelEndcapReadout',
        'ShortStripBarrelReadout':         'ShortStripBarrelReadout',
        'ShortStripEndcapReadout':         'ShortStripEndcapReadout',
        # 'ECalEndcapCollection':            'ECalEndcapCollection',
        # 'ECalBarrelCollection':            'ECalBarrelCollection',
        # 'HCalBarrelCollection':            'HCalBarrelCollection',
        # 'HCalEndcapCollection':            'HCalEndcapCollection',
        # 'HCalRingCollection':              'HCalRingCollection',
    }
    EDM4hep2Lcio.OutputLevel = DEBUG
    MyAIDAProcessor.EDM4hep2LcioTool = EDM4hep2Lcio


PixelBarrelDigitiser = MarlinProcessorWrapper("PixelBarrelDigitiser")
PixelBarrelDigitiser.OutputLevel = WARNING
PixelBarrelDigitiser.ProcessorType = "DDPlanarDigiProcessor"
PixelBarrelDigitiser.Parameters = {
                                 "IsStrip": ["false"],
                                 "ResolutionU": ["0.015"],
                                 "ResolutionV": ["0.015"],
                                 "SimTrackHitCollectionName": ["PixelBarrelReadout"],
                                 "SimTrkHitRelCollection": ["PixelBarrelTrackerHitRelations"],
                                 "SubDetectorName": ["Pixels"],
                                 "TrackerHitCollectionName": ["PixelBarrelTrackerHits"]
                                 }

ShortStripBarrelDigitiser = MarlinProcessorWrapper("ShortStripBarrelDigitiser")
ShortStripBarrelDigitiser.OutputLevel = WARNING
ShortStripBarrelDigitiser.ProcessorType = "DDPlanarDigiProcessor"
ShortStripBarrelDigitiser.Parameters = {
                                 "IsStrip": ["false"],
                                 "ResolutionU": ["0.043"],
                                 "ResolutionV": ["1.2"],
                                 "SimTrackHitCollectionName": ["ShortStripBarrelReadout"],
                                 "SimTrkHitRelCollection": ["ShortStripBarrelTrackerHitRelations"],
                                 "SubDetectorName": ["ShortStrips"],
                                 "TrackerHitCollectionName": ["ShortStripBarrelTrackerHits"]
                                 }

LongStripBarrelDigitiser = MarlinProcessorWrapper("LongStripBarrelDigitiser")
LongStripBarrelDigitiser.OutputLevel = WARNING
LongStripBarrelDigitiser.ProcessorType = "DDPlanarDigiProcessor"
LongStripBarrelDigitiser.Parameters = {
                                 "IsStrip": ["true"],
                                 "ResolutionU": ["0.072"],
                                 "ResolutionV": ["0.0"],
                                 "SimTrackHitCollectionName": ["LongStripBarrelReadout"],
                                 "SimTrkHitRelCollection": ["LongStripBarrelTrackerHitRelations"],
                                 "SubDetectorName": ["LongStrips"],
                                 "TrackerHitCollectionName": ["LongStripBarrelTrackerHits"]
                                 }

PixelEndcapDigitiser = MarlinProcessorWrapper("PixelEndcapDigitiser")
PixelEndcapDigitiser.OutputLevel = WARNING
PixelEndcapDigitiser.ProcessorType = "DDPlanarDigiProcessor"
PixelEndcapDigitiser.Parameters = {
                                 "IsStrip": ["false"],
                                 "ResolutionU": ["0.015"],
                                 "ResolutionV": ["0.015"],
                                 "SimTrackHitCollectionName": ["PixelEndcapReadout"],
                                 "SimTrkHitRelCollection": ["PixelEndcapTrackerHitRelations"],
                                 "SubDetectorName": ["Pixels"],
                                 "TrackerHitCollectionName": ["PixelEndcapTrackerHits"]
                                 }

ShortStripEndcapDigitiser = MarlinProcessorWrapper("ShortStripEndcapDigitiser")
ShortStripEndcapDigitiser.OutputLevel = WARNING
ShortStripEndcapDigitiser.ProcessorType = "DDPlanarDigiProcessor"
ShortStripEndcapDigitiser.Parameters = {
                                 "IsStrip": ["false"],
                                 "ResolutionU": ["0.043"],
                                 "ResolutionV": ["1.2"],
                                 "SimTrackHitCollectionName": ["ShortStripEndcapReadout"],
                                 "SimTrkHitRelCollection": ["ShortStripEndcapTrackerHitRelations"],
                                 "SubDetectorName": ["ShortStrips"],
                                 "TrackerHitCollectionName": ["ShortStripEndcapTrackerHits"]
                                 }

LongStripEndcapDigitiser = MarlinProcessorWrapper("LongStripEndcapDigitiser")
LongStripEndcapDigitiser.OutputLevel = WARNING
LongStripEndcapDigitiser.ProcessorType = "DDPlanarDigiProcessor"
LongStripEndcapDigitiser.Parameters = {
                                 "IsStrip": ["true"],
                                 "ResolutionU": ["0.072"],
                                 "ResolutionV": ["0.0"],
                                 "SimTrackHitCollectionName": ["LongStripEndcapReadout"],
                                 "SimTrkHitRelCollection": ["LongStripEndcapTrackerHitRelations"],
                                 "SubDetectorName": ["LongStrips"],
                                 "TrackerHitCollectionName": ["LongStripEndcapTrackerHits"]
                                 }

# LongStrips are real strips, need to build spacepoints before tracking
#   <processor name="FTDDDSpacePointBuilder" type="DDSpacePointBuilder">
#     <!--SpacePointBuilder combine si-strip measurements into 3D spacepoints (1TrackerHitPlanar+1TrackHitPlanar = 1 TrackerHit), that can be used by reconstruction-->
#     <!--Name of sub detector-->
#     <parameter name="SubDetectorName" type="string"> FTD </parameter>
#     <!--The length of the strips of the subdetector in mm-->
#     <parameter name="StripLength" type="Double"> 2.500000000e+02 </parameter>
#     <!--Name of the SpacePoint SimTrackerHit relation collection-->
#     <parameter name="SimHitSpacePointRelCollection" type="string" lcioOutType="LCRelation"> FTDSpacePointRelations </parameter>
#     <!--SpacePointsCollection-->
#     <parameter name="SpacePointsCollection" type="string" lcioOutType="TrackerHit"> FTDSpacePoints </parameter>
#     <!--TrackerHitCollection-->
#     <parameter name="TrackerHitCollection" type="string" lcioInType="TrackerHitPlane"> FTDStripTrackerHits </parameter>
#     <!--The name of the input collection of the relations of the TrackerHits to SimHits-->
#     <parameter name="TrackerHitSimHitRelCollection" type="string" lcioInType="LCRelation"> FTDStripTrackerHitRelations </parameter>
#     <!--verbosity level of this processor ("DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT")-->
#     <parameter name="Verbosity" type="string">MESSAGE </parameter>
#   </processor>
LongStripBarrelSpacePointBuilder = MarlinProcessorWrapper("LongStripBarrelSpacePointBuilder")
LongStripBarrelSpacePointBuilder.OutputLevel = WARNING
LongStripBarrelSpacePointBuilder.ProcessorType = "DDSpacePointBuilder"
LongStripBarrelSpacePointBuilder.Parameters = {
    "SubDetectorName": ["LongStrips"],
    "StripLength": ["78"],
    "SimHitSpacePointRelCollection": ["LongStripBarrelSpacePointRelations"],
    "SpacePointsCollection": ["LongStripBarrelSpacePoints"],
    "TrackerHitCollection": ["LongStripBarrelTrackerHits"],
    "TrackerHitSimHitRelCollection": ["LongStripBarrelTrackerHitRelations"],
}

LongStripEndcapSpacePointBuilder = MarlinProcessorWrapper("LongStripEndcapSpacePointBuilder")
LongStripEndcapSpacePointBuilder.OutputLevel = WARNING
LongStripEndcapSpacePointBuilder.ProcessorType = "DDSpacePointBuilder"
LongStripEndcapSpacePointBuilder.Parameters = {
    "SubDetectorName": ["LongStrips"],
    "StripLength": ["78"],
    "SimHitSpacePointRelCollection": ["LongStripEndcapSpacePointRelations"],
    "SpacePointsCollection": ["LongStripEndcapSpacePoints"],
    "TrackerHitCollection": ["LongStripEndcapTrackerHits"],
    "TrackerHitSimHitRelCollection": ["LongStripEndcapTrackerHitRelations"],
}

MyTruthTrackFinder = MarlinProcessorWrapper("MyTruthTrackFinder")
MyTruthTrackFinder.OutputLevel = WARNING
MyTruthTrackFinder.ProcessorType = "TruthTrackFinder"
MyTruthTrackFinder.Parameters = {
                                 "FitForward": ["true"],
                                 "MCParticleCollectionName": ["MCParticle"],
                                 "SiTrackCollectionName": ["SiTracks"],
                                 "SiTrackRelationCollectionName": ["SiTrackRelations"],
                                 "SimTrackerHitRelCollectionNames": [
                                     "PixelBarrelTrackerHitRelations",
                                    #  "ShortStripBarrelTrackerHitRelations",
                                    #  "LongStripBarrelSpacePointRelations",
                                     "PixelEndcapTrackerHitRelations",
                                    #  "ShortStripEndcapTrackerHitRelations",
                                    #  "LongStripEndcapSpacePointRelations",
                                     ],
                                 "TrackerHitCollectionNames": [
                                     "PixelBarrelTrackerHits",
                                     "ShortStripBarrelTrackerHits",
                                    #  "LongStripBarrelSpacePoints",
                                     "PixelEndcapTrackerHits",
                                     "ShortStripEndcapTrackerHits",
                                    #  "LongStripEndcapSpacePoints",
                                     ],
                                 "UseTruthInPrefit": ["false"]
                                 }

MyConformalTracking = MarlinProcessorWrapper("MyConformalTracking")
MyConformalTracking.OutputLevel = WARNING
MyConformalTracking.ProcessorType = "ConformalTrackingV2"
MyConformalTracking.Parameters = {
                                  "DebugHits": ["DebugHits"],
                                  "DebugPlots": ["false"],
                                  "DebugTiming": ["false"],
                                  "MCParticleCollectionName": ["MCParticle"],
                                  "MaxHitInvertedFit": ["0"],
                                  "MinClustersOnTrackAfterFit": ["3"],
                                  "RelationsNames": [
                                     "PixelBarrelTrackerHitRelations",
                                     "ShortStripBarrelTrackerHitRelations",
                                     "LongStripBarrelSpacePointRelations",
                                     "PixelEndcapTrackerHitRelations",
                                     "ShortStripEndcapTrackerHitRelations",
                                     "LongStripEndcapSpacePointRelations",
                                      ],
                                  "RetryTooManyTracks": ["false"],
                                  "SiTrackCollectionName": ["SiTracksCT"],
                                  "SortTreeResults": ["true"],
                                  "Steps":
                                  [
                                      "[VXDBarrel]",
                                      "@Collections", ":", "PixelBarrelTrackerHits",
                                      "@Parameters", ":", "MaxCellAngle", ":", "0.01;", "MaxCellAngleRZ", ":", "0.01;", "Chi2Cut", ":", "100;", "MinClustersOnTrack", ":", "4;", "MaxDistance", ":", "0.03;", "SlopeZRange:", "10.0;", "HighPTCut:", "10.0;",
                                      "@Flags", ":", "HighPTFit,", "VertexToTracker",
                                      "@Functions", ":", "CombineCollections,", "BuildNewTracks",
                                      "[VXDEncap]",
                                      "@Collections", ":", "PixelEndcapTrackerHits",
                                      "@Parameters", ":", "MaxCellAngle", ":", "0.01;", "MaxCellAngleRZ", ":", "0.01;", "Chi2Cut", ":", "100;", "MinClustersOnTrack", ":", "4;", "MaxDistance", ":", "0.03;", "SlopeZRange:", "10.0;", "HighPTCut:", "10.0;",
                                      "@Flags", ":", "HighPTFit,", "VertexToTracker",
                                      "@Functions", ":", "CombineCollections,", "ExtendTracks",
                                      "[LowerCellAngle1]",
                                      "@Collections", ":", "PixelBarrelTrackerHits,", "PixelEndcapTrackerHits",
                                      "@Parameters", ":", "MaxCellAngle", ":", "0.05;", "MaxCellAngleRZ", ":", "0.05;", "Chi2Cut", ":", "100;", "MinClustersOnTrack", ":", "4;", "MaxDistance", ":", "0.03;", "SlopeZRange:", "10.0;", "HighPTCut:", "10.0;",
                                      "@Flags", ":", "HighPTFit,", "VertexToTracker,", "RadialSearch",
                                      "@Functions", ":", "CombineCollections,", "BuildNewTracks",
                                      "[LowerCellAngle2]",
                                      "@Collections", ":",
                                      "@Parameters", ":", "MaxCellAngle", ":", "0.1;", "MaxCellAngleRZ", ":", "0.1;", "Chi2Cut", ":", "2000;", "MinClustersOnTrack", ":", "4;", "MaxDistance", ":", "0.03;", "SlopeZRange:", "10.0;", "HighPTCut:", "10.0;",
                                      "@Flags", ":", "HighPTFit,", "VertexToTracker,", "RadialSearch",
                                      "@Functions", ":", "BuildNewTracks,", "SortTracks",
                                      "[Tracker]",
                                      "@Collections", ":", "ShortStripBarrelTrackerHits,", "ShortStripEndcapTrackerHits,", "LongStripBarrelSpacePoints,", "LongStripEndcapSpacePoints,",
                                      "@Parameters", ":", "MaxCellAngle", ":", "0.1;", "MaxCellAngleRZ", ":", "0.1;", "Chi2Cut", ":", "2000;", "MinClustersOnTrack", ":", "4;", "MaxDistance", ":", "0.03;", "SlopeZRange:", "10.0;", "HighPTCut:", "1.0;",
                                      "@Flags", ":", "HighPTFit,", "VertexToTracker,", "RadialSearch",
                                      "@Functions", ":", "CombineCollections,", "ExtendTracks",
                                      "[Displaced]",
                                      "@Collections", ":", "PixelBarrelTrackerHits,", "PixelEndcapTrackerHits", "ShortStripBarrelTrackerHits,", "ShortStripEndcapTrackerHits,", "LongStripBarrelSpacePoints,", "LongStripEndcapSpacePoints,",
                                      "@Parameters", ":", "MaxCellAngle", ":", "0.1;", "MaxCellAngleRZ", ":", "0.1;", "Chi2Cut", ":", "1000;", "MinClustersOnTrack", ":", "5;", "MaxDistance", ":", "0.015;", "SlopeZRange:", "10.0;", "HighPTCut:", "10.0;",
                                      "@Flags", ":", "OnlyZSchi2cut,", "RadialSearch",
                                      "@Functions", ":", "CombineCollections,", "BuildNewTracks"
                                  ],
                                  "ThetaRange": ["0.05"],
                                  "TooManyTracks": ["100000"],
                                  "TrackerHitCollectionNames": [
                                     "PixelBarrelTrackerHits",
                                     "ShortStripBarrelTrackerHits",
                                     "LongStripBarrelSpacePoints",
                                     "PixelEndcapTrackerHits",
                                     "ShortStripEndcapTrackerHits",
                                     "LongStripEndcapSpacePoints",
                                      ],
                                  "trackPurity": ["0.7"]
                                  }

ClonesAndSplitTracksFinder = MarlinProcessorWrapper("ClonesAndSplitTracksFinder")
ClonesAndSplitTracksFinder.OutputLevel = WARNING
ClonesAndSplitTracksFinder.ProcessorType = "ClonesAndSplitTracksFinder"
ClonesAndSplitTracksFinder.Parameters = {
                                         "EnergyLossOn": ["true"],
                                         "InputTrackCollectionName": ["SiTracksCT"],
                                         "MultipleScatteringOn": ["true"],
                                         "OutputTrackCollectionName": ["SiTracks"],
                                         "SmoothOn": ["false"],
                                         "extrapolateForward": ["true"],
                                         "maxSignificancePhi": ["3"],
                                         "maxSignificancePt": ["2"],
                                         "maxSignificanceTheta": ["3"],
                                         "mergeSplitTracks": ["false"],
                                         "minTrackPt": ["1"]
                                         }

Refit = MarlinProcessorWrapper("Refit")
Refit.OutputLevel = WARNING
Refit.ProcessorType = "RefitFinal"
Refit.Parameters = {
                    "EnergyLossOn": ["true"],
                    "InputRelationCollectionName": ["SiTrackRelations"],
                    "InputTrackCollectionName": ["SiTracks"],
                    "Max_Chi2_Incr": ["1.79769e+30"],
                    "MinClustersOnTrackAfterFit": ["3"],
                    "MultipleScatteringOn": ["true"],
                    "OutputRelationCollectionName": ["SiTracks_Refitted_Relation"],
                    "OutputTrackCollectionName": ["SiTracks_Refitted"],
                    "ReferencePoint": ["-1"],
                    "SmoothOn": ["false"],
                    "extrapolateForward": ["true"]
                    }

MyClicEfficiencyCalculator = MarlinProcessorWrapper("MyClicEfficiencyCalculator")
MyClicEfficiencyCalculator.OutputLevel = WARNING
MyClicEfficiencyCalculator.ProcessorType = "ClicEfficiencyCalculator"
MyClicEfficiencyCalculator.Parameters = {
                                         "MCParticleCollectionName": ["MCParticle"],
                                         "MCParticleNotReco": ["MCParticleNotReco"],
                                         "MCPhysicsParticleCollectionName": ["MCPhysicsParticles"],
                                         "TrackCollectionName": ["SiTracks_Refitted"],
                                         "TrackerHitCollectionNames": [
                                             "PixelBarrelTrackerHits",
                                             "ShortStripBarrelTrackerHits",
                                             "LongStripBarrelSpacePoints",
                                             "PixelEndcapTrackerHits",
                                             "ShortStripEndcapTrackerHits",
                                             "LongStripEndcapSpacePoints",
                                             ],
                                         "TrackerHitRelCollectionNames": [
                                             "PixelBarrelTrackerHitRelations",
                                             "ShortStripBarrelTrackerHitRelations",
                                             "LongStripBarrelSpacePointRelations",
                                             "PixelEndcapTrackerHitRelations",
                                             "ShortStripEndcapTrackerHitRelations",
                                             "LongStripEndcapSpacePointRelations",
                                             ],
                                         "efficiencyTreeName": ["trktree"],
                                         "mcTreeName": ["mctree"],
                                         "morePlots": ["false"],
                                         "purityTreeName": ["puritytree"],
                                         "reconstructableDefinition": ["ILDLike"],
                                         "vertexBarrelID": ["1"]
                                         }

MyTrackChecker = MarlinProcessorWrapper("MyTrackChecker")
MyTrackChecker.OutputLevel = WARNING
MyTrackChecker.ProcessorType = "TrackChecker"
MyTrackChecker.Parameters = {
                             "MCParticleCollectionName": ["MCParticle"],
                             "TrackCollectionName": ["SiTracks_Refitted"],
                             "TrackRelationCollectionName": ["SiTracksMCTruthLink"],
                             "TreeName": ["checktree"],
                             "UseOnlyTree": ["true"]
                             }

MyStatusmonitor = MarlinProcessorWrapper("MyStatusmonitor")
MyStatusmonitor.OutputLevel = WARNING
MyStatusmonitor.ProcessorType = "Statusmonitor"
MyStatusmonitor.Parameters = {
                              "HowOften": ["100"]
                              }

# TODO: adapt
# MyRecoMCTruthLinker = MarlinProcessorWrapper("MyRecoMCTruthLinker")
# MyRecoMCTruthLinker.OutputLevel = WARNING
# MyRecoMCTruthLinker.ProcessorType = "RecoMCTruthLinker"
# MyRecoMCTruthLinker.Parameters = {
#                                   "BremsstrahlungEnergyCut": ["1"],
#                                   "CalohitMCTruthLinkName": ["CalohitMCTruthLink"],
#                                   "ClusterCollection": ["PandoraClusters"],
#                                   "ClusterMCTruthLinkName": ["ClusterMCTruthLink"],
#                                   "FullRecoRelation": ["true"],
#                                   "InvertedNonDestructiveInteractionLogic": ["false"],
#                                   "KeepDaughtersPDG": ["22", "111", "310", "13", "211", "321", "3120"],
#                                   "MCParticleCollection": ["MCPhysicsParticles"],
#                                   "MCParticlesSkimmedName": ["MCParticlesSkimmed"],
#                                   "MCTruthClusterLinkName": ["MCTruthClusterLink"],
#                                   "MCTruthRecoLinkName": ["MCTruthRecoLink"],
#                                   "MCTruthTrackLinkName": ["MCTruthSiTracksLink"],
#                                   "RecoMCTruthLinkName": ["RecoMCTruthLink"],
#                                   "RecoParticleCollection": ["PandoraPFOs"],
#                                   "SaveBremsstrahlungPhotons": ["true"],
#                                   "SimCaloHitCollections": ["ECalBarrelCollection", "ECalEndcapCollection", "HCalBarrelCollection", "HCalEndcapCollection", "HCalRingCollection", "YokeBarrelCollection", "YokeEndcapCollection", "LumiCalCollection"],
#                                   "SimCalorimeterHitRelationNames": ["RelationCaloHit", "RelationMuonHit"],
#                                   "SimTrackerHitCollections": ["VertexBarrelCollection", "VertexEndcapCollection", "InnerTrackerBarrelCollection", "OuterTrackerBarrelCollection", "InnerTrackerEndcapCollection", "OuterTrackerEndcapCollection"],
#                                   "TrackCollection": ["SiTracks_Refitted"],
#                                   "TrackMCTruthLinkName": ["SiTracksMCTruthLink"],
#                                   "TrackerHitsRelInputCollections": ["VXDTrackerHitRelations", "VXDEndcapTrackerHitRelations", "InnerTrackerBarrelHitsRelations", "OuterTrackerBarrelHitsRelations", "InnerTrackerEndcapHitsRelations", "OuterTrackerEndcapHitsRelations"],
#                                   "UseTrackerHitRelations": ["true"],
#                                   "UsingParticleGun": ["false"],
#                                   "daughtersECutMeV": ["10"]
#                                   }

MyHitResiduals = MarlinProcessorWrapper("MyHitResiduals")
MyHitResiduals.OutputLevel = WARNING
MyHitResiduals.ProcessorType = "HitResiduals"
MyHitResiduals.Parameters = {
                             "EnergyLossOn": ["true"],
                             "MaxChi2Increment": ["1000"],
                             "MultipleScatteringOn": ["true"],
                             "SmoothOn": ["false"],
                             "TrackCollectionName": ["SiTracks_Refitted"],
                             "outFileName": ["residuals.root"],
                             "treeName": ["restree"]
                             }


EventNumber = MarlinProcessorWrapper("EventNumber")
EventNumber.OutputLevel = WARNING
EventNumber.ProcessorType = "Statusmonitor"
EventNumber.Parameters = {
                          "HowOften": ["1"]
                          }

# TODO: put this somewhere else, needs to be in front of the output for now :(
# setup AIDA histogramming and add eventual background overlay
algList.append(MyAIDAProcessor)
# tracker hit digitisation
algList.append(PixelBarrelDigitiser)
algList.append(PixelEndcapDigitiser)
algList.append(ShortStripBarrelDigitiser)
algList.append(ShortStripEndcapDigitiser)
algList.append(LongStripBarrelDigitiser)
algList.append(LongStripEndcapDigitiser)
# spacepoint building for real strips
# algList.append(LongStripBarrelSpacePointBuilder)
# algList.append(LongStripEndcapSpacePointBuilder)

# tracking
# if CONFIG["Tracking"] == "Truth":
    # algList.append(MyTruthTrackFinder)
# elif CONFIG["Tracking"] == "Conformal":
#     algList.append(MyConformalTracking)
#     algList.append(ClonesAndSplitTracksFinder)

# algList.append(Refit)

# monitoring and Reco to MCTruth linking
# algList.append(MyClicEfficiencyCalculator)
# algList.append(MyRecoMCTruthLinker)
# algList.append(MyTrackChecker)
# event number processor, down here to attach the conversion back to edm4hep to it
algList.append(EventNumber)

if CONFIG["OutputMode"] == "LCIO":
    Output_REC = MarlinProcessorWrapper("Output_REC")
    Output_REC.OutputLevel = WARNING
    Output_REC.ProcessorType = "LCIOOutputProcessor"
    Output_REC.Parameters = {
                             "DropCollectionNames": [],
                             "DropCollectionTypes": [],
                             "FullSubsetCollections": ["EfficientMCParticles", "InefficientMCParticles"],
                             "KeepCollectionNames": [],
                             "LCIOOutputFile": [f"{output_basename}_REC.slcio"],
                             "LCIOWriteMode": ["WRITE_NEW"]
                             }

    Output_DST = MarlinProcessorWrapper("Output_DST")
    Output_DST.OutputLevel = WARNING
    Output_DST.ProcessorType = "LCIOOutputProcessor"
    Output_DST.Parameters = {
                             "DropCollectionNames": [],
                             "DropCollectionTypes": ["MCParticle", "LCRelation", "SimCalorimeterHit", "CalorimeterHit", "SimTrackerHit", "TrackerHit", "TrackerHitPlane", "Track", "ReconstructedParticle", "LCFloatVec"],
                             "FullSubsetCollections": ["EfficientMCParticles", "InefficientMCParticles", "MCPhysicsParticles"],
                             "KeepCollectionNames": ["MCParticlesSkimmed", "MCPhysicsParticles", "RecoMCTruthLink", "SiTracks", "SiTracks_Refitted", "PandoraClusters", "PandoraPFOs", "SelectedPandoraPFOs", "LooseSelectedPandoraPFOs", "TightSelectedPandoraPFOs", "RefinedVertexJets", "RefinedVertexJets_rel", "RefinedVertexJets_vtx", "RefinedVertexJets_vtx_RP", "BuildUpVertices", "BuildUpVertices_res", "BuildUpVertices_RP", "BuildUpVertices_res_RP", "BuildUpVertices_V0", "BuildUpVertices_V0_res", "BuildUpVertices_V0_RP", "BuildUpVertices_V0_res_RP", "PrimaryVertices", "PrimaryVertices_res", "PrimaryVertices_RP", "PrimaryVertices_res_RP", "RefinedVertices", "RefinedVertices_RP"],
                             "LCIOOutputFile": [f"{output_basename}_DST.slcio"],
                             "LCIOWriteMode": ["WRITE_NEW"]
                             }
    algList.append(Output_REC)
    algList.append(Output_DST)

if CONFIG["OutputMode"] == "EDM4Hep":
    from Configurables import Lcio2EDM4hepTool
    lcioConvTool = Lcio2EDM4hepTool("lcio2EDM4hep")
    lcioConvTool.convertAll = True
    lcioConvTool.collNameMapping = {
        "MCParticle": "MCParticles"
    }
    lcioConvTool.OutputLevel = DEBUG
# attach to the last non output processor
    EventNumber.Lcio2EDM4hepTool = lcioConvTool

    from Configurables import PodioOutput
    out = PodioOutput("PodioOutput", filename = f"{output_basename}_edm4hep.root")
    out.outputCommands = ["keep *"]
    algList.append(out)

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = algList,
                EvtSel = 'NONE',
                EvtMax = 3,
                ExtSvc = svcList,
                OutputLevel=WARNING
              )
