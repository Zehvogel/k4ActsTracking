import os
import sys
from pprint import pprint
from Gaudi.Configuration import *

from Configurables import ActsGeoSvc

algList = []

a = ActsGeoSvc("ActsGeoSvc")
a.detectors = ["/home/delitez/ACTS/acts/thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
a.debugGeometry = True
a.outputFileName = "MyObjFile"
from Configurables import ApplicationMgr

ApplicationMgr(TopAlg=algList, EvtSel="NONE", EvtMax=2, ExtSvc=[
a
], OutputLevel=INFO)
