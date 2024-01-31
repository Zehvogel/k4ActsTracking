#! /usr/bin/bash
set -e

COMPACT_FILE=../OpenDataDetector/xml/OpenDataDetector.xml
# COMPACT_FILE=$OPENDATADETECTOR/xml/OpenDataDetector.xml

ddsim --compactFile $COMPACT_FILE \
      --outputFile single_mu-_1GeV_80deg_fix-offset.edm4hep.root \
      --random.seed 0123456789 \
      --enableGun \
      --gun.particle mu- \
      --gun.energy "1*GeV" \
      --gun.distribution uniform \
      --gun.thetaMin 80 \
      --gun.thetaMax 80 \
      --numberOfEvents 1000
