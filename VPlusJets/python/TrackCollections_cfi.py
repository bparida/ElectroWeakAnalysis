import FWCore.ParameterSet.Config as cms

from ElectroWeakAnalysis.VPlusJets.AllPassFilter_cfi import AllPassFilter

#--------------------------
# Counter1: All read events
AllEventsStep = AllPassFilter.clone()

##-------- Scraping veto --------
noscraping = cms.EDFilter("FilterOutScraping",
   applyfilter = cms.untracked.bool(True),
   debugOn = cms.untracked.bool(False),
   numtrack = cms.untracked.uint32(10),
   thresh = cms.untracked.double(0.25)
)
noscrapingStep = AllPassFilter.clone()

##---------HBHE Noise Filter ------
from CommonTools.RecoAlgos.HBHENoiseFilter_cfi import HBHENoiseFilter
HBHENoiseFilter.minIsolatedNoiseSumE = cms.double(999999.)
HBHENoiseFilter.minNumIsolatedNoiseChannels = cms.int32(999999)
HBHENoiseFilter.minIsolatedNoiseSumEt = cms.double(999999.)
HBHENoiseStep = AllPassFilter.clone()



##-------- Primary vertex filter --------
primaryVertex = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),                                
    cut = cms.string("!isFake && ndof >= 4 && abs(z) <= 24 && position.Rho <= 2"), # tracksSize() > 3 for the older cut
    filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)
primaryVertexStep = AllPassFilter.clone()




TrackVtxPath = cms.Sequence(
    AllEventsStep +
    noscraping +
    noscrapingStep +
    HBHENoiseFilter + 
    HBHENoiseStep +
    primaryVertex +
    primaryVertexStep
)
