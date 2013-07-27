import FWCore.ParameterSet.Config as cms

from ElectroWeakAnalysis.VPlusJets.WenuCollections_cfi import looseElectrons
from ElectroWeakAnalysis.VPlusJets.WenuCollections_cfi import looseMuons

isQCD = False

isolationCutString = cms.string("")
if isQCD:
    isolationCutString = "(pfIsolationR04().sumChargedHadronPt+max(0.,pfIsolationR04().sumNeutralHadronEt+pfIsolationR04().sumPhotonEt-0.5*pfIsolationR04().sumPUPt))/pt> 0.12" 
else:
    isolationCutString = "(pfIsolationR04().sumChargedHadronPt+max(0.,pfIsolationR04().sumNeutralHadronEt+pfIsolationR04().sumPhotonEt-0.5*pfIsolationR04().sumPUPt))/pt< 0.12"

tightMuons = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("pt>20 && isGlobalMuon && isPFMuon && abs(eta)<2.4"
                     " && globalTrack().normalizedChi2<10"
                     " && globalTrack().hitPattern().numberOfValidMuonHits>0"
                     " && globalTrack().hitPattern().numberOfValidPixelHits>0"
                     " && numberOfMatchedStations>1"
                     " && globalTrack().hitPattern().trackerLayersWithMeasurement>5"
                     " && " + isolationCutString
                     )
)


WToMunu = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tightMuons pfMet"),
## Note: the 'mt()' method doesn't compute the transverse mass correctly, so we have to do it by hand.
    cut = cms.string('daughter(1).pt >20  && sqrt(2*daughter(0).pt*daughter(1).pt*(1-cos(daughter(0).phi-daughter(1).phi)))>30'), 
    checkCharge = cms.bool(False),
)



bestWmunu = cms.EDFilter("LargestPtCandViewSelector",
    maxNumber = cms.uint32(10),
    src = cms.InputTag("WToMunu")
)


electronFilter = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(0),
    src = cms.InputTag("looseElectrons")                     
)


muonFilter = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1),
    src = cms.InputTag("looseMuons")                     
)


WSequence = cms.Sequence(tightMuons *
                         WToMunu *
                         bestWmunu
                         )

VetoSequence = cms.Sequence( looseElectrons *
                             electronFilter *
                             looseMuons *
                             muonFilter
                             )

WPath = cms.Sequence(WSequence*VetoSequence)



