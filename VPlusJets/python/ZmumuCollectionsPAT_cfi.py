import FWCore.ParameterSet.Config as cms

from ElectroWeakAnalysis.VPlusJets.AllPassFilter_cfi import AllPassFilter

from ElectroWeakAnalysis.VPlusJets.WenuCollectionsPAT_cfi import looseElectrons
from ElectroWeakAnalysis.VPlusJets.WenuCollectionsPAT_cfi import looseMuons
isQCD = False

isolationCutString = cms.string("")
if isQCD:
    isolationCutString = "(pfIsolationR04().sumChargedHadronPt+max(0.,pfIsolationR04().sumNeutralHadronEt+pfIsolationR04().sumPhotonEt-0.5*pfIsolationR04().sumPUPt))/pt> 0.12" 
else:
    isolationCutString = "(pfIsolationR04().sumChargedHadronPt+max(0.,pfIsolationR04().sumNeutralHadronEt+pfIsolationR04().sumPhotonEt-0.5*pfIsolationR04().sumPUPt))/pt< 0.12"

tightMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsPFlow"),
    cut = cms.string("pt>20 && isGlobalMuon && isPFMuon && abs(eta)<2.4"
                     " && globalTrack().normalizedChi2<10"
                     " && globalTrack().hitPattern().numberOfValidMuonHits>0"
                     " && globalTrack().hitPattern().numberOfValidPixelHits>0"
                     " && numberOfMatchedStations>1"
                     " && globalTrack().hitPattern().trackerLayersWithMeasurement>5"
                     " && " + isolationCutString
                     )
)

## tight mu filter
tightMuonFilter = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("tightMuons")                     
)
tightLeptonStep = AllPassFilter.clone()

ZToMuMu = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tightMuons@+ tightMuons@-"),
    #cut = cms.string(' daughter(1).pt >20  && sqrt(2*daughter(0).pt*daughter(1).pt*(1-cos(daughter(0).phi-daughter(1).phi)))>30'), 
    cut = cms.string(' 60 < mass <120 '), 
    checkCharge = cms.bool(True),
)



bestZmumu = cms.EDFilter("LargestPtCandViewSelector",
    maxNumber = cms.uint32(10),
    src = cms.InputTag("ZToMuMu")
)
bestZToMuMuStep = AllPassFilter.clone()


electronFilter = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(0),
    src = cms.InputTag("looseElectrons")                     
)
looseElectronStep = AllPassFilter.clone()


muonFilter = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(2),
    src = cms.InputTag("looseMuons")                     
)
looseMuonStep = AllPassFilter.clone()


ZSequence = cms.Sequence(tightMuons *
                         tightMuonFilter *
                         tightLeptonStep *
                         ZToMuMu *
                         bestZmumu *
                         bestZToMuMuStep
                         )

VetoSequence = cms.Sequence( looseElectrons *
                             electronFilter *
                             looseElectronStep *
                             looseMuons *
                             muonFilter *
                             looseMuonStep
                             )

ZPath = cms.Sequence(ZSequence*VetoSequence)



