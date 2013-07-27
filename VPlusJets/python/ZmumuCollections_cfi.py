import FWCore.ParameterSet.Config as cms


selectMuons = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("pt>20 && isGlobalMuon && isTrackerMuon && abs(eta)<2.4"
                     " && globalTrack().normalizedChi2<10"
                     " && innerTrack().numberOfValidHits>10"
                     " && globalTrack().hitPattern().numberOfValidMuonHits>0"
                     " && globalTrack().hitPattern().numberOfValidPixelHits>0"
                     " && numberOfMatches>1"
                     " && (isolationR03().sumPt+isolationR03().emEt+isolationR03().hadEt)/pt< 0.3"
                     )
)



ZToMM = cms.EDProducer("NamedCandViewShallowCloneCombiner",
    cut = cms.string('60 < mass < 120'),
    name = cms.string('ZToMM'),
    roles = cms.vstring('muon1', 'muon2'),
    decay = cms.string('selectMuons@+ selectMuons@-'),
   checkCharge = cms.bool(True)                    
)

bestZmumu = cms.EDFilter("LargestPtCandViewSelector",
    maxNumber = cms.uint32(1),
    src = cms.InputTag("ZToMM")
)

ZPath = cms.Sequence(selectMuons*ZToMM*bestZmumu)

