import FWCore.ParameterSet.Config as cms

#WP80 electrons, only track iso, remove H/E cut


isQCD = False

tightEleIdLabel = "tight"
looseEleIdLabel = "loose"
if isQCD:
    tightEleIdLabel = "qcd"
    looseEleIdLabel = "qcd"

## modified WP70
tightElectrons = cms.EDProducer("GsfElectronIdSelector",
    src = cms.InputTag( "gsfElectrons" ),
    idLabel = cms.string("tight"),
    useMVAbasedID_   = cms.bool(False)                          
)



WToEnu = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tightElectrons pfMet"),
## Note: the 'mt()' method doesn't compute the transverse mass correctly, so we have to do it by hand.
    cut = cms.string('daughter(1).pt >20  && sqrt(2*daughter(0).pt*daughter(1).pt*(1-cos(daughter(0).phi-daughter(1).phi)))>30'),
    checkCharge = cms.bool(False),
)




bestWToEnu =cms.EDFilter("LargestPtCandViewSelector",
    maxNumber = cms.uint32(10),
    src = cms.InputTag("WToEnu")                 
)


##  Define loose electron selection for veto ######
## modified WP90
looseElectrons = tightElectrons.clone()
looseElectrons.idLabel = cms.string("loose")
looseElectronFilter = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1),
    src = cms.InputTag("looseElectrons")                     
)



##  Define loose muon selection for veto ######
looseMuons = cms.EDFilter("MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("pt>10 &&isPFMuon && (isGlobalMuon || isTrackerMuon) && abs(eta)<2.5"
                     " && (pfIsolationR04().sumChargedHadronPt+max(0.,pfIsolationR04().sumNeutralHadronEt+pfIsolationR04().sumPhotonEt-0.5*pfIsolationR04().sumPUPt))/pt< 0.2")     
)

looseMuonFilter = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(0),
    src = cms.InputTag("looseMuons")                     
)


WSequence = cms.Sequence(tightElectrons *
                         WToEnu *
                         bestWToEnu
                         )
VetoSequence = cms.Sequence( looseElectrons *
                             looseElectronFilter *
                             looseMuons *
                             looseMuonFilter
                             )

WPath = cms.Sequence(WSequence*VetoSequence)

