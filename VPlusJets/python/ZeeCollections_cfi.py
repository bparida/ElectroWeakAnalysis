import FWCore.ParameterSet.Config as cms

#WP80 electrons, only track iso, remove H/E cut
selectElectrons = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag( "gsfElectrons" ),
    cut = cms.string(
    "(ecalDrivenSeed==1) && (abs(superCluster.eta)<2.5)"
    " && !(1.4442<abs(superCluster.eta)<1.566)"
    " && (et>20.0)"
    " && (gsfTrack.trackerExpectedHitsInner.numberOfHits==0 && !(-0.02<convDist<0.02 && -0.02<convDcot<0.02))"
    " && (dr03TkSumPt/p4.Pt <0.1)"    
    " && ((isEB"
    " && (sigmaIetaIeta<0.01)"
    " && ( -0.06<deltaPhiSuperClusterTrackAtVtx<0.06 )"
    " && ( -0.004<deltaEtaSuperClusterTrackAtVtx<0.004 )"
    ")"
    " || (isEE"
    " && (sigmaIetaIeta<0.03)"
    " && ( -0.03<deltaPhiSuperClusterTrackAtVtx<0.03 )"
    " && ( -0.007<deltaEtaSuperClusterTrackAtVtx<0.007 )"
    "))"
    )
)





ZToEE = cms.EDProducer("NamedCandViewShallowCloneCombiner",
    cut = cms.string('60 < mass < 120'),
    name = cms.string('ZToEE'),
    roles = cms.vstring('electron1', 'electron2'),
    decay = cms.string('selectElectrons@+ selectElectrons@-')
)

bestZee = cms.EDFilter("LargestPtCandViewSelector",
    maxNumber = cms.uint32(1),
    src = cms.InputTag("ZToEE")
)

ZPath = cms.Sequence(selectElectrons+ZToEE+bestZee)

