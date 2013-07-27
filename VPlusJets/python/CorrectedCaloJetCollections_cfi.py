import FWCore.ParameterSet.Config as cms
from JetMETCorrections.Configuration.DefaultJEC_cff import *

ak5CaloL2Relative.useCondDB = False
ak5CaloL3Absolute.useCondDB = False
ak5CaloResidual.useCondDB = False

ic5CaloL2Relative.useCondDB = False
ic5CaloL3Absolute.useCondDB = False
ic5CaloResidual.useCondDB = False

kt4CaloL2Relative.useCondDB = False
kt4CaloL3Absolute.useCondDB = False
kt4CaloResidual.useCondDB = False

kt6CaloL2Relative.useCondDB = False
kt6CaloL3Absolute.useCondDB = False
kt6CaloResidual.useCondDB = False


#####################
ic5CaloJetsCorClean = cms.EDProducer("JetViewCleaner",
    srcJets = cms.InputTag("ic5CaloJetsL2L3"),
    module_label = cms.string(""),
    srcObjects = cms.VInputTag(cms.InputTag("gsfElectrons")),
    deltaRMin = cms.double(0.3)
)


#####################
kt4CaloJetsCorClean = cms.EDProducer("JetViewCleaner",
    srcJets = cms.InputTag("kt4CaloJetsL2L3"),
    module_label = cms.string(""),
    srcObjects = cms.VInputTag(cms.InputTag("gsfElectrons")),
    deltaRMin = cms.double(0.3)
)

#####################
kt6CaloJetsCorClean = cms.EDProducer("JetViewCleaner",
    srcJets = cms.InputTag("kt6CaloJetsL2L3"),
    module_label = cms.string(""),
    srcObjects = cms.VInputTag(cms.InputTag("gsfElectrons")),
    deltaRMin = cms.double(0.3)
)

#####################
ak5CaloJetsCorClean = cms.EDProducer("JetViewCleaner",
    srcJets = cms.InputTag("ak5CaloJetsL2L3"),
    module_label = cms.string(""),
    srcObjects = cms.VInputTag(cms.InputTag("gsfElectrons")),
    deltaRMin = cms.double(0.3)
)

##########################################
##########################################
CorJetPath = cms.Sequence(
    ic5CaloJetsL2L3 * ic5CaloJetsCorClean +
    kt4CaloJetsL2L3 * kt4CaloJetsCorClean +
    ak5CaloJetsL2L3 * ak5CaloJetsCorClean
    )

