import FWCore.ParameterSet.Config as cms

##########################################################################

# Apply loose PileUp PF jet ID
ak5PFnoPUJets = cms.EDProducer("PATPuJetIdSelector",
    src = cms.InputTag( "selectedPatJetsPFlow" ),
    idLabel = cms.string("loose"),
    valueMapLabel = cms.string("puJetMvaChs")
)

# Apply loose PF jet ID
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
ak5PFGoodJets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
     filterParams = pfJetIDSelector.clone(),
     src = cms.InputTag("ak5PFnoPUJets"),
     filter = cms.bool(True)
)

##-------- Remove electrons and muons from jet collection ----------------------
ak5PFJetsClean = cms.EDProducer("PFPATJetCleaner",
    srcJets = cms.InputTag("ak5PFGoodJets"),
    module_label = cms.string(""),
    srcObjects = cms.VInputTag(cms.InputTag("looseElectrons"),cms.InputTag("looseMuons")),
    deltaRMin = cms.double(0.3)
)


ak5PFJetsLooseIdAll = cms.EDFilter("PATJetRefSelector",
    src = cms.InputTag("ak5PFJetsClean"),
    cut = cms.string('pt > 20.0')
)

ak5PFJetsLooseId = cms.EDFilter("PATJetRefSelector",
    src = cms.InputTag("ak5PFJetsClean"),
    cut = cms.string('pt > 20.0 && abs(eta) < 2.4')
)

ak5PFJetsLooseIdVBFTag = cms.EDFilter("PATJetRefSelector",
    src = cms.InputTag("ak5PFJetsClean"),
    cut = cms.string('pt > 20.0 && abs(eta) > 2.4 && abs(eta) < 9.9')
)

##########################################
## Filter to require at least two jets in the event
RequireTwoJets = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(100),
    src = cms.InputTag("ak5PFJetsLooseId"),
)


############################################
PFJetPath = cms.Sequence( ak5PFnoPUJets + ak5PFGoodJets + ak5PFJetsClean + ak5PFJetsLooseId + ak5PFJetsLooseIdAll +
	ak5PFJetsLooseIdVBFTag )
#	ak5PFJetsLooseIdVBFTag + RequireTwoJets )
##########################################
##########################################################################
#############  Jets in Monte Carlo  #############
##########################################################################
# ak5GenJets are NOT there: First load the needed modules
from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5GenJets_cfi import *
from RecoMET.Configuration.GenMETParticles_cff import *
from RecoMET.METProducers.genMetTrue_cfi import *
 
GenJetPath = cms.Sequence( genParticlesForJets + ak5GenJets + genParticlesForMETAllVisible + genMetTrue)

##################### Tag jets: Needed for MC flavor matching
myPartons = cms.EDProducer("PartonSelector",
    src = cms.InputTag("genParticles"),
    withLeptons = cms.bool(False)
)

ak5flavourByRef = cms.EDProducer("JetPartonMatcher",
    jets = cms.InputTag("selectedPatJetsPFlow"),
    coneSizeToAssociate = cms.double(0.3),
    partons = cms.InputTag("myPartons")
)

ak5tagJet = cms.EDProducer("JetFlavourIdentifier",
    srcByReference = cms.InputTag("ak5flavourByRef"),
    physicsDefinition = cms.bool(False)
)

TagJetPath = cms.Sequence( myPartons +  ak5flavourByRef*ak5tagJet) 
#############################################

