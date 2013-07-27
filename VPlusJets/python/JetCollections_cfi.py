import FWCore.ParameterSet.Config as cms

mcFlag = False
                 

print "################### Reading Jet Collections ##################"
print "isMC = ",
if mcFlag:
    print "Yup"
else:
    print "Nope"
print "##############################################################"

from ElectroWeakAnalysis.VPlusJets.WenuCollections_cfi import looseElectrons
from ElectroWeakAnalysis.VPlusJets.WenuCollections_cfi import looseMuons
##########################################################################
##########################################################################
##########################################################################
##################### Clean PFJets
##-------------------- Import the Jet RECO modules -----------------------
##-------------------- Turn-on the FastJet density calculation -----------------------
from RecoJets.JetProducers.kt4PFJets_cfi import *
kt6PFJets = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True, voronoiRfact = 0.9)
kt6PFJets.Rho_EtaMax = cms.double(4.5)

# to compute FastJet rho to correct isolation (note: EtaMax restricted to 2.5)
kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True, voronoiRfact = 0.9 )
kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

##--------- Turn-on the FastJet jet area calculation for your favorite algorithm ----
from RecoJets.JetProducers.ak5PFJets_cfi import *
ak5PFJets.Rho_EtaMax = cms.double(4.5)
ak5PFJets.doAreaFastjet = True
##########################################################################
##-------- Remove electrons and muons from jet collection ----------------------
ak5PFJetsClean = cms.EDProducer("PFJetCleaner",
    srcJets = cms.InputTag("ak5PFJets"),
    module_label = cms.string(""),
    idLevel = cms.int32(1),
    etaMin  =  cms.double(0.0),
    etaMax  =  cms.double(2.4),
    ptMin   =  cms.double(0.0),                                 
    srcObjects = cms.VInputTag(cms.InputTag("looseElectrons"),cms.InputTag("looseMuons")),
    deltaRMin = cms.double(0.3)
)
ak5PFJetsCleanVBFTag = cms.EDProducer("PFJetCleaner",
    srcJets = cms.InputTag("ak5PFJets"),
    module_label = cms.string(""),
    idLevel = cms.int32(1),
    etaMin  =  cms.double(2.4),
    etaMax  =  cms.double(9.9),
    ptMin   =  cms.double(0.0),                                 
    srcObjects = cms.VInputTag(cms.InputTag("looseElectrons"),cms.InputTag("looseMuons")),
    deltaRMin = cms.double(0.3)
)
##-------- Apply JetId ----------------------
ak5PFJetsLooseId = cms.EDFilter("PFJetSelector",
    src     = cms.InputTag( "ak5PFJetsClean" ),                                     
    cut = cms.string(
    "neutralHadronEnergyFraction<0.99"
    " && neutralEmEnergyFraction<0.99"
    " && chargedMultiplicity>0"
    " && nConstituents>1"
    " && chargedHadronEnergyFraction>0.0"
    " && chargedEmEnergyFraction<0.99"
    )
)
ak5PFJetsLooseIdVBFTag = cms.EDFilter("PFJetSelector",
    src     = cms.InputTag( "ak5PFJetsCleanVBFTag" ),                                     
    cut = cms.string(
    "neutralHadronEnergyFraction<0.99"
    " && neutralEmEnergyFraction<0.99"
    ##" && chargedMultiplicity>0"
    " && nConstituents>1"
    ##" && chargedHadronEnergyFraction>0.0"
    ##" && chargedEmEnergyFraction<0.99"
    )
)
############################################
PFJetPath = cms.Sequence(
    looseElectrons +
    looseMuons +
    kt6PFJets +
    kt6PFJetsForIsolation +
    ak5PFJets +
    ak5PFJetsClean +
    ak5PFJetsLooseId +
    ak5PFJetsCleanVBFTag +
    ak5PFJetsLooseIdVBFTag
    )
##################### Corrected PFJets
from JetMETCorrections.Configuration.DefaultJEC_cff import *
ak5PFL1Fastjet.srcRho = cms.InputTag('kt6PFJets', 'rho')
if mcFlag:
    ak5PFJetsCor = ak5PFJetsL1FastL2L3.clone()
    ak5PFJetsCorVBFTag = ak5PFJetsL1FastL2L3.clone()
else:
    ak5PFJetsCor = ak5PFJetsL1FastL2L3Residual.clone()
    ak5PFJetsCorVBFTag = ak5PFJetsL1FastL2L3Residual.clone()
    

ak5PFJetsCor.src = "ak5PFJetsLooseId"
ak5PFJetsCorClean = cms.EDFilter("PtMinPFJetSelector",  
    src = cms.InputTag("ak5PFJetsCor"),
    ptMin = cms.double(20.0)
)
ak5PFJetsCorVBFTag.src = "ak5PFJetsLooseIdVBFTag"
ak5PFJetsCorCleanVBFTag = cms.EDFilter("PtMinPFJetSelector",  
    src = cms.InputTag("ak5PFJetsCorVBFTag"),
    ptMin = cms.double(0.0)
)
##########################################
## Filter to require at least two jets in the event
RequireTwoJets = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(100),
    src = cms.InputTag("ak5PFJetsCorClean"),                      
)

CorPFJetPath = cms.Sequence(
    looseElectrons +
    looseMuons + 
    ak5PFJetsClean +
    ak5PFJetsCleanVBFTag +
    ak5PFJetsL2L3 +
    ak5PFJetsL1FastL2L3 +    
    ak5PFJetsL1FastL2L3Residual +
    ak5PFJetsCor +
    ak5PFJetsCorClean +
    ak5PFJetsCorVBFTag +
    ak5PFJetsCorCleanVBFTag +
    RequireTwoJets
    )
if mcFlag:
    CorPFJetPath.remove( ak5PFJetsL1FastL2L3Residual )
##########################################################################
#############  Jets in Monte Carlo  #############
##########################################################################
# ak5GenJets are NOT there: First load the needed modules
from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5GenJets_cfi import *
GenJetPath = cms.Sequence( genParticlesForJets + ak5GenJets )


##################### Tag jets: Needed for MC flavor matching
myPartons = cms.EDProducer("PartonSelector",
    src = cms.InputTag("genParticles"),
    withLeptons = cms.bool(False)
)
###############
ak5flavourByRef = cms.EDProducer("JetPartonMatcher",
    jets = cms.InputTag("ak5PFJets"),
    coneSizeToAssociate = cms.double(0.3),
    partons = cms.InputTag("myPartons")
)

ak5tagJet = cms.EDProducer("JetFlavourIdentifier",
    srcByReference = cms.InputTag("ak5flavourByRef"),
    physicsDefinition = cms.bool(False)
)

#############################################
TagJetPath = cms.Sequence(
    myPartons + 
    ak5flavourByRef*ak5tagJet
    ) 
#############################################
#############################################
