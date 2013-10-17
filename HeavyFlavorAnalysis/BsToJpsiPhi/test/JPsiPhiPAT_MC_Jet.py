import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1), 

)

process.source = cms.Source("PoolSource",
	duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),			    
	fileNames = cms.untracked.vstring(
        #'rfio:/castor/cern.ch/user/g/gfedi/bs4eff/bs4eff100RECO.root',
	#'file:/tmp/terhi/step2_RAW2DIGI_L1Reco_RECO_219_1_HIi.root' #ei toimi 5_3_3_patch1
#	'file:/tmp/terhi/F895D87A-BDDC-E111-A4DA-00259073E346.root' #toimii 5_3_3_patch1

#	'file:/tmp/terhi/TagTestiBplusMC/2EADB8B3-63DF-E111-B174-008CFA0008C4.root' #Bplus toimii
#	'file:/tmp/terhi/58967E7B-47DD-E111-9AAF-E61F13191CAB.root' #Bplus toimii 5_3_3_patch1
#'file:/tmp/terhi/mcData2012/F4E74846-478E-E011-AF4C-0017A4770414.root'#ei toimi 5_3_3_patch1
	#'file:/tmp/gfedi/step2_RAW2DIGI_L1Reco_RECO_27_1_dcf.root'
        #'file:/tmp/gfedi/step2_RAW2DIGI_L1Reco_RECO_219_1_HIi.root',
#'file:/tmp/terhi/DC17936B-47DD-E111-9823-00215E21D690.root'
#'file:/tmp/terhi/C013D1DA-E5DE-E111-9246-001A645C0E8C.root'
'file:/tmp/terhi/BsMCTestFile.root'
	),
	skipEvents = cms.untracked.uint32(0)
)

#from myAnalyzers.JPsiKsPAT.RecoInput2_cfi import *

process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = cms.string('START311_V2::All')
process.GlobalTag.globaltag = cms.string('START53_V7C::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.coreTools import *

process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)

process.load("PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi")
#process.muonMatch.matched = cms.InputTag("genParticlesPlusSim")
process.muonMatch.matched = cms.InputTag("genParticles")
process.muonMatch.maxDeltaR = cms.double(0.02)
process.muonMatch.resolveByMatchQuality = cms.bool(True)

from PhysicsTools.PatAlgos.tools.trackTools import *
#makeTrackCandidates(process, 
#        label='TrackCandsK',                            # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
#        tracks=cms.InputTag('generalTracks'),           # input track collection
#        particleType='K+',                              # particle type (for assigning a mass)
#        preselection='pt>0.3 & abs(eta)<2.5',   # preselection cut on candidates. Only methods of 'reco::Candidate' are available
#        selection='',                           # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
#        isolation={},                                   # Isolations to use ('source':deltaR; set to {} for None)
#        isoDeposits=[],
#        mcAs= None                                     # Replicate MC match as the one used for Muons
#        );                                              #  you can specify more than one collection for this

#makeTrackCandidates(process,
#        label='TrackCandsPi',                             # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
#        tracks=cms.InputTag('generalTracks'),           # input track collection
#        particleType='pi+',                              # particle type (for assigning a mass)
#        preselection='pt>0.3 & abs(eta)<2.5',   # preselection cut on candidates. Only methods of 'reco::Candidate' are available
#        selection='',                           # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
#        isolation={},                                   # Isolations to use ('source':deltaR; set to {} for None)
#        isoDeposits=[],
#        mcAs= None                                     # Replicate MC match as the one used for Muons
#        );



process.allKTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
                                    src = cms.InputTag("generalTracks"),
                                    particleType = cms.string('K+')
                                    )

process.allPiTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
                                    src = cms.InputTag("generalTracks"),
                                    particleType = cms.string('pi+')
                                    )

process.kTracks = cms.EDFilter("CandViewRefSelector",
                               src = cms.InputTag("allKTracks"),
                               #cut = cms.string("pt > 0.3 & abs(eta) < 2.5")
                               cut = cms.string("pt > 0.59 & abs(eta) < 2.5")
                               )

process.piTracks = cms.EDFilter("CandViewRefSelector",
                               src = cms.InputTag("allPiTracks"),
                               #cut = cms.string("pt > 0.3 & abs(eta) < 2.5")
                               cut = cms.string("pt > 0.59 & abs(eta) < 2.5")
                               )


#process.genParticlesPlusSim = cms.EDProducer("GenPlusSimParticleProducer",
#                                             src              = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
#                                             setStatus     = cms.int32(8),
#                                             #   particleTypes = cms.vstring(""),
#                                             filter        = cms.vstring("pt > 0.0"),  # just for testing
#                                             genParticles   = cms.InputTag("genParticles")
#                                             )

#genParticles = cms.EDProducer("GenParticleProducer",
#                  src = cms.InputTag("source")
#)



#process.printTree = cms.EDAnalyzer("ParticleListDrawer",
#  maxEventsToPrint = cms.untracked.int32(-1),
#  printVertex = cms.untracked.bool(False),
#  src = cms.InputTag("genParticles")
#)


#process.printTree = cms.EDAnalyzer( "ParticleTreeDrawer",
# src = cms.InputTag("genParticles"),  
# printP4 = cms.untracked.bool(False),
# printPtEtaPhi = cms.untracked.bool(False),
# printVertex = cms.untracked.bool(False),
# printStatus = cms.untracked.bool(False),
# printIndex = cms.untracked.bool(False),
# status = cms.untracked.vint32( 3 )
#)

#process.p = cms.Path( process.printTree  )

###### PAT JETS ########
import PhysicsTools.PatAlgos.tools.jetTools as jetTools
jetTools.switchJetCollection(process, 
                    cms.InputTag('ak5PFJets'),   
                    doJTA            = True,            
                    doBTagging       = True,            
                    jetCorrLabel     = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute']),  
# data:                   jetCorrLabel     = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']), 
                    doType1MET       = False,            
                    genJetCollection = cms.InputTag("ak5GenJets"),
                    doJetID      = True,
                    jetIdLabel   = "ak5",
                    outputModules = [],
)


###### PAT JETS ########

process.bsVertexAnalysis = cms.EDAnalyzer("BsToJpsiPhiAnalysis",
                                          isMCstudy = cms.bool( True ),
                                          genParticlesLabel  = cms.InputTag("genParticles"),  
                                          TrackLabel_K = cms.InputTag("kTracks"),
                                          TrackLabel_pi = cms.InputTag("piTracks"),
                                          TriggerTag = cms.InputTag("TriggerResults::HLT"),
					  JetCollection = cms.InputTag("patJets"),
                                          MuonTag = cms.InputTag("patMuons"),
                                          StoreDeDxInfo = cms.bool( False ),
                                          JpsiMassWindowBeforeFit = cms.double(0.31), #leave this selection looser than the trigger one for the efficiency calculation
                                          JpsiMassWindowAfterFit = cms.double(0.150), 
                                          JpsiPtCut      = cms.double(6), 
                                          KaonTrackPtCut = cms.double(0.6),
                                          BdKaonTrackPtCut = cms.double(0.6),
                                          PhiMassWindowBeforeFit  = cms.double(0.03), 
                                          PhiMassWindowAfterFit  = cms.double(0.02), 
                                          BsLowerMassCutBeforeFit = cms.double(4.5),
                                          BsUpperMassCutBeforeFit = cms.double(6),
                                          BsLowerMassCutAfterFit  = cms.double(5),
                                          BsUpperMassCutAfterFit  = cms.double(6),
                                          KstarMassWindowBeforeFit =cms.double(0.2),
                                          KstarMassWindowAfterFit =cms.double(0.15),
                                          BdLowerMassCutBeforeFit = cms.double(4.5),
                                          BdUpperMassCutBeforeFit = cms.double(6),
                                          BdLowerMassCutAfterFit = cms.double(4.9),
                                          BdUpperMassCutAfterFit = cms.double(5.7),
                                          verbose                = cms.bool( False ), 
                                          outputFile = cms.untracked.string("/tmp/terhi/BsTests.root"),
                                         )				



##############



###################################################################
###################################################################
# New (easier) Onia2mumu trigger matching
#
#    # Make PAT Muons

process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi
process.patMuonsWithoutTrigger = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.clone()



from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import addMCinfo, changeRecoMuonInput, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution
useL1MatchingWindowForSinglets(process)
addMCinfo(process)
#changeTriggerProcessName(process, "REDIGI36X")
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
#process.muonMatchHLTL3.maxDeltaR = 0.1
#process.muonMatchHLTL3.maxDPtRel = 10.0
#process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
#process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
#process.muonMatchHLTTrackMu.maxDeltaR = 0.1
#process.muonMatchHLTTrackMu.maxDPtRel = 10.0

#
#################################################################
#################################################################


### ==== Apply some final selection (none by default) ====
process.patMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
   # cut = cms.string("p>2 && abs(eta)<2.4"), 
    cut = cms.string("p>0 && abs(eta)<1000"), 
)


#process.patDefaultSequence.remove(process.patJetCorrFactors)
#process.patDefaultSequence.remove(process.patJetCharge)
#process.patDefaultSequence.remove(process.patJetPartonMatch)
#process.patDefaultSequence.remove(process.patJetGenJetMatch)
#process.patDefaultSequence.remove(process.patJetPartons)
#process.patDefaultSequence.remove(process.patJetPartonAssociation)
#process.patDefaultSequence.remove(process.patJetFlavourAssociation)
#process.patDefaultSequence.remove(process.patJets)

#process.patDefaultSequence.remove(process.patMETs)
#process.patDefaultSequence.remove(process.selectedPatJets)
#process.patDefaultSequence.remove(process.cleanPatJets)
#process.patDefaultSequence.remove(process.countPatJets)

# can I do a replace of patMuons with the sequence that includes the trigger matching?
process.patDefaultSequence.replace(process.patMuons,process.patMuonsWithoutTrigger * process.patTriggerMatching * process.patMuons)


process.pat = cms.Path( process.patDefaultSequence )
#process.gp = cms.Path(process.genParticlesPlusSim )
process.ntup = cms.Path( process.allPiTracks * process.allKTracks * process.kTracks * process.piTracks * process.bsVertexAnalysis )

#process.schedule = cms.Schedule(process.gp, process.pat, process.ntup )
process.schedule = cms.Schedule( process.pat, process.ntup )
