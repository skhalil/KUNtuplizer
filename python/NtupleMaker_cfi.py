import FWCore.ParameterSet.Config as cms

makeNtuples = cms.EDAnalyzer('KUNtuplizer',
    genParts    = cms.InputTag("packedGenParticles"),    
    vertices    = cms.InputTag("offlineSlimmedPrimaryVertices"),
    puInfo      = cms.InputTag("slimmedAddPileupInfo"),
)
