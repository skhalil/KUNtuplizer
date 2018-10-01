import sys, os
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as opts

options = opts.VarParsing ('analysis')

options.register('sample',
     '/store/mc/RunIISummer16MiniAODv2/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/289B3CE4-89B8-E611-89BF-D8D385AE8B08.root',
     #'/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v3/000/284/036/00000/1ABD0A12-619F-E611-AAFC-02163E013674.root',
     opts.VarParsing.multiplicity.singleton,
     opts.VarParsing.varType.string,
     'Sample to tupalize')

options.register('outputLabel',
     'testNtuple.root',
     opts.VarParsing.multiplicity.singleton,
     opts.VarParsing.varType.string,
     'Output label')

options.register('wantSummary',
      True,
      opts.VarParsing.multiplicity.singleton,
      opts.VarParsing.varType.bool,
     'Want summary report')

options.setDefault('maxEvents', 100)

options.parseArguments()


process = cms.Process("Ntupalizer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        options.sample
        #'file:myfile.root'
    )
)

process.makeNtuples = cms.EDAnalyzer('KUNtuplizer')
process.load('Framework.KUNtuplizer.NtupleMaker_cfi')

# load the output root file service routine
process.TFileService = cms.Service("TFileService",
        fileName = cms.string(options.outputLabel)
)

process.p = cms.Path(process.makeNtuples)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(options.wantSummary) )

open('NtupleFileDump.py','w').write(process.dumpPython())
