import sys, os
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as opts

options = opts.VarParsing ('analysis')

options.register('sample',
     #'/store/mc/RunIISummer16MiniAODv3/SMS-T2bb_mSbot-1650to2600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSummer16v3Fast_pileupfromucsd_94X_mcRun2_asymptotic_v3-v1/80000/F668BBE0-DA1B-E911-8BE0-A0369FD0B3EC.root',
     #'/store/mc/RunIISummer16MiniAODv3/SMS-T2bb_mSbot-1275to2225_dM-25to150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1/90000/F0A575C4-5229-E911-ACE2-0CC47A4C8EA8.root',
     '/store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/FEC6063F-BE88-E811-9D5F-1866DAEB5C78.root', 
     #'/store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/D45302C2-FB87-E811-8C43-F01FAFE15857.root',
     #'/store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/D2A4CF1F-C292-E811-A4B6-008CFA111220.root',
     
     #'/store/mc/RunIIFall17MiniAOD/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/80000/AC04B21B-6212-E811-BC96-0242AC130002.root',
     #'/store/mc/RunIIFall17MiniAODv2/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/00ED846F-7D43-E811-9ED0-20CF305616D1.root',
     #'/store/mc/RunIISummer16MiniAODv2/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/289B3CE4-89B8-E611-89BF-D8D385AE8B08.root',
     #'/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v3/000/284/036/00000/1ABD0A12-619F-E611-AAFC-02163E013674.root',
     opts.VarParsing.multiplicity.singleton,
     opts.VarParsing.varType.string,
     'Sample to tupalize')

options.register('DataProcessing',
     'MC_Fall17MiniAOD',
     opts.VarParsing.multiplicity.singleton,
     opts.VarParsing.varType.string,
     'Data processing types. Options are:\
        Data_94x, MC_Fall17MiniAOD' )

options.register('outputLabel',
     'testNtuple.root',
     opts.VarParsing.multiplicity.singleton,
     opts.VarParsing.varType.string,
     'Output label')

options.register('globalTag',
     '',
     opts.VarParsing.multiplicity.singleton,
     opts.VarParsing.varType.string,
     'Global Tag')

options.register('wantSummary',
      True,
      opts.VarParsing.multiplicity.singleton,
      opts.VarParsing.varType.bool,
     'Want summary report')

options.setDefault('maxEvents', -1)

options.parseArguments()

#if options.DataProcessing == "":
#  sys.exit("!!!!ERROR: Enter 'DataProcessing' period. Options are: Data_94X and MC_Fall17MiniAOD.\n")


if options.globalTag != "": 
  print "!!!!WARNING: You have chosen globalTag as", options.globalTag,". Please check if this corresponds to your dataset."
else: 
  if options.DataProcessing=="Data_94x": options.globalTag="94X_dataRun2_v6"
  elif options.DataProcessing=="MC_Fall17MiniAOD": options.globalTag="94X_mc2017_realistic_v14"
  else: sys.exit("!!!!ERROR: Enter 'DataProcessing' period. Options are: \
      'Data_94X', 'MC_Fall17MiniAOD' \n")

process = cms.Process("Ntupalizer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #options.sample
        #'file:myfile.root'
         '/store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/FEC6063F-BE88-E811-9D5F-1866DAEB5C78.root',
         '/store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/D45302C2-FB87-E811-8C43-F01FAFE15857.root',
         '/store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/D2A4CF1F-C292-E811-A4B6-008CFA111220.root',
        #'/store/mc/RunIISummer16MiniAODv3/SMS-T2bb_mSbot-1650to2600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSummer16v3Fast_pileupfromucsd_94X_mcRun2_asymptotic_v3-v1/80000/F668BBE0-DA1B-E911-8BE0-A0369FD0B3EC.root',
    )
)

# set the global tag 
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = options.globalTag 

# load geometry and electron/photon ID related modules
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi')
process.load('RecoEgamma.ElectronIdentification.ElectronIDValueMapProducer_cfi')

# load the parameters of the ntuplizer
process.makeNtuples = cms.EDAnalyzer('KUNtuplizer')
process.load('Framework.KUNtuplizer.NtupleMaker_cfi')

# load the output root file service routine
process.TFileService = cms.Service("TFileService",
        fileName = cms.string(options.outputLabel)
)

process.p = cms.Path(process.makeNtuples)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(options.wantSummary) )

open('NtupleFileDump.py','w').write(process.dumpPython())
