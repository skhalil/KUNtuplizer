from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import os
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException
from multiprocessing import Process

PPNdir = os.environ['CMSSW_BASE']+'/src/Framework/KUNtuplizer/'


batch = "signal"
#batch = "background"

if "signal" in batch:
	mydatasets = 'crab_input_files_signal.txt'
	filesplit = 'FileBased'
	units = 1	
if "background" in batch:
        mydatasets = 'crab_input_files_background.txt'
        filesplit = 'FileBased'
	units = 1

config = config()

config.section_('General')
config.General.requestName = None
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = batch

config.section_('JobType')
config.JobType.psetName = PPNdir+'test/NtupleMaker_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.maxJobRuntimeMin = 2500#2000
config.JobType.maxMemoryMB = 4000 #2500
#config.JobType.allowUndistributedCMSSW = True
#config.JobType.disableAutomaticOutputCollection = True

config.section_('Data')
config.Data.inputDataset = None
config.Data.secondaryInputDataset = None
config.Data.inputDBS = 'global'
config.Data.splitting = filesplit
config.Data.unitsPerJob = units
#config.Data.outLFNDirBase = '/store/group/lpcbprime/noreplica/skhalil/Upgrade'
config.Data.outLFNDirBase = '/store/user/skhalil/susySVNtuples/'
#physical location:  /mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/group/skhalil/upgradeNtuples/
config.Data.allowNonValidInputDataset = True
#config.Data.ignoreLocality = False
config.Data.publication = False

config.section_('Site')
config.Site.storageSite = 'T2_US_Nebraska' #'T3_US_FNALLPC'

def submit(config):
    try:
        crabCommand('submit', config = config)
    except HTTPException as hte:
        print "Failed submitting task: %s" % (hte.headers)
    except ClientException as cle:
        print "Failed submitting task: %s" % (cle)

datasetsFile = open( mydatasets )
jobsLines = datasetsFile.readlines()
for ijob in jobsLines :
    s = ijob.rstrip()
    if (len(s)==0 or s[0][0]=='#'): continue
    #print 's: ', s
    #cdi = s + ''
    cgr = s.split('/')[1]
    #print 'cgr: ', cgr  
    if 'ext1' in s: cgr_v1 = cgr+'_ext1' 
    else: cgr_v1 = cgr 
    config.Data.inputDataset = s
    config.General.requestName = cgr_v1
    print "Submitting to Crab:"
    print "Inputdataset: ",s
    print "requestName: ",cgr_v1
    print 
    submit(config)
