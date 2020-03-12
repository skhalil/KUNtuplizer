universe = vanilla
executable = condor.sh
#inPath = /mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/skhalil/
#inPath = /home/t3-ku/skhalil/CMSSW_10_1_4_patch1/src/
inPath = /home/t3-ku/skhalil/CMSSW_10_1_4_patch1/src/Framework/KUNtuplizer/test/analysis/input/SAMPLE
outPath = /home/t3-ku/skhalil/CMSSW_10_1_4_patch1/src/Framework/KUNtuplizer/test/analysis/condor/out
logDir = /home/t3-ku/skhalil/CMSSW_10_1_4_patch1/src/Framework/KUNtuplizer/test/analysis/condor/out/log
arguments = $(cluster) $(process) $(inPath) $(outPath)
output = $(logDir)/batch_$(cluster)_$(process).stdout
error =  $(logDir)/batch_$(cluster)_$(process).stderr
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = analysis.py
getenv = True
queue QUEUE 
