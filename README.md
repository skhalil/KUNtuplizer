# KUNTuplizer
To check out code:
```
setenv SCRAM_ARCH slc6_amd64_gcc630 ; ###chs/ tcsh 
export SCRAM_ARCH=slc6_amd64_gcc630 ; ### bash
cmsrel CMSSW_9_4_4
cd CMSSW_9_4_4/src
cmsenv
git cms-init
git clone git@github.com:skhalil/KUNTuplizer.git Framework/KUNTuplizer
scram b -j4
```
To produce ntuple:
```
cd test
cmsRun NtuplerMaker_cfg.py
```
