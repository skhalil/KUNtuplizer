#!/usr/bin/env python
import subprocess, os, glob, sys, math, itertools
options = [
           ['TTJets', 10],
           ['WJets', 20],
           ['SMS-T2-4bd_genMET-80_mStop-500_mLSP-420', 5],
           ['SMS-T2-4bd_genMET-80_mStop-500_mLSP-490', 5],           
          ]

for opt in options:
    tmpDir = '../input/'+opt[0]
    os.system('mkdir -p '+ tmpDir)
    os.system('mkdir -p out/'+opt[0])
    f_in = open('../input/'+'input_'+opt[0]+'.txt', 'r')
    subFiles = f_in.readlines()       
    print len(subFiles)
    maxJobs = opt[1]
    inputsPerJob = math.ceil(len(subFiles) / maxJobs)
    
    
    for job in range(0, maxJobs):
        #print 'job:', job
        lines = []
        #print 'f_in: ', f_in
        for line in itertools.islice(open('../input/'+'input_'+opt[0]+'.txt', 'r'), job * inputsPerJob, (job + 1) * inputsPerJob):
            lines.append(line.strip())
        if lines:
            with open('../input/'+opt[0]+'/in_'+opt[0]+'_'+str(job)+'.txt', 'w') as output:
                #print '../input/'+opt[0]+'_'+str(job)+'.txt'
                for line in lines: output.write(line+'\n')
                    
    #inCondorFile = open('condor.sh', 'r')
    #readCondorFile = inCondorFile.read()
    #c1 = readCondorFile.replace('TTJets', opt[0])
    #c2 = c1.replace(1, inputsPerJob)
    #fcondor = open('condor_'+opt[0]+'.sh', 'w')
    #fcondor.write(c2)
    #fcondor.close()

    inSubmitFile = open('submit.sh', 'r')
    readSubmitFile = inSubmitFile.read()
    s1 = readSubmitFile.replace('SAMPLE', opt[0])
    s2 = s1.replace('QUEUE', str(maxJobs)) 
    fsubmit = open('submit_'+opt[0]+'.sh', 'w')
    fsubmit.write(s2)
    fsubmit.close()
    condorSubmit = "condor_submit submit_"+opt[0]+".sh"
    subprocess.call( ["echo %s"%condorSubmit, ""], shell=True)
    subprocess.call( [condorSubmit, ""], shell=True)
        
        
