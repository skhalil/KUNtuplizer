#!/bin/bash



echo "input parameters: cluster, process, inpath, outpath," $1 $2 $3 $4 #Don't need inpath because we hardcoded the path into the text file

echo

CLUSTER=$1

PROCESS=$2

INPATH=$3

OUTPATH=$4

START_TIME=`/bin/date`

echo "started at $START_TIME"

#cd ${_CONDOR_SCRATCH_DIR}

echo "executing ..."
cd /home/t3-ku/skhalil/CMSSW_10_1_4_patch1/src/
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

cd  ${_CONDOR_SCRATCH_DIR}

counter=0

for txt in `ls $INPATH`; do
    echo "going here: "$txt
    if test -f $INPATH/$txt; then
       echo "file   $INPATH/$txt   exists"
       name=`basename $txt .txt`
       name+=".root"
       echo "output file name: ${name/in/out}"
        if test $counter -eq $PROCESS; then
        echo "python analysis.py --files $INPATH/$txt --outFile $OUTPATH/${name/in/out}"
        python analysis.py --files $INPATH/$txt --outFile $OUTPATH/${name/in/out}
        fi
        let "counter+=1"
    fi
done


#echo "python analysis.py --files $INPATH/TTJets.txt --outFile $OUTPATH/TTJets_out.root" 
#python analysis.py --files $INPATH/TTJets.txt --outFile $OUTPATH/"TTJets_out.root"

exitcode=$?



echo ""

END_TIME=`/bin/date`

echo "finished at ${END_TIME}"

