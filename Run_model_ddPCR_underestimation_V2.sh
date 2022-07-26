#!/bin/bash -
#title          :Run_model_ddPCR_underestimation
#description    :Model the chance of getting multiple 16S PCR amplicons from one DNA fragment.
#                Wrapper for running the pipeline n times.
#author         :Roey Angel and Sebastian Vranjes
#date           :20200818
#version        :2
#usage          :./Run_model_ddPCR_underestimation_V2.sh ${nTIMES} ${nGENOMES} ${nFRAGS} ${OLIGOS}
#requires       : [ncbi_genome_download](https://pypi.org/project/ncbi-genome-download/); [SIPSim](https://github.com/nick-youngblut/SIPSim)
#notes          :
#bash_version   :4.4.20(1)-release
#============================================================================

nTIMES=${1:-100} # either parameter 1 or 1000
nGENOMES=${2:-1000} # either parameter 2 or 1000
nFRAGS=${3:-1000} # number of genome fragments to Generate
OLIGOS=${4:-"F338-R805.oligos"} # Oligo file for in silico PCR
RUNS="../Runs/"
RUNDIR=$PWD

# remove existing run folders
#if [ -d $RUNS ]
#then
# rm -rf $RUNS
#fi
#mkdir $RUNS

#if [ -f prokaryotes.txt ]
#then
#    rm -f prokaryotes.txt
#fi
#echo -e "Download list of genome ACC #"
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt

mkdir -p $RUNS
lastRun=$(find ../Runs/* -maxdepth 0 | sed 's/\.\.\/Runs\/Run_//g' | sort -n| tail -n 1)
if [ -z "$lastRun" ]
then
      lastRun=0
fi

# Run Model_ddPCR_underestimation_V2.sh n-times
for i in $(seq $(bc <<< "${lastRun}+1") $(bc <<< "$nTIMES+${lastRun}")); do
    mkdir ${RUNS}/Run_${i}
    cp Model_ddPCR_underestimation_V2.sh Grab_frags.py Plot_fragments.R ${OLIGOS} prokaryotes.txt ${RUNS}/Run_${i}
    cd ${RUNS}/Run_${i}
    bash ./Model_ddPCR_underestimation_V2.sh ${nGENOMES} ${nFRAGS} ${OLIGOS}&
    pids+=" $!"
    cd ${RUNDIR}
done
