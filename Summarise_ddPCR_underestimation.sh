#!/bin/bash -
#title          :Summarise_ddPCR_underestimation
#description    :Summaeise Model_ddPCR_underestimation reuslts
#author         :Roey Angel
#date           :20200921
#version        :1
#usage          :./Summarise_ddPCR_underestimation_V2.sh ${RUNS}
#requires       :
#notes          :
#bash_version   :4.4.20(1)-release
#============================================================================

RUNS=${1:-"../Runs/"} # Runs folder
nGENOMES=${2:-1000} # either parameter 2 or 1000
nFRAGS=${3:-1000} # number of genome fragments to Generate
OLIGOS=${4:-"F338-R805.oligos"} # Oligo file for in silico PCR
RUNS="../Runs/"
RUNDIR=$PWD

# Copy the summary headers
head -n 1 ${RUNS}/Run_1/Results/Summary.csv > ../Summary.csv

for RUN in ${RUNS}/*; do
    tail -n 1 ${RUN}/Results/Summary.csv >> ../Summary.csv
done
