#!/bin/bash -
#title          :Model ddPCR underestimation
#description    :Model the chance of getting multiple 16S PCR amplicons from one DNA fragment.
#author         :Roey Angel and Sebastian Vranjes
#date           :20200818
#version        :2
#usage          :./01_Download_genomes.sh ${nGENOMES} ${HMMs}
#requires       : [ncbi_genome_download](https://pypi.org/project/ncbi-genome-download/); [SIPSim](https://github.com/nick-youngblut/SIPSim)
#notes          :
#bash_version   :4.4.20(1)-release
#============================================================================

eval "$(conda shell.bash hook)" # this is needed to be able to use conda activate

nGENOMES=${1:-1000} # either parameter 2 or 1000
nFRAGS=${2:-1000} # number of genome fragments to Generate
OLIGOS=${3:-"F338-R805.oligos"} # Oligo file for in silico PCR
GeneLen=1600 # maximum length of the 16S gene (amplicons must be smaller)
GENOMES="Genomes/" # sotrage dir for downloaded genomes
FRAGMENTS="Fragments/" # sotrage dir for fragments
FINAL="Results/" # sotrage dir 16S matching results
RAW="Raw_results/"
RESOURCES="/proj/Resources/"
# Distrubtion parameters
Xi=2547
Omega=6631
Alpha=5 


# Functions
# download_genomes()
# {
# }

# remove and make folders
if [ -d $GENOMES ]
then
 rm -rf $GENOMES
fi
mkdir $GENOMES

if [ -d $FRAGMENTS ]
then
 rm -rf $FRAGMENTS
fi
mkdir $FRAGMENTS

if [ -d $RAW ]
then
 rm -rf $RAW
fi
mkdir $RAW
mkdir ${RAW}/Genome_hits

if [ -d $FINAL ]
then
 rm -rf $FINAL
fi
mkdir $FINAL

RUNDIR=$PWD
RUNLOG="${RUNDIR}/Run_model_ddPCR_underestimation_$(basename $PWD).log"

echo -e "01 Download genomes" >$RUNLOG
# Download list of prok genomes from NCBI

# 01 Download genomes from NCBI
if [ -f ${FINAL}Random_assembly_ACC.txt ]
then
    rm -f ${FINAL}Random_assembly_ACC.txt
fi
touch ${FINAL}Random_assembly_ACC.txt

GCOUNT=0
# Continue downloading more genomes until reaching nGENOMES
while [ $GCOUNT -lt $nGENOMES ]; do
    conda activate base
    moreGenomes=$(echo "$nGENOMES - $(find ${GENOMES} -name "*.fna" 2> /dev/null | wc -l)" | bc -l) # how many more genomes do we need? (2> /dev/null in case 0 genomes were downloaded so far)
    # Retain only Complete and Chromosome status; retain only unique names; get Assembly accession; grab random $nGENOMES genomes; replace GCA_ with GCF_ (for ncbi-genome-download)
    cat ./prokaryotes.txt  | grep "Complete Genome\|Chromosome" | awk -F '\t' $'!seen[$1]++' | awk -F $'\t' 'BEGIN {OFS = FS} {print $19,$6}' | sort -R | head -n $moreGenomes | sed 's/GCA_/GCF_/g' > ${FINAL}More_assembly_ACC.txt
    ACC=$(cat ${FINAL}More_assembly_ACC.txt | awk -v ORS="," '{print $1}') # store ACC # only
    echo -e "Run ncbi-genome-download" >>${RUNLOG}
    ncbi-genome-download -v --format fasta --assembly-accessions $ACC bacteria --parallel 10 --retries 20 --output-folder $GENOMES --flat-output >>${RUNLOG} 2>&1
    cat ${FINAL}More_assembly_ACC.txt >> ${FINAL}Random_assembly_ACC.txt # list of all the genomes it tried downloading
    downloadedGenomes=$(find ${GENOMES} -type f \( -iname \*.fna -o -iname \*.fna.gz \) -printf '%f\n' | perl -p -e 's/^(GCF_.*?)_.*/$1/' | sort)
    rm ${FINAL}More_assembly_ACC.txt

    # 02 Detect 16S in genomes
    ## run barrnap on each genome (for determining the number of 16S gene copies that could potentially be amplified)
    echo -e "\n02 Detect 16S rRNA genes in genomes using barrnap" >>${RUNLOG}
    cd $GENOMES
    for gz_genome in *.gz; do
      gunzip -f ${gz_genome}
      barrnap -k bac ${gz_genome%.gz} --reject 0.8 -t 10 > ../${RAW}/Genome_hits/${gz_genome%.gz}.barrnap 2>&1
    done
    rm -rf *.gz # for some reason some files retain their gz copy
    conda deactivate

    # # make sure there are no barrnap files without matching genome files (I'm not sure why this happens)
    # mismatch=$(comm -3 <(cat<<<$(find ../${RAW}/Genome_hits/ -type f -name "*.barrnap" |  perl -p -e 's/.*(GCF_.*?)_.*/$1/' | sort)) <(cat<<<$(find ./ -name "*.fna" | perl -p -e 's/.*(GCF_.*?)_.*/$1/' | sort)))
    # if [ ${#mismatch} -gt 0 ]; then
    #   cat <<< $mismatch | while read line; do
    #     echo -e "removing ../${RAW}/Genome_hits/${line}* \n" >>${RUNLOG}
    #     rm -f ../${RAW}/Genome_hits/${line}*
    #     echo -e "../${RAW}/Genome_hits/${line}* \n" >>${RUNLOG}
    #   done
    # fi

    # find and remove genomes without 16S (or 16S length < 0.8)
    nGenome16SHits=$(cat ../${RAW}/Genome_hits/*.barrnap | grep "Name=16S" | wc -l) # total 16S hits in all the genomes (incl. all paralogs)
    nGenomesw16S=$(find ../${RAW}/Genome_hits/ -type f -name "*.barrnap" -exec grep -m 1 "Name=16S" '{}' ';' | wc -l)
    echo -e "barrnap detected at least one 16S rRNA gene in $nGenomesw16S out of $(find ./ -name "*.fna" | wc -l) downloaded genomes" >>${RUNLOG}
    Genomeswo16S=$(find ../${RAW}/Genome_hits/ -type f -name "*.barrnap" -exec grep -HL "Name=16S" '{}' ';' | perl -p -e 's/.*(GCF_.*?)_.*/$1/')
    if [ ${#Genomeswo16S} -gt 0 ]; then
      cat <<< $Genomeswo16S | while read line; do
        # echo -e ${line} >>${RUNLOG}
        rm -f ${line}* 2> /dev/null
        rm -f ../${RAW}/Genome_hits/${line}* 2> /dev/null
      done
    fi

    # run in-silico PCR on the genomes
    echo -e "Run in-silico PCR on each genome" >>${RUNLOG}
    #if [ ! -f PCR_genome_hits.tsv ]; then
          # echo -e  "#seq_name\tbarrnap_V\tGene\tStart\tEnd\tevalue\tstrand\t?\tProduct\n" > PCR_genome_hits.tsv
            echo -e  "Seq name\tACC #\tF position\tF mismatches\tR position\tR mismatches\tAmplicon length" > PCR_genome_hits.tsv
    #fi
    for genome in *.fna; do
      primersearch ${genome} "../${OLIGOS}" 20 -scircular1 Y ${genome%.seq}.match # 10% mismatches from primer allowed
      AMPLIMERS=$(grep Amplimer\ length ${genome%.seq}.match | awk -F $' ' 'BEGIN {OFS = FS} $3<1600' | wc -l) # number of amplicons <1600 bp
      if [ ${AMPLIMERS} -gt 0 ]; then
        cat ${genome}.match \
        | awk -v fname=$(basename ${genome%.fna}) 'BEGIN{RS="Amplimer [0-9]+";FS="\n";};NR>1{printf fname"";for(i=2;i<=NF;i++){printf($i)};print("")}' \
        | awk 'BEGIN{FS="\t"}{$3=""; print}' \
        | awk 'BEGIN{FS=" ";OFS="\t"}{gsub(/[\[,\]]/,"",$18); print $1,$3,$9,$11,$18,$20,$24}1'\
        | awk -v GL="$GeneLen" 'BEGIN {OFS=FS} $7<GL' >> PCR_genome_hits.tsv
      fi
    done

    nGenomePCRHits=$(tail -n+2 PCR_genome_hits.tsv | wc -l) # total 16S-PCR hits in all the genomes (incl. all paralogs)
    nGenomeswPCR=$(awk 'NR==1 {next} !seen[$1]++' PCR_genome_hits.tsv | wc -l) # no. of genomes with at least one 16S PCR hit
    echo -e "primersearch detected at least one 16S rRNA amplicon in ${nGenomeswPCR} out of $(find ../${GENOMES} -name "*.fna" | wc -l) downloaded genomes" >>${RUNLOG}
    pAmplified=$(echo "$nGenomePCRHits/$nGenome16SHits" | bc -l)
    echo -e "${nGenomePCRHits} 16S PCR hits were detected in the genomes out of ${nGenome16SHits} total 16S hits ($(printf '%.3f' ${pAmplified}))" >>${RUNLOG}

    # remove genomes with no 16S match
    goodGenomes=$(cat PCR_genome_hits.tsv | awk '!(NR<=1) {print $1}' | perl -p -e 's/^(GCF_.*?)_.*/$1/' | sort | uniq)
    badGenomes=$(comm -3 <(cat<<<$downloadedGenomes) <(cat<<<$goodGenomes))

    if [ ${#badGenomes} -gt 0 ]; then
      cat <<< $badGenomes | while read line; do
        # echo -e ${line} >>${RUNLOG}
        rm -f ${line}* 2> /dev/null
        rm -f ../${RAW}/Genome_hits/${line}* 2> /dev/null
      done
    fi
    cd ../
    GCOUNT=$(find ${GENOMES} -name "*.fna" | wc -l)
done # end while

fgrep -f <(cat <<<$downloadedGenomes) ${FINAL}/Random_assembly_ACC.txt > ${FINAL}Downloaded_genomes.txt # list of all downloaded genomes
echo -e "$(find ${GENOMES} -name "*.fna"  | wc -l) genomes were downloaded successfully" >>${RUNLOG}
cp ${GENOMES}PCR_genome_hits.tsv ${FINAL}

# Again make sure there are no barrnap files without matching genome files (I'm not sure why this happens)
# mismatch=$(comm -3 <(cat<<<$(find ${RAW}/Genome_hits/ -type f -name "*.barrnap" |  perl -p -e 's/.*(GCF_.*?)_.*/$1/' | sort)) <(cat<<<$(find ./${GENOMES} -name "*.fna" | perl -p -e 's/.*(GCF_.*?)_.*/$1/' | sort)))
# if [ ${#mismatch} -gt 0 ]; then
#   cat <<< $mismatch | while read line; do
#   echo -e "removing ${RAW}/Genome_hits/${line}* \n" >>${RUNLOG}
#   rm -f ${RAW}/Genome_hits/${line}*
#   echo -e "${RAW}/Genome_hits/${line}* \n" >>${RUNLOG}
#   done
# fi

nGenome16SHits=$(cat ${RAW}/Genome_hits/*.barrnap | grep "Name=16S" | wc -l) # total 16S hits in all the genomes (incl. all paralogs)
nGenomesw16S=$(find ${RAW}/Genome_hits/ -type f -name "*.barrnap" -exec grep -m 1 "Name=16S" '{}' ';' | wc -l)
echo -e "barrnap detected at least one 16S rRNA gene in $nGenomesw16S out of $(find ./${GENOMES} -name "*.fna" | wc -l) downloaded genomes" >>${RUNLOG}

# 03 Fragment genomes
conda activate SIPSim
# Generate a list of genome fragments based on a given distribution
echo -e "\n03 Fragment genomes" >>${RUNLOG}
# awk create genome_index.txt (needed for SIPSim fragments)
find ${GENOMES}*.fna | awk 'sub( /.*\//,"",$0 )' | awk -v OFS='\t' '$1=$1"\t"$1' | awk -v OFS='\t' '{gsub(".fna$","",$1); print($1,$2)}' > ${GENOMES}genome_index.txt

# run SIPSim fragments
nFRAGS2grab=$(echo "scale=0; ($nFRAGS * 1.1)/1" | bc) # add 10% because Grab_frags.py fails on some fragments
# output to FRAGMENTS (skewed-normal params are empirically measured, see Fit_selm.html)
SIPSim fragments  \
    $GENOMES/genome_index.txt \
    --fp $GENOMES \
    --fld skewed-normal,${Xi},${Omega},${Alpha} \
    --flr None,None \
    --nf $nFRAGS2grab \
    --np=20 \
    --tbl \
    > ${FRAGMENTS}Frags_table.txt
cp ${FRAGMENTS}Frags_table.txt  ${FINAL}
conda deactivate
echo -e "SIPSim fragments generated $(tail -n+2 ${FRAGMENTS}Frags_table.txt | wc -l) DNA fragments (${nGENOMES} * ${nFRAGS} +10%)" >>${RUNLOG}

# Grab the fragments from the genomes using the table
echo -e "Grab the fragments from the genomes" >>${RUNLOG}
python Grab_frags.py $GENOMES $FRAGMENTS
dDownlodedFrags=$(find ${FRAGMENTS} -name "*.fa" | wc -l)
# delete the unneeded fragments
totalFrags=$(echo "scale=0; (${nFRAGS} * ${nGENOMES})/1" | bc) # total desired fragmnes
if [ ${dDownlodedFrags} -gt ${totalFrags} ]; then
  tooManyFrags=$(echo "scale=0; (${dDownlodedFrags} - ${totalFrags})" | bc)
  find ${FRAGMENTS} -name "*.fa" | sort -R | head -n ${tooManyFrags} | xargs rm
fi

realnFRAGS=$(find ${FRAGMENTS} -name "*.fa" | wc -l)
echo -e "Grab_frags.py grabbed $realnFRAGS DNA fragments" >>${RUNLOG}

# plot the fragment distribution
echo -e "Plot the fragment distribution" >>${RUNLOG}
Rscript Plot_fragments.R

# 04 Detect 16S rRNA genes in the fragments
echo -e "\n04 Detect 16S rRNA genes in the fragments" >>${RUNLOG}
# run barrnap on each fragment
conda activate base
cd ${FRAGMENTS}
for fasta in *.fa; do
    #echo -e "$fasta"
    barrnap -k bac ${fasta} --outseq ../${RAW}/${fasta} -t 5 > ../${RAW}/${fasta%.fa}.barrnap 2>&1
done
conda deactivate

mkdir ../${FINAL}/Fragment_hits
## copy frags with 16S hits only
for rawResult in ../${RAW}/*.fa; do
  if [ -s $rawResult ] && $(cat $rawResult | grep -q 16S); then
    cp ${rawResult} ${rawResult%.fa}.barrnap ../${FINAL}/Fragment_hits
  fi
done
echo -e "barrnap detected at least one 16S rRNA gene in $(for HIT in ../${FINAL}/Fragment_hits/*.barrnap; do cat ${HIT} | grep -m 1 "Name=16S"; done | wc -l) out of $realnFRAGS fragments" >>${RUNLOG}

# Output all hits to files and summarise in a table
cd ../${FINAL}
cat ./Fragment_hits/*.barrnap | grep "Name=16S" | awk 'BEGIN{printf "#seq_name\tbarrnap_V\tGene\tStart\tEnd\tevalue\tstrand\t?\tProduct\n"} {print $0}' > All_hits.tsv

# 05 in silico PCR on fragments
echo -e "\n05 Run in-silico PCR on each fragment" >>${RUNLOG}
echo -e  "Seq name\tFragment\tF position\tF mismatches\tR position\tR mismatches\tAmplicon length" | tee -a PCR_hits.tsv Mult_PCR_hits.tsv
for HIT in ./Fragment_hits/*.fa; do
  primersearch ${HIT} "../${OLIGOS}" 20 -scircular1 N ${HIT%.seq}.match # 20% mismatches from primer allowed
  AMPLIMERS=`grep Amplimer\ length ${HIT%.seq}.match | wc -l`
  if [ ${AMPLIMERS} -gt 0 ]; then
    cat ${HIT}.match \
    | awk -v fname=$(basename ${HIT%.fna}) 'BEGIN{RS="Amplimer [0-9]+";FS="\n";};NR>1{printf fname"";for(i=2;i<=NF;i++){printf($i)};print("")}' \
    | awk 'BEGIN{FS="\t"}{$3=""; print}' \
    | awk 'BEGIN{FS=" ";OFS="\t"}{print $1,$3,$9,$11,$18,$20,$24}1' \
    | awk -v GL="$GeneLen" 'BEGIN {OFS=FS} $7<GL' >> PCR_hits.tsv
  fi
done

# 06 Summarise
echo -e "\n06 Summarise results" >>${RUNLOG}
# fragments with 16S
mean16SperGenome=$(echo -e "$(cat PCR_genome_hits.tsv | wc -l)/${nGENOMES}" | bc -l)
echo -e "Average number of 16S genes per genome: $(printf '%.1f' ${mean16SperGenome})" >>${RUNLOG}
nHITs=$(echo -e "$(find ./ -name "*.fa" | wc -l)" | bc -l)  # count frags with (one or more) 16S  (i.e. once each hit fragment)
echo -e "${nHITs} frags with (one or more) 16S gene copies (i.e. once each hit fragment)" >>${RUNLOG}
pHITs=$(echo -e "${nHITs}/${realnFRAGS}" | bc -l)  # fraction of frags with (one or more) 16S per genome per fragment (i.e. once each hit fragment)
echo -e "Fraction of frags with (one or more) 16S per genome per fragment (i.e. once each hit fragment): $(printf '%.4f' ${pHITs})" >>${RUNLOG}
nMultHits=$(awk 'c[$1]++ && c[$1]==2' All_hits.tsv | wc -l) # count fragments with more than 1 16S hit (i.e. >1 coord lines)
echo -e "${nMultHits} fragments with more than 1 16S hit (i.e. >1 coord lines)" >>${RUNLOG}
pMultHits=$(echo "${nMultHits}/${nHITs}" | bc -l) # fraction of fragments with more than 1 16S hit (i.e. >1 coord lines) per fragments with 16S
echo -e "Fraction of fragments with more than 1 16S hit (i.e. >1 coord lines) per fragments with 16S: $(printf '%.3f' ${pMultHits})" >>${RUNLOG}
# fragments with PCR amplicon
awk 'NR==FNR {count[$1]++; next} count[$1]>1' PCR_hits.tsv PCR_hits.tsv >> Mult_PCR_hits.tsv # grab the fragments with more than one PCR hit per fragment
redundantHits=$(awk 'c[$1]++; c[$1]==2' PCR_hits.tsv | awk 'seen[$1]++' | wc -l) # count redundant hits (# multiple hits excl. first occurrence)
nPCRHits=$(echo -e "($(tail -n +2 PCR_hits.tsv | wc -l)-${redundantHits})" | bc -l) # count ddPCR hits (i.e only one hit per fragment is counted)
echo -e "No. of gene copies identified by ddPCR (PCR hits - hits on the same fragment): ${nPCRHits}" >>${RUNLOG}
pPCRHits=$(echo -e "${nPCRHits}/${nHITs}" | bc -l) # fraction of ddPCR hits per fragments with 16S
echo -e "Fraction of ddPCR hits per fragments with 16S (i.e. fraction of non-split genes): $(printf '%.2f' ${pPCRHits})" >>${RUNLOG}
nMultPCRHits=$(awk 'c[$1]++ && c[$1]==2' Mult_PCR_hits.tsv | wc -l) # count fragments multiple PCR hits (i.e. once per fragment!)
echo -e "${nMultPCRHits} fragments contain more than one 16S gene copy" >>${RUNLOG}
pMultPCRHits=$(echo -e "${nMultPCRHits}/${nPCRHits}" | bc -l) # fraction of multiple PCR hits per fragments with PCR hit
echo -e "Fraction of multiple PCR hits per fragments with a PCR hit: $(printf '%.2f' ${pMultPCRHits})" >>${RUNLOG}
pFragWPCRFragmentsW16S=$(echo -e "${nPCRHits}/${nHITs}" | bc -l) # fraction of PCR hits in total number of fragments with 16S (a measure of how many genes were split)
echo -e "Fraction of of PCR hits in total number of fragments with 16S (a measure of how many genes were split): $(printf '%.2f' ${pFragWPCRFragmentsW16S})" >>${RUNLOG}
redundant16S=$(awk 'c[$1]++; c[$1]==2' PCR_genome_hits.tsv | awk 'seen[$1]++' | wc -l) # count redundant 16S genes in genomes (# multiple 16S excl. first occurrence)
nGenomesWPCR=$(echo -e "($(tail -n +2 PCR_genome_hits.tsv | wc -l)-${redundant16S})" | bc -l) # count ddPCR hits for whole genomes (i.e only one hit per fragment is counted) if not = $nGENOMES then there were PCR mismatches
echo -e "${nGenomesWPCR} ddPCR hits for whole genomes (i.e only one hit per fragment is counted). If not =nGENOMES then there were PCR mismatches" >>${RUNLOG}

# Correct ddPCR value
ddPCRcFactor=$(echo -e "(1/${pFragWPCRFragmentsW16S} + ${redundantHits}/${nPCRHits})" | bc -l) # correct for fraction of split genes + for multiple amplicons per fragment
cddPCR=$(echo -e "${nPCRHits}*($ddPCRcFactor)" | bc -l)

# Create a summery table
echo -e "# Genomes,# Frags per genome,Total # frags,Mean 16S per genome,Frags w. 16S,P frags w. 16S,Frags w. PCR (ddPCR outcome), P frags w. PCR,Frags w. multi 16S,P frags w. multi 16S,Frags w. multi PCR,P frags w. multi PCR,ddPCR correction factor, Corrected ddPCR value" > Summary.csv
echo -e "${nGENOMES},`printf '%.0f' $(echo "$realnFRAGS / $nGENOMES" | bc -l)`,$(echo "$nGENOMES * $realnFRAGS" | bc -l),`printf '%.1f' ${mean16SperGenome}`,`printf '%.0f' ${nHITs}`,`printf '%.4f' ${pHITs}`, `printf '%.0f' ${nPCRHits}`,`printf '%.4f' ${pPCRHits}`, `printf '%.0f' ${nMultHits}`,`printf '%.4f' ${pMultHits}`,`printf '%.0f' ${nMultPCRHits}`,`printf '%.4f' ${pMultPCRHits}`,`printf '%.3f' ${ddPCRcFactor}`,`printf '%.0f' ${cddPCR}`" >> Summary.csv
if [ $(cat Summary.csv | wc -l | bc -l) == 2 ]; then
    echo -e "Results written to file successfully" >>${RUNLOG}
fi

cd ../
mv ${RUNLOG} ${FINAL}
rm -rf Fragments
zip -rm Genomes.zip Genomes
rm -rf ${RAW}

# [1] https://github.com/kblin/ncbi-genome-download
# [2] Seemann T. barrnap 0.9 : rapid ribosomal RNA prediction, https://github.com/tseemann/barrnap
# [3] SIPSim
