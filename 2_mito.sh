#!/bin/bash
###########################################################################
# Toussaint lab QC and assemble data on pyrgus
# Cody Raul Cardenas; v2026.02
# this script identifies mitochondrial genomes using mitofinder and 
# getorganelle; refernces are generated from NCBI data using a command
# like:
# esearch -db nuccore -query "" | efetch -format fasta > reference.fasta
# and/or 
# esearch -db nuccore -query "" efetch -format gbwithparts > reference.gb
# It expects your sequence files to have been assembled using the
# 0_QC_assembly.sh script
# ensure you have an adequate barcode database and close 
# mitochondrial reference genome for mitofinder
# run like
# nohup bash 2_mito.sh 6 > 2_mito.202602.out &!
#   !!!! Adjust cores to whats reasonable on pyrgus,
#   !!!! check using htop and look at load averages
#   !!!! be kind and leave resources available for others
#   !!!! be sure to check the path to the list of samples
###########################################################################
# KILL ME : [1] 4009
# need to sort out singularity

CORES=$1
source /local/anaconda3/bin/activate 
conda activate /home/cody/.conda/envs/singularity

# sample list should have the library and sample name. e.g., 
# Gaedephaga_lib3,CBX1687
# Gaedephaga_lib3,CBX1691
# Gaedephaga_lib3,CBX1710
# Gaedephaga_lib3,CBX1718
# Gaedephaga_lib3,CBX1730
# Gaedephaga_lib3,CBX1731
# ...

WORKING="/data/work/Toussaint_UCE/TEST_BARCODES" # singularity requires explicit paths, no symlinks!
while read LIST; do
    SAMP=$(echo ${LIST} | cut -d "," -f 2)
    LIB=$(echo ${LIST} | cut -d "," -f 1)
    DATAPATH=${WORKING}/1_ASSEMBLED/${LIB}/spades/${SAMP} # from list
    REFPATH=${WORKING} # where the .gb files live

    singularity run -B ${DATAPATH}:${DATAPATH},${REFPATH}:${REFPATH} \
        ${WORKING}/mitofinder_v1.4.2.sif \
        -j ${SAMP}_MF \
        -r ${REFPATH}/reference.gb \
        -a ${DATAPATH}/contigs.fasta \
        -o 5 \
        -p ${CORES} \
        --rename-contig no;
done < ${WORKING}/samples.list


conda deactivate
conda activate /home/cody/.conda/envs/assembly
# getorganelle
# basically the same process but comparing the two analyses will be useful
WORKING=/data/work/Toussaint_UCE
while read LIST; do
    SAMP=$(echo ${LIST} | cut -d "," -f 2)
    LIB=$(echo ${LIST} | cut -d "," -f 1)
    DATAPATH=${WORKING}/1_ASSEMBLED/${LIB}/spades/${SAMP} # from list
    REFPATH=${WORKING} # where the .gb files live
    
    # removed this flag, it seems useful for barcodes rather than whole mtDNA genomes
    #   --max-extending-len 100 \ 
    get_organelle_from_reads.py \
        -1 ${DATAPATH}/${SAMP}_r1.fastq.gz \
        -2 ${DATAPATH}/${SAMP}_r2.fastq.gz \
        -o ${SAMP}_GO/ \
        -P 0 \
        -F anonym \
        -s ${REFPATH}/reference.fasta \
        --genes ${REFPATH}/reference.fasta \
        --expected-max-size 20000 \
        -t ${CORES};
done < ${WORKING}/samples.list
