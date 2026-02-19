#!/bin/bash
######################################################
# Cody Raul Cardenas 2026.02.19 -- Toussaint Lab
# this script blasts the mitochondrial genes that have
# been recovered against the barcodes using blastn
# this sample is converted into a bed file
# and then them the barcode fasta file is generated 
# using bedtools.
#
# no positional arguments in this script; run like:
# nohup bash 3_blastn.sh > 3_blasnt_20260219.out &!
######################################################

source /local/anaconda3/bin/activate

# useful info thats more digestable than ncbi
# https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
for mtGEN in 3_BARCODES/*.fasta; do
    SAMPLE=$(echo ${mtGEN} | cut -d "/" -f 2 | cut -d "_" -f 1,2);
    conda activate /home/cody/.conda/envs/blastn;
    blastn -db /data/work/Toussaint_UCE/TEST_BARCODES/0_DATA/blast_db/barcode_clean.fasta \
        -query ${mtGEN} \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq salltitles" \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -evalue 1e-25 \
        -num_threads 4 \
        -out 3_BARCODES/${SAMPLE}.blastn;

    # use awk to get the longset barcode; if present
    # probably a quicker way to do this in awk; but
    # it works well enough 
    awk '$4>max { max=$4; out=$0} END {print out}' 3_BARCODES/${SAMPLE}.blastn | cut -f 1,7,8  > 3_BARCODES/${SAMPLE}.bed;
    conda deactivate;

    # use bedtools to extract 
    conda activate /home/cody/.conda/envs/sam-bam-bedtools;
    bedtools getfasta \
        -fi 3_BARCODES/${SAMPLE}_mtGEN.fasta \
        -bed 3_BARCODES/${SAMPLE}.bed \
        > 3_BARCODES/${SAMPLE}_barcode.fasta;
done

# identify empty barcode fasta files and remove them
find ./3_BARCODES/*_barcode.fasta -empty > 3_BARCODES/noBarcode.list
find ./3_BARCODES/*_barcode.fasta -empty -delete