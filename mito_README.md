# mitochondrial barcode recovery

## Data preperation
Get reference database for Carabidae

can install or use conda environment. To install, you will need X version
Alternatively, you can use an already installed version in your scripts
e.g.,

```
source /local/anaconda3/bin/activate 
conda activate /home/cody/.conda/envs/sratoolkit
```

Once you've installed or are using an already installed version you need to get your reference genomes and mitochondrial barcodes like so: 

```
conda activate sratoolkit
esearch -db nuccore -query "\"mitochondrion\"[All Fields] AND (\"Carabidae\"[Organism]) AND (refseq[filter] AND mitochondrion[filter] AND (\"12000\"[SLEN] : \"20000\"[SLEN]))" | efetch -format gbwithparts > reference.gb
esearch -db nuccore -query "\"mitochondrion\"[All Fields] AND (\"Carabidae\"[Organism]) AND (refseq[filter] AND mitochondrion[filter] AND (\"12000\"[SLEN] : \"20000\"[SLEN]))" | efetch -format fasta > reference.fasta
esearch -db nuccore -query "barcode (Carabidae[Organism]) 658[SLEN]" | efetch -format fasta > barcode.fasta
```
This can take a little while depending on the number of barcodes.

Clean up barcodes with any gaps using awk
```
# record which sequences are excluded
awk '/^>/{if(NR>1&&p)printf "%s",r; r=$0 ORS; p=0; next} $0~/[Nn-]/{p=1} {r=r $0 ORS} END{if(p)printf "%s",r}' barcode.fasta > barcode_gaps.fasta
# now exclude
awk '/^>/{if(NR>1&&p)printf "%s",r; r=$0 ORS; p=1; next} $0~/[Nn-]/{p=0} {r=r $0 ORS} END{if(p)printf "%s",r}' barcode.fasta > tmp.fasta
# identify taxnomic uncertainty
awk '/^>/{if(NR>1&&p)printf "%s",r; r=$0 ORS; p=($0~/(Carabidae|Coleoptera) sp\./); next} {r=r $0 ORS} END{if(p)printf "%s",r}' tmp.fasta > barcode_taxanomic_uncertainty.fasta
# final clean script
awk '/^>/{if(NR>1&&p)printf "%s",r; r=$0 ORS; p=($0!~/(Carabidae|Coleoptera) sp\./); next} {r=r $0 ORS} END{if(p)printf "%s",r}' tmp.fasta > barcode_clean.fasta
```

Your reference and barcode files should now be ready.

## process sample lists:


!! per library run not total

get all samples and their associated libraries:

```
awk 'FNR==1{lib=FILENAME;sub(/_IDs\.list$/,"",lib)}{print lib "," $1}' *_IDs.list > samples.list
```
run the `2_mito.sh script`

### Isolate mtDNA barcodes

once run is complete, move the data into the appropriate directory and create a symlink like so:
```
make symlink
mkdir 2_MTDATA/mtGENOME
mv CBX*GO/ 2_MTDNA/ && mv CBX*MF/ 2_MTDNA/

find ./2_MTDNA/*_MF/*_MitoFinder_mitfi_Final_Results/ -name '*_final_genes_NT.fasta' -print > MF_results.list
find ./2_MTDNA/*_GO/  -name '*.path_sequence.fasta' -print > GO_results.list

mkdir 3_BARCODES

for SOFTWARE in MF GO; do
    while read LIST; do
        # extract sample ID; adjust cut for subdirectory structure
        SAMPLE_ID=$(echo ${LIST} | cut -d "/" -f 3 | cut -d "_" -f 1)
        # clean up the leading fasta file path for use as a variable
        mtFASTA=$(echo ${LIST} | sed 's|.\/||')
        # always use an explicit path with symlinks
        WORKING="/data/work/Toussaint_UCE/TEST_BARCODES"
        ln -s ${WORKING}/${mtFASTA} ${WORKING}/3_BARCODES/${SAMPLE_ID}_${SOFTWARE}_mtGEN.fasta
    done < ${SOFTWARE}_results.list;
done
```

Once symlinked, blast the samples mtDNA generated with the respective methods using the barcode database.
```
nohup bash 3_blasnt.sh > 3_blastn_20260219.out &!
```



