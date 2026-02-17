# mitochondrial barcode recovery

## Data preperation
Get reference database for Carabidae

can install or use conda environment. To install, you will need X version
Alternatively, you can use an already installed version.
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

!! generate accession list script
```
grep ">" barcode_clean.fasta | cut -d " " -f 1 | sed 's/>//' > accession.list
```

now  generate a barcode genbank file; this one can take a few hours. I recommend writing a shell script and running it with nohup
```
for LIST in $(cat accessions.list); do
esearch -db nuccore -query ${LIST} | efetch -format gbwithparts > ./tmp/${LIST}.gb;
done
cat tmp/*.gb > barcode_clean.gb
```

Your reference barcode files should now be ready.

## process sample lists:


!! per library run not total

get all samples and their associated libraries:

```
awk 'FNR==1{lib=FILENAME;sub(/_IDs\.list$/,"",lib)}{print lib "," $1}' *_IDs.list > samples.list
```


In the genbank file, we need to ensure the gene name is consistent!

grep "/gene="COX1"
should be "/gene="CO1"
(or vise versa, either way all of the genes should be identically named!)



## concatenate your identified mitochondrial barcodes

mitofinder
```
find ./*_barcode/*_barcode_MitoFinder_mitfi_Final_Results/ -name '*_barcode_final_genes_NT.fasta' -print > MF_results.list
```

get organelle
```
find ./*_barcode/*_getorganelle  -name '*.path_sequence.fasta' -print > GO_results.list
```

### Create a symlink to the files

these files have a path to the directory of barcode results like:
```
./CBX2470_barcode/CBX2470_getorganelle/anonym.K85.scaffolds.graph1.1.path_sequence.fasta
./CBX0304_barcode/CBX0304_barcode_MitoFinder_mitfi_Final_Results/CBX0304_barcode_final_genes_NT.fasta
```
we can use these paths to create a symlink those results
```
mkdir barcode_results

for SOFTWARE in GO MF; do
    while read LIST; do
        # extract sample ID
        SAMPLE_ID=$(echo ${LIST} | cut -d "/" -f 2 | cut -d "_" -f 1)
        # clean up the leading fasta file path
        BARCODE=$(echo ${LIST} | sed 's|.\/||')
        # always use an explicit path with symlinks, its safer
        WORKING="/data/work/Toussaint_UCE/2_BARCODES"
        ln -s ${WORKING}/${BARCODE} ${WORKING}/barcode_results/${SAMPLE_ID}_${SOFTWARE}_barcode.fasta
    done < ${SOFTWARE}_results.list;
done
```


