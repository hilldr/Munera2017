
#! /bin/bash
GENES=/data/genomes/hg19_genes_refFlat.gtf
GENOME=/data/genomes/hg19.fa
INDEX=/data/genomes/hg19

#FASTQC quality control report generator
mkdir ../DATA/QC
for file in ../DATA/FASTQ/*.fastq* # will output filename as "$file"
do
    FILENAME="$file"       
    fastqc --outdir=../DATA/QC $FILENAME
done

#! /bin/bash
mkdir ../DATA/BAM
for file in ../DATA/FASTQ/*.fastq*
do
    SHORTNAME=$(basename "$file")
    NAME2="${SHORTNAME##*/}"
    DIRNAME="${NAME2%.*}"
    tophat2 -p 8 --b2-very-sensitive --no-coverage-search --no-novel-juncs --GTF $GENES -o ../DATA/BAM/$DIRNAME $INDEX $file
done

#! /bin/bash
for d in ../DATA/BAM/*/
do
    FILENAME="$file"       #set variable FILENAME equal to file from line 1
    SHORTNAME=$(basename "$file")
    NAME2="${SHORTNAME##*/}"
    DIRNAME="${d}"
    cufflinks -p 8 -o $DIRNAME --multi-read-correct --compatible-hits-norm --upper-quartile-norm --GTF $GENES ${d}*hits.bam
done

#! /bin/bash
cuffmerge -g $GENES -s $GENOME -p 8 -o ../DATA/merged_asm gtf_assembly.txt
for d in ../DATA/BAM/*/
do
    FILENAME="$file"       #set variable FILENAME equal to file from line 1
    SHORTNAME=$(basename "$file")
    NAME2="${SHORTNAME##*/}"
    DIRNAME="${d}"
    cuffquant -p 8 -o $DIRNAME --max-mle-iterations 100000 -v --multi-read-correct ../DATA/merged_asm/merged.gtf ${d}*hits.bam
done

mkdir ../RESULTS
# CUFFNORM
cuffnorm -o ../RESULTS/normout -p 8 -L Noggin-TXP,BMP-TXP,adult-SI,adult-colon ../DATA/merged_asm/merged.gtf \
../DATA/BAM/run1205_lane12_read1_indexD708-D501=NogginTXP9.fastq/abundances.cxb,../DATA/BAM/run1205_lane12_read1_indexD708-D503=NogginTXP5.fastq/abundances.cxb,../DATA/BAM/run1205_lane12_read1_indexD708-D504=NogginTXP4.fastq/abundances.cxb \
../DATA/BAM/run1205_lane12_read1_indexD709-D502=5331BMPTXP.fastq/abundances.cxb,../DATA/BAM/run1205_lane12_read1_indexD709-D503=5330BMPTXP.fastq/abundances.cxb,../DATA/BAM/run1205_lane12_read1_indexD709-D507=4904BMPTXP.fastq/abundances.cxb \
../DATA/BAM/ERR315344_merged/abundances.cxb,../DATA/BAM/ERR315381_merged/abundances.cxb,../DATA/BAM/ERR315409_merged/abundances.cxb,../DATA/BAM/ERR315442_merged/abundances.cxb,../DATA/BAM/ERR315461_merged/abundances.cxb \
../DATA/BAM/ERR315348_merged/abundances.cxb,../DATA/BAM/ERR315357_merged/abundances.cxb,../DATA/BAM/ERR315484_merged/abundances.cxb
# Run R script
Rscript *.R
