#!/bin/bash
# **************************************************
# Script to extract promoter sequence of genes or CDS
# Date    :  Dec. 27th, 2022
# Author  :  Xiaodong Li
# **************************************************

# stop the script when an error happens
set -ue

# 1 Get our script name
name=$(basename $0) 
echo 
echo The script name is: $name


# 2 Processing options & parameters with getopts
LENGTH=1500
OUTPREFIX="output"

while getopts hi:g:l:o: opt
do
    case "$opt" in 
    h) 
        echo "Found the -h option"
        echo "Usage: bash $name -g genome_fasta -i input_gff -l 1500 -o out_prefix"
        echo -e "\t-i: input gff file, required"
        echo -e "\t-g: genome sequence with fasta format, required"
        echo -e "\t-l: the upstream length of genes or CDS, default 1500bp"
        echo -e "\t-o: the prefix of output filename, default output"
        echo "`date '+Date: %D %T'`"
        echo
        exit 0
        ;;
    i) 
        echo "Found the -i option, with value $OPTARG"
        GFFFILE=$OPTARG
        ;;
    g) 
        echo "Found the -g option, with value $OPTARG"
        FASTAFILE=$OPTARG
        ;;
    l) 
        echo "Found the -l option, with value $OPTARG"
        LENGTH=$OPTARG
        ;;
    o) 
        echo "Found the -o option, with value $OPTARG"
        OUTPREFIX=$OPTARG
        ;;
    *) 
        echo "Unknown option: $opt"
        ;;
    esac
done

if [ -z $FASTAFILE ]
then
    echo "Not found the -g option, required!"
    exit 1
fi

if [ -z $GFFFILE ]
then
    echo "Not found the -i option, required!"
    exit 1
fi


# 3 Generate gene/cds bed6 file
## A BED file where each feature is described by chrom, start, end, name, score, and strand.
## Here cds coordinates is CDS_start(first cds):CDS_end(last cds)
gffread $GFFFILE --table @chr,@cds,@id,@numexons,@strand | \
awk '{
    gsub("-", ",", $2);
    n=split($2, arr, ",");
    $2=arr[1]"\t"arr[n];
    print $0;
}' OFS="\t" > $OUTPREFIX.cds.bed

gffread $GFFFILE --table @chr,@start,@end,@id,@numexons,@strand > $OUTPREFIX.genes.bed

# 4 Extract promoter location (Define your promoter as the 1.5kb upstream of your gene/cds)
seqkit faidx -f $FASTAFILE

cut -f1,2 $FASTAFILE.seqkit.fai > $OUTPREFIX.genome

bedtools flank \
-i $OUTPREFIX.cds.bed \
-g $OUTPREFIX.genome \
-l 1500 -r 0 -s > $OUTPREFIX.cds.promoters.bed

bedtools flank \
-i $OUTPREFIX.genes.bed \
-g $OUTPREFIX.genome \
-l 1500 -r 0 -s > $OUTPREFIX.genes.promoters.bed

# 5 Get 1500bp up-stream sequences of all genes or CDS, not including gene
seqkit subseq --bed $OUTPREFIX.cds.promoters.bed $FASTAFILE | \
seqkit replace -p "^\S+ " -r '' | \
seqkit sort -o $OUTPREFIX.cds.promoters.fa

seqkit subseq --bed $OUTPREFIX.genes.promoters.bed $FASTAFILE | \
seqkit replace -p "^\S+ " -r '' | \
seqkit sort -o $OUTPREFIX.genes.promoters.fa

exit 0
