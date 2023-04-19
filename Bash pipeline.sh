#!/bin/bash

#Have the directories set:
#dnaseq/
#  ├── results/
#  ├── data/
#        ├── untrimmed_fastq
#        ├── trimmed_fastq
#        ├── reference
#        ├── aligned_data
#annotation file in reference: wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed
#have the fastqs in untrimmed_fastq: wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
#                                    wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
#have hg19 reference genome in reference: wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz


echo "Running"

#Change to working directory
cd ~/dnaseq_pipeline

#Correct the extension on fastq files and decompress
cd data/untrimmed_fastq

for file in *.qz; do
    mv "$file" "${file/.qz/.gz}"
done

for file in *.gz; do
    gzip -d "$file"
    done

#Fastqc on raw data
for file in *.fastq; do
    fastqc "$file"
    done

for file in *html *zip; do 
    mv "$file" ~/dnaseq_pipeline/results/ 
    done

#Trimming
trimmomatic PE  \
-phred33 \
*R1.fastq *R2.fastq \
-baseout /home/ubuntu/dnaseq_pipeline/data/trimmed_fastq/NGS0001.trimmed.fastq \
ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
TRAILING:25 MINLEN:50

#Removing old files
rm *fastq

#Change working directory
cd ..
cd trimmed_fastq

#Post-trimming fastqc
for file in *P.fastq; do
    fastqc "$file"
    done

for file in *html *zip; do 
    mv "$file" ~/dnaseq_pipeline/results/ 
    done

#Change working directory
cd ..
cd reference

#Decompress and index hg19.fa
gzip -d hg19.fa.gz

bwa index -a bwtsw hg19.fa

#Change working directory
cd ..
cd aligned_data

#Alignment
bwa mem -t 4 -v 1 -R \
'@RG\tID:HiSeq2500\tPL:ILLUMINA\tSM:NGS0001\tLB:Nextera\tPU:111' -I 250,50 \
~/dnaseq_pipeline/data/reference/hg19.fa \
~/dnaseq_pipeline/data/trimmed_fastq/NGS0001.trimmed_1P.fastq \
~/dnaseq_pipeline/data/trimmed_fastq/NGS0001.trimmed_2P.fastq > \
~/dnaseq_pipeline/data/aligned_data/NGS0001.hg19.bwa.raw.sam

#Removing old files
cd ..
cd trimmed_fastq
rm *fastq
cd ..
cd aligned_data

#Sort and convert to bam
samtools sort -O bam -o NGS0001.hg19.bwa.sorted.bam NGS0001.hg19.bwa.raw.sam

#Index bam
samtools index NGS0001.hg19.bwa.sorted.bam

#Duplicate marking
picard MarkDuplicates I=NGS0001.hg19.bwa.sorted.bam \
O=NGS0001.hg19.bwa.marked_dup.bam \
CREATE_INDEX=true \
M=NGS0001.hg19.bwa.marked_dup_metrics.txt 

#Filter bam file
samtools view -F 1796  -q 20 -h -b \
NGS0001.hg19.bwa.marked_dup.bam > NGS0001.hg19.bwa.filtered.bam

#Index bam file
samtools index NGS0001.hg19.bwa.filtered.bam

#Variant calling
freebayes --bam NGS0001.hg19.bwa.filtered.bam \
--fasta-reference ~/dnaseq_pipeline/data/reference/hg19.fa \
--vcf ~/dnaseq_pipeline/results/NGS0001.hg19.bwa.freebayes.raw.vcf

#Change working directory
cd ../..
cd  results

#Compress and tabix
bgzip NGS0001.hg19.bwa.freebayes.raw.vcf
tabix -p vcf NGS0001.hg19.bwa.freebayes.raw.vcf.gz

#Variant filter
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
NGS0001.hg19.bwa.freebayes.raw.vcf.gz > NGS0001.hg19.bwa.freebayes.filtered.vcf

#Filter with bed file
bedtools intersect -header -wa -a NGS0001.hg19.bwa.freebayes.filtered.vcf -b ~/dnaseq_pipeline/data/reference/annotation.bed \
> NGS0001.hg19.bwa.freebayes.filtered.bed.vcf

#Compress and tabix
bgzip NGS0001.hg19.bwa.freebayes.filtered.bed.vcf
tabix -p vcf NGS0001.hg19.bwa.freebayes.filtered.bed.vcf.gz

#Variant annotation
~/annovar/./convert2annovar.pl -format vcf4 NGS0001.hg19.bwa.freebayes.filtered.bed.vcf.gz \
> NGS0001.hg19.bwa.freebayes.filtered.bed.avinput

~/annovar/./table_annovar.pl NGS0001.hg19.bwa.freebayes.filtered.bed.avinput ~/annovar/humandb/ -buildver hg19 \
-out ~/dnaseq_pipeline/results/NGS0001.hg19.bwa.freebayes.filtered.bed -remove \
-protocol  refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout

echo "Done"
