#!/bin/bash

#This example is for the case of paied end ChIP-seq

#This pipeline performs;
#1. adaptor trimming (trimmomatic)
#2. QC (fastqc)
#3. mapping (bowtie)
#4. duplicate removal (picard).

index="XXX" #Path to your index
adapter="XXX/Trimmomatic-0.39/adapters/TruSeq3-PE.fa" #Path to your adaptor used by trimmomatic

dname="H3K27Ac_TNFa_0min" #Your condition name

dout="XXX/${dname}" #Path to your output directory named ${dname}
nthreads=4 #Number of threads used by trimmomatic and fastqc

# Your fastq files 
forward_fq="XXX.fastq.gz"
reverse_fq="XXX.fastq.gz"


# trimming adapters
java -jar trimmomatic-0.39.jar PE \
-threads ${nthreads} \
-phred33 \
${forward_fq} ${reverse_fq} \
${dout}/${dname}_forward_paired.fq.gz ${dout}/${dname}_forward_unpaired.fq.gz \
${dout}/${dname}_reverse_paired.fq.gz ${dout}/${dname}_reverse_unpaired.fq.gz \
ILLUMINACLIP:${adapter}:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:30 


# quality check
fastqc --threads ${nthreads} --nogroup -o ${dout} ${dout}/${dname}_forward_paired.fq

fastqc --threads ${nthreads} --nogroup -o ${dout} ${dout}/${dname}_reverse_paired.fq


# mapping using Bowtie

bowtie -S --best -p 4 ${index} -1 ${dout}/${dname}_forward_paired.fq -2 ${dout}/${dname}_reverse_paired.fq ${dout}/${dname}.sam


# sorting
samtools view -b -o ${dout}/${dname}.bam ${dout}/${dname}.sam
samtools sort ${dout}/${dname}.bam -o ${dout}/${dname}_sorted.bam

# use only mapped reads
samtools view -@ 8 -b -F 4 ${dout}/${dname}_sorted.bam > ${dout}/${dname}_tmp.bam
# dupulicate removal
java -jar picard.jar MarkDuplicates \
I=${dout}/${dname}_tmp.bam \
O=${dout}/${dname}_rm_dups.bam \ #This is your final processed file
M=${dout}/${dname}_picard_report.txt \
REMOVE_DUPLICATES=true

samtools index ${dout}/${dname}_rm_dups.bam 






