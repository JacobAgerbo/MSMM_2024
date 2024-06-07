#!/bin/sh
#module load bbmap
echo "renaming fastq file for $1"
bbrename.sh in1=filt_"$1"_1.fastq.gz in2=filt_"$1"_2.fastq.gz out1=renamed_filt_"$1"_1.fastq.gz out2=renamed_filt_"$1"_2.fastq.gz