#!/bin/sh
#module load fastp
echo "starting QC on $1"
fastp -i "$1"_1.fastq.gz -I "$1"_2.fastq.gz -o filt_"$1"_1.fastq.gz -O filt_"$1"_2.fastq.gz
echo "Done with $1"