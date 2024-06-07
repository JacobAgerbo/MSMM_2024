# Using Anvi'o for Genome-Resolved Metagenomics and other 'omics

## Introduction
Anvi'o is a powerful tool for analyzing metagenomic data and performing genome-resolved metagenomics. In this guide, we will walk through the steps to utilize Anvi'o for a functional metagenomics analysis.

## Installation
To get started with Anvi'o, you need to install it on your system. Follow the installation instructions provided on the official Anvi'o website.

For more information and detailed documentation, refer to the [Anvi'o documentation](https://merenlab.org/software/anvio/).

## Data Preparation
Before running Anvi'o, make sure you have the necessary metagenomic data ready. In our case, we use an EU-based data portal for multi omics analysis of host-microbe interactions in domesticated chicken and Atlantic salmon. This portal includes raw sequencing reads, assembly files, and any other relevant data for your analysis.

For more information and detailed documentation, refer to the [Data portal](https://www.holofooddata.org/).

From the portal we can get ENA accessions and MAGs. Now lets download samples and the MAGs.

```bash
# Get samples
bash get_data.sh SAMPLES.txt

# Get MAGs (sorry for the ugly code)
while IFS= read -r genome; do
    # Process each genome here 
    echo "Downloading MAG: $genome"
    wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/non-model-fish-gut/v2.0/species_catalogue/MGYG0003075/"$genome"/genome/"$genome".fna
    wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/non-model-fish-gut/v2.0/species_catalogue/MGYG0002995/"$genome"/genome/"$genome".fna
    wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/non-model-fish-gut/v2.0/species_catalogue/MGYG0002996/"$genome"/genome/"$genome".fna
done < MAGs.txt
```

Now we will use anvi'o to get functional and taxonomic information, and microbial dynamics across samples. In short, we will do as below.
## Running Anvi'o
1. **Initialize Anvi'o:** Start by initializing Anvi'o and setting up your analysis environment.

    Use link to install anvi'o. [Anvi'o installation](https://anvio.org/install/).

2. **Import Data:** Import your metagenomic data into Anvi'o for analysis.

    ### Make contig database
    Lets start with merging all MAGs into one contig database (DB).
    ```bash
    cat *.fna > contigs-fixed.fa
    ```
    ...and index the contig database for re-mapping.
    ```bash
    bwa index contigs-fixed.fa
    ```
    Now lets make the anvi'o contig database from our HoloFood MAGs.
    ```bash
    anvi-gen-contigs-database -f contigs-fixed.fa -o ./CONTIGS.db -n 'CONTIG DB for HoloFood MAGs'
    ```
    Now we can annotated all genes in thee contig database, using NCBI COG and KEGG koFAMs.
    ```bash
    # NCBI COGs
    anvi-setup-ncbi-cogs
    anvi-run-ncbi-cogs -c CONTIGS.db

    # KEGG koFAM
    anvi-setup-kegg-data
    anvi-run-kegg-kofams -c CONTIGS.db
    ```

    ### Make profile database
    
    Now we need to make profiles for each sample to analyse MAG dynamics across samples.

    We need to ensure data quality of the downloaded data. Here we use fastp and bbmap. 

    ```bash
    mapfile -t FILES < SAMPLES.txt
    for file in "${FILES[@]}"
        do
        doQC.sh $file
        doRename.sh $file
    done
    ```

    Now we can map our sample files to the MAG catalogue. 

    ```bash    
    mapfile -t FILES < SAMPLES.txt
    for sample in "${FILES[@]}"
        do    
        # Set the input file names based on the user input
        read1="renamed_filt_$sample_1.fastq.gz"
        read2="renamed_filt_$sample_2.fastq.gz"

        # Print the starting message
        echo "Starting genome-resolved metagenomics on: $sample"

        # Index contigs for profiling
        # Align reads to contigs using bwa
        echo "Alligning reads to contigs using BWA"
        mkdir $sample
        bwa mem -t 2 contigs-fixed.fa "$read1" "$read2" > $sample/$sample_aln_pe.sam

        # Convert SAM to BAM
        echo "convert SAM to BAM"
        samtools view -@ 2 -bS "$sample/$sample_aln_pe.sam" > "$sample/$sample_aln_pe.bam"

        # Initialize BAM file for anvi'o
        # Profile the BAM file using anvi'o

        # Finally, because single profiles are rarely used for genome binning or visualization,
        # and since the clustering step increases the profiling runtime for no good reason,
        # the default behavior of profiling is to not cluster contigs automatically.
        # However, if you are planning to work with single profiles,
        # and if you would like to visualize them using the interactive interface without any merging,
        # you can use the --cluster-contigs flag to initiate clustering of contigs.
        echo "start profiling"
        anvi-init-bam "$sample/$sample_aln_pe.bam" -o "$sample/$sample.bam"
        anvi-profile -i "$sample/$sample.bam" -c CONTIGS.db" -o "PROFILES/$sample/" -T 2 --cluster-contigs

        # Print the completion message
        echo "Done with: $sample"
    ```

    Now all single profiles has been made and we only need to merge them. :neckbeard:

3. **Import MAGs from collections:** Use Anvi'o's binning tools to group contigs into metagenome-assembled genomes (MAGs).
4. **Summarize Data:** Explore and visualize the genomic and metabolic information using Anvi'o's interactive interface.
5. **Analyze and Interpret Results:** Analyze the MAGs to gain insights into the microbial community structure and function.

## Advanced Analysis
Anvi'o offers advanced features for in-depth metagenomic analysis, including:
- **Functional Annotation:** Annotate genes within MAGs to understand their functions.
- **Comparative Genomics:** Compare multiple MAGs to identify similarities and differences.
- **Metabolic Pathway Analysis:** Explore metabolic pathways encoded within the MAGs.

__Happy analyzing with Anvi'o!__


```bash
anvi-gen-genomes-storage -e contigs-db -o genomes-storage.db
anvi-map -t genomes-storage.db -p PROFILE.db -c CONTIGS.db -o BAM_FILE.bam
```