# Using Anvi'o for Genome-Resolved Metagenomics and other 'omics

## Introduction
Anvi'o is a powerful tool for analyzing metagenomic data and performing genome-resolved metagenomics. In this guide, we will walk through the steps to utilize Anvi'o for a functional metagenomics analysis.

## Installation
To get started with Anvi'o, you need to install it on your system. Follow the installation instructions provided on the official Anvi'o website.

For more information and detailed documentation, refer to the [Anvi'o documentation](https://merenlab.org/software/anvio/).

## Data Preparation
Before running Anvi'o, make sure you have the necessary metagenomic data ready. In our case, we use an EU-based data portal for multi omics analysis of host-microbe interactions in domesticated chicken and Atlantic salmon. This portal includes raw sequencing reads, assembly files, and any other relevant data for your analysis.

For more information and detailed documentation, refer to the [HoloFood Data portal](https://www.holofooddata.org/).

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

    # Annotation of marker genes with HMMs
    anvi-run-hmms -c CONTIGS.db
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

    Now we can map our samples to the MAG catalogue, using an old school loop. :yum:

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

    Which is an easy-peasy-one-liner. 

    ```bash
    anvi-merge -c CONTIGS.db PROFILES/*/PROFILE.db -o MERGED_PROFILES
    ```

    3. **Import MAGs from collections:** 

    Now we need to make and add a collection, to include our MAG information across contigs. The MAG information is stored in the header of the contigs, but we need them as splits. Therefore, we do some `anvi` and `awk` magic. 

    ```bash
    # Get splits and coverages from all contigs in your contig databases.
    anvi-export-splits-and-coverages -c CONTIGS.db -p MERGED_PROFILES/PROFILE.db
    # get split name and old mag name in to columns, using awk
    awk '{split($1, arr, "_"); print $1 "\t" arr[1]}' MERGED_PROFILES/MERGED_PROFILES-COVs.txt > MERGED_PROFILES/COLLECTION.txt
    ```
    Now include the collection to our CONTIG.db and PROFILE.db

    ```bash
    anvi-import-collection -p MERGED_PROFILES/PROFILE.db -c CONTIGS.db MERGED_PROFILES/COLLECTION.txt -C MAGs
    ```

    4. **Summarize and export Data:** Explore and visualize the genomic and metabolic information using Anvi'o's interactive interface.

    ```bash
    # Get functional annotations
    anvi-export-functions -c CONTIGS.db -o FUNCTIONS.txt
    # get coverage for each gene
    anvi-export-gene-coverage-and-detection -c CONTIGS.db -p PROFILE.db -O FUNCTIONS
    ```

    
    5. **Analyze and Interpret Results:** Analyze the MAGs to gain insights into the microbial community structure and function.

    Now we should be able to see an anvio illustatrin of our MAGs across samples. 
    ```bash
    anvi-interactive -c CONTIGS.db -p PROFILE.db -C MAGs
    ```

    Furthermore, we can add additional sample information, as in this case, we will add fatty acid measurements, and grouping (based on diet).
    ```bash
    anvi-import-misc-data view.txt \
                              -p PROFILE.db \
                              --target-data-table layers --just-do-it
    ```
    
    
    ...and here we get all other data related to the MAGs across samples.
    ```bash
    anvi-summarize -c CONTIGS.db -p PROFILE.db -C MAGs
    ``` 

    ## Advanced Analysis

    Now we will incrporate the information gotten from anvi'o with other data levels, including the metatranscriptome. 

    The extra data can be found in the integration folder [folder](https://www.holofooddata.org/).

    Anvi'o offers advanced features for in-depth metagenomic analysis, including:
    - **Functional Annotation:** Annotate genes within MAGs to understand their functions.
    - **Comparative Genomics:** Compare multiple MAGs to identify similarities and differences.
    - **Metabolic Pathway Analysis:** Explore metabolic pathways encoded within the MAGs.

    __Happy analyzing with Anvi'o!__


    ```bash
    anvi-gen-genomes-storage -e contigs-db -o genomes-storage.db
    anvi-map -t genomes-storage.db -p PROFILE.db -c CONTIGS.db -o BAM_FILE.bam
    ```