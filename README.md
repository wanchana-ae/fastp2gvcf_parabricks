# fastp2gvcf_parabricks
## System Requirements
- **GPU:** NVIDIA GeForce RTX 5070 Ti  
- **CPU:** AMD Ryzen 9 7950X 16-Core Processor  
- **Memory (RAM):** 32 GB

## Description
Pipeline to Call SNP/Indels (GVCF file) by Parabricks

## Tutorials
[NVIDIA Parabricks Tutorials](https://docs.nvidia.com/clara/parabricks/latest/tutorials.html)

## Sample data
```bash
wget -O parabricks_sample.tar.gz "https://s3.amazonaws.com/parabricks.sample/parabricks_sample.tar.gz" 
```

## Metadata format for multi-samples
The metadata file for multiple samples should have the following format:

Sequence_ID,Sample_ID,/PATH_to_data/sample_1.fq.gz,/PATH_to_data/sample_2.fq.gz

- `Sequence_ID` → The ID of the sequencing run  
- `Sample_ID` → The name of the sample  
- `/PATH_to_data/sample_1.fq.gz` → Path to paired-end read 1  
- `/PATH_to_data/sample_2.fq.gz` → Path to paired-end read 2

## Environment Setup

Before running the pipeline, make sure to load the required modules in your `run_parabricks.sh` script.

### Edit the script at the section:

```bash
# ==========================
# Environment setup
# ==========================
ml load fastp 
ml load clara-parabricks
ml load bwa
```

### Working Directories
Set your working directory and output paths in run_parabricks.sh:
```bash
# Set working dir
workdir="$HOME/Downloads/parabricks_sample/"

# Base directory paths
# Output dir fastp output file
OUT_TRIMMED="${workdir}Call_SNP/Trimed/"
# Output dir Mapping output file
OUT_MAPPED="${workdir}Call_SNP/Mapped/"
# Output dir GVCF output file
OUT_GVCF="${workdir}Call_SNP/GVCF/"
# Output dir log output file
OUT_LOG="${workdir}Call_SNP/log/"
# Path to reference genome
REF="${workdir}Ref/Homo_sapiens_assembly38.fasta"
```
### Sample Metadata File
The pipeline requires a metadata file to define your samples. Edit the following section in run_parabricks.sh
```bash
# Sample metadata file
# Columns: seq_id,sample_id,forward_read_path,reverse_read_path
# Example usage:
# SEQ001,Sample_A,/data/raw/SEQ001_R1.fastq.gz,/data/raw/SEQ001_R2.fastq.gz
# SEQ002,Sample_B,/data/raw/SEQ002_R1.fastq.gz,/data/raw/SEQ002_R2.fastq.gz
INFO_FILE="metadata_fix.txt"
```
### Main Loop Overview
The main loop in `run_parabricks.sh` processes each sample in the metadata file sequentially. For each sample, the pipeline performs the following steps:

```bash
while IFS=, read -r seq_id sample_id raw1 raw2; do
    # Step 0: check ref index
    # Ensures the reference genome is indexed for BWA. Creates the index if missing.
    check_bwa_index "$REF"

    # Step 1: QC + trimming
    # Uses fastp to perform quality control and trim adapters/low-quality bases. Output is saved in $OUT_TRIMMED.
    run_fastp "$sample_id" "$raw1" "$raw2"

    # Step 2: Mapping
    # Maps the trimmed reads to the reference genome using Parabricks/BWA. Output BAM files are saved in $OUT_MAPPED.
    # The fq2bam tool requires at least 38 GB of GPU memory by default; 
    # --low-memory option will reduce this to 16 GB of GPU memory at the cost of slower processing.
    run_fq2bam "$sample_id" "$in1" "$in2"
    
    # Step 3: Run HaplotypeCaller
    # Calls variants using HaplotypeCaller and generates GVCF files in $OUT_GVCF.
    # --htvc-low-memory option will reduce this to 16 GB of GPU memory
    run_haplotypecaller "$sample_id"

    echo "[DONE] Pipeline completed for $sample_id"
done < "$INFO_FILE"

