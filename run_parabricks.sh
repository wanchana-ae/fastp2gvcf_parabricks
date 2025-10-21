#!/usr/bin/env bash
set -euo pipefail

# ==========================
# Environment setup
# ==========================
ml load fastp 
ml load clara-parabricks
ml load bwa

#set working dir
workdir="$HOME/Downloads/parabricks_sample/"

# Base directory paths
OUT_TRIMMED="${workdir}Call_SNP/Trimed/"
OUT_MAPPED="${workdir}Call_SNP/Mapped/"
OUT_GVCF="${workdir}Call_SNP/GVCF/"
OUT_LOG="${workdir}Call_SNP/log/"
REF="${workdir}Ref/Homo_sapiens_assembly38.fasta"

# Sample metadata file
# Columns: seq_id,sample_id,forward_read_path,reverse_read_path
# Example usage:
#SEQ001,Sample_A,/data/raw/SEQ001_R1.fastq.gz,/data/raw/SEQ001_R2.fastq.gz
#SEQ002,Sample_B,/data/raw/SEQ002_R1.fastq.gz,/data/raw/SEQ002_R2.fastq.gz
INFO_FILE="metadata_fix.txt"

# Create directories if they don't exist
mkdir -p "${OUT_MAPPED}" "${OUT_GVCF}" "${OUT_TRIMMED}" "${OUT_LOG}"

# ==========================
# Functions
# ==========================
check_bwa_index() {
    local ref=$1
    local missing=false
    local exts=("amb" "ann" "bwt" "pac" "sa")

    echo "[INFO] Checking BWA index for $ref ..."
    for ext in "${exts[@]}"; do
        if [[ ! -s "${ref}.${ext}" ]]; then
            echo "[WARN] Missing ${ref}.${ext}"
            missing=true
        fi
    done

    if [[ "$missing" == true ]]; then
        echo "[INFO] Building BWA index for $ref ..."
        bwa index "$ref"
    else
        echo "[INFO] BWA index already exists."
    fi
}

run_fastp() {
    local sample_id=$1
    local raw1=$2
    local raw2=$3
    local out1="${OUT_TRIMMED}${sample_id}_forward_paired.fq.gz"
    local out2="${OUT_TRIMMED}${sample_id}_reverse_paired.fq.gz"
    local html="${OUT_LOG}${sample_id}.fastp.html"
    local json="${OUT_LOG}${sample_id}.fastp.json"
    local log_file="${OUT_LOG}${sample_id}.fastp.log"

    if [[ ! -s "$raw1" || ! -s "$raw2" ]]; then
        echo "[ERROR] Missing raw FASTQ for $sample_id" | tee -a "$log_file"
        return 1
    fi

    echo "[INFO] Running fastp for $sample_id ..."
    fastp \
        -i "$raw1" -I "$raw2" \
        -o "$out1" -O "$out2" \
        --thread 8 \
        --detect_adapter_for_pe \
        --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
        --length_required 16 \
        -h "$html"\
        -j "$json" 2>&1 | tee -a "$log_file"
}

run_fq2bam() {
    local sample_id=$1
    local in1=$2
    local in2=$3
    local bam_out="${OUT_MAPPED}${sample_id}.bam"
    local log_file="${OUT_LOG}${sample_id}.fq2bam.log"

    if [[ ! -s "$in1" || ! -s "$in2" ]]; then
        echo "[ERROR] Missing input FASTQ for $sample_id" | tee -a "$log_file"
        return 1
    fi

    echo "[INFO] Running fq2bam for $sample_id ..."
    #--low-memory for GPU low memory Example 16 GB
    pbrun fq2bam \
        --ref "$REF" \
        --in-fq "$in1" "$in2" \
        --read-group-sm ${sample_id} \
        --read-group-lb ${sample_id} \
        --read-group-pl ILLUMINA \
        --read-group-id-prefix ${sample_id}
        --low-memory \
        --out-bam "$bam_out" 2>&1 | tee -a "$log_file"
}

run_haplotypecaller() {
    local sample_id=$1
    local bam_in="${OUT_MAPPED}${sample_id}.bam"
    local vcf_out="${OUT_GVCF}${sample_id}.g.vcf.gz"
    local log_file="${OUT_LOG}${sample_id}.haplotypecaller.log"

    if [[ ! -s "$bam_in" ]]; then
        echo "[ERROR] BAM file missing for $sample_id" | tee -a "$log_file"
        return 1
    fi

    echo "[INFO] Running haplotypecaller for $sample_id ..."
    pbrun haplotypecaller \
        --ref "$REF" \
        --gvcf \
        --max-alternate-alleles 6 \
        --htvc-low-memory \
        --num-htvc-threads 16 \
        --in-bam "$bam_in" \
        --out-variants "$vcf_out" 2>&1 | tee -a "$log_file"
}

# ==========================
# Main loop
# ==========================
while IFS=, read -r seq_id sample_id raw1 raw2; do
    # Step 0: check ref index
    check_bwa_index "$REF"

    # Step 1: QC + trimming
    run_fastp "$sample_id" "$raw1" "$raw2"

    # Step 2: Mapping
    #in1="${OUT_TRIMMED}${sample_id}_forward_paired.fq.gz"
    #in2="${OUT_TRIMMED}${sample_id}_reverse_paired.fq.gz"
    run_fq2bam "$sample_id" "$in1" "$in2"
    
    #Step 3 : run haplotypecalle
    run_haplotypecaller "$sample_id"

    echo "[DONE] Pipeline completed for $sample_id"
done < "$INFO_FILE"