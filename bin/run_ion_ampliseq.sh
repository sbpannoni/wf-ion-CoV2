#!/bin/bash
set -euo pipefail

# run_ion_ampliseq.sh — Core Ion Torrent Ampliseq SARS-CoV-2 analysis
#
# Args:
#   1  sample_name
#   2  bam_file
#   3  scheme_name        (e.g. SARS-CoV-2)
#   4  scheme_dir         (path to primer scheme directory)
#   5  threads
#   6  ivar_min_trim_len  (min read length after primer trimming, default 30)
#   7  ivar_min_qual      (min base quality for pileup, default 20)
#   8  ivar_min_freq      (min allele frequency to call variant, default 0.25)
#   9  ivar_min_depth     (min depth to call variant/consensus, default 10)
#   10 ivar_consensus_freq (freq threshold for consensus base, default 0.0 = most common)

SAMPLE="$1"
BAM="$2"
SCHEME_NAME="$3"
SCHEME_DIR="$4"
THREADS="$5"
MIN_TRIM_LEN="${6:-30}"
MIN_QUAL="${7:-20}"
MIN_FREQ="${8:-0.25}"
MIN_DEPTH="${9:-10}"
CONSENSUS_FREQ="${10:-0.0}"

REFERENCE="${SCHEME_DIR}/${SCHEME_NAME}.reference.fasta"
PRIMERS="${SCHEME_DIR}/${SCHEME_NAME}.scheme.bed"

run_analysis() {
    # 1. Check alignment: does BAM header reference MN908947.3?
    if samtools view -H "${BAM}" | grep -q "SN:MN908947.3"; then
        echo "[ion] BAM is aligned to MN908947.3 — sorting in place"
        samtools sort -@ "${THREADS}" -o input.sorted.bam "${BAM}"
    else
        echo "[ion] BAM is unaligned or uses a different reference — extracting and re-aligning"
        samtools fastq -@ "${THREADS}" "${BAM}" \
            | minimap2 -a -x sr -t "${THREADS}" "${REFERENCE}" - \
            | samtools sort -@ "${THREADS}" -o input.sorted.bam
    fi
    samtools index input.sorted.bam

    # 2. Primer trimming with iVar
    # -e: include reads that don't match any primer (common with Ion Torrent edge amplicons)
    ivar trim \
        -i input.sorted.bam \
        -b "${PRIMERS}" \
        -p trimmed \
        -m "${MIN_TRIM_LEN}" \
        -q "${MIN_QUAL}" \
        -e
    samtools sort -@ "${THREADS}" -o primertrimmed.sorted.bam trimmed.bam
    samtools index primertrimmed.sorted.bam

    # Emit trimmed BAM under expected output names
    cp input.sorted.bam "${SAMPLE}.trimmed.rg.sorted.bam"
    cp input.sorted.bam.bai "${SAMPLE}.trimmed.rg.sorted.bam.bai"
    cp primertrimmed.sorted.bam "${SAMPLE}.primertrimmed.rg.sorted.bam"
    cp primertrimmed.sorted.bam.bai "${SAMPLE}.primertrimmed.rg.sorted.bam.bai"

    # 3. Variant calling — iVar outputs TSV; converted to VCF by ivar_tsv_to_vcf.py downstream
    # -d 600000: high depth cap to avoid pileup truncation on amplicons
    # -Q 0: pass all bases to ivar (ivar applies its own quality filter with -q)
    samtools mpileup \
        -aa -A \
        -d 600000 \
        -Q 0 \
        --reference "${REFERENCE}" \
        primertrimmed.sorted.bam \
        | ivar variants \
            -p "${SAMPLE}.ivar_variants" \
            -q "${MIN_QUAL}" \
            -t "${MIN_FREQ}" \
            -m "${MIN_DEPTH}" \
            -r "${REFERENCE}"

    # 4. Consensus generation
    # -t 0.0: use most frequent base at each position (appropriate for clonal virus)
    # -n N: mask low-coverage sites with N
    samtools mpileup \
        -aa -A \
        -d 600000 \
        -Q 0 \
        --reference "${REFERENCE}" \
        primertrimmed.sorted.bam \
        | ivar consensus \
            -p "${SAMPLE}.consensus_raw" \
            -q "${MIN_QUAL}" \
            -t "${CONSENSUS_FREQ}" \
            -m "${MIN_DEPTH}" \
            -n N

    # Rename FASTA header to sample name (downstream tools parse this)
    printf ">%s\n" "${SAMPLE}" > "${SAMPLE}.consensus.fasta"
    tail -n +2 "${SAMPLE}.consensus_raw.fa" >> "${SAMPLE}.consensus.fasta"

    # 5. Depth statistics by primer pool
    # Pool assignment by BED overlap — does not depend on read group tags
    awk 'BEGIN{OFS="\t"} $5==1' "${PRIMERS}" > pool1.bed
    awk 'BEGIN{OFS="\t"} $5==2' "${PRIMERS}" > pool2.bed

    # Write depth header
    printf "sample_name\tchrom\tpos\tdepth_fwd\tdepth_rev\tdepth\tprimer_set\n" \
        > "${SAMPLE}.depth.txt"

    for pool in 1 2; do
        samtools view -b -L "pool${pool}.bed" primertrimmed.sorted.bam \
            | samtools depth -aa -d 0 - \
            | awk -v s="${SAMPLE}" -v p="${pool}" \
                  'BEGIN{OFS="\t"} {print s, $1, $2, $3, $3, $3, p}' \
            >> "${SAMPLE}.depth.txt"
    done
}

# Run analysis; on any failure emit a clean mock output so the workflow continues
if run_analysis; then
    echo "[ion] Analysis complete for ${SAMPLE}"
else
    echo "[ion] WARNING: analysis failed for ${SAMPLE} — emitting placeholder outputs"

    # Consensus with sentinel string parsed by allConsensus process
    printf ">%s Artic-Fail\nN\n" "${SAMPLE}" > "${SAMPLE}.consensus.fasta"

    # Empty iVar TSV (header only)
    printf "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\n" \
        > "${SAMPLE}.ivar_variants.tsv"

    # Zero-depth file
    printf "sample_name\tchrom\tpos\tdepth_fwd\tdepth_rev\tdepth\tprimer_set\n" \
        > "${SAMPLE}.depth.txt"

    # Empty BAMs
    for suffix in trimmed.rg.sorted primertrimmed.rg.sorted; do
        samtools view -b -o "${SAMPLE}.${suffix}.bam" /dev/null 2>/dev/null || \
            touch "${SAMPLE}.${suffix}.bam"
        touch "${SAMPLE}.${suffix}.bam.bai"
    done
fi
