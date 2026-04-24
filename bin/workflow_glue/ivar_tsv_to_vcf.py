"""Convert iVar variants TSV to bgzipped, tabix-indexed VCF."""
import argparse
import sys
from datetime import date

import pandas as pd
import pysam


def argparser():
    """Return argument parser."""
    p = argparse.ArgumentParser(
        description="Convert iVar variants TSV to VCF (bgzipped + tabix-indexed).")
    p.add_argument("-t", "--tsv", required=True, help="iVar variants TSV")
    p.add_argument("-r", "--reference", required=True,
                   help="Reference FASTA (needed for deletion REF allele lookup)")
    p.add_argument("-s", "--sample", required=True, help="Sample name for VCF SAMPLE column")
    p.add_argument("-o", "--output", required=True, help="Output VCF path (.vcf.gz)")
    return p


def _decode_alt(ref_base, alt_str, fasta, chrom, pos_1based):
    """Return (vcf_ref, vcf_alt) from iVar ALT encoding.

    iVar encodes:
      SNP:         single base  e.g. 'A', 'T'
      Insertion:   +ACGT        e.g. '+ATG'  → REF=ref_base, ALT=ref_base+ATG
      Deletion:    -N           e.g. '-3'    → need to look up N deleted bases
    """
    if alt_str.startswith("+"):
        inserted = alt_str[1:]
        return ref_base, ref_base + inserted

    if alt_str.startswith("-"):
        try:
            del_len = int(alt_str[1:])
        except ValueError:
            return ref_base, ref_base  # malformed; treat as no-op
        # fetch deleted bases (0-based: pos to pos+del_len)
        pos_0 = pos_1based - 1
        try:
            deleted_seq = fasta.fetch(chrom, pos_0, pos_0 + del_len).upper()
        except Exception:
            deleted_seq = "N" * del_len
        vcf_ref = ref_base + deleted_seq
        vcf_alt = ref_base
        return vcf_ref, vcf_alt

    # SNP
    return ref_base, alt_str


def _phred_to_prob(phred):
    """Convert Phred score to probability; clamp to [0, 1]."""
    if phred <= 0:
        return 1.0
    return min(1.0, 10 ** (-phred / 10.0))


def main(argv=None):
    """Entry point."""
    args = argparser().parse_args(argv)

    df = pd.read_csv(args.tsv, sep="\t")

    # iVar TSV columns:
    # REGION POS REF ALT REF_DP REF_RV REF_QUAL ALT_DP ALT_RV ALT_QUAL
    # ALT_FREQ TOTAL_DP PVAL PASS
    required = {"REGION", "POS", "REF", "ALT", "ALT_DP", "ALT_QUAL",
                "ALT_FREQ", "TOTAL_DP", "PASS"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"iVar TSV missing columns: {missing}")

    out_vcf = args.output
    # Strip .gz for the raw VCF path; we will bgzip+tabix at the end
    if out_vcf.endswith(".gz"):
        raw_vcf = out_vcf[:-3]
    else:
        raw_vcf = out_vcf
        out_vcf = out_vcf + ".gz"

    fasta = pysam.FastaFile(args.reference)

    header_lines = [
        "##fileformat=VCFv4.2",
        f"##fileDate={date.today().strftime('%Y%m%d')}",
        f"##source=ivar_tsv_to_vcf wf-ion-CoV2",
        f"##reference={args.reference}",
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
        '##INFO=<ID=ALT_DP,Number=A,Type=Integer,Description="Alt allele read depth">',
        '##FILTER=<ID=PASS,Description="iVar PASS filter">',
        '##FILTER=<ID=FAIL,Description="Did not pass iVar quality filters">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}",
    ]

    records = []
    for _, row in df.iterrows():
        chrom = str(row["REGION"])
        pos = int(row["POS"])
        ref_base = str(row["REF"]).upper()
        alt_raw = str(row["ALT"]).upper()

        vcf_ref, vcf_alt = _decode_alt(ref_base, alt_raw, fasta, chrom, pos)

        qual = int(row.get("ALT_QUAL", 0))
        filt = "PASS" if str(row.get("PASS", "FALSE")).upper() == "TRUE" else "FAIL"
        dp = int(row.get("TOTAL_DP", 0))
        af = float(row.get("ALT_FREQ", 0.0))
        alt_dp = int(row.get("ALT_DP", 0))

        info = f"DP={dp};AF={af:.4f};ALT_DP={alt_dp}"
        fmt_val = f"1/1:{dp}:{af:.4f}"

        records.append((chrom, pos, ".", vcf_ref, vcf_alt, qual, filt, info, "GT:DP:AF", fmt_val))

    # Sort by chrom then pos
    records.sort(key=lambda r: (r[0], r[1]))

    with open(raw_vcf, "w") as fh:
        for line in header_lines:
            fh.write(line + "\n")
        for rec in records:
            fh.write("\t".join(str(x) for x in rec) + "\n")

    fasta.close()

    # bgzip and tabix-index
    pysam.tabix_compress(raw_vcf, out_vcf, force=True)
    pysam.tabix_index(out_vcf, preset="vcf", force=True)

    import os
    os.remove(raw_vcf)

    print(f"Written: {out_vcf} (+ .tbi)")
