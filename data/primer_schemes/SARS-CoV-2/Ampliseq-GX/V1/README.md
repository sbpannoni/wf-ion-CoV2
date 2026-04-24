# Ampliseq SARS-CoV-2 Insight Research Assay - GX Primer Scheme

## Required file: SARS-CoV-2.scheme.bed

This file must be obtained from your Ion Torrent system and placed in this directory.

### How to export from Ion Reporter / Torrent Suite

1. In Ion Reporter, open your SARS-CoV-2 workflow
2. Navigate to the panel manifest for "AmpliSeq SARS-CoV-2 Insight Research Panel GX"
3. Export the BED file (primer coordinates against reference MN908947.3)
ion_AmpliSeq_SARSCoV-2.2020323.Designed.bed

Alternatively, the BED can be extracted from the Torrent Suite server at:
  `/results/plugins/AmpliconCoverageAnalysis/results/<run_name>/amplicon.bed`
or from the panel directory on your Torrent Server.

### Required BED format (BED6)

The file must have exactly 6 tab-separated columns with NO header line:

```
MN908947.3	<START_0BASED>	<END>	<PRIMER_NAME>	<POOL_ID>	<STRAND>
```

- Column 5 (POOL_ID): integer 1 or 2 only — not "pool1", not "Pool_1"
- Column 6 (STRAND): + for LEFT primers, - for RIGHT primers
- Coordinates are 0-based start, 1-based end (standard BED)

### Example rows

```
MN908947.3	54	76	SARS-CoV-2_1_LEFT	1	+
MN908947.3	381	404	SARS-CoV-2_1_RIGHT	1	-
MN908947.3	320	342	SARS-CoV-2_2_LEFT	2	+
MN908947.3	704	726	SARS-CoV-2_2_RIGHT	2	-
```

### Conversion from Ion Reporter BED (if pool column is text)

```bash
# Convert "Pool_1"/"Pool_2" to integers 1/2
awk 'BEGIN{OFS="\t"} {
    gsub(/Pool_|pool_|POOL_/, "", $5)
    print
}' original.bed > SARS-CoV-2.scheme.bed
```
