# Example command to run within docker.  Typically, start docker first with 0_start_docker.sh

# Retain VCF_A variants within BED
# Retain VCF_B variants outside of BED

# Update header 
# A) get all lines starting with "##" common to both, write out
# B) get all lines starting with "##" which are unique to A; append `_A` to ID field, write out
# C) if VCF_B provided, get all lines starting with "##" unique to B, append `_B`, write out
# D) write out header for filter
    ##FILTER=<ID=hotspot,Description="Retaining calls where A intersects with BED and B does not intersect with BED.  A=..., B=..., BED=...">
#   or this, if VCF_B not provided
    ##FILTER=<ID=hotspot,Description="Retaining calls where A intersect BED.  A=..., BED=..">
# E) write out CHROM line
# Write out VCF entries using bedtools:
# * Sort output of:
#   A) bedtools intersect VCF_A BED 
#   B) if VCF_B provided, bedtools subtract VCF_B BED

VCF_A="/data/VCF_A.varscan_snv_vcf.vcf"
VCF_B="/data/VCF_B.varscan_indel_vcf.vcf"

BED="/data/test.bed"
OUTD="./out"
mkdir -p $OUTD
OUT="$OUTD/Hotspot1.vcf"

/bin/bash /opt/HotspotFilter/src/hotspot_filter.sh -A $VCF_A -B $VCF_B -D $BED -o $OUT

echo Written to $OUT
