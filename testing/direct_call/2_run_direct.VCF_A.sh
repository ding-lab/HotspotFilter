# Example command to run within docker.  Typically, start docker first with 0_start_docker.sh

# Variant - VCF_B is not specified.
# Retain VCF_A variants within BED

VCF_A="/data/VCF_A.varscan_snv_vcf.vcf"
BED="/data/test.bed"
OUTD="./out"
mkdir -p $OUTD
OUT="$OUTD/Hotspot2.vcf"

/bin/bash /opt/HotspotFilter/src/hotspot_filter.sh -A $VCF_A -D $BED -o $OUT

echo Written to $OUT
