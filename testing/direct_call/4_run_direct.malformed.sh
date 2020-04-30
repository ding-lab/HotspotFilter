# Example command to run within docker.  Typically, start docker first with 0_start_docker.sh

VCF_A="/data/VCF_A.varscan_snv_vcf.vcf"
VCF_B="/data/VCF_B.varscan_indel_vcf.vcf"

BED="/data/test-malformed.bed"
OUTD="./out"
mkdir -p $OUTD
OUT="$OUTD/Hotspot4.vcf"

/bin/bash /opt/HotspotFilter/src/hotspot_filter.sh -A $VCF_A -B $VCF_B -D $BED -o $OUT

echo Written to $OUT
