# Example command to run within docker.  Typically, start docker first with 0_start_docker.sh

VCF_A="/data/varscan_snv_vcf.vcf"
VCF_B="/data/varscan_indel_vcf.vcf"
BED="/data/test.bed"
OUTD="./out"
mkdir -p $OUTD
OUT="$OUTD/Hotspot1.vcf"

/bin/bash /opt/HotspotFilter/src/hotspot_filter.sh -A $VCF_A -B $VCF_B -D $BED -o $OUT

echo Written to $OUT
