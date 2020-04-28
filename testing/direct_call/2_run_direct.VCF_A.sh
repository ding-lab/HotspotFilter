# Example command to run within docker.  Typically, start docker first with 0_start_docker.sh

VCF_A="/data/varscan_snv_vcf.vcf"
BED="/data/test.bed"
OUTD="./out"
mkdir -p $OUTD
OUT="$OUTD/Hotspot2.vcf"

/bin/bash /opt/HotspotFilter/src/hotspot_filter.sh -A $VCF_A -D $BED -o $OUT

echo Written to $OUT
