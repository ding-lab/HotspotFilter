# Example command to run within docker.  Typically, start docker first with 0_start_docker.sh

# Variant - VCF_B is not specified.
# Retain VCF_A variants within BED

VCF_A="/cromwell-executions/TinJasmine-postcall.cwl/0a5b12e4-e40f-4d51-abcf-8b6b167638ef/call-roi_filter/inputs/731768474/filtered.vcf"
BED="/cromwell-executions/TinJasmine-postcall.cwl/0a5b12e4-e40f-4d51-abcf-8b6b167638ef/call-roi_filter/inputs/-1275953733/Homo_sapiens.GRCh38.95.allCDS.2bpFlanks.biomart.withCHR.bed"

OPT="-r"    # Retain-all

OUTD="./out"
mkdir -p $OUTD
OUT="$OUTD/Hotspot2b.vcf"
CMD="/bin/bash ../../src/hotspot_filter.sh $@ $OPT -A $VCF_A -D $BED -o $OUT"

>&2 echo Running $CMD
eval $CMD



echo Written to $OUT
