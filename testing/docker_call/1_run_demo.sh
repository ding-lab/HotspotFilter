IMAGE="mwyczalkowski/hotspot_filter:20200428"

VCF_A="/data/VCF_A.varscan_snv_vcf.vcf"
VCF_B="/data/VCF_B.varscan_indel_vcf.vcf"
BED="/data/test.bed"
OUT="/data/test_output/HotspotFiltered.vcf"

# This is what we want to run in docker
CMD_INNER="/bin/bash /opt/HotspotFilter/src/hotspot_filter.sh -A $VCF_A -B $VCF_B -D $BED -o $OUT"


SYSTEM=docker   # docker MGI or compute1
START_DOCKERD="~/Projects/WUDocker"  # https://github.com/ding-lab/WUDocker.git

VOLUME_MAPPING="../demo_data:/data"

>&2 echo Launching $IMAGE on $SYSTEM
CMD_OUTER="bash $START_DOCKERD/start_docker.sh -I $IMAGE -M $SYSTEM -c \"$CMD_INNER\" $@ $VOLUME_MAPPING "
echo Running: $CMD_OUTER
eval $CMD_OUTER
