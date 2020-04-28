# Launch docker environment for testing mutect-tool

SYSTEM=docker   # MGI and compute1
IMAGE="mwyczalkowski/hotspot_filter:20200428"
START_DOCKERD="~/Projects/WUDocker"  # https://github.com/ding-lab/WUDocker.git

VOLUME_MAPPING="../demo_data:/data"

# Also need: /storage1/fs1/dinglab/Active/CPTAC3/Common/CPTAC3.catalog
>&2 echo Launching $IMAGE on $SYSTEM
CMD="bash $START_DOCKERD/start_docker.sh -I $IMAGE -M $SYSTEM $@ $VOLUME_MAPPING"
echo Running: $CMD
eval $CMD

