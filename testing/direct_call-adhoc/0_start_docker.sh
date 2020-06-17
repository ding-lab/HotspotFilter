# Launch docker environment for testing mutect-tool

source ../../docker/docker_image.sh
SYSTEM=MGI   # MGI and compute1
START_DOCKERD="../../docker/WUDocker"  # https://github.com/ding-lab/WUDocker.git

VOLUME_MAPPING="/gscmnt/gc2541/cptac3_analysis/cromwell-workdir/cromwell-executions:/cromwell-executions"

# Also need: /storage1/fs1/dinglab/Active/CPTAC3/Common/CPTAC3.catalog
>&2 echo Launching $IMAGE on $SYSTEM
CMD="bash $START_DOCKERD/start_docker.sh -I $IMAGE -M $SYSTEM $@ $VOLUME_MAPPING"
echo Running: $CMD
eval $CMD

