IMAGE="mwyczalkowski/hotspot_filter:20200428"

cd ..
docker build -t $IMAGE -f docker/Dockerfile .
