IMAGE="mwyczalkowski/hotspot_filter:20200429"

cd ..
docker build -t $IMAGE -f docker/Dockerfile .
