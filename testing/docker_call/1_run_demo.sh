# TODO: follow model here: /home/mwyczalk_test/Projects/TinDaisy/mutect-tool/testing/docker_call


IMAGE="mwyczalkowski/varscan_vcf_remap"
DATD="../demo_data"

# Using python to get absolute path of DATD.  On Linux `readlink -f` works, but on Mac this is not always available
# see https://stackoverflow.com/questions/1055671/how-can-i-get-the-behavior-of-gnus-readlink-f-on-a-mac
ADATD=$(python -c 'import os,sys;print(os.path.realpath(sys.argv[1]))' $DATD)

VCF="/data/varscan_snv_vcf.vcf"
OUT="/data/out/varscan_snv_vcf-remapped.vcf"

CMD="python /opt/varscan_vcf_remap/src/varscan_vcf_remap.py --input $VCF --output $OUT"

docker run -v $ADATD:/data -it $IMAGE $CMD

