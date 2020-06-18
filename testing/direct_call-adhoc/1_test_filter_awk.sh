VCF="out/VCF-simple.vcf"

grep -v "^#" $VCF | awk 'BEGIN{FS="\t";OFS="\t"}{if ($7 == "." || $7 == "PASS") {$7 = "hotspot"} else {$7 = $7 ";hotspot"}; print }'
