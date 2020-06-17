Testing errors in TinJasmine run

```
Running: mkdir -p output
Processing VCF_A = /cromwell-executions/TinJasmine.cwl/739158e3-c905-4484-a2d8-dde3439c28d9/call-roi_filter/inputs/1204263469/filtered.vcf
BED = /cromwell-executions/TinJasmine.cwl/739158e3-c905-4484-a2d8-dde3439c28d9/call-roi_filter/inputs/-1541611278/Homo_sapiens.GRCh38.95.allCDS.2bpFlanks.biomart.bed
VCF_B not provided
Running /opt/conda/bin/bedtools intersect -a /cromwell-executions/TinJasmine.cwl/739158e3-c905-4484-a2d8-dde3439c28d9/call-roi_filter/inputs/1204263469/filtered.vcf -b /cromwell-executions/TinJasmine.cwl/739158e3-c905-4484-a2d8-dde3439c28d9/call-roi_filter/inputs/-1541611278/Homo_sapiens.GRCh38.95.allCDS.2bpFlanks.biomart.bed >> output/HotspotFiltered.vcf
***** WARNING: File /cromwell-executions/TinJasmine.cwl/739158e3-c905-4484-a2d8-dde3439c28d9/call-roi_filter/inputs/1204263469/filtered.vcf has inconsistent naming convention for record:
chr1    13273   .   G   C   1298.6  PASS    AC=1;AF=0.5;AN=2;DP=235;set=gatk_snv-varscan_snv    GT:AD:DP:GQ:PL  0/1:167,68:235:99:1306,0,4237

***** WARNING: File /cromwell-executions/TinJasmine.cwl/739158e3-c905-4484-a2d8-dde3439c28d9/call-roi_filter/inputs/1204263469/filtered.vcf has inconsistent naming convention for record:
chr1    13273   .   G   C   1298.6  PASS    AC=1;AF=0.5;AN=2;DP=235;set=gatk_snv-varscan_snv    GT:AD:DP:GQ:PL  0/1:167,68:235:99:1306,0,4237

Written to output/HotspotFiltered.vcf
```

Will run with docker mounting /gscmnt/gc2541/cptac3_analysis/cromwell-workdir/cromwell-executions
