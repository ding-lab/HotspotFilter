# Hotspot Filter

Retains variants according to genomic regions.  

Given a BED file, merges VCF_A and VCF_B where:
* Retain VCF_A variants within BED
* Retain VCF_B variants outside of BED

Optionally, excludes variants where FILTER field is not "PASS" or "."

Care is taken to retain header information.

IMAGE="mwyczalkowski/hotspot_filter:20200608"

## Details

* Update header 
  * get all lines starting with "##" common to both, write out
  * get all lines starting with "##" which are unique to A; append `_A` to ID field, write out
  * if VCF_B provided, get all lines starting with "##" unique to B, append `_B`, write out
  * write out header for filter
```
    ##FILTER=<ID=hotspot,Description="Retaining calls where A intersects with BED and B does not intersect with BED.  A=..., B=..., BED=...">
```
   or this, if VCF_B not provided
```
    ##FILTER=<ID=hotspot,Description="Retaining calls where A intersect BED.  A=..., BED=..">
```
  * write out CHROM line

* Write out VCF entries using bedtools:
 * Sort output of:
   A) bedtools intersect VCF_A BED 
   B) if VCF_B provided, bedtools subtract VCF_B BED

## hotspot_filter.sh documentation

From `src/hotspot_filter.sh`,
```
Retain calls from VCF_A which intersect with BED; optionally, retain calls from
VCF_B which do not intersect with BED

Usage:
  my_script.sh [options] ARGS

Options:
-h: Print this help message
-A VCF_A: required
-B VCF_B: optional
-D BED: required
-o OUTFN: Output file.  If not specified, write to stdout

Header of VCF_A and VCF_B is merged by retaining all common fields, appending _A, _B to ID field of
FILTER lines which come from VCF_A and VCF_B, respectively.
```

## INFO field

If VCF_B is specified, append an INFO field identifying whether this is from set A or B as a field,
    HOTSPOT=A / B
```
##INFO=<ID=HOTSPOT,Number=1,Type=Character,Description="Hotspot filter source">
```

Note that multiple hotspot filters may be used, this will just append their details

# Extra documentation details below

Starting project to merge VCFs according to BED file

Input:
  * A_VCF
  * B_VCF - optional
  * BED

Algorithm
  * write all variants in A which intersect with BED to output
  * If specified, write all variants in B not in BED to output
  * Update header appropriately


For the VCF part:
    cat <(bedtools intersect -a varscan_snv_vcf.vcf -b test.bed) <(bedtools subtract -a varscan_snv_vcf.vcf -b test.bed) | bedtools sort -i -

# Header

Example header:
```
##FILTER=<ID=read_depth,Description="Retain calls where read depth in tumor > 14 and normal > 8. Caller = varscan ">
```

## Merging headers

Generic filter field:
```
##FILTER=<ID=ID,Description="description">
```

I could add a suffix A and B to ID of every FILTER field which differs

This will work to add "-A" to ID field:

awk 'BEGIN{FS=",";OFS=","}{$1=$1"-A"; print}' filter.txt

## common files
    from : https://stackoverflow.com/questions/373810/unix-command-to-find-lines-common-in-two-files
    awk 'NR==FNR{arr[$0];next} $0 in arr' file1 file2

## Unique to one or the other

Retain lines unique to A:
$ awk 'FNR==NR {a[$0]++; next} !a[$0]' B.vcf A.vcf


# Summary of merging headers

A) get all lines starting with "##" common to both, write out
B) get all lines starting with "##" which are unique to A; append `_A` to ID field, write out
C) get all lines starting with "##" unique to B, append `_B`, write out
D) write out header for filter
E) write out CHROM line

Then write out sorted VCF

# Background reading

## GATK CombineVariants
Merging VCFs, good discussion of GATK CombineVariants
https://gatkforums.broadinstitute.org/gatk/discussion/53/combining-variants-from-different-files-into-one

## bedtools
https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html


# Contact

Matt Wyczalkowski (m.wyczalkowski@wustl.edu)


