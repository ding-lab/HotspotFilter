
Starting project to merge VCFs according to BED file

Input:
  * A_VCF
  * B_VCF - optional
  * BED

Algorithm
  * write all variants in A which intersect with BED to output
  * If specified, write all variants in B not in BED to output
  * Update header appropriately


# Background reading

## GATK CombineVariants
Merging VCFs, good discussion of GATK CombineVariants
https://gatkforums.broadinstitute.org/gatk/discussion/53/combining-variants-from-different-files-into-one

## bedtools
https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html


## Contact

Matt Wyczalkowski (m.wyczalkowski@wustl.edu)


