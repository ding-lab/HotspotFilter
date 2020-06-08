#!/bin/bash

# Matthew Wyczalkowski <m.wyczalkowski@wustl.edu>
# https://dinglab.wustl.edu/

read -r -d '' USAGE <<'EOF'
Retain calls from VCF_A which intersect with BED; optionally, retain calls from
VCF_B which do not intersect with BED

Usage:
  my_script.sh [options] ARGS

Options:
-h: Print this help message
-A VCF_A: required
-B VCF_B: optional
-D BED: required
-o OUTFN: required.  Output file.  
-R: Retain only variants wit PASS or . in VCF FILTER field

Header of VCF_A and VCF_B is merged by retaining all common fields, appending _A, _B to ID field of
FILTER lines which come from VCF_A and VCF_B, respectively.

If VCF_A or VCF_B end in .gz, they are uncompressed with /bin/gunzip before proceeding

EOF

# A) get all lines starting with "##" common to both, write out
# B) get all lines starting with "##" which are unique to A; append `_A` to ID field, write out
# C) get all lines starting with "##" unique to B, append `_B`, write out
# D) write out header for filter
# E) write out CHROM line
# F) write out sorted VCF, optionally excluding variants based on value of FILTER field

# specify path explicitly in CWL environment
BEDTOOLS="/opt/conda/bin/bedtools"
GUNZIP="/bin/gunzip"

function confirm {
    FN=$1
    WARN=$2
    if [ ! -s $FN ]; then
        if [ -z $WARN ]; then
            >&2 echo ERROR: $FN does not exist or is empty
            exit 1
        else
            >&2 echo WARNING: $FN does not exist or is empty.  Continuing
        fi
    fi
}

function test_exit_status {
    # Evaluate return value for chain of pipes; see https://stackoverflow.com/questions/90418/exit-shell-script-based-on-process-exit-code
    rcs=${PIPESTATUS[*]};
    for rc in ${rcs}; do
        if [[ $rc != 0 ]]; then
            >&2 echo Fatal error.  Exiting.
            exit $rc;
        fi;
    done
}

# http://wiki.bash-hackers.org/howto/getopts_tutorial
while getopts ":ho:A:B:D:R" opt; do
  case $opt in
    h)
      echo "$USAGE"
      exit 0
      ;;
    o) 
      OUTFN=$OPTARG
      ;;
    A) 
      VCF_A=$OPTARG
      confirm $VCF_A
      ;;
    B) 
      VCF_B=$OPTARG
      confirm $VCF_B
      ;;
    D) 
      BED=$OPTARG
      confirm $BED
      ;;
    R) 
      ONLY_PASS=1
      ;;
    \?)
      >&2 echo "Invalid option: -$OPTARG" 
      >&2 echo "$USAGE"
      exit 1
      ;;
    :)
      >&2 echo "Option -$OPTARG requires an argument." 
      >&2 echo "$USAGE"
      exit 1
      ;;
  esac
done
shift $((OPTIND-1))


if [ "$#" -ne 0 ]; then
    >&2 echo ERROR: Wrong number of arguments
    >&2 echo "$USAGE"
    exit 1
fi

if [ -z $VCF_A ]; then
    >&2 echo ERROR: VCF_A not specifid
    >&2 echo "$USAGE" 
    exit 1
fi
if [ -z $BED ]; then
    >&2 echo ERROR: BED not specifid
    >&2 echo "$USAGE" 
    exit 1
fi
if [ -z $OUTFN ]; then
    >&2 echo ERROR: OUTFN not specifid
    >&2 echo "$USAGE" 
    exit 1
fi

# Make the otuput directory and check for success
OUTD=$(dirname $OUTFN)
CMD="mkdir -p $OUTD"
>&2 echo Running: $CMD
eval $CMD
test_exit_status

# test if VCF_A ends in .gz and uncompress if it does.  Will write intermediate to OUTD
if [ ${VCF_A##*.} == "gz" ]; then
	VCF_A_NEW="$OUTD/VCF_A-uncompressed.vcf"
    CMD="$GUNZIP -c $VCF_A > $VCF_A_NEW"
    >&2 echo NOTE: $VCF_A is compressed.  
	>&2 echo "  " Writing to $VCF_A_NEW
	>&2 echo "  " Running: $CMD
    eval $CMD
    test_exit_status
    VCF_A=$VCF_A_NEW
fi
if [ ! -z $VCF_B ] && [ ${VCF_B##*.} == "gz" ]; then
	VCF_B_NEW="$OUTD/VCF_B-uncompressed.vcf"
    CMD="$GUNZIP -c $VCF_B > $VCF_B_NEW"
    >&2 echo NOTE: $VCF_B is compressed.  
	>&2 echo "  "Writing to $VCF_B_NEW
	>&2 echo "  "Running: $CMD
    eval $CMD
    test_exit_status
    VCF_B=$VCF_B_NEW
fi


# if VCF_B is not defined, then header is obtained just from $VCF_A
# Otherwise, write out common lines, then per-VCF lines

if [ ! -z $VCF_B ]; then

    >&2 echo Processing VCF_A = $VCF_A
    >&2 echo VCF_B = $VCF_B
    >&2 echo BED = $BED
    # First, write out common header lines
    # from : https://stackoverflow.com/questions/373810/unix-command-to-find-lines-common-in-two-files

    awk 'NR==FNR{arr[$0];next} $0 in arr' <( grep "^##" $VCF_A ) <( grep "^##" $VCF_B )  > $OUTFN
    test_exit_status 

    # Next, write out header lines unique to A, with _A appended to ID field
    ##FILTER=<ID=ID,Description="description">
    awk 'FNR==NR {a[$0]++; next} !a[$0]' <( grep "^##" $VCF_B ) <( grep "^##" $VCF_A ) | awk 'BEGIN{FS=",";OFS=","}{$1=$1"_A"; print}' >> $OUTFN
    test_exit_status 

    # Then, write out header lines unique to B, with _B appended to ID field
    awk 'FNR==NR {a[$0]++; next} !a[$0]' <( grep "^##" $VCF_A ) <( grep "^##" $VCF_B ) | awk 'BEGIN{FS=",";OFS=","}{$1=$1"_B"; print}' >> $OUTFN
    test_exit_status 

    # write out header for filter
    printf "##INFO=<ID=HOTSPOT,Number=1,Type=Character,Description=\"Hotspot filter source\">\n" >> $OUTFN
    printf "##FILTER=<ID=hotspot,Description=\"Retaining calls where A intersects with BED and B does not intersect with BED.  A=%s, B=%s, BED=%s\">\n" "$VCF_A" "$VCF_B" "$BED \n" >> $OUTFN
    test_exit_status 

    # Finally, write out CHROM header line
    grep "^#CHROM" $VCF_A >> $OUTFN
    test_exit_status 

    # now write out merged sorted VCF
	# the line below works, but fails silently when BEDTOOLS call fails (e.g., bad BED file)
	# Instead, we'll write to temp files
    # cat <($BEDTOOLS intersect -a $VCF_A -b $BED ) <($BEDTOOLS subtract -a $VCF_B -b $BED ) | $BEDTOOLS sort -i - >> $OUTFN
    # We retain intermediate files

    # adding INFO field here as ";HOTSPOT=A" or ";HOTSPOT=B"
    # INFO field is column 8 ; https://cseweb.ucsd.edu/classes/sp16/cse182-a/notes/VCFv4.2.pdf

	TMP_A="$OUTD/VCF_A.BED.tmp"

#   awk from: https://unix.stackexchange.com/questions/148114/how-to-add-words-to-an-existing-column
    CMD="$BEDTOOLS intersect -a $VCF_A -b $BED | awk 'BEGIN{FS=\"\\t\"; OFS=\"\\t\"}{\$8 = \$8 \";HOTSPOT=A\"; print}' > $TMP_A"
	>&2 echo Running $CMD
	eval $CMD
    test_exit_status 

	TMP_B="$OUTD/VCF_B.BED.tmp"
	CMD="$BEDTOOLS subtract -a $VCF_B -b $BED | awk 'BEGIN{FS=\"\\t\"; OFS=\"\\t\"}{\$8 = \$8 \";HOTSPOT=B\"; print}' > $TMP_B"
	>&2 echo Running $CMD
	eval $CMD
    test_exit_status 

	CMD="cat $TMP_A $TMP_B | $BEDTOOLS sort -i - "
    if [ $ONLY_PASS ]; then
    # command to filter out anything but PASS or . in FILTER column
        CMD="$CMD | awk 'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if (\$7 == \"PASS\" || \$7 == \".\") print}' "
    fi

	CMD="$CMD >> $OUTFN"
	>&2 echo Running $CMD
	eval $CMD
    test_exit_status 

else
    >&2 echo Processing VCF_A = $VCF_A
    >&2 echo BED = $BED
    >&2 echo VCF_B not provided

    # VCF_B is not specified.  Just write header from VCF_A, add out filter line, and VCF_A which intersects with BED
    grep "^##" $VCF_A > $OUTFN
    test_exit_status 

    printf "##FILTER=<ID=hotspot,Description=\"Retaining calls where A intersect BED.  A=%s, BED=%s\">\n" "$VCF_A" "$BED" >> $OUTFN
    test_exit_status 
    grep "^#CHROM" $VCF_A >> $OUTFN
    test_exit_status 

    CMD="$BEDTOOLS intersect -a $VCF_A -b $BED"
    if [ $ONLY_PASS ]; then
    # command to filter out anything but PASS or . in FILTER column
        CMD="$CMD | awk 'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if (\$7 == \"PASS\" || \$7 == \".\") print}' "
    fi

    CMD="$CMD >> $OUTFN"
	>&2 echo Running $CMD
	eval $CMD
    test_exit_status 
fi

>&2 echo Written to $OUTFN
     
