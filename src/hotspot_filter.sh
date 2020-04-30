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

Header of VCF_A and VCF_B is merged by retaining all common fields, appending _A, _B to ID field of
FILTER lines which come from VCF_A and VCF_B, respectively.

If VCF_A or VCF_B end in .gz, they are uncompressed with /bin/gunzip before proceeding

EOF

# A) get all lines starting with "##" common to both, write out
# B) get all lines starting with "##" which are unique to A; append `_A` to ID field, write out
# C) get all lines starting with "##" unique to B, append `_B`, write out
# D) write out header for filter
# E) write out CHROM line
# F) write out sorted VCF

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
while getopts ":ho:A:B:D:" opt; do
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
    CMD="$GUNZIP $VCF_A"
    >&2 echo VCF_A is compressed.  Running: $CMD
    eval $CMD
    test_exit_status
    VCF_A=${VCF_A%.*}
    >&2 echo New VCF_A = $VCF_A
fi
if [ ! -z $VCF_B ] && [${VCF_B##*.} == "gz" ]; then
    CMD="$GUNZIP $VCF_B"
    >&2 echo VCF_B is compressed.  Running: $CMD
    eval $CMD
    test_exit_status
    VCF_B=${VCF_B%.*}
    >&2 echo New VCF_B = $VCF_B
fi


# if VCF_B is not defined, then header is obtained just from $VCF_A
# Otherwise, write out common lines, then per-VCF lines

if [ ! -z $VCF_B ]; then

    >&2 echo Processing VCF_A = $VCF_A
    >&2 echo VCF_B = $VCF_B
    >&2 echo BED = $BED
    # First, write out common header lines
    # from : https://stackoverflow.com/questions/373810/unix-command-to-find-lines-common-in-two-files

    # Taking care to catch errors in process substitution with "|| kill $$"  See https://unix.stackexchange.com/questions/217605/bash-how-to-propagate-errors-in-process-substitution
    awk 'NR==FNR{arr[$0];next} $0 in arr' <( grep "^##" $VCF_A || kill $$ ) <( grep "^##" $VCF_B || kill $$ )  > $OUTFN
    test_exit_status 

    # Next, write out header lines unique to A, with _A appended to ID field
    ##FILTER=<ID=ID,Description="description">
    awk 'FNR==NR {a[$0]++; next} !a[$0]' <( grep "^##" $VCF_B || kill $$ ) <( grep "^##" $VCF_A || kill $$ ) | awk 'BEGIN{FS=",";OFS=","}{$1=$1"_A"; print}' >> $OUTFN
    test_exit_status 

    # Then, write out header lines unique to B, with _B appended to ID field
    awk 'FNR==NR {a[$0]++; next} !a[$0]' <( grep "^##" $VCF_A || kill $$ ) <( grep "^##" $VCF_B || kill $$ ) | awk 'BEGIN{FS=",";OFS=","}{$1=$1"_B"; print}' >> $OUTFN
    test_exit_status 

    # write out header for filter
    printf "##FILTER=<ID=hotspot,Description=\"Retaining calls where A intersects with BED and B does not intersect with BED.  A=%s, B=%s, BED=%s\">\n" "$VCF_A" "$VCF_B" "$BED" >> $OUTFN
    test_exit_status 

    # Finally, write out CHROM header line
    grep "^#CHROM" $VCF_A >> $OUTFN
    test_exit_status 

    # now write out merged sorted VCF
    cat <($BEDTOOLS intersect -a $VCF_A -b $BED || kill $$ ) <($BEDTOOLS subtract -a $VCF_B -b $BED || kill $$ ) | $BEDTOOLS sort -i - >> $OUTFN
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

    $BEDTOOLS intersect -a $VCF_A -b $BED >> $OUTFN
    test_exit_status 
fi

>&2 echo Written to $OUTFN
     
