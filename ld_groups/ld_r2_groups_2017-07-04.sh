#!/bin/bash
#
# Group variants by R^2 at R2 = 0.95 after filtering on MAF and missing
# data. Also keep track of the seed for each group and use that as the
# representative. (De novo assembly based CNVs)
#
# INPUT
#
# filter_variants_2017-07-04
#
#
# SETTINGS
#
WALLTIME="06:00:00"
QUEUE="small"
#
THRESHOLDS="1.00 0.95 0.90 0.85 0.80 0.75 0.70 0.65 0.60 0.55 0.50 0.45 0.40"
THRESHOLDS="0.95"
MIN_MAF=0.05
MAX_MISSING=0.8 # Really means 20% max missingness
N_CHECKED=1000  # Number of groups to check with plink

PROJDIR="${HOME}/project/metab_gwas"
MEL_153_LIST="${PROJDIR}/data/strain/lists/153meliloti.txt"
OUTDIR="${PROJDIR}/results/ld/r2_groups/2017-07-04"
GENOTYPES="${PROJDIR}/results/filter_variants/2017-07-04/genome.vcf.gz"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
ARRAYDIR="${OUTDIR}/arrayjobdata"
SCRIPTDIR="${OUTDIR}/script_copies"
SCRATCHDIR="/scratch.global/${USER}/gwas_ld_2017-07-04"
grouper="${PROJDIR}/script/bin/r2_groups_xtra_sort"
chooser="${PROJDIR}/script/bin/ld_grouping_choose_one.py"
vcf2tsv="${PROJDIR}/script/bin/vcf2tsv.py"
checker="${PROJDIR}/script/bin/check_r2_groups_onevar.py"

mkdir -p $OUTDIR
mkdir -p $LOGDIR
mkdir -p $SCRIPTDIR
mkdir -p $WORKDIR
mkdir -p $ARRAYDIR

if [[ $PBS_ENVIRONMENT == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load plink/1.90b
    module load vcftools/0.1.15

    set -euo pipefail

    R2=`cut -f 1 -d " " "${ARRAYDIR}/${PBS_ARRAYID}"`

    cd $OUTDIR
    mkdir -p $R2
    cd $R2

    mkdir -p $SCRATCHDIR
    cd $SCRATCHDIR
    mkdir -p $R2
    cd $R2

    # First filter the VCF file down to the variants that are most useful
    # for GWAS.
    vcftools --gzvcf $GENOTYPES --stdout --recode --keep $MEL_153_LIST | \
        vcftools --vcf - --stdout --recode --max-missing $MAX_MISSING  \
        --maf $MIN_MAF > \
        "tmp.vcf" \
        || { echo "filtering vcf failed"; exit 1; }

    cp $vcf2tsv .
    ./$(basename $vcf2tsv) --het-miss --rs --output "tmp.tsv" \
        "tmp.vcf" \
        || { echo "converting to tsv format failed"; exit 1; }

    cp $grouper .
    ./$(basename $grouper) $R2 1 "tmp.tsv" > "tmp.groups.tsv" \
       || { echo "grouping at ${R2} failed"; exit 1; }

   # Sort by group and by seed status
   head -n 1 "tmp.groups.tsv" > "output.tsv"
   tail -n +2 "tmp.groups.tsv" | sort -k 6n,6 -k 7rn,7 >> "output.tsv" \
       || { echo "sorting at ${R2} failed"; exit 1; }

    # Choose one variant to represent each group; given the way
    # the file is sorted, this will be the variant with the greatest
    # number of individuals genotyped, and if there are ties, the
    # greatest minor allele frequency.
    cp $chooser .
    ./$(basename $chooser) --output "one_variant.tsv" "output.tsv" \
        || { echo "choosing one variant for ${R2} failed"; exit 1; }

    # Check the results with plink - randomly choose N_CHECKED groups
    # to test that the LD between seed and other group members is
    # >= $R2. Will return an exit code if any fail. Checking all
    # groups is very time consuming, so only a random subset is
    # checked. Also, $N_CHECKED seeds are compared against each
    # other to make sure none of them should be grouped together. This
    # is not a complete check, but if a large enough number of groups
    # are checked, it should catch major issues. It will not catch
    # cases where a variant could have been placed in more than one
    # group - this can happen due to missing data issues, and is not
    # listed as an error.
    cp $checker .
    plink --make-bed --vcf "tmp.vcf" --allow-extra-chr --out "tmp" \
        || { echo "converting to bed format failed"; exit 1; }
    ./$(basename $checker) "tmp" "one_variant.tsv" $R2 $N_CHECKED \
        || { echo "one variant file wrong for ${R2}: failed"; exit 1; }

   # Eliminate temporary files
   rm tmp.*

    cd "${OUTDIR}/${R2}"
    rsync -av "${SCRATCHDIR}/${R2}/" "./" \
        || { echo "copying files failed"; exit 1; }

    rm "${ARRAYDIR}/${PBS_ARRAYID}"

else

    sfile="${SCRIPTDIR}/$(date '+%Y-%m-%d-%H%M')-$(basename $0)"
    cp $0 $sfile

    i=0
    for R2 in $THRESHOLDS; do
        echo $R2 > "${ARRAYDIR}/${i}"
        i=$(($i+1))
    done

    cd $WORKDIR
    qsub -q $QUEUE -l walltime=$WALLTIME -t 0-$(($i-1))%100 $sfile
    echo $WORKDIR
    echo "Submitted ${i} jobs"

fi
