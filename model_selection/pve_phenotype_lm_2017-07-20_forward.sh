#!/bin/bash
#
# Use forwards linear model selection based on R^2 to estimate PVE
# and choose the best variants out of the top LMM variants - using
# PAVs from the denovo assembly.
#
# INPUT
#
# filter_variants_2017-07-04
# gwas_gemma_2017-07-17 (_random)
#
# NOTES
#
# This script assumes the p-value we want to use is in the 13th column
# of the GWAS output files.
#
# UPDATE 26 Feb 2018: R^2 from residuals is now done correctly.
#
# UPDATE 15 May 2018: Includes raw bioclim traits now.
#
# SETTINGS
#
QUEUE="small"
WALLTIME="20:00:00"
#
PLANT_PHENOTYPES="nodule_A17 weight_A17 nodule_R108 weight_R108"
PHENOTYPES="Anti_Gent Anti_Spec Anti_Strept metal_Cd 2_Aminoethanol D_Trehalose Formic_Acid Growth_Rate L_Fucose N_Acetyl_D_Glucosamine PC1 PC2 PEG_max Putrescine Salt_max Temp_max $PLANT_PHENOTYPES bio_1 bio_12"
PHENOTYPES="bio_1 bio_12"
N_VARS_MAX=25
N_RAND_ITERS=100
MIN_MAF=0.05

RUN="2017-07-20_forward"
PROJDIR="${HOME}/project/metab_gwas"
GENOFILE="${PROJDIR}/results/filter_variants/2017-07-04/genome.variants.tsv"
LDFILE="${PROJDIR}/results/ld/r2_groups/2017-07-04/0.95/one_variant.tsv"
GWASDIR="${PROJDIR}/results/gwas/phenotype/gemma/2017-07-17"
GWASDIR_RND="${PROJDIR}/results/gwas/phenotype/gemma/2017-07-17_random"
OUTDIR="${PROJDIR}/results/pve/phenotype/lm/${RUN}"
SCRIPTDIR="${OUTDIR}/script_copies"
ARRAYJOBDIR="${OUTDIR}/arrayjobdata"
LOGDIR="${OUTDIR}/log"
WORKDIR="${OUTDIR}/working"
modsel="${PROJDIR}/script/bin/forward_r2_model_selection.r"

run_modsel() {
    # $1 GWAS summary file
    # $2 GWAS phenotype file
    # tail -n+2 followed by a head command causes an error, which
    # I try to avoid here
    awk '{print "group-" $5 "\\t" $4};' $LDFILE | \
        sort -k 1b,1 -t$'\\t' > "tmp.ld.groups.txt" \
        || return 1;
    zcat $1 | \
        tail -n +2 | \
        sort -k 13g,13 -t$'\\t' > \
        "tmp" || return 1;
    head -n $N_VARS_MAX "tmp" | \
        cut -f 2 | awk '{print $1 "\\t" NR};' | \
        sort -k 1b,1 -t$'\\t' > "tmp.groups.txt" \
        || return 1;
    join "tmp.groups.txt" "tmp.ld.groups.txt" | \
        sort -k 2n,2 -t" " | \
        cut -f 3 -d" " > "variants.txt" && \
    ../$(basename $modsel) $GENOFILE $2 "variants.txt" $N_VARS_MAX \
        $MIN_MAF "output" && \
    rm "tmp" "tmp.ld.groups.txt" "tmp.groups.txt" ||
    return 1;
}

mkdir -p $OUTDIR
mkdir -p $SCRIPTDIR
mkdir -p $LOGDIR
mkdir -p $WORKDIR
mkdir -p $ARRAYJOBDIR

if [[ $PBS_ENVIRONMENT == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load R/3.2.0_intel_mkl; module load ompi/intel

    set -euo pipefail

    cd $OUTDIR

    datafile="${ARRAYJOBDIR}/${PBS_ARRAYID}"
    PHENOTYPE=`cut -f 1 -d " " $datafile`

    mkdir -p $PHENOTYPE
    cd $PHENOTYPE
    cp $modsel .

    SECONDS=0

    # Do non-random first
    cd "${OUTDIR}/${PHENOTYPE}"
    mkdir -p "real"
    cd "real"
    run_modsel "${GWASDIR}/${PHENOTYPE}/output/output0.assoc.txt.gz" \
        "${GWASDIR}/${PHENOTYPE}/phenotype.tsv" \
        || { echo "running model selection on real ${PHENOTYPE} failed"; exit 1; }

    # Then random
    for i in $(seq 0 $(($N_RAND_ITERS-1))); do
        cd "${OUTDIR}/${PHENOTYPE}"
        mkdir -p "random-${i}"
        cd "random-${i}"
        run_modsel "${GWASDIR_RND}/${PHENOTYPE}/output/output.random.${i}.assoc.txt.gz" \
            "${GWASDIR_RND}/${PHENOTYPE}/phenotype.${i}.tsv" \
            || { echo "running model selection on random-${i} ${PHENOTYPE} failed"; exit 1; }
    done

    # Then combine
    cd "${OUTDIR}/${PHENOTYPE}"
    echo "run	$(head -n 1 real/output.stats.tsv)" > "output.stats.tsv" \
        || { echo "starting combined output for ${PHENOTYPE} failed"; exit 1; }
    for d in `find . -maxdepth 1 -type "d" | sed 's/\\.\\///g' | tail -n +2`; do
        tail -n +2 "${d}/output.stats.tsv" | \
            awk '{print "'$d'\\t" $0};' >> "output.stats.tsv" \
            || { echo "adding ${d} to combined output for ${PHENOTYPE} failed"; exit 1; }
    done

    echo "TIME TO RUN ${PHENOTYPE}: ${SECONDS}"

    rm $datafile
else

    sfile="${SCRIPTDIR}/$(date '+%Y-%m-%d-%H%M')-$(basename $0)"
    cp $0 $sfile

    i=0
    for pheno in $PHENOTYPES; do
        echo $pheno > "${ARRAYJOBDIR}/${i}"
        i=$(($i+1))
    done

    cd $WORKDIR
    qsub2 -A "tiffinp" -W group_list=tiffinp \
        -q $QUEUE -l walltime=$WALLTIME -t 0-$(($i-1)) $sfile
    echo $WORKDIR
    echo "Submitted ${i} jobs"

fi
