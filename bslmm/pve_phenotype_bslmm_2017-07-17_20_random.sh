#!/bin/bash
#
# Use GEMMA BSLMM (Bayesian Sparse Linear Mixed Models) to estimate the
# PVE and effect sizes for all variants and replicons combined. Runs on
# 15 phenotypes -- the 20 phenotype group minus the five focal phenotypes.
# Runs 1 very short chain per phenotype, so the BSLMM results should not
# be trusted.  Only includes one variant per LD group (R2 = 0.95,
# 2017-07-04 LD run). Uses de novo presence/absence variants.
#
# This script is a slightly modified version of the 2017-07-17_random
# script, which ran the permuted data for the five focal traits.
#
# INPUT
#
# phenotypes: process_phenotypes_2016-08-17_24hr,
#   Liana's phenotype data for the the non-biolog traits,
#   process_phenotypes_bioclim_30sec_2016-10-13,
#   process_phenotypes_composite_2017-01-31, and
#   growth_rate_20170210 for the growth rate phenotype ->
#   tabulate_phenotypes_2017-07-12
#
# genotypes: filter_variants_2017-07-04 ->
#   ld_r2_groups_2017-07-04
#
# NOTES
#
#   153 meliloti strains only
#   No outlier strains removed
#   Biolog phenotypes processed by subtracting water and dividing by
#       mean of three sugars
#   Binary phenotypes are run with the same model as continuous traits,
#       as recommended in the manual
#   Variant types are run together
#   RDVs from de novo assembly
#   Default priors for GEMMA
#   Runs the linear model MCMC version of BSLMM
#   MAF cutoff
#
# Uses the 0.94.1 version of GEMMA b/c 0.96 has proven challenging to
# install on this server.
#
# UPDATE 21 May 2018: Runs bio_1 and bio_12
#
# SETTINGS
#
WALLTIME="01:00:00"
QUEUE="small"
#
MIN_MAF=0.05                        # Min. MAF for variants
MAX_MISS=0.2                        # Max. missingness for variants
HWE=0                               # Don't filter on HWE
MODEL=1                             # Linear model using MCMC to fit parameters
BURNIN=2500                         # Burnin iterations
BINARY_CHAIN_LENGTH=6000            # Chain length for binary phenotypes
CONTINUOUS_CHAIN_LENGTH=6000        # Chain length for continuous phenotypes
RECORD_STEPS=5                      # Record every $RECORD_STEPS steps
NCHAINS=1                           # Number of chains per phenotype
N_RANDOM_RUNS=100
#
BINARY_PHENOTYPES="Anti_Gent Anti_Strept metal_Cd"
CONTINUOUS_PHENOTYPES="2_Aminoethanol D_Trehalose Formic_Acid Growth_Rate L_Fucose N_Acetyl_D_Glucosamine PC1 PC2 Salt_max Temp_max nodule_R108 weight_R108 bio_1 bio_12"
CONTINUOUS_PHENOTYPES="bio_1 bio_12"

RUN="2017-07-17_20_random"
PROJDIR="${HOME}/project/metab_gwas"
GENODIR="${PROJDIR}/results/filter_variants/2017-07-04"
RAND_PHENOS="${PROJDIR}/results/random_phenotype_tables/2017-09-21"
LDFILE="${PROJDIR}/results/ld/r2_groups/2017-07-04/0.95/one_variant.tsv"
MEL_153="${PROJDIR}/data/strain/lists/153meliloti.txt"

OUTDIR="${PROJDIR}/results/pve/phenotype/bslmm/${RUN}"
SCRIPTDIR="${OUTDIR}/script_copies"
ARRAYJOBDIR="${OUTDIR}/arrayjobdata"
LOGDIR="${OUTDIR}/log"
WORKDIR="${OUTDIR}/working"
SCRATCHDIR="/scratch.global/${USER}/${RUN}"

extract_stats="${PROJDIR}/script/bin/extract_random_bslmm_stats.r"
awkcols2='NR==1{for(i=1;i<=NF;i++){if($i==c1){n=i}; if($i==c2){nn=i}}}; NR>1 {print $n "\\t" $n "\\t" $nn};'

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load gemma/0.94.1
    module load plink/1.90b
    module load R/3.4.3
    module load ompi/intel

    set -euo pipefail

    TASK=`cut -f 1 -d " " "${ARRAYJOBDIR}/${PBS_ARRAYID}"`
    PHENOTYPE=`cut -f 2 -d " " "${ARRAYJOBDIR}/${PBS_ARRAYID}"`
    DATA_TYPE=`cut -f 3 -d " " "${ARRAYJOBDIR}/${PBS_ARRAYID}"`
    CHAIN=`cut -f 4 -d " " "${ARRAYJOBDIR}/${PBS_ARRAYID}"`
    CHAIN_LENGTH=`cut -f 5 -d " " "${ARRAYJOBDIR}/${PBS_ARRAYID}"`

    cd $OUTDIR
    mkdir -p $PHENOTYPE
    cd $PHENOTYPE
    mkdir -p "${SCRATCHDIR}/${PHENOTYPE}"


    # If using a binary phenotype, tell plink to treat 0 as control
    # instead of missing
    case $DATA_TYPE in
        "binary") ccphenos="--1" ;; # Tell Plink to treat 0 as control instead of missing
        "continuous") ccphenos="" ;;
    esac

    case $TASK in

        "prep")

            # Take only variant per LD group
            cp $LDFILE "one_variant.tsv" \
                || { echo "copying LD file failed"; exit 1; }
            awk '{print $4 "\\tgroup-" $5};' "one_variant.tsv" > \
                "rename.tsv" \
                || { echo "making renaming file failed"; exit 1; }
            cut -f 2 "rename.tsv" > "extract.txt" \
                || { echo "making list of variants to keep failed"; exit 1; }

            # Make sure only 153 is present
            awk '{print $1 "\\t" $1};' $MEL_153 > "keep.txt" \
                || { echo "making keep file failed"; exit 1; }

            ;;

        "run")

            SECONDS=0

            mkdir -p "output"

            for RRUN in $(seq 0 $(($N_RANDOM_RUNS-1))); do
                # Phenotype file
                awk '{ print $1 "\\t" $1 "\\t" $'$(($RRUN+2))'};' \
                    "${RAND_PHENOS}/${PHENOTYPE}.tsv" \
                    > "phenotype.${RRUN}.tsv" \
                    || { echo "making pheno file for ${PHENOTYPE}, ${RRUN} failed"; exit 1; }

                # Check that there are no real -9s in the phenotype data (used as
                # missing code by plink).
                awk '$3 == -9{print};' "phenotype.${RRUN}.tsv" | grep "." && \
                    { echo "overlap with missing data code: failed"; exit 1; }

                # Make plink binary file
                plink --vcf "${GENODIR}/genome.vcf.gz" \
                    --maf $MIN_MAF \
                    --geno $MAX_MISS \
                    --make-bed \
                    --keep "keep.txt" \
                    --pheno "phenotype.${RRUN}.tsv" \
                    --prune \
                    --extract "extract.txt" \
                    --update-map "rename.tsv" \
                    --update-name \
                    --allow-no-sex \
                    --allow-extra-chr \
                    $ccphenos \
                    --out "${SCRATCHDIR}/${PHENOTYPE}/input.${RRUN}" \
                    || { echo "making plink file for ${PHENOTYPE}, ${RRUN} failed"; exit 1; }

                seed=$RANDOM \
                    || { echo "getting random seed for ${PHENOTYPE}, ${CHAIN} failed"; exit 1; }
                echo $seed > "output/output_${CHAIN}.seed.txt" \
                    || { echo "writing random seed for ${PHENOTYPE}, ${CHAIN} failed"; exit 1; }

                # Run GEMMA
                gemma -bfile "${SCRATCHDIR}/${PHENOTYPE}/input.${RRUN}" -bslmm $MODEL \
                    -maf $MIN_MAF -miss $MAX_MISS -hwe $HWE -w $BURNIN -s $CHAIN_LENGTH \
                    -rpace $RECORD_STEPS -o "output_${RRUN}" -seed $seed | \
                    tee "output/output_${RRUN}.log" \
                    || { echo "running GEMMA for ${PHENOTYPE}, ${CHAIN} failed"; exit 1; }

                # Get rid of the trailing tab in the file with all the MCMC results
                # Then gzip
                mv "output/output_${RRUN}.hyp.txt" "output/hyp_${RRUN}.txt"
                sed 's/\\t$//g' "output/hyp_${RRUN}.txt" > \
                    "output/output_${RRUN}.hyp.txt" \
                    || { echo "reformatting hyp file for ${PHENOTYPE}, ${RRUN} failed"; exit 1; }
                rm "output/hyp_${RRUN}.txt"
                rm -f "output/output_${RRUN}.hyp.txt.gz"
                gzip "output/output_${RRUN}.hyp.txt" \
                    || { echo "compressing hyp output for ${PHENOTYPE}, ${RRUN} failed"; exit 1; }

                # Delete the gamma.txt file
                rm -f "output/output_${RRUN}.gamma.txt"
            done

            # See how long the script really took for future reference
            echo "BSLMM RUN TIME ${PHENOTYPE}, ${RRUN} = ${SECONDS} seconds or $(($SECONDS/60)) minutes"
            ;;

        "summary")
            cp "$extract_stats" .
            ./$(basename "$extract_stats") "output/output_" \
                "$RECORD_STEPS" > "stats.txt" \
                || { echo "calculating HPD intervals for ${PHENOTYPE} failed"; exit 1; }
            ;;

    esac

    rm -f "${ARRAYJOBDIR}/${PBS_ARRAYID}"

else

    mkdir -p $OUTDIR
    mkdir -p $SCRIPTDIR
    mkdir -p $LOGDIR
    mkdir -p $WORKDIR
    mkdir -p $ARRAYJOBDIR

    sfile="${SCRIPTDIR}/$(date '+%Y-%m-%d-%H%M')-$(basename $0)-${1}"
    cp "$0" "$sfile"

    i=0
    case $1 in

        "prep")
            WALLTIME="00:05:00"        
#            for phenotype in $BINARY_PHENOTYPES; do
#                echo $1 $phenotype "binary" "_" "_" > \
#                    "${ARRAYJOBDIR}/${i}"
#                i=$(($i+1))
#            done
            for phenotype in $CONTINUOUS_PHENOTYPES; do
                echo $1 $phenotype "continuous" "_" "_" > \
                    "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
            ;;

        "run")
            for chain in $(seq 0 $(($NCHAINS-1))); do
#                for phenotype in $BINARY_PHENOTYPES; do
#                    echo $1 $phenotype "binary" $chain \
#                        $BINARY_CHAIN_LENGTH > \
#                        "${ARRAYJOBDIR}/${i}"
#                    i=$(($i+1))
#                done
                for phenotype in $CONTINUOUS_PHENOTYPES; do
                    echo $1 $phenotype "continuous" $chain \
                        $CONTINUOUS_CHAIN_LENGTH > \
                        "${ARRAYJOBDIR}/${i}"
                    i=$(($i+1))
                done
            done
            ;;

        "summary")
            WALLTIME="00:05:00"
#            for phenotype in $BINARY_PHENOTYPES; do
#                echo $1 $phenotype "_" "_" "_" > \
#                    "${ARRAYJOBDIR}/${i}"
#                i=$(($i+1))
#            done
            for phenotype in $CONTINUOUS_PHENOTYPES; do
                echo $1 $phenotype "_" "_" "_" > \
                    "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
            ;;

        *)
            echo "${1} not recognized"
            exit 1
            ;;
    esac

    cd $WORKDIR
    qsub2 -A "tiffinp" -W group_list="tiffinp" \
        -q $QUEUE -l walltime=$WALLTIME -t 0-$(($i-1)) $sfile
    echo $WORKDIR
    echo "Submitted ${i} jobs"

fi
