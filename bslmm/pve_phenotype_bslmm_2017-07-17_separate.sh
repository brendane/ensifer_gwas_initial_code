#!/bin/bash
#
# Use GEMMA BSLMM (Bayesian Sparse Linear Mixed Models) to estimate the
# PVE and effect sizes for all variants and replicons combined. Uses
# 16 phenotypes and runs 5 chains per phenotype; only includes one
# variant per LD group (R2 = 0.95, 2017-05-04 LD run). This version is the
# same as the 2017-07-17 run except it runs SNP/MNPs and PAVs separately.
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
#   ld_r2_groups_2017-08-08
#
#
# Uses the 0.94.1 version of GEMMA b/c 0.96 has proven challenging to
# install on this server.
#
# UPDATE 15 May 2018: added bio_1 and bio_12
#
# SETTINGS
#
WALLTIME="12:00:00"
QUEUE="small"
#
R2=0.95                             # R2 threshold for variant grouping
MIN_MAF=0.05                        # Min. MAF for variants
MAX_MISS=0.2                        # Max. missingness for variants
HWE=0                               # Don't filter on HWE
MODEL=1                             # Linear model using MCMC to fit parameters
BURNIN=2500000                      # Burnin iterations
BINARY_CHAIN_LENGTH=25000000        # Chain length for binary phenotypes
METALCD_CHAIN_LENGTH=80000000       # Chain length for Cd tolerance - needed a bit longer to converge
CONTINUOUS_CHAIN_LENGTH=6000000     # Chain length for continuous phenotypes
RECORD_STEPS=500                    # Record every $RECORD_STEPS steps
NCHAINS=5                           # Number of chains per phenotype
#
PLANT_PHENOTYPES="nodule_A17 weight_A17 nodule_R108 weight_R108"
BINARY_PHENOTYPES="Anti_Gent Anti_Spec Anti_Strept metal_Cd"
CONTINUOUS_PHENOTYPES="2_Aminoethanol D_Trehalose Formic_Acid Growth_Rate L_Fucose N_Acetyl_D_Glucosamine PC1 PC2 PEG_max Putrescine Salt_max Temp_max $PLANT_PHENOTYPES bio_1 bio_12"

RUN="2017-07-17_separate"
PROJDIR="${HOME}/project/metab_gwas"
GENODIR="${PROJDIR}/results/filter_variants/2017-07-04"
PHENOFILE="${PROJDIR}/results/tabulate_phenotypes/2017-07-12/phenotypes.mel153.tsv"
LDDIR="${PROJDIR}/results/ld/r2_groups/2017-08-08"
MEL_153="${PROJDIR}/data/strain/lists/153meliloti.txt"

OUTDIR="${PROJDIR}/results/pve/phenotype/bslmm/${RUN}"
SCRIPTDIR="${OUTDIR}/script_copies"
ARRAYJOBDIR="${OUTDIR}/arrayjobdata"
LOGDIR="${OUTDIR}/log"
WORKDIR="${OUTDIR}/working"

merger="${PROJDIR}/script/bin/merge_bslmm_runs.py"
extract_stats="${PROJDIR}/script/bin/extract_bslmm_stats.r"
awkcols2='NR==1{for(i=1;i<=NF;i++){if($i==c1){n=i}; if($i==c2){nn=i}}}; NR>1 {print $n "\\t" $n "\\t" $nn};'

mkdir -p $OUTDIR
mkdir -p $SCRIPTDIR
mkdir -p $LOGDIR
mkdir -p $WORKDIR
mkdir -p $ARRAYJOBDIR

if [[ $PBS_ENVIRONMENT == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load gemma/0.94.1
    module load plink/1.90b

    set -euo pipefail

    TASK=`cut -f 1 -d " " "${ARRAYJOBDIR}/${PBS_ARRAYID}"`
    VARTYPE=`cut -f 2 -d " " "${ARRAYJOBDIR}/${PBS_ARRAYID}"`
    PHENOTYPE=`cut -f 3 -d " " "${ARRAYJOBDIR}/${PBS_ARRAYID}"`
    DATA_TYPE=`cut -f 4 -d " " "${ARRAYJOBDIR}/${PBS_ARRAYID}"`
    CHAIN=`cut -f 5 -d " " "${ARRAYJOBDIR}/${PBS_ARRAYID}"`
    CHAIN_LENGTH=`cut -f 6 -d " " "${ARRAYJOBDIR}/${PBS_ARRAYID}"`

    cd $OUTDIR
    mkdir -p $VARTYPE
    cd $VARTYPE
    mkdir -p $PHENOTYPE
    cd $PHENOTYPE

    LDFILE="${LDDIR}/${VARTYPE}/${R2}/one_variant.tsv"

    case $TASK in

        "prep")

            # If using a binary phenotype, tell plink to treat 0 as control
            # instead of missing
            case $DATA_TYPE in
                "binary") ccphenos="--1" ;; # Tell Plink to treat 0 as control instead of missing
                "continuous") ccphenos="" ;;
            esac

            # Phenotype file
            awk -v c1="strain" -v c2=$PHENOTYPE "$awkcols2" \
                $PHENOFILE > "phenotype.tsv" \
                || { echo "making phenotype file for ${PHENOTYPE} failed"; exit 1; }

            # Check that there are no real -9s in the phenotype data (used as
            # missing code by plink).
            awk '$3 == -9{print};' "phenotype.tsv" | grep "." && \
                { echo "overlap with missing data code: failed"; exit 1; }

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

            # Make plink binary file
            plink --vcf "${GENODIR}/genome.vcf.gz" \
                --maf $MIN_MAF \
                --geno $MAX_MISS \
                --make-bed \
                --keep "keep.txt" \
                --pheno "phenotype.tsv" \
                --prune \
                --extract "extract.txt" \
                --update-map "rename.tsv" \
                --update-name \
                --allow-no-sex \
                --allow-extra-chr \
                $ccphenos \
                --out "input" \
                || { echo "making plink file for ${PHENOTYPE} failed"; exit 1; }

            ;;

        "run")

            SECONDS=0

            mkdir -p "output"

            # I think that sometimes the different runs will start at the
            # same time and end up with the same random seed. This piece
            # of code is meant to prevent that from happening.
            while [[ 1 ]]; do
                { mkdir "lock" && break; } || sleep 10
            done
            rm -r "lock"
            seed=$RANDOM \
                || { echo "getting random seed for ${PHENOTYPE}, ${CHAIN} failed"; exit 1; }
            echo $seed > "output/output_${CHAIN}.seed.txt" \
                || { echo "writing random seed for ${PHENOTYPE}, ${CHAIN} failed"; exit 1; }

            # Run GEMMA
            gemma -bfile "input" -bslmm $MODEL -maf $MIN_MAF \
                -miss $MAX_MISS -hwe $HWE -w $BURNIN -s $CHAIN_LENGTH \
                -rpace $RECORD_STEPS -o "output_${CHAIN}" -seed $seed | \
                tee "output/output_${CHAIN}.log" \
                || { echo "running GEMMA for ${PHENOTYPE}, ${CHAIN} failed"; exit 1; }

            # Get rid of the trailing tab in the file with all the MCMC results
            # Then gzip
            mv "output/output_${CHAIN}.hyp.txt" "output/hyp_${CHAIN}.txt"
            sed 's/\\t$//g' "output/hyp_${CHAIN}.txt" > \
                "output/output_${CHAIN}.hyp.txt" \
                || { echo "reformatting hyp file for ${PHENOTYPE}, ${CHAIN} failed"; exit 1; }
            rm "output/hyp_${CHAIN}.txt"
            rm -f "output/output_${CHAIN}.hyp.txt.gz"
            gzip "output/output_${CHAIN}.hyp.txt" \
                || { echo "compressing hyp output for ${PHENOTYPE}, ${CHAIN} failed"; exit 1; }

            # Gzip the gamma.txt file
            rm -f "output/output_${CHAIN}.gamma.txt.gz"
            gzip "output/output_${CHAIN}.gamma.txt" \
                || { echo "gzipping failed for ${PHENOTYPE}, ${CHAIN}"; exit 1; }

            # See how long the script really took for future reference
            echo "BSLMM RUN TIME ${PHENOTYPE}, ${CHAIN} = ${SECONDS} seconds or $(($SECONDS/60)) minutes"
            ;;

        "combine")
            # Do not try to gzip files in this step - it ends up being
            # a problem. Do the gzipping manually once everything
            # seems okay.
            cp $merger .
            ./$(basename $merger) "output" \
                || { echo "merging output files failed"; exit 1; }

            cp "$extract_stats" .
            ./$(basename "$extract_stats") "output/output.hyp.merged.tsv" \
                "$RECORD_STEPS" "output/output_0.log.txt" > "stats.txt" \
                || { echo "calculating HPD intervals for ${PHENOTYPE} failed"; exit 1; }
            ;;

    esac

    rm -f "${ARRAYJOBDIR}/${PBS_ARRAYID}"

else

    sfile="${SCRIPTDIR}/$(date '+%Y-%m-%d-%H%M')-${1}-$(basename $0)"
    cp $0 $sfile

    i=0
    case $1 in

        "prep")
            WALLTIME="00:05:00"        
            for vartype in "snp" "rdv"; do
                for phenotype in $BINARY_PHENOTYPES; do
                    echo $1 $vartype $phenotype "binary" "_" "_" > \
                        "${ARRAYJOBDIR}/${i}"
                    i=$(($i+1))
                done
                for phenotype in $CONTINUOUS_PHENOTYPES; do
                    echo $1 $vartype $phenotype "continuous" "_" "_" > \
                        "${ARRAYJOBDIR}/${i}"
                    i=$(($i+1))
                done
            done
            ;;

        "run")
            for vartype in "snp" "rdv"; do
                for chain in $(seq 0 $(($NCHAINS-1))); do
                    for phenotype in $BINARY_PHENOTYPES; do
                        if [[ $phenotype == "metal_Cd" ]]; then
                            cl=$METALCD_CHAIN_LENGTH
                        else
                            cl=$BINARY_CHAIN_LENGTH
                        fi
                        echo $1 $vartype $phenotype "binary" $chain $cl > \
                            "${ARRAYJOBDIR}/${i}"
                        i=$(($i+1))
                    done
                    for phenotype in $CONTINUOUS_PHENOTYPES; do
                        echo $1 $vartype $phenotype "continuous" $chain $CONTINUOUS_CHAIN_LENGTH > \
                            "${ARRAYJOBDIR}/${i}"
                        i=$(($i+1))
                    done
                done
            done
            ;;

        "combine")
            WALLTIME="00:20:00"
            for vartype in "snp" "rdv"; do
                for phenotype in $BINARY_PHENOTYPES $CONTINUOUS_PHENOTYPES; do
                    echo $1 $vartype $phenotype "_" "_" "_" > "${ARRAYJOBDIR}/${i}"
                    i=$(($i+1))
                done
            done
            ;;

        *)
            echo "${1} not recognized"
            exit 1
            ;;
    esac

    cd $WORKDIR
    qsub -A "tiffinp" -W group_list="tiffinp" \
        -q $QUEUE -l walltime=$WALLTIME -t 0-$(($i-1)) $sfile
    echo $WORKDIR
    echo "Submitted ${i} jobs"

fi
