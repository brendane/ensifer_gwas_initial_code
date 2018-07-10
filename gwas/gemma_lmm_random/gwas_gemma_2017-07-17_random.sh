#!/bin/bash
#
# Run mixed model GWAS with a K matrix and initially no
# covariates on 16 phenotypes chosen by Mike and Brendan + symbiosis
# phenotypes. Randomize the phenotype values 100 times and run the
# analysis on each random set. Use # denovo-based presence/absence
# variants.
#
# This version uses only one variant per LD group (at R^2 = 0.95); the
# representative for each LD group is the seed for the group.
#
# This has three steps:
#   1. prep
#   2. gwas
#   3. annotate
#
# INPUT
#
# phenotypes: process_phenotypes_2016-08-17_24hr,
#   Liana's phenotype data for the the non-biolog traits,
#   process_phenotypes_bioclim_30sec_2016-10-13,
#   process_phenotypes_composite_2017-01-31, and
#   growth_rate_20170210 for the growth rate phenotype --->
#   tabulate_phenotypes_2017-07-12
#
# genotypes: filter_variants_2017-04-26 ->
#   ld_r2_groups_2017-05-04
#
# NOTES
#
#   - Just 153 meliloti strains
#   - Variants called by aligning to USDA1106
#   - SNP/MNPs + gene RDVs
#   - standardized K-matrix (preliminary analyses showed stronger
#     effects for low MAF variants than high MAF variants - see GEMMA
#     manual)
#   - Bi- and tri-allelic variants, split into bi-allelic components
#   - GEMMA does not do logistic regression, but the manual says that
#     running a regular linear model is fine
#   - K-matrix is calculated on all individuals, not separately for
#     each phenotype
#   - MAF >= 0.05 for all steps
#   - All variant types are used to construct the K-matrix and run the
#     association; all replicons are run together
#
# I'm using GEMMA 0.94.1 instead of 0.96 because 0.96 depends on a more
# recent version of GLIBC than is installed on this server, and I failed
# at installing a newer version.
#
# SETTINGS
#
# pbs
QUEUE="small"
WALLTIME="01:00:00"
#
# data processing
MIN_MAF=0.05      # Filter applied for all individuals and for each phenotype
MAX_MISS=0.2
#
# GWAS
N_ITERS=100      # Number of random iterations
PLANT_PHENOTYPES="nodule_A17 weight_A17 nodule_R108 weight_R108"
BINARY_PHENOTYPES="Anti_Gent Anti_Spec Anti_Strept metal_Cd"
CONTINUOUS_PHENOTYPES="2_Aminoethanol D_Trehalose Formic_Acid Growth_Rate L_Fucose N_Acetyl_D_Glucosamine PC1 PC2 PEG_max Putrescine Salt_max Temp_max $PLANT_PHENOTYPES bio_1 bio_12"
K_MAT_TYPE=2   # 1 = centered, 2 = standardized

RUN="2017-07-17_random"
PROJDIR="${HOME}/project/metab_gwas"
GENODIR="${PROJDIR}/results/filter_variants/2017-07-04"
GENESGFF="${PROJDIR}/data/genotype/annotation/USDA1106/2016-07-18/cds.gff3"
PHENOFILE="${PROJDIR}/results/tabulate_phenotypes/2017-07-12/phenotypes.mel153.tsv"
LDFILE="${PROJDIR}/results/ld/r2_groups/2017-07-04/0.95/one_variant.tsv"
MEL_153="${PROJDIR}/data/strain/lists/153meliloti.txt"

OUTDIR="${PROJDIR}/results/gwas/phenotype/gemma/${RUN}"
SCRIPTDIR="${OUTDIR}/script_copies"
ARRAYJOBDIR="${OUTDIR}/arrayjobdata"
LOGDIR="${OUTDIR}/log"
WORKDIR="${OUTDIR}/working"

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

    module load bedtools/2.26.0
    module load gemma/0.94.1
    module load plink/1.90b
    module load R/3.3.1

    set -euo pipefail

    cd $OUTDIR

    datafile="${ARRAYJOBDIR}/${PBS_ARRAYID}"
    TASK=`cut -f 1 -d " " $datafile`
    PHENOTYPE=`cut -f 2 -d " " $datafile`
    DATA_TYPE=`cut -f 3 -d " " $datafile`

    case $DATA_TYPE in
        "binary") ccphenos="--1" ;; # Tell Plink to treat 0 as control instead of missing
        "continuous") ccphenos="" ;;
    esac

    case $TASK in

        "prep")

            # Make sure only 153 is present
            awk '{print $1 "\\t" $1};' $MEL_153 > "keep.txt" \
                || { echo "making keep file failed"; exit 1; }

            # Need some phenotypes to stick in the plink file
            cat "keep.txt" | \
                Rscript -e 'x=read.table(file="stdin", as.is=TRUE); write.table(cbind(x, rnorm(153)), sep="\\t", col.names=FALSE, row.names=FALSE, quote=FALSE);' > \
                ph.txt \
                || { echo "getting random phenotypes failed"; exit 1; }
            
            # Take only variant per LD group
            cp $LDFILE "one_variant.tsv" \
                || { echo "copying LD file failed"; exit 1; }
            awk '{print $4 "\\tgroup-" $5};' "one_variant.tsv" > \
                "rename.tsv" \
                || { echo "making renaming file failed"; exit 1; }
            cut -f 2 "rename.tsv" > "extract.txt" \
                || { echo "making list of variants to keep failed"; exit 1; }

            # Convert to plink format
            plink --vcf "${GENODIR}/genome.vcf.gz" --make-bed --out "plink" \
                --allow-extra-chr --allow-no-sex \
               --maf $MIN_MAF --geno $MAX_MISS \
                --extract "extract.txt" --update-map "rename.tsv" \
                --update-name --pheno "ph.txt" --keep "keep.txt" \
                || { echo "making plink file failed"; exit 1; }

            # Calculate relatedness matrix
            gemma -bfile "plink" -gk $K_MAT_TYPE -o "K" -miss $MAX_MISS \
                || { echo "calculating K matrix failed"; exit 1; }

            rm -rf "K"
            mv "output" "K"

            rm "ph.txt"

            ;;

        "gwas")

            SECONDS=0

            mkdir -p $PHENOTYPE
            cd $PHENOTYPE
            rm -rf "output"
            mkdir -p "output"

            # Create a phenotype input file for plink.
            awk -v c1="strain" -v c2=$PHENOTYPE "$awkcols2" $PHENOFILE  \
                > "phenotype.tsv" \
                || { echo "making phenotype file for ${PHENOTYPE} failed"; exit 1; }

            # Check that there are no real -9s in the phenotype data (used as
            # missing code by plink).
            awk '$3 == -9{print};' "phenotype.tsv" | grep "." && \
                { echo "overlap with missing data code: failed"; exit 1; }

            for iter in $(seq 0 $(($N_ITERS-1))); do

                # Randomize the phenotype
                Rscript -e 'x=read.table("phenotype.tsv");xx=x[!is.na(x[,3]),];rnd=cbind(xx[,1:2],xx[sample(1:nrow(xx),nrow(xx),FALSE),3]);names(rnd)=names(x);write.table(rbind(rnd,x[is.na(x[,3]),]),col.names=FALSE,row.names=FALSE,quote=FALSE);' > \
                    "phenotype.${iter}.tsv" \
                    || { echo "Random iteration ${iter} for ${PHENOTYPE} failed"; exit 1; }

                # Make a binary plink file with this phenotype included
                # Note that the --prune option is not used, because all
                # the individuals used to construct the K-matrix must be
                # present. GEMMA will drop those with missing phenotype
                # values.
                plink --bfile "../plink" --allow-extra-chr \
                    --make-bed --allow-no-sex --out "plink.${iter}" \
                    --pheno "phenotype.${iter}.tsv" $ccphenos \
                    || { echo "making plink file for iter ${iter}, ${PHENOTYPE} failed"; exit 1; }

                # Run GEMMA's simplest LMM
                # -lmm 4 = all p-value tests - necessary to get the beta estimate
                # -miss = filter out sites with > $MAX_MISS missing
                # -maf = filter out sites with MAF < $MIN_MAF
                # -k = Use the k matrix made in the first step
                gemma -bfile "plink.${iter}" -k "../K/K.sXX.txt" -lmm 4 \
                    -o "output.random.${iter}" -miss $MAX_MISS -maf $MIN_MAF \
                    || { echo "gemma failed on iter ${iter}, ${PHENOTYPE}"; exit 1; }

                rm -f "output/output.random.${iter}.assoc.txt.gz"
                gzip "output/output.random.${iter}.assoc.txt" \
                    || { echo "gzipping iter ${iter}, ${PHENOTYPE} failed"; exit 1; }

                rm "plink.${iter}"*

            done

            # See how long the script really took for future reference
            echo "ASSOCIATION RUN TIME ${PHENOTYPE} = ${SECONDS} seconds or $(($SECONDS/60)) minutes"

            ;;

        "combine")
            cd $PHENOTYPE
            rm -f "output.combined.assoc.txt"
            rm -f "output.combined.assoc.txt.gz"
            
            echo "run" > "run"
            zcat "output/output.random.0.assoc.txt.gz" | \
                awk 'NR == 1 {print "run\\t" $0};' > \
                "output.combined.assoc.txt" \
                || { echo "starting combined file for ${PHENOTYPE} failed"; exit 1; }
            rm "run"

            for iter in $(seq 0 $(($N_ITERS-1))); do
                zcat "output/output.random.${iter}.assoc.txt.gz" | \
                    tail -n +2 | \
                    awk '{print "'$iter'\\t" $0};' >> "output.combined.assoc.txt" \
                    || { echo "adding ${iter} to combined output for ${PHENOTYPE} failed"; exit 1; }
            done

            gzip "output.combined.assoc.txt" \
                || { echo "gzipping combined output for ${PHENOTYPE} failed"; exit 1; }
            ;;
    esac

    rm -f $datafile

else

    sfile="${SCRIPTDIR}/$(date '+%Y-%m-%d-%H%M')-${1}-$(basename $0)"
    cp $0 $sfile

    # Run each phenotype in a separate job (of a single array job)
    i=0
    case $1 in
        "prep")
            WALLTIME="00:05:00"
            echo "prep" "_" "_" > "${ARRAYJOBDIR}/${i}"
            i=$(($i+1))
            ;;

        "gwas")
            for phenotype in $BINARY_PHENOTYPES; do
                echo "gwas" $phenotype "binary" > \
                    "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
            for phenotype in $CONTINUOUS_PHENOTYPES; do
                echo "gwas" $phenotype "continuous" > \
                    "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
            ;;

        "combine")
            WALLTIME="00:20:00"
            for phenotype in $BINARY_PHENOTYPES; do
                echo "combine" $phenotype "binary" > \
                    "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
            for phenotype in $CONTINUOUS_PHENOTYPES; do
                echo "combine" $phenotype "continuous" > \
                    "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
            ;;
    esac

    cd $WORKDIR
    qsub -A "tiffinp" -W group_list=tiffinp \
        -q $QUEUE -l walltime=$WALLTIME -t 0-$(($i-1)) $sfile
    echo $WORKDIR
    echo "Submitted ${i} jobs"

fi
