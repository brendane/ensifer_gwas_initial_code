#!/bin/bash
#
# Calculate basic popgen stats on the 153 strains using libsequence.
#
# SETTINGS
#
QUEUE="small"
WALLTIME="01:30:00"
MEM="16GB"
#
# Defaults (F84 distance method) are used for phylip
#
REPLICONS="Uni0 Uni1 Uni2"
MIN_COVERAGE=2

PROJDIR="${HOME}/project/metab_gwas"
INFILE="${PROJDIR}/results/filter_variants/2017-07-04/genome.primitives.vcf.gz"
REFERENCE="${PROJDIR}/data/genotype/annotation/USDA1106/2016-07-18/USDA1106.final.fasta"
GWASDIR="${PROJDIR}/results/gwas/phenotype/gemma/2017-07-17"
DEPTH="${PROJDIR}/../popgen/results/read_depth/153_meliloti/2016-08-24/depth.tsv"
STRAINS="${PROJDIR}/data/strain/lists/153meliloti.txt"
OUTDIR="${PROJDIR}/results/popgen/libseq/2018-03-05"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
SCRIPTDIR="${OUTDIR}/script_copies"

vcf2fasta="${PROJDIR}/script/bin/vcf2fasta2.py"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load analysis/0.8.2 ## libsequence
    module load htslib/1.3.1
    module load phylip/3.69
    module load vcflib/2016-10-12

    set -euo pipefail

    cd "$OUTDIR"

    ## Remove sneaky indels and the PAVs, convert to fasta, then run
    ## libsequence.
    cp "$vcf2fasta" .
    zcat "$INFILE" | vcfnoindels | grep -v "denovo" > \
       "snps.vcf" \
        || { echo "removing PAVs failed"; exit 1; }

    ./$(basename $vcf2fasta) --output "all.fasta" \
        --min-coverage "$MIN_COVERAGE" "snps.vcf" "$REFERENCE" "$DEPTH" \
        || { echo "making fasta file for all replicons failed"; exit 1; }
    compute -i "all.fasta" -o "all.tsv" \
        || { echo "calculating stats for all replicons failed"; exit 1; }
    for replicon in $REPLICONS; do
        ./$(basename $vcf2fasta) --output "${replicon}.fasta" \
            --replicon "$replicon" --min-coverage "$MIN_COVERAGE" \
            "snps.vcf" "$REFERENCE" "$DEPTH" \
            || { echo "making fasta file for ${replicon} failed"; exit 1; }
        compute -i "${replicon}.fasta" -o "${replicon}.tsv" \
            || { echo "calculating stats for ${replicon} failed"; exit 1; }
    done

    ## Delete intermediate files
    #rm *.fasta *.vcf

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$WORKDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y%m%d-%H%M')-$(basename $0)"
    cp "$0" "$sfile"
    cd "$WORKDIR"
    qsub2 -A "tiffinp" -W group_list="tiffinp" \
        -q "$QUEUE" -l "walltime=${WALLTIME},mem=${MEM}" $sfile
    echo "$WORKDIR"
fi
