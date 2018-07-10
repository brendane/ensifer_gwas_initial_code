#!/bin/bash
#
# Cluster denovo assembly genes - including the old genome assemblies -
# using CD-Hit.
#
# INPUT
#
# 188 de novo assemblies + USDA1106 PacBio
#
# SETTINGS
#
QUEUE="small"
WALLTIME="00:45:00"
#
MINID=0.90          # Minimum sequence identity
WORD_SIZE=10        # Word size for clustering - set high to make it run fast
GLOBAL_ALN=1        # Use global alignment to compute ID %
MIN_ALN_PROP=0.70   # alignment must cover >= 70% of both seqs; meant to prevent very long seqs from grabbing short seqs
NTHREADS=9
MAX_LENGTH=5000  # Maximum gene length clustered (some are very long)
PPN="nodes=1:ppn=${NTHREADS}"
MEM="4GB"

PROJDIR="${HOME}/project/metab_gwas"
OUTDIR="${PROJDIR}/results/ortho/cd-hit/2017-06-30"
DENOVO188="${HOME}/../shared/de_novo_assemblies/glimmer"
PACBIO_TRANS="${PROJDIR}/data/genotype/annotation/USDA1106/2016-07-18/cds.fasta"
PACBIO_GFF="${PROJDIR}/data/genotype/annotation/USDA1106/2016-07-18/cds.gff3"
SCRIPTDIR="${OUTDIR}/script_copies"
LOGDIR="${OUTDIR}/log"
WORKDIR="${OUTDIR}/working"
SCRATCHDIR="/scratch.local/${USER}/cdhit"
max_len="${PROJDIR}/script/bin/max_length_fasta.py"
cdhit2vcf="${PROJDIR}/script/bin/cdhit2vcf.py"
cdhit2tsv="${PROJDIR}/script/bin/cdhit2tsv.py"
annotate="${PROJDIR}/script/bin/annotate_cdhit_from_ipscan.py"

mkdir -p $OUTDIR
mkdir -p $SCRIPTDIR
mkdir -p $LOGDIR
mkdir -p $WORKDIR


if [[ $PBS_ENVIRONMENT == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load cdhit/4.6.8
    module load htslib/1.3.1

    set -euo pipefail

    cd $OUTDIR

    mkdir -p $SCRATCHDIR
    cd $SCRATCHDIR
    rm -f ./*
    
    SECONDS=0

    cp $max_len .

    # Gather all nucleotide sequences
    #
    # Sequence name conversions:
    #
    #   1020a -> 1020
    #   1667  -> 1676
    #   1024a -> 1024
    #   1588  -> 1580
    #   1836  -> 1830
    #   1160  -> 1161
    rm -f "seqs.fasta"
    for strain in `ls $DENOVO188`; do
        trans="${DENOVO188}/${strain}/${strain}_transcripts.fa"
        if [ -f $trans ]; then
            sed 's/^>/>'$strain'./g' $trans | \
                sed 's/1020a/1020/g' | \
                sed 's/1024a/1024/g' | \
                sed 's/1667/1676/g' | \
                sed 's/1588/1580/g' | \
                sed 's/1836/1830/g' | \
                sed 's/1160/1161/g' | \
                ./$(basename $max_len) $MAX_LENGTH >> \
                "seqs.fasta" \
                || { echo "gathering seqs failed on ${strain}"; exit 1; }
        fi
    done
    sed 's/^>/>1106PB./g' "${PACBIO_TRANS}" | \
        ./$(basename $max_len) $MAX_LENGTH >> \
        "seqs.fasta" \
        || { echo "adding PacBio USDA1106 failed"; exit 1; }

    # Run cd-hit
    cd-hit-est -i "seqs.fasta" -o "output" -c $MINID -M 3000 -n $WORD_SIZE \
        -T $NTHREADS -G 1 -d 50 -aL $MIN_ALN_PROP -AL $MIN_ALN_PROP \
        || { echo "running CD-Hit failed"; exit 1; }

    cp $cdhit2vcf .
    ./$(basename $cdhit2vcf) "output.clstr" "output.vcf" "denovo" \
        || { echo "making VCF failed"; exit 1; }

    rm "seqs.fasta"
    rm "output"
    gzip "output.clstr" \
        || { echo "gzipping output failed"; exit 1; }

    rm -f "output.vcf.gz" "output.vcf.gz.gzi"
    bgzip -i "output.vcf" \
        || { echo "bgzipping vcf output failed"; exit 1; }

    cp $cdhit2tsv .
    ./$(basename $cdhit2tsv) "output.clstr.gz" "output.tsv" "denovo" \
        || { echo "making tsv failed"; exit 1; }

    mkdir -p "annotation"
    for strain in `ls $DENOVO188`; do
        annot="${DENOVO188}/${strain}/${strain}_predicted_proteins.ipscan.tsv"
        if [ -f $annot ]; then
            f=$(echo $annot | sed 's/1020a/1020/g' | \
                sed 's/1024a/1024/g' | \
                sed 's/1667/1676/g' | \
                sed 's/1588/1580/g' | \
                sed 's/1836/1830/g' | \
                sed 's/1160/1161/g') \
                || { echo "altering name for ${annot} failed"; exit 1; }
            cp $annot "annotation/$(basename $f)" \
                || { echo "copying file for ${strain} failed"; exit 1; }
        fi
    done

    cp $annotate .
    ./$(basename $annotate) --ref-strain "1106PB" --ref-file "$PACBIO_GFF" \
        --output "output.gff3" "output.clstr.gz" "annotation" \
        || { echo "annotating failed"; exit 1; }

    rsync -av ./ $OUTDIR/ \
        || { echo "copying files failed"; exit 1; }

    echo "RUN TIME = ${SECONDS} ($(($SECONDS/60)) minutes)"
    
else

    sfile="${SCRIPTDIR}/$(date '+%Y-%m-%d-%H%M')-$(basename $0)"
    cp $0 $sfile

    cd $WORKDIR
    qsub -A "tiffinp" -W group_list="tiffinp" \
        -q $QUEUE \
        -l "walltime=${WALLTIME},mem=${MEM},${PPN}" $sfile
    echo $WORKDIR
fi
