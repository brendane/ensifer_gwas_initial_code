#!/bin/bash
#
# Combine SNPs with denovo assembly-based P/A variants. Other than the
# source of the P/A variants, this script does the same thing as the
# 2017-04-26 script.
#
# NOTES
#
# Multi-allelic sites are split into multiple lines using bcftools norm.
# It looks like sites with one alternate allele are listed as reference
# for the line with the other alternate allele.
#
# The variant IDs contain all the information needed for tracking down
# the variants after running GWAS:
#   - variant type (snp or rdv)
#   - chrom, pos, gene id for RDVs
#   - alternate allele - for multi-allelic variants
#
# UPDATE 21 July 2017: fixed minor allele frequency filter to allow
#   singletons.
#
# SETTINGS
#
MAX_MISSING=0.8
#
QUEUE="small"
WALLTIME="02:30:00"

PROJDIR="${HOME}/project/metab_gwas"
RAW_SNPCALLS="${PROJDIR}/data/genotype/snps/2016-07-11/all_smeliloti_filter_qual20.vcf.gz"
RAW_RDV_CALLS="${PROJDIR}/results/ortho/cd-hit/2017-06-30/output.vcf.gz"
MEL_153_LIST="${PROJDIR}/data/strain/lists/153meliloti.txt"
OUTDIR="${PROJDIR}/results/filter_variants/2017-07-04"
SCRIPTDIR="${OUTDIR}/script_copies"
LOGDIR="${OUTDIR}/log"
WORKDIR="${OUTDIR}/working"
vcf2tsv="${PROJDIR}/script/bin/vcf2tsv.py"
vcf2fasta="${PROJDIR}/script/bin/vcf2fasta_snps.py"
realcoords="${PROJDIR}/script/bin/mnp2snp.py"

mkdir -p $OUTDIR
mkdir -p $SCRIPTDIR
mkdir -p $LOGDIR
mkdir -p $WORKDIR

if [[ $PBS_ENVIRONMENT == "PBS_BATCH" ]]; then

    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"
    newgrp "tiffinp"

    module load bcftools/1.3.1
    module load htslib/1.3.1
    module load R/3.3.1b
    module load vcflib/2016-10-12
    module load vcftools/0.1.15

    set -euo pipefail

    cd $OUTDIR

    ######
    #
    # Steps:
    #   1. Filter small variant calls:
    #       a. Keep only SNPs and MNPs
    #       b. Split multi-allelic variants into multiple lines
    #       c. Keep only sites with variation
    #       d. Keep only the 153 individuals
    #       e. Keep only sites with >= 80% genotyped
    #       f. Set IDs to something informative
    #
    #   2. Filter CNVs using the same criteria as above
    #
    #   3. Combine and sort the files; CNVs are all on a contig called
    #      "denovo"
    #
    #   4. Set heterozygous calls to missing.
    #
    #   5. Delete intermediate files
    #

    # Keep only SNPs and MNPs
    zcat $RAW_SNPCALLS | vcfnoindels > "snps.type_filtered.vcf" \
        || { echo "filtering out indels failed"; exit 1; }

    # Get rid of some annotation info that upsets bcftools
    # Normalize (split multiple alleles into multiple lines)
    # Could also use vcflib::vcfbreakmulti here, but bcftools
    # seems okay now that the indels are really gone.
    bcftools annotate -x "FMT/DPR" "snps.type_filtered.vcf" | \
        bcftools norm -m "-any" > "snps.norm.vcf" \
        || { echo "normalizing SNPs failed"; exit 1; }

    # Filter on other criteria
    vcftools --vcf "snps.norm.vcf" --stdout \
        --recode --keep $MEL_153_LIST \
        --max-missing $MAX_MISSING | \
        vcftools --vcf - --stdout  --maf "$(Rscript -e 'cat(round(1/153,4))')" \
        --recode > \
        "filtered_snps.vcf" \
        || { echo "filtering SNPs by MAF and MISS. failed"; exit 1; }

    # Set variant IDs for SNPs
    bcftools annotate --set-id 'snp-%CHROM-%POS-%FIRST_ALT' \
        "filtered_snps.vcf" > "filtered_snps.id.vcf" \
        || { echo "adding IDs to SNPs failed"; exit 1; }

    # Normalize RDVs
    # NOTE: This step doesn't do anything useful
    bcftools norm -m "-any" $RAW_RDV_CALLS > "rdvs.norm.vcf" \
        || { echo "normalizing RDVs failed"; exit 1; }
    

    # Filter RDVs
    vcftools --gzvcf "$RAW_RDV_CALLS" --stdout \
        --recode --keep $MEL_153_LIST \
        --max-missing $MAX_MISSING | \
        vcftools --vcf - --stdout  --maf "$(Rscript -e 'cat(round(1/153,4))')" \
        --recode > \
        "filtered_rdvs.vcf" \
        || { echo "filtering RDVs by MAF and MISS. failed"; exit 1; }


    # Gzip the intermediate files because vcf-concat can't handle
    # plain text.
    rm -f "filtered_snps.id.vcf.gz"
    bgzip -i "filtered_snps.id.vcf" \
        || { echo "bgzipping SNPs failed"; exit 1; }
    rm -f "filtered_rdvs.vcf.gz"
    bgzip -i "filtered_rdvs.vcf" \
        || { echo "bgzipping RDVs failed"; exit 1; }

    # Make sure samples are ordered the same
    vcf-shuffle-cols -t "filtered_snps.id.vcf.gz" \
        "filtered_rdvs.vcf.gz" > "rdvs.ordered.vcf" \
        || { echo "re-ordering RDVs failed"; exit 1; }
    rm -f "rdvs.ordered.vcf.gz"
    bgzip -i "rdvs.ordered.vcf" \
        || { echo "bgzipping ordered RDVs failed"; exit 1; }

    # Combine. Also make genotypes coded as haploid missing
    # into diploid missing and change heterozygotes to just
    # missing.
    vcf-concat "filtered_snps.id.vcf.gz" \
        "rdvs.ordered.vcf.gz" | \
        vcf-sort | \
        sed 's/\\t\\.:/\\t.\\/.:/g' | \
        sed 's/\\t0\\/1:/\\t.\\/.:/g' | \
        sed 's/\\t1\\/0:/\\t.\\/.:/g' | \
        vcftools --vcf - --stdout  --maf "$(Rscript -e 'cat(round(1/153,4))')" \
        --recode | \
        vcftools --vcf - --stdout  --max-missing $MAX_MISSING \
        --recode \
        > "genome.vcf" \
        || { echo "merging failed"; exit 1; }

    rm -f "genome.vcf.gz"
    bgzip -i "genome.vcf" \
        || { echo "bgzipping combined file failed"; exit 1; }
    tabix "genome.vcf.gz" \
        || { echo "tabixxing combined file failed"; exit 1; }

    bcftools norm -m "+any" "genome.vcf.gz" > "genome.unnorm.vcf" \
        || { echo "merge multi-allelic sites failed"; exit 1; }
    rm -f "genome.unnorm.vcf.gz"
    bgzip -i "genome.unnorm.vcf" \
        || { echo "bgzipping non-norm. combined file failed"; exit 1; }
    tabix "genome.unnorm.vcf.gz" \
        || { echo "tabixxing non-norm. combined file failed"; exit 1; }

    # Delete intermediate files
    files=$(ls . | grep "vcf" | grep -v "genome")
    rm $files

    # Create a tsv file
   cp $vcf2tsv .
   ./$(basename $vcf2tsv) --het-miss --rs --output "genome.variants.tsv" \
       "genome.vcf.gz" \
       || { echo "making tsv file failed"; exit 1; }
   ./$(basename $vcf2tsv) --het-miss --rs --output "genome.unnorm.variants.tsv" \
       "genome.unnorm.vcf.gz" \
       || { echo "making tsv file failed"; exit 1; }

    # Make a fasta file - use only sites with single bp alleles
    #
    # This program sometimes adds indels back in because it realigns
    # variant alleles.
    vcfallelicprimitives "genome.unnorm.vcf.gz" > "genome.primitives.vcf" \
        || { echo "getting allelic primitives failed"; exit 1; }
    rm -f "genome.primitives.vcf.gz"
    bgzip -i "genome.primitives.vcf" \
        || { echo " gzipping primitives failed"; exit 1; }
    cp $vcf2fasta .
    ./$(basename $vcf2fasta) --output "genome.variants.fasta" \
        --exclude "denovo" "genome.primitives.vcf.gz" \
        || { echo "making fasta file failed"; exit 1; }

    # Take the un-normalized, primitives file and convert it to a ref/alt
    # tsv format. This is the best match to the allele frequency analyses
    # that were tested with a glm (in S&R project).
    ./$(basename $vcf2tsv) --het-miss \
        --output "genome.primitives.unnorm.variants.tsv" \
        "genome.primitives.vcf.gz" \
        || { echo "making tsv file failed"; exit 1; }

    # Also added later
    # Take MNPs left in the data and get the actual variant position
    cp $realcoords .
    ./$(basename $realcoords) --output "genome.realcoords.tsv" \
        "genome.vcf.gz" \
        || { echo "getting real coordinates failed"; exit 1; }

    # Allele frequency
    vcftools --gzvcf "genome.vcf.gz" --freq --out "genome" \
        || { echo "getting allele frequency for genome failed"; exit 1; }
    zcat "genome.primitives.vcf.gz" | \
        vcfnoindels | \
        vcftools --vcf - --freq --out "genome.primitives" \
        || { echo "getting allele frequency for primitives failed"; exit 1; }

else

    sfile="${SCRIPTDIR}/$(date '+%Y-%m-%d-%H%M')-$(basename $0)"
    cp $0 $sfile

    cd $WORKDIR
    qsub -A "tiffinp" -W group_list=tiffinp \
        -q $QUEUE -l walltime=$WALLTIME $sfile
    echo $WORKDIR

fi
