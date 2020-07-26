#!/bin/bash
#
# Test for phylogenetic signal in all genes.
#
# SETTINGS
#
QUEUE=mesabi
WALLTIME=03:00:00
NTHREADS=1
MEM=4GB
#
LAMBDA=0.1
N_SIM=10000  # Total number of MCMC iterations
THIN=10      # Keep every $THIN iterations after the burnin
BURN=1000    # Number of iterations ignored as burn-in
SE=0.5
N_BATCHES=10

PROJDIR="/home/tiffinp/epste051/project/symbiosis_gene_evolution"
SPP_TREE="${PROJDIR}/results/phylo/73strains_alpha_complete/iqtree/2019-09-26/single_copy_core_tree/sc.rooted.figtree.2020-04-24.nw" # Manually rooted by using figtree
ORTHO_FILE="${PROJDIR}/results/gene_tables/divergence/73strains_alpha_complete/2019-09-26/gene_distances.nopseudo.orthosets.tsv"     # This is actually from step 5f
PAIRWISE_DIV="${PROJDIR}/results/gene_tables/divergence/73strains_alpha_complete/2019-09-26/gene_distances.nopseudo.orthosets_pairwise_sequence.tsv" # Also from step 5f
GENE_DATA="${PROJDIR}/results/gene_tables/divergence/73strains_alpha_complete/2019-09-26/gene_distances.nopseudo.orthosets.tsv" # Likewise; and I don't know why I listed it twice!

OUTDIR="${PROJDIR}/results/phylo/73strains_alpha_complete/phylo_signal/2020-05-03"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
SCRIPTDIR="${OUTDIR}/script_copies"
ARRAYDIR="${OUTDIR}/arrayjobdata"

# Makes copy number table
copy_table="${PROJDIR}/script/bin/copy_number_table.py"
# Calculates delta
phylo_signal="${PROJDIR}/script/bin/delta_phylogenetic_signal.r"
# Calculates R^2 in species distance vs. gene distance regressions
r2="${PROJDIR}/script/bin/phylogenetic_distance_correlation.r"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    set -euo pipefail

    SECONDS=0

    TASK="$(cut -f 1 "${ARRAYDIR}/${PBS_ARRAYID}")"
    BATCH="$(cut -f 2 "${ARRAYDIR}/${PBS_ARRAYID}")"

    cd "$OUTDIR"

    case "$TASK" in
        prep)
            cp "$SPP_TREE" species_tree.nw

            cp "$copy_table" .
            "./$(basename "$copy_table")" "$ORTHO_FILE" genes subset > \
                copy_number.tsv \
                || { echo "making copy number table failed"; exit 1; }

            awk -F$'\t' \
                'NR==1 {for(i=2; i<=NF; i++) {header[i]=$i}}; NR>1 {for(i=2; i<=NF; i++) {print header[i] "\t" $i > $1".trait"}}' \
                copy_number.tsv \
                || { echo "making trait files failed"; exit 1; }

            cp "$phylo_signal" .
            ;;

        run)
            rm -f signal.$BATCH.tsv
            ls OG*.trait | split -n r/$BATCH/$N_BATCHES > files.$BATCH \
                || { echo "making list of files for ${BATCH} failed"; exit 1; }
            while read -r file; do
                ## NOTE: This was also run with the --pav option to calculate
                ## delta based on P/A only.
                "./$(basename "$phylo_signal")" species_tree.nw "$file" \
                    "$LAMBDA" "$SE" "$N_SIM" " $THIN" "$BURN" | \
                    awk '{ print "'"${file/.trait/}"'" "\t" $0 };' >> \
                    signal.$BATCH.tsv \
                    || { echo "testing ${file} failed"; exit 1; }
                rm "$file"
            done < files.$BATCH
            rm files.$BATCH
            ;;

        r2)
            cp "$r2" .
            "./$(basename "$r2")" "$PAIRWISE_DIV" "$GENE_DATA" > r2_scc_genes.tsv \
                || { echo "calculating correlation between scc and other genes failed"; exit 1; }
            ;;

        combine)
            cat signal.*.tsv > signal.tsv
            rm signal.*.tsv
            ;;
    esac

    rm "${ARRAYDIR}/${PBS_ARRAYID}"

    echo "RUN TIME $SECONDS ($(($SECONDS/60)) minutes ($(($SECONDS/3600))) hours)"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$ARRAYDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y%m%d-%H%M')-${1}-$(basename $0)"
    cp "$0" "$sfile"

    i=0
    case "$1" in
        "prep")
            WALLTIME="00:15:00"
            echo "$1"$'\t'"_" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
        "run")
            for j in $(seq 1 $N_BATCHES); do
                echo "$1"$'\t'"$j" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;
        "r2")
            WALLTIME="00:30:00"
            echo "$1"$'\t'"_" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
        "combine")
            WALLTIME="00:05:00"
            echo "$1"$'\t'"_" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
    esac

    if [[ "$i" == 1 ]]; then
        array="0"
    else
        array="0-$(($i-1))"
    fi

    cd "$WORKDIR"
    qsub -A "tiffinp" -W group_list="tiffinp" \
        -q "$QUEUE" -l "walltime=${WALLTIME}" \
        -t "$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
