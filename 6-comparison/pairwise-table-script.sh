#!/bin/bash
#
# Calculate distances between protein sequences in orthogroups
# identified by OrthoFinder. Runs on untrimmed alignments.
#
# Second OrthoFinder run that was done to include USDA1106.
#
# RERUN 2020-05-05 with some bug fixes.
#
# SETTINGS
#
QUEUE="small"
WALLTIME="15:00:00"
NTHREADS=8
MEM="32GB"

PROJDIR="/home/tiffinp/epste051/project/symbiosis_gene_evolution"

STRAINS="${PROJDIR}/data/strain_lists/73strains_complete_alpha"
ORTHODIR="${PROJDIR}/results/ortho/74strains/orthofinder/2019-09-24/Results_Directory/reconciliation_Sep25"
ORTHOSETS="${ORTHODIR}/orthosets.73strains_alpha_complete.tsv"
ORTHOSETS_BY_DISTANCE="${ORTHODIR}/orthosets.73strains_alpha_complete.dist_0.5.tsv"
PSEUDOGENES="${ORTHODIR}/../../pseudogenes.tsv"
ANNOTATION="${ORTHODIR}/annotation.73strains_alpha_complete.tsv"
SCC_DIST_MATRIX="${PROJDIR}/results/phylo/73strains_alpha_complete/iqtree/2019-09-26/single_copy_core_tree/sc.mldist"
KAKS_ALL_PAIRWISE="${PROJDIR}/results/phylo/73strains_alpha_complete/gestimator_dists/2019-09-25_orthosets/dists.tsv"
PROTEIN_DIVERGENCES_DIR="${PROJDIR}/results/phylo/73strains_alpha_complete/fastme/2019-09-25_orthogroups"
SPECIES="${PROJDIR}/notes/table/73strains_alpha_complete.species.tsv"

OUTDIR="${PROJDIR}/results/gene_tables/divergence/73strains_alpha_complete/2019-09-26"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
SCRIPTDIR="${OUTDIR}/script_copies"

tabulate="${PROJDIR}/script/bin/make_orthoset_divergence_table.r"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    set -euo pipefail

    SECONDS=0

    cd "$OUTDIR"

    "./$(basename "$tabulate")" --output "gene_distances.nopseudo" \
        --exclude "$PSEUDOGENES" \
        --strains "$STRAINS" \
        --single-copy-core-distances "$SCC_DIST_MATRIX" \
        --kaks "$KAKS_ALL_PAIRWISE" \
        --annotation "$ANNOTATION" \
        --protein-divergence-directory "$PROTEIN_DIVERGENCES_DIR" \
        --species "$SPECIES" \
        "$ORTHOSETS" "$ORTHOSETS_BY_DISTANCE" \
        || { echo "tabulating without pseudogenes failed"; exit 1; }

    echo "RUN TIME $SECONDS ($(($SECONDS/60)) minutes ($(($SECONDS/3600))) hours)"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$WORKDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y%m%d-%H%M')-$(basename $0)"
    cp "$0" "$sfile"

    cd "$WORKDIR"
    qsub2 -A "tiffinp" -W group_list="tiffinp" \
        -q "$QUEUE" -l "walltime=${WALLTIME},nodes=1:ppn=${NTHREADS},mem=${MEM}" \
        "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
