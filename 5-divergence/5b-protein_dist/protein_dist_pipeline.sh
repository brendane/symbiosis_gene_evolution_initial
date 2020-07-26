#!/bin/bash
#
# Calculate distances between protein sequences in orthogroups
# identified by OrthoFinder. Runs on untrimmed alignments.
#
# SETTINGS
#
QUEUE="small"
WALLTIME="04:00:00"
NTHREADS=1
MEM="4GB"
#
MODEL="JTT"             # Protein evolution model (most commonly chosen by IQ-tree)

PROJDIR="/home/tiffinp/epste051/project/symbiosis_gene_evolution"
PROT_ALN_DIR="${PROJDIR}/results/gene_aligns/73strains_alpha_complete/muscle/2019-09-25_orthogroups/aligned-aa"

OUTDIR="${PROJDIR}/results/phylo/73strains_alpha_complete/fastme/2019-09-25_orthogroups"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
SCRIPTDIR="${OUTDIR}/script_copies"

convert="${PROJDIR}/script/bin/fasta2fastmephylip.py"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load fastme/2.1.5

    set -euo pipefail

    SECONDS=0

    cd "$OUTDIR"

    cp "$convert" .
    for f in "$PROT_ALN_DIR/"*.fasta; do
        gene="$(basename "${f/.fasta/}")"
        "./$(basename "$convert")" "$f" > "${gene}.phy" \
            || { echo "convering ${f} to phylip format failed"; exit 1; }
        fastme \
            -i "${gene}.phy" \
            -p "$MODEL" \
            -c \
            -T "$NTHREADS" \
            -O "${gene}.dists" \
            || { echo "running fastme on ${f} failed"; exit 1; }
        rm "${gene}.phy"
    done

    echo "RUN TIME $SECONDS ($(($SECONDS/60)) minutes ($(($SECONDS/3600))) hours)"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$WORKDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y%m%d-%H%M')-$(basename $0)"
    cp "$0" "$sfile"

    cd "$WORKDIR"
    qsub -A "tiffinp" -W group_list="tiffinp" \
        -q "$QUEUE" -l "walltime=${WALLTIME},nodes=1:ppn=${NTHREADS},mem=${MEM}" \
        "$sfile"
    echo "$WORKDIR"

fi
