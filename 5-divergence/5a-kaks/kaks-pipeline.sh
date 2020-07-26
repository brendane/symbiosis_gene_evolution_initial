#!/bin/bash
#
# Calculate distances between protein sequences in orthogroups
# identified by OrthoFinder. Runs on untrimmed alignments.
#
# libsequence gestimator calculates Ka/Ks using Comeron's method. The
# software is no longer supported by Kevin Thornton, because it does
# not take into account high-throughput sequencing error rates, but I
# don't think there are any problems with the program.
#
# This script is similar to the fastme script, but it calculates
# produces a different output format, it calculates syn. and non-syn.
# dists, it runs on codon alignments instead of protein alignments, and
# it looks at orthosets instead of orthogroups, to save time.
# The codon alignments will have some weird sections in more divergent
# orthogroups.
#
# SETTINGS
#
QUEUE="small"
WALLTIME="01:30:00"
NTHREADS=1
MEM="4GB"
#
REMOVE_GAPS=""   # Don't remove gaps from the entire alignment (set to "-g" to remove gaps)
MAX_HITS=3       # Max hits per codon (default=3)

PROJDIR="/home/tiffinp/epste051/project/symbiosis_gene_evolution"
CODON_ALN_DIR="${PROJDIR}/results/gene_aligns/73strains_alpha_complete/muscle/2019-09-25_orthogroups/aligned-codon-orthosets"

OUTDIR="${PROJDIR}/results/phylo/73strains_alpha_complete/gestimator_dists/2019-09-25_orthosets"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
SCRIPTDIR="${OUTDIR}/script_copies"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load molpopgen

    set -euo pipefail

    SECONDS=0

    cd "$OUTDIR"

    for f in "$CODON_ALN_DIR/"*.fasta; do
        gene="$(basename "${f/.fasta/}")"
        gestimator \
            "$REMOVE_GAPS" \
            -m "$MAX_HITS" \
            -i "$f" \
            -o "${gene}.dists" \
            || { echo "running gestimator on ${f} failed"; exit 1; }
    done

    rm -f "dists.tsv"
    echo "orthoset	gene1	gene2	ka	ks	ka_ks" > "dists.tsv"
    for d in *.dists; do
        gene="${d/.dists/}"
        awk '{print "'$gene'" "\t" $0 };' "$d" | \
            sed 's/\t999\t/\tNaN\t/g' | \
            sed 's/\t999$/\tNaN/g' \
            >> "dists.tsv" \
            || { echo "adding dists from ${d} failed"; exit 1; }
    done

    #rm *.dists

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
