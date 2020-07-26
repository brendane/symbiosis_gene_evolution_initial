#!/bin/bash
#
# Make phylogenies for a concatenated single copy core alignment.
#
# Pseudogenes are not included.
#
# SETTINGS
#
QUEUE="small"
WALLTIME="07:00:00"
NTHREADS=16
MEM="32GB"
#
MODEL="MFP+MERGE"       # Automatically choose the best protein evolution model
MODEL_SC="JTT+G4"       # Model for concatenated alignment -- most common for single genes; chosen to avoid lengthy model selection step
N_BOOTS=1000            # Number of bootstraps
N_SCF=200               # Number of site concordance quartets

PROJDIR="/home/tiffinp/epste051/project/symbiosis_gene_evolution"
PROT_ALN_DIR="${PROJDIR}/results/gene_aligns/73strains_alpha_complete/muscle/2019-09-25_orthogroups/aligned-aa-orthosets"
SC_GENES="${PROJDIR}/notes/table/single_copy_core_genes.2019-09-24.nopseudo.tsv"

OUTDIR="${PROJDIR}/results/phylo/73strains_alpha_complete/iqtree/2019-09-26"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
ARRAYDIR="${OUTDIR}/arrayjobdata"
SCRIPTDIR="${OUTDIR}/script_copies"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load gcc/6.3.0 # For pre-compiled bljoin (higher versions work too)
    module load iqtree/1.6.10
    module load trimal/1.4.1

    set -euo pipefail

    TASK=$(cut -f 1 "${ARRAYDIR}/${PBS_ARRAYID}")

    SECONDS=0

    cd "$OUTDIR"

    case "$TASK" in
        "prep")
            ## Concatenate the single-copy core genes.
            ## Uses custom joining program.
            mkdir -p "tmp"
            tail -n +2 "$SC_GENES" | cut -f 2 > "single_copy_core_genes.txt" \
                || { echo "listing single-copy core genes failed"; exit 1; }
            while read -r gene; do
                trimal -in "${PROT_ALN_DIR}/${gene}.fasta" -automated1 -out "tmp/${gene}.fasta" -fasta \
                    || { echo "running trimAl on ${gene} failed"; exit 1; }
            done < "single_copy_core_genes.txt"
            bljoin -d "_" -f 1 tmp/*.fasta > "single_copy_core.concatenated.fasta" \
                || { echo "concatenating sequences failed"; exit 1; }
            rm -r "tmp"
            ;;

        "run-sc")
            mkdir -p "single_copy_core_tree"
            rm -f "sc.rooted_tree.figtree.nw"
            iqtree -s "single_copy_core.concatenated.fasta" \
                -pre "single_copy_core_tree/sc" \
                -m "$MODEL_SC" \
                -bb "$N_BOOTS" \
                -alrt "$N_BOOTS" \
                -nt "$NTHREADS" \
                || { echo "running iqtree on concatenated SC genes failed"; exit 1; }
                rm "single_copy_core_tree/sc.ckp.gz"
            
            set +euo pipefail
            module unload iqtree
            module load iqtree/1.7b12
            set -euo pipefail
            iqtree -s "single_copy_core.concatenated.fasta" \
                -t "single_copy_core_tree/sc.treefile" \
                --scf "$N_SCF" \
                -pre "single_copy_core_tree/sc" \
                -nt "$NTHREADS" \
                || { echo "estimating site concordance for species tree failed"; exit 1; }
            ;;

    esac

    echo "RUN TIME $SECONDS ($(($SECONDS/60)) minutes ($(($SECONDS/3600))) hours)"

    rm "${ARRAYDIR}/${PBS_ARRAYID}"

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
            NTHREADS=1
            MEM="2GB"
            echo "${1}" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
        "run-sc")
            NTHREADS=$NTHREADS
            echo "${1}" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
    esac


    if [[ "$i" == 1 ]]; then
        array="0"
    else
        array="0-$(($i-1))"
    fi

    cd "$WORKDIR"
    qsub2 -A "tiffinp" -W group_list="tiffinp" \
        -q "$QUEUE" -l "walltime=${WALLTIME},nodes=1:ppn=${NTHREADS},mem=${MEM}" \
        -t "$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
