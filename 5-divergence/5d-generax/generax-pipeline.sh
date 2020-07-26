#!/bin/bash
#
# Run GeneRax on all orthosets to infer duplication, transfer, and loss
# events. 
#
# GeneRax takes a species tree, gene trees, and gene sequence alignments,
# optimizes the gene trees using the species tree, and estimates rates of
# duplication, loss, and horizontal transfer. It is set up to estimate
# rates across all genes simultaneously, but I want rates for each
# gene, and running on all genes at once takes a very long
# time (even with --per-family-rates). So, I run GeneRax separately for
# each gene.
#
# SETTINGS
#
QUEUE=mesabi
WALLTIME=12:00:00
NTHREADS=1
MEM=4GB
#
MODEL=UndatedDTL  # reconciliation model; UndatedDTL looks for all three types of events (default=UndatedDL)
STRATEGY=SPR      # optimize gene tree topology (SPR, default) or not (EVAL)
SUBST_MODEL=GTR+G # substitution model; default is GTR, manual recommends adding "+G" to models
MAX_SPR_RADIUS=5  # Maximum radius for subtree pruning and regrafting moves (default=5)
MIN_TAXA=4        # Minimum number of species in a gene (orthoset) to run GeneRax
MIN_SEQS=6        # Minimum number of sequences in a gene (orthoset) to run GeneRax
#REC_OPT=Simplex   # Default option for this version for optimization algorithm
N_BATCHES=200

PROJDIR="/home/tiffinp/epste051/project/symbiosis_gene_evolution"
CODON_ALN_DIR="${PROJDIR}/results/gene_aligns/73strains_alpha_complete/muscle/2019-09-25_orthogroups/aligned-codon-orthosets"
SPP_TREE="${PROJDIR}/results/phylo/73strains_alpha_complete/iqtree/2019-09-26/single_copy_core_tree/sc.rooted.figtree.2020-04-24.nw" # Manually rooted by using figtree
GENE_TREES="${PROJDIR}/results/ortho/74strains/orthofinder/2019-09-24/Results_Directory/reconciliation_Sep25/RootedTrees.all.tsv" # See 2-orthologs/rooted_gene_trees_from_orthofinder.txt
PSEUDOGENES="${PROJDIR}/results/ortho/74strains/orthofinder/2019-09-24/pseudogenes.tsv" # See 2-orthologs/pseudogenes.txt

OUTDIR="${PROJDIR}/results/reconciliation/73strains_alpha_complete/generax/2020-04-23"
ARRAYDIR="${OUTDIR}/arrayjobdata"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
SCRIPTDIR="${OUTDIR}/script_copies"
SCRATCHDIR="/scratch.global/${USER}/generax_2020-04-23"

prune="${PROJDIR}/script/bin/prune_phylo_tips.py"
filter_fasta="${PROJDIR}/script/bin/filter_fasta_names.py"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    set -euo pipefail

    SECONDS=0

    TASK="$(cut -f 1 "${ARRAYDIR}/${PBS_ARRAYID}")"
    BATCH="$(cut -f 2 "${ARRAYDIR}/${PBS_ARRAYID}")"

    cd "$SCRATCHDIR"

    case "$TASK" in

        "prep")

            ## Make families file
            cp "$prune" .
            cp "$filter_fasta" .
            rm -rf output
            rm -f inputs.txt
            rm -rf family_files
            rm -rf alignments
            mkdir -p family_files
            mkdir -p mappings
            mkdir -p starting_trees
            mkdir -p alignments
            for f in "$CODON_ALN_DIR/"*.fasta; do
                gene="$(basename "${f/.fasta/}")"
                og="$(echo "$gene" | sed 's/\..\+//g')"

                if [[ -f "family_files/${gene}.txt" ]]; then continue; fi
                
                "./$(basename "$filter_fasta")" -x \
                    "$f" "$PSEUDOGENES" \
                    > "alignments/${gene}.fasta" \
                    || { echo "removing pseudogenes from ${gene} failed"; exit 1; }

                grep ">" "alignments/${gene}.fasta" | \
                    sed 's/^>//g' | \
                    awk '{ print $1 " " gensub(/_.+/, "", "g", $1)};' \
                    > "mappings/${gene}.link" \
                    || { echo "making gene-species mapping file for ${gene} was not happy, might have too few seqs"; }
                cut -f 1 -d" " "mappings/${gene}.link" > tmp2

                if [[ $(sort -u tmp2 | wc -l | cut -f 1 -d" ") -lt $MIN_SEQS ]]; then continue; fi
                if [[ $(cut -f 2 -d" " "mappings/${gene}.link" | sort -u | wc -l | cut -f 1 -d" ") -lt $MIN_TAXA ]]; then continue; fi

                grep "^${og}" "$GENE_TREES" | \
                    cut -f 2 > tmp \
                    || { echo "finding ${og} in gene trees file failed"; continue; }

                ## Prune the gene tree down to just the sequences of interest
                "./$(basename "$prune")" tmp tmp2 > "starting_trees/${gene}.nw" \
                    || { echo "making starting tree for ${gene} failed"; exit 1; }
                rm tmp tmp2

                echo "[FAMILIES]" > "family_files/${gene}.txt"
                echo "- ${gene}" >> "family_files/${gene}.txt"
                echo "starting_gene_tree = starting_trees/${gene}.nw" >> "family_files/${gene}.txt"
                echo "alignment = alignments/${gene}.fasta" >> "family_files/${gene}.txt"
                echo "mapping = mappings/${gene}.link" >> "family_files/${gene}.txt"
                echo "subst_model = ${SUBST_MODEL}" >> "family_files/${gene}.txt"
            done

            cp "$SPP_TREE" species_tree.nw

            ls family_files/* > files.txt

            rm -rf GeneRax_Output
            ;;

        "run")

            source activate /home/tiffinp/epste051/modules/modules/generax/1.2.0/generax-env
            mkdir -p output

            split -n l/$BATCH/$N_BATCHES files.txt | \
                while read -r ff; do
                ## Notes on other options
                ## --reconciliation-samples: currently undocumented feature
                ##      that might be meant to report multiple reconciliation
                ##      scenarios; ignoring for now b/c not in the documentation
                ## --support-threshold: Also undocumented. I think this has to
                ##      do with the minimum bootstrap support threshold for
                ##      tree nodes during the reconciliation step and the default
                ##      is not to use a threshold
                ## --rec-weight: Not clear what this does or if it currently does
                ##      anything at all.
                rm -rf "output/GeneRax_Output_$(basename "${ff/.txt/}")"
                generax \
                    --prefix "output/GeneRax_Output_$(basename "${ff/.txt/}")" \
                    --rec-model "$MODEL" \
                    --species-tree species_tree.nw \
                    --strategy "$STRATEGY" \
                    --families "$ff" \
                    --max-spr-radius "$MAX_SPR_RADIUS" \
                    --reconcile \
                    || { echo "running GeneRax on ${ff} failure, skipping"; continue; }
#                    || { echo "running GeneRax on ${ff} failed"; exit 1; }
                    #--rec-opt "$REC_OPT" \
            done

            source deactivate

            ;;

        summary)
            rm -f dtl.tsv
            echo "gene	duplication	loss	transfer" > dtl.tsv
            for d in output/GeneRax_Output_*; do
                gene="$(basename "${d/GeneRax_Output_/}")"
                statfile="${d}/results/${gene}/stats.txt"
                if [[ ! -s "$statfile" ]]; then continue; fi
                awk 'NR > 1 {print "'"$gene"'\t" $4 "\t" $5 "\t" $6};' \
                    "$statfile" \
                    >> dtl.tsv \
                    || { echo "getting DTL rates for ${gene} failed"; exit 1; }
            done

            cp -r dtl.tsv output "$OUTDIR" \
                || { echo "copying output failed"; exit 1; }

    esac

    echo "RUN TIME $SECONDS ($(($SECONDS/60)) minutes ($(($SECONDS/3600))) hours)"

    rm "${ARRAYDIR}/${PBS_ARRAYID}"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$ARRAYDIR"
    mkdir -p "$SCRATCHDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y%m%d-%H%M')${1}-$(basename $0)"
    cp "$0" "$sfile"

    i=0
    case "$1" in
        "prep")
            WALLTIME="07:00:00"
            MEM="2GB"
            echo "$1"$'\t'"_" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
        "run")
            for j in $(seq 1 $N_BATCHES); do
                echo "$1"$'\t'"$j" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;
        "summary")
            WALLTIME="02:00:00"
            MEM="1GB"
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
        -q "$QUEUE" -l "walltime=${WALLTIME},mem=${MEM}" \
        -t "$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
