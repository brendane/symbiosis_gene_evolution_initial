#!/bin/bash
#
# Create alignments of orthologous genes identified by OrthoFinder.
#
# SETTINGS
#
QUEUE="small"
WALLTIME="00:30:00"
MEM=4GB
#
N_BATCHES=50

PROJDIR="${HOME}/project/symbiosis_gene_evolution"
ORTHODIR="${PROJDIR}/results/ortho/74strains/orthofinder/2019-09-24"
ORTHOSETS="${ORTHODIR}/Results_Directory/reconciliation_Sep25/orthosets.73strains_alpha_complete.tsv"

OUTDIR="${PROJDIR}/results/gene_aligns/73strains_alpha_complete/muscle/2019-09-25_orthogroups"
ARRAYDIR="${OUTDIR}/arrayjobdata"
ASSMDIR="${PROJDIR}/data/assembly"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
SCRIPTDIR="${OUTDIR}/script_copies"

makeseqs="${PROJDIR}/script/bin/extract_orthosets_seqs2.py"
codonaligner="${PROJDIR}/script/bin/translated_alignment_muscle.py"
trim_gaps="${PROJDIR}/script/bin/trim_gaps_fasta.py"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load muscle/3.8.31

    set -euo pipefail

    SECONDS=0

    TASK="$(cut -f 1 "${ARRAYDIR}/${PBS_ARRAYID}")"
    BATCH="$(cut -f 2 "${ARRAYDIR}/${PBS_ARRAYID}")"

    cd "$OUTDIR"

    case "$TASK" in

        "fasta")
            ## Copy amino acid sequences
            mkdir -p "fasta-aa" && cd "fasta-aa"
            cp "${ORTHODIR}/fasta/"*.fasta . \
                || { echo "copying protein sequences failed"; exit 1; }
            cd ..

            ## Copy nucleotide sequences and modify if needed
            mkdir -p "fasta-nucl" && cd "fasta-nucl"
            for fname in "${ORTHODIR}/fasta/"*.fasta; do
                strain="$(basename "${fname/.fasta/}")"
                if [[ -f "${ASSMDIR}/ncbi/${strain}/genes.fasta" ]]; then
                    sed 's/>.\+\[locus_tag=/>/g' "${ASSMDIR}/ncbi/${strain}/genes.fasta" | \
                        sed 's/\].\+//g' > \
                        "${strain}.fasta" \
                        || { echo "copying nucleotide gene sequences for ${strain} failed"; exit 1; }
                else
                    sed 's/|.\+//g' "${ASSMDIR}/mage/${strain}/genes.fasta" > "${strain}.fasta" \
                        || { echo "copying nucleotide gene sequences for ${strain} failed"; exit 1; }
                fi
            done
            cd ..
            ;;

        "align-prep")
            ## Concatenate all fasta sequences for nucl. and all for aa,
            ## making sure that the strain prefix is added.
            ## Then, use the makeseqs script to extract unaligned files
            ## for each target gene.

            ## Put all the fasta sequences together into one big file
            rm -f "aa.fasta"
            for f in fasta-aa/*.fasta; do
                strain="$(basename "${f/.fasta/}")"
                sed 's/^>/>'"$strain"'_/g' "$f" >> "aa.fasta" \
                    || { echo "adding ${strain} to aa db failed"; exit 1; }
            done

            ## Put all the fasta sequences together into one big file
            rm -f "nucl.fasta"
            for f in fasta-nucl/*.fasta; do
                strain="$(basename "${f/.fasta/}")"
                sed 's/^>/>'"$strain"'_/g' "$f" >> "nucl.fasta" \
                    || { echo "adding ${strain} to nucleotide db failed"; exit 1; }
            done

            ## Create unaligned groups of sequences
            cp "$makeseqs" .
            for seqtype in "aa" "nucl"; do
                rm -rf "unaligned-${seqtype}" "aligned-${seqtype}"
                mkdir -p "unaligned-${seqtype}" "aligned-${seqtype}"
                "./$(basename "$makeseqs")" --orthogroups \
                    "$ORTHOSETS" "${seqtype}.fasta" "unaligned-${seqtype}" \
                    || { echo "extracting ${seqtype} failed"; exit 1; }
            done
            rm -rf "aligned-codon"
            mkdir -p "aligned-codon"

            for seqtype in "nucl" "aa"; do
                rm -rf "unaligned-${seqtype}-orthosets" "aligned-${seqtype}-orthosets"
                mkdir -p "unaligned-${seqtype}-orthosets" "aligned-${seqtype}-orthosets"
                "./$(basename "$makeseqs")" \
                    "$ORTHOSETS" "${seqtype}.fasta" "unaligned-${seqtype}-orthosets" \
                    || { echo "extracting ${seqtype} failed"; exit 1; }
            done
            rm -rf "aligned-codon-orthosets"
            mkdir -p "aligned-codon-orthosets"
            rm -rf "aligned-aa-orthosets"
            mkdir -p "aligned-aa-orthosets"
            ;;

        "align-aa")
            ls unaligned-aa | awk '{ print "unaligned-aa/" $1};' > "files.all.${BATCH}.txt"
            split -n "r/$BATCH/$N_BATCHES" "files.all.${BATCH}.txt" > \
                "files.only.${BATCH}.txt" \
                || { echo "splitting files for ${BATCH} in ${GROUP} failed"; exit 1; }
            while read -r f; do
                muscle -in "$f" > "aligned-aa/$(basename "$f")" \
                    || { echo "protein alignment on ${f} for ${GROUP} (batch ${BATCH}) failed"; exit 1; }
            done < "files.only.${BATCH}.txt"

            ## Delete the unaligned sequences (do this manually)
            #rm -r "unaligned-aa"
            #rm files.*.txt
            ;;

        "align-aa-orthosets")
            ls "unaligned-aa-orthosets" | awk '{ print "unaligned-aa-orthosets/" $1};' > "files.all.${BATCH}.txt"
            split -n "r/$BATCH/$N_BATCHES" "files.all.${BATCH}.txt" > \
                "files.only.${BATCH}.txt" \
                || { echo "splitting files for ${BATCH} in ${GROUP} failed"; exit 1; }
            while read -r f; do
                muscle -in "$f" > "aligned-aa-orthosets/$(basename "$f")" \
                    || { echo "protein alignment on ${f} for ${GROUP} (batch ${BATCH}) failed"; exit 1; }
            done < "files.only.${BATCH}.txt"
            ;;

        "align-codon")
            ls "unaligned-nucl" | awk '{ print "unaligned-nucl/" $1};' > "files.all.${BATCH}.txt"
            split -n "r/$BATCH/$N_BATCHES" "files.all.${BATCH}.txt" > \
                "files.only.${BATCH}.txt" \
                || { echo "splitting files for ${BATCH} failed"; exit 1; }
            while read -r f; do
                if [[ -s "$f" ]]; then
                    "$codonaligner" --protein-aln "aligned-aa/$(basename "$f")" "$f" > \
                        "aligned-codon/$(basename "$f")" \
                        || { echo "codon alignment on ${f} (batch ${BATCH}) failed"; exit 1; }
                fi
            done < "files.only.${BATCH}.txt"

            ## Delete the unaligned sequences (do this manually)
            #rm -r "unaligned-nucl"
            #rm "files."*".${BATCH}.txt"
            ;;

        "align-codon-orthosets")
            ls "unaligned-nucl-orthosets" | awk '{ print "unaligned-nucl-orthosets/" $1};' > "files.all.${BATCH}.txt"
            split -n "r/$BATCH/$N_BATCHES" "files.all.${BATCH}.txt" > \
                "files.only.${BATCH}.txt" \
                || { echo "splitting files for ${BATCH} failed"; exit 1; }
            while read -r f; do
                if [[ -s "$f" ]]; then
                    "$codonaligner" --protein-aln "aligned-aa/$(basename "$f" | sed 's/\..\+/.fasta/g')" \
                        "$f" | \
                        "$trim_gaps" > \
                        "aligned-codon-orthosets/$(basename "$f")" \
                        || { echo "codon alignment on ${f} (batch ${BATCH}) failed"; exit 1; }
                fi
            done < "files.only.${BATCH}.txt"

            ## Delete the unaligned sequences (do this manually)
            #rm -r "unaligned-nucl"
            #rm "files."*".${BATCH}.txt"
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
        "fasta")
            WALLTIME="00:15:00"
            MEM="2GB"
            NTHREADS=1
            echo "$1"$'\t'"_" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
        "align-prep")
            WALLTIME="00:30:00"
            echo "$1"$'\t'"_" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
        "align-aa")
            for batch in $(seq 1 $N_BATCHES); do
                echo "$1"$'\t'"$batch" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;
        "align-aa-orthosets")
            for batch in $(seq 1 $N_BATCHES); do
                echo "$1"$'\t'"$batch" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;
        "align-codon")
            WALLTIME="00:05:00"
            for batch in $(seq 1 $N_BATCHES); do
                echo "$1"$'\t'"$batch" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;
        "align-codon-orthosets")
            WALLTIME="01:05:00"
            for batch in $(seq 1 $N_BATCHES); do
                echo "$1"$'\t'"$batch" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;
    esac

    if [[ "$i" == 1 ]]; then
        array="0"
    else
        array="0-$(($i-1))"
    fi

    cd "$WORKDIR"
    qsub -A "tiffinp" -W group_list="tiffinp" \
        -q "$QUEUE" -l "walltime=${WALLTIME},mem=${MEM}" -t "$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
