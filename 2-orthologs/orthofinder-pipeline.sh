#!/bin/bash
#
# Run OrthoFinder on 74 rhizobial strains (73 + USDA1106).
#
# A few differences between this analysis and other OrthoFinder
# analyses that I've run on larger and more closely related datasets:
# - Default diamond search settings, except for lower evalue (1E-7)
# - No removal of identical duplicates
# - Protein coding genes were downloaded along with genome sequences,
#   so no need to run a conversion script.
#
# Originally, this script was run using OrthoFinder 2.3.3, but that
# version is not amenable to being run in small parts.
#
# SETTINGS
#
QUEUE="small"
WALLTIME="96:00:00"
MEM=62GB
NTHREADS=24
#
SEARCH="diamond"        # program to use for finding similar sequences
OTHREADS=4              # number of threads for non-blast steps; memory limited
ONE_WAY="-1"            # if -1 only search in one direction; to do 2-way set to ""
MCL_INFLATION=4         # inflation parameter for MCL (default=1.5); higher gives smaller clusters
N_BATCHES=50            # Number of tree reconcilation batches
#
ORTHO_RUN="Sep25"

PROJDIR="/home/tiffinp/epste051/project/symbiosis_gene_evolution"
ASSMDIR="${PROJDIR}/data/assembly"
STRAINS="${PROJDIR}/data/strain_lists/74strains" ## 1-genomes/all_strains.txt

declare -A STRAINSLISTS
STRAINSLISTS["73strains_alpha_complete"]="${PROJDIR}/data/strain_lists/73strains_complete_alpha" ## 1-genomes/27_chosen_strains.txt

OUTDIR="${PROJDIR}/results/ortho/74strains/orthofinder/2019-09-24"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
ARRAYDIR="${OUTDIR}/arrayjobdata"
SCRIPTDIR="${OUTDIR}/script_copies"

## Note the "$recon" script listed below has to be placed with the OrthoFinder
## sources b/c it imports a bunch of stuff from OrthoFinder. It is just a modification
## of the trees2ologs_dash.py script from OrthoFinder.
recon="${HOME}/devel/OrthoFinder-2.2.7_source/orthofinder/scripts/trees2ologs_modified_dash.py"
concat_recon="${PROJDIR}/script/bin/combine_orthofinder_duplication_batches.py"
orthotable="${PROJDIR}/script/bin/extract_orthosets_from_orthofinder_trees.py"
magegff="${PROJDIR}/script/bin/mage_gff_processor.py"
ncbigff="${PROJDIR}/script/bin/ncbi_gff_processor.py"
annot="${PROJDIR}/script/bin/combine_orthotable_and_annotations.py"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load orthofinder/2.2.7
    module load parallel

    set -euo pipefail

    SECONDS=0

    TASK="$(cut -f 1 "${ARRAYDIR}/${PBS_ARRAYID}")"
    BATCH="$(cut -f 2 "${ARRAYDIR}/${PBS_ARRAYID}")"
    GROUP="$(cut -f 3 "${ARRAYDIR}/${PBS_ARRAYID}")"

    cd "$OUTDIR"

    case "$TASK" in

        "prep")
            ## Copy protein coding gene files to working directory; also
            ## extract useful annotation information
            cp "$magegff" .
            cp "$ncbigff" .
            mkdir -p "annotation"
            mkdir -p "fasta" && cd "fasta"
            while read -r strain; do
                echo "$strain"
                f="${ASSMDIR}/ncbi/${strain}/cds.fasta"
                strain2="$(echo "$strain" | sed 's/_/-/g')"
                if [[ ! -f "$f" ]]; then
                    f="${ASSMDIR}/mage/${strain}/cds.fasta"
                    sed 's/|.\+//g' "$f" > "${strain2}.fasta" \
                        || { echo "copying protein gene sequences for ${strain} failed"; exit 1; }
                    "../$(basename "$magegff")" "${ASSMDIR}/mage/${strain}/genome.gff3" \
                        "${strain2}_" > "../annotation/${strain2}.tsv" \
                        || { echo "extracting annotation information for ${strain} failed"; exit 1; }
                else
                    sed 's/>.\+\[locus_tag=/>/g' "$f" | sed 's/\].\+//g' > \
                        "${strain2}.fasta" \
                        || { echo "copying protein gene sequences for ${strain} failed"; exit 1; }
                    "../$(basename "$ncbigff")" "${ASSMDIR}/ncbi/${strain}/genome.gff3" \
                        "${strain2}_" > "../annotation/${strain2}.tsv" \
                        || { echo "extracting annotation information for ${strain} failed"; exit 1; }
                fi
            done < "$STRAINS"
            
            rm -f "../annotation/all.tsv"
            cat "../annotation/"* > "../annotation/all.tsv"
            ;;

        "blast")
            ## Set up directory structure and DIAMOND databases
            orthofinder -f "fasta" -op -S "$SEARCH" \
                -t "$NTHREADS" "$ONE_WAY" \
                -n "Directory" | tee "blast_setup_log.txt" \
                || { echo "setting up for blast failed"; exit 1; }

            ## Run DIAMOND in parallel (commands from log file)
            ## Note that parallel requires a citation
            grep "^diamond blastp" "blast_setup_log.txt" > "diamond.cmds" \
                || { echo "getting diamond commands failed"; exit 1; }
            cat "diamond.cmds" | parallel -j "$NTHREADS" \
                || { echo "running diamond failed"; exit 1; }

            ## Move silly location of output to something more useful
            current_results="$(ls -t "fasta" | grep "^Results_Directory" | head -n 1)"
            mv "fasta/${current_results}" "Results_Directory" \
                || { echo "moving output directory failed"; exit 1; }
            ;;

        "orthogroup")
            orthofinder "$ONE_WAY" -a "$OTHREADS" \
                -I "$MCL_INFLATION" \
                -b "Results_Directory/WorkingDirectory" \
                -og -t "$NTHREADS" \
                || { echo "getting orthogroups failed"; exit 1; }
            ;;

        "tree")
            orthofinder "$ONE_WAY" \
                -fg "Results_Directory/WorkingDirectory" \
                -n "Directory" -a "$OTHREADS" -t "$OTHREADS" -ot \
                || { echo "constructing trees failed"; exit 1; }
            ;;

        "fix-trees")
            ## Apparently, OrthoFinder changes "." to "_" in gene trees (but not
            ## anywhere else). Really stupid and irritating.
            for f in "Results_Directory/WorkingDirectory/Orthologues_Directory_${ORTHO_RUN}/Gene_Trees/"*; do
                cat "$f" | \
                    sed 's/th_b2/th.b2/g' | \
                    sed 's/Cp5_3/Cp5.3/g' | \
                    sed 's/IPA7_2/IPA7.2/g' > \
                    "tmp" \
                    || { echo "fixing ${f} failed"; exit 1; }
                mv "tmp" "$f" \
                    || { echo "renaming ${f} failed"; exit 1; }
            done
            ;;

        "reconcile")
            species_tree="Results_Directory/WorkingDirectory/Orthologues_Directory_${ORTHO_RUN}/SpeciesTree_rooted.txt"
            mkdir -p "Results_Directory/reconciliation_${ORTHO_RUN}/"
            cp "$recon" .
            PYTHONPATH="$HOME/usr/lib/python2.7/site-packages.old" "$recon" \
                --output "Results_Directory/reconciliation_${ORTHO_RUN}/" \
                "Results_Directory/WorkingDirectory/Orthologues_Directory_${ORTHO_RUN}/Gene_Trees/" \
                "$species_tree" \
                "Results_Directory/WorkingDirectory/SpeciesIDs.txt" \
                "Results_Directory/WorkingDirectory/SequenceIDs.txt" \
                "$N_BATCHES" \
                "$BATCH" \
                || { echo "reconciliation on batch ${BATCH} failed"; exit 1; }
            ;;

        "combine-reconciled")
            cp "$concat_recon" .
            cd "Results_Directory/reconciliation_${ORTHO_RUN}"
            "../../$(basename "$concat_recon")" \
                --output "all_duplications.tsv" \
                "." \
                "../../fasta" \
                "../WorkingDirectory/Orthogroups.csv" \
                || { echo "combining reconciled files failed"; exit 1; }
            rm -f "RootedTrees.all.tsv"
            cat RootedTrees.*.tsv | sort -k 1b,1 -t$'\t' > \
                "RootedTrees.all.tsv" \
                || { echo "combining tree files failed"; exit 1; }
            ;;

        "orthosets")
            cd "Results_Directory/reconciliation_${ORTHO_RUN}"
            "$orthotable" \
                --output "orthosets.${GROUP}.tsv" \
                "${STRAINSLISTS[$GROUP]}" \
                "RootedTrees.all.tsv" \
                "../WorkingDirectory/Orthogroups.csv" \
                "../WorkingDirectory/Orthogroups_UnassignedGenes.csv" \
                || { echo "making orthosets for ${GROUP} failed"; exit 1; }
            ;;

        "annotation")
            "$annot" \
                "Results_Directory/reconciliation_${ORTHO_RUN}/orthosets.${GROUP}.tsv" \
                "annotation/all.tsv" \
                "1021,WSM419,USDA1106" \
                > "Results_Directory/reconciliation_${ORTHO_RUN}/annotation.${GROUP}.tsv" \
                || { echo "annotating ${GROUP} failed"; exit 1; }
            if [[ "$GROUP" == "73strains_alpha_complete" ]]; then
                awk -F$'\t' '$5 == 1 { print $1 };' "annotation/all.tsv" > \
                    "pseudogenes.tsv" \
                    || { echo "making list of pseudogenes failed"; exit 1; }
            fi
            ;;

    esac

    echo "RUN TIME $SECONDS ($(($SECONDS/60)) minutes ($(($SECONDS/3600))) hours)"
    echo "done"

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
            MEM="1GB"
            NTHREADS=1
            echo "$1"$'\t'"_" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
        "blast")
            WALLTIME="01:30:00"  # Actually took 60 min
            echo "$1"$'\t'"_" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
        "orthogroup")
            WALLTIME="01:00:00" # Actually took 40 min
            NTHREADS="$OTHREADS"
            echo "$1"$'\t'"_" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
        "tree")
            WALLTIME="01:00:00" # Actually took 30 min
            NTHREADS="$OTHREADS"
            echo "$1"$'\t'"_" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
        "fix-trees")
            WALLTIME="01:00:00"
            NTHREADS=1
            echo "$1"$'\t'"_" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
        "reconcile")
            NTHREADS="1"
            MEM="4GB"
            WALLTIME="00:05:00" # Actually took a few seconds
            for b in $(seq 0 $(($N_BATCHES-1))); do
                echo "$1"$'\t'"$b" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;
        "combine-reconciled")
            NTHREADS="1"
            MEM="4GB"
            WALLTIME="00:05:00" # Actually took a few seconds
            echo "$1"$'\t'"_" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
        "orthosets")
            NTHREADS="1"
            MEM="4GB"
            WALLTIME="00:10:00"
            for group in "73strains_alpha_complete"; do
                echo "$1"$'\t'"_"$'\t'"$group" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;
        "annotation")
            NTHREADS="1"
            MEM="4GB"
            WALLTIME="02:30:00"
            for group in "73strains_alpha_complete"; do
                echo "$1"$'\t'"_"$'\t'"$group" > "${ARRAYDIR}/${i}"
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
        -q "$QUEUE" -l "walltime=${WALLTIME},mem=${MEM},nodes=1:ppn=${NTHREADS}" \
        -t "$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
