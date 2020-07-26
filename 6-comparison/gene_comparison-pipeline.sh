#!/bin/bash
#
# Compare divergence summary stats for symbiosis and background genes.
#
# UPDATE 2020-07-06: Changed data merging code to fix missing columns.
#
# UPDATE 2020-07-17: Added PAV-only delta calculation. Also ran a
#    new calculation for Ka/Ks values that only includes those with
#    3 or more comparisons. In addition, stopped sorting the output,
#    because the sorting messes up correlation calculations. And,
#    ran the delta comparison on genes with at least 3 sequences only.
#
# SETTINGS
#
QUEUE="small"
WALLTIME="06:00:00"
NTHREADS=1
MEM="8GB"
#
N_RANDOM_DRAWS=1000
COPY_NUMBER_TOLERANCE=0.2

flip="""
import sys
data = []
with open(sys.argv[1], 'rt') as h:
    for line in h:
        row = line.split('\t')
        for i, field in enumerate(row):
            f = field.strip()
            if i + 1 > len(data):
                data.append([])
            data[i].append(f)
for d in data:
    sys.stdout.write('\t'.join(d) + '\n')
"""

PROJDIR="/home/tiffinp/epste051/project/symbiosis_gene_evolution"
STRAINS="${PROJDIR}/data/strain_lists/73strains_complete_alpha"
ORTHODIR="${PROJDIR}/results/ortho/74strains/orthofinder/2019-09-24"
PSEUDOGENES="${ORTHODIR}/pseudogenes.tsv"
ORTHOTABLE="${ORTHODIR}/Results_Directory/reconciliation_Sep25/orthosets.73strains_alpha_complete.tsv"
DIV_DIR="${PROJDIR}/results/gene_tables/divergence/73strains_alpha_complete/2019-09-26"
DIV_TABLE="${DIV_DIR}/gene_distances.nopseudo.orthosets.tsv"
ALL_PAIRWISE="${DIV_DIR}/gene_distances.nopseudo.orthosets_pairwise_sequence.tsv"
SPECIES="${PROJDIR}/notes/table/73strains_alpha_complete.species.tsv"
ORTHOSETS="${PROJDIR}/results/ortho/74strains/orthofinder/2019-09-24/Results_Directory/reconciliation_Sep25/orthosets.73strains_alpha_complete.tsv"
PHYLO_SIGNAL_DIR="${PROJDIR}/results/phylo/73strains_alpha_complete/phylo_signal/2020-05-03"
PHYLO_SIGNAL_DIR_PAONLY="${PROJDIR}/results/phylo/73strains_alpha_complete/phylo_signal/2020-07-16_pa_only"
DELTA_FP="${PHYLO_SIGNAL_DIR}/signal.tsv"
DELTA_FP_PAONLY="${PHYLO_SIGNAL_DIR_PAONLY}/signal.tsv"
R2_SCC="${PHYLO_SIGNAL_DIR}/r2_scc_genes.tsv"
DTL="${PROJDIR}/results/reconciliation/73strains_alpha_complete/generax/2020-04-23/dtl.tsv"

declare -A CANDIDATE_FILES
CANDIDATE_FILES[symbiosis]="${PROJDIR}/notes/table/symbiosis_genes.edited.2020-05-07.txt"
CANDIDATE_FILES[signaling]="${PROJDIR}/notes/table/symbiosis_signaling_genes.edited.2020-05-07.txt"
CANDIDATE_FILES[fixation]="${PROJDIR}/notes/table/symbiosis_fixation_genes.edited.2020-05-07.txt"
CANDIDATE_FILES[biomass_A17]="${PROJDIR}/notes/table/gwas_candidates/2020-05-11/plant_biomass.A17.top10.txt"
CANDIDATE_FILES[biomass_R108]="${PROJDIR}/notes/table/gwas_candidates/2020-05-11/plant_biomass.R108.top10.txt"
CANDIDATE_FILES[rhizobial_fitness_A17]="${PROJDIR}/notes/table/gwas_candidates/2020-05-11/rhizobial_fitness.A17.top10.txt"
CANDIDATE_FILES[rhizobial_fitness_R108]="${PROJDIR}/notes/table/gwas_candidates/2020-05-11/rhizobial_fitness.R108.top10.txt"

OUTDIR="${PROJDIR}/results/gene_comparison/73strains_alpha_complete/2020-05-11"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
SCRIPTDIR="${OUTDIR}/script_copies"
ARRAYDIR="${OUTDIR}/arrayjobdata"

BIN_DIR="${PROJDIR}/script/bin"
compare="${BIN_DIR}/compare_divergence_orthosets2.py"
table="${BIN_DIR}/tabulate_pairwise_divergence.py"
randomizer="${BIN_DIR}/randomize_divergence_orthosets.py"
jaccard="${BIN_DIR}/jaccard_index.py"
decile_plotter="${BIN_DIR}/decile_plots_from_full_table.r"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    set -euo pipefail

    SECONDS=0

    cd "$OUTDIR"

    TASK="$(cut -f 1 "${ARRAYDIR}/${PBS_ARRAYID}")"
    GENESET="$(cut -f 2 "${ARRAYDIR}/${PBS_ARRAYID}")"

    case "$TASK" in


        calculate-divergence)
            cp "$table" .
            "./$(basename "$table")" \
                --strains "$STRAINS" \
                --exclude "$PSEUDOGENES" \
                "$ALL_PAIRWISE" "$ORTHOTABLE" | \
               sort -k 1b,1 -t$'\t' > \
                divergence.tsv \
                || { echo "making gene divergence stats table failed"; exit 1; }

            joinr=''
            joinr+='argv = commandArgs(trailingOnly=TRUE); '
            joinr+='x = read.csv(argv[1], sep="\t", header=TRUE, as.is=TRUE, check.names=FALSE); '
            joinr+='y = read.csv(argv[2], sep="\t", header=TRUE, as.is=TRUE, check.names=FALSE); '
            joinr+='write.table(merge(x, y, all=TRUE, by="gene"), stdout(), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE); '

            sort -k 1b,1 -t$'\t' "$R2_SCC" > tmp \
                || { echo "sorting r2 file failed"; exit 1; }
            sort -k 1b,1 -t$'\t' "$DTL" | \
                awk 'NR == 1 { print $0 "\t" "transfer_minus_dup"}; NR > 1 { print $0 "\t" $4-$2};' \
                > tmp2 \
                || { echo "sorting dtl file failed"; exit 1; }
            sed -e '1i\gene	delta	fitz_purvis	fitz_purvis.p0	fitz_purvis.p1' \
                "$DELTA_FP" \
                | sort -k 1b,1 -t$'\t' \
                > tmp3a
            sed -e '1i\gene	delta_pav	fitz_purvis	fitz_purvis.p0	fitz_purvis.p1' \
                "$DELTA_FP_PAONLY" \
                | sort -k 1b,1 -t$'\t' \
                | cut -f 1,2 \
                > tmp3b

            Rscript -e "$joinr" tmp3a tmp3b > tmp3 \
                || { echo "joining phylogenetic signal datasets failed"; exit 1; }

            Rscript -e "$joinr" divergence.tsv tmp3 \
                > tmp4 \
                || { echo "joining divergence and phylo signal failed"; exit 1; }
            Rscript -e "$joinr" tmp4 tmp \
                > tmp5 \
                || { echo "joining data and r2 failed"; exit 1; }
            Rscript -e "$joinr" tmp5 tmp2 \
                > data.tsv \
                || { echo "joining data & reconciliation failed"; exit 1; }
            rm tmp*
            ;;

        gene-lists)
            mkdir -p genes
            mkdir -p randomized
            mkdir -p randomized/all
            mkdir -p randomized/matched_strain_count
            mkdir -p randomized/matched_count_and_copy
            mkdir -p randomized/phenotype_permutations

            for geneset in "${!CANDIDATE_FILES[@]}"; do
                cp "${CANDIDATE_FILES[$geneset]}" "genes/${geneset}" \
                    || { echo "copying ${geneset} failed"; exit 1; }
            done

            cat genes/biomass_A17 genes/biomass_R108 > \
                genes/biomass \
                || { echo "combining biomass across genotypes failed"; exit 1; }
            cat genes/rhizobial_fitness_A17 genes/rhizobial_fitness_R108 > \
                genes/rhizobial_fitness \
                || { echo "combining fitness across genotypes failed"; exit 1; }
            cat genes/rhizobial_fitness genes/biomass \
                > genes/gwas \
                || { echo "combining gwas across phenotypes failed"; exit 1; }

            for geneset in biomass_A17 biomass_R108 rhizobial_fitness_A17 rhizobial_fitness_R108; do
                fname="$(echo "${CANDIDATE_FILES[$geneset]}" | sed 's/top10/top10random/g')"
                python -c "$flip" "$fname" \
                    > "randomized/phenotype_permutations/${geneset}" \
                    || { echo "copying phenotype permutations for ${geneset} failed"; exit 1; }
            done

            paste -d$'\t' randomized/phenotype_permutations/biomass_A17 \
                randomized/phenotype_permutations/biomass_R108 > \
                randomized/phenotype_permutations/biomass \
                || { echo "combining biomass random across genotypes failed"; exit 1; }
            paste -d$'\t' randomized/phenotype_permutations/rhizobial_fitness_A17 \
                randomized/phenotype_permutations/rhizobial_fitness_R108 > \
                randomized/phenotype_permutations/rhizobial_fitness \
                || { echo "combining fitness random across genotypes failed"; exit 1; }
            paste -d$'\t' randomized/phenotype_permutations/rhizobial_fitness \
                randomized/phenotype_permutations/biomass\
                > randomized/phenotype_permutations/gwas
            ;;

        randomize)
            ## Pull out random samples of the same size as the real sample
            ## from all genes

            ## Pull out random samples that are matched in characteristics
            ## 1 -- number of strains (for copy number comparisons)
            ## 2 -- number of strains and number of copies 
            ## Note that earlier versions of this code also tried to match
            ## on phylogenetic distance between strains, but that was
            ## dropped.
            cp "$randomizer" .
            for genes_file in genes/*; do
                "./$(basename "$randomizer")" \
                    --copy-tol "$COPY_NUMBER_TOLERANCE" \
                    --output-all "randomized/all/$(basename "$genes_file")" \
                    --output-strain-count "randomized/matched_strain_count/$(basename "$genes_file")" \
                    --output-copy-count "randomized/matched_count_and_copy/$(basename "$genes_file")" \
                    divergence.tsv "$genes_file" "$N_RANDOM_DRAWS" \
                    || { echo "making randomization for ${genes_file} failed"; exit 1; }
            done
            ;;

        compare-divergence)
            ## Using the pre-calculated randomized datasets, compare
            ## strain count, mean copy number, pairwise ka/ks, and pairwise
            ## protein divergence. Also calculate the Jaccard index
            ## between pairs of strains based on gene sharing.

            mkdir -p "results/${GENESET}" && cd "results/${GENESET}"
            cp "$compare" .
            cp "$jaccard" .

            rand="all"
            rand_file="../../randomized/all/${GENESET}"
            "./$(basename "$compare")" \
                --output "table.${rand}" \
                --targets "../../genes/${GENESET}" \
                --strains "$STRAINS" \
                --strain-count-rand "$rand_file" \
                --copy-number-rand "$rand_file" \
                --ngenes-rand "$rand_file" \
                --kaks-rand "$rand_file" \
                --kaks-rand-3 "$rand_file" \
                --protein-divergence-rand "$rand_file" \
                --relative-protein-divergence-rand "$rand_file" \
                --delta-rand "$rand_file" \
                --delta-pav-rand-3 "$rand_file" \
                --delta-pav-rand "$rand_file" \
                --fitz-purvis-rand "$rand_file" \
                --r2-scc-rand "$rand_file" \
                --dtl-rand "$rand_file" \
                ../../data.tsv \
                || { echo "${GENESET}, ${rand} failed"; exit 1; }

            rand="matched"
            rand_file1="../../randomized/matched_strain_count/${GENESET}"
            rand_file2="../../randomized/matched_count_and_copy/${GENESET}"
            rand_file2_kaks3="../../randomized/matched_count_and_copy_kaks3/${GENESET}"
            "./$(basename "$compare")" \
                --output "table.${rand}" \
                --targets "../../genes/${GENESET}" \
                --strains "$STRAINS" \
                --copy-number-rand "$rand_file1" \
                --ngenes-rand "$rand_file1" \
                --kaks-rand "$rand_file2" \
                --kaks-rand-3 "$rand_file2" \
                --protein-divergence-rand "$rand_file2" \
                --relative-protein-divergence-rand "$rand_file2" \
                --delta-pav-rand-3 "$rand_file2" \
                --delta-pav-rand "$rand_file2" \
                --fitz-purvis-rand "$rand_file2" \
                --r2-scc-rand "$rand_file2" \
                --dtl-rand "$rand_file2" \
                ../../data.tsv \
                || { echo "${GENESET}, ${rand} failed"; exit 1; }

            if [[ -s "../../randomized/phenotype_permutations/${GENESET}" ]]; then
                rand="phenotype_permutations"
                rand_file="../../randomized/phenotype_permutations/${GENESET}"
                "./$(basename "$compare")" \
                    --output "table.${rand}" \
                    --targets "../../genes/${GENESET}" \
                    --strains "$STRAINS" \
                    --strain-count-rand "$rand_file" \
                    --copy-number-rand "$rand_file" \
                    --ngenes-rand "$rand_file" \
                    --kaks-rand "$rand_file" \
                    --protein-divergence-rand "$rand_file" \
                    --relative-protein-divergence-rand "$rand_file" \
                    --delta-rand "$rand_file" \
                    --delta-pav-rand "$rand_file" \
                    --fitz-purvis-rand "$rand_file" \
                    --r2-scc-rand "$rand_file" \
                    --dtl-rand "$rand_file" \
                    ../../data.tsv \
                    || { echo "${GENESET}, ${rand} failed"; exit 1; }
            fi

            "./$(basename "$jaccard")" "$DIV_TABLE" \
                "strains" "../../genes/${GENESET}" "subset" > \
                "jaccard_index.txt" \
                || { echo "calculating jaccard index without pseudogenes for ${GENSET} failed"; exit 1; }
            ;;

        plots)
            mkdir -p plots && cd plots
            mkdir -p "$GENESET" && cd "$GENESET"
            cp "$decile_plotter" .
            "./$(basename "$decile_plotter")" \
                decile_plots.pdf \
                "../../results/${GENESET}/table."* \
                || { echo "making plots for ${GENESET} failed"; exit 1; }

            if [[ -s "../../randomized/phenotype_permutations/${GENESET}" ]]; then
                "./$(basename "$decile_plotter")" \
                    key_decile_plots.pdf \
                    "../../results/${GENESET}/table.all.strain_count.tsv" \
                    "../../results/${GENESET}/table.all.copy_number.tsv" \
                    "../../results/${GENESET}/table.matched.copy_number.tsv" \
                    "../../results/${GENESET}/table.phenotype_permutations.copy_number.tsv" \
                    "../../results/${GENESET}/table.all.log2_median_kaks.tsv" \
                    "../../results/${GENESET}/table.matched.log2_median_kaks.tsv" \
                    "../../results/${GENESET}/table.phenotype_permutations.log2_median_kaks.tsv" \
                    "../../results/${GENESET}/table.all.r2_scc.tsv" \
                    "../../results/${GENESET}/table.matched.r2_scc.tsv" \
                    "../../results/${GENESET}/table.phenotype_permutations.r2_scc.tsv" \
                    "../../results/${GENESET}/table.all.delta.tsv" \
                    "../../results/${GENESET}/table.matched.delta.tsv" \
                    "../../results/${GENESET}/table.phenotype_permutations.delta.tsv" \
                    "../../results/${GENESET}/table.all.fitz_purvis.tsv" \
                    "../../results/${GENESET}/table.matched.fitz_purvis.tsv" \
                    "../../results/${GENESET}/table.phenotype_permutations.fitz_purvis.tsv" \
                    || { echo "making plots for ${GENESET} failed"; exit 1; }

                "./$(basename "$decile_plotter")" \
                    dup_transfer_decile_plots.pdf \
                    "../../results/${GENESET}/table.all.dup_rate.tsv" \
                    "../../results/${GENESET}/table.all.transfer_minus_dup_rate.tsv" \
                    "../../results/${GENESET}/table.all.transfer_rate.tsv" \
                    "../../results/${GENESET}/table.matched.dup_rate.tsv" \
                    "../../results/${GENESET}/table.matched.transfer_minus_dup_rate.tsv" \
                    "../../results/${GENESET}/table.matched.transfer_rate.tsv" \
                    "../../results/${GENESET}/table.phenotype_permutations.dup_rate.tsv" \
                    "../../results/${GENESET}/table.phenotype_permutations.transfer_minus_dup_rate.tsv" \
                    "../../results/${GENESET}/table.phenotype_permutations.transfer_rate.tsv" \
                    || { echo "making plots for ${GENESET} failed"; exit 1; }
            else
                "./$(basename "$decile_plotter")" \
                    key_decile_plots.pdf \
                    "../../results/${GENESET}/table.all.strain_count.tsv" \
                    "../../results/${GENESET}/table.all.copy_number.tsv" \
                    "../../results/${GENESET}/table.matched.copy_number.tsv" \
                    "../../results/${GENESET}/table.all.log2_median_kaks.tsv" \
                    "../../results/${GENESET}/table.matched.log2_median_kaks.tsv" \
                    "../../results/${GENESET}/table.all.r2_scc.tsv" \
                    "../../results/${GENESET}/table.matched.r2_scc.tsv" \
                    "../../results/${GENESET}/table.all.delta.tsv" \
                    "../../results/${GENESET}/table.matched.delta.tsv" \
                    "../../results/${GENESET}/table.all.fitz_purvis.tsv" \
                    "../../results/${GENESET}/table.matched.fitz_purvis.tsv" \
                    || { echo "making plots for ${GENESET} failed"; exit 1; }

                "./$(basename "$decile_plotter")" \
                    dup_transfer_decile_plots.pdf \
                    "../../results/${GENESET}/table.all.dup_rate.tsv" \
                    "../../results/${GENESET}/table.all.transfer_minus_dup_rate.tsv" \
                    "../../results/${GENESET}/table.all.transfer_rate.tsv" \
                    "../../results/${GENESET}/table.matched.dup_rate.tsv" \
                    "../../results/${GENESET}/table.matched.transfer_minus_dup_rate.tsv" \
                    "../../results/${GENESET}/table.matched.transfer_rate.tsv" \
                    || { echo "making plots for ${GENESET} failed"; exit 1; }
            fi
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
        compare-divergence|plots)
            for geneset in "${OUTDIR}/genes/"*; do
                echo "$1"$'\t'"$(basename "$geneset")" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;
        *)
            echo "$1" > "${ARRAYDIR}/${i}"
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
        -q "$QUEUE" -l "walltime=${WALLTIME},nodes=1:ppn=${NTHREADS},mem=${MEM}" \
        -t "$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
