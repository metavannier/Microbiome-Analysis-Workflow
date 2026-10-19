# ----------------------------------------------------------------------
# Final files of the common part and of each comparison
# ----------------------------------------------------------------------

def core_targets():
    """Common part, run once for all the comparisons."""
    targets = [
        # Quality control (FastQC/MultiQC before and after Trimmomatic)
        f"{QCDIR}/raw_multiqc.html",
        f"{QCDIR}/trimmed_multiqc.html",
        # QIIME 2 import and primer removal
        f"{IMPORTDIR}/{PROJ}-demux-paired-end.qzv",
        f"{IMPORTDIR}/{PROJ}-PE-demux-noprimer.qzv",
        f"{IMPORTDIR}/{PROJ}-cutadapt-stats.qzv",
        # DADA2
        f"{DADA2DIR}/{PROJ}-{GROUP}-stats.qzv",
        f"{DADA2DIR}/{PROJ}-{GROUP}-base-transitions.qzv",
        # Taxonomy
        f"{TAXDIR}/{PROJ}-rep-filtered-seqs-taxa-{GROUP}.qzv",
        f"{TAXDIR}/{PROJ}-taxa-table-filtered-{GROUP}.qzv",
        f"{TAXDIR}/taxonomy_filtered.tsv",
        f"{TAXDIR}/taxonomy_ranks.tsv",
        f"{TAXDIR}/{PROJ}-asv-table-with-taxonomy.tsv",
        f"{TAXDIR}/{PROJ}-taxa-bar-plots.qzv",
        f"{TAXDIR}/{PROJ}-level-{MA_LEVEL}.csv",
        # Phylogeny: table of the diversity analyses
        f"{PHYLODIR}/{PROJ}-taxa-table-filtered-tree-{PHYLO}-{GROUP}.qzv",
    ]
    if CHECK_CLASSIFIER:
        # Genera cross-checked with the uniform classifier
        targets.append(f"{TAXDIR}/taxonomy_genus_check.tsv")
    return targets


def comparison_targets(cid):
    """All the results of one comparison."""
    columns = COMPARISONS[cid]["group_columns"]
    patterns = [
        COMP_METADATA,
        # Diversity
        f"{DIVDIR}/{COMP_PREFIX}-table-{GROUP}.qzv",
        f"{DIVDIR}/{COMP_PREFIX}-alpha_rarefaction_curves.qzv",
        f"{DIVDIR}/{COMP_PREFIX}-taxa-bar-plots.qzv",
        *expand(f"{CORE}/{{beta}}_emperor.qzv", beta=CORE_BETA,
                allow_missing=True),
        *expand(f"{CORE}/{{alpha}}-group-significance.qzv",
                alpha=ALPHA_METRICS, allow_missing=True),
        *expand(f"{CORE}/{{column}}-{{beta}}-significance.qzv",
                column=columns, beta=BETA_METRICS, allow_missing=True),
        # Differential abundance
        f"{DADIR}/{DA_PREFIX}-table-filtered-{GROUP}.qzv",
        *expand(f"{DADIR}/{DA_PREFIX}-ancom-{{column}}-{GROUP}.qzv",
                column=columns, allow_missing=True),
        *expand(f"{DADIR}/{DA_PREFIX}-lefse-{{column}}.tsv",
                column=columns, allow_missing=True),
        # MetaboAnalyst (metaboanalyst: taxa_level, 8 = species)
        *expand(f"{MADIR}/{COMP_PREFIX}-level-{MA_LEVEL}-{{column}}-metaboanalyst.csv",
                column=columns, allow_missing=True),
    ]
    return expand(patterns, comparison=cid)


def all_targets():
    targets = core_targets()
    for cid in selected_comparisons():
        targets.extend(comparison_targets(cid))
    return targets
