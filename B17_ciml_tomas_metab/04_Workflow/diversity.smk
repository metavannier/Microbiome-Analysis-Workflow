# ----------------------------------------------------------------------
# Diversity analyses of each comparison
# ----------------------------------------------------------------------

DIVDIR = f"{COMPDIR}/01_diversity"
CORE = f"{DIVDIR}/core-metrics-results"

ALPHA_METRICS = config["diversity"]["alpha_metrics"]
BETA_METRICS = config["diversity"]["beta_metrics"]
# Outputs of core-metrics-phylogenetic, in the order of the command line
CORE_ALPHA = ["faith_pd", "observed_features", "shannon", "evenness"]
CORE_BETA = ["unweighted_unifrac", "weighted_unifrac", "jaccard", "bray_curtis"]


# Feature table restricted to the samples of the comparison
rule comparison_table:
    input:
        table=DIVERSITY_TABLE,
        metadata=COMP_METADATA,
    output:
        f"{DIVDIR}/{COMP_PREFIX}-table-{GROUP}.qza",
    log:
        f"{LOGDIR}/06_comparisons/{{comparison}}/table.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime feature-table filter-samples \
            --i-table {input.table} \
            --m-metadata-file {input.metadata} \
            --p-filter-empty-features \
            --o-filtered-table {output} > {log} 2>&1
        """


rule comparison_table_summary:
    input:
        table=f"{DIVDIR}/{COMP_PREFIX}-table-{GROUP}.qza",
        metadata=COMP_METADATA,
    output:
        summary=report(
            f"{DIVDIR}/{COMP_PREFIX}-table-{GROUP}.qzv",
            caption=caption("featuretable"),
            category="05 Diversity",
            subcategory="{comparison}",
        ),
        feature_frequencies=f"{DIVDIR}/{COMP_PREFIX}-table-{GROUP}-feature-frequencies.qza",
        sample_frequencies=f"{DIVDIR}/{COMP_PREFIX}-table-{GROUP}-sample-frequencies.qza",
    log:
        f"{LOGDIR}/06_comparisons/{{comparison}}/table_summary.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime feature-table summarize \
            --i-table {input.table} \
            --m-metadata-file {input.metadata} \
            --o-feature-frequencies {output.feature_frequencies} \
            --o-sample-frequencies {output.sample_frequencies} \
            --o-summary {output.summary} \
            --no-recycle > {log} 2>&1
        """


# The maximum depth of the curves is config maxdepth, lowered to the
# largest sample frequency of the comparison when needed (QIIME 2 refuses
# a max depth greater than the frequency of every sample).
rule alpha_rarefaction:
    input:
        table=f"{DIVDIR}/{COMP_PREFIX}-table-{GROUP}.qza",
        tree=TREE,
        metadata=COMP_METADATA,
    output:
        report(
            f"{DIVDIR}/{COMP_PREFIX}-alpha_rarefaction_curves.qzv",
            caption=caption("alphararefaction"),
            category="05 Diversity",
            subcategory="{comparison}",
        ),
    log:
        f"{LOGDIR}/06_comparisons/{{comparison}}/alpha_rarefaction.log",
    params:
        mindepth=config["diversity"]["mindepth"],
        maxdepth=config["diversity"]["maxdepth"],
        exportdir=f"{DIVDIR}/alpha_rarefaction_export",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime tools export --input-path {input.table} --output-path {params.exportdir} > {log} 2>&1
        maxfreq=$(biom summarize-table -i {params.exportdir}/feature-table.biom \
            | awk '$1 == "Max:" {{gsub(",", "", $2); print int($2)}}')
        rm -rf {params.exportdir}
        depth={params.maxdepth}
        if [ "$maxfreq" -lt "$depth" ]; then depth=$maxfreq; fi
        echo "Max depth of the rarefaction curves: $depth (largest sample: $maxfreq)" >> {log}
        qiime diversity alpha-rarefaction \
            --i-table {input.table} \
            --i-phylogeny {input.tree} \
            --m-metadata-file {input.metadata} \
            --p-min-depth {params.mindepth} \
            --p-max-depth $depth \
            --o-visualization {output} >> {log} 2>&1
        """


rule core_metrics:
    input:
        table=f"{DIVDIR}/{COMP_PREFIX}-table-{GROUP}.qza",
        tree=TREE,
        metadata=COMP_METADATA,
    output:
        rarefied=f"{CORE}/rarefied_table.qza",
        vectors=expand(f"{CORE}/{{alpha}}_vector.qza", alpha=CORE_ALPHA, allow_missing=True),
        matrices=expand(f"{CORE}/{{beta}}_distance_matrix.qza", beta=CORE_BETA, allow_missing=True),
        pcoa=expand(f"{CORE}/{{beta}}_pcoa_results.qza", beta=CORE_BETA, allow_missing=True),
        emperor=report(
            expand(f"{CORE}/{{beta}}_emperor.qzv", beta=CORE_BETA, allow_missing=True),
            caption=caption("emperor"),
            category="05 Diversity",
            subcategory="{comparison}",
        ),
    log:
        f"{LOGDIR}/06_comparisons/{{comparison}}/core_metrics.log",
    params:
        depth=config["diversity"]["samplingdepth"],
        seed=("" if config["diversity"].get("random_seed") is None
              else f"--p-random-seed {config['diversity']['random_seed']}"),
    threads: config["diversity"]["threads"]
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime diversity core-metrics-phylogenetic \
            --i-table {input.table} \
            --i-phylogeny {input.tree} \
            --m-metadata-file {input.metadata} \
            --p-sampling-depth {params.depth} {params.seed} \
            --p-n-jobs-or-threads {threads} \
            --o-rarefied-table {output.rarefied} \
            --o-faith-pd-vector {output.vectors[0]} \
            --o-observed-features-vector {output.vectors[1]} \
            --o-shannon-vector {output.vectors[2]} \
            --o-evenness-vector {output.vectors[3]} \
            --o-unweighted-unifrac-distance-matrix {output.matrices[0]} \
            --o-weighted-unifrac-distance-matrix {output.matrices[1]} \
            --o-jaccard-distance-matrix {output.matrices[2]} \
            --o-bray-curtis-distance-matrix {output.matrices[3]} \
            --o-unweighted-unifrac-pcoa-results {output.pcoa[0]} \
            --o-weighted-unifrac-pcoa-results {output.pcoa[1]} \
            --o-jaccard-pcoa-results {output.pcoa[2]} \
            --o-bray-curtis-pcoa-results {output.pcoa[3]} \
            --o-unweighted-unifrac-emperor {output.emperor[0]} \
            --o-weighted-unifrac-emperor {output.emperor[1]} \
            --o-jaccard-emperor {output.emperor[2]} \
            --o-bray-curtis-emperor {output.emperor[3]} \
            --no-recycle --verbose > {log} 2>&1
        """


# Alpha diversity: Kruskal-Wallis between the groups of every categorical
# column of the comparison metadata
rule alpha_significance:
    input:
        vector=f"{CORE}/{{alpha}}_vector.qza",
        metadata=COMP_METADATA,
    output:
        report(
            f"{CORE}/{{alpha}}-group-significance.qzv",
            caption=caption("alphasignificance"),
            category="05 Diversity",
            subcategory="{comparison}",
        ),
    log:
        f"{LOGDIR}/06_comparisons/{{comparison}}/{{alpha}}_significance.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime diversity alpha-group-significance \
            --i-alpha-diversity {input.vector} \
            --m-metadata-file {input.metadata} \
            --o-visualization {output} > {log} 2>&1
        """


# Beta diversity: PERMANOVA (global and pairwise) for each compared column
rule beta_significance:
    input:
        matrix=f"{CORE}/{{beta}}_distance_matrix.qza",
        metadata=COMP_METADATA,
    output:
        report(
            f"{CORE}/{{column}}-{{beta}}-significance.qzv",
            caption=caption("beta_significance"),
            category="05 Diversity",
            subcategory="{comparison}",
        ),
    log:
        f"{LOGDIR}/06_comparisons/{{comparison}}/{{column}}-{{beta}}_significance.log",
    params:
        permutations=config["diversity"]["permutations"],
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime diversity beta-group-significance \
            --i-distance-matrix {input.matrix} \
            --m-metadata-file {input.metadata} \
            --m-metadata-column {wildcards.column} \
            --p-method permanova \
            --p-permutations {params.permutations} \
            --p-pairwise \
            --o-visualization {output} > {log} 2>&1
        """


# Taxa bar plot of the samples of the comparison
rule comparison_barplot:
    input:
        table=f"{TAXDIR}/{PROJ}-rarefaction-table-filtered-{GROUP}.qza",
        taxonomy=f"{TAXDIR}/{PROJ}-tax_sklearn.qza",
        metadata=COMP_METADATA,
    output:
        table=temp(f"{DIVDIR}/{COMP_PREFIX}-rarefaction-table-{GROUP}.qza"),
        barplot=report(
            f"{DIVDIR}/{COMP_PREFIX}-taxa-bar-plots.qzv",
            caption=caption("taxbarplot"),
            category="05 Diversity",
            subcategory="{comparison}",
        ),
    log:
        f"{LOGDIR}/06_comparisons/{{comparison}}/taxa_barplot.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        (qiime feature-table filter-samples \
            --i-table {input.table} \
            --m-metadata-file {input.metadata} \
            --p-filter-empty-features \
            --o-filtered-table {output.table}
        qiime taxa barplot \
            --i-table {output.table} \
            --i-taxonomy {input.taxonomy} \
            --m-metadata-file {input.metadata} \
            --o-visualization {output.barplot}) > {log} 2>&1
        """
