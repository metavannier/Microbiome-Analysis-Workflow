# ----------------------------------------------------------------------
# Level 7 (deepest rank) taxon abundance tables for MetaboAnalyst
# ----------------------------------------------------------------------
# Source: the level-N.csv of the project taxa bar plot, i.e. the file of
# "Download CSV" in QIIME 2 View (rule barplot_csv). For each compared
# column, the samples of the comparison are kept, counts are converted to
# per-sample relative abundances and the Group column is filled from the
# comparison metadata.

MA = config["metaboanalyst"]
MA_LEVEL = MA["taxa_level"]
MADIR = f"{COMPDIR}/03_metaboanalyst"


rule metaboanalyst:
    input:
        csv=f"{TAXDIR}/{PROJ}-level-{MA_LEVEL}.csv",
        metadata=COMP_METADATA,
        script=script("format_taxa_table.py"),
    output:
        report(
            f"{MADIR}/{COMP_PREFIX}-level-{MA_LEVEL}-{{column}}-metaboanalyst.csv",
            caption=caption("metaboanalyst"),
            category="07 MetaboAnalyst",
            subcategory="{comparison}",
        ),
    log:
        f"{LOGDIR}/06_comparisons/{{comparison}}/metaboanalyst_{{column}}.log",
    params:
        sep=MA["sep"],
        decimal=MA["decimal"],
        min_group_prevalence=MA["min_group_prevalence"],
        counts="" if MA["relative_abundance"] else "--counts",
        zero_taxa="" if MA["drop_zero_taxa"] else "--keep-zero-taxa",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        python {input.script} \
            --taxa-csv {input.csv} \
            --metadata {input.metadata} \
            --group-column {wildcards.column} \
            --sep '{params.sep}' --decimal '{params.decimal}' \
            --min-group-prevalence {params.min_group_prevalence} \
            {params.counts} {params.zero_taxa} \
            --output {output} > {log} 2>&1
        """
