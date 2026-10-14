# ----------------------------------------------------------------------
# Differential abundance of each comparison (taxa collapsed at
# config[differential_abundance][taxa_level], species level by default)
# ----------------------------------------------------------------------

DADIR = f"{COMPDIR}/02_differential_abundance"
DA_LEVEL = config["differential_abundance"]["taxa_level"]
DA_PREFIX = f"{COMP_PREFIX}-L{DA_LEVEL}"


# Collapsed table restricted to the samples of the comparison
rule da_subset:
    input:
        table=f"{TAXDIR}/{PROJ}-collapse-table-L{DA_LEVEL}-{GROUP}.qza",
        metadata=COMP_METADATA,
    output:
        f"{DADIR}/{DA_PREFIX}-collapse-table-{GROUP}.qza",
    log:
        f"{LOGDIR}/06_comparisons/{{comparison}}/da_subset.log",
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


# Frequency filtering by group of samples: the table is split by the
# groups of the comparison 'filter_column', each part keeps the taxa with
# a frequency >= min_frequency found in >= min_samples samples, and the
# filtered parts are merged back.
rule da_filter_by_group:
    input:
        table=f"{DADIR}/{DA_PREFIX}-collapse-table-{GROUP}.qza",
        metadata=COMP_METADATA,
    output:
        table=f"{DADIR}/{DA_PREFIX}-table-filtered-{GROUP}.qza",
        splitdir=directory(f"{DADIR}/frequency_filtering"),
    log:
        f"{LOGDIR}/06_comparisons/{{comparison}}/da_filter_by_group.log",
    params:
        column=lambda wildcards: comparison(wildcards)["filter_column"],
        groups=lambda wildcards: " ".join(
            comparison(wildcards)["groups"][comparison(wildcards)["filter_column"]]),
        prefix=lambda wildcards: f"{PROJ}-{wildcards.comparison}-{GROUP}",
        min_frequency=config["differential_abundance"]["min_frequency"],
        min_samples=config["differential_abundance"]["min_samples"],
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        mkdir -p {output.splitdir}
        : > {log}
        filtered=()
        for group in {params.groups}; do
            split={output.splitdir}/{params.prefix}-table_split-$group.qza
            filt={output.splitdir}/{params.prefix}-table_filtered-$group.qza
            echo "== {params.column} = $group" >> {log}
            qiime feature-table filter-samples \
                --i-table {input.table} \
                --m-metadata-file {input.metadata} \
                --p-where "[{params.column}]='$group'" \
                --p-allow-empty-table \
                --o-filtered-table $split >> {log} 2>&1
            qiime feature-table filter-features \
                --i-table $split \
                --p-min-frequency {params.min_frequency} \
                --p-min-samples {params.min_samples} \
                --p-allow-empty-table \
                --o-filtered-table $filt >> {log} 2>&1
            filtered+=("$filt")
        done
        qiime feature-table merge \
            --i-tables "${{filtered[@]}}" \
            --o-merged-table {output.table} >> {log} 2>&1
        """


rule da_table_viz:
    input:
        f"{DADIR}/{DA_PREFIX}-table-filtered-{GROUP}.qza",
    output:
        report(
            f"{DADIR}/{DA_PREFIX}-table-filtered-{GROUP}.qzv",
            caption=caption("da_table"),
            category="06 Differential abundance",
            subcategory="{comparison}",
        ),
    log:
        f"{LOGDIR}/06_comparisons/{{comparison}}/da_table_viz.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime metadata tabulate --m-input-file {input} --o-visualization {output} > {log} 2>&1
        """


rule da_table_tsv:
    input:
        f"{DADIR}/{DA_PREFIX}-table-filtered-{GROUP}.qza",
    output:
        biom=f"{DADIR}/{DA_PREFIX}-table-filtered-{GROUP}.biom",
        tsv=f"{DADIR}/{DA_PREFIX}-table-filtered-{GROUP}.tsv",
    log:
        f"{LOGDIR}/06_comparisons/{{comparison}}/da_table_tsv.log",
    params:
        exportdir=f"{DADIR}/table_export",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime tools export --input-path {input} --output-path {params.exportdir} > {log} 2>&1
        mv {params.exportdir}/feature-table.biom {output.biom}
        rm -rf {params.exportdir}
        biom convert -i {output.biom} -o {output.tsv} --to-tsv >> {log} 2>&1
        """


# Table formatted for LEfSe (sample IDs row, class row, taxa with '|')
rule lefse:
    input:
        tsv=f"{DADIR}/{DA_PREFIX}-table-filtered-{GROUP}.tsv",
        metadata=COMP_METADATA,
        script=script("format_lefse_table.py"),
    output:
        report(
            f"{DADIR}/{DA_PREFIX}-lefse-{{column}}.tsv",
            caption=caption("lefse"),
            category="06 Differential abundance",
            subcategory="{comparison}",
        ),
    log:
        f"{LOGDIR}/06_comparisons/{{comparison}}/lefse_{{column}}.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        python {input.script} \
            --table-tsv {input.tsv} \
            --metadata {input.metadata} \
            --group-column {wildcards.column} \
            --output {output} > {log} 2>&1
        """


rule pseudocount:
    input:
        f"{DADIR}/{DA_PREFIX}-table-filtered-{GROUP}.qza",
    output:
        f"{DADIR}/{DA_PREFIX}-table-comp-{GROUP}.qza",
    log:
        f"{LOGDIR}/06_comparisons/{{comparison}}/pseudocount.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime composition add-pseudocount \
            --i-table {input} \
            --o-composition-table {output} > {log} 2>&1
        """


rule ancom:
    input:
        table=f"{DADIR}/{DA_PREFIX}-table-comp-{GROUP}.qza",
        metadata=COMP_METADATA,
    output:
        report(
            f"{DADIR}/{DA_PREFIX}-ancom-{{column}}-{GROUP}.qzv",
            caption=caption("ancom"),
            category="06 Differential abundance",
            subcategory="{comparison}",
        ),
    log:
        f"{LOGDIR}/06_comparisons/{{comparison}}/ancom_{{column}}.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime composition ancom \
            --i-table {input.table} \
            --m-metadata-file {input.metadata} \
            --m-metadata-column {wildcards.column} \
            --o-visualization {output} > {log} 2>&1
        """
