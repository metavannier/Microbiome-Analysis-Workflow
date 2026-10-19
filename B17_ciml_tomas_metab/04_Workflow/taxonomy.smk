# ----------------------------------------------------------------------
# Assigning taxonomy (naive Bayes classifier, SILVA)
# ----------------------------------------------------------------------

# Optional control classifier ("" = no genus check)
CHECK_CLASSIFIER = config["taxonomy"].get("check_classifier", "")


rule assign_tax:
    input:
        rep=f"{DADA2DIR}/{PROJ}-rep-filtered-seqs-{GROUP}.qza",
        classifier=config["taxonomy"]["classifier"],
    output:
        f"{TAXDIR}/{PROJ}-tax_sklearn.qza",
    log:
        f"{LOGDIR}/04_taxonomy/{PROJ}_sklearn_q2.log",
    params:
        confidence=config["taxonomy"]["confidence"],
        reads_per_batch=config["taxonomy"].get("reads_per_batch", "auto"),
    threads: config["taxonomy"]["njobs"]
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime feature-classifier classify-sklearn \
            --i-classifier {input.classifier} \
            --i-reads {input.rep} \
            --p-confidence {params.confidence} \
            --p-n-jobs {threads} \
            --p-reads-per-batch {params.reads_per_batch} \
            --o-classification {output} \
            --verbose > {log} 2>&1
        """


# Same classification with the control classifier (taxonomy:
# check_classifier, the SILVA uniform classifier by default)
use rule assign_tax as assign_tax_check with:
    input:
        rep=f"{DADA2DIR}/{PROJ}-rep-filtered-seqs-{GROUP}.qza",
        classifier=CHECK_CLASSIFIER,
    output:
        f"{TAXDIR}/{PROJ}-tax_sklearn-check.qza",
    log:
        f"{LOGDIR}/04_taxonomy/{PROJ}_sklearn-check_q2.log",


# Genera of the weighted (habitat-specific) classifier compared with the
# uniform classifier, as recommended by SILVA: the "different" ASVs are
# listed first
rule taxonomy_genus_check:
    input:
        taxonomy=f"{TAXDIR}/taxonomy.tsv",
        check=f"{TAXDIR}/{PROJ}-tax_sklearn-check.qza",
        script=script("compare_taxonomy.py"),
    output:
        report(
            f"{TAXDIR}/taxonomy_genus_check.tsv",
            caption=caption("taxonomy_check"),
            category="04 Taxonomy",
            labels={"result": "genus check (uniform classifier)"},
        ),
    log:
        f"{LOGDIR}/04_taxonomy/{PROJ}_taxonomy_genus_check.log",
    params:
        exportdir=f"{TAXDIR}/taxonomy_check_export",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        (rm -rf {params.exportdir}
        qiime tools export --input-path {input.check} --output-path {params.exportdir}
        python {input.script} --taxonomy {input.taxonomy} \
            --check {params.exportdir}/taxonomy.tsv --output {output}
        rm -rf {params.exportdir}) > {log} 2>&1
        """


# Remove the ASVs assigned to unwanted taxa (default: d__Eukaryota)
rule taxa_filter_table:
    input:
        table=f"{DADA2DIR}/{PROJ}-table-filtered-{GROUP}.qza",
        taxonomy=f"{TAXDIR}/{PROJ}-tax_sklearn.qza",
    output:
        f"{TAXDIR}/{PROJ}-taxa-table-filtered-{GROUP}.qza",
    log:
        f"{LOGDIR}/04_taxonomy/{PROJ}_taxa_filter_table.log",
    params:
        exclude=config["taxonomy"]["exclude"],
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime taxa filter-table \
            --i-table {input.table} \
            --i-taxonomy {input.taxonomy} \
            --p-exclude '{params.exclude}' \
            --o-filtered-table {output} > {log} 2>&1
        """


rule filter_seq:
    input:
        rep=f"{DADA2DIR}/{PROJ}-rep-filtered-seqs-{GROUP}.qza",
        table=f"{TAXDIR}/{PROJ}-taxa-table-filtered-{GROUP}.qza",
    output:
        f"{TAXDIR}/{PROJ}-rep-filtered-seqs-taxa-{GROUP}.qza",
    log:
        f"{LOGDIR}/04_taxonomy/{PROJ}_filter_seq.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime feature-table filter-seqs \
            --i-data {input.rep} \
            --i-table {input.table} \
            --o-filtered-data {output} > {log} 2>&1
        """


# Taxonomy table (all ASVs) and the same table without the excluded taxa
rule gen_tax:
    input:
        f"{TAXDIR}/{PROJ}-tax_sklearn.qza",
    output:
        table_tax=f"{TAXDIR}/taxonomy.tsv",
        table_tax_filtered=report(
            f"{TAXDIR}/taxonomy_filtered.tsv",
            caption=caption("tax"),
            category="04 Taxonomy",
            labels={"result": "taxonomy of the ASVs"},
        ),
    log:
        f"{LOGDIR}/04_taxonomy/{PROJ}_export_taxonomy.log",
    params:
        exportdir=f"{TAXDIR}/taxonomy_export",
        exclude=config["taxonomy"]["exclude"].replace(",", "|"),
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime tools export --input-path {input} --output-path {params.exportdir} > {log} 2>&1
        mv {params.exportdir}/taxonomy.tsv {output.table_tax}
        rm -rf {params.exportdir}
        grep -Eiv '{params.exclude}' {output.table_tax} > {output.table_tax_filtered}
        """


# Rank held by each taxonomic level (SILVA 144: level 7 = genus, level 8 =
# species with the weighted classifiers; SILVA 138: level 7 = species)
rule taxonomy_ranks:
    input:
        taxonomy=f"{TAXDIR}/taxonomy.tsv",
        script=script("taxonomy_ranks.py"),
    output:
        report(
            f"{TAXDIR}/taxonomy_ranks.tsv",
            caption=caption("taxonomy_ranks"),
            category="04 Taxonomy",
            labels={"result": "taxonomic levels"},
        ),
    log:
        f"{LOGDIR}/04_taxonomy/{PROJ}_taxonomy_ranks.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        python {input.script} --taxonomy {input.taxonomy} --output {output} > {log} 2>&1
        """


# ASV sequences with their taxonomy
rule rep_seqs_viz:
    input:
        rep=f"{TAXDIR}/{PROJ}-rep-filtered-seqs-taxa-{GROUP}.qza",
        taxonomy=f"{TAXDIR}/{PROJ}-tax_sklearn.qza",
    output:
        report(
            f"{TAXDIR}/{PROJ}-rep-filtered-seqs-taxa-{GROUP}.qzv",
            caption=caption("dada2seq"),
            category="03 DADA2",
            labels={"result": "ASV sequences"},
        ),
    log:
        f"{LOGDIR}/04_taxonomy/{PROJ}_rep_seqs_viz.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime feature-table tabulate-seqs \
            --i-data {input.rep} \
            --i-taxonomy {input.taxonomy} \
            --p-merge-method intersect \
            --o-visualization {output} > {log} 2>&1
        """


# Frequency per sample: use it to choose the sampling depth (config.yaml)
rule table_summary:
    input:
        table=f"{TAXDIR}/{PROJ}-taxa-table-filtered-{GROUP}.qza",
        metadata=METADATA_FILE,
    output:
        summary=report(
            f"{TAXDIR}/{PROJ}-taxa-table-filtered-{GROUP}.qzv",
            caption=caption("featuretable"),
            category="04 Taxonomy",
            labels={"result": "feature table summary"},
        ),
        feature_frequencies=f"{TAXDIR}/{PROJ}-taxa-table-filtered-{GROUP}-feature-frequencies.qza",
        sample_frequencies=f"{TAXDIR}/{PROJ}-taxa-table-filtered-{GROUP}-sample-frequencies.qza",
    log:
        f"{LOGDIR}/04_taxonomy/{PROJ}_table_summary.log",
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


# Remove the samples with fewer reads than the sampling depth
rule filter_rarefaction:
    input:
        f"{TAXDIR}/{PROJ}-taxa-table-filtered-{GROUP}.qza",
    output:
        f"{TAXDIR}/{PROJ}-rarefaction-table-filtered-{GROUP}.qza",
    log:
        f"{LOGDIR}/04_taxonomy/{PROJ}_filter_rarefaction.log",
    params:
        depth=config["diversity"]["samplingdepth"],
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime feature-table filter-samples \
            --i-table {input} \
            --p-min-frequency {params.depth} \
            --o-filtered-table {output} > {log} 2>&1
        """


# ASV table (counts per sample) with the taxonomy of each ASV
rule asv_table:
    input:
        table=f"{TAXDIR}/{PROJ}-rarefaction-table-filtered-{GROUP}.qza",
        taxonomy=f"{TAXDIR}/taxonomy.tsv",
    output:
        table_biom=f"{TAXDIR}/feature-table.biom",
        taxo_table_biom=f"{TAXDIR}/{PROJ}-asv-table-with-taxonomy.biom",
        taxo_table_tsv=report(
            f"{TAXDIR}/{PROJ}-asv-table-with-taxonomy.tsv",
            caption=caption("asv_table"),
            category="04 Taxonomy",
            labels={"result": "ASV table with taxonomy"},
        ),
    log:
        f"{LOGDIR}/04_taxonomy/{PROJ}_asv_table.log",
    params:
        exportdir=f"{TAXDIR}/table_export",
        biom_taxonomy=f"{TAXDIR}/taxonomy_biom_header.tsv",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime tools export --input-path {input.table} --output-path {params.exportdir} > {log} 2>&1
        mv {params.exportdir}/feature-table.biom {output.table_biom}
        rm -rf {params.exportdir}
        # biom expects the header '#OTUID taxonomy confidence'
        sed '1s/.*/#OTUID\\ttaxonomy\\tconfidence/' {input.taxonomy} > {params.biom_taxonomy}
        biom add-metadata -i {output.table_biom} -o {output.taxo_table_biom} \
            --observation-metadata-fp {params.biom_taxonomy} --sc-separated taxonomy >> {log} 2>&1
        biom convert -i {output.taxo_table_biom} -o {output.taxo_table_tsv} \
            --to-tsv --header-key taxonomy >> {log} 2>&1
        rm -f {params.biom_taxonomy}
        """


rule taxa_barplot:
    input:
        table=f"{TAXDIR}/{PROJ}-rarefaction-table-filtered-{GROUP}.qza",
        taxonomy=f"{TAXDIR}/{PROJ}-tax_sklearn.qza",
        metadata=METADATA_FILE,
    output:
        report(
            f"{TAXDIR}/{PROJ}-taxa-bar-plots.qzv",
            caption=caption("taxbarplot"),
            category="04 Taxonomy",
            labels={"result": "taxa bar plots"},
        ),
    log:
        f"{LOGDIR}/04_taxonomy/{PROJ}_taxa_barplot.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime taxa barplot \
            --i-table {input.table} \
            --i-taxonomy {input.taxonomy} \
            --m-metadata-file {input.metadata} \
            --o-visualization {output} > {log} 2>&1
        """


# The level-N.csv of the bar plot = "Download CSV" in QIIME 2 View
rule barplot_csv:
    input:
        qzv=f"{TAXDIR}/{PROJ}-taxa-bar-plots.qzv",
        script=script("extract_barplot_csv.py"),
    output:
        report(
            f"{TAXDIR}/{PROJ}-level-{{level}}.csv",
            caption=caption("level_csv"),
            category="04 Taxonomy",
            labels={"result": "taxa bar plot table, level {level}"},
        ),
    log:
        f"{LOGDIR}/04_taxonomy/{PROJ}_level-{{level}}_csv.log",
    wildcard_constraints:
        level=r"\d+",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        python {input.script} --qzv {input.qzv} \
            --level {wildcards.level} --output {output} > {log} 2>&1
        """


# Taxa collapsed at the level used for the differential abundance
rule taxa_collapse:
    input:
        table=f"{TAXDIR}/{PROJ}-rarefaction-table-filtered-{GROUP}.qza",
        taxonomy=f"{TAXDIR}/{PROJ}-tax_sklearn.qza",
    output:
        f"{TAXDIR}/{PROJ}-collapse-table-L{{level}}-{GROUP}.qza",
    log:
        f"{LOGDIR}/04_taxonomy/{PROJ}_collapse_L{{level}}.log",
    wildcard_constraints:
        level=r"\d+",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime taxa collapse \
            --i-table {input.table} \
            --i-taxonomy {input.taxonomy} \
            --p-level {wildcards.level} \
            --o-collapsed-table {output} > {log} 2>&1
        """
