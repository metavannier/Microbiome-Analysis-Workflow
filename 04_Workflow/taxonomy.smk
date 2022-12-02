# ----------------------------------------------
# Assigning taxonomy
# ----------------------------------------------

rule assign_tax:
    input:
        q2_repseq = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-rep-filtered-seqs-" + GROUP + ".qza",
        db_classified = REF + "/" + DB_classifier
    output:
        sklearn = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-tax_sklearn.qza"
    log:
        LOGDIR + "/03_taxonomy/" + PROJ + "_sklearn_q2.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"  
    shell:
        """
        qiime feature-classifier classify-sklearn \
        --i-classifier {input.db_classified} \
        --i-reads {input.q2_repseq} \
        --o-classification {output.sklearn} \
        --p-confidence {config[taxonomy][confidence]} \
        --p-n-jobs {config[taxonomy][njobs]}
        """

rule taxa_filter_table:
    input:
        filtertable = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-table-filtered-" + GROUP + ".qza",
        sklearn = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-tax_sklearn.qza"
    output:
        taxafiltertable = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-taxa-table-filtered-" + GROUP + ".qza"
    log:
        LOGDIR + "/03_taxonomy/" + PROJ + "_taxa_filter_table.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell: 
        """
        qiime taxa filter-table \
        --i-table {input.filtertable} \
        --i-taxonomy {input.sklearn} \
        --p-exclude {config[taxonomy][filter]} \
        --o-filtered-table {output.taxafiltertable}
        """

rule filter_seq:
    input:
        q2_repseq = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-rep-filtered-seqs-" + GROUP + ".qza",
        taxafiltertable = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-taxa-table-filtered-" + GROUP + ".qza"
    output:
        q2_repseq_filtered = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-rep-filtered-seqs-taxa-" + GROUP + ".qza"
    log:
        LOGDIR + "/03_taxonomy/" + PROJ + "_filter_seq.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell: 
        """
        qiime feature-table filter-seqs \
        --i-data {input.q2_repseq} \
        --i-table {input.taxafiltertable} \
        --o-filtered-data {output.q2_repseq_filtered}
        """

rule gen_tax:
    input:
        sklearn = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-tax_sklearn.qza"
    output:
        table_tax = OUTPUTDIR + "/04_taxonomy/taxonomy.tsv",
        table_tax_filtered = report(OUTPUTDIR + "/04_taxonomy/taxonomy_filtered.tsv", caption = ROOTDIR + "/07_Report/tax.rst", category="04 taxonomy"),
    log:
        LOGDIR + "/03_taxonomy/" + PROJ + "_exportTAXTSV_q2.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    params:
        directory(OUTPUTDIR + "/04_taxonomy/")
    shell:
        """
        qiime tools export --input-path {input.sklearn} --output-path {params}
        filter=({config[taxonomy][filter]})
        words="${{filter//,/|}}"
        grep -Ev $words {output.table_tax} > {output.table_tax_filtered}
        """

rule stats:
    input:
        filterrep = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-rep-filtered-seqs-taxa-" + GROUP + ".qza",
        stats = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-" + GROUP + "-stats.qza",
        taxafiltertable = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-taxa-table-filtered-" + GROUP + ".qza",
    output:
        rep_viz = report(OUTPUTDIR + "/04_taxonomy/" + PROJ + "-rep-filtered-seqs-taxa-" + GROUP + ".qzv", caption = ROOTDIR + "/07_Report/dada2seq.rst", category="03 dada2"),
        stats_viz = report(OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-" + GROUP + "-stats.qzv", caption = ROOTDIR + "/07_Report/dada2summary.rst", category="03 dada2"),
        featurestat = report(OUTPUTDIR + "/04_taxonomy/" + PROJ + "-taxa-table-filtered-" + GROUP + ".qzv", caption = ROOTDIR + "/07_Report/dada2summary.rst", category="04 taxonomy"),
    params:
        metadata = ROOTDIR + METADATA
    log:
        LOGDIR + "/03_" + GROUP + "/" + PROJ + "_" + GROUP + "-stat_q2.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime metadata tabulate --m-input-file {input.stats} --o-visualization {output.stats_viz}
        qiime feature-table tabulate-seqs --i-data {input.filterrep} --o-visualization {output.rep_viz}
        qiime feature-table summarize --i-table {input.taxafiltertable} --m-sample-metadata-file {params.metadata} --o-visualization {output.featurestat} 
        """

rule filter_rarefaction:
    input:
        taxafiltertable = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-taxa-table-filtered-" + GROUP + ".qza"
    output:
        rarefactionfiltertable = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-rarefaction-table-filtered-" + GROUP + ".qza"
    log:
        LOGDIR + "/04_taxonomy/" + PROJ + "_" + GROUP + "-filter_rarefaction.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime feature-table filter-samples \
        --i-table {input.taxafiltertable} \
        --p-min-frequency {config[diversity][samplingdepth]} \
        --o-filtered-table {output.rarefactionfiltertable}
        """

rule relative_frequency:
    input:
        rarefactionfiltertable = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-rarefaction-table-filtered-" + GROUP + ".qza"
    output:
        relativefreqtable = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-relative-frequency-" + GROUP + ".qza"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime feature-table relative-frequency \
        --i-table {input.rarefactionfiltertable} \
        --o-relative-frequency-table {output.relativefreqtable}
        """

rule gen_table:
    input:
        rarefactionfiltertable = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-rarefaction-table-filtered-" + GROUP + ".qza"
    output:
        table_biom = OUTPUTDIR + "/04_taxonomy/feature-table.biom",
    log:
        LOGDIR + "/03_" + GROUP + "/" + PROJ + "_" + GROUP + "-gen_table.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    params:
        directory(OUTPUTDIR + "/04_taxonomy/")
    shell:
        """
        qiime tools export --input-path {input.rarefactionfiltertable} --output-path {params}
        """

rule convert:
    input:
        table_biom = OUTPUTDIR + "/04_taxonomy/feature-table.biom",
        taxo_table = OUTPUTDIR + "/04_taxonomy/taxonomy.tsv",
    output:
        taxo_table_biom = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-asv-table-with-taxonomy.biom",
        taxo_table_tsv = report(OUTPUTDIR + "/04_taxonomy/" + PROJ + "-asv-table-with-taxonomy.tsv", caption = ROOTDIR + "/07_Report/asv_table.rst", category="03 dada2"),
    params:
        directory(OUTPUTDIR + "/04_taxonomy/")
    log:
        LOGDIR + "/03_" + GROUP + "/" + PROJ + "-asv-table.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        taxo_table={input.taxo_table}
        var="#OTUID\ttaxonomy\tconfidence"
        sed -i "1s/.*/$var/" ${{taxo_table}}
        biom add-metadata -i {input.table_biom} -o {output.taxo_table_biom} --observation-metadata-fp {input.taxo_table} --sc-separated taxonomy
        biom convert -i {output.taxo_table_biom} -o {output.taxo_table_tsv} --to-tsv --header-key taxonomy
        """

rule taxa_barplot:
    input:
        rarefactionfiltertable = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-rarefaction-table-filtered-" + GROUP + ".qza",
        sklearn = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-tax_sklearn.qza",
    output:
        taxabarplots = report(OUTPUTDIR + "/04_taxonomy/" + PROJ + "-taxa-bar-plots.qzv", caption = ROOTDIR + "/07_Report/taxbarplot.rst", category="04 taxonomy")
    log:
        LOGDIR + "/03_taxonomy/" + PROJ + "_taxa_barplot.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    params:
        metadata = ROOTDIR + "/sample-metadata.tsv"
    shell:
        """
        qiime taxa barplot \
        --i-table {input.rarefactionfiltertable} \
        --i-taxonomy {input.sklearn} \
        --m-metadata-file {params.metadata} \
        --o-visualization {output.taxabarplots}
        """
