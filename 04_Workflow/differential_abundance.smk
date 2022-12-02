rule taxa_collapse:
    input:
        rarefactionfiltertable = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-rarefaction-table-filtered-" + GROUP + ".qza",
        sklearn = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-tax_sklearn.qza",
    output:
        table_collapse = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-collapse-table-" + GROUP + ".qza",
    log:
        LOGDIR + "/07_differential_abundance/" + PROJ + "_" + GROUP + "-taxa_collapse.log"
    params:
        level = config["ancom"]["taxa_level"]
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime taxa collapse \
        --i-table {input.rarefactionfiltertable} \
        --i-taxonomy {input.sklearn} \
        --p-level {params.level} \
        --o-collapsed-table {output.table_collapse}
        """

rule filter_features:
    input:
        table_collapse = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-collapse-table-" + GROUP + ".qza"
    output:
        table_abond = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-table-abund-" + GROUP + ".qza"
    log:
        LOGDIR + "/07_differential_abundance/" + PROJ + "_" + GROUP + "-filter_features.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime feature-table filter-features \
        --i-table {input.table_collapse} \
        --p-min-frequency {config[ancom][min_frequency]} \
        --p-min-samples {config[ancom][min_sample]} \
        --o-filtered-table {output.table_abond}
        """

rule table_filter_features:
    input:
        table_abond = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-table-abund-" + GROUP + ".qza",
    output:
        table_abond_qzv = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-table-abund-" + GROUP + ".qzv",
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime metadata tabulate --m-input-file {input.table_abond} --o-visualization {output.table_abond_qzv}
        """

rule filter_sample:
    input:
        table_abond = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-table-abund-" + GROUP + ".qza"
    output:
        table_abond_selectedsample = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-table-abund-selectedsample-" + GROUP + ".qza",
        table_abond_qzv = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-table-abund-selectedsample-" + GROUP + ".qzv"
    log:
        LOGDIR + "/07_differential_abundance/" + PROJ + "_" + GROUP + "-filter_features.log"
    params:
        metadata = ROOTDIR + RMSAMPLE
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime feature-table filter-samples \
        --i-table {input.table_abond} \
        --m-metadata-file {params.metadata}\
        --p-exclude-ids True \
        --o-filtered-table {output.table_abond_selectedsample}
        qiime metadata tabulate --m-input-file {output.table_abond_selectedsample} --o-visualization {output.table_abond_qzv}
        """

rule pseudocount:
    input:
        table_abond_selectedsample = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-table-abund-selectedsample-" + GROUP + ".qza"
    output:
        table_abond_comp = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-table-abund-comp-" + GROUP + ".qza"
    log:
        LOGDIR + "/07_differential_abundance/" + PROJ + "_" + GROUP + "-pseudocount.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime composition add-pseudocount \
        --i-table {input.table_abond} \
        --o-composition-table {output.table_abond_comp}
        """

rule ancom:
    input:
        table_abond_comp = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-table-abund-comp-" + GROUP + ".qza"
    output:
        ancom = report(expand(OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-ancom-{column}-" + GROUP + ".qzv", column=COLUMN), caption = ROOTDIR + "/07_Report/ancom.rst", category="07 differential abundance")
    params:
        column = COLUMN,
        metadata = ROOTDIR + "/sample-metadata.tsv"
    log:
        LOGDIR + "/07_differential_abundance/" + PROJ + "_" + GROUP + "-ancom.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        conditions=({output.ancom})
        column=({params.column})
        len=${{#column[@]}}
        for (( i=0; i<$len; i=i+1 ))
        do qiime composition ancom \
        --i-table {input.table_abond_comp} \
        --m-metadata-file {params.metadata} \
        --m-metadata-column ${{column[$i]}} \
        --o-visualization ${{conditions[$i]}}
        done
        """
