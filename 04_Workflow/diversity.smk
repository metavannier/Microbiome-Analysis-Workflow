rule filter_taxa_table:
    input:
        filtertable = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-taxa-table-filtered-" + GROUP + ".qza"
    output:
        filtertable_selectedsample = OUTPUTDIR + "/06_diversity/" + PROJ + "-taxa-table-filtered-selectedsample-" + GROUP + ".qza",
        table_count_qzv = OUTPUTDIR + "/06_diversity/" + PROJ + "-taxa-table-filtered-selectedsample-" + GROUP + ".qzv"
    params:
        metadata = ROOTDIR + RMSAMPLE
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime feature-table filter-samples \
        --i-table {input.filtertable} \
        --m-metadata-file {params.metadata}\
        --p-exclude-ids True \
        --o-filtered-table {output.filtertable_selectedsample}        
        qiime metadata tabulate --m-input-file {output.filtertable_selectedsample} --o-visualization {output.table_count_qzv}
        """

rule alpha_rarefaction:
    input:
        filtertable_selectedsample = OUTPUTDIR + "/06_diversity/" + PROJ + "-taxa-table-filtered-selectedsample-" + GROUP + ".qza",
        metadata = ROOTDIR + METADATA
    output:
        rarefaction = report(OUTPUTDIR + "/06_diversity/" + PROJ + "-alpha_rarefaction_curves.qzv", caption = ROOTDIR + "/07_Report/alphararefaction.rst", category="06 diversity/alphadiversity")
    log:
        LOGDIR + "/06_diversity/" + PROJ + "_alpha_rarefaction.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime diversity alpha-rarefaction \
        --i-table {input.filtertable_selectedsample} \
        --m-metadata-file {input.metadata} \
        --o-visualization {output.rarefaction} \
        --p-min-depth {config[diversity][mindepth]} \
        --p-max-depth {config[diversity][maxdepth]}
        """

rule diversity_metrics:
    input:
    # If you want exclude some sample from the analyses :
        filtertable = OUTPUTDIR + "/06_diversity/" + PROJ + "-taxa-table-filtered-selectedsample-" + GROUP + ".qza",
        #filtertable = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-taxa-table-filtered-" + GROUP + ".qza",
        tree = OUTPUTDIR + "/05_phylogeny/" + PROJ + "-rooted-tree.qza",
    output:
        d1 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/rarefied_table.qza",
        d2 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/faith_pd_vector.qza",
        d3 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/observed_features_vector.qza",
        d4 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/shannon_vector.qza",
        d5 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/evenness_vector.qza",
        d6 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/unweighted_unifrac_distance_matrix.qza",
        d7 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/weighted_unifrac_distance_matrix.qza",
        d8 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/jaccard_distance_matrix.qza",
        d9 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/bray_curtis_distance_matrix.qza",
        d11 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/unweighted_unifrac_pcoa_results.qza",
        d12 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/weighted_unifrac_pcoa_results.qza",
        d13 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/jaccard_pcoa_results.qza",
        d14 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/bray_curtis_pcoa_results.qza",
        d15 = report(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/unweighted_unifrac_emperor.qzv", caption = ROOTDIR + "/07_Report/emperor.rst", category="06 diversity/betadiversity"),
        d16 = report(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/weighted_unifrac_emperor.qzv", caption = ROOTDIR + "/07_Report/emperor.rst", category="06 diversity/betadiversity"),
        d17 = report(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/jaccard_emperor.qzv", caption = ROOTDIR + "/07_Report/emperor.rst", category="06 diversity/betadiversity"),
        d18 = report(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/bray_curtis_emperor.qzv", caption = ROOTDIR + "/07_Report/emperor.rst", category="06 diversity/betadiversity"),
    params:
        metadata = ROOTDIR + METADATA
    log:
        LOGDIR + "/06_diversity/" + PROJ + "_alpha_rarefaction.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime diversity core-metrics-phylogenetic \
        --i-table {input.filtertable} \
        --i-phylogeny {input.tree} \
        --m-metadata-file {params.metadata} \
        --p-sampling-depth {config[diversity][samplingdepth]} \
        --o-rarefied-table {output.d1} \
        --o-faith-pd-vector {output.d2} \
        --o-observed-features-vector {output.d3} \
        --o-shannon-vector {output.d4} \
        --o-evenness-vector {output.d5} \
        --o-unweighted-unifrac-distance-matrix {output.d6} \
        --o-weighted-unifrac-distance-matrix {output.d7} \
        --o-jaccard-distance-matrix {output.d8} \
        --o-bray-curtis-distance-matrix {output.d9} \
        --o-unweighted-unifrac-pcoa-results {output.d11} \
        --o-weighted-unifrac-pcoa-results {output.d12} \
        --o-jaccard-pcoa-results {output.d13} \
        --o-bray-curtis-pcoa-results {output.d14} \
        --o-unweighted-unifrac-emperor {output.d15} \
        --o-weighted-unifrac-emperor {output.d16} \
        --o-jaccard-emperor {output.d17} \
        --o-bray-curtis-emperor {output.d18}
        """

rule alpha_significance:
    input:
        alphainput =  expand(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/{alphadiv}_vector.qza", alphadiv=ALPHADIV)
    output:
        alphasigni = report(expand(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/{alphadiv}-group-significance.qzv", alphadiv=ALPHADIV), caption = ROOTDIR + "/07_Report/phylodiv.rst", category="06 diversity/alphadiversity")
    params:
        metadata = ROOTDIR + "/sample-metadata.tsv"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        alphainput=({input.alphainput})
        alphasigni=({output.alphasigni})
        len=${{#alphainput[@]}}
        for (( i=0; i<$len; i=i+1 ))
        do 
        qiime diversity alpha-group-significance \
        --i-alpha-diversity ${{alphainput[$i]}} \
        --m-metadata-file {params.metadata} \
        --o-visualization ${{alphasigni[$i]}}
        done
        """

rule alphacor:
    input:
        alphainput = expand(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/{alphadiv}_vector.qza", alphadiv=ALPHADIV)
    output:
        coroutput = report(expand(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/{alphadiv}-correlation.qzv", alphadiv=ALPHADIV), caption = ROOTDIR + "/07_Report/alphacorrelation.rst", category="06 diversity/alphadiversity")
    params:
        metadata = ROOTDIR + "/sample-metadata.tsv"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        alphainput=({input.alphainput})
        coroutput=({output.coroutput})
        len=${{#alphainput[@]}}
        for (( i=0; i<$len; i=i+1 ))
        do 
        qiime diversity alpha-correlation \
        --i-alpha-diversity ${{alphainput[$i]}} \
        --m-metadata-file {params.metadata} \
        --o-visualization ${{coroutput[$i]}}
        done
        """

CONDITION = expand("{condition.condition}", condition=condition.itertuples())

rule beta_site:
    input:
        betainput =  expand(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/{betadiv}_distance_matrix.qza", betadiv=BETADIV)
    output:
        site = report(expand(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/{condition}-{betadiv}-significance.qzv", condition=CONDITION, betadiv=BETADIV), caption = ROOTDIR + "/07_Report/beta_significance.rst", category="06 diversity/betadiversity")
    params:
        metadata = ROOTDIR + "/sample-metadata.tsv"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        column=({CONDITION})
        betainput=({input.betainput})
        site=({output.site})
        len=${{#column[@]}}
        leninput=${{#betainput[@]}}
        x=0
        for (( j=0; j<$len; j=j+1 ))
        do 
        for (( i=0; i<$leninput; i=i+1 ))
        do 
        qiime diversity beta-group-significance \
        --i-distance-matrix ${{betainput[$i]}} \
        --m-metadata-file {params.metadata} \
        --m-metadata-column ${{column[$j]}} \
        --o-visualization ${{site[$i+$x]}} \
        --p-pairwise
        done
        x=$leninput
        done
        """

rule pcoa:
    input:
        pcoainput = expand(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/{pcoa}_pcoa_results.qza", pcoa=PCOA)
    output:
        pcoaoutput = report(expand(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/{pcoa}-emperor-days.qzv", pcoa=PCOA), caption = ROOTDIR + "/07_Report/emperor.rst", category="06 diversity/betadiversity")
    params:
        metadata = ROOTDIR + "/sample-metadata.tsv"
    log:
        LOGDIR + "/06_diversity/" + PROJ + "_pcoaunweighted.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        pcoainput=({input.pcoainput})
        output=({output.pcoaoutput})
        len=${{#pcoainput[@]}}
        for (( i=0; i<$len; i=i+1 ))
        do 
        qiime emperor plot \
        --i-pcoa ${{pcoainput[$i]}} \
        --m-metadata-file {params.metadata} \
        --p-custom-axes {config[diversity][pcoaaxes]} \
        --o-visualization ${{output[$i]}}
        done
        """
