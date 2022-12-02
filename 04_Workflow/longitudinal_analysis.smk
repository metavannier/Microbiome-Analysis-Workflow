rule pairwise_difference:
    input:
        relativefreqtable = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-relative-frequency-" + GROUP + ".qza"
    output:
        output_pairwise_difference = report(expand(OUTPUTDIR + "/08_longitudinal/" + PROJ + "-{feature}-pairwise-differences.qzv", feature=FEATURE), caption = ROOTDIR + "/07_Report/pairwise_difference.rst", category="08 longitudinal/pairwise difference")
    params:
        metadata = ROOTDIR + "/sample-metadata.tsv",
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        output=({output.output_pairwise_difference})
        len=${{#output[@]}}
        for (( i=0; i<$len; i=i+1 ))
        do 
        qiime longitudinal pairwise-differences \
        --m-metadata-file {params.metadata} \
        --i-table {input.relativefreqtable} \
        --p-metric "96f1df5356b17e2d4b6eefc878357fcb" \
        --p-group-column {config[volatility][default_group_column]} \
        --p-state-column {config[volatility][state_column]} \
        --p-state-1 {config[volatility][state1]} \
        --p-state-2 {config[volatility][state2]} \
        --p-individual-id-column {config[volatility][id_column]} \
        --p-replicate-handling random \
        --o-visualization ${{output[$i]}}
        done
        """

rule pairwise_distance:
    input:
        betainput =  expand(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/{betadiv}_distance_matrix.qza", betadiv=BETADIV)
    output:
        output_pairwise_distance = report(expand(OUTPUTDIR + "/08_longitudinal/" + PROJ + "-{betadiv}-pairwise-distances.qzv", betadiv=BETADIV), caption = ROOTDIR + "/07_Report/pairwise_distance.rst", category="08 longitudinal/pairwise distance")
    params:
        metadata = ROOTDIR + "/sample-metadata.tsv",
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        betainput=({input.betainput})
        outpairwisedist=({output.output_pairwise_distance})
        len=${{#betainput[@]}}
        for (( i=0; i<$len; i=i+1 ))
        do 
        qiime longitudinal pairwise-distances \
        --i-distance-matrix ${{betainput[$i]}} \
        --m-metadata-file {params.metadata} \
        --p-group-column {config[volatility][default_group_column]} \
        --p-state-column {config[volatility][state_column]} \
        --p-state-1 {config[volatility][state1]} \
        --p-state-2 {config[volatility][state2]} \
        --p-individual-id-column {config[volatility][id_column]} \
        --p-replicate-handling random \
        --o-visualization ${{outpairwisedist[$i]}}
        done
        """

rule volatility_alpha_diversity:
    input:
        d4 = expand(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/{alphadiv}_vector.qza", alphadiv=ALPHADIV),
    output:
        alphavolatility = report(expand(OUTPUTDIR + "/08_longitudinal/" + PROJ + "-{alphadiv}-volatility.qzv", alphadiv=ALPHADIV), caption = ROOTDIR + "/07_Report/volatility_alpha.rst", category="08 longitudinal/alpha volatility")
    params:
        metadata = ROOTDIR + "/sample-metadata.tsv",
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        alphadiv=({input.d4})
        output=({output.alphavolatility})
        len=${{#alphadiv[@]}}
        for (( i=0; i<$len; i=i+1 ))
        do 
        qiime longitudinal volatility \
        --m-metadata-file {params.metadata} \
        --m-metadata-file ${{alphadiv[$i]}} \
        --p-state-column {config[volatility][state_column]} \
        --p-individual-id-column {config[volatility][id_column]} \
        --p-default-group-column {config[volatility][default_group_column]} \
        --o-visualization ${{output[$i]}}
        done
        """

rule volatility_pcoa:
    input:
        pcoainput = expand(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/{pcoa}_pcoa_results.qza", pcoa=PCOA)
    output:
        pcoaoutput = report(expand(OUTPUTDIR + "/08_longitudinal/" + PROJ + "-{pcoa}_pcoa_results.qzv", pcoa=PCOA), caption = ROOTDIR + "/07_Report/volatility_pcoa.rst", category="08 longitudinal/pcoa volatility")
    params:
        metadata = ROOTDIR + "/sample-metadata.tsv",
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        pcoainput=({input.pcoainput})
        pcoaoutput=({output.pcoaoutput})
        len=${{#pcoainput[@]}}
        for (( i=0; i<$len; i=i+1 ))
        do 
        qiime longitudinal volatility \
        --m-metadata-file {params.metadata} \
        --m-metadata-file ${{pcoainput[$i]}} \
        --p-state-column {config[volatility][state_column]} \
        --p-individual-id-column {config[volatility][id_column]} \
        --p-default-group-column {config[volatility][default_group_column]} \
        --p-default-metric 'Axis 2' \
        --o-visualization ${{pcoaoutput[$i]}}
        done
        """

rule feature_volatility:
    input:
        table_collapse = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-rarefaction-table-filtered-" + GROUP + ".qza",
        sklearn = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-tax_sklearn.qza",
    output:
        featlong1 = OUTPUTDIR + "/08_longitudinal/feat_volatility/filtered_table.qza",
        featlong2 = OUTPUTDIR + "/08_longitudinal/feat_volatility/sample_estimator.qza",
        featlong3 = OUTPUTDIR + "/08_longitudinal/feat_volatility/feature_importance.qza",
        featlong4 = report(OUTPUTDIR + "/08_longitudinal/feat_volatility/accuracy_results.qzv", caption = ROOTDIR + "/07_Report/feature_volatility.rst", category="08 longitudinal/feature volatility"),
        featlong5 = report(OUTPUTDIR + "/08_longitudinal/feat_volatility/volatility_plot.qzv", caption = ROOTDIR + "/07_Report/feature_volatility.rst", category="08 longitudinal/feature volatility"),
        featlong6 = report(OUTPUTDIR + "/08_longitudinal/feat_volatility/feature_importance.qzv", caption = ROOTDIR + "/07_Report/feature_volatility.rst", category="08 longitudinal/feature volatility")
    params:
        metadata = ROOTDIR + "/sample-metadata.tsv",
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:  
        """
        qiime longitudinal feature-volatility \
        --i-table {input.table_collapse} \
        --m-metadata-file {params.metadata} \
        --p-state-column {config[volatility][state_column]} \
        --p-individual-id-column {config[volatility][id_column]} \
        --p-n-estimators {config[volatility][n_estimators]} \
        --p-random-state {config[volatility][random_state]} \
        --p-parameter-tuning \
        --p-estimator {config[volatility][p_estimator]} \
        --o-filtered-table {output.featlong1} \
        --o-feature-importance {output.featlong3} \
        --o-volatility-plot {output.featlong5} \
        --o-accuracy-results {output.featlong4} \
        --o-sample-estimator {output.featlong2} \
        --p-n-jobs {config[volatility][n_jobs]}
        qiime metadata tabulate \
        --m-input-file {output.featlong3} \
        --m-input-file {input.sklearn} \
        --o-visualization {output.featlong6}
        """

rule feature_heatmap:
    input:
        table_collapse = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-rarefaction-table-filtered-" + GROUP + ".qza",
        featlong3 = OUTPUTDIR + "/08_longitudinal/feat_volatility/feature_importance.qza"
    output:
        important_feature_table_top = expand(OUTPUTDIR + "/08_longitudinal/feat_volatility/important-feature-table-top-{heatmap}.qza", heatmap=HEATMAP),
        feature_heatmap = report(expand(OUTPUTDIR + "/08_longitudinal/feat_volatility/important-feature-{heatmap}-heatmap.qzv", heatmap=HEATMAP), caption = ROOTDIR + "/07_Report/feature_heatmap.rst", category="08 longitudinal08 longitudinal/feature volatility")
    params:
        metadata = ROOTDIR + "/sample-metadata.tsv",
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:  
        """
        heatmap=({HEATMAP})
        len=${{#heatmap[@]}}
        feature_heatmap=({output.feature_heatmap})
        important_feature_table_top=({output.important_feature_table_top})
        for (( i=0; i<$len; i=i+1 ))
        do 
        qiime sample-classifier heatmap \
        --i-table {input.table_collapse} \
        --i-importance {input.featlong3} \
        --m-sample-metadata-file {params.metadata} \
        --m-sample-metadata-column ${{heatmap[$i]}} \
        --p-group-samples \
        --p-metric correlation \
        --p-cluster features \
        --p-feature-count {config[volatility][feature_count]} \
        --o-filtered-table ${{important_feature_table_top[$i]}} \
        --o-heatmap ${{feature_heatmap[$i]}}
        done
        """

rule regress_samples:
    input:
        table_collapse = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-rarefaction-table-filtered-" + GROUP + ".qza",
        sklearn = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-tax_sklearn.qza",
    output:
        regressor1 = OUTPUTDIR + "/08_longitudinal/regressor/sample_estimator.qza",
        regressor2 = OUTPUTDIR + "/08_longitudinal/regressor/feature_importance.qza",
        regressor3 = OUTPUTDIR + "/08_longitudinal/regressor/predictions.qza",
        regressor4 = report(OUTPUTDIR + "/08_longitudinal/regressor/accuracy_results.qzv", caption = ROOTDIR + "/07_Report/regress_samples.rst", category="08 longitudinal/regression"),
        regressor5 = OUTPUTDIR + "/08_longitudinal/regressor/model_summary.qzv",
        regressor6 = report(OUTPUTDIR + "/08_longitudinal/regressor/feature_importance.qzv", caption = ROOTDIR + "/07_Report/regress_samples.rst", category="08 longitudinal/regression")
    params:
        metadata = ROOTDIR + "/sample-metadata.tsv",
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:  
        """
        qiime sample-classifier regress-samples \
        --i-table {input.table_collapse} \
        --m-metadata-file {params.metadata} \
        --m-metadata-column {config[volatility][state_column]} \
        --p-estimator {config[volatility][p_estimator]} \
        --p-n-estimators {config[volatility][n_estimators]} \
        --p-random-state {config[volatility][random_state]} \
        --p-parameter-tuning \
        --o-sample-estimator {output.regressor1} \
        --o-feature-importance {output.regressor2} \
        --o-predictions {output.regressor3} \
        --o-model-summary {output.regressor5} \
        --o-accuracy-results {output.regressor4}
        qiime metadata tabulate \
        --m-input-file {output.regressor2} \
        --m-input-file {input.sklearn} \
        --o-visualization {output.regressor6}
        """

rule maturity_index:
    input:
        table_collapse =  OUTPUTDIR + "/04_taxonomy/" + PROJ + "-rarefaction-table-filtered-" + GROUP + ".qza",
        sklearn = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-tax_sklearn.qza",
    output:
        maturity1 = OUTPUTDIR + "/08_longitudinal/maturity/maz_scores.qza",
        maturity2 = OUTPUTDIR + "/08_longitudinal/maturity/sample_estimator.qza",
        maturity3 = OUTPUTDIR + "/08_longitudinal/maturity/feature_importance.qza",
        maturity4 = OUTPUTDIR + "/08_longitudinal/maturity/predictions.qza",
        maturity5 = report(OUTPUTDIR + "/08_longitudinal/maturity/accuracy_results.qzv", caption = ROOTDIR + "/07_Report/maturity_index.rst", category="08 longitudinal/maturity index"),
        maturity6 = report(OUTPUTDIR + "/08_longitudinal/maturity/volatility_plots.qzv", caption = ROOTDIR + "/07_Report/maturity_index.rst", category="08 longitudinal/maturity index"),
        maturity7 = report(OUTPUTDIR + "/08_longitudinal/maturity/clustermap.qzv", caption = ROOTDIR + "/07_Report/maturity_index.rst", category="08 longitudinal/maturity index"),
        maturity8 = OUTPUTDIR + "/08_longitudinal/maturity/model_summary.qzv",
        maturity9 = report(OUTPUTDIR + "/08_longitudinal/maturity/feature_importance.qzv", caption = ROOTDIR + "/07_Report/maturity_index.rst", category="08 longitudinal/maturity index")
    params:
        metadata = ROOTDIR + "/sample-metadata.tsv",
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:  
        """
        qiime longitudinal maturity-index \
        --i-table {input.table_collapse} \
        --m-metadata-file {params.metadata} \
        --p-state-column {config[volatility][state_column]} \
        --p-individual-id-column {config[volatility][id_column]} \
        --p-control {config[volatility][control]} \
        --p-test-size {config[volatility][test_size]} \
        --p-stratify \
        --p-random-state 1010101 \
        --o-sample-estimator {output.maturity2} \
        --o-feature-importance {output.maturity3} \
        --o-predictions {output.maturity4} \
        --o-model-summary {output.maturity8} \
        --o-accuracy-results {output.maturity5} \
        --o-maz-scores {output.maturity1} \
        --o-clustermap {output.maturity7} \
        --o-volatility-plots {output.maturity6} \
        --p-group-by {config[volatility][default_group_column]}
        qiime metadata tabulate \
        --m-input-file {output.maturity3} \
        --m-input-file {input.sklearn} \
        --o-visualization {output.maturity9}
        """