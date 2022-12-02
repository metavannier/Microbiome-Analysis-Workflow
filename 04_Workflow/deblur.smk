#########################
# Denoising with deblur #
#########################

# Merging reads
rule vsearch:
  input:
    q2_primerRM = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux-noprimer.qza"
  output:
    q2_joined = OUTPUTDIR + "/03_deblur/" + PROJ + "-PE-demux-joined.qza"
  log:
    LOGDIR + "/03_deblur/" + PROJ + "_mergePE_q2.log"
  conda: 
    ROOTDIR + "/02_Container/qiime2.yaml"
  shell:
    "qiime vsearch join-pairs \
      --i-demultiplexed-seqs {input.q2_primerRM} \
      --o-joined-sequences {output.q2_joined} \
      --p-minovlen {config[merge][minoverlap]} \
      --p-maxdiffs {config[merge][maxdiff]} \
      --p-minmergelen {config[merge][minlength]}"

rule filter_qscore:
    input:
        q2_joined = OUTPUTDIR + "/03_deblur/" + PROJ + "-PE-demux-joined.qza"
    output:
        q2_filtered = OUTPUTDIR + "/03_deblur/" + PROJ + "-PE-demux-filtered.qza",
        q2_filterstats = OUTPUTDIR + "/03_deblur/" + PROJ + "-PE-demux-filtered-STATS.qza"
    log:
        LOGDIR + "/03_deblur/" + PROJ + "_filterPE_q2.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        "qiime quality-filter q-score \
            --i-demux {input.q2_joined} \
            --o-filtered-sequences {output.q2_filtered} \
            --p-min-quality {config[quality-filter][minphred]} \
            --p-quality-window {config[quality-filter][qualwindow]} \
            --o-filter-stats {output.q2_filterstats}"

rule deblur:
    input:
        q2_filtered = OUTPUTDIR + "/03_deblur/" + PROJ + "-PE-demux-filtered.qza"
    output:
        q2_repseq = OUTPUTDIR + "/03_deblur/" + PROJ + "-rep-seqs-deblur.qza",
        q2_deblurtab = OUTPUTDIR + "/03_deblur/" + PROJ + "-table-deblur.qza",
        q2_deblurstat = OUTPUTDIR + "/03_deblur/" + PROJ + "-deblur-stats.qza"
    log:
        LOGDIR + "/03_deblur/" + PROJ + "_deblur.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"        
    shell: 
        """
        qiime deblur denoise-16S \
        --i-demultiplexed-seqs {input.q2_filtered} \
        --p-trim-length {config[deblur][trim-length]} \
        --o-representative-sequences {output.q2_repseq} \
        --o-table {output.q2_deblurtab} \
        --p-sample-stats \
        --o-stats {output.q2_deblurstat}
        """        

rule get_stats:
    input:
        q2_import = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-demux-paired-end.qza",
        q2_primerRM = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux-noprimer.qza",
        q2_joined = OUTPUTDIR + "/03_deblur/" + PROJ + "-PE-demux-joined.qza",
        q2_filtered = OUTPUTDIR + "/03_deblur/" + PROJ + "-PE-demux-filtered-STATS.qza",
        q2_deblurstat = OUTPUTDIR + "/03_deblur/" + PROJ + "-deblur-stats.qza",
        q2_repseq = OUTPUTDIR + "/03_deblur/" + PROJ + "-rep-seqs-deblur.qza"
    output:
        raw = report(OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux.qzv", caption = ROOTDIR + "/07_Report/rawsum.rst", category="02 reads report"),
        primer = report(OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux-noprimer.qzv", caption = ROOTDIR + "/07_Report/rawsum.rst", category="02 reads report"),
        q2_joined_sum = report(OUTPUTDIR + "/03_deblur/" + PROJ + "-PE-demux-joined.qzv", caption = ROOTDIR + "/07_Report/mergesum.rst", category="03 merge report"),
        q2_filtered_sum = report(OUTPUTDIR + "/03_deblur/" + PROJ + "-PE-demux-filtered-STATS.qzv", caption = ROOTDIR + "/07_Report/filtersum.rst", category="04 qscore filtering report"),
        q2_deblurstat_sum = report(OUTPUTDIR + "/03_deblur/" + PROJ + "-deblur-stats.qzv", caption = ROOTDIR + "/07_Report/deblursummary.rst", category="05 grouping similar sequences report"),
        q2_repseq_sum = report(OUTPUTDIR + "/03_deblur/" + PROJ + "-rep-seqs-deblur.qzv", caption = ROOTDIR + "/07_Report/deblursummary.rst", category="05 grouping similar sequences report")
    log:
        LOGDIR + "/03_otu/" + PROJ + "_get_stats_q2.log"
    conda: 
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell: 
        """
        qiime demux summarize --i-data {input.q2_import} --o-visualization {output.raw}
        qiime demux summarize --i-data {input.q2_primerRM} --o-visualization {output.primer}
        qiime demux summarize --i-data {input.q2_joined} --o-visualization {output.q2_joined_sum}
        qiime metadata tabulate --m-input-file {input.q2_filtered} --o-visualization {output.q2_filtered_sum}
        qiime deblur visualize-stats --i-deblur-stats {input.q2_deblurstat} --o-visualization {output.q2_deblurstat_sum}
        qiime feature-table tabulate-seqs --i-data {input.q2_repseq} --o-visualization {output.q2_repseq_sum}
        """

rule gen_table:
    input:
        q2_deblurtab = OUTPUTDIR + "/03_deblur/" + PROJ + "-table-" GROUP + ".qza",
    output:
        table_otu = OUTPUTDIR + "/03_deblur/feature-table.biom"
    log:
        LOGDIR + "/03_deblur/" + PROJ + "_deblur-gen_table.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    params:
        directory(OUTPUTDIR + "/03_deblur/")
    shell:
        "qiime tools export --input-path {input.q2_deblurtab} --output-path {params}"

rule convert:
    input:
        table_otu = OUTPUTDIR + "/03_deblur/feature-table.biom"
    output:
        otu_table = OUTPUTDIR + "/03_deblur/" + PROJ + "-otu-table.tsv"
    log:
        LOGDIR + "/03_deblur/" + PROJ + "-otu-table.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        "biom convert -i {input} -o {output} --to-tsv"