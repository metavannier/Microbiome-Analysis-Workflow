#######
# OTU #
#######

# Merging reads
rule vsearch:
  input:
    q2_primerRM = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux-noprimer.qza"
  output:
    q2_joined = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined.qza"
  log:
    LOGDIR + "/03_otu/" + PROJ + "_mergePE_q2.log"
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
    q2_joined = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined.qza"
  output:
    q2_filtered = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined-filtered.qza",
    q2_filterstats = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined-filtered-STATS.qza"
  log:
    LOGDIR + "/03_otu/" + PROJ + "_filterPE_q2.log"
  conda:
    ROOTDIR + "/02_Container/qiime2.yaml"
  shell:
    "qiime quality-filter q-score \
      --i-demux {input.q2_joined} \
      --o-filtered-sequences {output.q2_filtered} \
      --p-min-quality {config[quality-filter][minphred]} \
      --p-quality-window {config[quality-filter][qualwindow]} \
      --o-filter-stats {output.q2_filterstats}"

rule get_stats:
  input:
    q2_import = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-demux-paired-end.qza",
    q2_primerRM = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux-noprimer.qza",
    q2_joined = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined.qza",
    q2_filtered = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined-filtered.qza",
  output:
    raw = report(OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux.qzv", caption = ROOTDIR + "/07_Report/rawsum.rst", category="02 reads report"),
    primer = report(OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux-noprimer.qzv", caption = ROOTDIR + "/07_Report/rawsum.rst", category="02 reads report"),
    q2_joined_sum = report(OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined.qzv", caption = ROOTDIR + "/07_Report/mergesum.rst", category="03 merging report"),
    q2_filtered_sum = report(OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined-filtered.qzv", caption = ROOTDIR + "/07_Report/filtersum.rst", category="04 qscore filtering report"),
  log:
    LOGDIR + "/03_otu/" + PROJ + "_get_stats_q2.log"
  conda: 
    ROOTDIR + "/02_Container/qiime2.yaml"
  shell: 
    """
    qiime demux summarize --i-data {input.q2_import} --o-visualization {output.raw}
    qiime demux summarize --i-data {input.q2_primerRM} --o-visualization {output.primer}
    qiime demux summarize --i-data {input.q2_joined} --o-visualization {output.q2_joined_sum}
    qiime demux summarize --i-data {input.q2_filtered} --o-visualization {output.q2_filtered_sum}
    """

rule dedup:
  input:
    q2_filtered = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined-filtered.qza",
  output:
    q2_dedup = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined-filtered-dedup_table.qza",
    q2_dedup_seq = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined-filtered-dedup_seqs.qza"
  log:
    LOGDIR + "/03_otu/" + PROJ + "_dedupPE_q2.log"
  conda:
    ROOTDIR + "/02_Container/qiime2.yaml"
  shell:
    "qiime vsearch dereplicate-sequences \
      --i-sequences {input.q2_filtered} \
      --o-dereplicated-table {output.q2_dedup} \
      --o-dereplicated-sequences {output.q2_dedup_seq}"

rule denovo:
  input:
    q2_dedup = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined-filtered-dedup_table.qza",
    q2_dedup_seq = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined-filtered-dedup_seqs.qza"
  output:
    q2_cluster_table = OUTPUTDIR + "/03_otu/" + PROJ + "-cluster-table.qza",
    q2_cluster_seqs = OUTPUTDIR + "/03_otu/" + PROJ + "-cluster-seqs.qza",
  log:
    LOGDIR + "/03_otu/" + PROJ + "_dedupPE_q2.log"
  conda:
    ROOTDIR + "/02_Container/qiime2.yaml"
  shell:
    "qiime vsearch cluster-features-de-novo \
      --i-table {input.q2_dedup} \
      --i-sequences {input.q2_dedup_seq} \
      --p-perc-identity {config[cluster-denovo][denovo_perc_id]} \
      --o-clustered-table {output.q2_cluster_table} \
      --o-clustered-sequences {output.q2_cluster_seqs}  \
      --p-threads {config[cluster-denovo][denovo_otu-thread]}"

rule chimera_find:
  input:
    q2_cluster_table = OUTPUTDIR + "/03_otu/" + PROJ + "-cluster-table.qza",
    q2_cluster_seqs = OUTPUTDIR + "/03_otu/" + PROJ + "-cluster-seqs.qza",
  output:
    q2_chimeras = OUTPUTDIR + "/03_otu/" + PROJ + "-chimeras.qza",
    q2_nonchimeras = OUTPUTDIR + "/03_otu/" + PROJ + "-nonchimeras.qza",
    q2_chimeras_stats = OUTPUTDIR + "/03_otu/" + PROJ + "-chimeras-STATS.qza"
  log:
    LOGDIR + "/03_otu/" + PROJ + "_dedupPE_q2.log"
  conda:
    ROOTDIR + "/02_Container/qiime2.yaml"
  shell:
    "qiime vsearch uchime-denovo \
	    --i-table {input.q2_cluster_table} \
	    --i-sequences  {input.q2_cluster_seqs} \
	    --o-chimeras {output.q2_chimeras} \
      --o-nonchimeras {output.q2_nonchimeras} \
      --o-stats {output.q2_chimeras_stats}"
  
rule chimera_rm_table:
  input:
    q2_cluster_table = OUTPUTDIR + "/03_otu/" + PROJ + "-cluster-table.qza",
    q2_nonchimeras = OUTPUTDIR + "/03_otu/" + PROJ + "-nonchimeras.qza"
  output:
    q2_table_nc = OUTPUTDIR + "/03_otu/" + PROJ + "-table-nc.qza"
  log:
    LOGDIR + "/03_otu/" + PROJ + "_chimerarm_table_q2.log"
  conda:
    ROOTDIR + "/02_Container/qiime2.yaml"
  shell:
    "qiime feature-table filter-features \
	    --i-table {input.q2_cluster_table} \
      --m-metadata-file {input.q2_nonchimeras} \
      --o-filtered-table {output.q2_table_nc}"

rule chimera_rm_seq:
  input:
    q2_cluster_seqs = OUTPUTDIR + "/03_otu/" + PROJ + "-cluster-seqs.qza",
    q2_nonchimeras = OUTPUTDIR + "/03_otu/" + PROJ + "-nonchimeras.qza",
    q2_table_nc = OUTPUTDIR + "/03_otu/" + PROJ + "-table-nc.qza"
  output:
    q2_seqs_nc = OUTPUTDIR + "/03_otu/" + PROJ + "-rep-seqs-otu.qza
  log:
    LOGDIR + "/03_otu/" + PROJ + "_chimerarm_seqs_q2.log"
  conda:
    ROOTDIR + "/02_Container/qiime2.yaml"
  shell:
    "qiime feature-table filter-seqs \
        --i-data {input.q2_cluster_seqs} \
        --m-metadata-file {input.q2_nonchimeras} \
        --o-filtered-data {output.q2_seqs_nc}"

rule summarize_chimera:
  input:
    q2_table_nc = OUTPUTDIR + "/03_otu/" + PROJ + "-table-nc.qza"
  output:
    q2_table_nc_qzv = OUTPUTDIR + "/03_otu/" + PROJ + "-table-nc.qzv"
  log:
    LOGDIR + "/03_otu/" + PROJ + "_chimerarm_seqs_sum_q2.log"
  conda:
    ROOTDIR + "/02_Container/qiime2.yaml"
  shell:
    "qiime feature-table summarize \
        --i-table {input.q2_table_nc} \
        --o-visualization {output.q2_table_nc_qzv}"

rule gen_table_otu:
  input:
    q2_table_nc = OUTPUTDIR + "/03_otu/" + PROJ + "-table-nc.qza"
  output:
    table_biom_otu = OUTPUTDIR + "/03_otu/feature-table.biom"    
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_otutable_BIOM.log"
  conda:
    ROOTDIR + "/02_Container/qiime2.yaml"
  params:
    directory(OUTPUTDIR + "/03_otu/")
  shell:
    "qiime tools export --input-path {input.q2_table_nc} --output-path {params}"

rule convert_otu:
  input:
    table_biom_otu = OUTPUTDIR + "/03_otu/" + PROJ + "-feature-table.biom"    
  output:
    table_tsv_otu = OUTPUTDIR + "/03_otu/" + PROJ + "-otu-table.tsv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_otuexportTSV_q2.log"
  conda:
    ROOTDIR + "/02_Container/qiime2.yaml"
  shell:
    "biom convert -i {input} -o {output} --to-tsv"

# Export fastq for fastqc
# rule qiime_export:
#   input:
#     inputqza = OUTPUTDIR + "/03_otu/" + PROJ + "-seqs-nc.qza"
#   output:
#     outputfastqpath = OUTPUTDIR + "/04_fastqfiltered/"
#   conda: 
#     ROOTDIR + "/02_Container/qiime2.yaml"
#   log:
#     LOGDIR + "/02_qiime/qiime_fastqfiltered.log"
#   shell:
#     """
#     qiime tools export \
#       --input-path {input.inputqza} \
#       --output-path {output.outputfastqpath}
#     """

# rule fastqc_filtered:
#   input:
#     OUTPUTDIR + "/00_trimmed/{sample}_L001_{num}_001.fastq.gz"
#   output:
#     html = OUTPUTDIR + "/01_fastqc/{sample}_L001_{num}_trimmed_fastqc.html",
#     zip = OUTPUTDIR + "/01_fastqc/{sample}_L001_{num}_trimmed_fastqc.zip"
#   params: ""
#   log:
#     LOGDIR + "/01_fastqc/{sample}_L001_{num}_trimmed.log"
#   wrapper:
#     "0.35.2/bio/fastqc"