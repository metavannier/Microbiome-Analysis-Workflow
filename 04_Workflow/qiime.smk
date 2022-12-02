# ----------------------------------------------
# Qiime2 import multiplexed paired-end FASTQ
# with barcodes in sequence
# ----------------------------------------------

rule qiime_import:
  input:
    rowdatapath = ROOTDIR + "/manifest.tsv"

  output:
    outputimport = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-demux-paired-end.qza"

  conda: 
    ROOTDIR + "/02_Container/qiime2.yaml"

  log:
    LOGDIR + "/02_qiime/qiime_import.log"

  shell:
    """
    qiime tools import \
      --type {config[project][type]} \
      --input-path {input.rowdatapath} \
      --input-format {config[project][format]} \
      --output-path {output.outputimport}
    qiime tools validate {output.outputimport}
    """

# Removing non-biological sequences
rule rm_primers:
  input:
    q2_import = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-demux-paired-end.qza"
  output:
    q2_primerRM = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux-noprimer.qza"
  log:
    LOGDIR + "/02_qiime/primer_q2.log"
  conda: 
    ROOTDIR + "/02_Container/qiime2.yaml"
  shell:
    "qiime cutadapt trim-paired \
      --i-demultiplexed-sequences {input.q2_import} \
      --p-front-f {config[primer][primerF]} \
      --p-front-r {config[primer][primerR]} \
      --p-error-rate {config[primer][primer_err]} \
      --p-overlap {config[primer][primer_overlap]} \
      --p-match-adapter-wildcards \
      --p-match-read-wildcards \
      --p-discard-untrimmed \
      --o-trimmed-sequences {output.q2_primerRM}"

rule get_stats:
  input:
    q2_primerRM = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux-noprimer.qza",
  output:
    primer = report(OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux-noprimer.qzv", caption = ROOTDIR + "/07_Report/rawsum.rst", category="02 reads report"),
  log:
    LOGDIR + "/02_qiime/" + PROJ + "_row_stats_q2.log"
  conda: 
    ROOTDIR + "/02_Container/qiime2.yaml"
  shell: 
    """
    qiime demux summarize --i-data {input.q2_primerRM} --o-visualization {output.primer}
    """