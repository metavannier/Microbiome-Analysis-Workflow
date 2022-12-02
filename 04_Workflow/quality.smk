# ----------------------------------------------
# FastQC to check the row reads quality
# ----------------------------------------------

rule fastqc:
  input:
    ROWDATA + "/{sample}_L001_{num}_001.fastq.gz"
  output:
    html = OUTPUTDIR + "/01_fastqc/{sample}_L001_{num}_fastqc.html",
    zip = OUTPUTDIR + "/01_fastqc/{sample}_L001_{num}_fastqc.zip"
  params: ""
  log:
    LOGDIR + "/01_fastqc/{sample}_L001_{num}.log"
  wrapper:
    "0.35.2/bio/fastqc"

# ----------------------------------------------
# Trimmomatic: trimming reads
# ----------------------------------------------

rule trimmomatic_pe:
  input:
    r1 = ROWDATA + "/{sample}_L001_" + R1_SUF + SUF,
    r2 = ROWDATA + "/{sample}_L001_" + R2_SUF + SUF,
  output:
    r1 = OUTPUTDIR + "/00_trimmed/{sample}_L001_" + R1_SUF + SUF,
    r2 = OUTPUTDIR + "/00_trimmed/{sample}_L001_" + R2_SUF + SUF,
    # reads where trimming entirely removed the mate
    r1_unpaired = OUTPUTDIR + "/00_untrimmed/{sample}_L001_" + R1_SUF + ".unpaired.fastq.gz",
    r2_unpaired = OUTPUTDIR + "/00_untrimmed/{sample}_L001_" + R2_SUF + ".unpaired.fastq.gz",
  log:
    log = LOGDIR + "/00_trimmomatic/{sample}.log"  
  params:
    trimmer = ["LEADING:2", "TRAILING:2", "SLIDINGWINDOW:4:2", "MINLEN:25"],
    extra = ""
  wrapper:
    "0.35.2/bio/trimmomatic/pe"

# ----------------------------------------------
# FastQC to check the reads trimmed quality
# ----------------------------------------------

rule fastqc_trimmed:
  input:
    OUTPUTDIR + "/00_trimmed/{sample}_L001_{num}_001.fastq.gz"
  output:
    html = OUTPUTDIR + "/01_fastqc/{sample}_L001_{num}_trimmed_fastqc.html",
    zip = OUTPUTDIR + "/01_fastqc/{sample}_L001_{num}_trimmed_fastqc.zip"
  params: ""
  log:
    LOGDIR + "/01_fastqc/{sample}_L001_{num}_trimmed.log"
  wrapper:
    "0.35.2/bio/fastqc"

# ----------------------------------------------
# MultiQC to check the reads trimmed quality
# ----------------------------------------------

rule multiqc_untrimmed:
  input:
    raw_qc = expand("{outputdir}/01_fastqc/{sample}_L001_{num}_fastqc.zip", outputdir = OUTPUTDIR, sample=SAMPLE_SET, num=SET_NUMS),
  output:
    raw_multi_html = report(OUTPUTDIR + "/01_fastqc/raw_multiqc.html", caption = ROOTDIR + "/07_Report/multiqc.rst", category="01 quality report"), 
  params:
    multiqc_output_raw = OUTPUTDIR + "/01_fastqc/raw_multiqc_data"
  conda:
    ROOTDIR + "/02_Container/multiqc.yaml"
  shell: 
    """
    multiqc -n {output.raw_multi_html} {input.raw_qc} #run multiqc
    rm -rf {params.multiqc_output_raw} #clean-up
    """ 

rule multiqc_trimmed:
  input:
    trim_qc = expand("{outputdir}/01_fastqc/{sample}_L001_{num}_trimmed_fastqc.zip", outputdir = OUTPUTDIR, sample=SAMPLE_SET, num=SET_NUMS)
  output:
    trim_multi_html = report(OUTPUTDIR + "/01_fastqc/trimmed_multiqc.html", caption = ROOTDIR + "/07_Report/multiqc.rst", category="01 quality report"), 
  params:
    multiqc_output_trim = OUTPUTDIR + "/01_fastqc/trimmed_multiqc_data"
  conda:
    ROOTDIR + "/02_Container/multiqc.yaml"
  shell: 
    """
    multiqc -n {output.trim_multi_html} {input.trim_qc} --force #run multiqc
    rm -rf {params.multiqc_output_trim} #clean-up
    """
