# ----------------------------------------------------------------------
# QIIME 2 import of the trimmed paired-end FASTQ (manifest.tsv)
# ----------------------------------------------------------------------

# manifest.tsv with '$PWD' (and other variables) replaced by absolute
# paths: QIIME 2 then also checks that R1 and R2 have the same number of
# reads for every sample (this check is skipped for '$PWD' paths).
rule expand_manifest:
    input:
        manifest=MANIFEST_FILE,
        script=script("make_manifest.py"),
    output:
        f"{IMPORTDIR}/manifest-absolute-paths.tsv",
    log:
        f"{LOGDIR}/02_qiime_import/expand_manifest.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        python {input.script} --expand {input.manifest} --output {output} > {log} 2>&1
        """


rule qiime_import:
    input:
        manifest=f"{IMPORTDIR}/manifest-absolute-paths.tsv",
        reads=expand(f"{TRIMDIR}/{{name}}", name=FASTQ_NAMES),
    output:
        f"{IMPORTDIR}/{PROJ}-demux-paired-end.qza",
    log:
        f"{LOGDIR}/02_qiime_import/{PROJ}_qiime_import.log",
    params:
        type=config["project"]["type"],
        format=config["project"]["format"],
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        (qiime tools import \
            --type '{params.type}' \
            --input-path {input.manifest} \
            --input-format {params.format} \
            --output-path {output}
        qiime tools validate {output}) > {log} 2>&1
        """


# ----------------------------------------------------------------------
# Removing non-biological sequences (primers) with cutadapt
# ----------------------------------------------------------------------

rule rm_primers:
    input:
        f"{IMPORTDIR}/{PROJ}-demux-paired-end.qza",
    output:
        trimmed=f"{IMPORTDIR}/{PROJ}-PE-demux-noprimer.qza",
        stats=f"{IMPORTDIR}/{PROJ}-cutadapt-stats.qza",
    log:
        f"{LOGDIR}/02_qiime_import/{PROJ}_primer_q2.log",
    params:
        primer_f=config["primer"]["primerF"],
        primer_r=config["primer"]["primerR"],
        error_rate=config["primer"]["primer_err"],
        overlap=config["primer"]["primer_overlap"],
    threads: config["primer"]["threads"]
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences {input} \
            --p-cores {threads} \
            --p-front-f {params.primer_f} \
            --p-front-r {params.primer_r} \
            --p-error-rate {params.error_rate} \
            --p-overlap {params.overlap} \
            --p-match-adapter-wildcards \
            --p-match-read-wildcards \
            --p-discard-untrimmed \
            --o-trimmed-sequences {output.trimmed} \
            --o-stats {output.stats} \
            --verbose > {log} 2>&1
        """


# ----------------------------------------------------------------------
# Reads summaries (counts per sample and interactive quality plots)
# ----------------------------------------------------------------------

rule demux_summarize:
    input:
        f"{IMPORTDIR}/{PROJ}-{{reads}}.qza",
    output:
        report(
            f"{IMPORTDIR}/{PROJ}-{{reads}}.qzv",
            caption=caption("rawsum"),
            category="02 Reads (QIIME 2)",
            labels={"reads": "{reads}"},
        ),
    log:
        f"{LOGDIR}/02_qiime_import/{PROJ}-{{reads}}_summarize.log",
    wildcard_constraints:
        reads="demux-paired-end|PE-demux-noprimer",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime demux summarize --i-data {input} --o-visualization {output} > {log} 2>&1
        """


# Per-sample cutadapt statistics (reads with primers found, written, ...)
rule cutadapt_stats:
    input:
        f"{IMPORTDIR}/{PROJ}-cutadapt-stats.qza",
    output:
        report(
            f"{IMPORTDIR}/{PROJ}-cutadapt-stats.qzv",
            caption=caption("cutadapt"),
            category="02 Reads (QIIME 2)",
            labels={"reads": "cutadapt statistics"},
        ),
    log:
        f"{LOGDIR}/02_qiime_import/{PROJ}_cutadapt_stats.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime metadata tabulate --m-input-file {input} --o-visualization {output} > {log} 2>&1
        """
