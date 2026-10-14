# ----------------------------------------------------------------------
# FastQC to check the raw reads quality
# ----------------------------------------------------------------------

rule fastqc_raw:
    input:
        f"{ROWDATA}/{{stem}}{FASTQ_EXT}",
    output:
        html=f"{QCDIR}/raw/{{stem}}_fastqc.html",
        zip=f"{QCDIR}/raw/{{stem}}_fastqc.zip",
    log:
        f"{LOGDIR}/01_fastqc/raw/{{stem}}.log",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.html),
    threads: 1
    conda:
        conda_env("fastqc.yaml")
    container:
        config["containers"]["fastqc"]
    shell:
        """
        fastqc --quiet --noextract --threads {threads} --outdir {params.outdir} {input} > {log} 2>&1
        """


# ----------------------------------------------------------------------
# Trimmomatic: trimming reads
# ----------------------------------------------------------------------

rule trimmomatic_pe:
    input:
        r1=f"{ROWDATA}/{{pair}}{R1_SUF}",
        r2=f"{ROWDATA}/{{pair}}{R2_SUF}",
    output:
        r1=f"{TRIMDIR}/{{pair}}{R1_SUF}",
        r2=f"{TRIMDIR}/{{pair}}{R2_SUF}",
        # reads where trimming entirely removed the mate
        r1_unpaired=f"{UNPAIREDDIR}/{{pair}}{R1_SUF[:-len(FASTQ_EXT)]}.unpaired{FASTQ_EXT}",
        r2_unpaired=f"{UNPAIREDDIR}/{{pair}}{R2_SUF[:-len(FASTQ_EXT)]}.unpaired{FASTQ_EXT}",
    log:
        f"{LOGDIR}/00_trimmomatic/{{pair}}.log",
    params:
        trimmer=config["trimmomatic"]["trimmer"],
        extra=config["trimmomatic"]["extra"],
    threads: config["trimmomatic"]["threads"]
    conda:
        conda_env("trimmomatic.yaml")
    container:
        config["containers"]["trimmomatic"]
    shell:
        """
        trimmomatic PE -threads {threads} {params.extra} \
            {input.r1} {input.r2} \
            {output.r1} {output.r1_unpaired} \
            {output.r2} {output.r2_unpaired} \
            {params.trimmer} > {log} 2>&1
        """


# ----------------------------------------------------------------------
# FastQC to check the reads trimmed quality
# ----------------------------------------------------------------------

rule fastqc_trimmed:
    input:
        f"{TRIMDIR}/{{stem}}{FASTQ_EXT}",
    output:
        html=f"{QCDIR}/trimmed/{{stem}}_fastqc.html",
        zip=f"{QCDIR}/trimmed/{{stem}}_fastqc.zip",
    log:
        f"{LOGDIR}/01_fastqc/trimmed/{{stem}}.log",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.html),
    threads: 1
    conda:
        conda_env("fastqc.yaml")
    container:
        config["containers"]["fastqc"]
    shell:
        """
        fastqc --quiet --noextract --threads {threads} --outdir {params.outdir} {input} > {log} 2>&1
        """


# ----------------------------------------------------------------------
# MultiQC to summarise the quality before and after trimming
# ----------------------------------------------------------------------

rule multiqc:
    input:
        lambda wildcards: expand(f"{QCDIR}/{wildcards.step}/{{stem}}_fastqc.zip",
                                 stem=FASTQ_STEMS),
    output:
        report(
            f"{QCDIR}/{{step}}_multiqc.html",
            caption=caption("multiqc"),
            category="01 Quality control",
            labels={"reads": "{step}"},
        ),
    log:
        f"{LOGDIR}/01_fastqc/multiqc_{{step}}.log",
    params:
        title=lambda wildcards: f"{PROJ} - {wildcards.step} reads",
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    wildcard_constraints:
        step="raw|trimmed",
    conda:
        conda_env("multiqc.yaml")
    container:
        config["containers"]["multiqc"]
    shell:
        """
        multiqc --force --no-data-dir --title '{params.title}' \
            --outdir {params.outdir} --filename $(basename {output}) \
            {input} > {log} 2>&1
        """
