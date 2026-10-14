# ----------------------------------------------------------------------
# Denoising with DADA2 (amplicon sequence variants, ASVs)
# ----------------------------------------------------------------------

rule dada2:
    input:
        f"{IMPORTDIR}/{PROJ}-PE-demux-noprimer.qza",
    output:
        table=f"{DADA2DIR}/{PROJ}-table-{GROUP}.qza",
        rep=f"{DADA2DIR}/{PROJ}-rep-seqs-{GROUP}.qza",
        stats=f"{DADA2DIR}/{PROJ}-{GROUP}-stats.qza",
        transitions=f"{DADA2DIR}/{PROJ}-{GROUP}-base-transition-stats.qza",
    log:
        f"{LOGDIR}/03_dada2/{PROJ}_dada2_q2.log",
    params:
        trunc_q=config["dada2"]["truncation_err"],
        trunc_len_f=config["dada2"]["truncation_len-f"],
        trunc_len_r=config["dada2"]["truncation_len-r"],
        trim_left_f=config["dada2"]["trim-left-f"],
        trim_left_r=config["dada2"]["trim-left-r"],
        max_ee_f=config["dada2"]["quality_err-f"],
        max_ee_r=config["dada2"]["quality_err-r"],
        n_reads_learn=config["dada2"]["training"],
        chimera=config["dada2"]["chimera"],
    threads: config["dada2"]["threads"]
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs {input} \
            --p-trunc-q {params.trunc_q} \
            --p-trunc-len-f {params.trunc_len_f} \
            --p-trunc-len-r {params.trunc_len_r} \
            --p-trim-left-f {params.trim_left_f} \
            --p-trim-left-r {params.trim_left_r} \
            --p-max-ee-f {params.max_ee_f} \
            --p-max-ee-r {params.max_ee_r} \
            --p-n-reads-learn {params.n_reads_learn} \
            --p-n-threads {threads} \
            --p-chimera-method {params.chimera} \
            --o-table {output.table} \
            --o-representative-sequences {output.rep} \
            --o-denoising-stats {output.stats} \
            --o-base-transition-stats {output.transitions} \
            --verbose > {log} 2>&1
        """


# Keep the merged sequences long enough to be V3-V4 ribosomal sequences
rule filterseq:
    input:
        f"{DADA2DIR}/{PROJ}-rep-seqs-{GROUP}.qza",
    output:
        f"{DADA2DIR}/{PROJ}-rep-filtered-seqs-{GROUP}.qza",
    log:
        f"{LOGDIR}/03_dada2/{PROJ}_filterseq.log",
    params:
        min_length=config["dada2"]["seqlenth"],
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime feature-table filter-seqs \
            --i-data {input} \
            --m-metadata-file {input} \
            --p-where 'length(sequence) >= {params.min_length}' \
            --o-filtered-data {output} > {log} 2>&1
        """


# Remove from the table the ASVs discarded by the length filter
rule filterfeature:
    input:
        table=f"{DADA2DIR}/{PROJ}-table-{GROUP}.qza",
        filterrep=f"{DADA2DIR}/{PROJ}-rep-filtered-seqs-{GROUP}.qza",
    output:
        f"{DADA2DIR}/{PROJ}-table-filtered-{GROUP}.qza",
    log:
        f"{LOGDIR}/03_dada2/{PROJ}_filterfeature.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime feature-table filter-features \
            --i-table {input.table} \
            --m-metadata-file {input.filterrep} \
            --o-filtered-table {output} > {log} 2>&1
        """


rule dada2_stats:
    input:
        f"{DADA2DIR}/{PROJ}-{GROUP}-stats.qza",
    output:
        report(
            f"{DADA2DIR}/{PROJ}-{GROUP}-stats.qzv",
            caption=caption("dada2summary"),
            category="03 DADA2",
            labels={"result": "denoising statistics"},
        ),
    log:
        f"{LOGDIR}/03_dada2/{PROJ}_dada2_stats.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime metadata tabulate --m-input-file {input} --o-visualization {output} > {log} 2>&1
        """


# Error model learned by DADA2 (observed vs estimated error rates)
rule dada2_error_model:
    input:
        f"{DADA2DIR}/{PROJ}-{GROUP}-base-transition-stats.qza",
    output:
        report(
            f"{DADA2DIR}/{PROJ}-{GROUP}-base-transitions.qzv",
            caption=caption("dada2errors"),
            category="03 DADA2",
            labels={"result": "error model"},
        ),
    log:
        f"{LOGDIR}/03_dada2/{PROJ}_dada2_base_transitions.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime dada2 plot-base-transitions \
            --i-base-transition-stats {input} \
            --o-visualization {output} > {log} 2>&1
        """
