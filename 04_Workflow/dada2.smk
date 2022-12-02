#########################
# Denoising with dada2  #
#########################

rule dada2:
    input:
        q2_primerRM = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux-noprimer.qza"
    output:
        table = OUTPUTDIR + "/03_dada2/" + PROJ + "-table-" + GROUP + ".qza",
        rep = OUTPUTDIR + "/03_dada2/" + PROJ + "-rep-seqs-dada2.qza",
        stats = OUTPUTDIR + "/03_dada2/" + PROJ + "-dada2-stats.qza"
    log:
        LOGDIR + "/03_dada2/" + PROJ + "_dada2_q2.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime dada2 denoise-paired \
        --i-demultiplexed-seqs {input.q2_primerRM} \
        --p-trunc-q {config[dada2][truncation_err]} \
        --p-trunc-len-f {config[dada2][truncation_len-f]} \
        --p-trunc-len-r {config[dada2][truncation_len-r]} \
        --p-trim-left-f {config[dada2][trim-left-f]} \
        --p-trim-left-r {config[dada2][trim-left-r]} \
        --p-max-ee-f {config[dada2][quality_err-f]} \
        --p-max-ee-r {config[dada2][quality_err-r]} \
        --p-n-reads-learn {config[dada2][training]} \
        --p-n-threads {config[dada2][threads]} \
        --p-chimera-method {config[dada2][chimera]} \
        --o-table {output.table} \
        --o-representative-sequences {output.rep} \
        --o-denoising-stats {output.stats}
        """

rule filterseq:
    input:
        rep = OUTPUTDIR + "/03_dada2/" + PROJ + "-rep-seqs-" + GROUP + ".qza",
    output:
        filterrep = OUTPUTDIR + "/03_dada2/" + PROJ + "-rep-filtered-seqs-" + GROUP + ".qza",
    log:
        LOGDIR + "/03_dada2/" + PROJ + "filterseq.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime feature-table filter-seqs \
        --i-data {input.rep} \
        --m-metadata-file {input.rep} \
        --p-where 'length(sequence) >= {config[dada2][seqlenth]}' \
        --o-filtered-data {output.filterrep} 
        """

rule filterfeature:
    input:
        table = OUTPUTDIR + "/03_dada2/" + PROJ + "-table-" + GROUP + ".qza",
    output:
        filtertable = OUTPUTDIR + "/03_dada2/" + PROJ + "-table-filtered-" + GROUP + ".qza",
    params:
        filterrep = OUTPUTDIR + "/03_dada2/" + PROJ + "-rep-filtered-seqs-" + GROUP + ".qza",
    log:
        LOGDIR + "/03_dada2/" + PROJ + "filterseq.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime feature-table filter-features \
        --i-table {input.table} \
        --m-metadata-file {params.filterrep} \
        --o-filtered-table {output.filtertable}        
        """
