rule sepp_phylogeny:
    input:
        q2_repseq = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-rep-filtered-seqs-taxa-" + GROUP + ".qza",
        seppdb =  REF + "/" + DB_sepp
    output:
        tree = OUTPUTDIR + "/05_phylogeny/" + PROJ + "-rooted-tree.qza",
        insertion = OUTPUTDIR + "/05_phylogeny/" + PROJ + "-insertion-placements.qza"
    log:
        LOGDIR + "/05_phylogeny/" + PROJ + "_phylogeny.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime fragment-insertion sepp \
        --i-representative-sequences {input.q2_repseq} \
        --i-reference-database {input.seppdb} \
        --o-tree {output.tree} \
        --o-placements {output.insertion} \
        --p-threads {config[phylogeny][sepp_thread]}
        """
