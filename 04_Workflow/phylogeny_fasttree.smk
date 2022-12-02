rule alignment:
    input:
        q2_repseq = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-rep-filtered-seqs-taxa-" + GROUP + ".qza",
    output:
        aligned_repseq = OUTPUTDIR + "/05_phylogeny/" + PROJ + "-aligned-rep-filtered-seqs.qza",
        masked_aligned_repseq = OUTPUTDIR + "/05_phylogeny/" + PROJ + "-masked-aligned-rep-filtered-seqs.qza",
    log:
        LOGDIR + "/05_phylogeny/" + PROJ + "_alignment.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime alignment mafft \
        --i-sequences {input.q2_repseq} \
        --o-alignment {output.aligned_repseq}
        qiime alignment mask \
        --i-alignment {output.aligned_repseq} \
        --o-masked-alignment {output.masked_aligned_repseq}
        """

rule phylogeny:
    input:
        masked_aligned_repseq = OUTPUTDIR + "/05_phylogeny/" + PROJ + "-masked-aligned-rep-filtered-seqs.qza",
    output:
        unrooted_tree = OUTPUTDIR + "/05_phylogeny/" + PROJ + "-unrooted-tree.qza",
        rooted_tree = OUTPUTDIR + "/05_phylogeny/" + PROJ + "-rooted-tree.qza",
    log:
        LOGDIR + "/05_phylogeny/" + PROJ + "_phylogeny.log"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime phylogeny fasttree \
        --i-alignment {input.masked_aligned_repseq} \
        --o-tree {output.unrooted_tree}
        qiime phylogeny midpoint-root \
        --i-tree {output.unrooted_tree} \
        --o-rooted-tree {output.rooted_tree}    
        """

#combine the taxa table, rooted phylogeny, taxonomy, and mapping file metadata objects into a single phyloseq object
# rule phyloseq:
#     input:
#         filtertable = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-taxa-table-filtered-" + GROUP + ".qza",
#         rooted_tree = OUTPUTDIR + "/05_phylogeny/" + PROJ + "-rooted-tree.qza",
#         sklearn = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-tax_sklearn.qza",
#         metadata = ROOTDIR + "/sample-metadata.tsv"
#     output:
#         output_phylseq = OUTPUTDIR + "/05_phylogeny/" + PROJ + "-phylseq.tsv",
#     log:
#         LOGDIR + "/05_phylogeny/" + PROJ + "_phylogeny.log"
#     conda:
#         ROOTDIR + "/02_Container/phyloseq.yaml"
#     script:
#         SCRIPTSDIR + "/phyloseq.R"


