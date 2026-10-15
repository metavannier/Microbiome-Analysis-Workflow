# ----------------------------------------------------------------------
# Rooted phylogenetic tree for the phylogenetic diversity metrics
# (Faith PD, UniFrac). Method chosen in config.yaml: sepp or fasttree.
# ----------------------------------------------------------------------

PHYLO = config["phylogeny"]["method"]
# The method is in the file names: changing it builds a new tree and reruns
# the diversity analyses
TREE = f"{PHYLODIR}/{PROJ}-rooted-tree-{PHYLO}.qza"
# Feature table used for ALL the diversity analyses: only the ASVs in the tree
DIVERSITY_TABLE = f"{PHYLODIR}/{PROJ}-taxa-table-filtered-tree-{PHYLO}-{GROUP}.qza"


if PHYLO == "sepp":

    # Insertion of the ASVs into the SEPP reference phylogeny
    rule sepp_phylogeny:
        input:
            rep=f"{TAXDIR}/{PROJ}-rep-filtered-seqs-taxa-{GROUP}.qza",
            seppdb=config["phylogeny"]["seppdb"],
        output:
            tree=TREE,
            placements=f"{PHYLODIR}/{PROJ}-insertion-placements.qza",
        log:
            f"{LOGDIR}/05_phylogeny/{PROJ}_sepp.log",
        threads: config["phylogeny"]["threads"]
        conda:
            QIIME_ENV
        container:
            QIIME_CONTAINER
        shell:
            """
            qiime fragment-insertion sepp \
                --i-representative-sequences {input.rep} \
                --i-reference-database {input.seppdb} \
                --p-threads {threads} \
                --o-tree {output.tree} \
                --o-placements {output.placements} \
                --verbose > {log} 2>&1
            """

    # ASVs too remote to be inserted by SEPP are removed from the table
    # (UniFrac and Faith PD require every ASV of the table in the tree)
    rule sepp_filter_features:
        input:
            table=f"{TAXDIR}/{PROJ}-taxa-table-filtered-{GROUP}.qza",
            tree=TREE,
        output:
            table=DIVERSITY_TABLE,
            removed=f"{PHYLODIR}/{PROJ}-table-not-in-tree-{GROUP}.qza",
        log:
            f"{LOGDIR}/05_phylogeny/{PROJ}_sepp_filter_features.log",
        conda:
            QIIME_ENV
        container:
            QIIME_CONTAINER
        shell:
            """
            qiime fragment-insertion filter-features \
                --i-table {input.table} \
                --i-tree {input.tree} \
                --o-filtered-table {output.table} \
                --o-removed-table {output.removed} > {log} 2>&1
            """

elif PHYLO == "fasttree":

    # MAFFT alignment, masking of the highly variable positions,
    # FastTree and midpoint rooting
    rule fasttree_phylogeny:
        input:
            f"{TAXDIR}/{PROJ}-rep-filtered-seqs-taxa-{GROUP}.qza",
        output:
            alignment=f"{PHYLODIR}/{PROJ}-aligned-rep-filtered-seqs.qza",
            masked=f"{PHYLODIR}/{PROJ}-masked-aligned-rep-filtered-seqs.qza",
            unrooted=f"{PHYLODIR}/{PROJ}-unrooted-tree.qza",
            tree=TREE,
        log:
            f"{LOGDIR}/05_phylogeny/{PROJ}_fasttree.log",
        threads: config["phylogeny"]["threads"]
        conda:
            QIIME_ENV
        container:
            QIIME_CONTAINER
        shell:
            """
            qiime phylogeny align-to-tree-mafft-fasttree \
                --i-sequences {input} \
                --p-n-threads {threads} \
                --o-alignment {output.alignment} \
                --o-masked-alignment {output.masked} \
                --o-tree {output.unrooted} \
                --o-rooted-tree {output.tree} \
                --no-recycle --verbose > {log} 2>&1
            """

    # With a de novo tree every ASV is in the tree: same table
    rule tree_table:
        input:
            f"{TAXDIR}/{PROJ}-taxa-table-filtered-{GROUP}.qza",
        output:
            DIVERSITY_TABLE,
        log:
            f"{LOGDIR}/05_phylogeny/{PROJ}_tree_table.log",
        conda:
            QIIME_ENV
        container:
            QIIME_CONTAINER
        shell:
            """
            cp {input} {output} 2> {log}
            """

else:
    raise WorkflowError(
        f"phylogeny: method must be 'sepp' or 'fasttree', not '{PHYLO}'")


# Reads per sample of the table rarefied by the diversity analyses: choose
# diversity:samplingdepth on this table (with SEPP, ASVs that could not be
# inserted in the tree are not counted)
rule diversity_table_summary:
    input:
        table=DIVERSITY_TABLE,
        metadata=METADATA_FILE,
    output:
        summary=report(
            f"{PHYLODIR}/{PROJ}-taxa-table-filtered-tree-{PHYLO}-{GROUP}.qzv",
            caption=caption("diversitytable"),
            category="04 Taxonomy",
            labels={"result": "feature table summary (ASVs in the tree)"},
        ),
        feature_frequencies=f"{PHYLODIR}/{PROJ}-taxa-table-filtered-tree-{PHYLO}-{GROUP}-feature-frequencies.qza",
        sample_frequencies=f"{PHYLODIR}/{PROJ}-taxa-table-filtered-tree-{PHYLO}-{GROUP}-sample-frequencies.qza",
    log:
        f"{LOGDIR}/05_phylogeny/{PROJ}_diversity_table_summary.log",
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime feature-table summarize \
            --i-table {input.table} \
            --m-metadata-file {input.metadata} \
            --o-feature-frequencies {output.feature_frequencies} \
            --o-sample-frequencies {output.sample_frequencies} \
            --o-summary {output.summary} \
            --no-recycle > {log} 2>&1
        """
