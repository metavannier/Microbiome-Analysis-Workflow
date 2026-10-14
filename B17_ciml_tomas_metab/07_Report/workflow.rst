--------------------------
Microbiome data processing
--------------------------

This workflow performs a 16S rRNA gene (V3-V4) marker-based analysis for the **{{ snakemake.config["project"]["name"] }}** project (B17_ciml_tomas_metab).

{% set comps = snakemake.config["report_comparisons"] %}
{% if comps %}
Comparison(s) presented in this report
======================================

{% for c in comps %}
- **{{ c.id }}** -- {{ c.description }}

  - samples: {{ c.subset if c.subset else "all" }} ({{ c.n_samples }} samples in the metadata/manifest)
  - compared column(s): {{ c.group_columns }} (frequency filtering by group: {{ c.filter_column }})
{% endfor %}

The results of each comparison are in the categories *05 Diversity*, *06 Differential abundance* and *07 MetaboAnalyst* (subcategory = comparison).
{% else %}
This report presents the common part of the analysis (quality, denoising, taxonomy) without comparison.
{% endif %}

---------------------
Materials and methods
---------------------

We follow the `QIIME 2 <https://qiime2.org>`_ (release 2026.7) recommendations for the sequence of workflow for examining amplicon sequence data (`Bolyen et al., 2019 <https://doi.org/10.1038/s41587-019-0209-9>`_).

Some steps produce a summary report in .qzv which can be visualized on `QIIME 2 View <https://view.qiime2.org/>`_ (download, drag and drop the file).

Quality control
===============

- Quality of the raw reads were assessed using `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ (Andrews, 2010) and summarised with `MultiQC <https://doi.org/10.1093/bioinformatics/btw354>`_. Low quality reads were trimmed using `Trimmomatic <https://academic.oup.com/bioinformatics/article/30/15/2114/2390096>`_ (Bolger et al., 2014) with the steps ``{{ snakemake.config["trimmomatic"]["trimmer"] }}``.

- We remove the V3-V4 primers ({{ snakemake.config["primer"]["primerF"] }} / {{ snakemake.config["primer"]["primerR"] }}) using `cutadapt <https://journal.embnet.org/index.php/embnetjournal/article/view/200/479>`_ (error rate {{ snakemake.config["primer"]["primer_err"] }}, minimum overlap {{ snakemake.config["primer"]["primer_overlap"] }}; read pairs without primers are discarded).

Denoising with DADA2
====================

- We use `DADA2 <https://www.nature.com/articles/nmeth.3869>`_ to model and correct Illumina-sequenced amplicon errors and produce the "amplicon sequence variants" (ASVs): reads truncated at {{ snakemake.config["dada2"]["truncation_len-f"] }} (forward) and {{ snakemake.config["dada2"]["truncation_len-r"] }} (reverse) nt, maximum expected errors {{ snakemake.config["dada2"]["quality_err-f"] }}/{{ snakemake.config["dada2"]["quality_err-r"] }}, chimera removal method "{{ snakemake.config["dada2"]["chimera"] }}".

- We keep the merged reads with a length >= {{ snakemake.config["dada2"]["seqlenth"] }} nt to avoid non V3-V4 ribosomal sequences.

Assigning taxonomy
==================

We use a machine learning classifier to assign taxonomies to the ASVs with the `classify-sklearn <https://doi.org/10.1186/s40168-018-0470-z>`_ method (confidence {{ snakemake.config["taxonomy"]["confidence"] }}).
The naive Bayes classifier pre-trained on the `SILVA <https://www.arb-silva.de/>`_ SSU database, V3-V4 region (341f-806r), for QIIME 2 2026.7 is used. ASVs assigned to {{ snakemake.config["taxonomy"]["exclude"] }} are removed.

Samples with fewer than {{ snakemake.config["diversity"]["samplingdepth"] }} reads (sampling depth) are removed from the taxonomic, diversity and abundance analyses.

Phylogeny
=========

We generate a rooted tree for the phylogenetic diversity analyses (method: {{ snakemake.config["phylogeny"]["method"] }}).

- sepp: `SEPP <https://doi.org/10.1128/mSystems.00021-18>`_ inserts the ASVs into the SILVA 128 reference phylogeny. ASVs that cannot be inserted are removed from the diversity analyses.

- fasttree: `MAFFT <https://academic.oup.com/nar/article/30/14/3059/2904316>`_ multiple sequence alignment, masking of the highly variable positions, `FastTree <https://academic.oup.com/mbe/article/26/7/1641/1128976>`_ and midpoint rooting.

Diversity
=========

For each comparison, the feature table is restricted to the samples of the comparison and rarefied at {{ snakemake.config["diversity"]["samplingdepth"] }} reads per sample.

- Alpha diversity (Shannon, observed features, Pielou's evenness, Faith's PD): rarefaction curves and Kruskal-Wallis tests between groups.
- Beta diversity (Jaccard, Bray-Curtis, unweighted and weighted UniFrac): PCoA (Emperor) and PERMANOVA ({{ snakemake.config["diversity"]["permutations"] }} permutations, global and pairwise) for each compared column.

Differential abundance
======================

- Taxa are collapsed at the level {{ snakemake.config["differential_abundance"]["taxa_level"] }} (species).
- We process the frequency filtering by group of samples: one taxon is kept in a group if its frequency is >= {{ snakemake.config["differential_abundance"]["min_frequency"] }} and if it is present in at least {{ snakemake.config["differential_abundance"]["min_samples"] }} samples of the group.
- `Analysis of Composition of Microbiomes (ANCOM) <https://pubmed.ncbi.nlm.nih.gov/26028277/>`_ identifies the taxa that are differentially abundant between the groups of each compared column.
- The same filtered table is formatted for `LEfSe <https://doi.org/10.1186/gb-2011-12-6-r60>`_.

MetaboAnalyst tables
====================

The species-level (level {{ snakemake.config["metaboanalyst"]["taxa_level"] }}) abundance table of the taxa bar plot (identical to the "Download CSV" of QIIME 2 View) is formatted for `MetaboAnalyst <https://www.metaboanalyst.ca/>`_: samples of the comparison in rows, a Group column (compared column of the sample metadata) and the relative abundance of each taxon.
