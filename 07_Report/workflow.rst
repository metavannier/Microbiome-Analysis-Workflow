--------------------------
Microbiome data processing
--------------------------

This workflow performs a gene marker-based analysis for the weaning project.

---------------------
Materials and methods
---------------------

We follow the `qiime2docs recommandations <https://docs.qiime2.org/2021.2/tutorials/qiime2-for-experienced-microbiome-researchers/#id8>`_ for the sequence of workflow for examining amplicon sequence data.

Some steps produce a summary report in .qzv which can be visualized on `quiime2view <https://view.qiime2.org/>`_

- Quality of the raw reads were assessed using `FastQC v0.11.9 toolkit <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ (Andrews, 2010). Low quality reads were trimmed using `Trimmomatic v0.39 <https://academic.oup.com/bioinformatics/article/30/15/2114/2390096>`_ (Bolger et al., 2014).

- We remove the V3_V4 primer with overhang adapters (FwOvAd_341F and RevOvAd_785R) using `cutadapt <https://journal.embnet.org/index.php/embnetjournal/article/view/200/479>`_ following the parameters from Parada et al. (Parada et al. 2016 MEP primers for 16S).

Grouping similar sequences
==========================

There are two main approaches for grouping similar sequences together: denoising and clustering.

OTU Clustering
--------------

This method perform de novo and closed reference clustering.

- The reads are merged and dereplicating using `VSEARCH <https://peerj.com/articles/2584/>`_ (Rognes et al. 2016) with the join-pairs method.

- Low-quality reads are discarded by filtering using the per-nucleotide phred quality scores using the q-score method with recommanded parameter from `Bokulich and al. <https://www.nature.com/articles/nmeth.2276>`_

- Sequences are clustered de novo based on their genetic similarity (threshold 0.97) with VSEARCH.

- Chimeric feature sequences are discarted using vsearch uchime_denovo method. 

Denoising with deblur
---------------------

- The reads are merged using `VSEARCH <https://peerj.com/articles/2584/>`_ (Rognes et al. 2016) with the join-pairs method.

- Low-quality reads are discarded by filtering using the per-nucleotide phred quality scores using the q-score method with recommanded parameter from `Bokulich and al. <https://www.nature.com/articles/nmeth.2276>`_

- We use `deblur <https://msystems.asm.org/content/2/2/e00191-16>`_ denoise-16S plugin to denoise sequences (it discards any reads which do not have a minimum 60% identity similarity to sequences from the 85% OTU GreenGenes database) and to extract the sub-OTUs.

Used in this workflow: Denoising with dada2
-------------------------------------------

- We use `DADA2 <https://www.nature.com/articles/nmeth.3869>`_ to modeling and correcting Illumina-sequenced amplicon errors and produce the "amplicon sequence variants” (ASVs).

- We keep the merged read with a length >= 400 nt to avoid non v3-v4 ribosomal sequence.

Assigning taxonomy
==================

We use a machine learning classifier to assign taxonomies to reads with the `classify-sklearn <https://www.jmlr.org/papers/volume12/pedregosa11a/pedregosa11a.pdf>`_ method.

A Naive Bayes classifiers trained on the Silva (release 138) 99% OTUs full-length sequences are applied to obtain the pre-trained taxonomy classifiers.

Phylogeny
=========

We generate a rooted tree for phylogenetic diversity analyses. Two method can be use :

- align-to-tree-mafft-fasttree : A pipeline using the `MAFFT <https://academic.oup.com/nar/article/30/14/3059/2904316>`_ program to perform a multiple sequence alignment, filters the alignment to remove positions that are highly variable. Then `FastTree <https://academic.oup.com/mbe/article/26/7/1641/1128976>`_ to generate a phylogenetic tree from the masked alignment. A midpoint rooting is applied to have the rooted tree.   

- USED IN THIS WORKFLOW: `SEPP <https://www.worldscientific.com/doi/abs/10.1142/9789814366496_0024>`_ whose address the problem of Phylogenetic Placement. The objective is to insert query sequences into an existing phylogenetic tree and alignment on full-length sequences for the same gene. 

Analysis
========

To visualize the figure you have to download, drag and drop the .qzv file on `quiime2view <https://view.qiime2.org/>`_ or open it in a browser when it is a html file.

Taxonomic and abundance analyses are performed with the feature table [frequency] obtained after quality filtering steps of the OTUs/ASV (sequence quality, sequence lenght, taxonomic validity, sampling depth threshold, minimal frequency, minimal sample, taxa collapsing). 

Diversity metrics are performed with collapsing at the species taxa. If you decide to not collapse you will stay at the "nucleic" resolution (due to the ASVs method). 

Quality reports
---------------

- The quality reports of the samples are in `01 quality report/raw_multiqc.html </media/thomas/data/ciml_tomas_metab_weaning/05_Output/01_fastqc/raw_multiqc.html>`_ for the fastq file before treaming and `01 quality report/trimmed_multiqc.html </media/thomas/data/ciml_tomas_metab_weaning/05_Output/01_fastqc/raw_multiqc.html>`_ for the fastq file after treaming.

- The demultiplexed sequences counts summary and the quality plot after removing the primers are in `/02 reads report/weaning-PE-demux-noprimer.qzv  </media/thomas/data/ciml_tomas_metab_weaning/05_Output/02_qiime_import/weaning-PE-demux-noprimer.qzv>`_ .

Grouping similar sequences (OTU, ASVs or sub-OTUs)
--------------------------------------------------

- For ASVs method, the ASV sequences are in `/03 dada2/weaning-rep-filtered-seqs-taxa-dada2.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/04_taxonomy/weaning-rep-filtered-seqs-taxa-dada2.qzv>`_ .

- Statistics of the dada2 process are in `/03 dada2/weaning-dada2-stats.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/03_dada2/weaning-dada2-stats.qzv>`_ .

Taxonomic table and visualization
---------------------------------

=> 3 samples are removed from the taxonomy and diversity analysis due to theyr low depth of sequencing (<7500) : M6946, M6938 and M6935

- Frequency per sample and sampling depth information are in `/04 taxonomy/weaning-taxa-table-filtered-dada2.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/04_taxonomy/weaning-taxa-table-filtered-dada2.qzv>`_ .

- The table of the taxonomic assignation for each feature (OTU or ASV) are in `/04 taxonomy/taxonomy_filtered.tsv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/04_taxonomy/taxonomy_filtered.tsv>`_ .

- You can find the interactive barplot visualization of taxonomies in `/04 taxonomy/weaning-taxa-bar-plots.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/04_taxonomy/weaning-taxa-bar-plots.qzv>`_ .

- The biom table with the feature count for each samples and the taxonomic affiliation are in `/03 dada2/table-with-taxonomy.tsv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/04_taxonomy/weaning-asv-table-with-taxonomy.tsv>`_ .

Alpha Diversity Analysis
------------------------

- Alpha Rarefaction for the selection of the rarefaction depth

`/06 diversity/alphadiversity/weaning-alpha_rarefaction_curves.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/06_diversity/weaning-alpha_rarefaction_curves.qzv>`_

- Files for visually and statistically compare groups of alpha diversity values calculated with different metrics.

`/06 diversity/alphadiversity/evenness-group-significance.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/06_diversity/weaning-core-metrics-results/evenness-group-significance.qzv>`_

`/06 diversity/alphadiversity/faith_pd-group-significance.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/06_diversity/weaning-core-metrics-results/faith_pd-group-significance.qzv>`_

`/06 diversity/alphadiversity/observed_features-group-significance.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/06_diversity/weaning-core-metrics-results/observed_features-group-significance.qzv>`_

`/06 diversity/alphadiversity/shannon-group-significance.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/06_diversity/weaning-core-metrics-results/shannon-group-significance.qzv>`_

- File representing the alpha-correlation to determine whether numeric sample metadata columns are correlated with alpha diversity.

`/06 diversity/alphadiversity/evenness-correlation.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/06_diversity/weaning-core-metrics-results/evenness-correlation.qzv>`_

`/06 diversity/alphadiversity/faith_pd-correlation.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/06_diversity/weaning-core-metrics-results/faith_pd-correlation.qzv>`_

`/06 diversity/alphadiversity/observed_features-correlation.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/06_diversity/weaning-core-metrics-results/observed_features-correlation.qzv>`_

`/06 diversity/alphadiversity/shannon-correlation.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/06_diversity/weaning-core-metrics-results/shannon-correlation.qzv>`_

Beta Diversity Analysis
-----------------------

- File representing the sample composition in the context of categorical metadata using PERMANOVA.

/06 diversity/betadiversity/"condition"-"betadiversity_distance"-significance.qzv 

Principal component analysis
----------------------------

- Files for the ACP visualisation of the beta diversity values calculated with different metrics.

`06 diversity/betadiversity/unweighted_unifrac_emperor.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/06_diversity/weaning-core-metrics-results/unweighted_unifrac_emperor.qzv>`_

`06 diversity/betadiversity/weighted_unifrac_emperor.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/06_diversity/weaning-core-metrics-results/weighted_unifrac_emperor.qzv>`_

`06 diversity/betadiversity/bray_curtis_emperor.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/06_diversity/weaning-core-metrics-results/bray_curtis_emperor.qzv>`_

`06 diversity/betadiversity/jaccard_emperor.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/06_diversity/weaning-core-metrics-results/jaccard_emperor.qzv>`_

- Principal component distribution of the different kinetic points with the betadiversity distance.

`06 diversity/betadiversity/unweighted_unifrac-emperor-days.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/06_diversity/weaning-core-metrics-results/unweighted_unifrac-emperor-days.qzv>`_

`06 diversity/betadiversity/weighted_unifrac-emperor-days.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/06_diversity/weaning-core-metrics-results/weighted_unifrac-emperor-days.qzv>`_

`06 diversity/betadiversity/bray_curtis-emperor-days.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/06_diversity/weaning-core-metrics-results/bray_curtis-emperor-days.qzv>`_

`06 diversity/betadiversity/jaccard-emperor-days.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/06_diversity/weaning-core-metrics-results/jaccard-emperor-days.qzv>`_

Differential abundance with ANCOM
---------------------------------

- We process the frequency filtering by group of samples : One taxon at the species level is keep if the frequency is > 50 and if it is present in at minimum 2 samples.

- `Analysis of Composition of Microbiomes (ANCOM) <https://pubmed.ncbi.nlm.nih.gov/26028277/>`_ to identify features that are differentially abundant across groups.

Differential abondance between Peyer's patch and Villus against other compartiment

`07 differential_abundance/weaning-ancom-PP_villi_vs_other-dada2.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/07_differential_abundance/weaning-ancom-PP_villi_vs_other-dada2.qzv>`_

Differential abondance between Peyer's patch against other compartiment

`07 differential_abundance/weaning-ancom-PP_vs_other-dada2.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/07_differential_abundance/weaning-ancom-PP_vs_other-dada2.qzv>`_

Relative abundance table for LEfSe tools
----------------------------------------

`07 differential_abundance/weaning-table_abond_lefse.tsv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/07_differential_abundance/weaning-table_abond_lefse.tsv>`_

Longitudinal analysis
---------------------

- Pairwise difference tests determine whether the value of an ASVs changed significantly between pairs of paired samples

`08 longitudinal/pairwise difference/weaning-96f1df5356b17e2d4b6eefc878357fcb-pairwise-differences.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/08_longitudinal/weaning-96f1df5356b17e2d4b6eefc878357fcb-pairwise-differences.qzv>`_

- The pairwise-distances visualizer

`08 longitudinal/pairwise distance/unweighted_unifrac-pairwise-distances.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/08_longitudinal/weaning-unweighted_unifrac-pairwise-distances.qzv>`_

`08 longitudinal/pairwise distance/weighted_unifrac-pairwise-distances.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/08_longitudinal/weaning-weighted_unifrac-pairwise-distances.qzv>`_

`08 longitudinal/pairwise distance/bray_curtis-pairwise-distances.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/08_longitudinal/weaning-bray_curtis-pairwise-distances.qzv>`_

`08 longitudinal/pairwise distance/jaccard-pairwise-distances.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/08_longitudinal/weaning-jaccard-pairwise-distances.qzv>`_

- Volatility visualizer

08_longitudinal/alpha volatility/weaning-{alphadiv}-volatility.qzv : Examine how alpha diversity and other metadata changes across time.

08_longitudinal/pcoa volatility/weaning-{pcoa}-volatility.qzv : A volatility plot will let us look at patterns of variation along principle coordinate axes.

Feature volatility analysis
---------------------------

- Plots relative frequencies of features across states (only important features are plotted). A supervised learning regressor is used to identify important features.

`08 longitudinal/feature volatility/volatility_plot.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/08_longitudinal/feat_volatility/volatility_plot.qzv>`_

- Identifies features that are predictive of a numeric metadata column.

`08 longitudinal/feature volatility/accuracy_results.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/08_longitudinal/feat_volatility/accuracy_results.qzv>`_

- Abundance heatmap of the most important features in each sample or group. 

`08 longitudinal/feature volatility/important-feature-Age-days-heatmap.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/08_longitudinal/feat_volatility/important-feature-Age-days-heatmap.qzv>`_

`08 longitudinal/feature volatility/important-feature-Body-site-heatmap.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/08_longitudinal/feat_volatility/important-feature-Body-site-heatmap.qzv>`_

- Table of the most important features with their taxonomic affiliation.

`08 longitudinal/feature volatility/feature_importance.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/08_longitudinal/feat_volatility/feature_importance.qzv>`_

Regression accuracy results
---------------------------

- Scatter plot showing predicted vs. true values for each test sample, accompanied by a linear regression line fitted to the data with 95% confidence intervals (grey shading).

`08 longitudinal/regression/accuracy_results.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/08_longitudinal/regressor/accuracy_results.qzv>`_

- Table of the most important features with their taxonomic affiliation.

`08 longitudinal/regression/feature_importance.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/08_longitudinal/regressor/feature_importance.qzv>`_

Maturity Index prediction (MAZ)
-------------------------------

- Contains a linear regression plot of predicted vs. expected values on all control test samples

`08 longitudinal/maturity index/accuracy_results.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/08_longitudinal/maturity/accuracy_results.qzv>`_

- Interactive volitility chart. This visualization can be useful for assessing how MAZ and other metrics change over time in each sample group

`08 longitudinal/maturity index/volatility_plots.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/08_longitudinal/maturity/volatility_plots.qzv>`_

- Heatmap showing the frequency of each important feature across time in each group. This plot is useful for visualizing how the frequency of important features changes over time in each group, demonstrating how different patterns of feature abundance (e.g., trajectories of development in the case of age or time-based models) may affect model predictions and MAZ scores.

`08 longitudinal/maturity index/clustermap.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/08_longitudinal/maturity/clustermap.qzv>`_

- Table of the most important features with their taxonomic affiliation.

`08 longitudinal/maturity index/feature_importance.qzv </media/thomas/data/ciml_tomas_metab_weaning/05_Output/08_longitudinal/maturity/feature_importance.qzv>`_