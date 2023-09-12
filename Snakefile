import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")

##### load config and sample sheets #####
configfile: "config.yaml"
validate(config, schema="06_Schemas/config.schema.yaml")

## Get the different group for the feature filtering by abundance
group = pd.read_table(config["group"]).set_index(["group"], drop=False)
GROUPCONDITION = expand("{group.group}", group=group.itertuples())

## Get the column name of the metadata file
condition = pd.read_table(config["condition"]).set_index(["condition"], drop=False)
CONDITION = expand("{condition.condition}", condition=condition.itertuples())

##### Set variables ####
ROWDATA = srcdir("00_RowData")
REF = srcdir("01_Reference")
SCRIPTSDIR = srcdir("03_Script")
ENVDIR = srcdir("04_Workflow")
ROOTDIR = os.getcwd()
LOGDIR = srcdir("08_Log")
OUTPUTDIR = srcdir("05_Output")

PROJ = config["project"]["name"]
SUF = config["project"]["suffix"]
R1_SUF = str(config["project"]["r1_suf"])
R2_SUF = str(config["project"]["r2_suf"])
METADATA = config["project"]["metadata"]
RMSAMPLE = config["project"]["rmsample"]


# Use glob statement to find all samples in 'raw_data' directory
## Wildcard '{num}' must be equivalent to 'R1' or '1', meaning the read pair designation.
SAMPLE_LIST,NUMS = glob_wildcards(ROWDATA + "/{sample}_L001_{num}" + SUF)
# # Unique the output variables from glob_wildcards
SAMPLE_SET = set(SAMPLE_LIST)
SET_NUMS = set(NUMS)

# Grouping method : name of the plugin use for denoising or clustering
GROUP = config["grouping_method"]["group"]

# Database information to assign taxonomy
DB_classifier = config["taxonomy"]["database_classified"]

# Method and database for the phylogeny
PHYLO = config["phylogeny"]["method"]
DB_sepp = config["phylogeny"]["seppdb"]

# ANCOM: column with the group to compare
COLUMN = config["ancom"]["metadata_column"].split(',')

# Alpha diversity metrics
ALPHADIV = config["volatility"]["alpha_diversity"].split(',')
BETADIV =  config["volatility"]["beta_diversity"].split(',')
PCOA = config["volatility"]["pcoa"].split(',')

# Column to present on the heatmap
HEATMAP = config["volatility"]["features_column"].split(',')

# Feature for the pairwise comparison
FEATURE = config["volatility"]["feature"].split(',')

# ----------------------------------------------
# Target rules
# ----------------------------------------------

rule all:
  input:
    # fastqc output before trimming
    # raw_html = expand("{outputdir}/01_fastqc/{sample}_L001_{num}_fastqc.html", outputdir=OUTPUTDIR, sample=SAMPLE_SET, num=SET_NUMS),
    # raw_zip = expand("{outputdir}/01_fastqc/{sample}_L001_{num}_fastqc.zip", outputdir=OUTPUTDIR, sample=SAMPLE_SET, num=SET_NUMS),
    # raw_multi_html = OUTPUTDIR + "/01_fastqc/raw_multiqc.html",
    # # Trimmed data output
    # trimmedData = expand("{outputdir}/00_trimmed/{sample}_L001_{num}_001.fastq.gz", outputdir=OUTPUTDIR, sample=SAMPLE_SET, num=SET_NUMS),
    # trim_html = expand("{outputdir}/01_fastqc/{sample}_L001_{num}_trimmed_fastqc.html", outputdir=OUTPUTDIR, sample=SAMPLE_SET, num=SET_NUMS),
    # trim_zip = expand("{outputdir}/01_fastqc/{sample}_L001_{num}_trimmed_fastqc.zip", outputdir=OUTPUTDIR, sample=SAMPLE_SET, num=SET_NUMS),
    # trim_multi_html = OUTPUTDIR + "/01_fastqc/trimmed_multiqc.html", #next change to include proj name
    # # ####Qiime import
    # q2_import = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-demux-paired-end.qza",
    # # Qiime remove primer
    # q2_primerRM = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux-noprimer.qza",
    # primer = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux-noprimer.qzv",
    # ########### Choose the method for grouping similar sequences ###########
    # ################### Clustering OTU ###################
    # # # Merge with vsearch
    # # q2_joined = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined.qza",
    # # # Quality filter
    # # q2_filtered = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined-filtered.qza",
    # # q2_filterstats = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined-filtered-STATS.qza",
    # # # Quality visualization
    # # raw = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux.qzv",
    # # primer = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux-noprimer.qzv",
    # # q2_joined_sum = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined.qzv",
    # # q2_filtered_sum = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined-filtered.qzv",
    # # # Dereplicate sequence
    # # q2_dedup = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined-filtered-dedup_table.qza",
    # # q2_dedup_seq = OUTPUTDIR + "/03_otu/" + PROJ + "-PE-demux-joined-filtered-dedup_seqs.qza",
    # # # De novo clustering
    # # q2_cluster_table = OUTPUTDIR + "/03_otu/" + PROJ + "-cluster-table.qza",
    # # q2_cluster_seqs = OUTPUTDIR + "/03_otu/" + PROJ + "-cluster-seqs.qza",
    # # # Detect chimera
    # # q2_chimeras = OUTPUTDIR + "/03_otu/" + PROJ + "-chimeras.qza",
    # # q2_nonchimeras = OUTPUTDIR + "/03_otu/" + PROJ + "-nonchimeras.qza",
    # # q2_chimeras_stats = OUTPUTDIR + "/03_otu/" + PROJ + "-chimeras-STATS.qza",
    # # # Filter chimera
    # # q2_table_nc = OUTPUTDIR + "/03_otu/" + PROJ + "-table-nc.qza",
    # # q2_seqs_nc = OUTPUTDIR + "/03_otu/" + PROJ + "-rep-seqs-otu.qza
    # # q2_table_nc_qzv = OUTPUTDIR + "/03_otu/" + PROJ + "-table-nc.qzv",
    # # table_biom_otu = OUTPUTDIR + "/03_otu/feature-table.biom",  
    # # table_tsv_otu = OUTPUTDIR + "/03_otu/" + PROJ + "-otu-table.tsv",
    # ################### Denoising with Deblur ###################
    # # # # Merge
    # # q2_joined = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-PE-demux-joined.qza",
    # # # Quality filter
    # # q2_filtered = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-PE-demux-filtered.qza",
    # # q2_filterstats = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-PE-demux-filtered-STATS.qza",
    # # # Deblur
    # # q2_repseq = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-rep-seqs-" + GROUP + ".qza",
    # # q2_deblurtab = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-table-" + GROUP + ".qza",
    # # q2_deblurstat = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-" + GROUP + "-stats.qza",
    # # # Quality
    # # raw = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux.qzv",
    # # primer = OUTPUTDIR + "/02_qiime_import/" + PROJ + "-PE-demux-noprimer.qzv",
    # # q2_joined_sum = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-PE-demux-joined.qzv",
    # # q2_filtered_sum = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-PE-demux-filtered-STATS.qzv",
    # # q2_deblurstat_sum = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-" + GROUP + "-stats.qzv",
    # # q2_repseq_sum = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-rep-seqs-" + GROUP + ".qzv",
    # # table_otu = OUTPUTDIR + "/03_" + GROUP + "/feature-table.biom",
    # # otu_table = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-otu-table.tsv",
    # ################### Denoising with Dada2 ###################
    # table = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-table-" + GROUP + ".qza",
    # rep = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-rep-seqs-" + GROUP + ".qza",
    # stats = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-dada2-stats.qza",
    # filterrep = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-rep-filtered-seqs-dada2.qza",
    # filtertable = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-table-filtered-" + GROUP + ".qza",
    # # # Taxonomy
    # sklearn = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-tax_sklearn.qza",
    # taxafiltertable = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-taxa-table-filtered-" + GROUP + ".qza",
    # q2_repseq_filtered = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-rep-filtered-seqs-taxa-" + GROUP + ".qza",
    # table_tax = OUTPUTDIR + "/04_taxonomy/taxonomy.tsv",
    # table_tax_filtered = OUTPUTDIR + "/04_taxonomy/taxonomy_filtered.tsv",
    # stats_viz = OUTPUTDIR + "/03_" + GROUP + "/" + PROJ + "-dada2-stats.qzv",
    # rep_viz = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-rep-filtered-seqs-taxa-" + GROUP + ".qzv",
    # featurestat = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-taxa-table-filtered-" + GROUP + ".qzv",
    # rarefactionfiltertable = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-rarefaction-table-filtered-" + GROUP + ".qza",
    # relativefreqtable = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-relative-frequency-" + GROUP + ".qza",
    # table_biom = OUTPUTDIR + "/04_taxonomy/feature-table.biom",
    # taxo_table_biom = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-asv-table-with-taxonomy.biom",
    # taxo_table_tsv = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-asv-table-with-taxonomy.tsv", 
    # taxabarplots = OUTPUTDIR + "/04_taxonomy/" + PROJ + "-taxa-bar-plots.qzv",
    # # Phylogeny with fasttree
    # # aligned_repseq = OUTPUTDIR + "/05_phylogeny/" + PROJ + "-aligned-rep-filtered-seqs.qza",
    # # masked_aligned_repseq = OUTPUTDIR + "/05_phylogeny/" + PROJ + "-masked-aligned-rep-filtered-seqs.qza",  
    # # unrooted_tree = OUTPUTDIR + "/05_phylogeny/" + PROJ + "-unrooted-tree.qza",
    # # rooted_tree = OUTPUTDIR + "/05_phylogeny/" + PROJ + "-rooted-tree.qza",
    # # output_phylseq = OUTPUTDIR + "/05_phylogeny/" + PROJ + "-phylseq.tsv",
    # # Phylogeny with SEPP
    # tree = OUTPUTDIR + "/05_phylogeny/" + PROJ + "-rooted-tree.qza",
    # insertion = OUTPUTDIR + "/05_phylogeny/" + PROJ + "-insertion-placements.qza",
    # ## Diversity
    # ### Only if you need to remove sample
    # filtertable_selectedsample = OUTPUTDIR + "/06_diversity/" + PROJ + "-taxa-table-filtered-selectedsample-" + GROUP + ".qza",
    # table_count_qzv = OUTPUTDIR + "/06_diversity/" + PROJ + "-taxa-table-filtered-selectedsample-" + GROUP + ".qzv",
    # rarefaction = OUTPUTDIR + "/06_diversity/" + PROJ + "-alpha_rarefaction_curves.qzv",
    # d1 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/rarefied_table.qza",
    # d2 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/faith_pd_vector.qza",
    # d3 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/observed_features_vector.qza",
    # d4 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/shannon_vector.qza",
    # d5 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/evenness_vector.qza",
    # d6 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/unweighted_unifrac_distance_matrix.qza",
    # d7 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/weighted_unifrac_distance_matrix.qza",
    # d8 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/jaccard_distance_matrix.qza",
    # d9 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/bray_curtis_distance_matrix.qza",
    # d11 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/unweighted_unifrac_pcoa_results.qza",
    # d12 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/weighted_unifrac_pcoa_results.qza",
    # d13 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/jaccard_pcoa_results.qza",
    # d14 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/bray_curtis_pcoa_results.qza",
    # d15 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/unweighted_unifrac_emperor.qzv",
    # d16 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/weighted_unifrac_emperor.qzv",
    # d17 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/jaccard_emperor.qzv",
    # d18 = OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/bray_curtis_emperor.qzv",    
    # alphasigni = expand(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/{alphadiv}-group-significance.qzv", alphadiv=ALPHADIV),
    # ## Correlation (need numeric values)
    # coroutput = expand(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/{alphadiv}-correlation.qzv", alphadiv=ALPHADIV),
    # site = expand(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/{condition}-{betadiv}-significance.qzv", condition=CONDITION, betadiv=BETADIV),
    # ## PCOa with cinetic in the axis
    # pcooutput = expand(OUTPUTDIR + "/06_diversity/" + PROJ + "-core-metrics-results/{pcoa}-emperor-days.qzv", pcoa=PCOA),
    # # Differential abundance
    table_collapse = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-collapse-table-" + GROUP + ".qza",
    table_collapse_viz = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-collapse-table-" + GROUP + ".qzv",
    output_table_split = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-table_split-" + GROUP + ".txt",
    output_table_filtered = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-table_filtered-" + GROUP + ".txt",
    output_table_merged = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-table_merged-" + GROUP + ".txt",
    output_visualization_feature_table = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-table_merged-" + GROUP + ".qzv",
    ## To do if you want remove samples for the differential analyses (with ANCOM)
    # table_abond_selectedsample = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-table-abund-selectedsample-" + GROUP + ".qza",
    # table_abond_qzv = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-table-abund-selectedsample-" + GROUP + ".qzv",
    # table_abond_comp = OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-table-abund-comp-" + GROUP + ".qza",
    # ancom = expand(OUTPUTDIR + "/07_differential_abundance/" + PROJ + "-ancom-{column}-" + GROUP + ".qzv", column=COLUMN),
    # #### Longitudinal analysis (NEED numerical data)
    # output_pairwise_difference = expand(OUTPUTDIR + "/08_longitudinal/" + PROJ + "-{feature}-pairwise-differences.qzv", feature=FEATURE),
    # output_pairwise_distance = expand(OUTPUTDIR + "/08_longitudinal/" + PROJ + "-{betadiv}-pairwise-distances.qzv", betadiv=BETADIV),
    # alphavolatility = expand(OUTPUTDIR + "/08_longitudinal/" + PROJ + "-{alphadiv}-volatility.qzv", alphadiv=ALPHADIV),
    # pcoaoutput = expand(OUTPUTDIR + "/08_longitudinal/" + PROJ + "-{pcoa}_pcoa_results.qzv", pcoa=PCOA),
    # featlong1 = OUTPUTDIR + "/08_longitudinal/feat_volatility/filtered_table.qza",
    # featlong2 = OUTPUTDIR + "/08_longitudinal/feat_volatility/sample_estimator.qza",
    # featlong3 = OUTPUTDIR + "/08_longitudinal/feat_volatility/feature_importance.qza",
    # featlong4 = OUTPUTDIR + "/08_longitudinal/feat_volatility/accuracy_results.qzv",
    # featlong5 = OUTPUTDIR + "/08_longitudinal/feat_volatility/volatility_plot.qzv",
    # featlong6 = OUTPUTDIR + "/08_longitudinal/feat_volatility/feature_importance.qza",
    # important_feature_table_top = expand(OUTPUTDIR + "/08_longitudinal/feat_volatility/important-feature-table-top-{heatmap}.qza", heatmap=HEATMAP),
    # feature_heatmap = expand(OUTPUTDIR + "/08_longitudinal/feat_volatility/important-feature-{heatmap}-heatmap.qzv", heatmap=HEATMAP), 
    # ## Predicting continuous (i.e., numerical) sample data
    # regressor1 = OUTPUTDIR + "/08_longitudinal/regressor/sample_estimator.qza",
    # regressor2 = OUTPUTDIR + "/08_longitudinal/regressor/feature_importance.qza",
    # regressor3 = OUTPUTDIR + "/08_longitudinal/regressor/predictions.qza",
    # regressor4 = OUTPUTDIR + "/08_longitudinal/regressor/accuracy_results.qzv",
    # regressor5 = OUTPUTDIR + "/08_longitudinal/regressor/model_summary.qzv",
    # regressor6 = OUTPUTDIR + "/08_longitudinal/regressor/feature_importance.qzv",
    # ## Maturity Index prediction
    # maturity1 = OUTPUTDIR + "/08_longitudinal/maturity/maz_scores.qza",
    # maturity2 = OUTPUTDIR + "/08_longitudinal/maturity/sample_estimator.qza",
    # maturity3 = OUTPUTDIR + "/08_longitudinal/maturity/feature_importance.qza",
    # maturity4 = OUTPUTDIR + "/08_longitudinal/maturity/predictions.qza",
    # maturity5 = OUTPUTDIR + "/08_longitudinal/maturity/accuracy_results.qzv",
    # maturity6 = OUTPUTDIR + "/08_longitudinal/maturity/volatility_plots.qzv",
    # maturity7 = OUTPUTDIR + "/08_longitudinal/maturity/clustermap.qzv",
    # maturity8 = OUTPUTDIR + "/08_longitudinal/maturity/model_summary.qzv",
    # maturity9 = OUTPUTDIR + "/08_longitudinal/maturity/feature_importance.qzv",

# ----------------------------------------------
# setup singularity 
# ----------------------------------------------

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
#singularity: "docker://continuumio/miniconda3"

# ----------------------------------------------
# setup report
# ----------------------------------------------

report: "07_Report/workflow.rst"

# ----------------------------------------------
# Impose rule order for the execution of the workflow 
# ----------------------------------------------

#ruleorder: fastqc > trimmomatic_pe > fastqc_trimmed > multiqc > qiime_import > rm_primers > dada2 > filterseq > filterfeature > dada2_stats > gen_table > convert > assign_tax > get_stats_tax > gen_tax

# ----------------------------------------------
# Load rules 
# ----------------------------------------------

# include: "04_Workflow/quality.smk"
# include: "04_Workflow/qiime.smk"
# include: "04_Workflow/" + GROUP + ".smk"
# include: "04_Workflow/taxonomy.smk"
# include: "04_Workflow/phylogeny_" + PHYLO + ".smk"
# include: "04_Workflow/diversity.smk"
include: "04_Workflow/differential_abundance.smk"
include: "04_Workflow/longitudinal_analysis.smk"