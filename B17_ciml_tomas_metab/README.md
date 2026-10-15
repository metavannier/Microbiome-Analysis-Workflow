# B17_ciml_tomas_metab — Microbiome analysis workflow (project *undernutrition*)

16S rRNA gene (V3-V4) metabarcoding analysis from paired-end FASTQ files to
one interactive report **per comparison**, with [Snakemake](https://snakemake.github.io)
(9.x) run through [Apptainer](https://apptainer.org) (1.4.4) and
[QIIME 2](https://qiime2.org) **2026.7**.

This project is a reworked, standalone version of
[Microbiome-Analysis-Workflow](https://github.com/metavannier/Microbiome-Analysis-Workflow)
(Thomas Vannier, @metavannier). It keeps the same methodology (FastQC/Trimmomatic,
cutadapt, DADA2, SILVA naive Bayes classifier, SEPP, diversity, frequency
filtering per group, ANCOM, LEfSe table) and adds:

- comparisons defined in a table (`comparisons.tsv`), each analysed on its own
  subset of samples and producing its own report;
- the level-7 taxon abundance tables (taxa bar plot "Download CSV" at
  level 7) formatted for **MetaboAnalyst**, generated automatically for every
  comparison;
- containers / conda environments for every step, a single launcher
  (`run_workflow.sh`), configuration and sample-sheet validation.

---------------------------------------------------------------------------

## 1. Project layout

| Path | Content |
|---|---|
| `Snakefile` | entry point of the workflow |
| `config.yaml` | all the parameters (paths, primers, DADA2, taxonomy, depth, ...) |
| `comparisons.tsv` | the comparisons (one report each) |
| `sample-metadata.tsv` | QIIME 2 sample metadata (`Sexe`, `Body-site`, `Diet`, ...) |
| `manifest.tsv` | QIIME 2 manifest (to create, see §3) — `manifest.example.tsv` shows the format |
| `run_workflow.sh` | launcher (Snakemake through Apptainer) |
| `00_RowData/` | raw FASTQ files (`M6969_S2_L001_R1_001.fastq.gz`, ...) |
| `01_Reference/` | SILVA classifier and SEPP reference (`.qza`) |
| `02_Container/` | conda environments (QIIME 2 2026.7, FastQC, Trimmomatic, MultiQC) |
| `03_Script/` | helper scripts (manifest, comparison metadata, MetaboAnalyst and LEfSe formatting) |
| `04_Workflow/` | Snakemake rules |
| `05_Output/` | results (created by the workflow) |
| `06_Schemas/` | validation schemas of `config.yaml` and `comparisons.tsv` |
| `07_Report/` | texts of the reports |
| `08_Log/` | logs of every step (created by the workflow) |
| `profiles/default/` | default Snakemake options of this workflow |
| `reports/` | the reports, one per comparison (created by `run_workflow.sh`) |
| `.snakemake/` | Snakemake metadata, conda environments (`conda/`) and their package cache (`conda-pkgs/`) |

## 2. Requirements

- **Apptainer 1.4.4** on the server, with Internet access to pull
  `docker://snakemake/snakemake:v9.27.0` (once, into
  `02_Container/snakemake_v9.27.0.sif`) and either the conda channels
  (default `SDM=conda`) or the container registries (`SDM=apptainer`).
- Disk, in the project directory:
  - `SDM=conda`: ~10 GB for the conda environments (`.snakemake/conda`,
    QIIME 2 2026.7 ≈ 6 GB) and their package cache (`.snakemake/conda-pkgs`,
    mostly hard links; set by the Snakefile so that nothing goes to
    `~/.conda/pkgs`);
  - `SDM=apptainer`: the images (QIIME 2 ≈ 3-5 GB); converting the QIIME 2
    image needs ~6 GB of temporary space in `APPTAINER_TMPDIR` (default
    `$TMPDIR` or `/tmp`) and the layers are cached in `APPTAINER_CACHEDIR`
    (default `~/.apptainer/cache`). If `/tmp` or the home quota is small:
    `export APPTAINER_TMPDIR=/big/disk/apptainer-tmp APPTAINER_CACHEDIR=/big/disk/apptainer-cache`
    (`run_workflow.sh` creates, binds and forwards them);
  - plus the results; QIIME 2 also writes large temporary files: set
    `tmpdir:` in `config.yaml` to a large directory if `/tmp` is small.

Nothing else has to be installed: Snakemake itself runs in the official
Snakemake container, and every step runs in its own software environment
(`./run_workflow.sh manifest` uses `python3` of the server if there is one,
otherwise the Python of the Snakemake container).

## 3. Set up the analysis

1. **Raw reads** — copy the paired FASTQ files in `00_RowData/`
   (`<sample>_S<n>_L001_R1_001.fastq.gz` / `..._R2_001.fastq.gz`).

2. **Manifest** — create `manifest.tsv` (sample IDs = part of the file name
   before the first `_`, paths to the reads trimmed by the workflow):

   ```bash
   ./run_workflow.sh manifest     # = python3 03_Script/make_manifest.py --metadata sample-metadata.tsv
   ```

   or write it yourself in the same format:

   ```
   sample-id	forward-absolute-filepath	reverse-absolute-filepath
   M6969	$PWD/05_Output/00_trimmed/M6969_S2_L001_R1_001.fastq.gz	$PWD/05_Output/00_trimmed/M6969_S2_L001_R2_001.fastq.gz
   ```

   The samples analysed are those of the manifest; they must all be in
   `sample-metadata.tsv` (the workflow checks it).

3. **References** (`01_Reference/`, see its README):
   - the SILVA classifier for QIIME 2 **2026.7** downloaded from
     <https://www.arb-silva.de/current-release/QIIME2/2026.7/SSU/V3V4-341f-806r/uniform>,
     expected as `01_Reference/SILVA_144_SSURef_NR99_uniform_classifier_V3V4-341f-806r.qza`
     (`taxonomy: classifier` in `config.yaml`; adapt it if your file name
     differs). A classifier only works with the scikit-learn version it was
     trained with (1.7.1 for QIIME 2 2026.4/2026.7); check with
     `unzip -p 01_Reference/<classifier>.qza '*/data/sklearn_version.json'`.
     The SILVA site is a JavaScript browser: download with `wget` from
     `https://www.arb-silva.de/fileadmin/silva_databases/current/QIIME2/2026.7/SSU/V3V4-341f-806r/uniform/<file>`
     and keep the `.md5`/`.doi` files next to it;
   - `sepp-refs-silva-128.qza` for the SEPP phylogeny (`phylogeny: seppdb`).

4. **Parameters** — check `config.yaml`, in particular the DADA2 truncation
   lengths and the sampling depth (see §4, step 1).

## 4. Run

```bash
./run_workflow.sh list                                   # the comparisons
CORES=16 ./run_workflow.sh dryrun L1_Feces_Female_CTL_vs_SAM   # what will be run
CORES=16 ./run_workflow.sh core                          # common part only
CORES=16 ./run_workflow.sh run L1_Feces_Female_CTL_vs_SAM     # one comparison + its report
CORES=16 ./run_workflow.sh run L3_Feces_SAM_Female_vs_Male,L3_Ileal_SAM_Female_vs_Male
CORES=16 ./run_workflow.sh all                           # every comparison + one report each
./run_workflow.sh report L1_Feces_Female_CTL_vs_SAM      # rebuild a report
```

Recommended order:

1. `./run_workflow.sh core` runs the common part once (quality control,
   primer removal, DADA2, taxonomy, phylogeny). Look, on
   <https://view.qiime2.org>, at
   - `05_Output/02_qiime_import/undernutrition-PE-demux-noprimer.qzv` (read
     quality → DADA2 truncation lengths; if you change them, rerun `core`);
   - `05_Output/05_phylogeny/undernutrition-taxa-table-filtered-tree-sepp-dada2.qzv`
     (reads per sample of the table rarefied by the diversity analyses, i.e.
     without the ASVs SEPP could not insert → `diversity: samplingdepth` and
     `maxdepth`); `05_Output/04_taxonomy/undernutrition-taxa-table-filtered-dada2.qzv`
     is the same table with every ASV (used for the taxonomy, differential
     abundance and MetaboAnalyst tables, also filtered at `samplingdepth`).
2. Run the comparisons (`run <ID>` or `all`). The common part is reused; only
   the comparison is computed. Each comparison gets its own report
   `reports/undernutrition_<ID>_report.zip` (unzip it and open
   `report.html`; set `REPORT_FORMAT=html` for a single HTML file). If a
   comparison fails, the reports of the other ones are still written and the
   failed ones are listed.

Snakemake only reruns what is needed: after a change of a parameter in
`config.yaml`, the steps using it (and the following ones) are rerun
(`./run_workflow.sh dryrun <ID>` shows them). `./run_workflow.sh report <ID>`
refuses to build a report over results that are not up to date with
`config.yaml` (the report text describes the current parameters); run
`./run_workflow.sh run <ID>` instead, or force it with `FORCE_REPORT=1`.

Keep the same software deployment (`SDM`) for a project: switching between
`conda` and `apptainer` changes the software environment of every step, and
Snakemake then recomputes everything (to avoid it after checking the
dry-run, add `--rerun-triggers mtime params input code`).

### Options of `run_workflow.sh` (environment variables)

| Variable | Default | |
|---|---|---|
| `CORES` | 8 | cores used by Snakemake |
| `SDM` | `conda` | `conda`: each step in a conda environment created by Snakemake (in `.snakemake/conda`, once); `apptainer`: each step in its official container (see *Software environments*); `none`: tools already installed (only with `SNAKEMAKE=...`) |
| `SNAKEMAKE_IMAGE` | `docker://snakemake/snakemake:v9.27.0` | Snakemake container, pulled once into `02_Container/snakemake_<tag>.sif` (or the path of a local `.sif`) |
| `SNAKEMAKE` | — | use this Snakemake instead of the container (e.g. `SNAKEMAKE=snakemake`) |
| `REPORT_FORMAT` | `zip` | `zip` or `html` |
| `BINDS` | — | directories outside the project to make visible in Apptainer (comma-separated), e.g. the target of symbolic links in `00_RowData`/`01_Reference` (`tmpdir`, `APPTAINER_TMPDIR` and `APPTAINER_CACHEDIR` are added automatically) |
| `FORCE_REPORT` | — | `1`: build the report even if the results are not up to date |

Any other Snakemake option can be added at the end, e.g.
`./run_workflow.sh run L2_Ileal_CTL_vs_SAM --forcerun core_metrics`
(`--config KEY=VALUE` given there is merged with the comparison selection).

### Equivalent Snakemake commands

```bash
apptainer pull 02_Container/snakemake_v9.27.0.sif docker://snakemake/snakemake:v9.27.0   # once
mkdir -p reports

# Snakemake 9 in its container, each step in a conda environment (default)
apptainer exec --bind "$PWD" 02_Container/snakemake_v9.27.0.sif \
    snakemake --cores 16 --sdm conda --config run_comparisons=L1_Feces_Female_CTL_vs_SAM
apptainer exec --bind "$PWD" 02_Container/snakemake_v9.27.0.sif \
    snakemake --cores 1 --config run_comparisons=L1_Feces_Female_CTL_vs_SAM \
    --report reports/undernutrition_L1_Feces_Female_CTL_vs_SAM_report.zip

# Snakemake installed on the server, each step in its Apptainer container
snakemake --cores 16 --sdm apptainer --apptainer-args "--bind $PWD" \
    --config run_comparisons=L1_Feces_Female_CTL_vs_SAM
```

`--config run_comparisons=` takes a comparison ID, several IDs separated by
commas, `all` (default) or `none` (common part only); a plain
`snakemake --cores N --sdm conda` therefore runs every comparison. The
Snakefile sets flexible conda channel priority, which the QIIME 2
environment needs: ignore the Snakemake warning recommending strict channel
priorities and do not set `channel_priority: strict` for this project.

### Software environments

Every step declares both a conda environment (`02_Container/`) and a
container (`config.yaml`, `containers:`):

| Steps | conda (`--sdm conda`) | container (`--sdm apptainer`) |
|---|---|---|
| QIIME 2 steps and helper scripts | `rachis-qiime2-2026.7-linux-64-conda.yml` (official QIIME 2 2026.7 environment) | `docker://quay.io/qiime2/qiime2:2026.7` |
| FastQC | `fastqc=0.12.1` | `docker://staphb/fastqc:0.12.1` |
| Trimmomatic | `trimmomatic=0.40` | `docker://staphb/trimmomatic:0.40` |
| MultiQC | `multiqc=1.35` | `docker://multiqc/multiqc:v1.35` |

`SDM=conda` (default) is the safest when Snakemake runs in its container:
the environments are created once by conda inside the project.

`SDM=apptainer` runs every step in the official images (no environment to
build). From the Snakemake container it uses nested Apptainer (the image
contains Apptainer), which works when the server allows unprivileged user
namespaces and `/dev/fuse`; test it with
`apptainer exec docker://snakemake/snakemake:v9.27.0 apptainer exec docker://alpine:3.20 true`.
If this fails, use `SDM=conda`, or Snakemake installed on the server
(`conda create -n snakemake -c conda-forge -c bioconda snakemake=9.27.0`, or
`pip install "snakemake[reports]==9.27.0" pandas` — the workflow needs pandas
and the reports need pygments) with
`SNAKEMAKE=snakemake SDM=apptainer ./run_workflow.sh ...`. The step images
are pulled once (in `.snakemake/singularity`, or in `APPTAINER_CACHEDIR`
when it is set).

Since release 2026.4, the QIIME 2 "amplicon" distribution is called
**qiime2** (framework *rachis*): image `quay.io/qiime2/qiime2:2026.7`, conda
environment `rachis-qiime2-2026.7`. The `qiime` commands are unchanged.

## 5. Comparisons

`comparisons.tsv` (tab-separated) defines one comparison per row:

| Column | Example | Meaning |
|---|---|---|
| `comparison` | `L1_Feces_Female_CTL_vs_SAM` | identifier (letters, digits, `_ . + -`) |
| `level` | `1` | free label |
| `subset` | `Body-site=Feces;Sexe=Female` | samples kept: `COLUMN=VALUE[,VALUE]` conditions separated by `;` (empty = all) |
| `group_columns` | `Diet` | metadata column(s) compared, separated by `,` |
| `filter_column` | `Diet` | column whose groups are used for the frequency filtering of the differential abundance (empty = first group column) |
| `description` | `Level 1 - Feces, females: CTL vs SAM (Diet)` | shown on the report |

Predefined comparisons (12):

| Level | Comparisons | Samples | Compared |
|---|---|---|---|
| 1 — CTL vs SAM for each sex and body site | `L1_<Feces/Caecum/Ileal>_<Female/Male>_CTL_vs_SAM` | one body site, one sex | `Diet` |
| 2 — CTL vs SAM, both sexes, for each body site | `L2_<Feces/Caecum/Ileal>_CTL_vs_SAM` | one body site, the four conditions | `Diet` (CTL vs SAM) and `Name-experiment` (the four sex × diet conditions); frequency filtering by the four conditions |
| 3 — female SAM vs male SAM for each body site | `L3_<Feces/Caecum/Ileal>_SAM_Female_vs_Male` | one body site, SAM diet | `Sexe` |

For every comparison and compared column:

- **diversity** (`01_diversity/`): feature table summary, alpha rarefaction
  curves, taxa bar plots, core metrics (rarefaction at `samplingdepth`),
  Kruskal-Wallis on Shannon / observed features / evenness / Faith PD,
  PERMANOVA (global + pairwise) on Jaccard / Bray-Curtis / unweighted and
  weighted UniFrac, PCoA (Emperor);
- **differential abundance** (`02_differential_abundance/`): taxa collapsed at
  level 7, frequency filtering by group, ANCOM, LEfSe input table;
- **MetaboAnalyst** (`03_metaboanalyst/`): level-7 abundance table.

Results are in `05_Output/06_comparisons/<comparison>/`; the files of the
common part in `05_Output/00_*` to `05_Output/05_phylogeny`. The tree and
the table of the diversity analyses are named after the phylogeny method
(`...-rooted-tree-sepp.qza`): changing `phylogeny: method` builds a new tree
and reruns the diversity analyses.

## 6. Level-7 tables for MetaboAnalyst

> **Level 7 is the genus with the SILVA 144 classifiers.** The SILVA
> classifiers published for QIIME 2 2026.7 (SILVA release 144) have 7 ranks:
> domain, **kingdom** (new rank), phylum, class, order, family, genus — and
> SILVA states that the *uniform* classifiers do not include species. With
> the previous SILVA 138 classifier (your former `level-7.csv`,
> `d__;p__;c__;o__;f__;g__;s__`), level 7 was the species. The workflow
> writes `05_Output/04_taxonomy/taxonomy_ranks.tsv` (also in the report) with
> the rank of each level for your classifier. The level is set in
> `config.yaml` (`metaboanalyst: taxa_level`, `differential_abundance:
> taxa_level`); a species level needs a classifier that includes species
> (e.g. a SILVA 144 *weighted* classifier, which is habitat-specific, or a
> classifier trained with species labels for scikit-learn 1.7.1); check its
> levels in `taxonomy_ranks.tsv` and set `taxa_level` accordingly.

For each comparison and compared column, the workflow writes
`05_Output/06_comparisons/<ID>/03_metaboanalyst/undernutrition-<ID>-level-7-<column>-metaboanalyst.csv`:

```
Sample,Group,Unassigned;__;__;__;__;__;__,d__Bacteria;__;__;__;__;__;__,...
M6995,CTL,0,0,...
```

- source: `05_Output/04_taxonomy/undernutrition-level-7.csv`, extracted from
  `undernutrition-taxa-bar-plots.qzv` — the very file obtained with
  *Download CSV* at taxonomic level 7 in QIIME 2 View;
- samples of the comparison only, in rows; `Group` = value of the compared
  metadata column in `sample-metadata.tsv` (it replaces the second column of
  the former "Ileum-cleaned-F vs M.csv");
- values = relative abundance of each taxon in the sample (counts divided by
  the total of the sample over all taxa), as in the previous tables
  (default; `relative_abundance: false` gives the read counts);
- taxa absent from all the samples of the comparison are removed (default;
  `drop_zero_taxa: false` keeps them);
- comma-separated with `.` decimals (format read by MetaboAnalyst,
  *Statistical Analysis [one factor]*, samples in rows, unpaired).

Options in `config.yaml` (`metaboanalyst:`): taxonomic level, counts instead
of relative abundances, keep the absent taxa, prevalence filter
(`min_group_prevalence`), separator and decimal mark (`sep: ";"` and
`decimal: ","` reproduce the French Excel format of the old files).

The previous "cleaned" table (Ileum, family level) also lacked 48
low-abundance taxa; no simple rule reproduces exactly that selection (the
closest, "present in ≥ 3 samples of at least one group", matches 111 of 113
taxa), so no such filter is applied by default (`min_group_prevalence: 0`);
MetaboAnalyst offers its own data filtering.

The script can also be used by hand (Python 3 with pandas):

```bash
# with the metadata (new procedure)
python3 03_Script/format_taxa_table.py --taxa-csv level-7.csv \
    --metadata sample-metadata.tsv --group-column Name-experiment \
    --subset "Body-site=Ileal" --output Ileal-level-7.csv
# with an existing MetaboAnalyst table giving Sample/Group (old procedure)
python3 03_Script/format_taxa_table.py --taxa-csv level-7.csv \
    --groups-from "Ileum-cleaned-F vs M.csv" --output Ileal-level-7.csv
```

The taxa columns are recognised by their structure (same number of ranks as
the first taxon), so the metadata columns appended by QIIME 2 at the end of
`level-N.csv`, numeric or not, are never taken as taxa.

(`python3 03_Script/format_taxa_table.py --help` for all the options.)

## 7. Changes compared with the previous workflow

Methodology kept: Trimmomatic (`LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:25`),
cutadapt (same primers and options), DADA2 (same truncation, trimming,
expected errors, training reads), length filter ≥ 400 nt, naive Bayes
classification (confidence 0.7), exclusion of `d__Eukaryota`, removal of the
samples below the sampling depth, SEPP tree, core metrics, alpha/beta group
significance, collapse at level 7, frequency filtering by group
(≥ 50 reads in ≥ 2 samples), pseudocount + ANCOM, LEfSe table.

Changes required by the current software or by the new data:

- **DADA2 chimera removal: `consensus`** (DADA2 and QIIME 2 default) instead
  of `pooled`. The `pooled` option was removed from QIIME 2 in release 2025.4
  ([q2-dada2 #171](https://github.com/qiime2/q2-dada2/pull/171)); it is not
  accepted by any release that can use the 2026.7 SILVA classifier.
- **Classifier**: SILVA 138 99 % full-length naive Bayes classifier →
  SILVA 144 NR99 uniform classifier trained on the V3-V4 region (341f-806r)
  for QIIME 2 2026.7 (scikit-learn 1.7.1). Same confidence (0.7), but a new
  database release and region-specific training change the assignments, and
  the 7 ranks now include the kingdom: level 7 is the genus (see §6). Keep it
  in mind when comparing with earlier projects.
- QIIME 2 2026.7 commands: new required outputs (cutadapt statistics, DADA2
  base-transition statistics, feature-table summarize frequencies), shown in
  the reports.
- **SEPP**: the ASVs that cannot be inserted in the reference tree are removed
  (`fragment-insertion filter-features`) from the table used for **all** the
  diversity analyses (rarefaction curves, rarefaction at `samplingdepth`,
  every alpha and beta metric, also Shannon, observed features, evenness,
  Jaccard and Bray-Curtis): UniFrac and Faith PD fail otherwise, and
  core-metrics rarefies the table once for all the metrics. Taxa bar plots,
  ASV table, ANCOM, LEfSe and MetaboAnalyst tables keep every ASV. The removed
  ASVs are in `05_Output/05_phylogeny/undernutrition-table-not-in-tree-dada2.qza`.
  A sample just above `samplingdepth` in the taxonomy table can therefore fall
  below it in the diversity table: choose the depth with
  `undernutrition-taxa-table-filtered-tree-sepp-dada2.qzv` (§4).
- Alpha rarefaction: as each comparison is a subset of the samples, the
  maximum depth is lowered to the largest sample of the comparison when
  `maxdepth` is higher (QIIME 2, old and new, refuses a maximum depth above
  the largest sample).
- Core metrics use a fixed random seed for the rarefaction
  (`diversity: random_seed`, `null` = random as before).
- Tool versions: FastQC 0.12.1, Trimmomatic 0.40, MultiQC 1.35 (snakemake
  wrappers 0.35.2 replaced by explicit commands).

Not carried over (not applicable to this design): longitudinal analyses,
alpha correlation and PCoA with a "days" axis (the metadata has no numeric
or time column); OTU clustering and Deblur alternatives (only DADA2 was
supported); the relative-frequency table and the visualization of the
collapsed table before filtering (only used by those analyses).

Other minor differences:

- the LEfSe formatting was ported from R to Python (same table; the class
  row is named after the compared column instead of `day`);
- the merge of the frequency-filtered tables, which listed the weaning files
  by hand, now uses the groups of each comparison;
- samples are selected per comparison (and `exclude_samples` applied)
  *before* the frequency filtering by group, whereas the previous sample
  filter (`RMSAMPLE`) came after the merge;
- PERMANOVA is run on the compared columns of each comparison (the previous
  `conditions.tsv` targets were not used);
- the alpha rarefaction curves include Faith PD (the tree is given);
- added: summary of the imported reads before primer removal, cutadapt and
  DADA2 error-model visualizations, taxonomy in the ASV sequence table,
  per-comparison taxa bar plots, `taxonomy_ranks.tsv`, summary of the
  diversity table;
- `taxonomy_filtered.tsv` excludes the taxa case-insensitively, like
  `qiime taxa filter-table`;
- `maxdepth` guidance: about the median (not the mean) frequency per sample.

## 8. Troubleshooting

- **`No space left on device` in a QIIME 2 step**: set `tmpdir:` in
  `config.yaml` to a large directory (`run_workflow.sh` binds it).
- **`no space left on device` / `packer failed to pack` while pulling an
  image** (`SDM=apptainer`): set `APPTAINER_TMPDIR` and `APPTAINER_CACHEDIR`
  to a large directory (§2); `tmpdir:` does not apply to image pulls.
- **Disk quota exceeded in `$HOME`**: the conda package cache is
  `.snakemake/conda-pkgs` (set by the Snakefile unless `CONDA_PKGS_DIRS` is
  set); for Apptainer, set `APPTAINER_CACHEDIR` outside `$HOME`.
- **Classifier error "scikit-learn version ... does not match"**: the
  classifier was not trained for QIIME 2 2026.7 (see §3).
- **Files not found inside the containers** (symbolic links, data elsewhere):
  `BINDS=/path/to/data,/path/to/references ./run_workflow.sh ...`.
- **Interrupted run**: `./run_workflow.sh unlock`, then run again (incomplete
  files are recomputed).
- **conda environment creation fails to solve**: the QIIME 2 environment
  needs flexible channel priority; the Snakefile sets
  `CONDA_CHANNEL_PRIORITY=flexible` (do not override it).
- **Nested Apptainer fails** (`SDM=apptainer` from the Snakemake container):
  use `SDM=conda`, or Snakemake installed on the server (§4).
- **`apptainer pull` of `quay.io/qiime2/qiime2:2026.7` fails**: check the
  available tags on <https://quay.io/repository/qiime2/qiime2> and update
  `containers: qiime2` in `config.yaml`.

## Author

Thomas Vannier (@metavannier), <https://centuri-livingsystems.org/t-vannier/>
