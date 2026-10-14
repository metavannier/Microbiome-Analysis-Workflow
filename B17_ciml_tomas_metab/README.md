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
- the species-level (level 7) taxon abundance tables formatted for
  **MetaboAnalyst**, generated automatically for every comparison;
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

## 2. Requirements

- **Apptainer 1.4.4** on the server, with Internet access to pull
  `docker://snakemake/snakemake:v9.27.0` and either the conda channels
  (default mode) or the container registries (`--sdm apptainer` mode).
- Disk: ~5 GB for the QIIME 2 environment/image, plus the results.
  QIIME 2 writes large temporary files: set `tmpdir:` in `config.yaml` to a
  large directory if `/tmp` is small.

Nothing else has to be installed: Snakemake itself runs in the official
Snakemake container, and every step runs in its own software environment.

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
     <https://www.arb-silva.de/current-release/QIIME2/2026.7/SSU/V3V4-341f-806r/uniform>:
     copy the `.qza` and write its **exact file name** in `config.yaml`
     (`taxonomy: classifier`). A classifier only works with the scikit-learn
     version it was trained with (1.7.1 for QIIME 2 2026.4/2026.7); check with
     `unzip -p 01_Reference/<classifier>.qza '*/data/sklearn_version.json'`;
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
   primer removal, DADA2, taxonomy, phylogeny). Look at
   `05_Output/02_qiime_import/undernutrition-PE-demux-noprimer.qzv` (read
   quality → DADA2 truncation lengths) and
   `05_Output/04_taxonomy/undernutrition-taxa-table-filtered-dada2.qzv`
   (reads per sample → `diversity: samplingdepth` and `maxdepth`), on
   <https://view.qiime2.org>. If you change the DADA2 parameters, rerun `core`.
2. Run the comparisons (`run <ID>` or `all`). The common part is reused; only
   the comparison is computed. Each comparison gets its own report
   `reports/undernutrition_<ID>_report.zip` (unzip it and open
   `report.html`; set `REPORT_FORMAT=html` for a single HTML file).

Snakemake only reruns what is needed: after a change of a parameter in
`config.yaml`, the steps using it (and the following ones) are rerun.

### Options of `run_workflow.sh` (environment variables)

| Variable | Default | |
|---|---|---|
| `CORES` | 8 | cores used by Snakemake |
| `SDM` | `conda` | `conda`: each step in a conda environment created by Snakemake (in `.snakemake/conda`, once); `apptainer`: each step in its container (needs Snakemake on the server, see below); `none`: tools already installed |
| `SNAKEMAKE_IMAGE` | `docker://snakemake/snakemake:v9.27.0` | Snakemake container |
| `SNAKEMAKE` | — | use this Snakemake instead of the container (e.g. `SNAKEMAKE=snakemake`) |
| `REPORT_FORMAT` | `zip` | `zip` or `html` |
| `BINDS` | — | directories outside the project to make visible in Apptainer (comma-separated), e.g. the target of symbolic links in `00_RowData`/`01_Reference`, or `tmpdir` |

Any other Snakemake option can be added at the end, e.g.
`./run_workflow.sh run L2_Ileal_CTL_vs_SAM --forcerun core_metrics`.

### Equivalent Snakemake commands

```bash
# Snakemake 9 in its container, each step in a conda environment (default)
apptainer exec --bind "$PWD" docker://snakemake/snakemake:v9.27.0 \
    snakemake --cores 16 --sdm conda --config run_comparisons=L1_Feces_Female_CTL_vs_SAM
apptainer exec --bind "$PWD" docker://snakemake/snakemake:v9.27.0 \
    snakemake --cores 1 --config run_comparisons=L1_Feces_Female_CTL_vs_SAM \
    --report reports/undernutrition_L1_Feces_Female_CTL_vs_SAM_report.zip

# Snakemake installed on the server, each step in its Apptainer container
snakemake --cores 16 --sdm apptainer --apptainer-args "--bind $PWD" \
    --config run_comparisons=L1_Feces_Female_CTL_vs_SAM
```

`--config run_comparisons=` takes a comparison ID, several IDs separated by
commas, `all` (default) or `none` (common part only).

### Software environments

Every step declares both a conda environment (`02_Container/`) and a
container (`config.yaml`, `containers:`):

| Steps | conda (`--sdm conda`) | container (`--sdm apptainer`) |
|---|---|---|
| QIIME 2 steps and helper scripts | `rachis-qiime2-2026.7-linux-64-conda.yml` (official QIIME 2 2026.7 environment) | `docker://quay.io/qiime2/qiime2:2026.7` |
| FastQC | `fastqc=0.12.1` | `docker://staphb/fastqc:0.12.1` |
| Trimmomatic | `trimmomatic=0.40` | `docker://staphb/trimmomatic:0.40` |
| MultiQC | `multiqc=1.35` | `docker://multiqc/multiqc:v1.35` |

When Snakemake itself runs inside the Snakemake container (default of
`run_workflow.sh`), use `SDM=conda`: running the step containers would need
Apptainer inside Apptainer. With Snakemake installed on the server
(`pip install snakemake` or `conda create -n snakemake -c conda-forge -c bioconda snakemake`),
`SNAKEMAKE=snakemake SDM=apptainer ./run_workflow.sh ...` uses the official
images instead (pulled once in `.snakemake/singularity`).

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
  the species level, frequency filtering by group, ANCOM, LEfSe input table;
- **MetaboAnalyst** (`03_metaboanalyst/`): species-level abundance table.

Results are in `05_Output/06_comparisons/<comparison>/`; the files of the
common part in `05_Output/00_*` to `05_Output/05_phylogeny`.

## 6. Species-level tables for MetaboAnalyst

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
  the total of the sample over all taxa), as in the previous tables;
- taxa absent from all the samples of the comparison are removed;
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

(`python3 03_Script/format_taxa_table.py --help` for all the options.)

## 7. Changes compared with the previous workflow

Methodology kept: Trimmomatic (`LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:25`),
cutadapt (same primers and options), DADA2 (same truncation, trimming,
expected errors, training reads), length filter ≥ 400 nt, naive Bayes
classifier (confidence 0.7), exclusion of `d__Eukaryota`, removal of the
samples below the sampling depth, SEPP tree, core metrics, alpha/beta group
significance, collapse at level 7, frequency filtering by group
(≥ 50 reads in ≥ 2 samples), pseudocount + ANCOM, LEfSe table.

Changes required by the current software:

- **DADA2 chimera removal: `consensus`** (DADA2 and QIIME 2 default) instead
  of `pooled`. The `pooled` option was removed from QIIME 2 in release 2025.4
  ([q2-dada2 #171](https://github.com/qiime2/q2-dada2/pull/171)); it is not
  accepted by any release that can use the 2026.7 SILVA classifier.
- QIIME 2 2026.7 commands: new required outputs (cutadapt statistics, DADA2
  base-transition statistics, feature-table summarize frequencies), shown in
  the reports.
- Alpha rarefaction: the maximum depth is lowered to the largest sample of a
  comparison when `maxdepth` is higher (QIIME 2 refuses it otherwise).
- SEPP: ASVs that cannot be inserted in the reference tree are removed before
  the phylogenetic metrics (`fragment-insertion filter-features`); UniFrac and
  Faith PD fail otherwise.
- Core metrics use a fixed random seed for the rarefaction
  (`diversity: random_seed`, `null` = random as before).
- Tool versions: FastQC 0.12.1, Trimmomatic 0.40, MultiQC 1.35 (snakemake
  wrappers 0.35.2 replaced by explicit commands).

Not carried over (not applicable to this design): longitudinal analyses,
alpha correlation and PCoA with a "days" axis (the metadata has no numeric
or time column); OTU clustering and Deblur alternatives (only DADA2 was
supported). The LEfSe formatting was ported from R to Python (same table;
the class row is named after the compared column instead of `day`). The
merge of the frequency-filtered tables, which listed the weaning files by
hand, now uses the groups of each comparison.

## 8. Troubleshooting

- **`No space left on device` in a QIIME 2 step**: set `tmpdir:` in
  `config.yaml` to a large directory (and `BINDS=<that directory>`).
- **Classifier error "scikit-learn version ... does not match"**: the
  classifier was not trained for QIIME 2 2026.7 (see §3).
- **Files not found inside the containers** (symbolic links, data elsewhere):
  `BINDS=/path/to/data,/path/to/references ./run_workflow.sh ...`.
- **Interrupted run**: `./run_workflow.sh unlock`, then run again (incomplete
  files are recomputed).
- **conda environment creation fails**: the QIIME 2 environment needs flexible
  channel priority (`run_workflow.sh` sets `CONDA_CHANNEL_PRIORITY=flexible`).
- **`apptainer pull` of `quay.io/qiime2/qiime2:2026.7` fails**: check the
  available tags on <https://quay.io/repository/qiime2/qiime2> and update
  `containers: qiime2` in `config.yaml`.

## Author

Thomas Vannier (@metavannier), <https://centuri-livingsystems.org/t-vannier/>
