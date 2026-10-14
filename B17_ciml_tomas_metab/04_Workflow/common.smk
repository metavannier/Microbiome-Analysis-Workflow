# ----------------------------------------------------------------------
# Common definitions: paths, sample sheets, comparisons and helpers
# ----------------------------------------------------------------------

import mimetypes
import os
import re
import sys

import pandas as pd
from snakemake.exceptions import WorkflowError
from snakemake.utils import validate


# ---- Directories (relative to the project directory) -------------------

ROWDATA = config["project"]["raw_dir"]
ENVDIR = "02_Container"
SCRIPTSDIR = "03_Script"
OUTPUTDIR = "05_Output"
REPORTDIR = "07_Report"
LOGDIR = "08_Log"

TRIMDIR = f"{OUTPUTDIR}/00_trimmed"
UNPAIREDDIR = f"{OUTPUTDIR}/00_unpaired"
QCDIR = f"{OUTPUTDIR}/01_fastqc"
IMPORTDIR = f"{OUTPUTDIR}/02_qiime_import"
DADA2DIR = f"{OUTPUTDIR}/03_dada2"
TAXDIR = f"{OUTPUTDIR}/04_taxonomy"
PHYLODIR = f"{OUTPUTDIR}/05_phylogeny"
COMPDIR = f"{OUTPUTDIR}/06_comparisons/{{comparison}}"

PROJ = config["project"]["name"]
GROUP = "dada2"  # denoising method, kept in the file names
R1_SUF = config["project"]["r1_suffix"]
R2_SUF = config["project"]["r2_suffix"]
FASTQ_EXT = ".fastq.gz"


# ---- Software environments --------------------------------------------
# Each step declares a conda environment and a container image; the one
# used depends on --software-deployment-method (--sdm) conda | apptainer.

QIIME_CONTAINER = config["containers"]["qiime2"]


# conda environment files and report captions are given relative to the
# rule files (04_Workflow/), hence the '../'
def conda_env(name):
    return f"../{ENVDIR}/{name}"


def caption(name):
    return f"../{REPORTDIR}/{name}.rst"


# QIIME 2 results are zip archives: offered for download in the report
mimetypes.add_type("application/zip", ".qzv")
mimetypes.add_type("application/zip", ".qza")


QIIME_ENV = conda_env("rachis-qiime2-2026.7-linux-64-conda.yml")


def script(name):
    return f"{SCRIPTSDIR}/{name}"


# Shell prefix: strict mode, UTF-8 locale (required by the QIIME 2 CLI)
# and an optional temporary directory with enough space for QIIME 2.
_prefix = "set -euo pipefail; export LC_ALL=C.UTF-8 LANG=C.UTF-8; "
if config.get("tmpdir"):
    _prefix += "mkdir -p {0}; export TMPDIR={0}; ".format(config["tmpdir"])
shell.prefix(_prefix)


# ---- Sample metadata ---------------------------------------------------

def read_qiime_metadata(path):
    """QIIME 2 metadata as strings, indexed by sample ID (directive rows
    such as '#q2:types' and blank rows removed)."""
    if not os.path.exists(path):
        raise WorkflowError(f"Metadata file '{path}' not found.")
    md = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    id_col = md.columns[0]
    md[id_col] = md[id_col].str.strip()
    md = md[(md[id_col] != "") & (~md[id_col].str.startswith("#"))]
    if md[id_col].duplicated().any():
        dup = md[id_col][md[id_col].duplicated()].tolist()
        raise WorkflowError(f"Duplicated sample IDs in {path}: {dup}")
    return md.set_index(id_col)


METADATA_FILE = config["metadata"]
METADATA = read_qiime_metadata(METADATA_FILE)


# ---- Manifest (samples to analyse) -------------------------------------

def read_manifest(path):
    """Read manifest.tsv. The FASTQ paths must point to the reads trimmed
    by the workflow (05_Output/00_trimmed/<raw file name>); '$PWD' and
    other environment variables are expanded as QIIME 2 does."""
    if not os.path.exists(path):
        raise WorkflowError(
            f"Manifest '{path}' not found. Put the FASTQ files in "
            f"{ROWDATA}/ and create it with:\n"
            f"    python3 {SCRIPTSDIR}/make_manifest.py "
            f"--metadata {METADATA_FILE}\n"
            "(see manifest.example.tsv for the format)")
    mf = pd.read_csv(path, sep="\t", dtype=str, comment="#",
                     keep_default_na=False)
    expected = ["sample-id", "forward-absolute-filepath",
                "reverse-absolute-filepath"]
    if list(mf.columns[:3]) != expected:
        raise WorkflowError(f"{path} must have the header: "
                            + "\t".join(expected))
    mf = mf[mf["sample-id"].str.strip() != ""]
    if mf["sample-id"].duplicated().any():
        raise WorkflowError(f"Duplicated sample IDs in {path}")
    trimdir = os.path.abspath(TRIMDIR)
    fastqs = {}
    for _, row in mf.iterrows():
        sample = row["sample-id"].strip()
        for col in expected[1:]:
            fq = os.path.expandvars(row[col].strip().replace(
                "$PWD", os.getcwd()))
            if os.path.dirname(os.path.abspath(fq)) != trimdir:
                raise WorkflowError(
                    f"{path}: '{row[col]}' (sample {sample}) must be in "
                    f"$PWD/{TRIMDIR}/ (reads trimmed by the workflow).")
            name = os.path.basename(fq)
            if not name.endswith(FASTQ_EXT):
                raise WorkflowError(f"{path}: {name} is not a .fastq.gz")
            fastqs.setdefault(sample, []).append(name)
    return fastqs


MANIFEST_FILE = config["manifest"]
SAMPLE_FASTQS = read_manifest(MANIFEST_FILE)
SAMPLES = list(SAMPLE_FASTQS)

_no_metadata = [s for s in SAMPLES if s not in METADATA.index]
if _no_metadata:
    raise WorkflowError(
        f"Samples of {MANIFEST_FILE} missing from {METADATA_FILE}: "
        + ", ".join(_no_metadata))
_not_in_manifest = [s for s in METADATA.index if s not in SAMPLE_FASTQS]
if _not_in_manifest:
    print(f"WARNING: {len(_not_in_manifest)} samples of {METADATA_FILE} "
          f"are not in {MANIFEST_FILE} and will not be analysed: "
          + ", ".join(_not_in_manifest[:20])
          + (" ..." if len(_not_in_manifest) > 20 else ""),
          file=sys.stderr)

# Read pairs: R1 file name without the R1 suffix -> (R1 name, R2 name)
READ_PAIRS = {}
for _sample, (_r1, _r2) in SAMPLE_FASTQS.items():
    if not (_r1.endswith(R1_SUF) and _r2.endswith(R2_SUF)
            and _r1[: -len(R1_SUF)] == _r2[: -len(R2_SUF)]):
        raise WorkflowError(
            f"Sample {_sample}: '{_r1}' / '{_r2}' do not follow the "
            f"read suffixes '{R1_SUF}' / '{R2_SUF}' of config.yaml")
    READ_PAIRS[_r1[: -len(R1_SUF)]] = (_r1, _r2)
FASTQ_NAMES = [name for pair in READ_PAIRS.values() for name in pair]
FASTQ_STEMS = [name[: -len(FASTQ_EXT)] for name in FASTQ_NAMES]


# ---- Comparisons -------------------------------------------------------

SAFE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.+-]*$")


def parse_subset(text):
    """'Body-site=Feces;Sexe=Male' -> {'Body-site': ['Feces'], ...}"""
    conditions = {}
    for cond in filter(None, (c.strip() for c in str(text).split(";"))):
        if "=" not in cond:
            raise WorkflowError(f"Malformed subset condition '{cond}' "
                                "(expected COLUMN=VALUE[,VALUE])")
        col, values = cond.split("=", 1)
        conditions[col.strip()] = [v.strip() for v in values.split(",")]
    return conditions


def load_comparisons(path):
    """Read and check comparisons.tsv against the analysed samples."""
    sheet = pd.read_csv(path, sep="\t", dtype=str, comment="#",
                        keep_default_na=False)
    validate(sheet, schema="../06_Schemas/comparisons.schema.yaml")
    sheet["comparison"] = sheet["comparison"].str.strip()
    if sheet["comparison"].duplicated().any():
        raise WorkflowError(f"Duplicated comparison IDs in {path}")
    excluded = set(config["project"].get("exclude_samples") or [])
    analysed = METADATA.loc[[s for s in SAMPLES if s not in excluded]]
    comparisons = {}
    for _, row in sheet.iterrows():
        cid = row["comparison"]
        subset = parse_subset(row["subset"])
        group_columns = [c.strip() for c in row["group_columns"].split(",")
                         if c.strip()]
        filter_column = row["filter_column"].strip() or group_columns[0]
        for col in list(subset) + group_columns + [filter_column]:
            if col not in METADATA.columns:
                raise WorkflowError(
                    f"Comparison {cid}: column '{col}' not in "
                    f"{METADATA_FILE} ({', '.join(METADATA.columns)})")
            if not SAFE.match(col):
                raise WorkflowError(
                    f"Comparison {cid}: column name '{col}' may only "
                    "contain letters, digits and _ . + -")
        keep = pd.Series(True, index=analysed.index)
        for col, values in subset.items():
            keep &= analysed[col].str.strip().isin(values)
        samples = analysed[keep]
        if samples.empty:
            raise WorkflowError(
                f"Comparison {cid}: no sample of {MANIFEST_FILE} matches "
                f"'{row['subset']}'")
        groups = {}
        for col in dict.fromkeys(group_columns + [filter_column]):
            values = sorted(v for v in samples[col].str.strip().unique() if v)
            if len(values) < 2:
                raise WorkflowError(
                    f"Comparison {cid}: column '{col}' has fewer than 2 "
                    f"groups among the samples of {MANIFEST_FILE} matching "
                    f"'{row['subset']}': {values}")
            bad = [v for v in values if not SAFE.match(v)]
            if bad:
                raise WorkflowError(
                    f"Comparison {cid}: values of '{col}' must only "
                    f"contain letters, digits and _ . + - ({bad})")
            groups[col] = values
        comparisons[cid] = {
            "level": row["level"],
            "description": row["description"],
            "subset": row["subset"],
            "group_columns": group_columns,
            "filter_column": filter_column,
            "groups": groups,
            "n_samples": len(samples),
        }
    return comparisons


COMPARISONS = load_comparisons(config["comparisons"])


def selected_comparisons():
    """config 'run_comparisons': 'all', 'none' or IDs (comma-separated
    string or YAML list), e.g. --config run_comparisons=L1_Feces_Male"""
    sel = config.get("run_comparisons", "all")
    if isinstance(sel, (list, tuple)):
        ids = [str(s).strip() for s in sel]
    else:
        sel = str(sel).strip()
        if sel.lower() == "all":
            return list(COMPARISONS)
        if sel.lower() in ("none", ""):
            return []
        ids = [s.strip() for s in sel.split(",") if s.strip()]
    unknown = [s for s in ids if s not in COMPARISONS]
    if unknown:
        raise WorkflowError(
            f"Unknown comparison(s) {unknown}. Defined in "
            f"{config['comparisons']}: {', '.join(COMPARISONS)}")
    return ids


def comparison(wildcards):
    return COMPARISONS[wildcards.comparison]


def alternatives(values):
    """Regular expression matching exactly one of the values."""
    return "|".join(re.escape(v) for v in values) or "NO_MATCH"


COMPARED_COLUMNS = sorted({c for comp in COMPARISONS.values()
                           for c in comp["group_columns"]})


# Comparisons described on the first page of the report (07_Report/workflow.rst)
config["report_comparisons"] = [
    {
        "id": cid,
        "description": COMPARISONS[cid]["description"],
        "subset": COMPARISONS[cid]["subset"],
        "group_columns": ", ".join(COMPARISONS[cid]["group_columns"]),
        "filter_column": COMPARISONS[cid]["filter_column"],
        "n_samples": COMPARISONS[cid]["n_samples"],
    }
    for cid in selected_comparisons()
]
