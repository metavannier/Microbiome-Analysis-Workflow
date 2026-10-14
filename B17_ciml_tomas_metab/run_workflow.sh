#!/usr/bin/env bash
# ======================================================================
#  Run the microbiome workflow (B17_ciml_tomas_metab - undernutrition)
#  with the latest Snakemake executed through Apptainer.
#
#  Usage: ./run_workflow.sh <command> [comparison] [extra snakemake options]
#
#    list                  comparisons defined in comparisons.tsv
#    manifest              create manifest.tsv from the FASTQ of 00_RowData
#    dryrun [ID|all|none]  show the jobs that would be run (default: all)
#    core                  common part only (quality, DADA2, taxonomy, tree)
#    run ID [ID2,...]      run comparison(s) (+ the common part if needed),
#                          then write one report per comparison
#    all                   run every comparison, then one report each
#    report ID             (re)build the report of a comparison
#    unlock                remove the lock left by an interrupted run
#
#  Environment variables (defaults in brackets):
#    CORES            number of cores [8]
#    SDM              software deployment: conda | apptainer | none [conda]
#                     (none = the tools are already on the PATH)
#    SNAKEMAKE_IMAGE  Snakemake container [docker://snakemake/snakemake:v9.27.0]
#    SNAKEMAKE        run this Snakemake command instead of the container
#                     (e.g. SNAKEMAKE=snakemake if installed on the server)
#    REPORT_FORMAT    zip | html [zip]
#    BINDS            extra directories to bind in Apptainer, comma-separated
#                     (data, references or TMPDIR outside the project)
#
#  Examples:
#    CORES=16 ./run_workflow.sh run L1_Feces_Female_CTL_vs_SAM
#    CORES=16 ./run_workflow.sh all
#    SNAKEMAKE=snakemake SDM=apptainer ./run_workflow.sh run L3_Ileal_SAM_Female_vs_Male
# ======================================================================

set -euo pipefail

cd "$(dirname "$(readlink -f "$0")")"

CORES=${CORES:-8}
SDM=${SDM:-conda}
SNAKEMAKE_IMAGE=${SNAKEMAKE_IMAGE:-docker://snakemake/snakemake:v9.27.0}
REPORT_FORMAT=${REPORT_FORMAT:-zip}
BINDS=${BINDS:-}
PROJECT=$(awk '/^project:/ {p=1} p && $1 == "name:" {gsub("\"", "", $2); print $2; exit}' config.yaml)
REPORTDIR=reports

# The QIIME 2 conda environment is solved with flexible channel priority
export CONDA_CHANNEL_PRIORITY=flexible

bind_args() {
    local binds="$PWD"
    [ -n "$BINDS" ] && binds="$binds,$BINDS"
    echo "$binds"
}

snakemake_cmd() {
    if [ -n "${SNAKEMAKE:-}" ]; then
        # Snakemake installed on the server
        $SNAKEMAKE "$@"
    else
        if [ "$SDM" = "apptainer" ]; then
            echo "WARNING: SDM=apptainer from the Snakemake container needs" \
                 "Apptainer inside Apptainer; prefer SDM=conda, or install" \
                 "Snakemake on the server and set SNAKEMAKE=snakemake." >&2
        fi
        apptainer exec --bind "$(bind_args)" "$SNAKEMAKE_IMAGE" snakemake "$@"
    fi
}

deployment_args() {
    case "$SDM" in
        apptainer) echo "--sdm apptainer --apptainer-args=--bind=$(bind_args)" ;;
        conda) echo "--sdm conda" ;;
        none) echo "" ;;
        *) echo "SDM must be conda, apptainer or none" >&2; exit 1 ;;
    esac
}

comparisons() {
    awk -F'\t' 'NR > 1 && $1 !~ /^#/ && $1 != "" {print $1}' comparisons.tsv
}

check_comparison() {
    local id
    for id in ${1//,/ }; do
        if ! comparisons | grep -qx "$id"; then
            echo "Unknown comparison '$id'. Run: $0 list" >&2
            exit 1
        fi
    done
}

report_path() {
    echo "$REPORTDIR/${PROJECT}_$1_report.$REPORT_FORMAT"
}

make_report() {
    mkdir -p "$REPORTDIR"
    snakemake_cmd --cores 1 --config run_comparisons="$1" \
        --report "$(report_path "$1")" "${@:2}"
    echo "Report: $(report_path "$1")"
}

command=${1:-help}
shift || true

case "$command" in
    list)
        awk -F'\t' 'NR > 1 && $1 !~ /^#/ && $1 != "" {printf "%-32s %s\n", $1, $6}' comparisons.tsv
        ;;
    manifest)
        if [ -e manifest.tsv ]; then
            echo "manifest.tsv already exists (remove it to recreate it)" >&2
            exit 1
        fi
        python3 03_Script/make_manifest.py --metadata sample-metadata.tsv "$@"
        ;;
    dryrun)
        target=${1:-all}
        [ "$#" -gt 0 ] && shift
        snakemake_cmd --dry-run --cores "$CORES" $(deployment_args) \
            --config run_comparisons="$target" "$@"
        ;;
    core)
        snakemake_cmd --cores "$CORES" $(deployment_args) core "$@"
        ;;
    run)
        [ "$#" -ge 1 ] || { echo "Usage: $0 run <comparison ID>" >&2; exit 1; }
        ids=$1
        shift
        check_comparison "$ids"
        snakemake_cmd --cores "$CORES" $(deployment_args) \
            --config run_comparisons="$ids" "$@"
        for id in ${ids//,/ }; do make_report "$id"; done
        ;;
    all)
        snakemake_cmd --cores "$CORES" $(deployment_args) \
            --config run_comparisons=all "$@"
        for id in $(comparisons); do make_report "$id"; done
        ;;
    report)
        [ "$#" -ge 1 ] || { echo "Usage: $0 report <comparison ID>" >&2; exit 1; }
        check_comparison "$1"
        make_report "$@"
        ;;
    unlock)
        snakemake_cmd --cores 1 --unlock "$@"
        ;;
    help|-h|--help)
        sed -n '2,34p' "$0" | sed 's/^# \{0,1\}//'
        ;;
    *)
        echo "Unknown command '$command'" >&2
        sed -n '2,34p' "$0" | sed 's/^# \{0,1\}//'
        exit 1
        ;;
esac
