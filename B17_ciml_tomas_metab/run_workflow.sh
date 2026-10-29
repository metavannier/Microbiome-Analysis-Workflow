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
#    run ID[,ID2,...]      run comparison(s) (+ the common part if needed),
#                          then write one report per comparison
#    all                   run every comparison, then one report each
#    report ID             (re)build the report of an up-to-date comparison
#    unlock                remove the lock left by an interrupted run
#
#  Environment variables (defaults in brackets):
#    CORES            number of cores [8]
#    SDM              software deployment: conda | apptainer | none [conda]
#                     apptainer: official images (from the Snakemake
#                     container: nested Apptainer, see README)
#                     none: tools already installed (needs SNAKEMAKE=...)
#    SNAKEMAKE_IMAGE  Snakemake container [docker://snakemake/snakemake:v9.27.0]
#                     pulled once into 02_Container/snakemake_<tag>.sif
#    SNAKEMAKE        run this Snakemake command instead of the container
#                     (e.g. SNAKEMAKE=snakemake if installed on the server)
#    REPORT_FORMAT    zip | html [zip]
#    BINDS            extra directories to bind in Apptainer, comma-separated
#                     (data or references outside the project)
#    FORCE_REPORT=1   build a report even if the results are not up to date
#    APPTAINER_TMPDIR, APPTAINER_CACHEDIR  (optional) Apptainer temporary and
#                     cache directories, also used inside the container
#
#  Extra Snakemake options can be added at the end of the command, e.g.
#    ./run_workflow.sh run L2_Ileal_CTL_vs_SAM --forcerun core_metrics
#
#  Examples:
#    CORES=16 ./run_workflow.sh run L1_Feces_Female_CTL_vs_SAM
#    CORES=16 ./run_workflow.sh all
#    SDM=apptainer ./run_workflow.sh run L3_Ileal_SAM_Female_vs_Male
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

case "$SDM" in
    conda) DEPLOY=(--sdm conda) ;;
    apptainer) DEPLOY=(--sdm apptainer) ;;
    none)
        DEPLOY=()
        if [ -z "${SNAKEMAKE:-}" ]; then
            echo "SDM=none needs Snakemake installed on the server" \
                 "(SNAKEMAKE=snakemake): the Snakemake container does not" \
                 "see the tools of the server." >&2
            exit 1
        fi ;;
    *) echo "SDM must be conda, apptainer or none (not '$SDM')" >&2; exit 1 ;;
esac

# Directories visible in the containers: the project, BINDS, the
# Apptainer temporary/cache directories and the 'tmpdir' of config.yaml
bind_args() {
    local binds="$PWD" dir
    [ -n "$BINDS" ] && binds="$binds,$BINDS"
    for dir in "${APPTAINER_TMPDIR:-}" "${APPTAINER_CACHEDIR:-}" \
        "$(awk -F'"' '/^tmpdir:/ {print $2}' config.yaml)"; do
        case "$dir" in
            /*) mkdir -p "$dir"; binds="$binds,$dir" ;;
        esac
    done
    echo "$binds"
}
if [ "$SDM" = "apptainer" ]; then
    DEPLOY+=("--apptainer-args=--bind $(bind_args)")
fi

# Snakemake container image, pulled once (no registry access afterwards)
snakemake_image() {
    case "$SNAKEMAKE_IMAGE" in
        *://*)
            local sif="02_Container/snakemake_${SNAKEMAKE_IMAGE##*:}.sif"
            if [ ! -e "$sif" ]; then
                echo "Pulling $SNAKEMAKE_IMAGE into $sif" >&2
                apptainer pull "$sif" "$SNAKEMAKE_IMAGE" >&2
            fi
            echo "$sif" ;;
        *) echo "$SNAKEMAKE_IMAGE" ;;
    esac
}

snakemake_cmd() {
    if [ -n "${SNAKEMAKE:-}" ]; then
        # Snakemake installed on the server
        $SNAKEMAKE "$@"
    else
        # Apptainer does not pass APPTAINER_* variables into the container
        local envs=() var
        for var in APPTAINER_TMPDIR APPTAINER_CACHEDIR; do
            if [ -n "${!var:-}" ]; then envs+=(--env "$var=${!var}"); fi
        done
        apptainer exec ${envs[@]+"${envs[@]}"} --bind "$(bind_args)" \
            "$(snakemake_image)" snakemake "$@"
    fi
}

# Snakemake keeps only the last --config given: the --config KEY=VALUE
# of the extra options are merged with run_comparisons.
EXTRA_CONFIG=()
EXTRA_ARGS=()
split_extra() {
    local arg in_config=false
    for arg in "$@"; do
        case "$arg" in
            --config|-C) in_config=true ;;
            --config=*) EXTRA_CONFIG+=("${arg#--config=}"); in_config=false ;;
            -*) in_config=false; EXTRA_ARGS+=("$arg") ;;
            *) if $in_config; then EXTRA_CONFIG+=("$arg"); else EXTRA_ARGS+=("$arg"); fi ;;
        esac
    done
}

# snakemake_sel <run_comparisons> <snakemake options...>
snakemake_sel() {
    local selection=$1
    shift
    snakemake_cmd "$@" ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"} \
        --config run_comparisons="$selection" ${EXTRA_CONFIG[@]+"${EXTRA_CONFIG[@]}"}
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
    snakemake_cmd --cores 1 --report "$(report_path "$1")" \
        --config run_comparisons="$1" ${EXTRA_CONFIG[@]+"${EXTRA_CONFIG[@]}"} &&
        echo "Report: $(report_path "$1")"
}

# Results of a comparison complete and up to date (nothing left to run)
is_complete() {
    snakemake_sel "$1" --dry-run --cores 1 ${DEPLOY[@]+"${DEPLOY[@]}"} 2>&1 |
        grep -q "Nothing to be done"
}

# Reports of the comparisons of a run, even if some of them failed: after
# a failure, only the comparisons whose results are complete get a report
make_reports() {
    local status=$1 id missing=()
    shift
    for id in "$@"; do
        if [ "$status" -ne 0 ] && ! is_complete "$id"; then
            missing+=("$id")
            continue
        fi
        make_report "$id" || missing+=("$id")
    done
    if [ "${#missing[@]}" -gt 0 ]; then
        echo "No report for: ${missing[*]} (failed or incomplete results;" \
             "see the errors above, fix them and run again)." >&2
        [ "$status" -ne 0 ] || status=1
    fi
    return "$status"
}

# The report text is rendered from the current config.yaml: refuse to
# build it over results computed with other parameters, code or inputs.
check_up_to_date() {
    local out
    [ "${FORCE_REPORT:-0}" = 1 ] && return 0
    out=$(snakemake_sel "$1" --dry-run --cores 1 ${DEPLOY[@]+"${DEPLOY[@]}"} 2>&1) ||
        { echo "$out" >&2; exit 1; }
    if ! grep -q "Nothing to be done" <<< "$out"; then
        echo "$out" | sed -n '/^Job stats:/,/^total/p' >&2
        echo "ERROR: the results of '$1' are not up to date (config.yaml," \
             "code or input changed). Run: $0 run $1" \
             "(or FORCE_REPORT=1 $0 report $1 to build it anyway)." >&2
        exit 1
    fi
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
        if command -v python3 > /dev/null || [ -n "${SNAKEMAKE:-}" ]; then
            python3 03_Script/make_manifest.py --metadata sample-metadata.tsv "$@"
        else
            apptainer exec --bind "$(bind_args)" "$(snakemake_image)" \
                python3 03_Script/make_manifest.py --metadata sample-metadata.tsv "$@"
        fi
        ;;
    dryrun)
        target=all
        if [ "$#" -gt 0 ] && [ "${1#-}" = "$1" ]; then target=$1; shift; fi
        split_extra "$@"
        snakemake_sel "$target" --dry-run --cores "$CORES" ${DEPLOY[@]+"${DEPLOY[@]}"}
        ;;
    core)
        split_extra "$@"
        snakemake_sel none --cores "$CORES" ${DEPLOY[@]+"${DEPLOY[@]}"} core
        ;;
    run)
        [ "$#" -ge 1 ] || { echo "Usage: $0 run <comparison ID>" >&2; exit 1; }
        ids=$1
        shift
        check_comparison "$ids"
        split_extra "$@"
        status=0
        snakemake_sel "$ids" --cores "$CORES" ${DEPLOY[@]+"${DEPLOY[@]}"} || status=$?
        make_reports "$status" ${ids//,/ }
        ;;
    all)
        split_extra "$@"
        status=0
        snakemake_sel all --cores "$CORES" ${DEPLOY[@]+"${DEPLOY[@]}"} || status=$?
        make_reports "$status" $(comparisons)
        ;;
    report)
        [ "$#" -ge 1 ] || { echo "Usage: $0 report <comparison ID>" >&2; exit 1; }
        id=$1
        shift
        check_comparison "$id"
        split_extra "$@"
        check_up_to_date "$id"
        make_report "$id"
        ;;
    unlock)
        snakemake_cmd --cores 1 --unlock "$@"
        ;;
    help|-h|--help)
        sed -n '2,42p' "$0" | sed 's/^# \{0,1\}//'
        ;;
    *)
        echo "Unknown command '$command'" >&2
        sed -n '2,42p' "$0" | sed 's/^# \{0,1\}//'
        exit 1
        ;;
esac
