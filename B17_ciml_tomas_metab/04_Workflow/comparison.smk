# ----------------------------------------------------------------------
# Comparisons (comparisons.tsv): samples of each comparison
# ----------------------------------------------------------------------
# Each comparison is a subset of the samples (e.g. Feces + Male) and one
# or several metadata columns to compare (e.g. Diet). All its results
# are written in 05_Output/06_comparisons/<comparison>/.

COMP_METADATA = f"{COMPDIR}/00_metadata/{{comparison}}-metadata.tsv"
COMP_PREFIX = f"{PROJ}-{{comparison}}"


rule comparison_metadata:
    input:
        metadata=METADATA_FILE,
        script=script("make_comparison_metadata.py"),
    output:
        COMP_METADATA,
    log:
        f"{LOGDIR}/06_comparisons/{{comparison}}/metadata.log",
    params:
        subset=lambda wildcards: comparison(wildcards)["subset"],
        columns=lambda wildcards: ",".join(dict.fromkeys(
            comparison(wildcards)["group_columns"]
            + [comparison(wildcards)["filter_column"]])),
        exclude=",".join(config["project"].get("exclude_samples") or []),
    conda:
        QIIME_ENV
    container:
        QIIME_CONTAINER
    shell:
        """
        python {input.script} \
            --metadata {input.metadata} \
            --subset '{params.subset}' \
            --group-columns '{params.columns}' \
            --exclude-samples '{params.exclude}' \
            --output {output} > {log} 2>&1
        """
