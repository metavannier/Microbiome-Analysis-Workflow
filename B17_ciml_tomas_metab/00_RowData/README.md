# 00_RowData

Put here the raw paired-end reads of the project, e.g.

    M6969_S2_L001_R1_001.fastq.gz
    M6969_S2_L001_R2_001.fastq.gz

then create `manifest.tsv` at the root of the project:

    python3 03_Script/make_manifest.py --metadata sample-metadata.tsv
