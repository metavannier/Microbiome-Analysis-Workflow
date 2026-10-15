# 01_Reference

Reference databases (paths set in `config.yaml`):

- **Taxonomy classifier** (`taxonomy: classifier`): the SILVA (release 144)
  naive Bayes classifier pre-trained for QIIME 2 2026.7, V3-V4 region
  (341f-806r), uniform weights, expected as
  `SILVA_144_SSURef_NR99_uniform_classifier_V3V4-341f-806r.qza`, downloaded from
  <https://www.arb-silva.de/current-release/QIIME2/2026.7/SSU/V3V4-341f-806r/uniform>.
  Copy the `.qza` here and write its exact file name in `config.yaml`.
  A classifier must be used with the QIIME 2 release it was trained for
  (2026.7 here, the release used by the workflow).
- **SEPP reference** (`phylogeny: seppdb`, only with `method: "sepp"`):
  `sepp-refs-silva-128.qza` (QIIME 2 data resources).

If you use symbolic links pointing outside the project directory, bind
their target directory in Apptainer (README.md §4, variable `BINDS`, and §8).
