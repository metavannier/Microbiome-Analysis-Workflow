# 01_Reference

Reference databases (paths set in `config.yaml`):

- **Taxonomy classifier** (`taxonomy: classifier`): the SILVA (release 144)
  naive Bayes classifier pre-trained for QIIME 2 2026.7, V3-V4 region
  (341f-806r), habitat-weighted "animal distal gut", 8 ranks down to the
  species:
  `SILVA_144_SSURef_NR99_weighted_classifier_V3V4-341f-806r_animal-distal-gut.qza`
  (394,721,014 bytes, md5 `570fc11cb9ecdfb2b64e87bf48fc5c98`, DOI
  10.82364/144/SSU/Ref-NR99/QIIME2/2026.7/V3V4-341f-806r/weighted/animal-distal-gut),
  from <https://www.arb-silva.de/current-release/QIIME2/2026.7/SSU/V3V4-341f-806r/weighted/animal-distal-gut>.
- **Control classifier** (`taxonomy: check_classifier`, optional): the SILVA
  144 uniform classifier of the same region, down to the genus, used to
  cross-check the genera of the weighted classifier (SILVA recommendation):
  `SILVA_144_SSURef_NR99_uniform_classifier_V3V4-341f-806r.qza`
  (112,685,152 bytes, DOI 10.82364/144/SSU/Ref-NR99/QIIME2/2026.7/V3V4-341f-806r/uniform),
  from <https://www.arb-silva.de/current-release/QIIME2/2026.7/SSU/V3V4-341f-806r/uniform>.
- **SEPP reference** (`phylogeny: seppdb`, only with `method: "sepp"`):
  `sepp-refs-silva-128.qza` (QIIME 2 data resources).

Download (the SILVA site is a JavaScript browser; the files themselves are
under `fileadmin`):

```bash
cd 01_Reference
B=https://www.arb-silva.de/fileadmin/silva_databases/current/QIIME2/2026.7/SSU/V3V4-341f-806r
wget $B/weighted/animal-distal-gut/SILVA_144_SSURef_NR99_weighted_classifier_V3V4-341f-806r_animal-distal-gut.qza{,.md5,.doi}
wget $B/uniform/SILVA_144_SSURef_NR99_uniform_classifier_V3V4-341f-806r.qza{,.md5,.doi}
md5sum *.qza; cat *.md5      # the checksums must be identical
```

Keep the SILVA file names and the `.doi` files (DOI to cite). A classifier
must be used with the QIIME 2 release it was trained for (2026.7 here, the
release used by the workflow; scikit-learn 1.7.1). The weighted classifier
needs ~24-32 GB of RAM and ~40 GB of temporary space (README.md §2). Why
these files and what to expect from the species: README.md §6.1.

If you use symbolic links pointing outside the project directory, bind
their target directory in Apptainer (README.md §4, variable `BINDS`, and §8).
