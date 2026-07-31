# ADNI data access

The repository does not contain ADNI participant data. Reproducing the ADNI
analyses requires independent approval from ADNI and acceptance of its data-use
terms.

## Local preparation

1. Sign in to the LONI Image and Data Archive with an approved account.
2. Download the current ADNI source tables and their data dictionaries.
3. Record the download date, table versions and SHA-256 checksums in a local
   manifest.
4. Store the files under a dated local directory such as
   `data/adni_raw/YYYY-MM-DD/`.
5. Check the measurement dates and units against the source documentation.
6. Build the analytic tables locally.

For the MRI-AV45 data:

```bash
Rscript code/adni/prepare_adni_redesign_data.R \
  data/adni_raw/YYYY-MM-DD 90 180
```

For the MRI-ADAS13 data:

```bash
Rscript code/adni/prepare_adni_mri_adas13.R \
  data/adni_raw/YYYY-MM-DD 30 5
```

The analysis tables produced by these scripts remain under `data/` and are
ignored by Git. Serialized model objects and full result directories are also
ignored.

## Source documentation

- [ADNI getting started guide](https://adni.loni.usc.edu/quick-start-guide-asset/getting_started.html)
- [ADNI Study Files guide](https://adni.loni.usc.edu/quick-start-guide-asset/study_files.html)
- [ADNI data dictionary guide](https://adni.loni.usc.edu/quick-start-guide-asset101625/datadic.html)

The historical scripts use an older `adnimerge` object. They are retained for
reproduction of the thesis-era analysis, not as a substitute for the current
source tables.

Only the reviewed aggregate tables under `results_public/` are intended for
public distribution.
