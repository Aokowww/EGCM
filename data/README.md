# Local ADNI data

Keep each authorized snapshot in a dated directory such as:

`data/adni_raw/2026-07-26/`

Raw files and derived CSVs are ignored by Git and must not be committed,
redistributed, or uploaded to an external service. Record download dates,
versions and SHA-256 checksums in a local manifest without copying subject
data into that record.

Build the local 90- and 180-day aligned tables with:

```sh
Rscript code/adni/prepare_adni_redesign_data.R \
  data/adni_raw/2026-07-26 90 180
```

The analysis runner expects one row per PET-anchored `RID` and `EXAMDATE`.
Duplicate keys, unparseable dates, missing required columns, out-of-range
values, and date gaps outside the frozen window stop the analysis.
