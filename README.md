# massmatcher

`massmatcher` is an R package for annotating untargeted metabolomics features by mass matching, isotopologue-aware matching, and adduct/isotope clustering.

It is designed for large metabolite reference databases and uses Arrow-backed parquet data for scalable filtering before precise matching.

## What the package does

- Matches unknown feature `m/z` values to metabolite candidates by:
  - Monoisotopic mass (`mass_match()`)
  - Most abundant isotopologue mass (`isotopologue_mass_match()`)
- Supports two reference databases for database-backed workflows:
  - `"metorigindb"` (default)
  - `"pubchem"`
- Adds optional chemical classification (`get_chemical_classification()`)
- Provides a higher-level clustering workflow that combines:
  - m/z matching
  - adduct correlation
  - isotope correlation
  (`mz_match_clustering()`, `mz_match_clustering_pubchem()`)

## Installation

Install from a local checkout:

```r
install.packages("remotes")
remotes::install_local("path/to/massmatcher")
```

Or from inside the package directory:

```r
install.packages("devtools")
devtools::install(".")
```

Load:

```r
library(massmatcher)
library(dplyr)
```

## Quick Start

Run a minimal monoisotopic mass-matching workflow with built-in example data:

```r
library(massmatcher)
library(dplyr)

data("feature_table_exp_hilicpos", package = "massmatcher")

# keep first two columns as mz/time
unknown_feature <- tibble(mz = c(160.1332, 176.1282), time = c(30, 35))

hits <- mass_match(
  unknown_feature = unknown_feature,
  mz_ppm = 10,
  adduct = "M+H",
  database = "metorigindb"
)

hits |> select(mz, time, theoretical_mz, mz_ppm, Name, Adduct) |> head(10)
```

Switch to PubChem by adding `database = "pubchem"` in `mass_match()`.

## Shiny app

`massmatcher` includes a bundled Shiny app with three major modules:
- Search (`mass_match()`, `isotopologue_mass_match()`, `target_mass_search()`)
- Clustering (`mz_match_clustering()`, `mz_match_clustering_pubchem()`)
- Additional Meta Data (`annotate_match_results()`, `annotate_mz_match_clustering_results()`)

In the Clustering module, selecting database `"pubchem"` automatically uses the lightweight
`mz_match_clustering_pubchem()` path.

Run locally:

```r
library(massmatcher)
run_massmatcher_app()
```

Search module input tips:
- Select `Ion mode` (positive/negative), then choose supported adduct(s) from the adduct dropdown.
- For `m/z` search and isotopologue search, enter one `mz` per line.
- For `target_mass_search`, enter target masses one per line. This searches monoisotopic masses within an uploaded feature table.
- For `target_mass_search`, upload a feature table where first column is `mz` and second column is `time`.
- For m/z or isotopologue search uploads, use a CSV/TSV with exactly one column: `mz`.

Additional Meta Data module:
- Upload a result table (CSV/TSV/RDS) to add classification and concentration metadata.
- The uploaded table should include `InChIKey` and/or `HMDB_ID` columns.

For deployment, the app source is under:
- `inst/shiny/massmatcher-app/app.R`

## Built-in data and resources

- Metabolite databases (parquet):
  - `inst/metorigindb/metorigindb_database_major_enriched.parquet`
  - `inst/pubchem/pubchem_clean-0.parquet`
- ClassyFire cache:
  - `inst/classyfire/pubchem_classyfire-0.parquet`
- Concentration resources:
  - `inst/concentration/hmdb_concentrations_normal_condition.csv`
  - `inst/concentration/custom_concentrations.csv`
- Example feature table:
  - `data/feature_table_exp_hilicpos.rda`

## Database selection

For database-backed functions, use `database = "metorigindb"` or `database = "pubchem"`.

Default is `"metorigindb"`:

```r
mass_match(
  unknown_feature = tibble(mz = c(160.1332, 176.1282), time = c(30, 35))
)
```

Switch to PubChem:

```r
mass_match(
  unknown_feature = tibble(mz = c(160.1332, 176.1282), time = c(30, 35)),
  database = "pubchem"
)
```

## Add classification and concentration without recomputation

Use enrichment helpers to add chemical classification and concentration to
outputs you already computed:

```r
mass_hits <- mass_match(
  unknown_feature = tibble(mz = c(160.1332, 176.1282), time = c(30, 35)),
  mz_ppm = 10,
  adduct = c("M+H", "M+Na"),
  database = "metorigindb"
)
mass_hits_enriched <- annotate_match_results(
  result_table = mass_hits,
  database = "metorigindb"
)

iso_hits <- isotopologue_mass_match(
  unknown_feature = tibble(mz = c(161.1366, 177.1315), time = c(30, 35)),
  mz_ppm = 10,
  adduct = c("M+H", "M+Na"),
  database = "metorigindb"
)
iso_hits_enriched <- annotate_match_results(
  result_table = iso_hits,
  database = "metorigindb"
)

data("feature_table_exp_hilicpos", package = "massmatcher")
clustering_results <- mz_match_clustering(
  met_raw_wide = feature_table_exp_hilicpos[1:500, 1:12],
  database = "metorigindb",
  write_output = FALSE
)
clustering_results_enriched <- annotate_mz_match_clustering_results(
  clustering_output = clustering_results,
  database = "metorigindb"
)
```

The enrichment columns include:
- `Classification_Kingdom`, `Classification_Superclass`, `Classification_Class`,
  `Classification_Subclass`, `Classification_Direct_parent`,
  `Classification_Alternative_parent`
- `Concentration_average`, `Concentration_units`

## Function overview and examples

### 1) Precise m/z matching by monoisotopic mass

`mass_match()` returns matched feature-metabolite pairs and `mz_ppm` error:

```r
## the first feature is Delta-Valerobetaine and the second feature is Homocarnitine
unknown_feature <- tibble(mz = c(160.1332, 176.1282), time = c(30, 35))

mass_hits <- mass_match(
  unknown_feature = unknown_feature,
  mz_ppm = 10,
  adduct = c("M+H"),
  database = "metorigindb"
)

head(mass_hits)
```

### 2) Isotopologue-aware m/z matching

`isotopologue_mass_match()` matches features using most abundant isotopologue masses:

```r
iso_hits <- isotopologue_mass_match(
  unknown_feature = tibble(mz = c(161.1366, 177.1315), time = c(30, 35)),
  mz_ppm = 10,
  adduct = c("M+H", "M+Na"),
  database = "metorigindb"
)
```

### 3) Targeted search for known masses

`target_mass_search()` finds expected adduct m/z values for one or more target masses in a feature table:

```r
feature_table <- tibble(
  mz = c(181.0707, 203.0526, 500.2000),
  time = c(100, 200, 300),
  sample_1 = c(10000, 5000, 3000)
)

target_hits <- target_mass_search(
  mass = c(180.063388, 202.045322),
  feature_table = feature_table,
  adduct = "M+H",
  mz_ppm = 10
)

target_hits
```

### 5) End-to-end m/z matching + adduct/isotope clustering

`mz_match_clustering()` is the high-level workflow that can write output CSV files and also returns result tables.

```r
data("feature_table_exp_hilicpos", package = "massmatcher")

# use a subset for a fast demo
feature_subset <- feature_table_exp_hilicpos[1:1000, 1:12]

clustering_results <- mz_match_clustering(
  met_raw_wide = feature_subset,
  database = "metorigindb",
  mz_threshold = 10,
  ion_mode = "positive",
  write_output = FALSE
)

names(clustering_results)
# "mz_only" "mz_only_isotope" "clustering_all" "clustering_main"
```

For PubChem, use the lightweight variant to reduce load time and memory while extracting PubChem results:

```r
feature_subset <- feature_table_exp_hilicpos[1:500, 1:12]
clustering_results_pubchem <- mz_match_clustering_pubchem(
  met_raw_wide = feature_subset,
  mz_threshold = 10,
  ion_mode = "positive",
  write_output = FALSE
)
```

### 6) Low-level overlap and enrichment utilities

`find.Overlapping.mzs()` is a general overlap helper used internally and can be used directly:

```r
find.Overlapping.mzs(
  dataA = data.frame(mz = c(100.0000, 200.0000), time = c(10, 20)),
  dataB = data.frame(mz = c(100.0002, 300.0000), time = c(12, 40)),
  mz.thresh = 5,
  time.thresh = 30
)
```

`annotate_match_results()` enriches a result table (for example from
`mass_match()` or `isotopologue_mass_match()`) with classification and
concentration columns:

```r
mass_hits_enriched <- annotate_match_results(
  result_table = mass_hits,
  database = "metorigindb",
  query_missing_classification = TRUE
)
```

`annotate_mz_match_clustering_results()` enriches existing
`mz_match_clustering()` outputs without rerunning clustering:

```r
clustering_results_enriched <- annotate_mz_match_clustering_results(
  clustering_output = clustering_results,
  database = "metorigindb",
  query_missing_classification = TRUE
)
```

`get_chemical_classification()` remains available when you only need
classification by InChIKey (without concentration joins).

## Typical workflow

1. Prepare feature table with first columns as `mz`, `time` and the rest as sample intensities.
2. Start with `mass_match()` and/or `isotopologue_mass_match()` for candidate annotation.
3. Add classification and concentration with `annotate_match_results()` and/or `annotate_mz_match_clustering_results()`.
4. Run `mz_match_clustering()` for metorigindb or `mz_match_clustering_pubchem()` for PubChem.
5. Use `target_mass_search()` for focused checks on specific metabolites.

## Notes

- Database-backed functions use `"metorigindb"` by default.
- Multiple adducts can be supplied as vectors (for example `c("M+H", "M+Na")`).
- `mass_match()` and `isotopologue_mass_match()` use a broad exact-mass database prefilter followed by exact ppm overlap matching.
- Matching quality depends on ion mode, adduct assumptions, and ppm threshold.
- The package bundles large parquet reference files; keep them inside the package structure for standalone use.

## Maintainer

Jiada (James) Zhan  
Email: `jzha832@emory.edu` or `jamesjiadazhan@gmail.com`
