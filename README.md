# Project Fledgling: Data management and analysis

The code in this repository processes the radio-telemetry detections from the Oxbow Farm node array, stores the data in a duck databes, localizes detections, generates daily KDEs, calculates the daily usage of each field, and determines movement states (stationary vs. moving). The analysis is sequential (each step relies on the previous) and lives in `src/`.

---

## Before you start

**Install the companion package, `oxbowR`.** To keep things simple, I created an R package with all of the functions. Install the package from my git repo with:

```r
pak::pak("evansbrians/project_fledgling_package")
```

... and attach the library with:

```r
library(oxbowR)
```

**Database note**: The detection timestamps are timezone-aware, so
DuckDB needs its `icu` extension for the database steps. The first time you run
them you may hit an "`icu` extension not found" error — fix it once, with any
open DuckDB connection `con`:

```r
DBI::dbExecute(con, "INSTALL icu")
```

It downloads and caches the extension -- you only have to do this once (
it loads automatically in future sessions).

---

## Process

Run them in this order. "Uses" lists the main `oxbowR` functions each step
calls; the rest of each script is ordinary tidyverse/spatial R.

### 1. Build the database

*Purpose:** Downloads detections from the CTT (Cellular Tracking
Technologies) API, builds the local DuckDB database, adds the tag
deployment/recovery dates, and shortens the long tag IDs. **Warning!** This step takes a long time and will eat your internet! **Do not run this if you already have the database!**

* **Run:** `1.1_build_database_setup.R` → `1.2_build_database_updating_db.R` → `1.3_build_database_read_deployment_recoveries.R` → `1.4_build_database_shorten_ids.R`.
* **Uses:** `download_data()`.
* **Makes:** `data/oxbow.duckdb` (the database everything else reads from).

### 2. Calibration

Learns the relationship between radio signal strength (RSSI)
and distance, using Blair's transect walk where the true tag-to-node distances
are known. This is what lets us turn localize based on signal strength.

* **Run:** `2.1_calibration_deriving_data_blair.R` (prepares the input) → `2.2_calibration_blair.R`.
* **Uses:** `calculate_node_locations()`, `get_calibration_error()`,
`get_summary_stats()`, `haversine()`.
* **Makes:** the RSSI→distance calibration table,
`data/calibration/cal_rssi_v_dist_blair.csv` (read by steps 3 and 5).

### 3. Validation

Checks the calibrated method against known ground-truth GPS points.

* **Run:** `3.1_validation_groundtruth_setup.R` (prepares the input) → `3.2_validation_groundtruth.R`.
* **Uses:** `calculate_node_locations()`.

### 4. Filtering

Trims the raw detections down to what's usable. It keeps only detections **inside each tag's deployment window**, groups them into short time windows, keeps the strongest signal per node, and drops nodes too far from the strongest one.

* **Run:** `4_build_detection_windows.R`.
* **Uses:** `build_detection_windows()`.
* **Makes:** `data/localization/detection_windows.parquet`.

### 5. Localization

Turns each time window into a single location (longitude /
latitude) by multilateration, and attaches a confidence and an error estimate to
every fix.

* **Run:** `5_localize_windows.R`.
* **Makes:** `data/localization/localizations.parquet`.

### 6. KDE (space-use surfaces)

For each bird and each day, turns that day's locations into a
kernel density estimate — a smooth "heat map" of where the bird spent its time —
plus 50/80/95% home-range outlines.

* **Run:** `6_daily_kde.R`.
* **Makes:** daily rasters and polygons in `results/localization/kde_output/`.

### 7. KDE → fields

Overlays each day's heat map on the farm's field boundaries and
computes the share of use that fell in each field — the biology we're after
(which fields the birds actually used).

* **Run:** `7_field_use.R`.
* **Makes:** per-field use tables in `results/localization/field_use/`.

### 8. Movement states

Classifies each window as **stationary** or **moving** from how
much the node signals vary, using a Bayesian hidden Markov model.

* **Run:** `8_detect_movement_states.R`.
* **Uses:** `fit_rssi_state_hmm()`.
* **Note:** the Stan model ships inside the `oxbowR` package and is fit from R via the `cmdstanr` package.

---

## Repository structure

| Folder | What's in it |
|---|---|
| `src/` | The eight workflow steps above, as flat numbered files (steps 1–3 have sub-steps, e.g. `1.1`, `1.2`). |
| `data/` | Inputs: the DuckDB database, the tag table, field boundaries. |
| `results/` | Outputs, written by the steps (locations, KDEs, field use). |
| `docs/` | An interactive field-usage map and its published web version. |

---
