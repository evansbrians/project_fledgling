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

### 1. Build the database  ·  `src/01_build_database/`

*Purpose:** Downloads detections from the CTT (Cellular Tracking
Technologies) API, builds the local DuckDB database, adds the tag
deployment/recovery dates, and shortens the long tag IDs.
**Run:** `setup.R` → `updating_db.R` → `read_deployment_recoveries.R` →
`shorten_ids.R`  (`functions.R` is a helper these load).

* **Uses:** `download_data()`.
* **Makes:** `data/oxbow.duckdb` (the database everything else reads from).

### 2. Calibration  ·  `src/02_calibration/`

Learns the relationship between radio signal strength (RSSI)
and distance, using Blair's transect walk where the true tag-to-node distances
are known. This is what lets us turn signal strengths into locations.

* **Run:** `calibration_blair.R`  (`deriving_data_blair.R` prepares its input).
* **Uses:** `calculate_node_locations()`, `get_calibration_error()`,
`get_summary_stats()`, `haversine()`.

### 3. Validation  ·  `src/03_validation/`

Checks the calibrated method against known ground-truth GPS points.

* **Run:** `groundtruth.R`  (`groundtruth_setup.R` prepares its input).
* **Uses:** `calculate_node_locations()`.

### 4. Filtering  ·  `src/04_filtering/`

Trims the raw detections down to what's usable. It keeps only detections **inside each tag's deployment window**, groups them into short time windows, keeps the strongest signal per node, and drops nodes too far from the strongest one.

* **Run:** `build_detection_windows.R`.
* **Uses:** `build_detection_windows()`.
* **Makes:** `data/localization/detection_windows.parquet`.

### 5. Localization  ·  `src/05_localization/`

Turns each time window into a single location (longitude /
latitude) by multilateration, and attaches a confidence and an error estimate to
every fix.

* **Run:** `localize_windows.R`.
* **Makes:** `data/localization/localizations.parquet`.

### 6. KDE (space-use surfaces)  ·  `src/06_kde/`

For each bird and each day, turns that day's locations into a
kernel density estimate — a smooth "heat map" of where the bird spent its time —
plus 50/80/95% home-range outlines.

* **Run:** `daily_kde.R`.
* **Makes:** daily rasters and polygons in `results/localization/kde_output/`.

### 7. KDE → fields  ·  `src/07_kde_to_fields/`

Overlays each day's heat map on the farm's field boundaries and
computes the share of use that fell in each field — the biology we're after
(which fields the birds actually used).

* **Run:** `field_use.R`.
* **Makes:** per-field use tables in `results/localization/field_use/`.

### 8. Movement states  ·  `src/08_movement_states/`

Classifies each window as **stationary** or **moving** from how
much the node signals vary, using a Bayesian hidden Markov model.

* **Run:** `detect_movement_states.R`  (`functions_movement_states.R` and the model
file `rssi_state_hmm.stan` sit alongside it).
* **Uses:** `fit_rssi_state_hmm()`.
* **Note:** the `.stan` file defines the model for detecting movment states and is fit in R via the `cmdstanr` package.

---

## Repository structure

| Folder | What's in it |
|---|---|
| `src/` | The eight workflow steps above. |
| `data/` | Inputs: the DuckDB database, the tag table, field boundaries. |
| `results/` | Outputs, written by the steps (locations, KDEs, field use). |
| `docs/` | An interactive field-usage map and its published web version. |

---
