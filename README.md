
```markdown
# TCR Analysis Pipeline

This repository contains scripts and notebooks for analyzing T-cell receptor (TCR) sequencing data, including distance calculations, clustering, diversity, and overlap analyses for both alpha and beta chains. The pipeline is organized to allow stepwise processing from raw data to summary outputs.

---

## 1. `setup_results_directory.sh`

### Purpose
Creates the directory structure required to store TCR analysis outputs for both **alpha** and **beta chains**.

### Directories Created

**Alpha chain:**
```

outputs/alpha/tcrdist/
├── alignments/img           # Plots of sequence alignments
├── cluster\_overlaps         # Heatmaps and Jaccard index outputs
├── cluster\_sizes            # Cluster summary statistics
└── diversity                # Diversity metrics (Shannon, Simpson, etc.)

```

**Beta chain:**
```

outputs/beta/tcrdist/
├── alignments/img
├── cluster\_overlaps
├── cluster\_sizes
└── diversity

````

### Usage
```bash
chmod +x setup_results_directory.sh
./setup_results_directory.sh
````

* `mkdir -p` ensures the script doesn’t fail if directories already exist.
* Must be run **before running any R or Python TCR analysis scripts**.

---

## 2. `TCR analysis.R`

### Purpose

Prepares the TCR data for downstream analysis.

* Reads in processed TCR sequencing CSV files.
* Cleans and standardizes the column names for V/J genes and CDR3 sequences.
* Outputs tidy datasets ready for distance and clustering analyses.

---

## 3. `TCRdist.ipynb`

### Purpose

Generates TCR-specific distance matrices using the [`tcrdist3`](https://tcrdist3.readthedocs.io/en/latest/tcrdistances.html) library in Python.

### Features

* Handles both **alpha** and **beta chains**.
* Computes generation probabilities (`pgen`) using OLGA models.
* Saves distance matrices and probability tables for each input file.
* Outputs compressed results for download (`tcrdist.tar.gz`).

### Usage

1. Upload CSV files for the chain of interest.
2. Create a folder for output:

   ```bash
   ! mkdir tcrdist
   ```
3. Install dependencies:

   ```python
   !pip install tcrdist3
   ```
4. Run the notebook cells for **alpha** or **beta chains**.
5. Download the `tcrdist.tar.gz` results file.

---

## 4. `TCR clustering.R`

### Purpose

Performs clustering of TCR sequences based on distance metrics.

* Uses `TCRdist` or other similarity measures to group TCRs.
* Outputs cluster summaries including sizes and members.

---

## 5. `TCR diversity.R`

### Purpose

Calculates diversity metrics for TCR repertoires.

* Computes Shannon and Simpson diversity indices.
* Can normalize by read depth or other size factors.
* Produces plots summarizing diversity changes across samples or timepoints.

---

## 6. `TCR overlap.R`

### Purpose

Analyzes overlap of TCR clusters across samples or mice.

* Creates **presence/absence heatmaps** for clusters and amino acid clonotypes.
* Computes **Jaccard indices** to quantify similarity between repertoires.
* Saves both raw data and plots to the `cluster_overlaps` directories.

---

## Notes

* Scripts are designed to work sequentially. Recommended workflow:

  1. Run `setup_results_directory.sh` to create directories.
  2. Prepare and clean data with `TCR analysis.R`.
  3. Generate TCR distances with `TCRdist.ipynb`.
  4. Cluster TCRs with `TCR clustering.R`.
  5. Calculate diversity with `TCR diversity.R`.
  6. Analyze overlaps with `TCR overlap.R`.

* All scripts support both **alpha** and **beta chains**, with output stored under `outputs/{chain}/tcrdist/`.

```

I can also make a **slightly fancier version** with a **diagram of the workflow** showing how data flows through each script if you want—this is often very helpful for new users.  

Do you want me to add that diagram?
```
