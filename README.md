[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19371503.svg)](https://doi.org/10.5281/zenodo.19371503)

# Data and code associated with the paper "We might not notice a 'mass' extinction" by Giovanni Strona and Corey J.A. Bradshaw

## Repository structure

```
.
├── main.R                    # Main analysis script
├── extinction_timeline.csv   # Extinction timeline data (see Data section below)
├── COL_2026-02-13_XR/
│   └── NameUsage.tsv         # Catalogue of Life snapshot (see Data section below)
└── iucn_species_data/
    ├── <group_1>/
    │   ├── assessments.csv
    │   └── taxonomy.csv
    ├── <group_2>/
    │   ├── assessments.csv
    │   └── taxonomy.csv
    └── ...
```

The script will automatically create two output folders:
- `FIGURES/` — PDF figures
- `NUMBERS/` — numerical summaries and model outputs

## Requirements

The analysis was run in R. The following packages are required:

```r
install.packages(c("viridis", "data.table", "dplyr", "readr",
                   "ggplot2", "tidyr", "patchwork", "zoo", "mgcv", "scam", "scales"))
```

## Data

### Catalogue of Life
The file `COL_2026-02-13_XR/NameUsage.tsv` corresponds to the **Catalogue of Life snapshot of 13 February 2026**.
It can be downloaded from the [Catalogue of Life download portal](https://download.catalogueoflife.org/col/annual/).
After downloading, place `NameUsage.tsv` inside a folder named `COL_2026-02-13_XR/`.

### IUCN Red List
The `iucn_species_data/` folder contains species assessments downloaded from the [IUCN Red List](https://www.iucnredlist.org/).
Each subfolder corresponds to a taxonomic group and must contain two files: `assessments.csv` and `taxonomy.csv`.
A free account is required to access and download IUCN data.

### Extinction timeline
The file `extinction_timeline.csv` is provided in the Zenodo data repository associated with this paper (see DOI below).

## Zenodo data repository

The complete dataset (including `extinction_timeline.csv` and all input data) is archived at:

> [INSERT ZENODO DOI HERE]

## Citation

If you use this code or data, please cite:

> Strona G. and Bradshaw C.J.A. (2026). We might not notice a 'mass' extinction. *[Journal name]*.

## License

Code released under the [MIT License](LICENSE).
