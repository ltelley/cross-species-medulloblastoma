# Downloading the large data files

Two files exceed GitHub's 100 MB limit and are **not** tracked in this repository.
Download them once into the paths below, then `Rscript run_analysis.R` will run.

| File (expected path)                 | Size   | Source |
|--------------------------------------|--------|--------|
| `data/human/GSE85217_raw.rds`        | 118 MB | GEO **GSE85217** (Cavalli et al. 2017) |
| `data/mouse_gtml/GPL6096.soft.gz`    | 576 MB | GEO platform **GPL6096** (Affymetrix Mouse Exon 1.0 ST) |

> A third file, `data/mouse/GPL1261.soft.gz` (63 MB), is **kept** in the repo
> (under the 100 MB limit), so no action is needed for it.

## Option A — regenerate with GEOquery (recommended)

```r
# install.packages("BiocManager"); BiocManager::install("GEOquery")
library(GEOquery)

# 1) Human expression series (GSE85217) -> data/human/GSE85217_raw.rds
#    The script loads this object with readRDS(); we cache the getGEO result.
dir.create("data/human", recursive = TRUE, showWarnings = FALSE)
human_list <- getGEO("GSE85217", GSEMatrix = TRUE, getGPL = FALSE)
saveRDS(human_list, "data/human/GSE85217_raw.rds")

# 2) Mouse Exon 1.0 ST platform annotation (GPL6096) -> data/mouse_gtml/GPL6096.soft.gz
dir.create("data/mouse_gtml", recursive = TRUE, showWarnings = FALSE)
# getGEO downloads and caches the SOFT file; copy it to the expected path:
gpl6096 <- getGEO("GPL6096", destdir = "data/mouse_gtml")
# getGEO saves it as data/mouse_gtml/GPL6096.soft(.gz). Ensure the name matches:
#   file.rename(list.files("data/mouse_gtml", "GPL6096.*soft.*", full.names = TRUE)[1],
#               "data/mouse_gtml/GPL6096.soft.gz")
```

## Option B — direct download from NCBI FTP

```bash
# GPL6096 platform SOFT file (~576 MB)
mkdir -p data/mouse_gtml
curl -L -o data/mouse_gtml/GPL6096.soft.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6096/soft/GPL6096_family.soft.gz"
```

GEO pages:
- GSE85217 — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85217
- GPL6096  — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6096
