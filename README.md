# WGCNA Analysis of Dengue Infection in Pediatric PBMCs (GSE38246)

Weighted Gene Co-expression Network Analysis (WGCNA) of peripheral blood
mononuclear cell (PBMC) transcriptomes from Nicaraguan children with dengue
infection, using the public dataset **GSE38246**.

## Dataset

- **GEO accession:** [GSE38246](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38246)
- **Title:** Transcriptional response to dengue infection in peripheral blood
  mononuclear cells of Nicaraguan children
- **Organism:** *Homo sapiens*
- **Platform:** GPL15615
- **Samples:** 105 dengue-infected PBMC samples (DF1 / DF2 / DHF / DSS) + 8
  healthy controls; clinical group and fever day encoded in GEO sample titles

## Objective

1. Build gene co-expression modules from PBMC expression profiles.
2. Relate modules to dengue clinical category (DF1, DF2, DHF, DSS, Healthy)
   and fever day.
3. Identify biologically meaningful modules and candidate hub genes.
4. Run GO Biological Process and KEGG functional enrichment on key modules.
5. Assess robustness of the strongest module–trait association.

## Repository structure

```
.
├── scripts/
│   └── GSE38246_WGCNA_Final_Script.R   # complete, reproducible analysis script
├── results/
│   ├── figures/                        # output plots (heatmaps, boxplots, etc.)
│   └── tables/                         # output CSVs (module sizes, correlations, hub genes...)
└── README.md
```

## How to run

1. Download the GEO series matrix file for GSE38246 and the GPL15615
   platform table (the script fetches GPL15615 automatically via
   `GEOquery::getGEO()`).
2. Set `setwd()` at the top of `scripts/GSE38246_WGCNA_Final_Script.R` to
   the folder containing `GSE38246_series_matrix.txt`.
3. Run the script top to bottom in R / RStudio. Required packages
   (WGCNA, GEOquery, clusterProfiler, org.Hs.eg.db, enrichplot, ggplot2,
   ggpubr, pheatmap) are installed automatically if missing.
4. All figures and tables are written to the working directory; move them
   into `results/figures` and `results/tables` as needed.

## Methods summary

- **QC:** Samples with >30% missing values removed (10 samples); probes
  with >20% missing values removed after sample QC.
- **Annotation:** Probes mapped to gene symbols via the GPL15615 platform
  table; multiple probes per gene collapsed to the highest-variance probe
  (21,488 unique genes).
- **Outlier removal:** Six samples with mean sample–sample correlation
  < 0.45 removed (pragmatic cutoff), leaving 97 samples.
- **Network construction:** Top 8,000 most-variable genes, KNN-imputed to
  zero missingness (97 × 7,915 final matrix). Signed network, biweight
  midcorrelation (`bicor`), soft-thresholding power = 12 (selected
  pragmatically — higher powers produced unusably sparse networks despite
  higher R²).
- **Modules:** 8 co-expression modules + grey (unassigned); analysis
  focuses on **Yellow** (ME4), **Green** (ME5), and **Pink** (ME8).
- **Trait association:** Module eigengenes correlated with 5 binary disease
  indicators + fever day (48 tests), BH-FDR corrected jointly.
- **Hub genes:** Operational criteria — `|kME| ≥ 0.80`, `|GS| ≥ 0.20`,
  genome-wide BH-FDR < 0.05 for gene significance.
- **Enrichment:** `clusterProfiler::enrichGO()` (BP) and `enrichKEGG()`,
  background = all genes in the WGCNA input matrix.

## Key results

| Module | Trait | r | FDR |
|---|---|---|---|
| Yellow (ME4) | DF1 | 0.391 | 0.0036 |
| Green (ME5) | FeverDay | −0.338 | 0.0179 |
| Green (ME5) | Healthy | −0.305 | 0.0287 |
| Pink (ME8) | FeverDay | −0.377 | 0.0060 |

- **Yellow:** platelet activation / coagulation / hemostasis (GO + KEGG
  concordant); candidate hub genes for DF1: `HIST1H2BN`, `HIST1H2AE`,
  `HIST1H2BJ`, `TNNC2`.
- **Green:** mitotic / cell-cycle processes; associated with fever day and
  the Healthy indicator.
- **Pink:** antiviral / innate-immune (interferon-stimulated gene)
  signature; associated with fever day.

## Limitations

See the **Limitations / Interpretation Caveats** section at the end of
`scripts/GSE38246_WGCNA_Final_Script.R` for a full list, including: the
pragmatic (non-R²-optimal) choice of soft-thresholding power, the pragmatic
sample-outlier cutoff, the small Healthy control group (n = 7), the
compositional nature of the disease indicator variables, and the caution
that KEGG "Hepatitis C" / "Influenza A" hits in the Pink module reflect
shared host antiviral gene sets and are not evidence of co-infection.

## License

MIT (see `LICENSE`).
