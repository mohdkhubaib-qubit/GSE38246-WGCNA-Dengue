# ================================================================================================
# GSE38246 — WGCNA ANALYSIS OF DENGUE-INFECTED PBMC SAMPLES
# "Transcriptional response to dengue infection in peripheral blood mononuclear
#  cells of Nicaraguan children"
# Organism: Homo sapiens | Platform: GPL15615
#
# FINAL, COMPLETE, REPRODUCIBLE SCRIPT
#
# This script reconstructs an already-completed analysis. It does NOT redesign the
# methodology. All thresholds, sample removals, power selection, module sizes,
# correlations, p-values and FDR values below are the ACTUAL results obtained during
# the original analysis (as reproduced from the original R console history). Nothing
# in this script invents new numbers. Where an object could not be unambiguously
# reconstructed from the session history, it is flagged with "# TODO: verify this step".
#
# Working directory used in the original analysis: E:/GSE38246_WGCNA
# ================================================================================================


# ================================================================================================
# SECTION 1 — SETUP
# ================================================================================================
# Why: Load every package actually used in the analysis. Bioconductor packages are
# installed conditionally so re-running this script does not reinstall packages
# unnecessarily. WGCNA depends on 'impute' and 'preprocessCore', which are Bioconductor
# packages and are installed separately because install.packages() alone cannot resolve
# them (this caused the original "there is no package called 'impute'" / 'preprocessCore'
# errors in the console history).

setwd("E:/GSE38246_WGCNA")   # TODO: verify this step (adjust to your local path)
getwd()
options(stringsAsFactors = FALSE)

# ---- CRAN / Bioconductor manager ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# ---- GEOquery (used to download the GPL15615 platform annotation table) ----
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")

# ---- WGCNA and its Bioconductor dependencies ----
if (!requireNamespace("impute", quietly = TRUE))          BiocManager::install("impute")
if (!requireNamespace("preprocessCore", quietly = TRUE))  BiocManager::install("preprocessCore")
if (!requireNamespace("WGCNA", quietly = TRUE))            install.packages("WGCNA")

# ---- Functional enrichment ----
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))    BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("enrichplot", quietly = TRUE))      BiocManager::install("enrichplot")

# ---- Plotting / heatmaps / stats helpers ----
if (!requireNamespace("ggplot2", quietly = TRUE))  install.packages("ggplot2")
if (!requireNamespace("ggpubr", quietly = TRUE))    install.packages("ggpubr")
if (!requireNamespace("pheatmap", quietly = TRUE))  install.packages("pheatmap")

library(GEOquery)
library(preprocessCore)
library(impute)
library(WGCNA)              # masks stats::cor -- intentional; WGCNA::cor is used throughout
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(ggpubr)
library(pheatmap)

options(stringsAsFactors = FALSE)
# allowWGCNAThreads()  # TODO: verify this step (enable multithreading if desired; not used in original run)


# ================================================================================================
# SECTION 2 — DATA IMPORT
# ================================================================================================
# Why: The GEO series matrix file stores the processed expression matrix as a
# tab-delimited block beginning after the sample/series metadata header lines and
# ending at "!series_matrix_table_end". In the original file, this block began at
# line 71 (skip = 70) and contained 45,696 data rows (probes).
#
# ID_REF          = platform probe ID (integer), matches GPL15615 "ID" column
# GSM937xxx cols  = one column per sample (Cy3-channel log-ratio expression values)
# "null" strings  = missing values in the raw text file -> converted to NA via na.strings

geo_file <- "GSE38246_series_matrix.txt"
file.exists(geo_file)

lines <- readLines(geo_file)
which(grepl("series_matrix_table_begin", lines))   # 70
which(grepl("series_matrix_table_end", lines))      # 45768

expr_raw <- read.delim(
  geo_file,
  header = TRUE,
  sep = "\t",
  skip = 70,
  nrows = 45696,
  na.strings = "null",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

dim(expr_raw)                 # 45696 x 114  (ID_REF + 113 samples)
sum(is.na(expr_raw))          # 728,903 missing values in the raw matrix

# Extract sample titles (contain disease/fever-day info, e.g. "63.3 DF1 D6")
sample_line   <- lines[grep("^!Sample_title", lines)]
sample_titles <- strsplit(sample_line, "\t")[[1]]
sample_titles <- gsub('"', "", sample_titles)
sample_titles <- sample_titles[-1]
length(sample_titles)         # 113

# Interpretation:
# The raw matrix contains 45,696 probes across 113 GSM samples, with 728,903 missing
# expression values encoded as the string "null" in the original file. Quality control
# is required before any downstream network analysis.


# ================================================================================================
# SECTION 3 — SAMPLE QUALITY CONTROL (missingness-based removal)
# ================================================================================================
# Why: Samples with a very high proportion of missing expression values are unreliable
# and were removed BEFORE probe-level filtering. Removal was based purely on the
# percentage of missing values per sample (>30%), NOT on disease status -- the
# disease-state breakdown of the removed samples was inspected only afterward, for
# transparency, not as a removal criterion.

sample_missing     <- colSums(is.na(expr_raw[, -1]))
sample_missing_pct <- sample_missing / nrow(expr_raw) * 100

# Samples with >30% missing values (the 10 poor-quality samples)
bad_samples <- names(sample_missing_pct[sample_missing_pct > 30])
bad_samples
# [1] "GSM937331" "GSM937337" "GSM937350" "GSM937357" "GSM937369" "GSM937372"
# [7] "GSM937406" "GSM937411" "GSM937420" "GSM937438"

bad_indices <- match(bad_samples, colnames(expr_raw)[-1])
data.frame(
  Sample = bad_samples,
  Title  = sample_titles[bad_indices],
  Missing_Percent = sample_missing_pct[bad_samples]
)
# GSM937331  68.2 DHF D5     40.15%
# GSM937337 136.1 DHF D5     36.74%
# GSM937350  37.2 DF2 D4     49.09%
# GSM937357  55.1 DF1 D3     35.63%
# GSM937369 1922.4 Healthy   36.38%
# GSM937372  50.9 DF2 D6     34.85%
# GSM937406  29.9 DHF D7     43.23%
# GSM937411 189.2 DHF D4     31.48%
# GSM937420  37.3 DF2 D5     43.91%
# GSM937438  63.1 DF1 D4     45.93%

# Remove the 10 poor-quality samples
keep_samples <- sample_missing_pct <= 30
expr_qc <- expr_raw[, c(TRUE, keep_samples)]
dim(expr_qc)   # 45696 x 104  (ID_REF + 103 remaining samples)

# Interpretation:
# Removal of these 10 samples was based solely on missing-data burden (>30% missing),
# not on clinical group. The QC matrix retains 103 samples for probe-level filtering.


# ================================================================================================
# SECTION 4 — PROBE QUALITY CONTROL
# ================================================================================================
# Why: After sample QC, probe-level missingness was recalculated across the retained
# 103 samples. Several thresholds were evaluated before selecting <=20% missing.

probe_missing_qc  <- rowSums(is.na(expr_qc[, -1]))
probe_missing_pct <- probe_missing_qc / ncol(expr_qc[, -1]) * 100

c(
  `<=10%` = sum(probe_missing_pct <= 10),   # 25,571
  `<=20%` = sum(probe_missing_pct <= 20),   # 35,088   <- threshold used
  `<=30%` = sum(probe_missing_pct <= 30),   # 41,385
  `<=50%` = sum(probe_missing_pct <= 50)    # 45,322
)

# The <=20% missingness threshold was selected as a practical balance between
# retaining probe coverage and limiting per-probe missing-data burden.
keep_probes   <- probe_missing_pct <= 20
expr_filtered <- expr_qc[keep_probes, ]
dim(expr_filtered)      # 35,088 probes x 104 (ID_REF + 103 samples)

# Interpretation:
# expr_filtered contains 35,088 probes x 103 samples after combined sample- and
# probe-level QC, matching the pre-specified <=20% probe-missingness threshold.


# ================================================================================================
# SECTION 5 — EXPRESSION DISTRIBUTION
# ================================================================================================
# Why: Before any downstream transformation, the expression value distribution was
# inspected to confirm the data were already normalized/log-ratio transformed by the
# platform pipeline (values ranged roughly -19 to +17, centered near 0), so NO further
# log2 transformation was applied. Applying log2 to already log-ratio data would be
# an error and was deliberately avoided.

summary(as.numeric(unlist(expr_filtered[, -1])))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.       NAs
# -19.03300  -1.04600   0.15200   0.09595   1.29700  17.29200    222530

quantile(
  as.numeric(unlist(expr_filtered[, -1])),
  probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1),
  na.rm = TRUE
)
#      0%      1%      5%     25%     50%     75%     95%     99%    100%
# -19.033  -6.288  -3.595  -1.046   0.152   1.297   3.544   6.085  17.292

hist(
  as.numeric(unlist(expr_filtered[, -1])),
  breaks = 100,
  main = "Expression value distribution",
  xlab = "Expression value"
)

# Interpretation:
# The distribution is roughly symmetric around zero with a long tail, consistent with
# already log-ratio-transformed two-color array data. No further log2 transform applied.


# ================================================================================================
# SECTION 6 — PROBE ANNOTATION (mapping probes to gene symbols/ORFs via GPL15615)
# ================================================================================================
# Why: Probe-to-gene mapping is required before collapsing multiple probes per gene
# and running downstream network/enrichment analysis. Only "Experimental" probes with
# a non-empty ORF (gene symbol) field were retained as annotated experimental probes.

gpl <- getGEO("GPL15615", AnnotGPL = TRUE)   # "Annotation GPL not available, so will use submitter GPL instead"
gpl_table <- Table(gpl)
dim(gpl_table)          # 45696 x 9
table(gpl_table$`Reporter Group[Role]`)
#      Control Experimental
#         3916        41780

experimental <- gpl_table$`Reporter Group[Role]` == "Experimental"
sum(experimental & gpl_table$ORF != "")   # 40,666 experimental probes have an ORF (gene symbol)
sum(experimental & gpl_table$ORF == "")   # 1,114 experimental probes have no ORF
length(unique(gpl_table$ORF[experimental & gpl_table$ORF != ""]))  # 24,391 unique gene symbols

annot <- gpl_table[
  gpl_table$`Reporter Group[Role]` == "Experimental" &
    gpl_table$ORF != "",
  c("ID", "ORF")
]
annot$ID <- as.integer(annot$ID)
dim(annot)                                    # 40,666 x 2

sum(expr_filtered$ID_REF %in% annot$ID)       # 32,330 probes matched genes (Experimental + ORF)
sum(!expr_filtered$ID_REF %in% annot$ID)      # 2,758 probes did not match

# Breakdown of the 35,088 retained (post-QC) probes:
#   Experimental + gene symbol   32,330
#   Experimental + no ORF           466
#   Control + gene symbol         1,873
#   Control + no ORF                419
#   Total                        35,088

keep_annotated <- expr_filtered$ID_REF %in% annot$ID
expr_annotated <- expr_filtered[keep_annotated, ]

probe_annot <- annot[match(expr_annotated$ID_REF, annot$ID), ]
all(expr_annotated$ID_REF == probe_annot$ID)   # TRUE
dim(expr_annotated)     # 32,330 x 104
dim(probe_annot)        # 32,330 x 2

# Interpretation:
# 32,330 of the 35,088 QC-passing probes are annotated "Experimental" probes with a
# gene symbol and were carried forward for gene-level collapsing.


# ================================================================================================
# SECTION 7 — MULTIPLE PROBES PER GENE (collapse to one probe per gene by max SD)
# ================================================================================================
# Why: 6,776 genes were represented by more than one probe. For each such gene, the
# single probe with the HIGHEST standard deviation across the 103 samples was retained
# as the representative probe (maximizes biological signal / variance captured).

expr_values <- expr_annotated[, -1]
probe_sd <- apply(expr_values, 1, sd, na.rm = TRUE)
probe_annot$SD <- probe_sd

summary(probe_annot$SD)
sum(probe_annot$SD == 0, na.rm = TRUE)          # 0 zero-variance probes

probe_counts <- table(probe_annot$ORF)
sum(probe_counts > 1)                           # 6,776 genes with multiple probes

# Order by gene, then by decreasing SD; keep first (highest-SD) row per gene
probe_annot_ordered  <- probe_annot[order(probe_annot$ORF, -probe_annot$SD), ]
selected_probes      <- !duplicated(probe_annot_ordered$ORF)
probe_annot_selected <- probe_annot_ordered[selected_probes, ]

nrow(probe_annot_selected)                       # 21,488 genes
length(unique(probe_annot_selected$ORF))          # 21,488 (all unique)

selected_ids   <- probe_annot_selected$ID
selected_rows  <- match(selected_ids, expr_annotated$ID_REF)
expr_gene      <- expr_annotated[selected_rows, -1]
rownames(expr_gene) <- probe_annot_selected$ORF

dim(expr_gene)                # 21,488 genes x 103 samples
sum(duplicated(rownames(expr_gene)))   # 0

# Interpretation:
# The final gene-level expression matrix contains 21,488 unique genes across 103
# samples, each gene represented by its single most-variable probe.


# ================================================================================================
# SECTION 8 — SAMPLE CLUSTERING (sample-sample correlation, NOT gene-gene)
# ================================================================================================
# Why: Sample outlier detection requires a SAMPLE x SAMPLE correlation/distance matrix.
# IMPORTANT METHODOLOGICAL NOTE (preserved intentionally):
# The first attempt mistakenly computed cor(datExpr) directly on the samples-as-rows
# matrix, which WGCNA/R interprets as a GENE-GENE correlation (21,488 x 21,488) -- this
# is enormous and caused "Error: cannot allocate vector of size 1.7 Gb" during hclust().
# The CORRECT approach is to transpose appropriately so that cor() is computed on
# SAMPLES (i.e., cor(t(datExpr)) when datExpr has samples as rows / genes as columns),
# producing a 103 x 103 sample-sample correlation matrix. This distinction matters:
# gene-gene correlation describes co-expression between genes (used later in WGCNA's
# network construction itself), while sample-sample correlation describes similarity
# between arrays/patients and is used here purely for outlier detection.

datExpr <- t(expr_gene)     # samples as rows, genes as columns (103 x 21,488)

# INCORRECT (caused memory error) -- kept here only to document the mistake, NOT run:
# sample_cor  <- cor(datExpr, use = "pairwise.complete.obs", method = "pearson")  # gene-gene!
# sample_dist <- as.dist(1 - sample_cor)
# sample_tree <- hclust(sample_dist, method = "average")   # Error: cannot allocate vector of size 1.7 Gb

# CORRECT: sample-sample correlation
sample_cor <- cor(
  t(datExpr),
  use = "pairwise.complete.obs",
  method = "pearson"
)
dim(sample_cor)             # 103 x 103

sample_dist <- as.dist(1 - sample_cor)
sample_tree <- hclust(sample_dist, method = "average")

plot(
  sample_tree,
  main = "Sample clustering",
  xlab = "", sub = "", cex = 0.6
)

# Interpretation:
# Gene-gene correlation (cor(datExpr) on a samples-as-rows matrix) would compute
# correlations BETWEEN GENES and is inappropriate/impractical for sample outlier
# detection at this scale. Sample-sample correlation (cor(t(datExpr))) is the correct
# and computationally tractable approach used here.


# ================================================================================================
# SECTION 9 — SAMPLE OUTLIER REMOVAL (mean sample-sample correlation)
# ================================================================================================
# Why: Samples with low average correlation to all other samples are candidate
# technical/biological outliers. A pragmatic cutoff of mean correlation < 0.45 was
# used to flag six clear outliers. This cutoff is NOT presented as an objectively
# optimal threshold -- it was a pragmatic decision, and a sensitivity analysis of this
# choice would ideally be performed (see Section 29 / Limitations).

mean_sample_cor <- rowMeans(sample_cor, na.rm = TRUE)

outlier_candidates <- data.frame(
  Sample = names(mean_sample_cor),
  MeanCorrelation = mean_sample_cor
)
outlier_candidates <- outlier_candidates[order(outlier_candidates$MeanCorrelation), ]
head(outlier_candidates, 15)

# Match sample titles for context (NOT used as a removal criterion)
sample_ids <- colnames(expr_gene)
sample_meta <- data.frame(
  Sample = sample_ids,
  Title  = sample_titles[match(sample_ids, colnames(expr_raw)[-1])],
  stringsAsFactors = FALSE
)
outlier_candidates$Title <- sample_meta$Title[match(outlier_candidates$Sample, sample_meta$Sample)]
head(outlier_candidates, 15)

# The six samples with mean correlation < 0.45 (removed):
remove_samples <- c(
  "GSM937397",  # 0.352, 191.1 DHF D4
  "GSM937359",  # 0.382,  31.9 DHF D8
  "GSM937366",  # 0.386,  50.2 DF2 D5
  "GSM937396",  # 0.434,  44.1 DSS D4
  "GSM937426",  # 0.437,  21.3 DHF D6
  "GSM937424"   # 0.442, 174.2 DHF D5
)

# Interpretation:
# This cutoff (mean correlation < 0.45) was a pragmatic decision, not a formal
# statistical threshold. Removal of these six samples reduces the working sample set
# to 97, used for all downstream WGCNA network construction.


# ================================================================================================
# SECTION 10 — WGCNA PREPROCESSING (gene filtering, missingness, KNN imputation)
# ================================================================================================
# Why: WGCNA requires a complete (no-NA) expression matrix. The workflow here:
#   (1) goodSamplesGenes() standard WGCNA QC check
#   (2) restrict to the 8,000 most variable genes (pragmatic dimensionality reduction)
#   (3) remove the 6 low-correlation outlier samples identified in Section 9
#   (4) remove any remaining genes with >20% missingness
#   (5) KNN-impute the small remaining fraction of missing values
# KNN imputation is performed AFTER variance-based filtering and outlier removal so
# that imputation is not wasted computing values for genes/samples later discarded,
# and so that the neighbor structure used by KNN reflects the final, cleaned sample set.

gsg <- goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK    # TRUE

# Select the 8,000 most variable genes (pragmatic feature-selection choice)
gene_variance <- apply(datExpr, 2, var, na.rm = TRUE)
n_genes <- 8000
top_genes <- names(sort(gene_variance, decreasing = TRUE))[1:n_genes]
datExpr_wgcna <- datExpr[, top_genes]
dim(datExpr_wgcna)     # 103 x 8000

# Remove the six low-correlation outlier samples
datExpr_wgcna_clean <- datExpr_wgcna[!rownames(datExpr_wgcna) %in% remove_samples, ]
dim(datExpr_wgcna_clean)   # 97 x 8000

# Remove genes with >20% missingness after sample QC
missing_gene_clean     <- colSums(is.na(datExpr_wgcna_clean))
missing_gene_clean_pct <- missing_gene_clean / nrow(datExpr_wgcna_clean) * 100
sum(missing_gene_clean_pct > 20)   # 85 genes exceed 20% missingness

keep_genes_clean   <- missing_gene_clean_pct <= 20
datExpr_wgcna_final <- datExpr_wgcna_clean[, keep_genes_clean]
dim(datExpr_wgcna_final)   # 97 x 7,915

# KNN imputation (performed on the FINAL, filtered matrix)
imputed <- impute.knn(t(datExpr_wgcna_final))
datExpr_wgcna_imputed <- t(imputed$data)
dim(datExpr_wgcna_imputed)     # 97 x 7,915
sum(is.na(datExpr_wgcna_imputed))   # 0 -- no missing values remain

# Interpretation:
# The final pre-network input matrix contains 97 samples x 7,915 genes with zero
# missing values. Both the 8,000-gene variance filter and the 97-sample set reflect
# pragmatic, explicitly-flagged analytical choices rather than data-driven optima.


# ================================================================================================
# SECTION 11 — SOFT THRESHOLD ANALYSIS
# ================================================================================================
# Why: pickSoftThreshold() was run using a signed network with biweight midcorrelation
# (bicor), across a wide range of candidate powers, to evaluate the scale-free
# topology fit (R^2) at each power. A signed network + bicor combination was used
# throughout for robustness to outliers and to preserve the direction of correlation
# (up- vs down-regulation) in the resulting network.

powers <- c(6:20, seq(22, 30, by = 2))

sft_bicor_signed <- pickSoftThreshold(
  datExpr_wgcna_imputed,
  powerVector = powers,
  networkType = "signed",
  corFnc = "bicor",
  verbose = 5
)

data.frame(
  Power = sft_bicor_signed$fitIndices$Power,
  R2 = round(sft_bicor_signed$fitIndices$SFT.R.sq, 3),
  MeanConnectivity = round(sft_bicor_signed$fitIndices$mean.k., 1)
)
# Selected values from the actual run:
#  Power 12 -> R2 = 0.383
#  Power 14 -> R2 = 0.505
#  Power 18 -> R2 = 0.659
#  Power 20 -> R2 = 0.709
#  Power 26 -> R2 = 0.814
#  Power 30 -> R2 = 0.847

# --- Diagnostic network-construction trials at high R^2 powers (documented, NOT final) ---
# Power 26: blockwiseModules() detected NO modules at all in either block
#           (mergeCloseModules failed: "Color levels are empty ... only grey").
# Power 18: table(net$colors) -> 0:7752, 1:125, 2:38  -> only 163/7915 genes (~2%) non-grey,
#           a very small number of non-grey genes.
# Power 14: table(net$colors) -> 0:7444,1:306,2:112,3:53 -> 471/7915 = 5.95% non-grey genes.
# Power 12: table(net$colors) -> 0:5747,1:1516,2:431,3:84,4:81,5:56 -> 2168/7915 = 27.39%
#           non-grey genes assigned across 5 (pre-merge) modules -- a practically usable network.

# IMPORTANT: Power 12 was NOT selected because it was the mathematically "optimal" power
# by R^2 alone (R^2 = 0.383 at power 12 is well below the classic >=0.80 rule-of-thumb).
# It was selected PRAGMATICALLY: powers achieving higher R^2 (18, 26, etc.) produced
# excessively sparse networks with almost no genes assigned to any non-grey module, or
# no modules at all, and were therefore not usable for downstream module-trait analysis.

soft_power <- 12


# ================================================================================================
# SECTION 12 — FINAL WGCNA NETWORK
# ================================================================================================
# Why: The final network was constructed with the exact parameters below. corType =
# "bicor" (robust correlation) and networkType = "signed" / TOMType = "signed" were
# deliberately used (NOT pearson / unsigned) -- these were corrections made during the
# original analysis and must not be silently reverted.

net_final <- blockwiseModules(
  datExpr_wgcna_imputed,
  power = 12,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  maxPOutliers = 1,
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  pamRespectsDendro = FALSE,
  numericLabels = TRUE,
  verbose = 3
)

table(labels2colors(net_final$colors))
#    black      blue     brown     green      grey      pink       red turquoise    yellow
#       43      1303       234        89       146        36        63      5874       127

length(unique(labels2colors(net_final$colors)))   # 9 (8 colored modules + grey)

# Visualize module structure (both pre-clustering blocks)
plotDendroAndColors(
  net_final$dendrograms[[1]],
  labels2colors(net_final$colors)[net_final$blockGenes[[1]]],
  "Module",
  main = "Gene dendrogram and WGCNA modules - Block 1",
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05
)
plotDendroAndColors(
  net_final$dendrograms[[2]],
  labels2colors(net_final$colors)[net_final$blockGenes[[2]]],
  "Module",
  main = "Gene dendrogram and WGCNA modules - Block 2",
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05
)

# Interpretation:
# The final signed, bicor-based network at power 12 assigns 2,168+ genes to 8
# non-grey modules; the 5,874-gene turquoise module and 146-gene grey (unassigned)
# module are the largest groups, with several smaller, more biologically tractable
# modules (yellow, green, pink, etc.) used for downstream trait association.


# ================================================================================================
# SECTION 13 — FINAL MODULES (sizes and module-number-to-color mapping)
# ================================================================================================
# Why: blockwiseModules() returns numeric module labels (net_final$colors); these were
# converted to WGCNA standard color names via labels2colors(), and mapped explicitly
# to the module-eigengene column names (MEx) used throughout the rest of the script.

module_colors <- labels2colors(net_final$colors)
module_sizes <- sort(table(module_colors), decreasing = TRUE)
module_sizes
# turquoise      blue     brown      grey    yellow     green       red     black      pink
#      5874      1303       234       146       127        89        63        43        36

module_color_map <- data.frame(
  Module = names(table(net_final$colors)),
  Color  = labels2colors(as.numeric(names(table(net_final$colors)))),
  Genes  = as.numeric(table(net_final$colors))
)
module_color_map
#   0      grey   146
#   1 turquoise  5874
#   2      blue  1303
#   3     brown   234
#   4    yellow   127   <- ME4 = Yellow
#   5     green    89   <- ME5 = Green
#   6       red    63
#   7     black    43
#   8      pink    36   <- ME8 = Pink

write.csv(as.data.frame(module_sizes), "WGCNA_Module_Sizes.csv", row.names = FALSE)

# Interpretation:
# There are 8 biological (colored) modules plus the grey (unassigned) "module".
# ME4 = Yellow, ME5 = Green, ME8 = Pink -- these are the three modules on which the
# biological interpretation of this project centers (Sections 18-29).


# ================================================================================================
# SECTION 14 — MODULE EIGENGENES
# ================================================================================================
# Why: A module eigengene (ME) is the first principal component of a module's
# standardized expression profile across samples -- it summarizes the dominant
# co-expression pattern of that module in a single value per sample, and is the
# quantity correlated against clinical traits below.

MEs_final <- orderMEs(net_final$MEs)
MEs_nonGrey <- MEs_final[, !grepl("^ME0$", names(MEs_final))]

dim(MEs_final)     # 97 x 9
dim(MEs_nonGrey)    # 97 x 8

# Interpretation:
# MEs_nonGrey excludes the grey (unassigned-genes) pseudo-module, since grey does not
# represent a biologically coherent co-expression module.


# ================================================================================================
# SECTION 15 — CLINICAL METADATA
# ================================================================================================
# Why: Clinical/disease information is embedded in the GEO sample title field (e.g.
# "63.3 DF1 D6" = patient 63.3, primary dengue fever, fever day 6; "932.4 Healthy" =
# healthy control). Disease category and fever day were extracted directly from the
# title strings -- NEVER inferred from the GSM accession number itself.

traitData <- sample_meta[match(rownames(MEs_final), sample_meta$Sample), ]

traitData$Disease <- ifelse(
  grepl("Healthy", traitData$Title),
  "Healthy",
  sub("^[^ ]+ ([A-Za-z0-9]+) D[0-9]+$", "\\1", traitData$Title)
)

traitData$FeverDay <- NA_real_
dengue_idx <- traitData$Disease != "Healthy"
traitData$FeverDay[dengue_idx] <- as.numeric(
  sub(".* D([0-9]+)$", "\\1", traitData$Title[dengue_idx])
)

table(traitData$Disease, useNA = "ifany")
#     DF1     DF2     DHF     DSS Healthy
#      21      26      24      19       7

table(traitData$FeverDay, useNA = "ifany")   # 90 dengue samples have a FeverDay value; 7 Healthy are NA

# Binary disease-group traits
traitData$DF1     <- as.numeric(traitData$Disease == "DF1")
traitData$DF2     <- as.numeric(traitData$Disease == "DF2")
traitData$DHF     <- as.numeric(traitData$Disease == "DHF")
traitData$DSS     <- as.numeric(traitData$Disease == "DSS")
traitData$Healthy <- as.numeric(traitData$Disease == "Healthy")

all(rownames(MEs_final) == traitData$Sample)   # TRUE

# Interpretation:
# Final sample counts (n = 97) after outlier removal: DF1 = 21, DF2 = 26, DHF = 24,
# DSS = 19, Healthy = 7. FeverDay is available for the 90 dengue (non-Healthy) samples.
# IMPORTANT: the Healthy group (n = 7) is small, limiting statistical power for any
# comparison involving Healthy controls -- interpret Healthy-related results cautiously.
# The five disease binary variables (DF1/DF2/DHF/DSS/Healthy) are mutually exclusive
# and compositional (they sum to 1 per sample), which should be kept in mind when
# interpreting their correlations with module eigengenes.


# ================================================================================================
# SECTION 16 — MODULE-TRAIT CORRELATIONS
# ================================================================================================
# Why: Pearson correlation between each module eigengene and each clinical trait
# (5 binary disease indicators + continuous FeverDay), with Student's t-based p-values
# via corPvalueStudent() (using the actual number of complete-case samples per trait --
# 97 for the disease indicators, 90 for FeverDay), followed by Benjamini-Hochberg FDR
# correction applied ACROSS ALL 48 module-trait tests jointly (8 modules x 6 traits).

traitData_final <- traitData[match(rownames(MEs_nonGrey), traitData$Sample), ]
all(rownames(MEs_nonGrey) == traitData_final$Sample)  # TRUE

trait_matrix <- traitData_final[, c("DF1", "DF2", "DHF", "DSS", "Healthy", "FeverDay")]

moduleTraitCor_final <- cor(
  MEs_nonGrey, trait_matrix,
  use = "pairwise.complete.obs", method = "pearson"
)

moduleTraitPvalue_final <- matrix(
  NA, nrow = nrow(moduleTraitCor_final), ncol = ncol(moduleTraitCor_final),
  dimnames = dimnames(moduleTraitCor_final)
)
for (j in 1:ncol(trait_matrix)) {
  complete_n <- sum(complete.cases(MEs_nonGrey, trait_matrix[, j]))
  moduleTraitPvalue_final[, j] <- corPvalueStudent(moduleTraitCor_final[, j], nSamples = complete_n)
  cat(colnames(trait_matrix)[j], ":", complete_n, "complete samples\n")
}
# DF1 : 97 | DF2 : 97 | DHF : 97 | DSS : 97 | Healthy : 97 | FeverDay : 90

moduleTraitFDR_final <- matrix(
  p.adjust(as.vector(moduleTraitPvalue_final), method = "BH"),
  nrow = nrow(moduleTraitPvalue_final), ncol = ncol(moduleTraitPvalue_final),
  dimnames = dimnames(moduleTraitPvalue_final)
)

round(moduleTraitCor_final, 3)
#        DF1    DF2    DHF    DSS Healthy FeverDay
# ME4  0.391  0.069 -0.266 -0.094  -0.154    0.021
# ME6  0.009  0.059  0.041 -0.044  -0.115    0.083
# ME1  0.163 -0.032 -0.170  0.063  -0.017   -0.011
# ME2  0.017  0.086  0.016 -0.069  -0.094    0.176
# ME5 -0.242  0.230  0.055  0.134  -0.305   -0.338
# ME8  0.076  0.065 -0.042  0.076  -0.277   -0.377
# ME3  0.181 -0.143 -0.177  0.083   0.125   -0.209
# ME7  0.035 -0.017 -0.069  0.077  -0.030   -0.146

round(moduleTraitFDR_final, 4)
#        DF1    DF2    DHF    DSS Healthy FeverDay
# ME4 0.0036 0.7800 0.0683 0.7800  0.4260   0.9171
# ME6 0.9333 0.8018 0.8701 0.8701  0.6270   0.7800
# ME1 0.3817 0.9015 0.3610 0.7846  0.9171   0.9333
# ME2 0.9171 0.7800 0.9171 0.7800  0.7800   0.3610
# ME5 0.1156 0.1405 0.8172 0.5073  0.0287   0.0179
# ME8 0.7800 0.7846 0.8701 0.7800  0.0573   0.0060
# ME3 0.3610 0.4793 0.3610 0.7800  0.5575   0.2580
# ME7 0.9015 0.9171 0.7800 0.7800  0.9015   0.4793

write.csv(moduleTraitCor_final, "Module_Trait_Correlations.csv")
write.csv(moduleTraitFDR_final, "Module_Trait_FDR.csv")

# Interpretation:
# Significant results (FDR < 0.05, out of 48 module-trait tests):
#   ME4 (Yellow) / DF1      : r =  0.391, FDR = 0.0036
#   ME5 (Green)  / Healthy  : r = -0.305, FDR = 0.0287
#   ME5 (Green)  / FeverDay : r = -0.338, FDR = 0.0179
#   ME8 (Pink)   / FeverDay : r = -0.377, FDR = 0.0060
# Near-significant:
#   ME8 (Pink)   / Healthy  : FDR = 0.0573
#   ME4 (Yellow) / DHF      : FDR = 0.0683
# All other module-trait associations do not survive FDR correction.


# ================================================================================================
# SECTION 17 — MODULE-TRAIT HEATMAP
# ================================================================================================
# Why: Standard WGCNA labeledHeatmap() visualization of the module-trait correlation
# matrix, annotated with correlation and FDR in each cell, using a blue-white-red
# diverging palette (blueWhiteRed(50)) and a fixed color scale (zlim = c(-1, 1)) so
# color intensity is comparable across all modules/traits.

textMatrix <- paste(
  signif(moduleTraitCor_final, 2), "\n(", signif(moduleTraitFDR_final, 2), ")", sep = ""
)

png("FINAL_Module_Trait_Heatmap.png", width = 2200, height = 1800, res = 300)
par(mar = c(7, 7, 3, 2))
labeledHeatmap(
  Matrix = moduleTraitCor_final,
  xLabels = colnames(trait_matrix),
  yLabels = rownames(moduleTraitCor_final),
  ySymbols = rownames(moduleTraitCor_final),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.65,
  cex.lab = 1.0,
  zlim = c(-1, 1),
  main = "Module\u2013trait relationships"
)
dev.off()


# ================================================================================================
# SECTION 18 — CLINICAL GROUP ANALYSIS (Kruskal-Wallis + pairwise Wilcoxon)
# ================================================================================================
# Why: For the three biologically-focal modules (Yellow, Green, Pink), non-parametric
# Kruskal-Wallis tests across the five clinical groups were followed by pairwise
# Wilcoxon rank-sum tests with Benjamini-Hochberg correction, to characterize which
# specific group comparisons drive the overall group difference.

# ---- YELLOW (ME4) ----
yellow_group_df <- data.frame(
  Sample = rownames(MEs_nonGrey),
  Disease = traitData_final$Disease,
  Yellow_ME = as.numeric(MEs_nonGrey[, "ME4"])
)
yellow_group_df$Disease <- factor(yellow_group_df$Disease, levels = c("Healthy","DF1","DF2","DHF","DSS"))

ggplot(yellow_group_df, aes(x = Disease, y = Yellow_ME)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.12, size = 2.5) +
  labs(title = "Yellow module eigengene across clinical groups",
       x = "Clinical group", y = "Yellow module eigengene") +
  theme_classic(base_size = 14)
ggsave("Yellow_Module_Clinical_Groups.png", width = 8, height = 6, dpi = 300)

yellow_kw <- kruskal.test(Yellow_ME ~ Disease, data = yellow_group_df)
yellow_kw
# Kruskal-Wallis chi-squared = 18.052, df = 4, p-value = 0.001206

yellow_pairwise <- pairwise.wilcox.test(
  yellow_group_df$Yellow_ME, yellow_group_df$Disease, p.adjust.method = "BH", exact = FALSE
)
yellow_pairwise
#             Healthy    DF1     DF2     DHF
# DF1          0.0318      -       -       -
# DF2          0.1490   0.2014     -       -
# DHF          0.9812   0.0014  0.0318     -
# DSS          0.6258   0.0318  0.3786  0.5514

# ---- GREEN (ME5) ----
green_group_df <- data.frame(
  Sample = rownames(MEs_nonGrey),
  Disease = traitData_final$Disease,
  Green_ME = as.numeric(MEs_nonGrey[, "ME5"])
)
green_group_df$Disease <- factor(green_group_df$Disease, levels = c("Healthy","DF1","DF2","DHF","DSS"))

ggplot(green_group_df, aes(x = Disease, y = Green_ME)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.12, size = 2.2) +
  labs(title = "Green module eigengene across clinical groups",
       x = "Clinical group", y = "Green module eigengene") +
  theme_classic(base_size = 14)
ggsave("Green_Module_Clinical_Groups.png", width = 8, height = 6, dpi = 300)

green_kw <- kruskal.test(Green_ME ~ Disease, data = green_group_df)
green_kw
# Kruskal-Wallis chi-squared = 19.21, df = 4, p-value = 0.0007148

green_pairwise <- pairwise.wilcox.test(
  green_group_df$Green_ME, green_group_df$Disease, p.adjust.method = "BH", exact = FALSE
)
green_pairwise
#             Healthy    DF1     DF2     DHF
# DF1          0.083       -       -       -
# DF2          0.007    0.007      -       -
# DHF          0.016    0.091   0.421      -
# DSS          0.013    0.061   0.638   0.592

# ---- PINK (ME8) ----
pink_group_df <- data.frame(
  Sample = rownames(MEs_nonGrey),
  Disease = traitData_final$Disease,
  Pink_ME = as.numeric(MEs_nonGrey[, "ME8"])
)
pink_group_df$Disease <- factor(pink_group_df$Disease, levels = c("Healthy","DF1","DF2","DHF","DSS"))

ggplot(pink_group_df, aes(x = Disease, y = Pink_ME)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.12, size = 2.2) +
  labs(title = "Pink module eigengene across clinical groups",
       x = "Clinical group", y = "Pink module eigengene") +
  theme_classic(base_size = 14)

pink_kw <- kruskal.test(Pink_ME ~ Disease, data = pink_group_df)
pink_kw
# Kruskal-Wallis chi-squared = 10.031, df = 4, p-value = 0.03991

pink_pairwise <- pairwise.wilcox.test(
  pink_group_df$Pink_ME, pink_group_df$Disease, p.adjust.method = "BH", exact = FALSE
)
pink_pairwise
#             Healthy    DF1     DF2     DHF
# DF1          0.032       -       -       -
# DF2          0.032    0.957      -       -
# DHF          0.018    0.743   0.829      -
# DSS          0.018    0.906   0.957   0.743

# Interpretation:
# Yellow differs significantly overall (p = 0.0012), driven mainly by DF1 vs Healthy/
# DHF/DSS contrasts. Green differs significantly overall (p = 0.0007), driven mainly by
# Healthy vs DF2/DHF/DSS. Pink differs overall (p = 0.040), driven mainly by Healthy vs
# every dengue subgroup. Interpret cautiously given the small Healthy group (n = 7).


# ================================================================================================
# SECTION 19 — HUB GENE ANALYSIS (module membership kME + gene significance GS)
# ================================================================================================
# Why: Module membership (kME, via signedKME) quantifies how strongly each gene's
# expression correlates with its own module eigengene. Gene significance (GS)
# quantifies each gene's correlation with a specific clinical trait. Candidate hub
# genes require BOTH high kME AND high |GS| AND a significant FDR-corrected GS p-value.
# CRITICAL: BH-FDR correction for GS p-values is applied ACROSS ALL 7,915 genes in the
# network for each trait (not only within a candidate shortlist), avoiding the
# statistical error of correcting only among pre-selected candidates.

expr_final  <- datExpr_wgcna_imputed
gene_colors <- labels2colors(net_final$colors)
gene_names  <- colnames(expr_final)

kME <- signedKME(expr_final, MEs_nonGrey, outputColumnName = "kME")
# colnames(kME): kME4 kME6 kME1 kME2 kME5 kME8 kME3 kME7

GS_DF1     <- cor(expr_final, traitData_final$DF1,     use = "pairwise.complete.obs")
GS_FeverDay<- cor(expr_final, traitData_final$FeverDay,use = "pairwise.complete.obs")
GS_Healthy <- cor(expr_final, traitData_final$Healthy, use = "pairwise.complete.obs")

GS_DF1_p      <- corPvalueStudent(GS_DF1,      nSamples = 97)
GS_FeverDay_p <- corPvalueStudent(GS_FeverDay, nSamples = 90)
GS_Healthy_p  <- corPvalueStudent(GS_Healthy,  nSamples = 97)

# Global BH-FDR correction ACROSS ALL 7,915 genes for each trait
GS_DF1_FDR      <- p.adjust(GS_DF1_p,      method = "BH")
GS_FeverDay_FDR <- p.adjust(GS_FeverDay_p, method = "BH")
GS_Healthy_FDR  <- p.adjust(GS_Healthy_p,  method = "BH")

# Strict candidate-hub-gene criteria: |kME| >= 0.80, |GS| >= 0.20, GS_FDR < 0.05
yellow_idx <- gene_colors == "yellow"
green_idx  <- gene_colors == "green"
pink_idx   <- gene_colors == "pink"

yellow_hubs_final <- data.frame(
  Gene = gene_names[yellow_idx], kME = kME[yellow_idx, "kME4"],
  GS = GS_DF1[yellow_idx], GS_p = GS_DF1_p[yellow_idx], GS_FDR = GS_DF1_FDR[yellow_idx]
)
yellow_hubs_final <- subset(yellow_hubs_final, abs(kME) >= 0.80 & abs(GS) >= 0.20 & GS_FDR < 0.05)
yellow_hubs_final <- yellow_hubs_final[order(-abs(yellow_hubs_final$kME)), ]

green_hubs_fever_final <- data.frame(
  Gene = gene_names[green_idx], kME = kME[green_idx, "kME5"],
  GS = GS_FeverDay[green_idx], GS_p = GS_FeverDay_p[green_idx], GS_FDR = GS_FeverDay_FDR[green_idx]
)
green_hubs_fever_final <- subset(green_hubs_fever_final, abs(kME) >= 0.80 & abs(GS) >= 0.20 & GS_FDR < 0.05)

pink_hubs_fever_final <- data.frame(
  Gene = gene_names[pink_idx], kME = kME[pink_idx, "kME8"],
  GS = GS_FeverDay[pink_idx], GS_p = GS_FeverDay_p[pink_idx], GS_FDR = GS_FeverDay_FDR[pink_idx]
)
pink_hubs_fever_final <- subset(pink_hubs_fever_final, abs(kME) >= 0.80 & abs(GS) >= 0.20 & GS_FDR < 0.05)

nrow(yellow_hubs_final)          # 4
nrow(green_hubs_fever_final)     # 0
nrow(pink_hubs_fever_final)      # 0

yellow_hubs_final
#          Gene       kME        GS      GS_FDR
#      HIST1H2BN 0.8187457 0.4210180 0.01387969
#      HIST1H2AE 0.8223197 0.4114770 0.01863722
#      HIST1H2BJ 0.8251694 0.4026897 0.02379173
#          TNNC2 0.8742395 0.3850752 0.04323535

# Consolidated hub table
hub_results <- rbind(
  data.frame(yellow_hubs_final[, c("Gene","kME","GS","GS_FDR")], Module = "Yellow"),
  data.frame(green_hubs_fever_final[, c("Gene","kME","GS","GS_FDR")], Module = "Green")[0, ],
  data.frame(pink_hubs_fever_final[, c("Gene","kME","GS","GS_FDR")], Module = "Pink")[0, ]
)
hub_results <- hub_results[order(hub_results$Module, hub_results$GS_FDR), ]
write.csv(hub_results, "Hub_Genes.csv", row.names = FALSE)

# Interpretation:
# Under the strict operational hub criteria, four Yellow/DF1 candidate hub genes were
# identified: HIST1H2BN, HIST1H2AE, HIST1H2BJ, TNNC2 (kME 0.82-0.87, GS 0.39-0.42,
# FDR 0.014-0.043). These are "candidate hub genes under the defined criteria" and NOT
# universally-defined WGCNA hub genes. No genes passed all strict criteria for
# Green/FeverDay or Pink/FeverDay, even though the modules themselves are significantly
# associated with FeverDay at the eigengene level (Section 16) -- module-level
# significance does not guarantee individual-gene-level significance after genome-wide
# FDR correction.


# ================================================================================================
# SECTION 20 — PINK ANTIVIRAL/INNATE-IMMUNE GENES
# ================================================================================================
# Why: The Pink module was enriched for antiviral/innate-immune GO terms (Section 22).
# A curated set of antiviral/interferon-stimulated genes within Pink was examined for
# module membership (kME8) and gene significance for FeverDay (with global FDR).

antiviral_genes <- c("ETV7","XAF1","IFIT1","RSAD2","LAMP3","IFIT2","SIGLEC1","CXCL10","CCL2","OAS1")

pink_genes_df <- data.frame(
  Gene = gene_names[pink_idx], kME = kME[pink_idx, "kME8"],
  GS_FeverDay = GS_FeverDay[pink_idx], GS_FDR = GS_FeverDay_FDR[pink_idx]
)
pink_genes_df$Priority <- abs(pink_genes_df$kME) * abs(pink_genes_df$GS_FeverDay)
pink_genes_ranked <- pink_genes_df[order(-pink_genes_df$Priority), ]

pink_antiviral_final <- pink_genes_ranked[
  pink_genes_ranked$Gene %in% antiviral_genes,
  c("Gene", "kME", "GS_FeverDay", "GS_FDR", "Priority")
]
pink_antiviral_final
#     Gene       kME GS_FeverDay      GS_FDR
#     ETV7 0.7960955  -0.4140780  0.04370434
#    IFIT1 0.7807453  -0.3473550  0.18776637
#     XAF1 0.6115395  -0.4273184  0.04370434
#    RSAD2 0.7716852  -0.3197458  0.29701576
#    LAMP3 0.8057392  -0.2957388  0.37761288
#    IFIT2 0.8148799  -0.2584996  0.64865130
#  SIGLEC1 0.6507684  -0.3138758  0.31064169
#   CXCL10 0.8003792  -0.2281400  0.77368855
#     CCL2 0.8058115  -0.2227972  0.77807453
#     OAS1 0.5213891  -0.2824824  0.46069218

write.csv(pink_antiviral_final, "Pink_Antiviral_Genes.csv", row.names = FALSE)

# Antiviral gene expression heatmap ordered by fever day
dengue_samples <- traitData_final$Sample[!is.na(traitData_final$FeverDay)]
pink_heatmap_df <- expr_final[dengue_samples, antiviral_genes, drop = FALSE]
fever_order <- order(traitData_final$FeverDay[match(rownames(pink_heatmap_df), traitData_final$Sample)])
pink_heatmap_df <- pink_heatmap_df[fever_order, , drop = FALSE]
pink_heatmap_scaled <- scale(pink_heatmap_df)

png("Pink_Antiviral_Gene_Heatmap.png", width = 1800, height = 900, res = 200)
pheatmap(
  t(pink_heatmap_scaled), cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = FALSE,
  fontsize_row = 11, main = "Antiviral and innate-immune genes in the Pink module"
)
dev.off()

# Interpretation:
# ETV7 and XAF1 both reach genome-wide GS-FDR < 0.05 for FeverDay (0.0437 each), but
# neither satisfies the strict kME >= 0.80 criterion (ETV7 = 0.796, XAF1 = 0.612), so
# neither qualifies as a strict candidate hub gene under the Section 19 criteria even
# though their gene-level trait association is itself statistically significant.


# ================================================================================================
# SECTION 21 — FUNCTIONAL ENRICHMENT (GO Biological Process)
# ================================================================================================
# Why: Module gene symbols were mapped to Entrez IDs via bitr() and tested for GO
# Biological Process enrichment with enrichGO(), using the full 7,915-gene WGCNA input
# set (mapped to 5,467 unique Entrez IDs) as the statistical background/universe.

yellow_genes <- gene_names[gene_colors == "yellow"]
green_genes  <- gene_names[gene_colors == "green"]
pink_genes   <- gene_names[gene_colors == "pink"]
length(yellow_genes)   # 127
length(green_genes)    # 89
length(pink_genes)     # 36

yellow_map <- bitr(yellow_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
green_map  <- bitr(green_genes,  fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
pink_map   <- bitr(pink_genes,   fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

yellow_entrez <- unique(yellow_map$ENTREZID)   # 104 unique Entrez IDs / 127 genes
green_entrez  <- unique(green_map$ENTREZID)    #  67 unique Entrez IDs /  89 genes
pink_entrez   <- unique(pink_map$ENTREZID)     #  34 unique Entrez IDs /  36 genes
# Mapping failures: ~18.9% of Yellow genes, ~24.7% of Green genes, ~5.6% of Pink genes
# failed to map to an Entrez ID (mostly histone-cluster variants, BCR/TCR segments,
# LOC placeholders, and non-coding RNA probe annotations).

background_map    <- bitr(gene_names, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
background_entrez <- unique(background_map$ENTREZID)
length(background_entrez)   # 5,467

yellow_GO_BP <- enrichGO(
  gene = yellow_entrez, universe = background_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE
)
green_GO_BP <- enrichGO(
  gene = green_entrez, universe = background_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE
)
pink_GO_BP <- enrichGO(
  gene = pink_entrez, universe = background_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE
)

nrow(as.data.frame(yellow_GO_BP))   # 41 significant BP terms
nrow(as.data.frame(green_GO_BP))    # 80
nrow(as.data.frame(pink_GO_BP))     # 31

# Redundant/highly-overlapping GO terms were collapsed via simplify()
yellow_GO_simplified <- simplify(yellow_GO_BP, cutoff = 0.7, by = "p.adjust", select_fun = min)
green_GO_simplified  <- simplify(green_GO_BP,  cutoff = 0.7, by = "p.adjust", select_fun = min)
pink_GO_simplified   <- simplify(pink_GO_BP,   cutoff = 0.7, by = "p.adjust", select_fun = min)

dotplot(yellow_GO_simplified, showCategory = 10, title = "Yellow Module \u2013 GO BP")
dotplot(green_GO_simplified,  showCategory = 10, title = "Green Module \u2013 GO BP")
dotplot(pink_GO_simplified,   showCategory = 10, title = "Pink Module \u2013 GO BP")

write.csv(as.data.frame(yellow_GO_BP), "Yellow_GO_BP.csv", row.names = FALSE)
write.csv(as.data.frame(green_GO_BP),  "Green_GO_BP.csv",  row.names = FALSE)
write.csv(as.data.frame(pink_GO_BP),   "Pink_GO_BP.csv",   row.names = FALSE)


# ================================================================================================
# SECTION 22 — GO INTERPRETATION
# ================================================================================================
# YELLOW: strong enrichment for platelet activation, wound healing, blood coagulation,
#   hemostasis, response to wounding, regulation of body fluid levels, platelet
#   aggregation, myeloid cell differentiation.
#
# GREEN: strong enrichment for chromosome localization, mitotic sister chromatid
#   segregation, mitotic nuclear division, metaphase chromosome alignment, nuclear
#   chromosome segregation, sister chromatid segregation, mitotic cell cycle,
#   chromosome segregation.
#
# PINK: enrichment for viral life cycle, viral genome replication, viral process,
#   antiviral innate immune response, defense response to other organism, response to
#   cytokine, response to virus.
#
# CAUTION: GO "viral process"/"viral life cycle" enrichment in the Pink module reflects
# HOST gene signatures/pathways that respond to or interact with viral infection -- it
# does NOT by itself constitute direct evidence of active viral replication in these
# samples.


# ================================================================================================
# SECTION 23 — KEGG PATHWAY ENRICHMENT
# ================================================================================================
yellow_KEGG <- enrichKEGG(
  gene = yellow_entrez, universe = background_entrez, organism = "hsa",
  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05
)
green_KEGG <- enrichKEGG(
  gene = green_entrez, universe = background_entrez, organism = "hsa",
  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05
)
pink_KEGG <- enrichKEGG(
  gene = pink_entrez, universe = background_entrez, organism = "hsa",
  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05
)

as.data.frame(yellow_KEGG)[, c("Description","Count","p.adjust")]
# Platelet activation                       10  3.3715e-05
# ECM-receptor interaction                   8  4.318846e-04
# Focal adhesion                             8  7.488105e-03
# Neutrophil extracellular trap formation    5  4.340564e-02
# Integrin signaling                         6  4.340564e-02
# Hematopoietic cell lineage                 5  4.340564e-02

as.data.frame(green_KEGG)[, c("Description","Count","p.adjust")]
# p53 signaling pathway                4  0.003516593
# Protein processing in ER             4  0.004888743

as.data.frame(pink_KEGG)[, c("Description","Count","p.adjust")]
# Hepatitis C     4  0.02696312
# Influenza A     4  0.02849452

dotplot(yellow_KEGG, showCategory = 6, title = "Yellow Module \u2013 KEGG Pathways")
dotplot(green_KEGG,  showCategory = 2, title = "Green Module \u2013 KEGG Pathways")
dotplot(pink_KEGG,   showCategory = 2, title = "Pink Module \u2013 KEGG Pathways")

write.csv(as.data.frame(yellow_KEGG), "Yellow_KEGG.csv", row.names = FALSE)
write.csv(as.data.frame(green_KEGG),  "Green_KEGG.csv",  row.names = FALSE)
write.csv(as.data.frame(pink_KEGG),   "Pink_KEGG.csv",   row.names = FALSE)

# Interpretation:
# IMPORTANT: The Pink module's KEGG "Hepatitis C" and "Influenza A" pathway hits should
# NOT be interpreted as evidence that these patients had hepatitis C or influenza
# infection. These KEGG pathways are curated collections of shared HOST antiviral/
# pathogen-response genes (e.g. IFIT1, CXCL10, RSAD2, OAS1) that are broadly reused
# across many viral-response contexts, not dengue-specific or disease-specific markers.


# ================================================================================================
# SECTION 24 — YELLOW PLATELET GENES
# ================================================================================================
# Why: Given the strong platelet-activation/coagulation GO and KEGG signal in Yellow,
# a curated panel of platelet-associated genes was examined individually for
# correlation with DF1 status.

platelet_genes_check <- c("ITGB3","PTGS1","VWF","GP1BA","GP6","GP1BB","PTGIR","MYLK","ITGA2B","GP9")

platelet_expr <- data.frame(
  Sample = rownames(expr_final), Disease = traitData_final$Disease, DF1 = traitData_final$DF1,
  expr_final[, platelet_genes_check, drop = FALSE]
)

platelet_cor <- sapply(
  platelet_genes_check,
  function(g) cor(platelet_expr[[g]], platelet_expr$DF1, use = "pairwise.complete.obs")
)
platelet_cor
#    ITGB3     PTGS1       VWF     GP1BA       GP6     GP1BB     PTGIR      MYLK    ITGA2B       GP9
# 0.3225429 0.2859228 0.3279254 0.3089126 0.2930505 0.2516369 0.2791534 0.2219912 0.1748644 0.1110497

platelet_plot_df <- data.frame(Gene = names(platelet_cor), Correlation = as.numeric(platelet_cor))
platelet_plot_df <- platelet_plot_df[order(platelet_plot_df$Correlation), ]

ggplot(platelet_plot_df, aes(x = Correlation, y = reorder(Gene, Correlation))) +
  geom_point(size = 4) + geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Correlation of platelet-associated genes with DF1",
       x = "Pearson correlation with DF1", y = "Gene") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(10, 20, 10, 10))
ggsave("Platelet_Genes_DF1_Correlation.png", width = 8, height = 6, dpi = 300)

# Interpretation:
# All ten platelet genes show POSITIVE correlation with DF1 (range 0.11-0.33),
# consistent with the module-level GO/KEGG platelet-activation signal. However, NONE
# of these genes individually reaches genome-wide GS-FDR < 0.05 for DF1. The correct
# interpretation is that Yellow shows a coherent platelet-associated transcriptional
# program and a significant MODULE-level association with DF1, not that every
# individual platelet gene is independently significant.


# ================================================================================================
# SECTION 25 — PINK MODULE VS FEVER DAY
# ================================================================================================
# Final: Pink ME8 vs FeverDay: r = -0.377, FDR = 0.0060

pink_fever_df <- data.frame(
  Sample = rownames(MEs_nonGrey), FeverDay = traitData_final$FeverDay,
  Pink_ME = as.numeric(MEs_nonGrey[, "ME8"])
)
pink_fever_df <- pink_fever_df[!is.na(pink_fever_df$FeverDay), ]

ggplot(pink_fever_df, aes(x = FeverDay, y = Pink_ME)) +
  geom_point(size = 3) + geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Pink module eigengene across dengue fever days",
       x = "Fever day", y = "Pink module eigengene") +
  theme_classic(base_size = 14)
ggsave("Pink_Module_FeverDay.png", width = 8, height = 6, dpi = 300)

# Interpretation:
# The negative Pink-FeverDay relationship indicates that Pink module eigengene values
# tend to decrease as fever day (illness duration) increases in these dengue samples.
# This is an association observed at a single (cross-sectional) time point per patient,
# not a longitudinal trajectory, and no causal claim is made.


# ================================================================================================
# SECTION 26 — GREEN MODULE VS FEVER DAY
# ================================================================================================
# Final: Green ME5 vs FeverDay: r = -0.338, FDR = 0.0179

green_fever_df <- data.frame(
  Sample = rownames(MEs_nonGrey), FeverDay = traitData_final$FeverDay,
  Green_ME = as.numeric(MEs_nonGrey[, "ME5"])
)
green_fever_df <- green_fever_df[!is.na(green_fever_df$FeverDay), ]

ggplot(green_fever_df, aes(x = FeverDay, y = Green_ME)) +
  geom_point(size = 3) + geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Green module eigengene across dengue fever days",
       x = "Fever day", y = "Green module eigengene") +
  theme_classic(base_size = 14)
ggsave("Green_Module_FeverDay.png", width = 8, height = 6, dpi = 300)

# Interpretation:
# As with Pink, this negative association should be interpreted cautiously as a
# cross-sectional pattern (not a validated time-course trend), with no causal claim.


# ================================================================================================
# SECTION 27 — CELL-CYCLE GENES IN GREEN
# ================================================================================================
# Why: Given the strong mitotic/cell-cycle GO signal in Green, a curated panel of
# cell-cycle genes was examined for kME (green module) and gene significance (FeverDay).

cell_cycle_genes <- c("BUB1","RRM2","CDK1","DLGAP5","KIF18A","SPC25","NUF2","SKA3")

green_cellcycle <- data.frame(
  Gene = cell_cycle_genes,
  kME = kME[match(cell_cycle_genes, gene_names), "kME5"],
  GS_FeverDay = GS_FeverDay[match(cell_cycle_genes, gene_names)],
  GS_FDR = GS_FeverDay_FDR[match(cell_cycle_genes, gene_names)]
)
green_cellcycle$Priority <- abs(green_cellcycle$kME) * abs(green_cellcycle$GS_FeverDay)
green_cellcycle <- green_cellcycle[order(-green_cellcycle$Priority), ]
green_cellcycle
#    Gene       kME GS_FeverDay    GS_FDR
#    BUB1 0.8532406  -0.3830390 0.1055778   <- highest priority, but NOT individually significant
#    RRM2 0.7749403  -0.3363854 0.2239647
#    CDK1 0.5444101  -0.1981073 0.8369782
#  DLGAP5 0.5491903  -0.1848712 0.8816061
#  KIF18A 0.3696880  -0.2561869 0.6496758
#   SPC25 0.4165289  -0.2204460 0.7833529
#    NUF2 0.3909220  -0.2329757 0.7376008
#    SKA3 0.3761703  -0.0721616 0.9995016

# Heatmap of cell-cycle genes ordered by fever day
green_heatmap_df <- expr_final[dengue_samples, cell_cycle_genes, drop = FALSE]
fever_order2 <- order(traitData_final$FeverDay[match(rownames(green_heatmap_df), traitData_final$Sample)])
green_heatmap_df <- green_heatmap_df[fever_order2, , drop = FALSE]
green_heatmap_scaled <- scale(green_heatmap_df)

pheatmap(
  t(green_heatmap_scaled), cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = FALSE,
  fontsize_row = 11, main = "Cell-cycle genes in the Green module"
)

# Interpretation:
# BUB1 (kME = 0.853, GS = -0.383) is NOT individually significant after genome-wide FDR
# correction (FDR = 0.106) and should NOT be called individually significant. The Green
# module ITSELF is significantly associated with FeverDay at the eigengene level
# (Section 16) and is enriched for cell-cycle/mitotic processes (Section 22), but this
# module-level signal is not fully reducible to any single gene tested here.


# ================================================================================================
# SECTION 28 — MODULE-MODULE RELATIONSHIPS (Yellow, Green, Pink eigengene correlations)
# ================================================================================================
key_MEs <- MEs_nonGrey[, c("ME4", "ME5", "ME8")]
key_ME_cor <- cor(key_MEs, use = "pairwise.complete.obs", method = "pearson")
round(key_ME_cor, 3)
#        ME4    ME5   ME8
# ME4  1.000 -0.373 0.405
# ME5 -0.373  1.000 0.052
# ME8  0.405  0.052 1.000

key_ME_p <- matrix(NA, nrow = 3, ncol = 3, dimnames = list(colnames(key_MEs), colnames(key_MEs)))
for (i in 1:3) for (j in 1:3) {
  key_ME_p[i, j] <- corPvalueStudent(key_ME_cor[i, j], nSamples = nrow(key_MEs))
}
round(key_ME_p, 4)
# Yellow-Green p = 0.0002 | Yellow-Pink p ~ 0 | Green-Pink p = 0.6102

# Gene overlap check between modules (WGCNA modules are, by construction, disjoint
# gene sets -- each gene belongs to exactly one module)
yellow_pink_overlap <- intersect(gene_names[gene_colors == "yellow"], gene_names[gene_colors == "pink"])
length(yellow_pink_overlap)   # 0

# Interpretation:
# Yellow-Green: r = -0.373, p = 0.0002 (significant negative correlation).
# Yellow-Pink:  r =  0.405, p \u2248 0 (significant positive correlation).
# Green-Pink:   r =  0.052, p = 0.6102 (not significant).
# WGCNA modules contain MUTUALLY EXCLUSIVE genes by construction (confirmed: 0 gene
# overlap between Yellow and Pink), so zero gene overlap between modules is expected
# and does NOT preclude their eigengenes from being correlated -- module eigengene
# correlation reflects shared sample-level transcriptional dynamics, not shared genes.


# ================================================================================================
# SECTION 29 — YELLOW ROBUSTNESS ANALYSIS (leave-one-sample-out)
# ================================================================================================
# Why: To assess whether the Yellow-DF1 association is driven by one or two influential
# samples, the Yellow-DF1 Pearson correlation was recomputed 97 times, each time
# leaving out one sample.

yellow_ME <- MEs_nonGrey[, "ME4"]
df1_trait <- traitData_final$DF1
full_cor <- cor(yellow_ME, df1_trait, method = "pearson")   # 0.3911 (full-sample)

loo_cor <- sapply(seq_along(yellow_ME), function(i) {
  cor(yellow_ME[-i], df1_trait[-i], method = "pearson")
})

robustness_df <- data.frame(
  Sample = rownames(MEs_nonGrey), Disease = traitData_final$Disease,
  DF1 = df1_trait, LOO_Correlation = loo_cor
)

cat("Full Yellow-DF1 correlation:", round(full_cor, 4), "\n")     # 0.3911
cat("Minimum LOO correlation:",     round(min(loo_cor), 4), "\n") # 0.3555
cat("Maximum LOO correlation:",     round(max(loo_cor), 4), "\n") # 0.4238

write.csv(robustness_df, "Yellow_DF1_LOO_Robustness.csv", row.names = FALSE)

robustness_df$Order <- rank(robustness_df$LOO_Correlation, ties.method = "first")
ggplot(robustness_df, aes(x = Order, y = LOO_Correlation)) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = full_cor, linetype = "dashed") +
  labs(title = "Leave-one-sample-out robustness of Yellow\u2013DF1 association",
       x = "Samples ordered by leave-one-out correlation",
       y = "Yellow\u2013DF1 Pearson correlation") +
  theme_classic(base_size = 14)
ggsave("Yellow_DF1_LOO_Robustness.png", width = 8, height = 5.5, dpi = 300)

# Interpretation:
# The Yellow-DF1 association remains positive across all 97 leave-one-sample-out
# iterations within this final 97-sample dataset (range 0.3555-0.4238 vs full-sample
# 0.3911), indicating it is not driven by one or two outlier samples WITHIN this
# dataset. IMPORTANT: this robustness check does NOT validate the original six-sample
# outlier-removal decision itself (Section 9) -- it only shows stability conditional on
# that earlier, pragmatic removal having already been made.


# ================================================================================================
# SECTION 30 — OUTPUT FILES
# ================================================================================================
# All key outputs generated by this script (already written at the relevant steps
# above, consolidated here for reference):
#
#   WGCNA_Module_Sizes.csv
#   Module_Trait_Correlations.csv
#   Module_Trait_FDR.csv
#   FINAL_Module_Trait_Heatmap.png
#   Yellow_Module_Clinical_Groups.png / Green_Module_Clinical_Groups.png (Pink not saved to file in original run)
#   Hub_Genes.csv
#   Pink_Antiviral_Genes.csv
#   Pink_Antiviral_Gene_Heatmap.png
#   Yellow_GO_BP.csv / Green_GO_BP.csv / Pink_GO_BP.csv
#   Yellow_KEGG.csv / Green_KEGG.csv / Pink_KEGG.csv
#   Platelet_Genes_DF1_Correlation.png
#   Pink_Module_FeverDay.png
#   Green_Module_FeverDay.png
#   Yellow_DF1_LOO_Robustness.png / Yellow_DF1_LOO_Robustness.csv

write.csv(data.frame(Gene = yellow_genes), "Yellow_Module_Genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = green_genes),  "Green_Module_Genes.csv",  row.names = FALSE)
write.csv(data.frame(Gene = pink_genes),   "Pink_Module_Genes.csv",   row.names = FALSE)


# ================================================================================================
# SECTION 31 — FINAL BIOLOGICAL INTERPRETATION
# ================================================================================================
# Summary of findings, stated without exaggeration:
#
# 1. Yellow module (ME4): positively associated with DF1 status (r = 0.391,
#    FDR = 0.0036) and enriched for platelet activation/coagulation, wound healing,
#    and related hemostasis pathways (GO BP and KEGG concordant).
#
# 2. Green module (ME5): negatively associated with fever day (r = -0.338,
#    FDR = 0.0179) and with the Healthy indicator (r = -0.305, FDR = 0.0287); enriched
#    for mitotic/cell-cycle processes (chromosome segregation, mitotic nuclear
#    division).
#
# 3. Pink module (ME8): negatively associated with fever day (r = -0.377,
#    FDR = 0.0060); enriched for antiviral/innate-immune processes (viral process,
#    antiviral innate immune response, response to cytokine/virus).
#
# 4. Yellow and Pink module eigengenes are positively correlated (r = 0.405,
#    p \u2248 0), despite containing entirely disjoint gene sets (0 shared genes).
#
# 5. Yellow and Green module eigengenes are negatively correlated (r = -0.373,
#    p = 0.0002).
#
# 6. Green and Pink module eigengenes are not significantly correlated (r = 0.052,
#    p = 0.610).
#
# 7. Strict hub-gene analysis (|kME| >= 0.80, |GS| >= 0.20, genome-wide GS-FDR < 0.05)
#    identified four Yellow/DF1 candidate hub genes: HIST1H2BN, HIST1H2AE, HIST1H2BJ,
#    TNNC2.
#
# 8. No genes passed the strict hub criteria for Green/FeverDay or Pink/FeverDay, even
#    though both modules are significantly associated with FeverDay at the eigengene
#    (module) level.
#
# 9. Pink antiviral/innate-immune genes (ETV7, XAF1, IFIT1, RSAD2, LAMP3, IFIT2,
#    SIGLEC1, CXCL10, CCL2, OAS1) show coordinated negative directional relationships
#    with fever day, but individual-gene statistical significance varies substantially
#    (only ETV7 and XAF1 reach genome-wide GS-FDR < 0.05, and neither satisfies the
#    strict kME >= 0.80 hub criterion).


# ================================================================================================
# LIMITATIONS / INTERPRETATION CAVEATS
# ================================================================================================
# - Soft-thresholding power 12 was selected PRAGMATICALLY, not because it maximized
#   scale-free-topology R^2 -- higher powers (14, 18, 26) produced excessively sparse
#   networks (as low as 0 or ~6% non-grey genes) or no detectable modules at all.
# - Restricting the network to the 8,000 most variable genes was a pragmatic
#   dimensionality-reduction choice, not a data-driven optimum.
# - The sample-outlier cutoff of mean sample-sample correlation < 0.45 (six samples
#   removed) was a pragmatic decision; a formal sensitivity analysis of this threshold
#   (e.g., testing alternative cutoffs) would strengthen confidence in the final
#   sample set, and was not performed here. The leave-one-out robustness check in
#   Section 29 assesses stability CONDITIONAL ON this prior removal -- it does not
#   independently validate the removal decision itself.
# - The Healthy control group is small (n = 7), limiting statistical power for any
#   comparison involving Healthy samples; Healthy-related p-values/FDRs should be
#   interpreted with this limitation in mind.
# - FeverDay analyses use dengue samples only (n = 90); Healthy samples have no
#   FeverDay value by definition and are excluded from those analyses.
# - The five disease-group binary variables (DF1, DF2, DHF, DSS, Healthy) are mutually
#   exclusive and compositional (sum to 1 per sample); correlations among them and with
#   module eigengenes should be interpreted with this structure in mind.
# - Module-trait (and module-module) correlation does NOT imply causality in either
#   direction.
# - Hub-gene thresholds (|kME| >= 0.80, |GS| >= 0.20, GS-FDR < 0.05) are OPERATIONAL,
#   pre-specified criteria for this analysis, not a universal or field-standard
#   definition of "hub gene."
# - The Pink module's KEGG "Hepatitis C" and "Influenza A" pathway enrichments reflect
#   shared HOST antiviral-response gene sets and must NOT be interpreted as evidence
#   that these patients were infected with hepatitis C or influenza.
# - Platelet-associated genes in Yellow can show a coherent, biologically meaningful
#   module-level/pathway-level signal (consistent positive correlations with DF1,
#   strong GO/KEGG enrichment) without any single platelet gene individually passing
#   genome-wide FDR correction for gene significance -- module-level and gene-level
#   significance are distinct claims and should not be conflated.
# ================================================================================================
