## 5-Minute Setup

### 1. Install Required Packages (2 minutes)
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "curatedTCGAData",
  "MultiAssayExperiment", 
  "TCGAutils"
))
```

### 2. Load the Tool (30 seconds)
```r
source("curatedTCGAData_tools.R")
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
```

### 3. Your First Analysis (2 minutes)
```r
# Step 1: See what's available
cancers <- list_cancer_types()
head(cancers)

# Step 2: Pick a cancer and see its data
list_assays_for_disease("BRCA")

# Step 3: Get the data
mae <- build_tcga_mae(
  diseaseCode = "BRCA",
  assays = "RNASeq2Gene",
  dry.run = FALSE
)

# Step 4: Filter to primary tumors
mae <- extract_primary_tumors(mae)

# Step 5: Get clinical data
clinical <- get_clinical_data(mae, "BRCA")
```

🎉 **You're ready to analyze TCGA data!**

---

## Common Use Cases

### Use Case 1: Differential Expression
**Goal**: Compare gene expression between cancer subtypes

```r
# Get RNA-seq data
mae <- build_tcga_mae("BRCA", "RNASeq2Gene")
mae <- extract_primary_tumors(mae)

# Get subtypes
subtypes <- get_subtype_info(mae)

# Extract expression matrix
expr <- assay(experiments(mae)[[1]])

# Now use DESeq2, edgeR, or limma for DE analysis
```

### Use Case 2: Survival Analysis
**Goal**: Correlate gene expression with survival

```r
# Get data
mae <- build_tcga_mae("OV", "RNASeq2Gene")

# Get clinical variables
clinical <- get_clinical_data(mae, "OV")

# Key survival variables
survival_vars <- c("days_to_death", "days_to_last_followup", "vital_status")

# Now use survival package
```

### Use Case 3: Multi-Omic Integration
**Goal**: Integrate RNA-seq, copy number, and methylation

```r
# Get multiple assays
mae <- build_tcga_mae(
  diseaseCode = "BRCA",
  assays = c("RNASeq2Gene", "GISTIC_AllByGene", "Methylation_methyl450")
)

mae <- extract_primary_tumors(mae)

# Find common samples
common <- Reduce(intersect, lapply(experiments(mae), colnames))
mae <- mae[, common, ]

# Now use MOFA2, iCluster, or mixOmics
```

### Use Case 4: Pan-Cancer Analysis
**Goal**: Compare multiple cancer types

```r
# Get data for multiple cancers
mae <- build_tcga_mae(
  diseaseCode = c("BRCA", "OV", "UCEC"),  # 3 gynecologic cancers
  assays = "RNASeq2Gene"
)

# Get clinical data
clinical <- get_clinical_data(mae)

# Compare by disease code
table(clinical$admin.disease_code)
```

---

## Troubleshooting

### Problem: Package installation fails
**Solution**:
```r
# Update BiocManager
install.packages("BiocManager")

# Try installing one package at a time
BiocManager::install("curatedTCGAData")
BiocManager::install("MultiAssayExperiment")
BiocManager::install("TCGAutils")
```

### Problem: "Disease code not found"
**Solution**:
```r
# Check valid codes
cancers <- list_cancer_types()
print(cancers$Code)
```

### Problem: "Assay not available"
**Solution**:
```r
# Check what's available
list_assays_for_disease("YOUR_CODE")
```

### Problem: Download is taking forever
**Solution**:
```r
# Preview first
build_tcga_mae("BRCA", "*", dry.run = TRUE)

# Download one assay at a time
mae1 <- build_tcga_mae("BRCA", "RNASeq2Gene")
mae2 <- build_tcga_mae("BRCA", "Mutation")
```

---

## Natural Language Examples

### Example 1: Starting from Scratch
```
User: "I want to analyze breast cancer"
→ list_cancer_types() 
→ Show BRCA is available

User: "What data can I get?"
→ list_assays_for_disease("BRCA")
→ Show available assays

User: "Get me gene expression and mutations"
→ build_tcga_mae("BRCA", c("RNASeq2Gene", "Mutation"))
→ Build and show summary

User: "Only tumor samples please"
→ extract_primary_tumors(mae)
→ Filter and report
```


### Example 2: Specific Request
```
User: "Get lung cancer RNA-seq data"
→ Clarify: LUAD (adenocarcinoma) or LUSC (squamous)?

User: "Both"
→ build_tcga_mae(c("LUAD", "LUSC"), "RNASeq2Gene")
→ Build combined dataset

User: "Show me the clinical variables"
→ get_clinical_data(mae)
→ Display clinical data
```

### Example 3: Exploring
```
User: "What cancers are available?"
→ list_cancer_types()
→ Show all 33 cancer types

User: "I'm interested in kidney cancer"
→ list_assays_for_disease("KIRC")
→ Show KIRC assays

User: "Get everything"
→ Warn about size, suggest specific assays
→ Or use dry.run first
```

---

## Next Steps After Loading Data

Once you have your MultiAssayExperiment:

### 1. Explore the Data
```r
# Describe experiments
names(experiments(mae))

# Describe specific assay
describe_assay(mae, "BRCA_RNASeq2Gene-20160128")

# Check sample overlap
upsetSamples(mae)
```

### 2. Extract Data for Analysis
```r
# Get expression matrix
expr <- assay(experiments(mae)[[1]])

# Get clinical data
clinical <- colData(mae)

# Get mutation data as GRangesList
mutations <- as(experiments(mae)[["Mutation"]], "GRangesList")
```

### 3. Filter and Subset
```r
# Primary tumors only
mae <- extract_primary_tumors(mae)

# Specific patients
mae_subset <- mae[, mae$patientID %in% patients_of_interest, ]

# Specific genes
mae_genes <- mae[c("TP53", "BRCA1", "BRCA2"), , ]
```

### 4. Perform Analysis
```r
# Differential expression (use DESeq2, edgeR, limma)
# Survival analysis (use survival package)
# Pathway analysis (use fgsea, clusterProfiler)
# Integration (use MOFA2, mixOmics)
```

---

## Cheat Sheet

| Want to... | Use this function |
|-----------|------------------|
| See available cancers | `list_cancer_types()` |
| Check assays for a cancer | `list_assays_for_disease("CODE")` |
| Get the data | `build_tcga_mae("CODE", "ASSAYS")` |
| Preview before download | `build_tcga_mae(..., dry.run=TRUE)` |
| Get only tumors | `extract_primary_tumors(mae)` |
| Get clinical data | `get_clinical_data(mae, "CODE")` |
| Get subtypes | `get_subtype_info(mae)` |
| Explore an assay | `describe_assay(mae, "NAME")` |
| Get help | `tcga_help("TOPIC")` |

## Key Reminders

✅ **DO**:
- Start with `dry.run = TRUE`
- Use wildcards (`CN*`, `RNA*`) for groups of assays
- Filter to primary tumors for most analyses
- Check sample overlap with `upsetSamples()`
- Use standard clinical variables via `getClinicalNames()`

❌ **DON'T**:
- Download all assays at once without previewing
- Forget to filter sample types
- Assume all samples overlap across experiments
- Mix data versions

---

## Getting Help

1. **Built-in help**: `tcga_help()` or `tcga_help("topic")`
2. **Examples**: See `examples.R`
3. **Documentation**: See `SKILL.md` and `README.md`
4. **Package help**: `?curatedTCGAData` in R
5. **Bioconductor**: https://bioconductor.org/packages/curatedTCGAData

---

## You're Ready! 🚀

Start with:
```r
source("curatedTCGAData_tools.R")
list_cancer_types()
```

Then follow your research question!
`
