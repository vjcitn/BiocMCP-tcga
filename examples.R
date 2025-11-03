# This demonstrates conversational workflows for accessing TCGA data

# Load the tools
source("curatedTCGAData_tools.R")

# Load required packages
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)

# ============================================================================
# Example 1: Exploring Available Data
# ============================================================================

cat("\
=== Example 1: What cancer types are available? ===\
")
cancers <- list_cancer_types()
head(cancers, 10)

cat("\
=== What assays are available for breast cancer? ===\
")
brca_assays <- list_assays_for_disease("BRCA")

# ============================================================================
# Example 2: Building a Simple MultiAssayExperiment
# ============================================================================

cat("\
=== Example 2: Get breast cancer RNA-seq and mutation data ===\
")

# First, do a dry run to see what will be downloaded
cat("\
Dry run (preview):\
")
preview <- build_tcga_mae(
  diseaseCode = "BRCA",
  assays = c("RNASeq2Gene", "Mutation"),
  version = "2.1.1",
  dry.run = TRUE
)
print(preview)

# Actually build the MultiAssayExperiment
cat("\
Building MultiAssayExperiment...\
")
brca_mae <- build_tcga_mae(
  diseaseCode = "BRCA",
  assays = c("RNASeq2Gene", "Mutation"),
  version = "2.1.1",
  dry.run = FALSE
)

# ============================================================================
# Example 3: Filtering to Primary Tumors
# ============================================================================

cat("\
=== Example 3: Filter to primary tumor samples ===\
")
brca_primary <- extract_primary_tumors(brca_mae)

# ============================================================================
# Example 4: Extracting Clinical Data
# ============================================================================

cat("\
=== Example 4: Get clinical data ===\
")

# Get all clinical data
all_clinical <- get_clinical_data(brca_primary)
cat(sprintf("\
Total clinical variables: %d\
", ncol(all_clinical)))

# Get standard clinical variables for BRCA
standard_clinical <- get_clinical_data(brca_primary, diseaseCode = "BRCA")
head(standard_clinical)

# ============================================================================
# Example 5: Exploring Subtype Information
# ============================================================================

cat("\
=== Example 5: Get breast cancer subtypes ===\
")
subtype_info <- get_subtype_info(brca_primary)

if (!is.null(subtype_info)) {
  # Show distribution of PAM50 subtypes
  if ("PAM50" %in% names(subtype_info$data)) {
    cat("\
PAM50 subtype distribution:\
")
    print(table(subtype_info$data$PAM50))
  }
}

# ============================================================================
# Example 6: Describing Individual Assays
# ============================================================================

cat("\
=== Example 6: Examine the RNA-seq assay ===\
")
describe_assay(brca_primary, "BRCA_RNASeq2Gene-20160128")

cat("\
=== Examine the mutation assay ===\
")
describe_assay(brca_primary, "BRCA_Mutation-20160128")

# ============================================================================
# Example 7: Multi-Cancer Comparison
# ============================================================================

cat("\
=== Example 7: Compare lung adenocarcinoma and squamous cell ===\
")

# Build combined dataset
lung_mae <- build_tcga_mae(
  diseaseCode = c("LUAD", "LUSC"),
  assays = "RNASeq2Gene",
  version = "2.1.1",
  dry.run = FALSE
)

# Get clinical data to compare
lung_clinical <- get_clinical_data(lung_mae)

# Check sample sizes
cat("\
Sample counts by disease:\
")
print(table(lung_clinical$admin.disease_code))

# ============================================================================
# Example 8: Using Wildcards for Assay Selection
# ============================================================================

cat("\
=== Example 8: Get all copy number data for colorectal cancer ===\
")

# Use wildcard to get all copy number assays
coad_cn <- build_tcga_mae(
  diseaseCode = "COAD",
  assays = "CN*",  # Matches CNASNP, CNVSNP, CNASeq
  version = "2.1.1",
  dry.run = FALSE
)

cat("\
Copy number experiments:\
")
print(names(experiments(coad_cn)))

# ============================================================================
# Example 9: Splitting by Sample Type
# ============================================================================

cat("\
=== Example 9: Separate tumor and normal samples ===\
")

# Get both tumor and normal samples for COAD
coad_split <- split_by_sample_type(
  coad_cn,
  sample_codes = c("01", "10", "11")  # Primary tumor, blood normal, solid tissue normal
)

cat("\
Experiments after splitting:\
")
print(names(experiments(coad_split)))

# ============================================================================
# Example 10: Getting Help
# ============================================================================

cat("\
=== Example 10: Using the help system ===\
")

# General help
tcga_help()

# Specific topics
cat("\
")
tcga_help("assays")

cat("\
")
tcga_help("workflow")

# ============================================================================
# Conversational Workflow Example
# ============================================================================

cat("\
=== Conversational Workflow: Ovarian Cancer Analysis ===\
")

# "I want to analyze ovarian cancer"
cat("\
Step 1: Check available assays for OV\
")
ov_assays <- list_assays_for_disease("OV")

# "Get me RNA-seq, copy number, and mutation data"
cat("\
Step 2: Build MultiAssayExperiment\
")
ov_mae <- build_tcga_mae(
  diseaseCode = "OV",
  assays = c("RNASeq2Gene", "GISTIC_AllByGene", "Mutation"),
  version = "2.1.1",
  dry.run = FALSE
)

# "Show me only primary tumors"
cat("\
Step 3: Filter to primary tumors\
")
ov_primary <- extract_primary_tumors(ov_mae)

# "What subtypes are available?"
cat("\
Step 4: Check for subtype information\
")
ov_subtypes <- get_subtype_info(ov_primary)

# "Get the clinical data"
cat("\
Step 5: Extract clinical variables\
")
ov_clinical <- get_clinical_data(ov_primary, diseaseCode = "OV")

cat("\
Ready for downstream analysis!\
")
cat("Now you could perform:\
")
cat("  - Differential expression analysis between subtypes\
")
cat("  - Survival analysis using clinical data\
")
cat("  - Integration of RNA-seq and copy number\
")
cat("  - Mutation enrichment analysis\
")

# 

# ============================================================================
# Tips for MCP Integration
# ============================================================================

cat("\
=== Tips for MCP Tool Integration ===\
")
cat("
1. Natural Language Mapping:
   User says: 'breast cancer' → Use diseaseCode = 'BRCA'
   User says: 'gene expression' → Suggest 'RNASeq2Gene'
   User says: 'mutations' → Use 'Mutation'

2. Progressive Disclosure:
   - Start with list_cancer_types() to show options
   - Then list_assays_for_disease() for specific cancer
   - Use dry.run=TRUE to preview before downloading
   - Build the actual dataset last

3. Context Awareness:
   - Remember which MAE object is currently active
   - Suggest logical next steps based on workflow
   - Offer to filter, extract, or describe as needed

4. Error Prevention:
   - Always validate disease codes against list_cancer_types()
   - Check assay availability before building
   - Warn about large downloads

5. Helpful Responses:
   - Explain what each assay type contains
   - Clarify sample type codes (01, 10, 11, etc.)
   - Suggest relevant TCGAutils functions
   - Link to documentation when appropriate
")


