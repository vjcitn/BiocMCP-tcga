## Overview
This skill provides conversational tools for working with TCGA (The Cancer Genome Atlas) data through the curatedTCGAData Bioconductor package. It simplifies the process of building MultiAssayExperiment objects by allowing users to conversationally specify tumor types and assay types.

## Core Functions

### Tool 1: List Available Cancer Types
**Function**: `list_cancer_types`
**Purpose**: Show available TCGA disease codes and their descriptions

**Implementation**:
```r
list_cancer_types <- function() {
  data('diseaseCodes', package = "TCGAutils")
  available_codes <- diseaseCodes[diseaseCodes$Available == "Yes", ]
  
  # Format output nicely
  result <- data.frame(
    Code = available_codes$Study.Abbreviation,
    Name = available_codes$Study.Name,
    Subtype_Data = available_codes$SubtypeData,
    stringsAsFactors = FALSE
  )
  
  return(result)
}
```

**Usage**: "What cancer types are available?" or "Show me all TCGA disease codes"

---

### Tool 2: List Available Assays
**Function**: `list_assays_for_disease`
**Purpose**: Show what data types (assays) are available for a specific disease code

**Parameters**:
- `diseaseCode`: TCGA disease code (e.g., "BRCA", "COAD")
- `version`: Data version (default "2.1.1")

**Implementation**:
```r
list_assays_for_disease <- function(diseaseCode, version = "2.1.1") {
  library(curatedTCGAData)
  
  # Use dry.run to get available assays without downloading
  result <- curatedTCGAData(
    diseaseCode = diseaseCode,
    assays = "*",
    version = version,
    dry.run = TRUE,
    verbose = FALSE
  )
  
  # Extract and format assay names
  assay_names <- unique(gsub(paste0("^", diseaseCode, "_"), "", result$title))
  assay_names <- gsub("-\\\\d{8}.*$", "", assay_names)
  
  return(list(
    diseaseCode = diseaseCode,
    version = version,
    available_assays = sort(unique(assay_names)),
    count = length(unique(assay_names))
  ))
}
```

**Usage**: "What assays are available for breast cancer?" or "Show me data types for BRCA"

---

### Tool 3: Build MultiAssayExperiment
**Function**: `build_tcga_mae`
**Purpose**: Create a MultiAssayExperiment with specified disease codes and assays

**Parameters**:
- `diseaseCode`: Character vector of TCGA disease codes
- `assays`: Character vector of assay types (supports wildcards like "CN*", "RNA*")
- `version`: Data version (default "2.1.1")
- `dry.run`: Whether to preview before downloading (default FALSE)

**Implementation**:
```r
build_tcga_mae <- function(diseaseCode, assays, version = "2.1.1", dry.run = FALSE) {
  library(curatedTCGAData)
  library(MultiAssayExperiment)
  
  message(sprintf("Building MultiAssayExperiment for %s with assays: %s",
                  paste(diseaseCode, collapse = ", "),
                  paste(assays, collapse = ", ")))
  
  mae <- curatedTCGAData(
    diseaseCode = diseaseCode,
    assays = assays,
    version = version,
    dry.run = dry.run,
    verbose = TRUE
  )
  
  if (!dry.run) {
    # Provide summary information
    summary_info <- list(
      diseases = diseaseCode,
      n_experiments = length(experiments(mae)),
      experiment_names = names(experiments(mae)),
      n_samples = ncol(mae),
      n_features = nrow(mae)
    )
    
    message("\
MultiAssayExperiment Summary:")
    message(sprintf("  Diseases: %s", paste(diseaseCode, collapse = ", ")))
    message(sprintf("  Experiments: %d", summary_info$n_experiments))
    message(sprintf("  Total samples: %d", summary_info$n_samples))
    
    return(mae)
  } else {
    return(mae)
  }
}
```

**Usage**: "Get breast cancer RNA-seq and mutation data" or "Build a MAE with BRCA and COAD with all copy number assays"

---

### Tool 4: Get Primary Tumors Only
**Function**: `extract_primary_tumors`
**Purpose**: Filter MultiAssayExperiment to include only primary tumor samples

**Parameters**:
- `mae`: A MultiAssayExperiment object

**Implementation**:
```r
extract_primary_tumors <- function(mae) {
  library(TCGAutils)
  
  message("Extracting primary tumor samples...")
  primary_mae <- TCGAprimaryTumors(mae)
  
  # Show sample counts before and after
  orig_samples <- sampleTables(mae)
  new_samples <- sampleTables(primary_mae)
  
  message("\
Sample counts:")
  message("Original:")
  print(orig_samples)
  message("\
Primary tumors only:")
  print(new_samples)
  
  return(primary_mae)
}
```

**Usage**: "Filter to primary tumors only" or "Keep only the primary tumor samples"

---

### Tool 5: Get Clinical Variables
**Function**: `get_clinical_data`
**Purpose**: Extract key clinical variables from the MultiAssayExperiment

**Parameters**:
- `mae`: A MultiAssayExperiment object
- `diseaseCode`: TCGA disease code to get standard clinical variables

**Implementation**:
```r
get_clinical_data <- function(mae, diseaseCode = NULL) {
  library(TCGAutils)
  
  clinical_data <- colData(mae)
  
  # If disease code provided, get standard clinical variables
  if (!is.null(diseaseCode)) {
    standard_vars <- getClinicalNames(diseaseCode)
    
    # Check which are available
    available_vars <- standard_vars[standard_vars %in% names(clinical_data)]
    
    message(sprintf("Standard clinical variables for %s:", diseaseCode))
    message(sprintf("  Available: %s", paste(available_vars, collapse = ", ")))
    
    if (length(available_vars) > 0) {
      clinical_subset <- clinical_data[, available_vars, drop = FALSE]
      return(as.data.frame(clinical_subset))
    }
  }
  
  return(as.data.frame(clinical_data))
}
```

**Usage**: "Get clinical data" or "Show me the clinical variables for BRCA"

---

### Tool 6: Get Subtype Information
**Function**: `get_subtype_info`
**Purpose**: Extract cancer subtype information from the MultiAssayExperiment

**Parameters**:
- `mae`: A MultiAssayExperiment object

**Implementation**:
```r
get_subtype_info <- function(mae) {
  library(TCGAutils)
  
  subtype_map <- getSubtypeMap(mae)
  
  if (nrow(subtype_map) > 0) {
    message("Available subtype information:")
    print(subtype_map)
    
    # Extract subtype columns from colData
    subtype_cols <- unique(subtype_map[[2]])
    clinical_data <- colData(mae)
    
    available_subtype_cols <- subtype_cols[subtype_cols %in% names(clinical_data)]
    
    if (length(available_subtype_cols) > 0) {
      subtype_data <- clinical_data[, available_subtype_cols, drop = FALSE]
      return(list(
        map = subtype_map,
        data = as.data.frame(subtype_data)
      ))
    }
  } else {
    message("No subtype information available for this dataset.")
    return(NULL)
  }
}
```

**Usage**: "What subtypes are available?" or "Show me the subtype classifications"

---

### Tool 7: Describe Assay
**Function**: `describe_assay`
**Purpose**: Get detailed information about a specific assay in the MultiAssayExperiment

**Parameters**:
- `mae`: A MultiAssayExperiment object
- `assay_name`: Name of the assay to describe

**Implementation**:
```r
describe_assay <- function(mae, assay_name) {
  if (!assay_name %in% names(experiments(mae))) {
    stop(sprintf("Assay '%s' not found. Available assays: %s",
                 assay_name, paste(names(experiments(mae)), collapse = ", ")))
  }
  
  assay_data <- experiments(mae)[[assay_name]]
  
  info <- list(
    name = assay_name,
    class = class(assay_data)[1],
    dimensions = dim(assay_data),
    n_features = nrow(assay_data),
    n_samples = ncol(assay_data),
    feature_names_preview = head(rownames(assay_data), 10),
    sample_names_preview = head(colnames(assay_data), 10)
  )
  
  # Additional info for specific classes
  if (is(assay_data, "RaggedExperiment")) {
    info$metadata_columns <- names(mcols(assay_data))
  } else if (is(assay_data, "SummarizedExperiment")) {
    info$assay_names <- assayNames(assay_data)
    info$col_data_columns <- names(colData(assay_data))
  }
  
  message(sprintf("\
Assay: %s", info$name))
  message(sprintf("Class: %s", info$class))
  message(sprintf("Dimensions: %d features × %d samples", info$n_features, info$n_samples))
  
  return(info)
}
```

**Usage**: "Describe the mutation assay" or "What's in the RNASeq2Gene experiment?"

---

## Conversation Flow Examples

### Example 1: Basic Workflow
```
User: "I want to work with breast cancer data"
Assistant: [calls list_assays_for_disease("BRCA")]
         "Breast cancer (BRCA) has the following assays available: [lists assays]
          Which assays would you like to include?"

User: "Get me RNA-seq and mutation data"
Assistant: [calls build_tcga_mae(diseaseCode="BRCA", assays=c("RNASeq2Gene", "Mutation"))]
         "I've built your MultiAssayExperiment with 2 experiments..."

User: "Filter to primary tumors only"
Assistant: [calls extract_primary_tumors(mae)]
         "Filtered to primary tumors. Original had X samples, now Y samples..."
```

### Example 2: Multi-Cancer Comparison
```
User: "Compare lung adenocarcinoma and lung squamous cell carcinoma"
Assistant: [calls list_assays_for_disease for both LUAD and LUSC]
         "Both LUAD and LUSC have these common assays: [lists]
          Would you like to build a combined dataset?"

User: "Yes, get RNA-seq for both"
Assistant: [calls build_tcga_mae(diseaseCode=c("LUAD", "LUSC"), assays="RNASeq2Gene")]
         "Created MultiAssayExperiment with RNA-seq data for both cancer types..."
```

### Example 3: Exploring Available Data
```
User: "What cancer types can I analyze?"
Assistant: [calls list_cancer_types()]
         "Here are all 33 available TCGA cancer types: [shows table]"

User: "What data is available for colorectal cancer?"
Assistant: [calls list_assays_for_disease("COAD")]
         "COAD (Colon adenocarcinoma) has these assays: [lists]"
```

## Common Assay Patterns

Users often request data using natural language. Here are mappings:

- **"RNA-seq"** → `RNASeq2Gene` or `RNASeq2GeneNorm`
- **"gene expression"** → `RNASeq2Gene`, `mRNAArray`
- **"mutations"** → `Mutation`
- **"copy number"** → `CNASNP`, `CNVSNP`, `CNASeq`, or `CN*` (wildcard)
- **"methylation"** → `Methylation_methyl450`, `Methylation_methyl27`
- **"protein"** → `RPPAArray`
- **"miRNA"** → `miRNASeqGene`
- **"GISTIC"** → `GISTIC_AllByGene`, `GISTIC_ThresholdedByGene`

## Best Practices

1. **Always start with dry.run**: Show users what will be downloaded before actually downloading
2. **Offer guidance**: When users specify a disease, show available assays
3. **Use wildcards intelligently**: Suggest `CN*` for all copy number data, `RNA*` for all RNA data
4. **Explain data types**: When returning an MAE, explain what each experiment contains
5. **Suggest filtering**: Remind users they can filter to primary tumors or specific sample types
6. **Reference TCGAutils**: Mention that TCGAutils has additional helper functions for working with the data

## Error Handling

Common issues and solutions:

1. **Invalid disease code**: Check against `diseaseCodes$Study.Abbreviation`
2. **Assay not available**: Use `list_assays_for_disease` first
3. **Version mismatch**: Default to latest version (2.1.1) unless specified
4. **Large downloads**: Warn users about data size for multiple diseases/assays
5. **Memory constraints**: Suggest downloading one assay at a time for large datasets

## Integration with Analysis

After building a MAE, suggest next steps:
- "Would you like to filter to primary tumors?"
- "Shall I extract the clinical data?"
- "Would you like to see the subtype information?"
- "Ready to perform differential expression analysis?"

## Version Information

Default version: **2.1.1** (latest)
- Includes log2 RPM miRNA expression values
- Includes log2 normalized RSEM for RNA-seq
- Corrected subtype data for OV and SKCM

Previous versions available: 2.1.0, 2.0.1, 1.1.38

## Required Packages

- `curatedTCGAData`
- `MultiAssayExperiment`
- `TCGAutils`
- `SummarizedExperiment`
- `RaggedExperiment` (for mutation and copy number data)

