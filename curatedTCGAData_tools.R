# Implementation of conversational tools for TCGA data access

#' List Available TCGA Cancer Types
#'
#' @return A data.frame with disease codes, names, and subtype availability
#' @export
list_cancer_types <- function() {
  data('diseaseCodes', package = "TCGAutils")
  available_codes <- TCGAutils::diseaseCodes[TCGAutils::diseaseCodes$Available == "Yes", ]
  
  result <- data.frame(
    Code = available_codes$Study.Abbreviation,
    Name = available_codes$Study.Name,
    Subtype_Data = available_codes$SubtypeData,
    stringsAsFactors = FALSE
  )
  
  message(sprintf("Found %d available TCGA cancer types", nrow(result)))
  return(result)
}


#' List Available Assays for a Disease
#'
#' @param diseaseCode TCGA disease code (e.g., "BRCA", "COAD")
#' @param version Data version (default "2.1.1")
#' @return List with disease code, version, and available assays
#' @export
list_assays_for_disease <- function(diseaseCode, version = "2.1.1") {
  if (!requireNamespace("curatedTCGAData", quietly = TRUE)) {
    stop("Package 'curatedTCGAData' is required but not installed.")
  }
  
  # Use dry.run to get available assays without downloading
  result <- curatedTCGAData::curatedTCGAData(
    diseaseCode = diseaseCode,
    assays = "*",
    version = version,
    dry.run = TRUE,
    verbose = FALSE
  )
  
  # Extract and format assay names
  assay_names <- unique(gsub(paste0("^", diseaseCode, "_"), "", result$title))
  assay_names <- gsub("-\\\\d{8}.*$", "", assay_names)
  assay_names <- sort(unique(assay_names))
  
  message(sprintf("Disease: %s (version %s)", diseaseCode, version))
  message(sprintf("Available assays: %d", length(assay_names)))
  message(paste(assay_names, collapse = ", "))
  
  return(list(
    diseaseCode = diseaseCode,
    version = version,
    available_assays = assay_names,
    count = length(assay_names)
  ))
}


#' Build MultiAssayExperiment from TCGA Data
#'
#' @param diseaseCode Character vector of TCGA disease codes
#' @param assays Character vector of assay types (supports wildcards)
#' @param version Data version (default "2.1.1")
#' @param dry.run Whether to preview before downloading (default FALSE)
#' @return MultiAssayExperiment object or preview data.frame
#' @export
build_tcga_mae <- function(diseaseCode, assays, version = "2.1.1", dry.run = FALSE) {
  if (!requireNamespace("curatedTCGAData", quietly = TRUE)) {
    stop("Package 'curatedTCGAData' is required but not installed.")
  }
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop("Package 'MultiAssayExperiment' is required but not installed.")
  }
  
  message(sprintf("Building MultiAssayExperiment for %s with assays: %s",
                  paste(diseaseCode, collapse = ", "),
                  paste(assays, collapse = ", ")))
  
  mae <- curatedTCGAData::curatedTCGAData(
    diseaseCode = diseaseCode,
    assays = assays,
    version = version,
    dry.run = dry.run,
    verbose = TRUE
  )
  
  if (!dry.run) {
    # Provide summary information
    message("\
=== MultiAssayExperiment Summary ===")
    message(sprintf("Diseases: %s", paste(diseaseCode, collapse = ", ")))
    message(sprintf("Experiments: %d", length(MultiAssayExperiment::experiments(mae))))
    message(sprintf("Experiment names: %s", 
                    paste(names(MultiAssayExperiment::experiments(mae)), collapse = ", ")))
    message(sprintf("Total samples: %d", ncol(mae)))
    message(sprintf("Patients: %d", length(unique(mae$patientID))))
    
    # Show sample types present
    message("\
=== Sample Types Present ===")
    sample_tabs <- TCGAutils::sampleTables(mae)
    for (exp_name in names(sample_tabs)) {
      message(sprintf("%s: %s", 
                      exp_name, 
                      paste(names(sample_tabs[[exp_name]]), "=", sample_tabs[[exp_name]], collapse = ", ")))
    }
  }
  
  return(mae)
}

#' Extract Primary Tumor Samples
#'
#' @param mae A MultiAssayExperiment object
#' @return MultiAssayExperiment with only primary tumor samples
#' @export
extract_primary_tumors <- function(mae) {
  if (!requireNamespace("TCGAutils", quietly = TRUE)) {
    stop("Package 'TCGAutils' is required but not installed.")
  }
  
  message("Extracting primary tumor samples...")
  
  # Show original sample counts
  orig_samples <- TCGAutils::sampleTables(mae)
  message("\
Original sample counts:")
  for (exp_name in names(orig_samples)) {
    message(sprintf("  %s: %s", 
                    exp_name, 
                    paste(names(orig_samples[[exp_name]]), "=", orig_samples[[exp_name]], collapse = ", ")))
  }
  
  # Extract primary tumors
  primary_mae <- TCGAutils::TCGAprimaryTumors(mae)
  
  # Show new sample counts
  new_samples <- TCGAutils::sampleTables(primary_mae)
  message("\
Primary tumor sample counts:")
  for (exp_name in names(new_samples)) {
    message(sprintf("  %s: %s", 
                    exp_name, 
                    paste(names(new_samples[[exp_name]]), "=", new_samples[[exp_name]], collapse = ", ")))
  }
  
  message(sprintf("\
Reduced from %d to %d samples", ncol(mae), ncol(primary_mae)))
  
  return(primary_mae)
}


#' Get Clinical Data
#'
#' @param mae A MultiAssayExperiment object
#' @param diseaseCode Optional TCGA disease code to get standard clinical variables
#' @return data.frame with clinical data
#' @export
get_clinical_data <- function(mae, diseaseCode = NULL) {
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop("Package 'MultiAssayExperiment' is required but not installed.")
  }
  
  clinical_data <- MultiAssayExperiment::colData(mae)
  
  # If disease code provided, get standard clinical variables
  if (!is.null(diseaseCode)) {
    if (!requireNamespace("TCGAutils", quietly = TRUE)) {
      warning("Package 'TCGAutils' needed for standard clinical variables")
    } else {
      standard_vars <- TCGAutils::getClinicalNames(diseaseCode)
      
      # Check which are available
      available_vars <- standard_vars[standard_vars %in% names(clinical_data)]
      
      message(sprintf("Standard clinical variables for %s:", diseaseCode))
      if (length(available_vars) > 0) {
        message(sprintf("  Available (%d): %s", length(available_vars), 
                        paste(available_vars, collapse = ", ")))
        clinical_subset <- clinical_data[, available_vars, drop = FALSE]
        return(as.data.frame(clinical_subset))
      } else {
        message("  None of the standard variables found in this dataset")
      }
    }
  }
  
  message(sprintf("Total clinical variables: %d", ncol(clinical_data)))
  message(sprintf("Total patients: %d", nrow(clinical_data)))
  
  return(as.data.frame(clinical_data))
}


#' Get Subtype Information
#'
#' @param mae A MultiAssayExperiment object
#' @return List with subtype mapping and data, or NULL if no subtypes available
#' @export
get_subtype_info <- function(mae) {
  if (!requireNamespace("TCGAutils", quietly = TRUE)) {
    stop("Package 'TCGAutils' is required but not installed.")
  }
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop("Package 'MultiAssayExperiment' is required but not installed.")
  }
  
  subtype_map <- TCGAutils::getSubtypeMap(mae)
  
  if (nrow(subtype_map) > 0) {
    message("Available subtype information:")
    print(subtype_map)
    
    # Extract subtype columns from colData
    subtype_cols <- unique(subtype_map[[2]])
    clinical_data <- MultiAssayExperiment::colData(mae)
    
    available_subtype_cols <- subtype_cols[subtype_cols %in% names(clinical_data)]
    
    if (length(available_subtype_cols) > 0) {
      message(sprintf("\
Subtype columns in data: %s", 
                      paste(available_subtype_cols, collapse = ", ")))
      subtype_data <- clinical_data[, available_subtype_cols, drop = FALSE]
      
      # Show summary of subtype values
      for (col in available_subtype_cols) {
        vals <- table(subtype_data[[col]], useNA = "ifany")
        message(sprintf("\
%s:", col))
        print(vals)
      }
      
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

#' Describe Assay
#'
#' @param mae A MultiAssayExperiment object
#' @param assay_name Name of the assay to describe
#' @return List with assay information
#' @export
describe_assay <- function(mae, assay_name) {
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop("Package 'MultiAssayExperiment' is required but not installed.")
  }
  
  if (!assay_name %in% names(MultiAssayExperiment::experiments(mae))) {
    stop(sprintf("Assay '%s' not found.\
Available assays: %s",
                 assay_name, 
                 paste(names(MultiAssayExperiment::experiments(mae)), collapse = ", ")))
  }
  
  assay_data <- MultiAssayExperiment::experiments(mae)[[assay_name]]
  
  info <- list(
    name = assay_name,
    class = class(assay_data)[1],
    dimensions = dim(assay_data),
    n_features = nrow(assay_data),
    n_samples = ncol(assay_data)
  )
  
  message(sprintf("\
=== Assay: %s ===", info$name))
  message(sprintf("Class: %s", info$class))
  message(sprintf("Dimensions: %d features × %d samples", info$n_features, info$n_samples))
  
  # Show feature names
  if (info$n_features > 0) {
    n_preview <- min(10, info$n_features)
    info$feature_names_preview <- head(rownames(assay_data), n_preview)
    message(sprintf("\
First %d features:", n_preview))
    message(paste("  ", info$feature_names_preview, collapse = "\
"))
  }
  
  # Show sample names
  if (info$n_samples > 0) {
    n_preview <- min(5, info$n_samples)
    info$sample_names_preview <- head(colnames(assay_data), n_preview)
    message(sprintf("\
First %d samples:", n_preview))
    message(paste("  ", info$sample_names_preview, collapse = "\
"))
  }
  
  # Additional info for specific classes
  if (methods::is(assay_data, "RaggedExperiment")) {
    info$metadata_columns <- names(S4Vectors::mcols(assay_data))
    message(sprintf("\
Metadata columns (%d): %s", 
                    length(info$metadata_columns),
                    paste(head(info$metadata_columns, 10), collapse = ", ")))
  } else if (methods::is(assay_data, "SummarizedExperiment")) {
    info$assay_names <- SummarizedExperiment::assayNames(assay_data)
    info$col_data_columns <- names(SummarizedExperiment::colData(assay_data))
    message(sprintf("\
Assay matrices: %s", paste(info$assay_names, collapse = ", ")))
    message(sprintf("Column data variables (%d): %s", 
                    length(info$col_data_columns),
                    paste(head(info$col_data_columns, 10), collapse = ", ")))
  }
  
  return(invisible(info))
}


#' Split Assays by Sample Type
#'
#' @param mae A MultiAssayExperiment object
#' @param sample_codes Character vector of TCGA sample codes (e.g., "01" for primary tumor)
#' @return MultiAssayExperiment with assays split by sample type
#' @export
split_by_sample_type <- function(mae, sample_codes = c("01", "10", "11")) {
  if (!requireNamespace("TCGAutils", quietly = TRUE)) {
    stop("Package 'TCGAutils' is required but not installed.")
  }
  
  # Load sample type reference
  data('sampleTypes', package = "TCGAutils")
  
  message("Sample codes requested:")
  selected_types <- TCGAutils::sampleTypes[TCGAutils::sampleTypes$Code %in% sample_codes, ]
  print(selected_types[, c("Code", "Definition")])
  
  message("\
Splitting assays by sample type...")
  result <- TCGAutils::TCGAsplitAssays(mae, sample_codes)
  
  message(sprintf("\
Original experiments: %d", length(MultiAssayExperiment::experiments(mae))))
  message(sprintf("Split experiments: %d", length(MultiAssayExperiment::experiments(result))))
  
  return(result)
}


#' Provide Conversational Help
#'
#' @param topic Optional topic for help (e.g., "assays", "diseases", "workflow")
#' @return Prints helpful information
#' @export
tcga_help <- function(topic = NULL) {
  if (is.null(topic)) {
    message("=== curatedTCGAData Help ===\
")
    message("Available help topics:")
    message("  'diseases' - List available cancer types")
    message("  'assays' - Explain different data types")
    message("  'workflow' - Basic analysis workflow")
    message("  'samples' - Understanding TCGA sample types")
    message("  'versions' - Data version information")
    message("\
Usage: tcga_help('topic')")
    return(invisible(NULL))
  }
  
  topic <- tolower(topic)
  
  if (topic == "diseases") {
    message("=== TCGA Cancer Types ===\
")
    message("Use list_cancer_types() to see all 33 available cancer types.")
    message("\
Common abbreviations:")
    message("  BRCA - Breast invasive carcinoma")
    message("  LUAD - Lung adenocarcinoma")
    message("  COAD - Colon adenocarcinoma")
    message("  OV - Ovarian serous cystadenocarcinoma")
    message("  GBM - Glioblastoma multiforme")
    
  } else if (topic == "assays") {
    message("=== TCGA Assay Types ===\
")
    message("Gene Expression:")
    message("  RNASeq2Gene - RSEM TPM gene expression")
    message("  RNASeq2GeneNorm - Normalized RSEM expression")
    message("  mRNAArray - Microarray gene expression")
    message("\
Mutations:")
    message("  Mutation - Somatic mutation calls")
    message("\
Copy Number:")
    message("  CNASNP - Copy number from SNP array")
    message("  CNVSNP - Germline copy number variants")
    message("  GISTIC_AllByGene - Gene-level copy number")
    message("\
Epigenetics:")
    message("  Methylation_methyl450 - 450K methylation array")
    message("\
Proteins:")
    message("  RPPAArray - Reverse phase protein array")
    message("\
Tip: Use wildcards like 'CN*' for all copy number data")
    
  } else if (topic == "workflow") {
    message("=== Basic Workflow ===\
")
    message("1. List available cancers:")
    message("   cancers <- list_cancer_types()")
    message("\
2. Check available assays:")
    message("   list_assays_for_disease('BRCA')")
    message("\
3. Build MultiAssayExperiment:")
    message("   mae <- build_tcga_mae('BRCA', c('RNASeq2Gene', 'Mutation'))")
    message("\
4. Extract primary tumors:")
    message("   mae_primary <- extract_primary_tumors(mae)")
    message("\
5. Get clinical data:")
    message("   clinical <- get_clinical_data(mae_primary, 'BRCA')")
    
  } else if (topic == "samples") {
    message("=== TCGA Sample Types ===\
")
    message("TCGA barcodes contain sample type codes:")
    message("  01 - Primary Solid Tumor")
    message("  10 - Blood Derived Normal")
    message("  11 - Solid Tissue Normal")
    message("  02 - Recurrent Solid Tumor")
    message("  06 - Metastatic")
    message("\
Use TCGAprimaryTumors() to filter to primary tumors")
    message("Use split_by_sample_type() to separate sample types")
    
  } else if (topic == "versions") {
    message("=== Data Versions ===\
")
    message("Current version: 2.1.1 (recommended)")
    message("  - Log2 RPM miRNA expression")
    message("  - Log2 normalized RSEM RNA-seq")
    message("  - Corrected subtype data")
    message("\
Previous versions: 2.1.0, 2.0.1, 1.1.38")
    message("Specify version in build_tcga_mae(version = '2.1.1')")
    
  } else {
    message(sprintf("Unknown topic: %s", topic))
    message("Available topics: diseases, assays, workflow, samples, versions")
  }
  
  return(invisible(NULL))
}


