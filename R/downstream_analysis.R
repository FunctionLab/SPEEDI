
#' Perform downstream analyses (RNA)
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param reference_tissue Reference tissue
#' @param species Species (human or mouse)
#' @param metadata_df Data frame containing metadata for samples. Rownames should be sample names and column names should be metadata attributes with two classes (e.g., condition: disease and control)
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]). Directory will be created if it doesn't already exist.
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A list containing differential expression analyses
#' @examples
#' \dontrun{differential_expression_results <- RunDE_RNA(sc_obj = sc_obj, metadata_df = metadata_df)}
#' @export
run_downstream_analyses_RNA <- function(sc_obj, reference_tissue, species, metadata_df, output_dir = getwd(), log_flag = log_flag) {
  # We run differential expression on each metadata attribute provided by the user
  differential_expression_results <- RunDE_RNA(sc_obj, metadata_df = metadata_df, output_dir = output_dir, log_flag = log_flag)
  # Next, if species is human, we will run functional module discovery using the DE results
  # TODO: Pick appropriate networks depending on cell type
  FMD_results <- list()
  if(species == "human") {
    index <- 1
    for(current_de_result in differential_expression_results) {
      for(current_cell_type in unique(current_de_result$Cell_Type)) {
        # Subset to DE results for our current cell type
        cell_specific_de_result <- current_de_result[current_de_result$Cell_Type == current_cell_type,]
        # Grab HB networks associated with cell type and reference tissue
        hb_networks <- grab_hb_networks(current_cell_type, reference_tissue)
        # We grab high FC genes (positive FC) - if there are at least 20 genes, we can do FMD
        high_genes <- cell_specific_de_result[cell_specific_de_result$sc_log2FC > 0.1,]$Gene_Name
        # We also grab highly negative FC genes (opposite direction) - if there are at least 20 genes, we can do FMD
        low_genes <- cell_specific_de_result[cell_specific_de_result$sc_log2FC < -0.1,]$Gene_Name
        # Run FMD for each network
        for(network in hb_networks) {
          # Run FMD for high genes
          FMD_result_high <- run_fmd_wrapper(high_genes, network, output_dir, current_cell_type, unique(current_de_result$metadata_attribute), "high", log_flag = log_flag)
          if(!is.null(FMD_result_high)) {
            FMD_results[[index]] <- FMD_result_high
            index <- index + 1
          }
          # Run FMD for low genes
          FMD_result_low <- run_fmd_wrapper(low_genes, network, output_dir, current_cell_type, unique(current_de_result$metadata_attribute), "low", log_flag = log_flag)
          if(!is.null(FMD_result_low)) {
            FMD_results[[index]] <- FMD_result_low
            index <- index + 1
          }
        }
      }
    }
  }
  return(list(differential_expression_results, FMD_results))
}



#' Perform differential expression analysis (RNA)
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param metadata_df Data frame containing metadata for samples. Rownames should be sample names and column names should be metadata attributes with two classes (e.g., condition: disease and control)
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]). Directory will be created if it doesn't already exist.
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A list containing differential expression analyses
#' @examples
#' \dontrun{differential_expression_results <- RunDE_RNA(sc_obj = sc_obj, metadata_df = metadata_df)}
#' @export
RunDE_RNA <- function(sc_obj, metadata_df, output_dir = getwd(), log_flag = FALSE) {
  de_results <- list()
  index <- 1
  # Normalize paths (in case user provides relative paths)
  output_dir <- normalize_dir_path(output_dir)
  print_SPEEDI("\n", log_flag, silence_time = TRUE)
  print_SPEEDI("Running cell-type level differential expression analysis (RNA)", log_flag)
  for(metadata_attribute in colnames(metadata_df)) {
    print_SPEEDI(paste0("Currently analyzing metadata attribute: ", metadata_attribute), log_flag)
    print_SPEEDI(paste0("First metadata attribute value: ", unique(sc_obj[[metadata_attribute]])[,1][1]), log_flag)
    print_SPEEDI(paste0("Second metadata attribute (base level) value: ", unique(sc_obj[[metadata_attribute]])[,1][2]), log_flag)
    final_current_de <- data.frame(Cell_Type = character(), Gene_Name = character(), sc_pval_adj = character(), sc_log2FC = character(), pseudo_bulk_pval = character(),
                                   pseudo_bulk_log2FC = character())
    for(current_cell_type in unique(sc_obj$predicted_celltype_majority_vote)) {
      print_SPEEDI(paste0("Currently analyzing cell type: ", current_cell_type), log_flag)
      # Run FindMarkers to find DEGs
      idxPass <- which(sc_obj$predicted_celltype_majority_vote %in% current_cell_type)
      cellsPass <- names(sc_obj$orig.ident[idxPass])
      # Dummy declaration to avoid check() complaining
      cell_name <- NULL
      sc_obj_cell_type_subset <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
      Seurat::DefaultAssay(sc_obj_cell_type_subset) <- "SCT"
      Seurat::Idents(sc_obj_cell_type_subset) <- metadata_attribute
      current_de <- Seurat::FindMarkers(sc_obj_cell_type_subset, ident.1 = unique(sc_obj_cell_type_subset[[metadata_attribute]])[,1][1], ident.2 = unique(sc_obj_cell_type_subset[[metadata_attribute]])[,1][2],
                                logfc.threshold = 0.1, min.pct = 0.1, assay = "SCT", recorrect_umi = FALSE)
      current_de <- current_de[current_de$p_val_adj < 0.05,]
      # Run DESeq2 for pseudobulk filtering
      Seurat::DefaultAssay(sc_obj_cell_type_subset) <- "RNA"
      pseudobulk_counts <- create_pseudobulk_counts(sc_obj_cell_type_subset, log_flag)
      pseudobulk_metadata <- metadata_df
      pseudobulk_metadata$aliquots <- rownames(pseudobulk_metadata)
      pseudobulk_metadata <- pseudobulk_metadata[match(colnames(pseudobulk_counts), pseudobulk_metadata$aliquots),]
      # Dummy declaration to avoid check() complaining
      aliquots <- NULL
      pseudobulk_metadata <- subset(pseudobulk_metadata, select = -c(aliquots))
      pseudobulk_analysis <- DESeq2::DESeqDataSetFromMatrix(countData = pseudobulk_counts, colData = pseudobulk_metadata, design = stats::formula(paste("~",metadata_attribute)))
      pseudobulk_analysis <- DESeq2::DESeq(pseudobulk_analysis)
      pseudobulk_analysis_results_contrast <- utils::tail(DESeq2::resultsNames(pseudobulk_analysis), n=1)
      pseudobulk_analysis_results <- DESeq2::results(pseudobulk_analysis, name=pseudobulk_analysis_results_contrast)
      pseudobulk_analysis_results <- pseudobulk_analysis_results[rowSums(is.na(pseudobulk_analysis_results)) == 0, ] # Remove NAs
      pseudobulk_analysis_results <- pseudobulk_analysis_results[pseudobulk_analysis_results$pvalue < 0.05,]
      pseudobulk_analysis_results <- pseudobulk_analysis_results[pseudobulk_analysis_results$log2FoldChange < -0.3 | pseudobulk_analysis_results$log2FoldChange > 0.3,]
      # Filter genes from single cell based on DESeq2 pseudobulk results
      final_genes <- intersect(rownames(current_de), rownames(pseudobulk_analysis_results))
      # Record information about remaining genes in final_current_de
      for(current_gene in final_genes) {
        current_sc_pval_adj <- current_de[rownames(current_de) == current_gene,]$p_val_adj
        current_sc_log2FC <- current_de[rownames(current_de) == current_gene,]$avg_log2FC
        current_pseudo_bulk_pval <- pseudobulk_analysis_results[rownames(pseudobulk_analysis_results) == current_gene,]$pvalue
        current_pseudo_bulk_log2FC <- pseudobulk_analysis_results[rownames(pseudobulk_analysis_results) == current_gene,]$log2FoldChange
        current_pseudo_bulk_log2FC <- current_pseudo_bulk_log2FC * -1 # Pseudobulk test has condition 1 / condition 2 flipped relative to single cell test, so we just multiply by -1 so FC is in consistent direction
        current_row <- data.frame(current_cell_type, current_gene, current_sc_pval_adj, current_sc_log2FC, current_pseudo_bulk_pval, current_pseudo_bulk_log2FC)
        names(current_row) <- c("Cell_Type", "Gene_Name", "sc_pval_adj", "sc_log2FC", "pseudo_bulk_pval", "pseudo_bulk_log2FC")
        final_current_de <- rbind(final_current_de, current_row)
      }
    }
    print_SPEEDI(paste0("Writing results for metadata attribute ", metadata_attribute, "to file"), log_flag)
    metadata_attribute_no_spaces <- sub(" ", "_", metadata_attribute) # Remove spaces for file name
    utils::write.table(final_current_de, file = paste0(output_dir, metadata_attribute_no_spaces, ".DE.tsv"), sep = "\t", quote = FALSE)
    final_current_de$metadata_attribute <- metadata_attribute
    de_results[[index]] <- final_current_de
    index <- index + 1
  }
  Seurat::DefaultAssay(sc_obj) <- "SCT"
  print_SPEEDI("Differential expression analysis complete", log_flag)
  gc()
  return(de_results)
}

#' Perform functional module discovery (RNA)
#'
#' @param gene_list List of genes
#' @param network Background network
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A list (first element containing URL to see full results in web browser and second element containing enrichment results)
#' @examples
#' \dontrun{fmd_results <- RunFMD_RNA(gene_list = c("ADNP", "CHD8"), network = "blood")}
#' @export
RunFMD_RNA <- function(gene_list, network = "global", log_flag = FALSE) {
  print_SPEEDI("\n", log_flag, silence_time = TRUE)
  print_SPEEDI("Running functional module discovery analysis (RNA)", log_flag)
  print_SPEEDI(paste0("Length of initial gene list: ", length(gene_list)), log_flag)
  print_SPEEDI(paste0("Initial gene list: ", paste(gene_list, collapse = ' ')), log_flag)
  print_SPEEDI(paste0("Background network: ", network), log_flag)
  # HumanBase expects gene list concatenated with +
  concatenated_gene_list <- paste(gene_list, collapse = '+')
  # Grab relevant HumanBase URLs
  gene_check_url <- get_gene_check_url()
  fmd_submission_url <- get_fmd_submission_url()

  # Create payload for checking genes against genes present in HumanBase
  gene_check_payload <- list(
    q = concatenated_gene_list
  )
  # Submit POST request to gene_check_url with our list of genes
  gene_check_post_request = httr::POST(gene_check_url, body = gene_check_payload, encode = "form", httr::verbose())
  gene_check_post_response <- jsonlite::fromJSON(httr::content(gene_check_post_request, as = "text"))
  # We can parse the response to grab our final list of entrez IDs that we will use for our FMD job
  final_gene_list <- c()
  final_gene_list_entrez <- c()
  for(i in 1:nrow(gene_check_post_response)) {
    current_gene_row <- gene_check_post_response[i,]
    if(current_gene_row$num_matches > 0) {
      final_gene_list <- c(final_gene_list, current_gene_row$query)
      current_entrez_id <- current_gene_row$matches
      final_gene_list_entrez <- c(final_gene_list_entrez, current_gene_row$matches[[1]][1,]$entrez)
    }
  }
  final_gene_list_entrez <- as.character(final_gene_list_entrez)
  print_SPEEDI(paste0("Length of final gene list for HumanBase: ", length(final_gene_list)), log_flag)
  print_SPEEDI(paste0("Final gene list: ", paste(final_gene_list, collapse = ' ')), log_flag)

  if(length(final_gene_list) < 20) {
    print_SPEEDI("Fewer than 20 input genes were present in HumanBase, so we cannot run FMD", log_flag)
    gc()
    return(NULL)
  }

  # Each FMD job has an associated cache key (generated via SHA1) in the URL so we can refer to the results later
  # This cache key is based on the sorted gene list (entrez) and the chosen network
  # We generate the cache key below using the same algorithm as HumanBase (that way, you'll get the
  # same cache key regardless of whether you use the UI or this script)
  final_gene_list_entrez <- sort(final_gene_list_entrez)
  fmd_payload <- list(
    entrez = final_gene_list_entrez
  )
  fmd_hash <- as.character(jsonlite::toJSON(fmd_payload))
  fmd_hash <- paste(substr(fmd_hash, 1, 1), " ", substr(fmd_hash, 2, nchar(fmd_hash)), sep = "")
  fmd_hash <- paste(substr(fmd_hash, 1, 11), " ", substr(fmd_hash, 12, nchar(fmd_hash)), sep = "")
  fmd_hash <- paste(substr(fmd_hash, 1, nchar(fmd_hash) - 1), " ", substr(fmd_hash, nchar(fmd_hash), nchar(fmd_hash)), sep = "")
  fmd_hash <- paste0(fmd_hash, network)
  fmd_hash <- digest::digest(fmd_hash, algo = "sha1", serialize = F)
  # Create final FMD submission URL (based on network and hash)
  current_fmd_url <- paste0(fmd_submission_url, network, "&body_tag=", fmd_hash)
  # Submit FMD POST request (payload is our sorted list of entrez IDs)
  fmd_submission_post_request <- httr::POST(current_fmd_url, body = fmd_payload, encode = "json", httr::verbose())
  fmd_submission_post_response <- jsonlite::fromJSON(httr::content(fmd_submission_post_request, as = "text"))
  if(length(fmd_submission_post_response$enrichment) == 0) {
    print_SPEEDI("No enrichment found for input genes", log_flag)
    gc()
    return(NULL)
  }
  cached_url = paste0("https://hb.flatironinstitute.org/module/overview/?body_tag=", fmd_hash)
  print_SPEEDI(paste0("Done submitting FMD job! Associated URL is: ", cached_url), log_flag)
  print_SPEEDI("Functional module discovery analysis complete", log_flag)
  gc()
  return(list(cached_url, fmd_submission_post_response))
}

#' Grab associated HumanBase networks for current cell type and reference tissue
#'
#' @param current_cell_type Current cell type
#' @param reference_tissue Reference tissue
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A vector containing a list of HumanBase networks
#' @examples
#' \dontrun{hb_networks <- grab_hb_networks(current_cell_type = "B naive",
#' reference_tissue = "pbmc")}
#' @export
grab_hb_networks <- function(current_cell_type, reference_tissue, log_flag = FALSE) {
  print_SPEEDI(paste0("Grabbing HumanBase networks associated with cell type ", current_cell_type, " and reference tissue ", reference_tissue), log_flag)
  hb_networks <- c("global")
  # Change reference_tissue to all lowercase to prevent any issues with casing
  reference_tissue <- tolower(reference_tissue)
  # TODO: Add cell-type granularity to PBMC networks
  # TODO: Look into cortex (cerebral cortex tissue OK?)
  if(reference_tissue == "pbmc" | reference_tissue == "pbmc_full") {
    hb_networks <- c(hb_networks, "blood")
  } else if (reference_tissue == "adipose") {
    hb_networks <- c(hb_networks, "adipose-tissue")
  } else if (reference_tissue == "bone_marrow") {
    hb_networks <- c(hb_networks, "bone-marrow")
  } else if (reference_tissue == "fetus") {
    hb_networks <- c(hb_networks, "fetus")
  } else if (reference_tissue == "heart") {
    hb_networks <- c(hb_networks, "heart")
  } else if (reference_tissue == "kidney") {
    hb_networks <- c(hb_networks, "kidney")
  } else if (reference_tissue == "lung") {
    hb_networks <- c(hb_networks, "lung")
  } else if (reference_tissue == "pancreas") {
    hb_networks <- c(hb_networks, "pancreas")
  }else if (reference_tissue == "tonsil") {
    hb_networks <- c(hb_networks, "tonsil")
  }
  return(hb_networks)
}

#' Wrapper for running FMD job in the context of the SPEEDI wrapper function
#'
#' @param gene_list Input list of genes
#' @param network Background network
#' @param RNA_output_dir Output directory (RNA)
#' @param cell_type Cell type
#' @param metadata_attribute Metadata attribute
#' @param fc_flag Flag to indicate direction of fold change ("high" or "low")
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A list containing FMD results
#' @examples
#' \dontrun{hb_networks <- grab_hb_networks(current_cell_type = "B naive",
#' reference_tissue = "pbmc")}
#' @export
run_fmd_wrapper <- function(gene_list, network, RNA_output_dir, cell_type, metadata_attribute, fc_flag, log_flag = FALSE) {
  print_SPEEDI(paste0("Running FMD for ", fc_flag, " genes in cell type ", cell_type, " for metadata attribute ", metadata_attribute), log_flag)
  FMD_result <- NULL
  if(length(gene_list) >= 20) {
    # Run FMD and create output file where header line (starting with #) is a URL to see full results in web browser
    # The table below contains enrichment results from HumanBase
    FMD_result <- RunFMD_RNA(gene_list = gene_list, network = network, log_flag = log_flag)
    if(!is.null(FMD_result)) {
      print_SPEEDI("Writing FMD results to file", log_flag)
      cell_type <- sub(" ", "_", cell_type) # Remove spaces for file name
      output_file <- paste0(RNA_output_dir, "FMD_", fc_flag, "_", cell_type, "_", network, "_", metadata_attribute, ".csv")
      cat(paste0("# ", FMD_result[[1]], "\n"), file=output_file)
      current_enrichment_table <- FMD_result[[2]]$enrichment
      current_enrichment_table <- current_enrichment_table[,!names(current_enrichment_table) %in% c("genes", "edges")]
      utils::write.table(current_enrichment_table, file = output_file, append=TRUE, sep = ",", quote = FALSE)
    }
  } else {
    print_SPEEDI("Not enough genes in gene list for FMD", log_flag)
  }
  return(FMD_result)
}

#' Create pseudobulk counts
#' @param sc_obj Seurat object Input list of genes
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @examples
#' \dontrun{pseudobulk_counts <- create_pseudobulk_counts(sc_obj)}
#' @export
create_pseudobulk_counts <- function(sc_obj, log_flag) {
  print_SPEEDI("Computing pseudobulk counts", log_flag)
  cells_pseudobulk <- list()
  for (sample_name in unique(sc_obj$sample)) {
    idxMatch <- which(stringr::str_detect(as.character(sc_obj$sample), sample_name))
    if(length(idxMatch)>=1) {
      samples_subset <- subset(x = sc_obj, subset = sample %in% sample_name)
      samples_data <- samples_subset@assays$RNA@counts
      samples_data <- rowSums(as.matrix(samples_data))
      cells_pseudobulk[[sample_name]] <- samples_data
    } else {
      cells_pseudobulk[[sample_name]] <- numeric(nrow(sc_obj@assays$RNA))
    }
  }
  final_cells_pseudobulk_df <- dplyr::bind_cols(cells_pseudobulk[1])
  for (idx in 2:length(unique(sc_obj$sample))) {
    final_cells_pseudobulk_df <- dplyr::bind_cols(final_cells_pseudobulk_df, cells_pseudobulk[idx])
  }
  final_cells_pseudobulk_df <- as.data.frame(final_cells_pseudobulk_df)
  rownames(final_cells_pseudobulk_df) <- names(cells_pseudobulk[[1]])
  return(final_cells_pseudobulk_df)
}
