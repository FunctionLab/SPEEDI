#' Perform differential expression analysis (RNA)
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param metadata_df Dataframe containing metadata for samples. Rownames should be sample names and column names should be metadata attributes with two classes (e.g., condition: disease and control)
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]). Directory will be created if it doesn't already exist.
#' @param de_family Differential expression family (see [Libra::run_de()]).
#' @param de_method Differential expression method (see [Libra::run_de()]).
#' @param de_type Differential expression type (see [Libra::run_de()]).
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A list containing differential expression analyses
#' @examples
#' \dontrun{differential_expression_results <- RunDE_RNA(sc_obj = sc_obj, metadata_df = metadata_df)}
#' @export
RunDE_RNA <- function(sc_obj, metadata_df, de_family = "pseudobulk", de_method = "edgeR", de_type = "LRT", output_dir = getwd(), log_flag = FALSE) {
  de_results <- list()
  index <- 1
  # Normalize paths (in case user provides relative paths)
  output_dir <- normalize_dir_path(output_dir)
  print_SPEEDI("\n", log_flag, silence_time = TRUE)
  print_SPEEDI("Running differential expression analysis (RNA)", log_flag)
  Seurat::DefaultAssay(sc_obj) <- "RNA"
  for(metadata_attribute in colnames(metadata_df)) {
    current_de <- Libra::run_de(sc_obj, replicate_col = "sample",
                         cell_type_col = "predicted_celltype_majority_vote", label_col = metadata_attribute,
                         de_family = de_family, de_method = de_method, de_type = de_type)
    current_de$metadata_attribute <- metadata_attribute
    utils::write.table(current_de, file = paste0(output_dir, metadata_attribute, ".DE.tsv"), sep = "\t", quote = FALSE)
    de_results[[index]] <- current_de
    index <- index + 1
  }
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
      final_gene_list_entrez <- c(final_gene_list_entrez, current_gene_row$matches[[1]]$entrez)
    }
  }
  final_gene_list_entrez <- as.character(final_gene_list_entrez)
  print_SPEEDI(paste0("Length of final gene list for HumanBase: ", length(final_gene_list)), log_flag)
  print_SPEEDI(paste0("Final gene list: ", paste(final_gene_list, collapse = ' ')), log_flag)

  if(length(final_gene_list) < 20) {
    print_SPEEDI("Fewer than 20 input genes were present in HumanBase, so we cannot run FMD", log_flag)
    gc()
    return(FALSE)
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
  cached_url = paste0("https://hb.flatironinstitute.org/module/overview/?body_tag=", fmd_hash)
  print_SPEEDI(paste0("Done submitting FMD job! Associated URL is: ", cached_url), log_flag)
  print_SPEEDI("Functional module discovery analysis complete", log_flag)
  gc()
  return(list(cached_url, fmd_submission_post_response))
}
