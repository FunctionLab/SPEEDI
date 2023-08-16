#' Create QC output files for prefiltered data (RNA)
#'
#' @param sc_obj Seurat object containing cells for all samples (pre-filtering)
#' @param output_dir Directory where QC output will be saved
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return TRUE
#' @examples
#' \dontrun{Create_QC_Output_Prefiltered_RNA(sc_obj = sc_obj, output_dir = getwd()}
#' @export
Create_QC_Output_Prefiltered_RNA <- function(sc_obj, output_dir = getwd(), log_flag = FALSE) {
  # Normalize paths (in case user provides relative paths)
  output_dir <- normalize_dir_path(output_dir)
  print_SPEEDI("Printing QC plots on pre-filtered data (RNA)", log_flag)
  p <- Seurat::VlnPlot(sc_obj, features = c("nFeature_RNA"), split.by = "sample", group.by = "sample", raster = FALSE) + ggplot2::xlab("Sample")
  ggplot2::ggsave(paste0(output_dir, "pre-filtered_nFeature_violin_plot.png"), plot = p, device = "png", width = 10, height = 10, units = "in")
  p <- Seurat::VlnPlot(sc_obj, features = c("nCount_RNA"), split.by = "sample", group.by = "sample", raster = FALSE) + ggplot2::xlab("Sample")
  ggplot2::ggsave(paste0(output_dir, "pre-filtered_nCount_violin_plots.png"), plot = p, device = "png", width = 10, height = 10, units = "in")
  p <- Seurat::VlnPlot(sc_obj, features = c("percent.mt"), split.by = "sample", group.by = "sample", raster = FALSE) + ggplot2::xlab("Sample")
  ggplot2::ggsave(paste0(output_dir, "pre-filtered_percentMT_violin_plots.png"), plot = p, device = "png", width = 10, height = 10, units = "in")
  p <- Seurat::VlnPlot(sc_obj, features = c("percent.hb"), split.by = "sample", group.by = "sample", raster = FALSE) + ggplot2::xlab("Sample")
  ggplot2::ggsave(paste0(output_dir, "pre-filtered_percentHB_violin_plots.png"), plot = p, device = "png", width = 10, height = 10, units = "in")
  p <- Seurat::VlnPlot(sc_obj, features = c("percent.rp"), split.by = "sample", group.by = "sample", raster = FALSE) + ggplot2::xlab("Sample")
  ggplot2::ggsave(paste0(output_dir, "pre-filtered_percentRP_violin_plots.png"), plot = p, device = "png", width = 10, height = 10, units = "in")
  # Create nCount vs nFeature scatter plot for each sample
  individual_samples <- Seurat::SplitObject(sc_obj, split.by = "sample")
  # Visualize nCount vs nFeature via scatter plot for each sample
  nCount_vs_nFeature_plots <- list()
  for (i in 1:length(individual_samples)) {
    p <- Seurat::FeatureScatter(individual_samples[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")
    nCount_vs_nFeature_plots[[i]] <- p
  }
  n <- length(nCount_vs_nFeature_plots)
  nCol <- floor(sqrt(n))
  nCount_vs_nFeature_plots <- do.call(parse_namespace_and_function("gridExtra::grid.arrange"), c(nCount_vs_nFeature_plots, ncol=nCol))
  ggplot2::ggsave(paste0(output_dir, "pre-filtered_nCount_vs_nFeature_plots.png"), plot = nCount_vs_nFeature_plots, device = "png", width = 20, height = 20, units = "in")
  rm(individual_samples)
  gc()
  return(TRUE)
}

#' Create QC output files for prefiltered data (ATAC)
#'
#' @param proj ArchR project associated with data
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]).
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return TRUE
#' @examples
#' \dontrun{Create_QC_Output_Prefiltered_ATAC(proj = proj}
#' @export
Create_QC_Output_Prefiltered_ATAC <- function(proj, output_dir = getwd(), log_flag = FALSE) {
  num_cells <- length(proj$cellNames)
  num_samples <- length(unique(proj$Sample))
  sample_text <- ""
  if(num_samples == 1) {
    sample_text <- paste0("(1 Sample, ", num_cells, " Cells)")
  } else {
    sample_text <- paste0("(", num_samples, " Samples, ", num_cells, " Cells)")
  }
  print_SPEEDI("Printing QC plots on pre-filtered data (ATAC)", log_flag)
  # Plot out TSS Enrichment / Doublet Enrichment / Nucleosome Ratio for each sample to help us decide filtering thresholds
  p1 <- ArchR::plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "ridges") +
    ggplot2::ggtitle(paste0("TSS Enrichment\n", sample_text)) + ggplot2::theme(plot.title = ggplot2::element_text(size=18), axis.text=ggplot2::element_text(size=10), axis.title=ggplot2::element_text(size=14))
  ggplot2::ggsave(filename = paste0(output_dir, "pre-filtered_TSSEnrichment.png"), plot = p1, device = "png", width = 8, height = 8, units = "in")
  p2 <- ArchR::plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "DoubletEnrichment", plotAs = "ridges") +
    ggplot2::ggtitle(paste0("Doublet Enrichment\n", sample_text)) + ggplot2::theme(plot.title = ggplot2::element_text(size=18), axis.text=ggplot2::element_text(size=10), axis.title=ggplot2::element_text(size=14))
  ggplot2::ggsave(filename = paste0(output_dir, "pre-filtered_DoubletEnrichment.png"), plot = p2, device = "png", width = 8, height = 8, units = "in")
  p3 <- ArchR::plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "NucleosomeRatio", plotAs = "ridges") +
    ggplot2::ggtitle(paste0("Nucleosome Ratio\n", sample_text)) + ggplot2::theme(plot.title = ggplot2::element_text(size=18), axis.text=ggplot2::element_text(size=10), axis.title=ggplot2::element_text(size=14))
  ggplot2::ggsave(filename = paste0(output_dir, "pre-filtered_NucleosomeRatio.png"), plot = p3, device = "png", width = 8, height = 8, units = "in")
  ArchR::plotPDF(p1,p2,p3, name = "pre-filtered_QC_metrics_plots", ArchRProj = proj, addDOC = FALSE, width = 7, height = 5)
  gc()
  return(TRUE)
}

