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
  print_SPEEDI("Printing QC plots on pre-filtered data (RNA)", log_flag)
  p <- Seurat::VlnPlot(sc_obj, features = c("nFeature_RNA"), split.by = "sample", group.by = "sample", raster = FALSE)
  ggplot2::ggsave(paste0(output_dir, "pre-filtered_nFeature_violin_plot.png"), plot = p, device = "png", width = 10, height = 10, units = "in")
  p <- Seurat::VlnPlot(sc_obj, features = c("nCount_RNA"), split.by = "sample", group.by = "sample", raster = FALSE)
  ggplot2::ggsave(paste0(output_dir, "pre-filtered_nCount_violin_plots.png"), plot = p, device = "png", width = 10, height = 10, units = "in")
  p <- Seurat::VlnPlot(sc_obj, features = c("percent.mt"), split.by = "sample", group.by = "sample", raster = FALSE)
  ggplot2::ggsave(paste0(output_dir, "pre-filtered_percentMT_violin_plots.png"), plot = p, device = "png", width = 10, height = 10, units = "in")
  p <- Seurat::VlnPlot(sc_obj, features = c("percent.hb"), split.by = "sample", group.by = "sample", raster = FALSE)
  ggplot2::ggsave(paste0(output_dir, "pre-filtered_percentHB_violin_plots.png"), plot = p, device = "png", width = 10, height = 10, units = "in")
  p <- Seurat::VlnPlot(sc_obj, features = c("percent.rp"), split.by = "sample", group.by = "sample", raster = FALSE)
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
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return TRUE
#' @examples
#' \dontrun{Create_QC_Output_Prefiltered_ATAC(proj = proj}
#' @export
Create_QC_Output_Prefiltered_ATAC <- function(proj, log_flag = FALSE) {
  print_SPEEDI("Printing QC plots on pre-filtered data (ATAC)", log_flag)
  # Plot out TSS Enrichment / Doublet Enrichment / Nucleosome Ratio for each sample to help us decide filtering thresholds
  p1 <- ArchR::plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "ridges")
  p2 <- ArchR::plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "DoubletEnrichment", plotAs = "ridges")
  p3 <- ArchR::plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "NucleosomeRatio", plotAs = "ridges")
  ArchR::plotPDF(p1,p2,p3, name = "pre-filtered_QC_metrics.pdf", ArchRProj = proj, addDOC = FALSE, width = 7, height = 5)
  gc()
  return(TRUE)
}

