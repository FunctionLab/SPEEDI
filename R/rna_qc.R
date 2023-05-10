#' Filter raw data for RNA
#'
#' @param sc_obj Seurat object containing cells for all samples (pre-filtering)
#' @param output_dir Directory where QC output will be saved
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A Seurat object which contains filtered data from all input data
#' @examples
#' \dontrun{sc_obj <- Create_RNA_QC_Output(sc_obj = sc_obj, output_dir = getwd()}
#' @export
Create_RNA_QC_Output <- function(sc_obj, output_dir = getwd(), log_flag = FALSE) {
  print_SPEEDI("Printing QC plots on pre-filtered data (RNA)", log_flag)
  p <- Seurat::VlnPlot(sc_obj, features = c("nFeature_RNA"), split.by = "sample", group.by = "sample", raster = FALSE)
  ggplot2::ggsave(paste0(output_dir, "nFeature_violin_plot.png"), plot = p, device = "png", width = 10, height = 10, units = "in")
  p <- Seurat::VlnPlot(sc_obj, features = c("nCount_RNA"), split.by = "sample", group.by = "sample", raster = FALSE)
  ggplot2::ggsave(paste0(output_dir, "nCount_violin_plots.png"), plot = p, device = "png", width = 10, height = 10, units = "in")
  p <- Seurat::VlnPlot(sc_obj, features = c("percent.mt"), split.by = "sample", group.by = "sample", raster = FALSE)
  ggplot2::ggsave(paste0(output_dir, "percentMT_violin_plots.png"), plot = p, device = "png", width = 10, height = 10, units = "in")
  p <- Seurat::VlnPlot(sc_obj, features = c("percent.hb"), split.by = "sample", group.by = "sample", raster = FALSE)
  ggplot2::ggsave(paste0(output_dir, "percentHB_violin_plots.png"), plot = p, device = "png", width = 10, height = 10, units = "in")
  p <- Seurat::VlnPlot(sc_obj, features = c("percent.rp"), split.by = "sample", group.by = "sample", raster = FALSE)
  ggplot2::ggsave(paste0(output_dir, "percentRP_violin_plots.png"), plot = p, device = "png", width = 10, height = 10, units = "in")
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
  ggplot2::ggsave(paste0(output_dir, "nCount_vs_nFeature_plots.png"), plot = nCount_vs_nFeature_plots, device = "png", width = 20, height = 20, units = "in")
  rm(individual_samples)
  gc()
  return(sc_obj)
}
