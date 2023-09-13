#'
NULL

#' downsample_seurat: downsample the Seurat object to the minimum number of cells per
#' metadata category (sample, batch etc.)
#' 
#' @param seurat_obj A Seurat object
#' @param metadata_slot The name of the metadata slot 
#' @param target_cells The minimum number of cells per sample
#' @return A Seurat object with the minimum number of cells per sample
#' 
#' @export
#' 
downsample_seurat <- function(seurat_obj, metadata_slot, target_cells) {
  # Get the unique sample IDs from the metadata
  ids <- as.character(unique(seurat_obj@meta.data[[metadata_slot]]))
  
  # Define a helper function to downsample cells for a given sample
  downsample_cells <- function(sam) {
    # Get the cells for this sample
    cells <- Cells(seurat_obj)[which(seurat_obj[[metadata_slot]] == sam)]
    if (length(cells) > target_cells) {
      # If there are more cells than the target, downsample
      keep_cells <- sample(cells, target_cells)
    } else {
      # If there are fewer cells than the target, keep all
      keep_cells <- cells
    }
    return(keep_cells)
  }
  
  # Use lapply to iterate over sample IDs and downsample cells
  cells_to_keep <- unlist(lapply(ids, downsample_cells))
  
  # Subset the seurat_obj using the cells to keep
  seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
  
  print("Downsampled Seurat object to minimum number of cells per sample, cells per sample:")
  print(table(seurat_obj@meta.data[[metadata_slot]]))
  
  return(seurat_obj)
}

#' refine_metadata_levels: removes unused metadata levels from a Seurat object,
#' useful for removing unused levels from the metadata after subsetting. Credit to
#' github user michael-kotliar for this function 
#' (https://github.com/satijalab/seurat/issues/5069#issuecomment-1372599846).
#' 
#' @param seurat_obj A Seurat object
#' @return A Seurat object with unused metadata levels removed
#' 
#' @export
#' 
refine_metadata_levels <- function(seurat_data){
    for (i in base::colnames(seurat_data@meta.data)){
        if (base::is.factor(seurat_data@meta.data[[i]])){
            base::print(base::paste("Re-evaluating levels for a factor column", i))
            base::print(
                base::paste(
                    "before:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
                )
            )
            seurat_data@meta.data[[i]] <- base::droplevels(seurat_data@meta.data[[i]])  # need to drop levels of the removed values
            base::print(
                base::paste(
                    "after:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
                )
            )
        }
    }
    return (seurat_data)
}

#' plot_metadata: Plots all of the categorical metadata in grouped
#' and split UMAPs, autosizes based on the number of unique values
#' in the metadata. Plots are output as PDFs to a directory.
#' 
#' @param seurat_obj A Seurat object
#' @param output_dir Output directory
#' @return Writes out plots to directory
#' 
#' @export
#' 
plot_metadata <- function(seurat_obj, output_dir) {
  ## Plots all of the categorical metadata in grouped and split UMAPs,
  ## autosizes based on the number of unique values in the metadata

  # get the metadata columns that not numeric
  meta_cols <- names(seurat_obj@meta.data)[!sapply(seurat_obj@meta.data, is.numeric)]
  # Exclude "orig.ident"
  meta_cols <- meta_cols[meta_cols != "orig.ident"]
  for (meta in meta_cols) {
    print(paste0("***** Plotting ", meta, " *****"))
    pdf(file.path(paste0(output_dir, "/", meta, "_UMAP.pdf")),
      width = 10,
      height = 10
    )
    print(DimPlot(seurat_obj,
      raster = FALSE,
      order = TRUE,
      label = FALSE,
      group.by = meta
    ))
    dev.off()
    h <- ifelse(length(unique(seurat_obj@meta.data[[meta]])) == 2,
      10,
      3 * ceiling(length(unique(seurat_obj@meta.data[[meta]])) / 2)
    )
    pdf(file.path(paste0(output_dir, "/", meta, "_split_UMAP.pdf")),
      width = 15,
      height = h
    )
    print(DimPlot(seurat_obj,
      raster = FALSE,
      order = TRUE,
      label = FALSE,
      group.by = meta,
      split.by = meta,
      ncol = 2
    ) + NoLegend())
    dev.off()
  }
}
#' plot_metadata_numeric: Plots all of the numeric metadata
#' using feature plot UMAPs. Plots are output as PDFs to a 
#' directory.
#' 
#' @param seurat_obj A Seurat object
#' @param output_dir Output directory
#' @return Writes out plots to directory
#' 
#' @export
#' 
plot_metadata_numeric <- function(seurat_obj, output_dir) {
  # Plot UMAPs for all numeric metadata columns
  nonchar_cols <- names(seurat_obj@meta.data)[sapply(seurat_obj@meta.data, is.numeric)]
  for (meta in nonchar_cols) {
    print(paste0("***** Plotting ", meta, " *****"))
    pdf(file.path(paste0(output_dir, "/", meta, "_UMAP.pdf")),
      width = 10,
      height = 10
    )
    print(FeaturePlot(seurat_obj,
      features = meta,
      raster = FALSE,
      order = TRUE
    ))
    dev.off()
  }
}